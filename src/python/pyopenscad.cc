/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <Python.h>
#include "pyopenscad.h"
#include "CsgOpNode.h"
#include "Value.h"
#include "Expression.h"
#include "PlatformUtils.h"
#include <Context.h>
#include <Selection.h>
#include "PlatformUtils.h"

// #define HAVE_PYTHON_YIELD
static PyObject *PyInit_openscad(void);

// https://docs.python.org/3.10/extending/newtypes.html

void PyObjectDeleter (PyObject *pObject) { Py_XDECREF(pObject); };

PyObjectUniquePtr pythonInitDict(nullptr, PyObjectDeleter) ;
PyObjectUniquePtr pythonMainModule(nullptr, PyObjectDeleter) ;
std::list<std::string> pythonInventory;
AssignmentList customizer_parameters;
AssignmentList customizer_parameters_finished;
bool pythonDryRun=false;
std::shared_ptr<AbstractNode> python_result_node = nullptr; /* global result veriable containing the python created result */
std::vector<SelectedObject> python_result_handle;
bool python_active;  /* if python is actually used during evaluation */
bool python_trusted; /* global Python trust flag */
bool python_runipython = false;
bool pythonMainModuleInitialized = false;
bool pythonRuntimeInitialized = false;

std::vector<std::string> mapping_name;
std::vector<std::string> mapping_code;
std::vector<int> mapping_level;
std::shared_ptr<const FileContext> osinclude_context = nullptr;


void PyOpenSCADObject_dealloc(PyOpenSCADObject *self)
{
  Py_XDECREF(self->dict);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyOpenSCADObject_alloc(PyTypeObject *cls, Py_ssize_t nitems)
{
  PyObject *self = PyType_GenericAlloc(cls, nitems);
  ((PyOpenSCADObject *)self)->dict = PyDict_New();
  PyObject *origin=PyList_New(4);
  for(int i=0;i<4;i++) {
  	PyObject *row=PyList_New(4);
	for(int j=0;j<4;j++)
		PyList_SetItem(row,j,PyFloat_FromDouble(i==j?1.0:0.0));
	PyList_SetItem(origin,i,row);
//  	Py_XDECREF(row);
  }
  PyDict_SetItemString(((PyOpenSCADObject *)self)->dict,"origin",origin);
  Py_XDECREF(origin);
  return self;
}

/*
 *  allocates a new PyOpenSCAD Object including its internal dictionary
 */

static PyObject *PyOpenSCADObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  return PyOpenSCADObject_alloc(&PyOpenSCADType, 0);
}

/*
 *  allocates a new PyOpenSCAD to store an existing OpenSCAD Abstract Node
 */

PyObject *PyOpenSCADObjectFromNode(PyTypeObject *type, const std::shared_ptr<AbstractNode> &node)
{
  PyOpenSCADObject *self;
  self = (PyOpenSCADObject *)  type->tp_alloc(type, 0);
  if (self != nullptr) {
    self->node = node;
    return (PyObject *)self;
  }
  return nullptr;
}

//PyGILState_STATE gstate=PyGILState_LOCKED;
PyThreadState *tstate=nullptr;

void python_lock(void){
//#ifndef _WIN32	
  if(tstate != nullptr && pythonInitDict != nullptr) PyEval_RestoreThread(tstate);
//#endif  
}

void python_unlock(void) {
//#ifndef _WIN32	
  if(pythonInitDict != nullptr)	tstate = PyEval_SaveThread();
//#endif  
}

/*
 *  parses either a PyOpenSCAD Object or an List of PyOpenScad Object and adds it to the list of supplied children, returns 1 on success
 */

int python_more_obj(std::vector<std::shared_ptr<AbstractNode>>& children, PyObject *more_obj) {
  int i, n;
  PyObject *obj;
  PyObject *dummy_dict;
  std::shared_ptr<AbstractNode> child;
  if (PyList_Check(more_obj)) {
    n = PyList_Size(more_obj);
    for (i = 0; i < n; i++) {

      obj = PyList_GetItem(more_obj, i);
      child = PyOpenSCADObjectToNode(obj, &dummy_dict);
      children.push_back(child);
    }
  } else if (Py_TYPE(more_obj) == &PyOpenSCADType) {
    child = PyOpenSCADObjectToNode(more_obj, &dummy_dict);
    children.push_back(child);
  } else return 1;
  return 0;
}

/*
 *  extracts Absrtract Node from PyOpenSCAD Object
 */

std::shared_ptr<AbstractNode> PyOpenSCADObjectToNode(PyObject *obj, PyObject **dict)
{
  std::shared_ptr<AbstractNode> result = ((PyOpenSCADObject *) obj)->node;
  if(result.use_count() > 2) {
    result = result->clone();
  }
  *dict =  ((PyOpenSCADObject *) obj)->dict;
  return result;
}


/*
 * same as  python_more_obj but always returns only one AbstractNode by creating an UNION operation
 */

std::shared_ptr<AbstractNode> PyOpenSCADObjectToNodeMulti(PyObject *objs,PyObject **dict)
{
  std::shared_ptr<AbstractNode> result;
  if (Py_TYPE(objs) == &PyOpenSCADType) {
    result = ((PyOpenSCADObject *) objs)->node;
    if(result.use_count() > 2) {
	    result = result->clone();
    }
    *dict =  ((PyOpenSCADObject *) objs)->dict;
  } else if (PyList_Check(objs)) {
    DECLARE_INSTANCE
    auto node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);

    int n = PyList_Size(objs);
    for (int i = 0; i < n; i++) {
      PyObject *obj = PyList_GetItem(objs, i);
      if(Py_TYPE(obj) ==  &PyOpenSCADType) {
        std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNode(obj,dict);
        node->children.push_back(child);
      } else return nullptr;
    }
    result=node;
    *dict = nullptr; // TODO improve
  } else result=nullptr;
  return result;
}

void  python_hierdump(std::ostringstream &stream, const std::shared_ptr<AbstractNode> &node)
{
  if(node == nullptr) return;	
  stream <<  node->toString();	
  auto children = node->getChildren();
  if(children.size() < 1) stream << ";";
  else if(children.size() < 2)  python_hierdump(stream, children[0]);
  else {
   stream <<  "{ ";	  
   for(unsigned int i=0;i<children.size();i++) python_hierdump(stream, children[i]);
   stream << " }";	  
  }
}
void python_build_hashmap(const std::shared_ptr<AbstractNode> &node, int level)
{
	
  PyObject *maindict = PyModule_GetDict(pythonMainModule.get());
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  std::ostringstream stream;
//  python_hierdump(stream, node);
  std::string code = stream.str();
  while (PyDict_Next(maindict, &pos, &key, &value)) {
    if(value->ob_type != &PyOpenSCADType) continue;
    std::shared_ptr<AbstractNode> testnode = ((PyOpenSCADObject *) value)->node;
    if(testnode != node) continue;
    PyObject* key1 = PyUnicode_AsEncodedString(key, "utf-8", "~");
    if(key1 == nullptr) continue;
    const char *key_str =  PyBytes_AS_STRING(key1);
    if(key_str == nullptr) continue;
    mapping_name.push_back(key_str);
    mapping_code.push_back(code);
    mapping_level.push_back(pos);
  }
  if(level < 5) { // no  many level are unclear and error prone(overwrites memory)
    for(const auto &child:node->getChildren()) {
      python_build_hashmap(child,level+1);	  
    }
  }  
}

void python_retrieve_pyname(const std::shared_ptr<AbstractNode> &node)
{
  std::string name;	
  int level=-1;
  std::ostringstream stream;
  python_hierdump(stream, node);
  std::string my_code = stream.str();
  for(unsigned int i=0;i<mapping_code.size();i++) {
    if(mapping_code[i] == my_code) {
      if(level == -1 || level > mapping_level[i]) {	    
        name=mapping_name[i];	    
        level=mapping_level[i];	    
      }	
    }	    
  }
  node->setPyName(name);
}
/*
 * converts a python obejct into an integer by all means
 */

int python_numberval(PyObject *number, double *result)
{
  if(number == nullptr) return 1;
  if(number == Py_False) return 1;
  if(number == Py_True) return 1;
  if(number == Py_None) return 1;
  if (PyFloat_Check(number)) {
    *result = PyFloat_AsDouble(number);
    return 0;
  }
  if (PyLong_Check(number)) {
    *result = PyLong_AsLong(number);
    return 0;
  }
  return 1;
}

std::vector<int>  python_intlistval(PyObject *list)
{
  std::vector<int> result;	
  PyObject *item;
  if( PyLong_Check(list)) {
    result.push_back(PyLong_AsLong(list));
  }
  if( PyList_Check(list)) {
    for(int i=0;i<PyList_Size(list); i++) {
      item = PyList_GetItem(list, i);	    
      if( PyLong_Check(item)) {
        result.push_back(PyLong_AsLong(item));
      }
    }	    
  }
  return result;
}
/*
 * Tries to extract an 3D vector out of a python list
 */

int python_vectorval(PyObject *vec, int minval, int maxval, double *x, double *y, double *z, double *w)
{
  if(w != NULL ) *w = 0;
  if (PyList_Check(vec)) {
    if(PyList_Size(vec) < minval || PyList_Size(vec) > maxval) return 1;
    	  
    if (PyList_Size(vec) >= 1) {
      if (python_numberval(PyList_GetItem(vec, 0), x)) return 1;
    }
    if (PyList_Size(vec) >= 2) {
      if (python_numberval(PyList_GetItem(vec, 1), y)) return 1;
    }
    if (PyList_Size(vec) >= 3) {
      if (python_numberval(PyList_GetItem(vec, 2), z)) return 1;
    }
    if (PyList_Size(vec) >= 4 && w != NULL) {
      if (python_numberval(PyList_GetItem(vec, 3), w)) return 1;
    }
    return 0;
  }
  if (!python_numberval(vec, x)) {
    *y = *x;
    *z = *x;
    if(w != NULL) *w = *x;
    return 0;
  }
  return 1;
}

std::vector<Vector3d> python_vectors(PyObject *vec, int mindim, int maxdim) 
{
  std::vector<Vector3d> results;	
  Vector3d result;	
  if (PyList_Check(vec)) {
    // check if its a valid vec<Vector3d>
    int valid=1;
    for(int i=0;valid && i<PyList_Size(vec);i++) {
      PyObject *item = PyList_GetItem(vec,i);
      if(!PyList_Check(item)) valid=0;
    }	    
    if(valid) {
      for(int j=0;valid && j<PyList_Size(vec);j++) {
        PyObject *item = PyList_GetItem(vec,j);
        if(PyList_Size(item) >= mindim && PyList_Size(item) <= maxdim) {	  
          for(int i=0;i<PyList_Size(item);i++) {
            if (PyList_Size(item) > i) {
              if (python_numberval(PyList_GetItem(item, i), &result[i])) return results; // Error
            }
          }	
        }  
	results.push_back(result);
      }	
      return results;
    }
    if(PyList_Size(vec) >= mindim && PyList_Size(vec) <= maxdim) {	  
      for(int i=0;i<PyList_Size(vec);i++) {
        if (PyList_Size(vec) > i) {
          if (python_numberval(PyList_GetItem(vec, i), &result[i])) return results; // Error
        }
      }	
    }  
    results.push_back(result);
  }
  if (!python_numberval(vec, &result[0])) {
    result[1] = result[0];
    result[2] = result[1];
    results.push_back(result);
  }
  return results; // Error
}

/*
 * Helper function to extract actual values for fn, fa and fs
 */

void get_fnas(double& fn, double& fa, double& fs) {
  PyObject *mainModule = PyImport_AddModule("__main__");
  if (mainModule == nullptr) return;
  PyObjectUniquePtr varFn(PyObject_GetAttrString(mainModule, "fn"),PyObjectDeleter);
  PyObjectUniquePtr varFa(PyObject_GetAttrString(mainModule, "fa"),PyObjectDeleter);
  PyObjectUniquePtr varFs(PyObject_GetAttrString(mainModule, "fs"),PyObjectDeleter);
  if (varFn.get() != nullptr){
    fn = PyFloat_AsDouble(varFn.get());
  }
  if (varFa.get() != nullptr){
    fa = PyFloat_AsDouble(varFa.get());
  }
  if (varFs.get() != nullptr){
    fs = PyFloat_AsDouble(varFs.get());
  }
}

/*
 * Type specific init function. nothing special here
 */

static int PyOpenSCADInit(PyOpenSCADObject *self, PyObject *arfs, PyObject *kwds)
{
  (void)self;
  (void)arfs;
  (void)kwds;
  return 0;
}
Outline2d python_getprofile(void *v_cbfunc, int fn, double arg)
{
	PyObject *cbfunc = (PyObject *) v_cbfunc;
	Outline2d result;
	if(pythonInitDict == NULL)  initPython(0.0);
	PyObject* args = PyTuple_Pack(1,PyFloat_FromDouble(arg));
	PyObject* polygon = PyObject_CallObject(cbfunc, args);
        Py_XDECREF(args);
	if(polygon == NULL) { // TODO fix
		for(unsigned int i=0;i < (unsigned int) fn;i++) {
			double ang=360.0*(i/(double) fn);
			PyObject* args = PyTuple_Pack(2,PyFloat_FromDouble(arg),PyFloat_FromDouble(ang));
			Py_XINCREF(args);
			PyObject* pypt = PyObject_CallObject(cbfunc, args);
			double r=PyFloat_AsDouble(pypt);
			if(r < 0) r=-r;  // TODO who the hell knows, why this is needed
			double ang1=ang*3.1415/180.0;
			double x=r*cos(ang1);
			double y=r*sin(ang1);
			result.vertices.push_back(Vector2d(x,y));
		}
	} else if(polygon && PyList_Check(polygon)) {
		unsigned int n=PyList_Size(polygon);
		for(unsigned int i=0;i < n;i++) {
			PyObject *pypt = PyList_GetItem(polygon, i);
			if(PyList_Check(pypt) && PyList_Size(pypt) == 2) {
				double x=PyFloat_AsDouble(PyList_GetItem(pypt, 0));
				double y=PyFloat_AsDouble(PyList_GetItem(pypt, 1));
				result.vertices.push_back(Vector2d(x,y));
			}
		}
	}
	if(result.vertices.size() < 3)
	{
		Outline2d err;
		err.vertices.push_back(Vector2d(0,0));
		err.vertices.push_back(Vector2d(10,0));
		err.vertices.push_back(Vector2d(10,10));
		return err;
	}
	return result;
}

double python_doublefunc(void *v_cbfunc, double arg)
{
	PyObject *cbfunc = (PyObject *) v_cbfunc;
	double result=0;
	PyObject* args = PyTuple_Pack(1,PyFloat_FromDouble(arg));
	PyObject* funcresult = PyObject_CallObject(cbfunc, args);
	Py_XDECREF(args);
	if(funcresult)
		result=PyFloat_AsDouble(funcresult);
	return result;
}

/*
 * Try to call a python function by name using OpenSCAD module childs and OpenSCAD function arguments: argument order is childs, arguments
 */
PyObject *python_fromopenscad(const Value &val)
{	
    switch(val.type())
    {
      case Value::Type::UNDEFINED:
	return Py_None;
      case Value::Type::BOOL:
	return val.toBool()?Py_True:Py_False;
      case Value::Type::NUMBER:
	return  PyFloat_FromDouble(val.toDouble());
      case Value::Type::STRING:
	return PyUnicode_FromString(val.toString().c_str());
      case Value::Type::VECTOR:
	{
	  const VectorType& vec = val.toVector();
  	  PyObject *result=PyList_New(vec.size());
	  for(int j=0;j<vec.size();j++)
		PyList_SetItem(result,j,python_fromopenscad(vec[j]));
	  return result;
	}
//TODO  more types RANGE, OBJECT, FUNCTION
      default:
	return Py_None;
    }
    return Py_None;
}
void python_catch_error(std::string &errorstr)
{
    PyObject *pyExcType;
    PyObject *pyExcValue;
    PyObject *pyExcTraceback;
    PyErr_Fetch(&pyExcType, &pyExcValue, &pyExcTraceback);
    PyErr_NormalizeException(&pyExcType, &pyExcValue, &pyExcTraceback);
    if(pyExcType != nullptr) Py_XDECREF(pyExcType);

    if(pyExcValue != nullptr){
      PyObjectUniquePtr str_exc_value( PyObject_Repr(pyExcValue), PyObjectDeleter);
      PyObjectUniquePtr pyExcValueStr( PyUnicode_AsEncodedString(str_exc_value.get(), "utf-8", "~"), PyObjectDeleter);
      char *suberror = PyBytes_AS_STRING(pyExcValueStr.get());
      if(suberror != nullptr) errorstr +=  suberror;
      Py_XDECREF(pyExcValue);
    }
    if(pyExcTraceback != nullptr) {
      auto *tb_o = (PyTracebackObject *)pyExcTraceback;
      int line_num = tb_o->tb_lineno;
      errorstr += " in line ";
      errorstr += std::to_string(line_num);
      Py_XDECREF(pyExcTraceback);
    }
}

PyObject *python_callfunction(const std::shared_ptr<const Context> &cxt , const std::string &name, const std::vector<std::shared_ptr<Assignment> > &op_args, std::string &errorstr)
{
  PyObject *pFunc = nullptr;
  if(!pythonMainModule){
    return nullptr;
  }

  int dotpos=name.find(".");
  if(dotpos >= 0) {
    std::string varname = name.substr(1,dotpos-1); // assume its always in paranthesis
    std::string methodname = name.substr(dotpos+1,name.size()-dotpos-2);
    Value var = cxt->lookup_variable(varname, Location::NONE).clone();
    if(var.type() == Value::Type::PYTHONCLASS) {
      const PythonClassType &python_class = var.toPythonClass();
      PyObject* methodobj = PyUnicode_FromString(methodname.c_str());
      pFunc = PyObject_GenericGetAttr((PyObject *) python_class.ptr,methodobj);
    }
  }
  if(!pFunc) {
    PyObject *maindict = PyModule_GetDict(pythonMainModule.get());

    // search the function in all modules
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(maindict, &pos, &key, &value)) {
      PyObject *module = PyObject_GetAttrString(pythonMainModule.get(), PyUnicode_AsUTF8(key));
      if(module != nullptr){
        PyObject *moduledict = PyModule_GetDict(module);
        Py_DECREF(module);
        if(moduledict != nullptr) {
          pFunc = PyDict_GetItemString(moduledict, name.c_str());
         if(pFunc != nullptr) break;
        } 
      }
    }  
  }
  if (!pFunc) {
    errorstr="Function not found";    	  
    return nullptr;
  }
  if (!PyCallable_Check(pFunc)) {
    errorstr="Function not callable";    	  
    return nullptr;
  }

  PyObject *args = PyTuple_New(op_args.size());
  for(unsigned int i=0;i<op_args.size();i++)
  {
    Assignment *op_arg=op_args[i].get();

    std::shared_ptr<Expression> expr=op_arg->getExpr();
    PyObject *value= python_fromopenscad( expr.get()->evaluate(cxt));
    if(value != nullptr) {
      PyTuple_SetItem(args, i, value);
    }

  }
  PyObject* funcresult = PyObject_CallObject(pFunc, args);
  Py_XDECREF(args);

  errorstr="";
  if(funcresult == nullptr) {
    python_catch_error(errorstr);	  
    PyErr_SetString(PyExc_TypeError, errorstr.c_str());

    return nullptr;
  }
  return funcresult;
}

/*
 * Actually trying use python to evaluate a OpenSCAD Module
 */

std::shared_ptr<AbstractNode> python_modulefunc(const ModuleInstantiation *op_module,const std::shared_ptr<const Context> &cxt, std::string &error) // null & error: error, else: None
{
  PyObject *dummydict;   
  std::shared_ptr<AbstractNode> result=nullptr;
  std::string errorstr = "";
  {
    PyObject *funcresult = python_callfunction(cxt,op_module->name(),op_module->arguments, errorstr);
    if (errorstr.size() > 0){
      error = errorstr;
      return nullptr;
    }
    if(funcresult == nullptr) {error="function not executed"; return nullptr; }

    if(funcresult->ob_type == &PyOpenSCADType) result=PyOpenSCADObjectToNode(funcresult, &dummydict);
    Py_XDECREF(funcresult);
    errorstr="";
  }
  return result;
}

/*
 * Converting a python result to an openscad result. extra function required as it might call itself hierarchically
 */

Value python_convertresult(PyObject *arg, int &error)
{
  error=0;	
  if(arg == nullptr) return Value::undefined.clone();
  if(PyList_Check(arg)) {
    VectorType vec(nullptr);
    for(int i=0;i<PyList_Size(arg);i++) {
      PyObject *item=PyList_GetItem(arg,i);
      int suberror;
      vec.emplace_back(python_convertresult(item,suberror));
      error |= suberror;
    }
    return std::move(vec);
  } else if(PyFloat_Check(arg)) { return { PyFloat_AsDouble(arg) }; 
  } else if(arg == Py_False) { return false;
  } else if(arg == Py_True) {  return true; 
  } else if(PyLong_Check(arg))  { return { (double) PyLong_AsLong(arg) }; }
  else if(PyUnicode_Check(arg)) {
    auto str = std::string(PyUnicode_AsUTF8(arg));
    return { str } ;
  } else if(arg == Py_None) { return Value::undefined.clone(); 
  } else if(arg->ob_type->tp_base == &PyBaseObject_Type) {
	  Py_INCREF(arg);
	  return PythonClassType(arg);
  } else {
    printf("unsupported type %s\n",arg->ob_type->tp_base->tp_name);	  
    PyErr_SetString(PyExc_TypeError, "Unsupported function result\n");
    error=1;
  }
  return Value::undefined.clone();
}

/*
 * Actually trying use python to evaluate a OpenSCAD Function
 */

Value python_functionfunc(const FunctionCall *call,const std::shared_ptr<const Context> &cxt, int &error  )
{
  std::string errorstr = "";
  PyObject *funcresult = python_callfunction(cxt,call->name, call->arguments, errorstr);
  if (errorstr.size() > 0)
  {
    PyErr_SetString(PyExc_TypeError, errorstr.c_str());
    return Value::undefined.clone();
  }
  if(funcresult == nullptr) return Value::undefined.clone();

  Value res = python_convertresult(funcresult, error);
  Py_XDECREF(funcresult);
  return res;
}

extern PyObject *PyInit_data(void);
PyMODINIT_FUNC PyInit_PyData(void);
/*
 * Main python evaluation entry
 */
#ifdef HAVE_PYTHON_YIELD
std::vector<PyObject *> python_orphan_objs;
extern "C" {
	void set_object_callback(void (*object_capture_callback)(PyObject *));
}
void openscad_object_callback(PyObject *obj) {
	if(obj->ob_type == &PyOpenSCADType) {
  		Py_INCREF(obj);
		python_orphan_objs.push_back(obj);
	}
}
#endif
void initPython(double time)
{
  if(pythonInitDict) { /* If already initialized, undo to reinitialize after */
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *maindict = PyModule_GetDict(pythonMainModule.get());
    while (PyDict_Next(maindict, &pos, &key, &value)) {
      PyObjectUniquePtr key_(PyUnicode_AsEncodedString(key, "utf-8", "~"), PyObjectDeleter);
      if(key_ == nullptr) continue;
      const char *key_str =  PyBytes_AS_STRING(key_.get());
      if(key_str == nullptr) continue;
      if (std::find(std::begin(pythonInventory), std::end(pythonInventory), key_str) == std::end(pythonInventory))
      {
        PyDict_DelItemString(maindict, key_str);
      }
      // bug in  PyDict_GetItemString, thus iterating
      if(strcmp(key_str,"sys") == 0) {
        PyObject *sysdict = PyModule_GetDict(value);
	if(sysdict == nullptr) continue;
	// get builtin_module_names
        PyObject *key1, *value1;
        Py_ssize_t pos1 = 0;
        while (PyDict_Next(sysdict, &pos1, &key1, &value1)) {
          PyObjectUniquePtr key1_(PyUnicode_AsEncodedString(key1, "utf-8", "~"), PyObjectDeleter);
          if(key1_ == nullptr) continue;
          const char *key1_str =  PyBytes_AS_STRING(key1_.get());
          if(strcmp(key1_str,"modules") == 0) {
            PyObject *key2, *value2;
            Py_ssize_t pos2 = 0;
            while (PyDict_Next(value1, &pos2, &key2, &value2)) {
              PyObjectUniquePtr key2_(PyUnicode_AsEncodedString(key2, "utf-8", "~"), PyObjectDeleter);
              if(key2_ == nullptr) continue;
              const char *key2_str =  PyBytes_AS_STRING(key2_.get());
	      if(key2_str == nullptr) continue;
	      if(!PyModule_Check(value2)) continue;

	      PyObject *modrepr = PyObject_Repr(value2);
	      PyObject* modreprobj = PyUnicode_AsEncodedString(modrepr, "utf-8", "~");
              const char *modreprstr = PyBytes_AS_STRING(modreprobj);
	      if(modreprstr == nullptr) continue;
	      if(strstr(modreprstr,"(frozen)") != nullptr) continue;
	      if(strstr(modreprstr,"(built-in)") != nullptr) continue;
	      if(strstr(modreprstr,"/encodings/") != nullptr) continue;
	      if(strstr(modreprstr,"_frozen_") != nullptr) continue;
              PyDict_DelItem(value1, key2);

	    }
          }
        }
      }
    }
  } else {
    PyPreConfig preconfig;
    PyPreConfig_InitPythonConfig(&preconfig);
    Py_PreInitialize(&preconfig);
//    PyEval_InitThreads(); // https://stackoverflow.com/questions/47167251/pygilstate-ensure-causing-deadlock

#ifdef HAVE_PYTHON_YIELD
    set_object_callback(openscad_object_callback);
#endif
    PyImport_AppendInittab("openscad", &PyInit_openscad);
    PyImport_AppendInittab("libfive", &PyInit_data);
    PyConfig config;
    PyConfig_InitPythonConfig(&config);
    std::string libdir;
    std::ostringstream stream;
#ifdef _WIN32
    char sepchar = ';';
    stream << PlatformUtils::applicationPath() << "\\..\\libraries\\python";
#else
    char sepchar = ':';
    stream << PlatformUtils::applicationPath() << "/../libraries/python";
    stream << sepchar + PlatformUtils::applicationPath() << "/../lib/python"  <<  PY_MAJOR_VERSION  <<  "."  <<  PY_MINOR_VERSION ; // find it where linuxdeply put it
#endif   
    stream << sepchar << PlatformUtils::userLibraryPath() << sepchar << ".";
    PyConfig_SetBytesString(&config, &config.pythonpath_env, stream.str().c_str());
    PyStatus status = Py_InitializeFromConfig(&config);
    if (PyStatus_Exception(status)) {
      LOG( message_group::Error, "Python not found. Is it installed ?");
      return;
    }
    PyConfig_Clear(&config);

    pythonMainModule.reset(PyImport_AddModule("__main__"));
    pythonMainModuleInitialized = pythonMainModule != nullptr;
    pythonInitDict.reset(PyModule_GetDict(pythonMainModule.get()));
    pythonRuntimeInitialized = pythonInitDict != nullptr;
    PyInit_PyOpenSCAD();
    PyInit_PyData();
    PyRun_String("from builtins import *\n", Py_file_input, pythonInitDict.get(), pythonInitDict.get());
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *maindict = PyModule_GetDict(pythonMainModule.get());
    while (PyDict_Next(maindict, &pos, &key, &value)) {
      PyObjectUniquePtr key1(PyUnicode_AsEncodedString(key, "utf-8", "~"), PyObjectDeleter);
      const char *key_str =  PyBytes_AsString(key1.get());
      if(key_str != NULL) pythonInventory.push_back(key_str);
    }

  }
  std::ostringstream stream;
  stream << "fa=12.0\nfn=0.0\nfs=2.0\nt=" << time << "\nphi=" << 2*G_PI*time;
  PyRun_String(stream.str().c_str(), Py_file_input, pythonInitDict.get(), pythonInitDict.get());
  customizer_parameters_finished = customizer_parameters;
  customizer_parameters.clear();
}

void finishPython(void)
{
#ifdef HAVE_PYTHON_YIELD
      set_object_callback(NULL);
      if(python_result_node == nullptr) {
        if(python_orphan_objs.size() == 1) {
  	  python_result_node = PyOpenSCADObjectToNode(python_orphan_objs[0]);
        } else if(python_orphan_objs.size() > 1) {
          DECLARE_INSTANCE
	  auto node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);
          int n = python_orphan_objs.size();
          for (int i = 0; i < n; i++) {
            std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNode(python_orphan_objs[i]);
            node->children.push_back(child);
          } 
          python_result_node=node;
        } 
      }
#endif
}

std::string evaluatePython(const std::string & code, bool dry_run)
{
  std::string error;
  python_result_node = nullptr;
  python_result_handle.clear();
  PyObject *pyExcType = nullptr;
  PyObjectUniquePtr pyExcValue (nullptr, PyObjectDeleter);
  PyObjectUniquePtr pyExcTraceback (nullptr, PyObjectDeleter);
  /* special python code to catch errors from stdout and stderr and make them available in OpenSCAD console */
  for(ModuleInstantiation *mi : modinsts_list){
    delete mi; // best time to delete it
  }
  modinsts_list.clear();
  pythonDryRun=dry_run;
  if(!pythonMainModuleInitialized)
	  return "Python not initialized";
  const char *python_init_code="\
import sys\n\
class InputCatcher:\n\
   def __init__(self):\n\
      self.data = \"modules\"\n\
   def read(self):\n\
      return self.data\n\
   def readline(self):\n\
      return self.data\n\
class OutputCatcher:\n\
   def __init__(self):\n\
      self.data = ''\n\
   def write(self, stuff):\n\
      self.data = self.data + stuff\n\
   def flush(self):\n\
      pass\n\
catcher_in = InputCatcher()\n\
catcher_out = OutputCatcher()\n\
catcher_err = OutputCatcher()\n\
stdin_bak=sys.stdin\n\
stdout_bak=sys.stdout\n\
stderr_bak=sys.stderr\n\
sys.stdin = catcher_in\n\
sys.stdout = catcher_out\n\
sys.stderr = catcher_err\n\
";
  const char *python_exit_code="\
sys.stdin = stdin_bak\n\
sys.stdout = stdout_bak\n\
sys.stderr = stderr_bak\n\
";

    PyRun_SimpleString(python_init_code);
#ifdef HAVE_PYTHON_YIELD
    for(auto obj : python_orphan_objs) {
        Py_DECREF(obj);
    }
    python_orphan_objs.clear();
#endif
    PyObjectUniquePtr result(nullptr, PyObjectDeleter);
    result.reset(PyRun_String(code.c_str(), Py_file_input, pythonInitDict.get(), pythonInitDict.get())); /* actual code is run here */


    if(result  == nullptr) {
      PyErr_Print();
      error = ""; 
      python_catch_error(error);
    } 
    for(int i=0;i<2;i++)
    {
      PyObjectUniquePtr catcher(nullptr, PyObjectDeleter);
      catcher.reset( PyObject_GetAttrString(pythonMainModule.get(), i==1?"catcher_err":"catcher_out"));
      if(catcher == nullptr) continue;
      PyObjectUniquePtr command_output(nullptr, PyObjectDeleter);
      command_output.reset(PyObject_GetAttrString(catcher.get(), "data"));

      PyObjectUniquePtr command_output_value(nullptr,  PyObjectDeleter);
      command_output_value.reset(PyUnicode_AsEncodedString(command_output.get(), "utf-8", "~"));
      const char *command_output_bytes =  PyBytes_AS_STRING(command_output_value.get());
      if(command_output_bytes != nullptr && *command_output_bytes != '\0')
      {
        if(i ==1) error += command_output_bytes; /* output to console */
        else LOG(command_output_bytes); /* error to LOG */
      }
    }
    PyRun_SimpleString(python_exit_code);
    return error;
}
/*
 * the magical Python Type descriptor for an OpenSCAD Object. Adding more fields makes the type more powerful
 */


int python__setitem__(PyObject *dict, PyObject *key, PyObject *v);
PyObject *python__getitem__(PyObject *dict, PyObject *key);

PyObject *python__getattro__(PyObject *dict, PyObject *key)
{
	PyObject *result=python__getitem__(dict,key);
	if(result == Py_None || result == nullptr)  result = PyObject_GenericGetAttr(dict,key);
	return result;
}

int python__setattro__(PyObject *dict, PyObject *key, PyObject *v)
{
	return python__setitem__(dict, key, v);
}


PyTypeObject PyOpenSCADType = {
    PyVarObject_HEAD_INIT(nullptr, 0)
    "PyOpenSCAD",             			/* tp_name */
    sizeof(PyOpenSCADObject), 			/* tp_basicsize */
    0,                         			/* tp_itemsize */
    (destructor) PyOpenSCADObject_dealloc,	/* tp_dealloc */
    0,                         			/* vectorcall_offset */
    0,                         			/* tp_getattr */
    0,                         			/* tp_setattr */
    0,                         			/* tp_as_async */
    python_str,               			/* tp_repr */
    &PyOpenSCADNumbers,        			/* tp_as_number */
    0,                         			/* tp_as_sequence */
    &PyOpenSCADMapping,        			/* tp_as_mapping */
    0,                         			/* tp_hash  */
    0,                         			/* tp_call */
    python_str,                			/* tp_str */
    python__getattro__,      			/* tp_getattro */
    python__setattro__,  			/* tp_setattro */
    0,                         			/* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/* tp_flags */
    "PyOpenSCAD Object",          		/* tp_doc */
    0,                         			/* tp_traverse */
    0,                         			/* tp_clear */
    0,                         			/* tp_richcompare */
    0,                         			/* tp_weaklistoffset */
    0,                         			/* tp_iter */
    0,                         			/* tp_iternext */
    PyOpenSCADMethods,             		/* tp_methods */
    0,             				/* tp_members */
    0,                         			/* tp_getset */
    0,                         			/* tp_base */
    0,                         			/* tp_dict */
    0,                         			/* tp_descr_get */
    0,                         			/* tp_descr_set */
    0,                         			/* tp_dictoffset */
    (initproc) PyOpenSCADInit,      		/* tp_init */
    PyOpenSCADObject_alloc,    			/* tp_alloc */
    PyOpenSCADObject_new,                	/* tp_new */
};



static PyModuleDef OpenSCADModule = {
  PyModuleDef_HEAD_INIT,
  "openscad",
  "OpenSCAD Python Module",
  -1,
  PyOpenSCADFunctions,
  NULL, NULL, NULL, NULL
};

static PyObject *PyInit_openscad(void)
{
  return PyModule_Create(&OpenSCADModule);
}

PyMODINIT_FUNC PyInit_PyOpenSCAD(void)
{
  PyObject *m;

  if (PyType_Ready(&PyOpenSCADType) < 0) return NULL;
  m = PyModule_Create(&OpenSCADModule);
  if (m == NULL) return NULL;

  Py_INCREF(&PyOpenSCADType);
  PyModule_AddObject(m, "openscad", (PyObject *)&PyOpenSCADType);
  return m;
}

// ----------------------------------------------
// IPython Interpreter side
// ----------------------------------------------

static PyStatus
pymain_init(void)
{
    PyStatus status;

//    if (_PyStatus_EXCEPTION(status)) {
//        return status;
//    }

    PyPreConfig preconfig;
    PyPreConfig_InitPythonConfig(&preconfig);
//    status = _Py_PreInitializeFromPyArgv(&preconfig, args);
//    if (_PyStatus_EXCEPTION(status)) {
//        return status;
//    }

    PyConfig config;
    PyConfig_InitPythonConfig(&config);

//    if (args->use_bytes_argv) {
//        status = PyConfig_SetBytesArgv(&config, args->argc, args->bytes_argv);
//    }
//    else {
//        status = PyConfig_SetArgv(&config, args->argc, args->wchar_argv);
//    }
//    if (_PyStatus_EXCEPTION(status)) {
//        goto done;
//    }

    status = Py_InitializeFromConfig(&config);
//    if (_PyStatus_EXCEPTION(status)) {
//        goto done;
//    }
//    status = 0; // PyStatus_Ok;

done:
    PyConfig_Clear(&config);
    return status;
}


/* Write an exitcode into *exitcode and return 1 if we have to exit Python.
   Return 0 otherwise. */
static int
pymain_run_interactive_hook(int *exitcode)
{
    PyObject *sys, *hook, *result;
    sys = PyImport_ImportModule("sys");
    if (sys == NULL) {
        goto error;
    }

    hook = PyObject_GetAttrString(sys, "__interactivehook__");
    Py_DECREF(sys);
    if (hook == NULL) {
        PyErr_Clear();
        return 0;
    }

    if (PySys_Audit("cpython.run_interactivehook", "O", hook) < 0) {
        goto error;
    }
    result  =  PyObject_CallNoArgs(hook);
    Py_DECREF(hook);
    if (result == NULL) {
        goto error;
    }
    Py_DECREF(result);
    return 0;

error:
    PySys_WriteStderr("Failed calling sys.__interactivehook__\n");
//    return pymain_err_print(exitcode);
    return 0;
}






static void
pymain_repl(PyConfig *config, int *exitcode)
{
//    if (!config->inspect && _Py_GetEnv(config->use_environment, "PYTHONINSPECT")) {
//        pymain_set_inspect(config, 1);
//    }

//    if (!(config->inspect && stdin_is_interactive(config) && config_run_code(config))) {
//        return;
//    }

//    pymain_set_inspect(config, 0);
    if (pymain_run_interactive_hook(exitcode)) {
        return;
    }
    PyCompilerFlags cf = _PyCompilerFlags_INIT;

    int res = PyRun_AnyFileFlags(stdin, "<stdin>", &cf);
//    *exitcode = (res != 0);
}


static void
pymain_run_python(int *exitcode)
{
    PyObject *main_importer_path = NULL;
    PyInterpreterState *interp = PyInterpreterState_Get();
    /* pymain_run_stdin() modify the config */
    PyConfig *config = (PyConfig*)_PyInterpreterState_GetConfig(interp);

//    if (_PyStatus_EXCEPTION(_PyPathConfig_UpdateGlobal(config))) {
//        goto error;
//    }

//    if (config->run_filename != NULL) {
//        if (pymain_get_importer(config->run_filename, &main_importer_path,
//                                exitcode)) {
//            return;
//        }
//    }
    // import readline and rlcompleter before script dir is added to sys.path
//    pymain_import_readline(config);

//    PyObject *path0 = NULL;
//    if (main_importer_path != NULL) {
//        path0 = Py_NewRef(main_importer_path);
//    }
//    else if (!config->safe_path) {
//        int res = _PyPathConfig_ComputeSysPath0(&config->argv, &path0);
//        if (res < 0) {
//            goto error;
//        }
//        else if (res == 0) {
//            Py_CLEAR(path0);
//        }
//    }
//    if (path0 != NULL) {
//        wchar_t *wstr = PyUnicode_AsWideCharString(path0, NULL);
//        if (wstr == NULL) {
//            Py_DECREF(path0);
//            goto error;
//        }
//        config->sys_path_0 = _PyMem_RawWcsdup(wstr);
//        PyMem_Free(wstr);
//        if (config->sys_path_0 == NULL) {
//            Py_DECREF(path0);
//            goto error;
//        }
//        int res = pymain_sys_path_add_path0(interp, path0);
//        Py_DECREF(path0);
//        if (res < 0) {
//            goto error;
//        }
//    }
//
//    pymain_header(config);
//
//    _PyInterpreterState_SetRunningMain(interp);
//    assert(!PyErr_Occurred());
//
//    if (config->run_command) {
//        *exitcode = pymain_run_command(config->run_command);
//    }
//    else if (config->run_module) {
//        *exitcode = pymain_run_module(config->run_module, 1);
//    }
//    else if (main_importer_path != NULL) {
//        *exitcode = pymain_run_module(L"__main__", 0);
//    }
//    else if (config->run_filename != NULL) {
//        *exitcode = pymain_run_file(config);
//    }
//    else {
//        *exitcode = pymain_run_stdin(config);
//    }

    pymain_repl(config, exitcode);
    goto done;

error:
//    *exitcode = pymain_exit_err_print();

done:
//    _PyInterpreterState_SetNotRunningMain(interp);
    Py_XDECREF(main_importer_path);
}

int
Py_RunMain(void)
{
    int exitcode = 0;

    pymain_run_python(&exitcode);

//    if (Py_FinalizeEx() < 0) {
//        exitcode = 120;
//    }

//    pymain_free();

//    if (_PyRuntime.signals.unhandled_keyboard_interrupt) {
//        exitcode = exit_sigint();
//    }

    return exitcode;
}


void ipython(void) {
/*	
    _PyArgv args = {
        .argc = argc,
        .use_bytes_argv = 1,
        .bytes_argv = argv,
        .wchar_argv = NULL};
*/
//    PyStatus status = pymain_init();
    initPython(0.0);
/*    
    if (_PyStatus_IS_EXIT(status)) {
        pymain_free();
        return status.exitcode;
    }
    if (_PyStatus_EXCEPTION(status)) {
        pymain_exit_error(status);
    }

*/
    Py_RunMain();
}

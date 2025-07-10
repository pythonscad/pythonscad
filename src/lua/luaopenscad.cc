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

/* MuJS is Copyright © 2013-2020 Artifex Software, Inc */

#include "luaopenscad.h"
#include "CsgOpNode.h"
#include "Value.h"
#include "Expression.h"
#include "PlatformUtils.h"
#include <Context.h>
#include <Selection.h>


std::shared_ptr<AbstractNode> lua_result_node = nullptr; /* global result veriable containing the perl created result */

#if 0
// #define HAVE_PYTHON_YIELD
static PyObject *PyInit_openscad(void);


PyObject *pythonInitDict = nullptr;
PyObject *pythonMainModule = nullptr ;
std::list<std::string> pythonInventory;
AssignmentList customizer_parameters;
AssignmentList customizer_parameters_finished;
std::vector<SelectedObject> python_result_handle;
bool python_active;  /* if python is actually used during evaluation */
bool python_trusted; /* global Python trust flag */
#include "PlatformUtils.h"
bool pythonMainModuleInitialized = false;
bool pythonRuntimeInitialized = false;

std::vector<std::string> mapping_name;
std::vector<std::string> mapping_code;
std::vector<int> mapping_level;


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

std::string python_hierdump(const std::shared_ptr<AbstractNode> &node)
{
  std::string dump = node->toString();	
  auto children = node->getChildren();
  if(children.size() < 1) dump += ";";
  else if(children.size() < 2) dump += python_hierdump(children[0]);
  else {
   dump += "{ ";	  
   for(int i=0;i<children.size();i++) dump += python_hierdump(children[i]);
   dump += " }";	  
  }

  return dump;
}
void python_build_hashmap(const std::shared_ptr<AbstractNode> &node)
{
  PyObject *maindict = PyModule_GetDict(pythonMainModule);
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  std::string code=python_hierdump(node);
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
  for(const auto &child:node->getChildren()) {
    python_build_hashmap(child);	  
  }
}

void python_retrieve_pyname(const std::shared_ptr<AbstractNode> &node)
{
  std::string name;	
  int level=-1;
  std::string my_code = python_hierdump(node);
  for(int i=0;i<mapping_code.size();i++) {
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

/*
 * Tries to extract an 3D vector out of a python list
 */

int python_vectorval(PyObject *vec, double *x, double *y, double *z, double *w)
{
  *x = 1.0;
  *y = 1.0;
  if(w != NULL ) *w = 0;
  if (PyList_Check(vec)) {
    if (PyList_Size(vec) >= 1) {
      if (python_numberval(PyList_GetItem(vec, 0), x)) return 1;
    }
    if (PyList_Size(vec) >= 2) {
      if (python_numberval(PyList_GetItem(vec, 1), y)) return 1;
    }
    if (PyList_Check(vec) && PyList_Size(vec) >= 3) {
      if (python_numberval(PyList_GetItem(vec, 2), z)) return 1;
    }
    if (PyList_Check(vec) && PyList_Size(vec) >= 4 && w != NULL) {
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

/*
 * Helper function to extract actual values for fn, fa and fs
 */

void get_fnas(double& fn, double& fa, double& fs) {
  PyObject *mainModule = PyImport_AddModule("__main__");
  if (mainModule == nullptr) return;
  PyObject *varFn = PyObject_GetAttrString(mainModule, "fn");
  PyObject *varFa = PyObject_GetAttrString(mainModule, "fa");
  PyObject *varFs = PyObject_GetAttrString(mainModule, "fs");
  if (varFn != nullptr){
    fn = PyFloat_AsDouble(varFn);
    Py_XDECREF(varFn);
  }
  if (varFa != nullptr){
    fa = PyFloat_AsDouble(varFa);
    Py_XDECREF(varFa);
  }
  if (varFs != nullptr){
    fs = PyFloat_AsDouble(varFs);
    Py_XDECREF(varFs);
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
		for(unsigned int i=0;i < fn;i++) {
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

PyObject *python_callfunction(const std::shared_ptr<const Context> &cxt , const std::string &name, const std::vector<std::shared_ptr<Assignment> > &op_args, const char *&errorstr)
{
  PyObject *pFunc = nullptr;
  if(!pythonMainModule){
    return nullptr;
  }
  PyObject *maindict = PyModule_GetDict(pythonMainModule);

  // search the function in all modules
  PyObject *key, *value;
  Py_ssize_t pos = 0;

  while (PyDict_Next(maindict, &pos, &key, &value)) {
    PyObject *module = PyObject_GetAttrString(pythonMainModule, PyUnicode_AsUTF8(key));
    if(module != nullptr){
      PyObject *moduledict = PyModule_GetDict(module);
      Py_DECREF(module);
      if(moduledict != nullptr) {
        pFunc = PyDict_GetItemString(moduledict, name.c_str());
        if(pFunc != nullptr) break;
      }
    }
  }
  if (!pFunc) {
    return nullptr;
  }
  if (!PyCallable_Check(pFunc)) {
    return nullptr;
  }

  PyObject *args = PyTuple_New(op_args.size());
  for(int i=0;i<op_args.size();i++)
  {
    Assignment *op_arg=op_args[i].get();

    std::shared_ptr<Expression> expr=op_arg->getExpr();
    Value val = expr.get()->evaluate(cxt);
    PyObject *value=NULL;
    switch(val.type())
    {
      case Value::Type::NUMBER:
	value =  PyFloat_FromDouble(val.toDouble());
        break;
      case Value::Type::STRING:
	value = PyUnicode_FromString(val.toString().c_str());
        break;
//TODO  more types RANGE, VECTOR, OBEJCT, FUNCTION
      default:
	value= PyLong_FromLong(-1);
        break;
    }
    if(value != nullptr) {
      PyTuple_SetItem(args, i, value);
      Py_XDECREF(value);
    }

  }
  PyObject* funcresult = PyObject_CallObject(pFunc, args);
  Py_XDECREF(args);

  if(funcresult == nullptr) {
    PyObject *pyExcType;
    PyObject *pyExcValue;
    PyObject *pyExcTraceback;
    PyErr_Fetch(&pyExcType, &pyExcValue, &pyExcTraceback);
    PyErr_NormalizeException(&pyExcType, &pyExcValue, &pyExcTraceback);
    if(pyExcType != nullptr) Py_XDECREF(pyExcType);
    if(pyExcTraceback != nullptr) Py_XDECREF(pyExcTraceback);

    if(pyExcValue != nullptr){
      PyObject* str_exc_value = PyObject_Repr(pyExcValue);
      PyObject* pyExcValueStr = PyUnicode_AsEncodedString(str_exc_value, "utf-8", "~");
      Py_XDECREF(str_exc_value);
      errorstr =  PyBytes_AS_STRING(pyExcValueStr);
      PyErr_SetString(PyExc_TypeError, errorstr);
      Py_XDECREF(pyExcValueStr);
      Py_XDECREF(pyExcValue);
    }

    return nullptr;
  }
  return funcresult;
}

/*
 * Actually trying use python to evaluate a OpenSCAD Module
 */

std::shared_ptr<AbstractNode> python_modulefunc(const ModuleInstantiation *op_module,const std::shared_ptr<const Context> &cxt, int *modulefound)
{
   *modulefound=0;
  PyObject *dummydict;   
  std::shared_ptr<AbstractNode> result=nullptr;
  const char *errorstr = nullptr;
  {
    PyObject *funcresult = python_callfunction(cxt,op_module->name(),op_module->arguments, errorstr);
    if (errorstr != nullptr){
      PyErr_SetString(PyExc_TypeError, errorstr);
      return nullptr;
    }
    *modulefound=1;
    if(funcresult == nullptr) return nullptr;

    if(funcresult->ob_type == &PyOpenSCADType) result=PyOpenSCADObjectToNode(funcresult, &dummydict);
    Py_XDECREF(funcresult);
  }
  return result;
}

/*
 * Converting a python result to an openscad result. extra function required as it might call itself hierarchically
 */

Value python_convertresult(PyObject *arg)
{
  if(arg == nullptr) return Value::undefined.clone();
  if(PyList_Check(arg)) {
    VectorType vec(nullptr);
    for(int i=0;i<PyList_Size(arg);i++) {
      PyObject *item=PyList_GetItem(arg,i);
      vec.emplace_back(python_convertresult(item));
    }
    return std::move(vec);
  } else if(PyFloat_Check(arg)) { return { PyFloat_AsDouble(arg) }; }
  else if(PyUnicode_Check(arg)) {
    PyObject* repr = PyObject_Repr(arg);
    PyObject* strobj = PyUnicode_AsEncodedString(repr, "utf-8", "~");
    Py_XDECREF(repr);
    const char *chars =  PyBytes_AS_STRING(strobj);
    auto str = std::string(chars);
    Py_XDECREF(strobj);
    return { str } ;
  } else {
    PyErr_SetString(PyExc_TypeError, "Unsupported function result\n");
    return Value::undefined.clone();
  }
}

/*
 * Actually trying use python to evaluate a OpenSCAD Function
 */

Value python_functionfunc(const FunctionCall *call,const std::shared_ptr<const Context> &cxt  )
{
  const char *errorstr = nullptr;
  PyObject *funcresult = python_callfunction(cxt,call->name, call->arguments, errorstr);
  if (errorstr != nullptr)
  {
    PyErr_SetString(PyExc_TypeError, errorstr);
    return Value::undefined.clone();
  }
  if(funcresult == nullptr) return Value::undefined.clone();

  Value res = python_convertresult(funcresult);
  Py_XDECREF(funcresult);
  return res;
}

#ifdef ENABLE_LIBFIVE
extern PyObject *PyInit_libfive(void);
PyMODINIT_FUNC PyInit_PyLibFive(void);
#endif
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
    PyObject *maindict = PyModule_GetDict(pythonMainModule);
    while (PyDict_Next(maindict, &pos, &key, &value)) {
      PyObject* key1 = PyUnicode_AsEncodedString(key, "utf-8", "~");
      if(key1 != nullptr) {
        const char *key_str =  PyBytes_AS_STRING(key1);
        if(key_str == nullptr) continue;
        if (std::find(std::begin(pythonInventory), std::end(pythonInventory), key_str) == std::end(pythonInventory))
        {
          PyDict_DelItemString(maindict, key_str);
        }
        Py_XDECREF(key1);
      }
    }
  } else {
    PyPreConfig preconfig;
    PyPreConfig_InitPythonConfig(&preconfig);
    Py_PreInitialize(&preconfig);

#ifdef HAVE_PYTHON_YIELD
    set_object_callback(openscad_object_callback);
#endif
    PyImport_AppendInittab("openscad", &PyInit_openscad);
#ifdef ENABLE_LIBFIVE	    
    PyImport_AppendInittab("libfive", &PyInit_libfive);
#endif	    
    PyConfig config;
    PyConfig_InitPythonConfig(&config);
    char libdir[256];
#ifdef _WIN32
    snprintf(libdir, 256, "%s\\..\\libraries\\python\\:.",PlatformUtils::applicationPath().c_str()); /* add libraries/python to python search path */
#else
    snprintf(libdir, 256, "%s/../libraries/python/:.",PlatformUtils::applicationPath().c_str()); /* add libraries/python to python search path */
#endif   
    PyConfig_SetBytesString(&config, &config.pythonpath_env, libdir);
    PyStatus status = Py_InitializeFromConfig(&config);
    if (PyStatus_Exception(status)) {
      LOG( message_group::Error, "Python not found. Is it installed ?");
      return;
    }
    PyConfig_Clear(&config);

    pythonMainModule =  PyImport_AddModule("__main__");
    Py_XINCREF(pythonMainModule);
    pythonMainModuleInitialized = pythonMainModule != nullptr;
    pythonInitDict = PyModule_GetDict(pythonMainModule);
    Py_XINCREF(pythonInitDict);
    pythonRuntimeInitialized = pythonInitDict != nullptr;
    PyInit_PyOpenSCAD();
#ifdef ENABLE_LIBFIVE	    
    PyInit_PyLibFive();
#endif	    
    PyRun_String("from builtins import *\n", Py_file_input, pythonInitDict, pythonInitDict);
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *maindict = PyModule_GetDict(pythonMainModule);
    while (PyDict_Next(maindict, &pos, &key, &value)) {
      PyObject* key1 = PyUnicode_AsEncodedString(key, "utf-8", "~");
      const char *key_str =  PyBytes_AsString(key1);
      if(key_str != NULL) pythonInventory.push_back(key_str);
      Py_XDECREF(key1);
    }
  }
  char run_str[250];
  sprintf(run_str,"fa=12.0\nfn=0.0\nfs=2.0\nt=%g\nphi=%g",time,2*G_PI*time);
  PyRun_String(run_str, Py_file_input, pythonInitDict, pythonInitDict);
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

std::string evaluatePython(const std::string & code)
{
  std::string error;
  python_result_node = nullptr;
  python_result_handle.clear();
  PyObject *pyExcType = nullptr;
  PyObject *pyExcValue = nullptr;
  PyObject *pyExcTraceback = nullptr;
  /* special python code to catch errors from stdout and stderr and make them available in OpenSCAD console */
  if(!pythonMainModuleInitialized)
	  return "Python not initialized";
  const char *python_init_code="\
import sys\n\
class OutputCatcher:\n\
   def __init__(self):\n\
      self.data = ''\n\
   def write(self, stuff):\n\
      self.data = self.data + stuff\n\
   def flush(self):\n\
      pass\n\
catcher_out = OutputCatcher()\n\
catcher_err = OutputCatcher()\n\
stdout_bak=sys.stdout\n\
stderr_bak=sys.stderr\n\
sys.stdout = catcher_out\n\
sys.stderr = catcher_err\n\
";
  const char *python_exit_code="\
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
    PyObject *result = PyRun_String(code.c_str(), Py_file_input, pythonInitDict, pythonInitDict); /* actual code is run here */

    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *pFunc;

    if(result  == nullptr) PyErr_Print();
    PyRun_SimpleString(python_exit_code);

    for(int i=0;i<2;i++)
    {
      PyObject* catcher = PyObject_GetAttrString(pythonMainModule, i==1?"catcher_err":"catcher_out");
	if(catcher == nullptr) continue;
      PyObject* command_output = PyObject_GetAttrString(catcher, "data");
      Py_XDECREF(catcher);
      PyObject* command_output_value = PyUnicode_AsEncodedString(command_output, "utf-8", "~");
      Py_XDECREF(command_output);
      const char *command_output_bytes =  PyBytes_AS_STRING(command_output_value);
      if(command_output_bytes != nullptr && *command_output_bytes != '\0')
      {
        if(i ==1) error += command_output_bytes; /* output to console */
        else LOG(command_output_bytes); /* error to LOG */
      }
      Py_XDECREF(command_output_value);
    }

    PyErr_Fetch(&pyExcType, &pyExcValue, &pyExcTraceback); /* extract actual python stack trace in case of an expception and return the error string to the caller */
    if(pyExcType != nullptr) Py_XDECREF(pyExcType);
    if(pyExcValue != nullptr) {
      PyObject* str_exc_value = PyObject_Repr(pyExcValue);
      if(pyExcValue != nullptr) Py_XDECREF(pyExcValue);
      PyObject* pyExcValueStr = PyUnicode_AsEncodedString(str_exc_value, "utf-8", "~");
      const char *strExcValue =  PyBytes_AS_STRING(pyExcValueStr);
      if(strExcValue != nullptr && strcmp(strExcValue,"<NULL>") != 0) error += strExcValue;
      if(str_exc_value != nullptr) Py_XDECREF(str_exc_value);
    }
    if(pyExcTraceback != nullptr) {
      auto *tb_o = (PyTracebackObject *)pyExcTraceback;
      int line_num = tb_o->tb_lineno;
      error += " in line ";
      error += std::to_string(line_num);
      Py_XDECREF(pyExcTraceback);
    }

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

#endif

#include <stdio.h>


void LuaOpenSCADObjectFromNode(const std::shared_ptr<AbstractNode> &node)
{
}

std::shared_ptr<AbstractNode> LuaOpenSCADObjectToNode(void *data)
{
	return nullptr;	
}

//js_State *js_interp;

void initLua(double time)
{
	printf("initLua\n");
//  js_interp = js_newstate(NULL, NULL, JS_STRICT);
//  registerLuaFunctions();
}  
std::string evaluateLua(const std::string & code)
{
  char outbuf[BUFSIZ]="";
  printf("evaluateLuia %s\n",code.c_str());
#if 0	
  char errbuf[BUFSIZ]="";
  setbuf(stdout, outbuf);
  setbuf(stderr, errbuf);
  js_dostring(js_interp, code.c_str());
  fflush(stdout);
  fflush(stderr);
  setbuf(stdout, NULL);
  setbuf(stderr, NULL);
  if(*errbuf != '\0') LOG( message_group::Error, std::string(errbuf));
#endif  
  return outbuf;
}
void finishLua(void)
{
	printf("finishLua\n");
//  js_freestate(js_interp);
}

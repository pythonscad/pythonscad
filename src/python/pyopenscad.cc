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
#include <dlfcn.h>
#include <filesystem>

#include "pyopenscad.h"
#include "pydata.h"
#include "core/CsgOpNode.h"
#include "Value.h"
#include "Expression.h"
#include "PlatformUtils.h"
#include <Context.h>
#include <Selection.h>
#include "platform/PlatformUtils.h"
#include "primitives.h"
namespace fs = std::filesystem;

// #define HAVE_PYTHON_YIELD
extern "C" PyObject *PyInit_openscad(void);

bool python_active;
bool python_trusted;
fs::path python_scriptpath;
// https://docs.python.org/3.10/extending/newtypes.html

void PyObjectDeleter (PyObject *pObject) { Py_XDECREF_(pObject); };

PyObjectUniquePtr pythonInitDict(nullptr, PyObjectDeleter) ;
PyObjectUniquePtr pythonMainModule(nullptr, PyObjectDeleter) ;
std::list<std::string> pythonInventory;
AssignmentList customizer_parameters;
AssignmentList customizer_parameters_finished;
bool pythonDryRun=false;
std::shared_ptr<AbstractNode> python_result_node = nullptr; /* global result veriable containing the python created result */
PyObject *python_result_obj = nullptr;
std::vector<SelectedObject> python_result_handle;
std::vector<std::shared_ptr<AbstractNode> > shows;
bool python_runipython = false;
bool pythonMainModuleInitialized = false;
bool pythonRuntimeInitialized = false;

std::vector<std::string> mapping_name;
std::vector<std::string> mapping_code;
std::vector<int> mapping_level;
std::vector<std::shared_ptr<AbstractNode>> nodes_hold; // make sure, that these nodes are not yet freed
std::shared_ptr<AbstractNode> void_node, full_node;				       

void *libpython_handle=nullptr;

#define PYTHON_LOAD_FUNC(name, type) \
	pf.name= (type) dlsym(libpython_handle, #name);\
	printf("Function %s is %p\n", #name, pf.name);

struct PyKernel pf;

void initDynamic(void)
{
  if(libpython_handle) return;	
  libpython_handle = dlopen("/usr/lib/libpython3.13.so", RTLD_LAZY);	
  printf("handle=%p\n",libpython_handle);
  if(!libpython_handle) return;

  PYTHON_LOAD_FUNC(PyConfig_InitPythonConfig, void(*)(PyConfig *config) );
  PYTHON_LOAD_FUNC(PyPreConfig_InitPythonConfig, void(*)(PyPreConfig *config) );
  PYTHON_LOAD_FUNC(PyConfig_Clear, void(*)(PyConfig *config) );
  PYTHON_LOAD_FUNC(PyConfig_Read, PyStatus(*)(PyConfig *config) );
  PYTHON_LOAD_FUNC(PyConfig_SetBytesArgv, PyStatus(*)(PyConfig *config, Py_ssize_t argc, char * const *argv) );
  PYTHON_LOAD_FUNC(PyConfig_SetBytesString, PyStatus(*)( PyConfig *config, wchar_t **config_str, const char *str) );
  PYTHON_LOAD_FUNC(Py_PreInitialize, PyStatus(*)( const PyPreConfig *src_config) );
  PYTHON_LOAD_FUNC(Py_InitializeFromConfig, PyStatus(*)( const PyConfig *src_config) );

  PYTHON_LOAD_FUNC(PyType_GenericAlloc, PyObject *(*)(PyTypeObject *, Py_ssize_t) );
  PYTHON_LOAD_FUNC(_Py_Dealloc, void (*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyType_Ready, int (*)(PyTypeObject *) );

  PYTHON_LOAD_FUNC(PyEval_RestoreThread, void (*)(PyThreadState*) );
  PYTHON_LOAD_FUNC(PyEval_SaveThread, PyThreadState *(*)(void) );

  PYTHON_LOAD_FUNC(PyModule_GetDict, PyObject *(*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyBytes_AsString, const char *(*)(PyObject *) );

  PYTHON_LOAD_FUNC(PyImport_AddModule, PyObject *(*)(const char *) );
  PYTHON_LOAD_FUNC(PyImport_ImportModule, PyObject *(*)(const char *) );
  PYTHON_LOAD_FUNC(PyModule_Create2, PyObject *(*)(PyModuleDef *, int apiver) );

  PYTHON_LOAD_FUNC(PyArg_ParseTupleAndKeywords, int (*)(PyObject *, PyObject *, const char *, char **, ...) );
  PYTHON_LOAD_FUNC(PyErr_SetString, void (*)(PyObject *exception, const char *string) );

  PYTHON_LOAD_FUNC(PyRun_AnyFileExFlags, int (*)(FILE *fp, const char *filename,        int closeit, PyCompilerFlags *flags) );
  PYTHON_LOAD_FUNC(PyRun_SimpleStringFlags, int (*)(const char *filename, PyCompilerFlags *flags) );
  PYTHON_LOAD_FUNC(PyRun_StringFlags, PyObject *(*)(const char *, int, PyObject *, PyObject *, PyCompilerFlags *) );
  PYTHON_LOAD_FUNC(PyObject_CallObject, PyObject *(*)(PyObject *func, PyObject *args) );

  PYTHON_LOAD_FUNC(PyList_GetItem, PyObject *(*)(PyObject *, Py_ssize_t) );
  PYTHON_LOAD_FUNC(PyList_SetItem, void (*)(PyObject *, Py_ssize_t ,PyObject *) );
  PYTHON_LOAD_FUNC(PyList_New, PyObject *(*)(Py_ssize_t) );
  PYTHON_LOAD_FUNC(PyList_Size, Py_ssize_t (*)(PyObject *) );

  PYTHON_LOAD_FUNC(PyDict_New, PyObject *(*)(void) );
  PYTHON_LOAD_FUNC(PyDict_SetItem, PyObject *(*)(PyObject *, PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyDict_SetItemString, PyObject *(*)(PyObject *, const char *, PyObject *) );
  PYTHON_LOAD_FUNC(PyDict_GetItem, PyObject *(*)(PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyDict_GetItemString, PyObject *(*)(PyObject *, const char *) );
  PYTHON_LOAD_FUNC(PyDict_DelItem, void (*)(PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyDict_DelItemString, void (*)(PyObject *, const char *) );
  PYTHON_LOAD_FUNC(PyDict_Next, PyObject *(*)(PyObject *,Py_ssize_t *, PyObject **, PyObject ** ) );

  PYTHON_LOAD_FUNC(PyTuple_New, PyObject *(*)(Py_ssize_t size) );
  PYTHON_LOAD_FUNC(PyTuple_Pack, PyObject *(*)(Py_ssize_t, ...) );
  PYTHON_LOAD_FUNC(PyTuple_Size, Py_ssize_t (*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyTuple_GetItem, PyObject *(*)(PyObject *,Py_ssize_t ) );
  PYTHON_LOAD_FUNC(PyTuple_SetItem, void (*)(PyObject *, Py_ssize_t ,PyObject *) );

  PYTHON_LOAD_FUNC(PyFloat_FromDouble, PyObject *(*)(double) );
  PYTHON_LOAD_FUNC(PyFloat_AsDouble, double(*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyLong_AsLong, long(*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyLong_FromLong, PyObject *(*)(long) );

  PYTHON_LOAD_FUNC(PyUnicode_AsEncodedString, PyObject *(*)(PyObject *,const char *enc, const char *errors) );
  PYTHON_LOAD_FUNC(PyUnicode_FromString, PyObject *(*)(const char *u) );
  PYTHON_LOAD_FUNC(PyUnicode_FromStringAndSize, PyObject *(*)(const char *u, Py_ssize_t size) );
  PYTHON_LOAD_FUNC(PyUnicode_AsUTF8, const char *(*)(PyObject *) );

  PYTHON_LOAD_FUNC(PyErr_Fetch, void (*)(PyObject **, PyObject **, PyObject **) );
  PYTHON_LOAD_FUNC(PyType_IsSubtype, int (*)(PyTypeObject *, PyTypeObject *) );

  PYTHON_LOAD_FUNC(Py_ExitStatusException, void (*)(PyStatus err) );
  PYTHON_LOAD_FUNC(PyStatus_Exception, int (*)(PyStatus err) );
  PYTHON_LOAD_FUNC(PyStatus_IsExit, int (*)(PyStatus err) );
  PYTHON_LOAD_FUNC(PyWideStringList_Append, PyStatus (*)(PyWideStringList *list, const wchar_t *item) );

  PYTHON_LOAD_FUNC(PyModule_AddObject, int (*)(PyObject *mod, const char *, PyObject *value) );

  PYTHON_LOAD_FUNC(PyCallable_Check, int (*)(PyObject *) );
  PYTHON_LOAD_FUNC(PyErr_Clear, void (*)(void) );
  PYTHON_LOAD_FUNC(PyErr_NormalizeException, void (*)(PyObject**, PyObject**, PyObject**) );
  PYTHON_LOAD_FUNC(PyErr_Print, void (*)(void) );
  PYTHON_LOAD_FUNC(PyImport_AppendInittab, int (*)(const char *name, PyObject* (*initfunc)(void)) );
  PYTHON_LOAD_FUNC(PyObject_CallNoArgs, PyObject *(*)(PyObject *func) );
  PYTHON_LOAD_FUNC(PyObject_GenericGetAttr, PyObject *(*)(PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyObject_GetAttr, PyObject *(*)(PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyObject_GetAttrString, PyObject *(*)(PyObject *, const char *) );
  PYTHON_LOAD_FUNC(PyObject_HasAttr, int (*)(PyObject *, PyObject *) );
  PYTHON_LOAD_FUNC(PyObject_HasAttrString, int (*)(PyObject *, const char *) );
  PYTHON_LOAD_FUNC(PyObject_Repr, PyObject *(*)(PyObject *) );
  PYTHON_LOAD_FUNC(PySys_Audit, int (*)(const char *event, const char *format, ...) );
  PYTHON_LOAD_FUNC(PySys_WriteStdout, void (*)(const char *format, ...) );

  PYTHON_LOAD_FUNC(_Py_TrueStruct, PyObject *);
  PYTHON_LOAD_FUNC(_Py_FalseStruct, PyObject *);
  PYTHON_LOAD_FUNC(_Py_NoneStruct, PyObject *);
  PYTHON_LOAD_FUNC(PyExc_TypeError, PyObject *);
  PYTHON_LOAD_FUNC(PyFunction_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyList_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyFloat_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyLong_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyUnicode_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyTuple_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyModule_Type, PyTypeObject *);
  PYTHON_LOAD_FUNC(PyBaseObject_Type, PyTypeObject *);
}

void finishDynamic(void)
{
  if(!libpython_handle) return;	
  dlclose(libpython_handle);
  libpython_handle = nullptr;
}

void PyOpenSCADObject_dealloc(PyOpenSCADObject *self)
{
  Py_XDECREF_(self->dict);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyOpenSCADObject_alloc(PyTypeObject *cls, Py_ssize_t nitems)
{
  PyObject *self = pf.PyType_GenericAlloc(cls, nitems);
  ((PyOpenSCADObject *)self)->dict = pf.PyDict_New();
  PyObject *origin=pf.PyList_New(4);
  for(int i=0;i<4;i++) {
  	PyObject *row=pf.PyList_New(4);
	for(int j=0;j<4;j++)
		pf.PyList_SetItem(row,j,pf.PyFloat_FromDouble(i==j?1.0:0.0));
	pf.PyList_SetItem(origin,i,row);
  }
  pf.PyDict_SetItemString(((PyOpenSCADObject *)self)->dict,"origin",origin);
  Py_XDECREF_(origin);
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
    Py_XINCREF_((PyObject *) self);
    self->node = node;
    return (PyObject *)self;
  }
  return nullptr;
}

//PyGILState_STATE gstate=PyGILState_LOCKED;
PyThreadState *tstate=nullptr;

void python_lock(void){
//#ifndef _WIN32	
  if(tstate != nullptr && pythonInitDict != nullptr) pf.PyEval_RestoreThread(tstate);
//#endif  
}

void python_unlock(void) {
//#ifndef _WIN32	
  if(pythonInitDict != nullptr)	tstate = pf.PyEval_SaveThread();
//#endif  
}
/*
 *  extracts Absrtract Node from PyOpenSCAD Object
 */

/*
 *  parses either a PyOpenSCAD Object or an List of PyOpenScad Object and adds it to the list of supplied children, returns 1 on success
 */

int python_more_obj(std::vector<std::shared_ptr<AbstractNode>>& children, PyObject *more_obj) {
  int i, n;
  PyObject *obj;
  PyObject *dummy_dict;
  std::shared_ptr<AbstractNode> child;
  if (PyList_CHECK(more_obj)) {
    n = pf.PyList_Size(more_obj);
    for (i = 0; i < n; i++) {

      obj = pf.PyList_GetItem(more_obj, i);
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

std::string python_version(void)
{
  std::ostringstream stream;
  stream << "Python " <<  PY_MAJOR_VERSION  <<  "."  <<  PY_MINOR_VERSION  << "." << PY_MICRO_VERSION ;
  return stream.str();
}

/*
 * same as  python_more_obj but always returns only one AbstractNode by creating an UNION operation
 */

std::shared_ptr<AbstractNode> PyOpenSCADObjectToNodeMulti(PyObject *objs,PyObject **dict)
{
  std::shared_ptr<AbstractNode> result = nullptr;
  if (Py_TYPE(objs) == &PyOpenSCADType) {
    result = ((PyOpenSCADObject *) objs)->node;
    if(result.use_count() > 2) {
	    result = result->clone();
    }
    *dict =  ((PyOpenSCADObject *) objs)->dict;
  } else if (PyList_CHECK(objs)) {
    DECLARE_INSTANCE
    auto node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION);

    int n = pf.PyList_Size(objs);
    for (int i = 0; i < n; i++) {
      PyObject *obj = pf.PyList_GetItem(objs, i);
      if(Py_TYPE(obj) ==  &PyOpenSCADType) {
        std::shared_ptr<AbstractNode> child = PyOpenSCADObjectToNode(obj,dict);
        node->children.push_back(child);
      } else return nullptr;
    }
    result=node;
    *dict = nullptr; // TODO improve
  } else if(objs == Py_NONE || objs == Py_FALSE){
    result = void_node;	  
    *dict = nullptr; // TODO improve
  } else if(objs == Py_TRUE){
    result = full_node;	  
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
	
  PyObject *maindict = pf.PyModule_GetDict(pythonMainModule.get());
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  std::ostringstream stream;
//  python_hierdump(stream, node);
  std::string code = stream.str();
  while (pf.PyDict_Next(maindict, &pos, &key, &value)) {
    if(value->ob_type != &PyOpenSCADType) continue;
    std::shared_ptr<AbstractNode> testnode = ((PyOpenSCADObject *) value)->node;
    if(testnode != node) continue;
    PyObject* key1 = pf.PyUnicode_AsEncodedString(key, "utf-8", "~");
    if(key1 == nullptr) continue;
    const char *key_str =  pf.PyBytes_AsString(key1);
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

int python_numberval(PyObject *number, double *result, int *flags, int flagor)
{
  if(number == nullptr) return 1;
  if(number == Py_FALSE) return 1;
  if(number == Py_TRUE) return 1;
  if(number == Py_NONE) return 1;
  if (PyFloat_CHECK(number)) {
    *result = pf.PyFloat_AsDouble(number);
    return 0;
  }
  if (PyLong_CHECK(number)) {
    *result = pf.PyLong_AsLong(number);
    return 0;
  }
  if( number->ob_type == &PyDataType  && flags != nullptr) {
    *flags |= flagor;
     *result = PyDataObjectToValue(number);    
     return 0;
  }
  if (PyUnicode_CHECK(number) && flags != nullptr) {
    PyObjectUniquePtr str( pf.PyUnicode_AsEncodedString(number, "utf-8", "~"), PyObjectDeleter);
    const char *str1 = pf.PyBytes_AsString(str.get());
    sscanf(str1,"%lf",result);
    if(flags != nullptr) *flags |= flagor;
    return 0;
  }
  return 1;
}

std::vector<int>  python_intlistval(PyObject *list)
{
  std::vector<int> result;	
  PyObject *item;
  if( PyLong_CHECK(list)) {
    result.push_back(pf.PyLong_AsLong(list));
  }
  if( PyList_CHECK(list)) {
    for(int i=0;i<pf.PyList_Size(list); i++) {
      item = pf.PyList_GetItem(list, i);	    
      if( PyLong_CHECK(item)) {
        result.push_back(pf.PyLong_AsLong(item));
      }
    }	    
  }
  return result;
}
/*
 * Tries to extract an 3D vector out of a python list
 */

int python_vectorval(PyObject *vec, int minval, int maxval, double *x, double *y, double *z, double *w, int *flags)
{
  if(flags != nullptr) *flags = 0;
  if (PyList_CHECK(vec)) {
    if(pf.PyList_Size(vec) < minval || pf.PyList_Size(vec) > maxval) return 1;
    	  
    if (pf.PyList_Size(vec) >= 1) {
      if (python_numberval(pf.PyList_GetItem(vec, 0), x, flags, 1)) return 1;
    }
    if (pf.PyList_Size(vec) >= 2) {
      if (python_numberval(pf.PyList_GetItem(vec, 1), y, flags, 2)) return 1;
    }
    if (pf.PyList_Size(vec) >= 3) {
      if (python_numberval(pf.PyList_GetItem(vec, 2), z, flags, 4)) return 1;
    }
    if (pf.PyList_Size(vec) >= 4 && w != NULL) {
      if (python_numberval(pf.PyList_GetItem(vec, 3), w, flags, 8)) return 1;
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

std::vector<Vector3d> python_vectors(PyObject *vec, int mindim, int maxdim, int *dragflags) 
{
  std::vector<Vector3d> results;	
  if (PyList_CHECK(vec)) {
    // check if its a valid vec<Vector3d>
    int valid=1;
    for(int i=0;valid && i<pf.PyList_Size(vec);i++) {
      PyObject *item = pf.PyList_GetItem(vec,i);
      if(!PyList_CHECK(item)) valid=0;
    }	    
    if(valid) {
      for(int j=0;valid && j<pf.PyList_Size(vec);j++) {
        Vector3d result(0,0,0);	
        PyObject *item = pf.PyList_GetItem(vec,j);
        if(pf.PyList_Size(item) >= mindim && pf.PyList_Size(item) <= maxdim) {	  
          for(int i=0;i<pf.PyList_Size(item);i++) {
            if (pf.PyList_Size(item) > i) {
              if (python_numberval(pf.PyList_GetItem(item, i), &result[i])) return results; // Error
            }
          }	
        }  
	results.push_back(result);
      }	
      return results;
    }
    Vector3d result(0,0,0);	
    if(pf.PyList_Size(vec) >= mindim && pf.PyList_Size(vec) <= maxdim) {	  
      for(int i=0;i<pf.PyList_Size(vec);i++) {
        if (pf.PyList_Size(vec) > i) {
          if (python_numberval(pf.PyList_GetItem(vec, i), &result[i],dragflags, 1<<i)) return results; // Error
        }
      }	
    }  
    results.push_back(result);
  }
  Vector3d result(0,0,0);	
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
  PyObject *mainModule = pf.PyImport_AddModule("__main__");
  if (mainModule == nullptr) return;
  fn=0;
  fa=12;
  fs=2;

  if(pf.PyObject_HasAttrString(mainModule,"fn")) {
    PyObjectUniquePtr varFn(pf.PyObject_GetAttrString(mainModule, "fn"),PyObjectDeleter);
    if (varFn.get() != nullptr){
      fn = pf.PyFloat_AsDouble(varFn.get());
    }
  }  

  if(pf.PyObject_HasAttrString(mainModule,"fa")) {
    PyObjectUniquePtr varFa(pf.PyObject_GetAttrString(mainModule, "fa"),PyObjectDeleter);
    if (varFa.get() != nullptr){
      fa = pf.PyFloat_AsDouble(varFa.get());
    }
  }

  if(pf.PyObject_HasAttrString(mainModule,"fs")) {
    PyObjectUniquePtr varFs(pf.PyObject_GetAttrString(mainModule, "fs"),PyObjectDeleter);
    if (varFs.get() != nullptr){
      fs = pf.PyFloat_AsDouble(varFs.get());
    }
  }  
}

/*
 * Type specific init function. nothing special here
 */

static int PyOpenSCADInit(PyOpenSCADObject *self, PyObject *args, PyObject *kwds)
{
  (void)self;
  (void)args;
  (void)kwds;
  return 0;
}
Outline2d python_getprofile(void *v_cbfunc, int fn, double arg)
{
	PyObject *cbfunc = (PyObject *) v_cbfunc;
	Outline2d result;
	if(pythonInitDict == NULL)  initPython(PlatformUtils::applicationPath(),"", 0.0);
	PyObject* args = pf.PyTuple_Pack(1,pf.PyFloat_FromDouble(arg));
	PyObject* polygon = pf.PyObject_CallObject(cbfunc, args);
        Py_XDECREF_(args);
	if(polygon == NULL) { // TODO fix
		for(unsigned int i=0;i < (unsigned int) fn;i++) {
			double ang=360.0*(i/(double) fn);
			PyObject* args = pf.PyTuple_Pack(2,pf.PyFloat_FromDouble(arg),pf.PyFloat_FromDouble(ang));
			Py_XINCREF_(args);
			PyObject* pypt = pf.PyObject_CallObject(cbfunc, args);
			double r=pf.PyFloat_AsDouble(pypt);
			if(r < 0) r=-r;  // TODO who the hell knows, why this is needed
			double ang1=ang*3.1415/180.0;
			double x=r*cos(ang1);
			double y=r*sin(ang1);
			result.vertices.push_back(Vector2d(x,y));
		}
	} else if(polygon && PyList_CHECK(polygon)) {
		unsigned int n=pf.PyList_Size(polygon);
		for(unsigned int i=0;i < n;i++) {
			PyObject *pypt = pf.PyList_GetItem(polygon, i);
			if(PyList_CHECK(pypt) && pf.PyList_Size(pypt) == 2) {
				double x=pf.PyFloat_AsDouble(pf.PyList_GetItem(pypt, 0));
				double y=pf.PyFloat_AsDouble(pf.PyList_GetItem(pypt, 1));
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
	PyObject* args = pf.PyTuple_Pack(1,pf.PyFloat_FromDouble(arg));
	PyObject* funcresult = pf.PyObject_CallObject(cbfunc, args);
	Py_XDECREF_(args);
	if(funcresult)
		result=pf.PyFloat_AsDouble(funcresult);
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
	return Py_NONE;
      case Value::Type::BOOL:
	return val.toBool()?Py_TRUE:Py_FALSE;
      case Value::Type::NUMBER:
	return  pf.PyFloat_FromDouble(val.toDouble());
      case Value::Type::STRING:
	return pf.PyUnicode_FromString(val.toString().c_str());
      case Value::Type::VECTOR:
	{
	  const VectorType& vec = val.toVector();
  	  PyObject *result=pf.PyList_New(vec.size());
	  for(size_t j=0;j<vec.size();j++)
		pf.PyList_SetItem(result,j,python_fromopenscad(vec[j]));
	  return result;
	}
//TODO  more types RANGE, OBJECT, FUNCTION
      default:
	return Py_NONE;
    }
    return Py_NONE;
}
void python_catch_error(std::string &errorstr)
{
    PyObject *pyExcType;
    PyObject *pyExcValue;
    PyObject *pyExcTraceback;
    pf.PyErr_Fetch(&pyExcType, &pyExcValue, &pyExcTraceback);
    pf.PyErr_NormalizeException(&pyExcType, &pyExcValue, &pyExcTraceback);
    if(pyExcType != nullptr) Py_XDECREF_(pyExcType);

    if(pyExcValue != nullptr){
      PyObjectUniquePtr str_exc_value( pf.PyObject_Repr(pyExcValue), PyObjectDeleter);
      PyObjectUniquePtr pyExcValueStr( pf.PyUnicode_AsEncodedString(str_exc_value.get(), "utf-8", "~"), PyObjectDeleter);
      const char *suberror = pf.PyBytes_AsString(pyExcValueStr.get());
      if(suberror != nullptr) errorstr +=  suberror;
      Py_XDECREF_(pyExcValue);
    }
    if(pyExcTraceback != nullptr) {
      auto *tb_o = (PyTracebackObject *)pyExcTraceback;
      int line_num = tb_o->tb_lineno;
      errorstr += " in line ";
      errorstr += std::to_string(line_num);
      Py_XDECREF_(pyExcTraceback);
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
      PyObject* methodobj = pf.PyUnicode_FromString(methodname.c_str());
      pFunc = pf.PyObject_GenericGetAttr((PyObject *) python_class.ptr,methodobj);
    }
  }
  if(!pFunc) {
    PyObject *maindict = pf.PyModule_GetDict(pythonMainModule.get());

    // search the function in all modules
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (pf.PyDict_Next(maindict, &pos, &key, &value)) {
      PyObject *module = pf.PyObject_GetAttrString(pythonMainModule.get(), pf.PyUnicode_AsUTF8(key));
      if(module != nullptr){
        PyObject *moduledict = pf.PyModule_GetDict(module);
        Py_XDECREF_(module);
        if(moduledict != nullptr) {
          pFunc = pf.PyDict_GetItemString(moduledict, name.c_str());
         if(pFunc != nullptr) break;
        } 
      }
    }  
  }
  if (!pFunc) {
    errorstr="Function not found";    	  
    return nullptr;
  }
  if (!pf.PyCallable_Check(pFunc)) {
    errorstr="Function not callable";    	  
    return nullptr;
  }

  PyObject *args = pf.PyTuple_New(op_args.size());
  for(unsigned int i=0;i<op_args.size();i++)
  {
    Assignment *op_arg=op_args[i].get();

    std::shared_ptr<Expression> expr=op_arg->getExpr();
    PyObject *value= python_fromopenscad( expr.get()->evaluate(cxt));
    if(value != nullptr) {
      pf.PyTuple_SetItem(args, i, value);
    }

  }
  PyObject* funcresult = pf.PyObject_CallObject(pFunc, args);
  Py_XDECREF_(args);

  errorstr="";
  if(funcresult == nullptr) {
    python_catch_error(errorstr);	  
    pf.PyErr_SetString(pf.PyExc_TypeError, errorstr.c_str());

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
    Py_XDECREF_(funcresult);
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
  if(PyList_CHECK(arg)) {
    VectorType vec(nullptr);
    for(int i=0;i<pf.PyList_Size(arg);i++) {
      PyObject *item=pf.PyList_GetItem(arg,i);
      int suberror;
      vec.emplace_back(python_convertresult(item,suberror));
      error |= suberror;
    }
    return std::move(vec);
  } else if(PyFloat_CHECK(arg)) { return { pf.PyFloat_AsDouble(arg) }; 
  } else if(arg == Py_FALSE) { return false;
  } else if(arg == Py_TRUE) {  return true; 
  } else if(PyLong_CHECK(arg))  { return { (double) pf.PyLong_AsLong(arg) }; }
  else if(PyUnicode_CHECK(arg)) {
    auto str = std::string(pf.PyUnicode_AsUTF8(arg));
    return { str } ;
  } else if(arg == Py_NONE) { return Value::undefined.clone(); 
  } else if(arg->ob_type->tp_base == pf.PyBaseObject_Type) {
	  Py_XINCREF_(arg);
	  return PythonClassType(arg);
  } else {
    printf("unsupported type %s\n",arg->ob_type->tp_base->tp_name);	  
    pf.PyErr_SetString(pf.PyExc_TypeError, "Unsupported function result\n");
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
    pf.PyErr_SetString(pf.PyExc_TypeError, errorstr.c_str());
    return Value::undefined.clone();
  }
  if(funcresult == nullptr) return Value::undefined.clone();

  Value res = python_convertresult(funcresult, error);
  Py_XDECREF_(funcresult);
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
  		Py_XINCREF_(obj);
		python_orphan_objs.push_back(obj);
	}
}
#endif


void initPython(const std::string& binDir, const std::string &scriptpath, double time)
{
  static bool alreadyTried=false;
  if(alreadyTried) return;  
  const auto name = "openscad-python";
  const auto exe = binDir + "/" + name;
  initDynamic();
  if(scriptpath.size() > 0) python_scriptpath = scriptpath;	
  if(pythonInitDict) { /* If already initialized, undo to reinitialize after */
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *maindict = pf.PyModule_GetDict(pythonMainModule.get());
    while (pf.PyDict_Next(maindict, &pos, &key, &value)) {
      PyObjectUniquePtr key_(pf.PyUnicode_AsEncodedString(key, "utf-8", "~"), PyObjectDeleter);
      if(key_ == nullptr) continue;
      const char *key_str =  pf.PyBytes_AsString(key_.get());
      if(key_str == nullptr) continue;
      if (std::find(std::begin(pythonInventory), std::end(pythonInventory), key_str) == std::end(pythonInventory))
      {
        if(strlen(key_str) < 4 || strncmp(key_str,"stat",4) != 0){	      
          pf.PyDict_DelItemString(maindict, key_str);
	}  
      }
      // bug in  PyDict_GetItemString, thus iterating
      if(strcmp(key_str,"sys") == 0) {
        PyObject *sysdict = pf.PyModule_GetDict(value);
	if(sysdict == nullptr) continue;
	// get builtin_module_names
        PyObject *key1, *value1;
        Py_ssize_t pos1 = 0;
        while (pf.PyDict_Next(sysdict, &pos1, &key1, &value1)) {
          PyObject *key1_ = pf.PyUnicode_AsEncodedString(key1, "utf-8", "~");
          if(key1_ == nullptr) continue;
          const char *key1_str =  pf.PyBytes_AsString(key1_);
          if(strcmp(key1_str,"modules") == 0) {
            PyObject *key2, *value2;
            Py_ssize_t pos2 = 0;
            while (pf.PyDict_Next(value1, &pos2, &key2, &value2)) {
              PyObject *key2_ =pf.PyUnicode_AsEncodedString(key2, "utf-8", "~");
              if(key2_ == nullptr) continue;
              const char *key2_str =  pf.PyBytes_AsString(key2_);
	      if(key2_str == nullptr) continue;
	      if(!PyModule_CHECK(value2)) continue;

	      PyObject *modrepr = pf.PyObject_Repr(value2);
	      PyObject* modreprobj = pf.PyUnicode_AsEncodedString(modrepr, "utf-8", "~");
              const char *modreprstr = pf.PyBytes_AsString(modreprobj);
	      if(modreprstr == nullptr) continue;
	      if(strstr(modreprstr,"(frozen)") != nullptr) continue;
	      if(strstr(modreprstr,"(built-in)") != nullptr) continue;
	      if(strstr(modreprstr,"/encodings/") != nullptr) continue;
	      if(strstr(modreprstr,"_frozen_") != nullptr) continue;
	      if(strstr(modreprstr,"site-packages") != nullptr) continue;
	      if(strstr(modreprstr,"usr/lib") != nullptr) continue;

//  PyObject *mod_dict = pf.PyModule_GetDict(value2);
//  PyObject *loader = pf.PyDict_GetItemString(mod_dict,"__loader__");
//  PyObject *loaderrepr = PyObject_Repr(loader);
//  PyObject* loaderreprobj = pf.PyUnicode_AsEncodedString(loaderrepr, "utf-8", "~");
//  const char *loaderreprstr = pf.PyBytes_AsString(loaderreprobj);
//  if(strstr(loaderreprstr, "ExtensionFileLoader") != nullptr) continue; // dont delete extension files

              pf.PyDict_DelItem(value1, key2);

	    }
          }
        }
      }
    }
  } else {
    PyPreConfig preconfig;
    pf.PyPreConfig_InitPythonConfig(&preconfig);
    pf.Py_PreInitialize(&preconfig);
//    PyEval_InitThreads(); // https://stackoverflow.com/questions/47167251/pygilstate-ensure-causing-deadlock

#ifdef HAVE_PYTHON_YIELD
    set_object_callback(openscad_object_callback);
#endif
    pf.PyImport_AppendInittab("openscad", &PyInit_openscad);
    pf.PyImport_AppendInittab("libfive", &PyInit_data);
    PyConfig config;
    pf.PyConfig_InitPythonConfig(&config);

    std::string sep = "";
    std::ostringstream stream;
#ifdef _WIN32
    char sepchar = ';';
    sep = sepchar;
    stream << PlatformUtils::applicationPath() << "\\..\\libraries\\python";
#else
    char sepchar = ':';
    const auto pythonXY = "python" + std::to_string(PY_MAJOR_VERSION) + "." + std::to_string(PY_MINOR_VERSION);
    const std::array<std::string, 5> paths = {
        "../libraries/python",
        "../lib/" + pythonXY,
        "../python/lib/" + pythonXY,
        "../Frameworks/" + pythonXY,
        "../Frameworks/" + pythonXY + "/site-packages",
    };
    for (const auto& path : paths) {
        const auto p = fs::path(PlatformUtils::applicationPath() + fs::path::preferred_separator + path);
        if (fs::is_directory(p)) {
            stream << sep << fs::absolute(p).generic_string();
            sep = sepchar;
        }
    }
#endif
    fs::path scriptfile(python_scriptpath);
    stream << sep << PlatformUtils::userLibraryPath();
    stream << sep << scriptfile.parent_path().string();
    stream << sepchar << ".";
    pf.PyConfig_SetBytesString(&config, &config.pythonpath_env, stream.str().c_str());
    
    pf.PyConfig_SetBytesString(&config, &config.program_name, name);
    pf.PyConfig_SetBytesString(&config, &config.executable, exe.c_str());

    PyStatus status = pf.Py_InitializeFromConfig(&config);
    if (pf.PyStatus_Exception(status)) {
      alreadyTried=true;	    
      LOG( message_group::Error, "Python %1$lu.%2$lu.%3$lu not found. Is it installed ?",PY_MAJOR_VERSION, PY_MINOR_VERSION, PY_MICRO_VERSION);
      return;
    }
    pf.PyConfig_Clear(&config);

    pythonMainModule.reset(pf.PyImport_AddModule("__main__"));
    pythonMainModuleInitialized = pythonMainModule != nullptr;
    pythonInitDict.reset(pf.PyModule_GetDict(pythonMainModule.get()));
    pythonRuntimeInitialized = pythonInitDict != nullptr;
    PyInit_PyOpenSCAD();
    PyInit_PyData();
    pf.PyRun_String("from builtins import *\n", Py_file_input, pythonInitDict.get(), pythonInitDict.get());
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    PyObject *maindict = pf.PyModule_GetDict(pythonMainModule.get());
    while (pf.PyDict_Next(maindict, &pos, &key, &value)) {
      PyObjectUniquePtr key1(pf.PyUnicode_AsEncodedString(key, "utf-8", "~"), PyObjectDeleter);
      const char *key_str =  pf.PyBytes_AsString(key1.get());
      if(key_str != NULL) pythonInventory.push_back(key_str);
    }

  }
  std::ostringstream stream;
  stream << "t=" << time << "\nphi=" << 2*G_PI*time;
  pf.PyRun_StringFlags(stream.str().c_str(), Py_file_input, pythonInitDict.get(), pythonInitDict.get(), nullptr);
  customizer_parameters_finished = customizer_parameters;
  customizer_parameters.clear();
  python_result_handle.clear();
  nodes_hold.clear();
  DECLARE_INSTANCE
  void_node = std::make_shared<CubeNode>(instance); // just placeholders
  full_node = std::make_shared<CubeNode>(instance); // just placeholders
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
      python_show_final();
      finishDynamic();
}

std::string evaluatePython(const std::string & code, bool dry_run)
{
  std::string error;
  python_result_node = nullptr;
  python_result_handle.clear();
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
   def isatty(self):\n\
      return False\n\
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
stdin_bak = None\n\
stdout_bak = None\n\
stderr_bak = None\n\
";

#ifndef OPENSCAD_NOGUI  
    pf.PyRun_SimpleString(python_init_code);
#endif    
#ifdef HAVE_PYTHON_YIELD
    for(auto obj : python_orphan_objs) {
        Py_DECREF(obj);
    }
    python_orphan_objs.clear();
#endif
    PyObjectUniquePtr result(nullptr, PyObjectDeleter);
    result.reset(pf.PyRun_String(code.c_str(), Py_file_input, pythonInitDict.get(), pythonInitDict.get())); /* actual code is run here */


#ifndef OPENSCAD_NOGUI
    if(result  == nullptr) {
      pf.PyErr_Print();
      error = ""; 
      python_catch_error(error);
    } 
    for(int i=0;i<2;i++)
    {
      PyObjectUniquePtr catcher(nullptr, PyObjectDeleter);
      catcher.reset( pf.PyObject_GetAttrString(pythonMainModule.get(), i==1?"catcher_err":"catcher_out"));
      if(catcher == nullptr) continue;
      PyObjectUniquePtr command_output(nullptr, PyObjectDeleter);
      command_output.reset(pf.PyObject_GetAttrString(catcher.get(), "data"));

      PyObjectUniquePtr command_output_value(nullptr,  PyObjectDeleter);
      command_output_value.reset(pf.PyUnicode_AsEncodedString(command_output.get(), "utf-8", "~"));
      const char *command_output_bytes =  pf.PyBytes_AsString(command_output_value.get());
      if(command_output_bytes != nullptr && *command_output_bytes != '\0')
      {
        if(i ==1) error += command_output_bytes; /* output to console */
        else LOG(command_output_bytes); /* error to LOG */
      }
    }
    pf.PyRun_SimpleString(python_exit_code);
#endif    
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
	if(result == Py_NONE || result == nullptr)  result = pf.PyObject_GenericGetAttr(dict,key);
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

extern "C" PyObject *PyInit_openscad(void)
{
  return pf.PyModule_Create(&OpenSCADModule);
}

PyMODINIT_FUNC PyInit_PyOpenSCAD(void)
{
  PyObject *m;

  if (pf.PyType_Ready(&PyOpenSCADType) < 0) return NULL;
  
  m =PyInit_openscad();
  if (m == NULL) return NULL;

  Py_XINCREF_((PyObject *) &PyOpenSCADType);
  pf.PyModule_AddObject(m, "openscad", (PyObject *)&PyOpenSCADType);
  return m;
}

// ----------------------------------------------
// IPython Interpreter side
// ----------------------------------------------

static PyStatus
pymain_init(void)
{
    PyStatus status;

//    if (_pf.PyStatus_EXCEPTION(status)) {
//        return status;
//    }

    PyPreConfig preconfig;
    pf.PyPreConfig_InitPythonConfig(&preconfig);
//    status = _Py_PreInitializeFromPyArgv(&preconfig, args);
//    if (_pf.PyStatus_EXCEPTION(status)) {
//        return status;
//    }

    PyConfig config;
    pf.PyConfig_InitPythonConfig(&config);

//    if (args->use_bytes_argv) {
//        status = pf.PyConfig_SetBytesArgv(&config, args->argc, args->bytes_argv);
//    }
//    else {
//        status = pf.PyConfig_SetArgv(&config, args->argc, args->wchar_argv);
//    }
//    if (_pf.PyStatus_EXCEPTION(status)) {
//        goto done;
//    }

    status = pf.Py_InitializeFromConfig(&config);
//    if (_pf.PyStatus_EXCEPTION(status)) {
//        goto done;
//    }
//    status = 0; // pf.PyStatus_Ok;

    pf.PyConfig_Clear(&config);
    return status;
}


/* Write an exitcode into *exitcode and return 1 if we have to exit Python.
   Return 0 otherwise. */
static int
pymain_run_interactive_hook(int *exitcode)
{
    PyObject *sys, *hook, *result;
    sys = pf.PyImport_ImportModule("sys");
    if (sys == NULL) {
        goto error;
    }

    hook = pf.PyObject_GetAttrString(sys, "__interactivehook__");
    Py_XDECREF_(sys);
    if (hook == NULL) {
        pf.PyErr_Clear();
        return 0;
    }

    if (pf.PySys_Audit("cpython.run_interactivehook", "O", hook) < 0) {
        goto error;
    }
    result  =  pf.PyObject_CallNoArgs(hook);
    Py_XDECREF_(hook);
    if (result == NULL) {
        goto error;
    }
    Py_XDECREF_(result);
    return 0;

error:
    pf.PySys_WriteStdout("Failed calling sys.__interactivehook__\n");
//    return pymain_err_print(exitcode);
    return 0;
}






static void
pymain_repl(int *exitcode)
{
    if (pymain_run_interactive_hook(exitcode)) {
        return;
    }
    PyCompilerFlags cf = _PyCompilerFlags_INIT;

    pf.PyRun_AnyFileFlags(stdin, "<stdin>", &cf);
}


static void
pymain_run_python(int *exitcode)
{
    PyObject *main_importer_path = NULL;
//    PyInterpreterState *interp = PyInterpreterState_Get();

    pymain_repl(exitcode);
    goto done;

//    *exitcode = pymain_exit_err_print();

done:
//    _PyInterpreterState_SetNotRunningMain(interp);
    Py_XDECREF_(main_importer_path);
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
    initPython(PlatformUtils::applicationPath(),"", 0.0);
    Py_RunMain();
    return ;
}
// -------------------------

#ifdef ENABLE_JUPYTER
void python_startjupyter(void)
{
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
    initPython(0.0);
    pf.PyRun_SimpleString(python_init_code);

    auto logger = xeus::make_console_logger(xeus::xlogger::msg_type,
                                            xeus::make_file_logger(xeus::xlogger::full, "my_log_file.log"));
    try{	
	xeus::xconfiguration config = xeus::load_configuration(python_jupyterconfig);
	std::unique_ptr<xeus::xcontext> context = xeus::make_zmq_context();
	
	// Create interpreter instance
	using interpreter_ptr = std::unique_ptr<openscad_jupyter::interpreter>;
	interpreter_ptr interpreter = interpreter_ptr(new openscad_jupyter::interpreter());
		
	// Create kernel instance and start it
	xeus::xkernel kernel(config,
                         xeus::get_user_name(),
                         std::move(context),
                         std::move(interpreter),
                         xeus::make_xserver_shell_main,
			 xeus::make_in_memory_history_manager(),
			 std::move(logger));
	
	kernel.start();
    } catch(std::exception &e) {
	printf("Exception %s during startup of jupyter\n",e.what());	    
    }
}
#endif

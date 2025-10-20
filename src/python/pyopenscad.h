#include <Python.h>
#include <memory>
#include "python_public.h"
#include "geometry/Polygon2d.h"
#include "core/node.h"
#include "core/function.h"
#include "core/ScopeContext.h"
#include "core/UserModule.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

struct PyKernel
{
  void (*PyConfig_InitPythonConfig)(PyConfig *config);
  void (*PyPreConfig_InitPythonConfig)(PyPreConfig *config);
  void (*PyConfig_Clear)(PyConfig *);

  PyStatus (*PyConfig_Read)(PyConfig *config);
  PyStatus (*PyConfig_SetBytesArgv)( PyConfig *config, Py_ssize_t argc, char * const *argv);
  PyStatus (*PyConfig_SetBytesString)( PyConfig *config, wchar_t **config_str, const char *str);
  PyStatus (*Py_PreInitialize)( const PyPreConfig *src_config);
  PyStatus (*Py_InitializeFromConfig)( const PyConfig *config);

  PyObject *(*PyType_GenericAlloc)(PyTypeObject *, Py_ssize_t);
  void (*_Py_Dealloc)(PyObject *);	
  int (*PyType_Ready)(PyTypeObject *);	

  void (*PyEval_RestoreThread)(PyThreadState *);
  PyThreadState * (*PyEval_SaveThread)(void);

  PyObject *(*PyModule_GetDict)(PyObject *);	
  const char *(*PyBytes_AsString)(PyObject *);	

  PyObject *(*PyImport_AddModule)(const char *);	
  PyObject *(*PyImport_ImportModule)(const char *);	
  PyObject *(*PyModule_Create2)(PyModuleDef*, int apiver);	

  int (*PyArg_ParseTupleAndKeywords)(PyObject *, PyObject *, const char *, char **, ...);
  void (*PyErr_SetString)(PyObject *exception, const char *string);

  int (*PyRun_AnyFileExFlags)( FILE *fp, const char *filename,        int closeit, PyCompilerFlags *flags);
  int (*PyRun_SimpleStringFlags)(const char *, PyCompilerFlags *);
  PyObject *(*PyRun_StringFlags)(const char *, int, PyObject *, PyObject *, PyCompilerFlags *);
  PyObject *(*PyObject_CallObject)(PyObject *callable, PyObject *args);

  PyObject *(*PyList_GetItem)(PyObject *, Py_ssize_t);
  void (*PyList_SetItem)(PyObject *, Py_ssize_t, PyObject *);
  PyObject *(*PyList_New)(Py_ssize_t size);
  Py_ssize_t(*PyList_Size)(PyObject *);

  PyObject *(*PyDict_New)(void);	
  PyObject *(*PyDict_SetItem)(PyObject *,PyObject *,PyObject *);	
  PyObject *(*PyDict_SetItemString)(PyObject *,const char *,PyObject *);	
  PyObject *(*PyDict_GetItem)(PyObject *,PyObject *);	
  PyObject *(*PyDict_GetItemString)(PyObject *,const char *);	
  void (*PyDict_DelItem)(PyObject *,PyObject *);	
  void (*PyDict_DelItemString)(PyObject *,const char *);	
  PyObject *(*PyDict_Next)(PyObject *, Py_ssize_t *, PyObject **,PyObject **);	

  PyObject *(*PyTuple_New)(Py_ssize_t size);
  PyObject *(*PyTuple_Pack)(Py_ssize_t, ...);

  Py_ssize_t(*PyTuple_Size)(PyObject *);
  PyObject * (*PyTuple_GetItem)(PyObject *, Py_ssize_t);
  void (*PyTuple_SetItem)(PyObject *, Py_ssize_t, PyObject *);

  PyObject * (*PyFloat_FromDouble)(double);
  double (*PyFloat_AsDouble)(PyObject*);
  long (*PyLong_AsLong)(PyObject *);
  PyObject *(*PyLong_FromLong)(long );

  PyObject* (*PyUnicode_AsEncodedString)( PyObject *unicode, const char *encoding, const char *errors   );
  PyObject* (*PyUnicode_FromString)( const char *u);
  PyObject* (*PyUnicode_FromStringAndSize)( const char *u,  Py_ssize_t size );
  const char *(*PyUnicode_AsUTF8)(PyObject *unicode);

  void (*PyErr_Fetch)(PyObject **, PyObject **, PyObject **);
  int  (*PyType_IsSubtype)(PyTypeObject *, PyTypeObject *);

  void (*Py_ExitStatusException)(PyStatus err);
  int (*PyStatus_Exception)(PyStatus err);
  int (*PyStatus_IsExit)(PyStatus err);
  PyStatus (*PyWideStringList_Append)(PyWideStringList *list, const wchar_t *item);

  int (*PyModule_AddObject)(PyObject *mod, const char *, PyObject *value);

  int (*PyCallable_Check)(PyObject *); 
  void (*PyErr_Clear)(void);
  void (*PyErr_NormalizeException)(PyObject**, PyObject**, PyObject**);
  void (*PyErr_Print)(void);
  int (*PyImport_AppendInittab)( const char *name,    PyObject* (*initfunc)(void));
  PyObject * (*PyObject_CallNoArgs)(PyObject *func);
  PyObject * (*PyObject_GenericGetAttr)(PyObject *, PyObject *);
  PyObject * (*PyObject_GetAttr)(PyObject *, PyObject *);
  PyObject * (*PyObject_GetAttrString)(PyObject *, const char *);
  int (*PyObject_HasAttr)(PyObject *, PyObject *);
  int (*PyObject_HasAttrString)(PyObject *, const char *);
  PyObject * (*PyObject_Repr)(PyObject *);
  int (*PySys_Audit)( const char *event, const char *format, ...);
  void (*PySys_WriteStdout)(const char *format, ...);

  PyObject *_Py_TrueStruct;
  PyObject *_Py_FalseStruct;
  PyObject *_Py_NoneStruct;
  PyObject *PyExc_TypeError;
  PyTypeObject *PyFunction_Type;
  PyTypeObject *PyList_Type;
  PyTypeObject *PyFloat_Type;
  PyTypeObject *PyLong_Type;
  PyTypeObject *PyUnicode_Type;
  PyTypeObject *PyTuple_Type;
  PyTypeObject *PyModule_Type;
  PyTypeObject *PyBaseObject_Type;

};

extern struct PyKernel pf;

#define Py_NONE pf._Py_NoneStruct
#define Py_TRUE pf._Py_TrueStruct
#define Py_FALSE pf._Py_FalseStruct

#define PyList_CHECK(op) (op->ob_type == pf.PyList_Type)
#define PyFloat_CHECK(op) (op->ob_type == pf.PyFloat_Type)
#define PyLong_CHECK(op) (op->ob_type == pf.PyLong_Type)
#define PyUnicode_CHECK(op) (op->ob_type == pf.PyUnicode_Type)
#define PyTuple_CHECK(op) (op->ob_type == pf.PyTuple_Type)
#define PyModule_CHECK(op) (op->ob_type == pf.PyModule_Type)

static inline void Py_XDECREF_(PyObject *op)
{
    if (op != _Py_NULL) {
   // Non-limited C API and limited C API for Python 3.9 and older access
    // directly PyObject.ob_refcnt.
      if (_Py_IsImmortal(op)) return;
      _Py_DECREF_STAT_INC();
      if (--op->ob_refcnt == 0) {
          pf._Py_Dealloc(op);
      }
	    
    }
}

static inline void Py_XINCREF_(PyObject *op)
{
    // Non-limited C API and limited C API for Python 3.9 and older access
    // directly PyObject.ob_refcnt.
    op->ob_refcnt++;
}




typedef struct {
  PyObject_HEAD std::shared_ptr<AbstractNode> node;
  PyObject *dict;
  /* Type-specific fields go here. */
} PyOpenSCADObject;

void PyObjectDeleter(PyObject *pObject);
using PyObjectUniquePtr = std::unique_ptr<PyObject, decltype(PyObjectDeleter)&>;

PyMODINIT_FUNC PyInit_PyOpenSCAD(void);

extern PyTypeObject PyOpenSCADType;

extern PyObject *python_result_obj;
extern std::vector<SelectedObject> python_result_handle;
extern void python_catch_error(std::string& errorstr);

extern bool python_active;
extern fs::path python_scriptpath;
extern std::string trusted_edit_document_name;
extern std::string untrusted_edit_document_name;
extern std::vector<std::shared_ptr<AbstractNode>> nodes_hold;
extern std::shared_ptr<AbstractNode> void_node, full_node;
bool trust_python_file(const std::string& file, const std::string& content);
PyObject *PyOpenSCADObjectFromNode(PyTypeObject *type, const std::shared_ptr<AbstractNode>& node);
std::shared_ptr<AbstractNode> PyOpenSCADObjectToNode(PyObject *object, PyObject **dict);
std::shared_ptr<AbstractNode> PyOpenSCADObjectToNodeMulti(PyObject *object, PyObject **dict);
PyTypeObject *PyOpenSCADObjectType(PyObject *objs);
int python_more_obj(std::vector<std::shared_ptr<AbstractNode>>& children, PyObject *more_obj);
Outline2d python_getprofile(void *cbfunc, int fn, double arg);
double python_doublefunc(void *cbfunc, double arg);
std::shared_ptr<AbstractNode> python_modulefunc(const ModuleInstantiation *module,
                                                const std::shared_ptr<const Context>& context,
                                                std::string& error);
std::vector<int> python_intlistval(PyObject *list);

Value python_functionfunc(const FunctionCall *call, const std::shared_ptr<const Context>& context,
                          int& error);
int python_vectorval(PyObject *vec, int minarg, int maxarg, double *x, double *y, double *z,
                     double *w = NULL, int *flags = nullptr);
std::vector<Vector3d> python_vectors(PyObject *vec, int mindim, int maxdim, int *dragflags);
int python_numberval(PyObject *number, double *result, int *flags = nullptr, int flagor = 0);
void get_fnas(double& fn, double& fa, double& fs);
void python_retrieve_pyname(const std::shared_ptr<AbstractNode>& node);
void python_build_hashmap(const std::shared_ptr<AbstractNode>& node, int level);
PyObject *python_fromopenscad(const Value& val);

extern SourceFile *osinclude_source;

PyObject *python_str(PyObject *self);

extern PyNumberMethods PyOpenSCADNumbers;
extern PyMappingMethods PyOpenSCADMapping;
extern PyMethodDef PyOpenSCADFunctions[];
extern PyMethodDef PyOpenSCADMethods[];

extern PyObjectUniquePtr pythonInitDict;
extern PyObjectUniquePtr pythonMainModule;
extern int debug_num, debug_cnt;

/***************************************************************************
* Copyright (c) 2019, Sylvain Corlay, Johan Mabille, Wolf Vollprecht       *
* Copyright (c) 2019, QuantStack                                           *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#include <Python.h>
#include "pyopenscad.h"
#include "openscad_jupyter.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <stack>
#include <cctype>

#include "xeus/xinterpreter.hpp"


void PyObjectDeleter (PyObject *pObject);
using PyObjectUniquePtr = std::unique_ptr<PyObject, const decltype(PyObjectDeleter)&>;

namespace openscad_jupyter
{
    void interpreter::configure_impl()
    {
    }

    void interpreter::execute_request_impl(send_reply_callback cb,
                                           int execution_counter,
                                           const std::string& code,
                                           xeus::execute_request_config config,
                                           nl::json user_expression)
    {
	PyObject *emptystr = PyUnicode_FromString("");
        nl::json jresult;
        try
        {
	    static int python_firstrun=1;
//            if(python_firstrun) PyRun_SimpleString("print(\"\")\n");
//	    python_firstrun=0;
            nl::json pub_data;

            std::ostringstream stream;
	    int skip_eq=0;
	    if(std::strstr(code.c_str(),"import ") != nullptr) skip_eq=1;
	    if(std::strstr(code.c_str(),"=") != nullptr) skip_eq=1;

	    if(!skip_eq) stream << "__res__=";
	    stream << code;
//	    printf("tot %s\n",stream.str().c_str());
            int status = PyRun_SimpleString(stream.str().c_str());
//	    printf("status=%d\n",status);
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
                   publish_stream(i==0?"stdout":"stderr", command_output_bytes);
                }
                PyObject_SetAttrString(catcher.get(), "data", emptystr);
            }
	    if(status == 0) {
                PyObjectUniquePtr cmdresult(nullptr, PyObjectDeleter);
                cmdresult.reset( PyObject_GetAttrString(pythonMainModule.get(), "__res__"));
		if(!skip_eq && cmdresult.get() != nullptr) {
 //                   printf("result found\n");
                    if(Py_TYPE(cmdresult.get()) == &PyOpenSCADType){
                        pub_data["text/plain"] = "3D display coming soon";
		    } else if(cmdresult.get() != Py_None) {
                        PyObjectUniquePtr repr( PyObject_Repr(cmdresult.get()), PyObjectDeleter);
                        PyObjectUniquePtr reprstr( PyUnicode_AsEncodedString(repr.get(), "utf-8", "~"), PyObjectDeleter);
                        char *charstr = PyBytes_AS_STRING(reprstr.get());
			if(charstr != nullptr) {
                            pub_data["text/plain"] = charstr;
			}
		    }
                }		   
            }
            publish_execution_result(execution_counter, std::move(pub_data), nl::json::object());
            jresult["status"] = "ok";
            jresult["payload"] = nl::json::array();
            jresult["user_expressions"] = nl::json::object();
            cb(jresult);
        }
        catch (const std::runtime_error& err)
        {
            publish_stream("stderr", err.what());
            jresult["status"] = "error";
            cb(jresult);
        }
    }

    nl::json interpreter::complete_request_impl(const std::string& /*code*/, int /*cursor_pos*/)
    {
        nl::json jresult;
        jresult["status"] = "ok";
        return jresult;
    };

    nl::json interpreter::inspect_request_impl(const std::string& /*code*/,
                                               int /*cursor_pos*/,
                                               int /*detail_level*/)
    {
        nl::json jresult;
        jresult["status"] = "ok";
        return jresult;
    };

    nl::json interpreter::is_complete_request_impl(const std::string& /*code*/)
    {
        nl::json jresult;
        jresult["status"] = "complete";
        return jresult;
    };

    nl::json interpreter::kernel_info_request_impl()
    {
        nl::json result;
        result["implementation"] = "openscad";
        result["implementation_version"] = "0.1.0";
        std::string banner = "PythonSCAD\n"
        " goes Jupyter";
        result["banner"] = banner;
        result["language_info"]["name"] = "openscad";
        result["language_info"]["version"] = "";
        result["language_info"]["mimetype"] = "";
        result["language_info"]["file_extension"] = "py";
        return result;
// https://blog.jupyter.org/authoring-custom-jupyter-widgets-2884a462e724#:~:text=The%20idea%20behind%20Jupyter%20widgets,in%20the%20JavaScript%20front%2Dend.	
    }

    void interpreter::shutdown_request_impl()
    {
    }
}

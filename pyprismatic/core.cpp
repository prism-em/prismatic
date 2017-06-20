#include <Python.h>
#include "params.h"
#include "configure.h"
#include "parseInput.h"


static PyObject* pyprismatic_core_go(PyObject *self, PyObject *args){
	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
	meta.filename_atoms = "/home/aj/hdd1/clion/PRISM/SI100.XYZ";
	meta.filename_output = "/home/aj/hdd1/clion/PRISM/output_python.mrc";
	meta.algorithm = Prismatic::Algorithm::Multislice;
	// print metadata
    meta.toString();

	// configure simulation behavior
	Prismatic::configure(meta);

	// execute simulation
	Prismatic::execute_plan(meta);
	Py_RETURN_NONE;
	// if (!PyArg_ParseTuple(args, "i", &result)){
		// return NULL;
	// } else{
		// return Py_BuildValue("i",add4_local(result));	
		// return Py_BuildValue("i",add5_cuda(result));	
	// }
}

static PyMethodDef pyprismatic_core_methods[] = {
	{"go",(PyCFunction)pyprismatic_core_go, METH_VARARGS, "Execute Prismatic calculation"},
   	{NULL, NULL, 0, NULL}
};


static struct PyModuleDef module_def = {
	PyModuleDef_HEAD_INIT,"pypristmatic.core","Python wrapper for Prismatic\
	 package for fast image simulation using the PRISM and multislice\
	 algorithms in Scanning Transmission Electron Microscopy (STEM)",-1,pyprismatic_core_methods
};

PyMODINIT_FUNC PyInit_core(){
	return PyModule_Create(&module_def);
}

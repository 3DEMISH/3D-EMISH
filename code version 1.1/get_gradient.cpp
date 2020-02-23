//Python.h has all the required function definitions to manipulate the Python objects
#include <iostream>

#include "Python.h"
#include "numpy/arrayobject.h"



 //This is the function that is called from your python code
static PyObject* get_gradient_gradient(PyObject* self, PyObject* args){

    PyObject *arg1=NULL,*arg2=NULL ,*out=NULL;
    int axis;
    PyArrayObject *arr1=NULL, *arr2=NULL,*oarr=NULL;

    if (!PyArg_ParseTuple(args, "OOO!", &arg1, &arg2,
        &PyArray_Type, &out)) return NULL;

    arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
    oarr = (PyArrayObject*)PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);

    int dim0= arr1->dimensions[0]-2;
    int dim1= arr1->dimensions[1]-2;
    int dim2= arr1->dimensions[2]-2;
   
    for (int i=2; i<dim0; ++i) {
        for (int j=2; j<dim1; ++j) {
		for (int k=2; k<dim2 ;++k) {
            double *arr1_c = (double*)PyArray_GETPTR3(arr1, i, j,k);
            double *arr2_c = (double*)PyArray_GETPTR3(arr2, i, j,k);
            double *arr1_sip = (double*)PyArray_GETPTR3(arr1, i+2, j,k);
            double *arr1_sin = (double*)PyArray_GETPTR3(arr1, i-2, j,k);
	    double *arr1_sjp = (double*)PyArray_GETPTR3(arr1, i, j+2,k);
            double *arr1_sjn = (double*)PyArray_GETPTR3(arr1, i, j-2,k);
	    double *arr1_skp = (double*)PyArray_GETPTR3(arr1, i, j,k+2);
            double *arr1_skn = (double*)PyArray_GETPTR3(arr1, i, j,k-2);
	    double *arr1_ip = (double*)PyArray_GETPTR3(arr1, i+1, j,k);
            double *arr1_in = (double*)PyArray_GETPTR3(arr1, i-1, j,k);
	    double *arr1_jp = (double*)PyArray_GETPTR3(arr1, i, j+1,k);
            double *arr1_jn = (double*)PyArray_GETPTR3(arr1, i, j-1,k);
	    double *arr1_kp = (double*)PyArray_GETPTR3(arr1, i, j,k+1);
            double *arr1_kn = (double*)PyArray_GETPTR3(arr1, i, j,k-1);
	    double *arr2_ip = (double*)PyArray_GETPTR3(arr2, i+1, j,k);
            double *arr2_in = (double*)PyArray_GETPTR3(arr2, i-1, j,k);
	    double *arr2_jp = (double*)PyArray_GETPTR3(arr2, i, j+1,k);
            double *arr2_jn = (double*)PyArray_GETPTR3(arr2, i, j-1,k);
	    double *arr2_kp = (double*)PyArray_GETPTR3(arr2, i, j,k+1);
            double *arr2_kn = (double*)PyArray_GETPTR3(arr2, i, j,k-1);
            double *v = (double*)PyArray_GETPTR3(oarr, i, j,k);
            double der_i=*arr2_ip*(*arr1_sip-*arr1_c)-*arr2_in*(*arr1_c-*arr1_sin);
            double der_j=*arr2_jp*(*arr1_sjp-*arr1_c)-*arr2_jn*(*arr1_c-*arr1_sjn);
            double der_k=*arr2_kp*(*arr1_skp-*arr1_c)-*arr2_kn*(*arr1_c-*arr1_skn);
      
	    //*v =*arr2_c *0.25*(der_i+der_j+der_k)+0.25*((*arr1_ip-*arr1_in)*(*arr2_ip-*arr2_in)+(*arr1_jp-*arr1_jn)*(*arr2_jp-*arr2_jn)+(*arr1_kp-*arr1_kn)*(*arr2_kp-*arr2_kn));
	    *v =0.25*(der_i+der_j+der_k);
        
    }
    }
    }
  //value returned back to python code - another python object
  //build value here converts the C long to a python integer
 Py_DECREF(arr1);
    Py_DECREF(arr2);
    Py_DECREF(oarr);
    Py_INCREF(Py_None);
  return Py_None;
}

//This is the docstring that corresponds to our 'gradient' function.
static char get_gradient_docs[] =
    "gradient( ): gradient all elements of the list\n";

/* This table contains the relavent info mapping -
  <function-name in python module>, <actual-function>,
  <type-of-args the function expects>, <docstring associated with the function>
*/
static PyMethodDef get_gradient_funcs[] = {
    {"gradient", (PyCFunction)get_gradient_gradient, METH_VARARGS, get_gradient_docs},
    {NULL, NULL, 0, NULL}
};

/*
get_gradient is the module name, and this is the initialization block of the module.
<desired module name>, <the-info-table>, <module's-docstring>
*/
PyMODINIT_FUNC initget_gradient(void){
    Py_InitModule3("get_gradient", get_gradient_funcs,
                   "gradient all ze lists");
	import_array();
}
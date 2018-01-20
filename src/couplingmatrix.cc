#define NPY_NO_DEPRECATED_API  NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

static PyObject *module_scipy_linalg_lapack = NULL;
static PyObject *scipy_linalg_lapack_zsysv = NULL;

static PyObject* CM2S(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject *arg_matrix = NULL;
    PyObject *arg_normalizedFreq = NULL;
    static char *kwlist[] = {"matrix", "normalizedFreq", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO", kwlist, &arg_matrix, &arg_normalizedFreq))
        return NULL;

    PyObject *matrix = PyArray_FROM_OTF(arg_matrix, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (matrix == NULL){
        return NULL;
    }
    npy_double *matrix_ptr = (npy_double *)PyArray_DATA((PyArrayObject *)matrix);
    npy_intp size_matrix = PyArray_SIZE((PyArrayObject *)matrix);
    int ndim_matrix = PyArray_NDIM((PyArrayObject *)matrix);
    npy_intp *dims_matrix = PyArray_DIMS((PyArrayObject *)matrix);

    PyObject *normalizedFreq = PyArray_FROM_OTF(arg_normalizedFreq, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (normalizedFreq == NULL){
        Py_DECREF(matrix);
        return NULL;
    }
    npy_double *normalizedFreq_ptr = (npy_double *)PyArray_DATA((PyArrayObject *)normalizedFreq);
    npy_intp size_normalizedFreq = PyArray_SIZE((PyArrayObject *)normalizedFreq);
    npy_intp *dims_normalizedFreq = PyArray_DIMS((PyArrayObject *)normalizedFreq);

    npy_double R1 = 1.0, RN = 1.0, sqrt_R1_RN = npy_sqrt(R1 * RN);
    int N = (int)(*dims_matrix) - 2;

    PyObject *S11 = PyArray_ZEROS(1, dims_normalizedFreq, NPY_CDOUBLE, NPY_CORDER);
    PyObject *S21 = PyArray_ZEROS(1, dims_normalizedFreq, NPY_CDOUBLE, NPY_CORDER);
    npy_cdouble *S11_ptr = (npy_cdouble *)PyArray_DATA((PyArrayObject *)S11);
    npy_cdouble *S21_ptr = (npy_cdouble *)PyArray_DATA((PyArrayObject *)S21);
    
    PyObject *Z = PyArray_ZEROS(ndim_matrix, dims_matrix, NPY_CDOUBLE, NPY_CORDER);
    npy_cdouble *Z_ptr = (npy_cdouble *)PyArray_DATA((PyArrayObject *)Z);

    while(size_matrix--){
      (Z_ptr++)->real = *matrix_ptr++;
    }
    Z_ptr = (npy_cdouble *)PyArray_GETPTR2((PyArrayObject *)Z, 0, 0);
    Z_ptr->imag = -R1;
    Z_ptr = (npy_cdouble *)PyArray_GETPTR2((PyArrayObject *)Z, N + 1, N + 1);
    Z_ptr->imag = -RN;

    PyObject *T = NULL;
    PyObject *Y = NULL;
    npy_cdouble *Y_ptr;
    npy_intp dims_b[1] = {N + 2};
    PyObject *b = PyArray_ZEROS(1, dims_b, NPY_CDOUBLE, NPY_CORDER);
    npy_cdouble *b_ptr = (npy_cdouble *)PyArray_DATA((PyArrayObject *)b);
    b_ptr->real = 1.0;

    int i;
    while(size_normalizedFreq--){
      for (i = 0; i < N; i++){
        Z_ptr = (npy_cdouble *)PyArray_GETPTR2((PyArrayObject *)Z, i + 1, i + 1);
        matrix_ptr = (npy_double *)PyArray_GETPTR2((PyArrayObject *)matrix, i + 1, i + 1);
        Z_ptr->real = (*matrix_ptr) + (*normalizedFreq_ptr);
      }

      Py_XDECREF(T);
      T = PyObject_CallFunctionObjArgs(scipy_linalg_lapack_zsysv, Z, b, NULL);
      Y = PyTuple_GET_ITEM(T, (Py_ssize_t)2);

      Y_ptr = (npy_cdouble *)PyArray_GETPTR1((PyArrayObject *)Y, 0);
      S11_ptr->real = 1.0 - 2 * R1 * Y_ptr->imag;
      S11_ptr->imag = 2 * R1 * Y_ptr->real;

      Y_ptr = (npy_cdouble *)PyArray_GETPTR1((PyArrayObject *)Y, N + 1);
      S21_ptr->real = 2 * sqrt_R1_RN * Y_ptr->imag;
      S21_ptr->imag = -2 * sqrt_R1_RN * Y_ptr->real;

      normalizedFreq_ptr++;
      S11_ptr++;
      S21_ptr++;
    }
    
    Py_DECREF(Z);
    Py_DECREF(T);
    Py_DECREF(b);
    Py_DECREF(matrix);
    Py_DECREF(normalizedFreq);

    return Py_BuildValue("(O, O)", S11, S21);
}

static struct PyMethodDef MylibMethods[] = {
    {"CM2S", (PyCFunction)CM2S, METH_VARARGS | METH_KEYWORDS, "convert N+2 Matrix to S"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef CouplingmatrixModule = {
    PyModuleDef_HEAD_INIT,
    "couplingmatrix",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    MylibMethods
};

#ifdef __cplusplus
extern "C" {
#endif

PyMODINIT_FUNC
PyInit_couplingmatrix(void)
{
    import_array();

    module_scipy_linalg_lapack = PyImport_ImportModule("scipy.linalg.lapack");

    scipy_linalg_lapack_zsysv = PyObject_GetAttrString(module_scipy_linalg_lapack, "zsysv");
    if (!PyCallable_Check(scipy_linalg_lapack_zsysv)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be callable");
        return NULL;
    }

    Py_DECREF(module_scipy_linalg_lapack);

    return PyModule_Create(&CouplingmatrixModule);
}

#ifdef __cplusplus
}
#endif
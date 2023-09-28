#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "symnmf.h"
#define MAX_ITER 300



PyObject* get_py_list_from_matrix(double**, int, int);
static PyObject* symnmf(PyObject*, PyObject*);
static PyObject* sym(PyObject*, PyObject*);
static PyObject* ddg(PyObject*, PyObject*);
static PyObject* norm(PyObject*, PyObject*);

/* C matrix into python matrix*/

PyObject* get_py_list_from_matrix(double** A, int n, int k) {
    Py_ssize_t t;
    Py_ssize_t f;
    int i, j;
    PyObject* py_list;
    PyObject* temp_py_list;
    i = 0;
    j = 0;
    py_list = PyList_New(n);
    for (t = 0; t < n; t++) {
        temp_py_list = PyList_New(k);
        for (f = 0; f < k; f++) {
            PyList_SetItem(temp_py_list,f,PyFloat_FromDouble(A[i][j]));
            j = j + 1;
        }
        j = 0;
        PyList_SetItem(py_list,t,temp_py_list);
        i = i + 1;
    }
    return py_list;
}

/* C-Python wrapping fuctions */

static PyObject* symnmf(PyObject *self, PyObject *args) {
    Py_ssize_t i;
    Py_ssize_t j;
    PyObject* datapoints_py;
    PyObject* H_py;
    PyObject* final_H_py;
    int n, d, k, iter;
    double** datapoints;
    double** sym;
    double** ddg;
    double** norm;
    double** H;
    double** old_H;

    if(!PyArg_ParseTuple(args, "OOiii", &datapoints_py, &H_py, &n, &d, &k)) {
        printf("An Error Has Occurred");
        return Py_BuildValue("");
    }

    datapoints = create_matrix(n,d);
    H = create_matrix(n, k);
    old_H = create_matrix(n, k);


    for (i = 0; i < n; i++) { 
        for (j = 0; j < d; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(datapoints_py, i),j));
            
        }
    }
    for (i = 0; i < n; i++) { 
        for (j = 0; j < k; j++) {
            H[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(H_py, i),j));

        }
    }
    iter = 0;
    sym = get_similarity_matrix(datapoints, n, d);
    ddg = get_diag_deg_matrix(sym, n);
    norm = get_W(ddg, sym, n);
    copy_matrix_inplace(old_H, H, n, k);
    update_H(H, n, k, norm);
    while ((check_convergence(old_H, H, n, k) == 0) && iter < MAX_ITER ){
        copy_matrix_inplace(old_H, H, n, k);
        update_H(H, n, k, norm);
        iter += 1; 
    }
    final_H_py = get_py_list_from_matrix(H, n, k);
    free_matrix(old_H, n);
    free_matrix(H, n);
    free_matrix(norm, n);
    free_matrix(sym, n);
    free_matrix(ddg, n);
    free_matrix(datapoints, n);
    return final_H_py;

}

static PyObject* sym(PyObject *self, PyObject *args) {
    Py_ssize_t i;
    Py_ssize_t j;
    PyObject* datapoints_py;
    int n, d, k;
    double** datapoints;
    double** sym;

    if(!PyArg_ParseTuple(args, "Oiii", &datapoints_py, &n, &d, &k)) {
        printf("An Error Has Occurred");
        return Py_BuildValue("");
    }

    datapoints = create_matrix(n,d);

    for (i = 0; i < n; i++) { 
        for (j = 0; j < d; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(datapoints_py, i),j));
        }
    }

    sym = get_similarity_matrix(datapoints, n, d);
    print_matrix(sym, n, n);
    free_matrix(sym,n);
    return Py_BuildValue("");
}

static PyObject* ddg(PyObject *self, PyObject *args) {
    Py_ssize_t i;
    Py_ssize_t j;
    PyObject* datapoints_py;
    int n, d, k;
    double** datapoints;
    double** sym;
    double** ddg;


    if(!PyArg_ParseTuple(args, "Oiii", &datapoints_py, &n, &d, &k)) {
        printf("An Error Has Occurred");
        return Py_BuildValue("");
    }

    datapoints = create_matrix(n,d);

    for (i = 0; i < n; i++) { 
        for (j = 0; j < d; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(datapoints_py, i),j));
        }
    }


    sym = get_similarity_matrix(datapoints, n, d);
    ddg = get_diag_deg_matrix(sym, n);
    print_matrix(ddg, n, n);
    free_matrix(sym, n);
    free_matrix(ddg, n);
    return Py_BuildValue("");
}

static PyObject* norm(PyObject *self, PyObject *args) {
    Py_ssize_t i;
    Py_ssize_t j;
    PyObject* datapoints_py;
    int n, d, k;
    double** datapoints;
    double** sym;
    double** ddg;
    double** norm;
    PyObject* final_norm_py;


    if(!PyArg_ParseTuple(args, "Oiii", &datapoints_py, &n, &d, &k)) {
        printf("An Error Has Occurred");
        return Py_BuildValue("");
    }

    datapoints = create_matrix(n,d);

    for (i = 0; i < n; i++) { 
        for (j = 0; j < d; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(datapoints_py, i),j));
        }
    }

    sym = get_similarity_matrix(datapoints, n, d);
    ddg = get_diag_deg_matrix(sym, n);
    norm = get_W(ddg, sym, n);
    final_norm_py = get_py_list_from_matrix(norm, n, n);
    free_matrix(sym, n);
    free_matrix(ddg, n);
    free_matrix(norm, n);
    return final_norm_py;
}


static PyMethodDef symnmf_Methods[] = {
    {"symnmf",
    (PyCFunction) symnmf,
    METH_VARARGS,
    PyDoc_STR("symnmf wrapper")},
    {"sym",
    (PyCFunction) sym,
    METH_VARARGS,
    PyDoc_STR("sym wrapper")},
    {"ddg",
    (PyCFunction) ddg,
    METH_VARARGS,
    PyDoc_STR("ddg wrapper")},
    {"norm",
    (PyCFunction) norm,
    METH_VARARGS,
    PyDoc_STR("norm wrapper")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf",
    "symnmf Python wrapper for custom C extension library",
    -1,
    symnmf_Methods
};

PyMODINIT_FUNC PyInit_mysymnmf(void) {
    return PyModule_Create(&moduledef);
}

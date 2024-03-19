#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stddef.h>     // NULL, size_t

#include "../include/edit.h"    // varalg_edit


static PyObject*
edit(PyObject* self, PyObject* args, PyObject* kwargs)
{
    (void) self;

    size_t lhs = 0;
    size_t rhs = 40;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i|i", (char*[3]) {"lhs", "rhs", NULL}, &lhs, &rhs))
    {
        return NULL;
    } // if

    return Py_BuildValue("i", varalg_edit(lhs, rhs));
} // edit


static PyMethodDef methods[] =
{
    {"edit", (PyCFunction) (void(*)(void)) edit, METH_VARARGS | METH_KEYWORDS, "Calculate edit distance"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef algebra =
{
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "algebra",
    .m_doc = "Variant Algebra",
    .m_size = -1,
    .m_methods = methods,
};


static struct PyModuleDef lcs =
{
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "algebra.lcs",
    .m_doc = "LCS",
    .m_size = -1,
    .m_methods = methods,
};


PyMODINIT_FUNC
PyInit_algebra(void)
{
    PyObject* algebra_module = PyModule_Create(&algebra);
    if (algebra_module == NULL)
    {
        PyErr_SetString(PyExc_ImportError, "");
        return NULL;
    } // if

    PyObject* lcs_module = PyModule_Create(&lcs);
    if (lcs_module == NULL)
    {
        PyErr_SetString(PyExc_ImportError, "");
        return NULL;
    } // if

    PyDict_SetItemString(PyImport_GetModuleDict(), lcs.m_name, lcs_module);
    if (PyModule_AddObject(algebra_module, "lcs", lcs_module) < 0)
    {
        PyErr_SetString(PyExc_ImportError, "");
        return NULL;
    } // if

    return algebra_module;
} // PyInit_algebra

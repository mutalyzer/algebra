#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "../include/edit.h"    // varalg_edit


static PyObject*
edit(PyObject* self, PyObject* args, PyObject* kwargs)
{
    (void) self;

    size_t lhs = 0;
    size_t rhs = 40;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i|i:edit", (char*[3]) {"lhs", "rhs", NULL}, &lhs, &rhs))
    {
        return NULL;
    } // if

    return Py_BuildValue("i", varalg_edit(lhs, rhs));
} // edit


static PyModuleDef def =
{
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "lcs_ext",
    .m_doc = "LCS C extension",
    .m_size = -1,
    .m_methods = (PyMethodDef[])
    {
        {"edit", (PyCFunction) (void(*)(void)) edit, METH_VARARGS | METH_KEYWORDS, "Calculate the simple edit distance between two strings."},
        {NULL, NULL, 0, NULL}
    },
};


PyMODINIT_FUNC
PyInit_lcs_ext(void)
{

    return PyModule_Create(&def);
} // PyInit_lcs_ext

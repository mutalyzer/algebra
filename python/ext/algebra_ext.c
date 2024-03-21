#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t


static PyObject*
hello(PyObject* self, PyObject* args, PyObject* kwargs)
{
    (void) self;

    size_t lhs = 0;
    size_t rhs = 40;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i|i:hello", (char*[3]) {"lhs", "rhs", NULL}, &lhs, &rhs))
    {
        return NULL;
    } // if

    return Py_BuildValue("i", lhs * rhs);
} // hello


static struct PyModuleDef def =
{
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_doc = "Algebra C extension",
    .m_size = -1,
    .m_methods = (PyMethodDef[]) {
        {"hello", (PyCFunction) (void(*)(void)) hello, METH_VARARGS | METH_KEYWORDS, "Hello from parent"},
        {NULL, NULL, 0, NULL}
    },
};


PyMODINIT_FUNC
PyInit_algebra_ext(void)
{
    return PyModule_Create(&def);
} // PyInit_algebra_ext

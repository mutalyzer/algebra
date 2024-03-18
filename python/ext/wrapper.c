#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "../include/edit.h" // varalg_edit


static PyObject* edit(PyObject *self, PyObject* args, PyObject * kwargs){
    (void) self;
    size_t lhs = 0;
    size_t rhs = 40;
    static char const * keywords[] = {"lhs", "rhs", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i|i", (char **) keywords, &lhs, &rhs)){
        return NULL;
    }

    return Py_BuildValue("i", varalg_edit(lhs, rhs));
}


static PyMethodDef methods[] = {
    {"edit", (PyCFunction)(void(*)(void)) edit, METH_VARARGS | METH_KEYWORDS, "Calculate edit distance"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef algebra = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "algebra",
    .m_doc = "Variant Algebra",
    .m_size = -1,
    .m_methods = methods,
};

PyMODINIT_FUNC PyInit_algebra(void){
    return PyModule_Create(&algebra);
}

#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "../include/parser.h"      // va_parse_*
#include "../include/variant.h"     // VA_Variant


static PyObject*
parse_hgvs(PyObject* const self, PyObject* const args, PyObject* const kwargs)
{
    (void) self;

    char const* expression = NULL;
    size_t len = 0;
    char const* reference = NULL;
    size_t len_ref = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|z#:parse_hgvs", (char*[3]) {"expression", "reference", NULL}, &expression, &len, &reference, &len_ref))
    {
        return NULL;
    } // if

    VA_Variant variants[256] = {{0}};
    size_t const count = va_parse_hgvs(len, expression, 256, variants);

    return Py_BuildValue("i", count);
} // parse_hgvs


static PyObject*
parse_spdi(PyObject* const self, PyObject* const args, PyObject* const kwargs)
{
    (void) self;

    char const* expression = NULL;
    size_t len = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#:parse_spdi", (char*[2]) {"expression", NULL}, &expression, &len))
    {
        return NULL;
    } // if

    VA_Variant variants[1] = {{0}};
    size_t const count = va_parse_spdi(len, expression, variants);

    return Py_BuildValue("i", count);
} // parse_spdi


static struct PyModuleDef const def =
{
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_doc = "Algebra C extension",
    .m_size = -1,
    .m_methods = (PyMethodDef[]) {
        {"parse_hgvs", (PyCFunction) (void(*)(void)) parse_hgvs, METH_VARARGS | METH_KEYWORDS, "Parse HGVS expression"},
        {"parse_spdi", (PyCFunction) (void(*)(void)) parse_spdi, METH_VARARGS | METH_KEYWORDS, "Parse SPDI expression"},
        {NULL, NULL, 0, NULL}  // sentinel
    },
};


PyMODINIT_FUNC
PyInit_algebra_ext(void)
{
    return PyModule_Create(&def);
} // PyInit_algebra_ext

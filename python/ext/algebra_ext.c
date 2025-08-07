#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "lcsgraph.h"    // LCSgraph, LCSgraph_*
#include "variant.h"     // Variant, Variant_*

#include "../include/extract.h"     // gva_extract
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/variant.h"     // GVA_Variant
#include "../src/array.h"           // ARRAY_DESTROY, array_length


static PyObject*
canonical(PyObject* self, PyObject* args, PyObject* kwargs)
{
    LCSgraph* graph = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|:canonical",
                                     (char*[]) {"graph", NULL},
                                     &LCSgraph_Type, &graph))
    {
        return NULL;
    } // if

    GVA_Variant* variants = gva_extract(gva_std_allocator, graph->graph);

    PyObject* result = PyList_New(array_length(variants));
    if (result == NULL)
    {
        variants = ARRAY_DESTROY(gva_std_allocator, variants);
        return NULL;
    } // if

    for (size_t i = 0; i < array_length(variants); ++i)
    {
        PyObject* str = PyUnicode_FromFormat("%.*s", (int) variants[i].sequence.len, variants[i].sequence.str);
        if (str == NULL)
        {
            Py_DECREF(result);
            variants = ARRAY_DESTROY(gva_std_allocator, variants);
            return NULL;
        } // if
        PyObject* item = PyObject_CallFunction((PyObject*) &Variant_Type, "nnO", variants[i].start, variants[i].end, str);
        if (item == NULL)
        {
            Py_DECREF(str);
            Py_DECREF(result);
            variants = ARRAY_DESTROY(gva_std_allocator, variants);
            return NULL;
        } // if
        PyList_SET_ITEM(result, i, item);
    } // for
    variants = ARRAY_DESTROY(gva_std_allocator, variants);

    return result;
} // canonical


static PyModuleDef algebra_ext_module =
{
    PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_doc = PyDoc_STR("Algebra C extension"),
    .m_size = -1,
    .m_methods = (PyMethodDef[])
    {
        {
            "canonical",
            (PyCFunction) canonical,
            METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("Extract the canonical variant representation together with its supremal representation for an allele."),
        },
        {NULL},  // sentinel
    },
};


PyMODINIT_FUNC
PyInit_algebra_ext(void)
{
    PyObject* const module = PyModule_Create(&algebra_ext_module);
    if (module == NULL)
    {
        return NULL;
    } // if

    if (PyType_Ready(&LCSgraph_Type) < 0)
    {
        return NULL;
    } // if
    if (PyType_Ready(&Variant_Type) < 0)
    {
        return NULL;
    } // if

    if (PyModule_AddObjectRef(module, "LCSgraph", (PyObject*) &LCSgraph_Type) < 0)
    {
        return NULL;
    } // if
    if (PyModule_AddObjectRef(module, "Variant", (PyObject*) &Variant_Type) < 0)
    {
        return NULL;
    } // if

    return module;
} // PyInit_algebra_ext

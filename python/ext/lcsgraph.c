#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL

#include "lcsgraph.h"   // LCSgraph, LCSgraph_*
#include "variant.h"    // Variant, Variant_*

#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/types.h"       // GVA_NULL


PyObject*
LCSgraph_new(PyTypeObject* type, PyObject* args, PyObject* kwargs)
{
    (void) args;
    (void) kwargs;

    LCSgraph* const self = (LCSgraph*) type->tp_alloc(type, 0);
    if (self == NULL)
    {
        return NULL;
    } // if

    self->graph = (GVA_LCS_Graph) {NULL, NULL, NULL, {0, 0, {0, NULL}}, {0, NULL}, GVA_NULL, 0};
    return (PyObject*) self;
} // LCSgraph_new


inline void
LCSgraph_dealloc(LCSgraph* self)
{
    gva_lcs_graph_destroy(gva_std_allocator, self->graph);
    Py_TYPE(self)->tp_free((PyObject*) self);
} // LCSgraph_dealloc


PyObject*
LCSgraph_from_variants(PyObject* cls, PyObject* args, PyObject* kwargs)
{
    char const* reference = NULL;
    Py_ssize_t len_ref = 0;
    PyObject* variant_list = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#O!|:from_variants",
                                     (char* []) {"reference", "variants", NULL},
                                     &reference, &len_ref, &PyList_Type, &variant_list))
    {
        return NULL;
    } // if


    GVA_Variant* variants = gva_std_allocator.allocate(gva_std_allocator.context, NULL, 0, PyList_GET_SIZE(variant_list) * sizeof(*variants));
    if (variants == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "OOM");
        return NULL;
    } // if

    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(variant_list); ++i)
    {
        Variant* item = NULL;
        if (!PyArg_Parse(PyList_GET_ITEM(variant_list, i), "O!", &Variant_Type, &item))
        {
            variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, PyList_GET_SIZE(variant_list) * sizeof(*variants), 0);
            return NULL;
        } // if
    } // for

    LCSgraph* self = (LCSgraph*) PyObject_CallNoArgs(cls);
    if (self == NULL)
    {
        variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, PyList_GET_SIZE(variant_list) * sizeof(*variants), 0);
        return NULL;
    } // if

    self->graph = gva_lcs_graph_from_variants(gva_std_allocator, len_ref, reference, PyList_GET_SIZE(variant_list), variants);

    PySys_FormatStderr("supremal: " GVA_VARIANT_FMT_SPDI "\n", GVA_VARIANT_PRINT_SPDI("reference", self->graph.supremal));

    variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, PyList_GET_SIZE(variant_list) * sizeof(*variants), 0);
    return (PyObject*) self;
} // LCSgraph_from_variants


PyTypeObject LCSgraph_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algebra_ext.LCSgraph",
    .tp_basicsize = sizeof(LCSgraph),
    .tp_doc = PyDoc_STR("LCS graph class for storing all minimal alignments."),
    .tp_new = LCSgraph_new,
    .tp_dealloc = (destructor) LCSgraph_dealloc,
    .tp_methods = (PyMethodDef[])
    {
        {
            "from_variants",
            (PyCFunction) LCSgraph_from_variants,
            METH_CLASS | METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("Iteratively find the supremal LCS graph for an allele by repeatedly widening a range of influence."),
        },
        {NULL}  // sentinel
    },
};

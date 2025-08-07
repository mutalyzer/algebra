#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "lcsgraph.h"   // LCSgraph, LCSgraph_*
#include "variant.h"    // Variant, Variant_*

#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*, gva_edges
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // gva_string_destroy
#include "../include/types.h"       // GVA_NULL
#include "../include/variant.h"     // GVA_Variant
#include "../src/array.h"           // array_length


// DEBUG: PySys_FormatStderr


static PyObject*
LCSgraph_new(PyTypeObject* subtype, PyObject* args, PyObject* kwargs)
{
    (void) args;
    (void) kwargs;

    LCSgraph* const self = (LCSgraph*) subtype->tp_alloc(subtype, 0);
    if (self == NULL)
    {
        return NULL;
    } // if
    return (PyObject*) self;
} // LCSgraph_new


static inline void
LCSgraph_dealloc(LCSgraph* self)
{
    gva_string_destroy(gva_std_allocator, self->graph.observed);
    gva_lcs_graph_destroy(gva_std_allocator, self->graph);
    Py_TYPE(self)->tp_free((PyObject*) self);
} // LCSgraph_dealloc


static PyObject*
from_variants(PyObject* cls, PyObject* args, PyObject* kwargs)
{
    char const* reference = NULL;
    Py_ssize_t len_ref = 0;
    PyObject* variant_list = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#O!|:from_variants",
                                     (char*[]) {"reference", "variants", NULL},
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
        variants[i] = (GVA_Variant) {item->start, item->end, {item->len, item->sequence}};
    } // for

    LCSgraph* self = (LCSgraph*) PyObject_CallNoArgs(cls);
    if (self == NULL)
    {
        variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, PyList_GET_SIZE(variant_list) * sizeof(*variants), 0);
        return NULL;
    } // if

    self->graph = gva_lcs_graph_from_variants(gva_std_allocator, len_ref, reference, PyList_GET_SIZE(variant_list), variants);

    variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, PyList_GET_SIZE(variant_list) * sizeof(*variants), 0);
    return (PyObject*) self;
} // from_variants


static PyObject*
local_supremal(LCSgraph* self)
{
    PyObject* local = PyList_New(array_length(self->graph.local_supremal) - 1);
    if (local == NULL)
    {
        return NULL;
    } // if

    for (size_t i = 0; i < array_length(self->graph.local_supremal) - 1; ++i)
    {
        GVA_Variant variant;
        gva_edges(self->graph.observed.str,
            self->graph.local_supremal[i], self->graph.local_supremal[i + 1],
            i == 0, i == array_length(self->graph.local_supremal) - 2,
            &variant);

        PyObject* str = PyUnicode_FromFormat("%.*s", (int) variant.sequence.len, variant.sequence.str);
        if (str == NULL)
        {
            Py_DECREF(local);
            return NULL;
        } // if
        PyObject* item = PyObject_CallFunction((PyObject*) &Variant_Type, "nnO", variant.start, variant.end, str);
        if (item == NULL)
        {
            Py_DECREF(str);
            Py_DECREF(local);
            return NULL;
        } // if
        PyList_SET_ITEM(local, i, item);
    } // for
    return local;
} // local_supremal


static PyObject*
supremal(LCSgraph* self)
{
    PyObject* str = PyUnicode_FromFormat("%.*s", (int) self->graph.supremal.sequence.len, self->graph.supremal.sequence.str);
    if (str == NULL)
    {
        return NULL;
    } // if
    return PyObject_CallFunction((PyObject*) &Variant_Type, "nnO", self->graph.supremal.start, self->graph.supremal.end, str);
} // supremal


PyTypeObject LCSgraph_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algebra_ext.LCSgraph",
    .tp_basicsize = sizeof(LCSgraph),
    .tp_doc = PyDoc_STR("LCS graph class for storing all minimal alignments."),
    .tp_new = (newfunc) LCSgraph_new,
    .tp_dealloc = (destructor) LCSgraph_dealloc,
    .tp_methods = (PyMethodDef[])
    {
        {
            "from_variants",
            (PyCFunction) from_variants,
            METH_CLASS | METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("Iteratively find the supremal LCS graph for an allele by repeatedly widening a range of influence."),
        },
        {
            "local_supremal",
            (PyCFunction) local_supremal,
            METH_NOARGS,
            PyDoc_STR("The local supremal variant."),
        },
        {
            "supremal",
            (PyCFunction) supremal,
            METH_NOARGS,
            PyDoc_STR("The supremal variant."),
        },
        {NULL}  // sentinel
    },
};

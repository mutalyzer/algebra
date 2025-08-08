#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "lcsgraph.h"   // LCSgraph, LCSgraph_*
#include "variant.h"    // Variant, Variant_*

#include "../include/extractor.h"   // gva_canonical
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, gva_lcs_graph_*, gva_edges
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // gva_string_destroy
#include "../include/types.h"       // GVA_NULL
#include "../include/variant.h"     // GVA_Variant
#include "../src/array.h"           // ARRAY_DESTROY, array_length


static PyObject*
LCSgraph_new(PyTypeObject* subtype, PyObject* Py_UNUSED(args), PyObject* Py_UNUSED(kwargs))
{
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
LCSgraph_from_variants(PyObject* cls, PyObject* args, PyObject* kwargs)
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
} // LCSgraph_from_variants


static PyObject*
LCSgraph_canonical(LCSgraph* self)
{
    GVA_Variant* variants = gva_canonical(gva_std_allocator, self->graph);

    PyObject* result = PyList_New(array_length(variants));
    if (result == NULL)
    {
        variants = ARRAY_DESTROY(gva_std_allocator, variants);
        return NULL;
    } // if

    for (size_t i = 0; i < array_length(variants); ++i)
    {
        PyObject* str = PyUnicode_FromStringAndSize(variants[i].sequence.str, variants[i].sequence.len);
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
} // LCSgraph_canonical


static PyObject*
LCSgraph_local_supremal(LCSgraph* self)
{
    PyObject* result = PyList_New(array_length(self->graph.local_supremal) - 1);
    if (result == NULL)
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

        PyObject* str = PyUnicode_FromStringAndSize(variant.sequence.str, variant.sequence.len);
        if (str == NULL)
        {
            Py_DECREF(result);
            return NULL;
        } // if
        PyObject* item = PyObject_CallFunction((PyObject*) &Variant_Type, "nnO", variant.start, variant.end, str);
        if (item == NULL)
        {
            Py_DECREF(str);
            Py_DECREF(result);
            return NULL;
        } // if
        PyList_SET_ITEM(result, i, item);
    } // for
    return result;
} // LCSgraph_local_supremal


static PyObject*
LCSgraph_supremal(LCSgraph* self)
{
    PyObject* str = PyUnicode_FromStringAndSize(self->graph.supremal.sequence.str, self->graph.supremal.sequence.len);
    if (str == NULL)
    {
        return NULL;
    } // if
    PyObject* result = PyObject_CallFunction((PyObject*) &Variant_Type, "nnO", self->graph.supremal.start, self->graph.supremal.end, str);
    if (result == NULL)
    {
        Py_DECREF(str);
        return NULL;
    } // if
    return result;
} // LCSgraph_supremal


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
            (PyCFunction) LCSgraph_from_variants,
            METH_CLASS | METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("Iteratively find the supremal LCS graph for an allele by repeatedly widening a range of influence."),
        },
        {
            "canonical",
            (PyCFunction) LCSgraph_canonical,
            METH_NOARGS,
            PyDoc_STR("The canonical variant."),
        },
        {
            "local_supremal",
            (PyCFunction) LCSgraph_local_supremal,
            METH_NOARGS,
            PyDoc_STR("The local supremal variant."),
        },
        {
            "supremal",
            (PyCFunction) LCSgraph_supremal,
            METH_NOARGS,
            PyDoc_STR("The supremal variant."),
        },
        {NULL}  // sentinel
    },
};

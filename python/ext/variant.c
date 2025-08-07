#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, offsetof

#include "variant.h"    // Variant, Variant_*

#include "../include/std_alloc.h"   // gva_std_allocator


static PyObject*
Variant_new(PyTypeObject* subtype, PyObject* args, PyObject* kwargs)
{
    Py_ssize_t start = 0;
    Py_ssize_t end = 0;
    char const* sequence = NULL;
    Py_ssize_t len = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "nns#|:Variant",
                                     (char*[]) {"start", "end", "sequence", NULL},
                                     &start, &end, &sequence, &len))
    {
        return NULL;
    } // if

    Variant* const self = (Variant*) subtype->tp_alloc(subtype, 0);
    if (self == NULL)
    {
        return NULL;
    } // if

    self->start = start;
    self->end = end;
    self->len = len;
    self->sequence = sequence;
    return (PyObject*) self;
} // Variant_new


static inline void
Variant_dealloc(Variant* self)
{
    Py_TYPE(self)->tp_free((PyObject*) self);
} // Variant_dealloc


static inline PyObject*
Variant_repr(Variant* self)
{
    return PyUnicode_FromFormat("%zu:%zu/%.*s", self->start, self->end, (int) self->len, self->sequence);
} // Variant_repr


PyTypeObject Variant_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algebra_ext.Variant",
    .tp_basicsize = sizeof(Variant),
    .tp_doc = PyDoc_STR("Variant class for deletion/insertions."),
    .tp_new = (newfunc) Variant_new,
    .tp_dealloc = (destructor) Variant_dealloc,
    .tp_repr = (reprfunc) Variant_repr,
    .tp_members = (PyMemberDef[])
    {
        {
            "start",
            Py_T_PYSSIZET,
            offsetof(Variant, start),
            Py_READONLY,
            PyDoc_STR("The start position (included) of the deleted part (zero-based)."),
        },
        {
            "end",
            Py_T_PYSSIZET,
            offsetof(Variant, end),
            Py_READONLY,
            PyDoc_STR("The end position (not included) of the deleted part."),
        },
        {
            "sequence",
            Py_T_STRING,
            offsetof(Variant, sequence),
            Py_READONLY,
            PyDoc_STR("The inserted sequence."),
        },
        {NULL},  // sentinel
    },
};

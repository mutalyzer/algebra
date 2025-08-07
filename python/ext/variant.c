#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL

#include "variant.h"    // Variant, Varaint_*


PyObject*
Variant_new(PyTypeObject* type, PyObject* args, PyObject* kwargs)
{
    Py_ssize_t start = 0;
    Py_ssize_t end = 0;
    char const* sequence = NULL;
    Py_ssize_t len = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "nns#:Variant",
                                     (char* []) {"start", "end", "sequence", NULL},
                                     &start, &end, &sequence, &len))
    {
        return NULL;
    } // if

    Variant* const self = (Variant*) type->tp_alloc(type, 0);
    if (self == NULL)
    {
        return NULL;
    } // if

    self->variant = (GVA_Variant) {start, end, {len, sequence}};
    return (PyObject*) self;
} // Variant_new


inline void
Variant_dealloc(Variant* self)
{
    Py_TYPE(self)->tp_free((PyObject*) self);
} // Variant_dealloc


PyTypeObject Variant_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algebra_ext.Variant",
    .tp_basicsize = sizeof(Variant),
    .tp_doc = PyDoc_STR("Variant class for deletion/insertions."),
    .tp_new = Variant_new,
    .tp_dealloc = (destructor) Variant_dealloc,
};

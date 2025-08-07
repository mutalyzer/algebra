#ifndef EXT_VARIANT_H
#define EXT_VARIANT_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include "../include/variant.h"     // GVA_Variant


typedef struct
{
    PyObject_HEAD
    GVA_Variant variant;
} Variant;


PyObject*
Variant_new(PyTypeObject* type, PyObject* args, PyObject* kwargs);


void
Variant_dealloc(Variant* self);


extern PyTypeObject Variant_Type;


#endif // EXT_VARIANT_H

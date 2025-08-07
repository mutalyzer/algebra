#ifndef EXT_VARIANT_H
#define EXT_VARIANT_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*


typedef struct
{
    PyObject_HEAD
    Py_ssize_t  start;
    Py_ssize_t  end;
    Py_ssize_t  len;
    char const* sequence;
} Variant;


extern PyTypeObject Variant_Type;


#endif // EXT_VARIANT_H

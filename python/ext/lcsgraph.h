#ifndef EXT_LCSGRAPH_H
#define EXT_LCSGRAPH_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include "../include/lcs_graph.h"   // GVA_LCS_Graph


typedef struct
{
    PyObject_HEAD
    GVA_LCS_Graph graph;
} LCSgraph;


PyObject*
LCSgraph_new(PyTypeObject* type, PyObject* args, PyObject* kwargs);


void
LCSgraph_dealloc(LCSgraph* self);


PyObject*
LCSgraph_from_variants(PyObject* cls, PyObject* args, PyObject* kwargs);


extern PyTypeObject LCSgraph_Type;


#endif // EXT_LCSGRAPH_H

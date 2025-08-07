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


extern PyTypeObject LCSgraph_Type;


#endif // EXT_LCSGRAPH_H

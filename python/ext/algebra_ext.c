#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL

#include "lcsgraph.h"   // LCSgraph, LCSgraph_*
#include "variant.h"    // Variant, Variant_*

#include "../include/compare.h"     // gva_compare_graphs
#include "../include/relations.h"   // GVA_RELATION_LABELS, GVA_Relation
#include "../include/std_alloc.h"   // gva_std_allocator


// DEBUG: PySys_FormatStderr


static PyObject*
algebra_ext_compare(PyObject* self, PyObject* args, PyObject* kwargs)
{
    char const* reference = NULL;
    Py_ssize_t len_ref = 0;
    LCSgraph* lhs = NULL;
    LCSgraph* rhs = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#O!O!|:compare",
                                     (char*[]) {"reference", "lhs", "rhs", NULL},
                                     &reference, &len_ref, &LCSgraph_Type, &lhs, &LCSgraph_Type, &rhs))
    {
        return NULL;
    } // if

    GVA_Relation const relation = gva_compare_graphs(gva_std_allocator, len_ref, reference, lhs->graph, rhs->graph);

    PyObject* relations = PyObject_GetAttrString(self, "Relation");
    if (relations == NULL)
    {
        return NULL;
    } // if

    PyObject* result = PyObject_GetAttrString(relations, GVA_RELATION_LABELS[relation]);
    Py_DECREF(relations);
    return result;
} // algebra_ext_compare


static PyModuleDef algebra_ext_module =
{
    PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_doc = PyDoc_STR("Algebra C extension"),
    .m_size = -1,
    .m_methods = (PyMethodDef[])
    {
        {
            "compare",
            (PyCFunction) algebra_ext_compare,
            METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("Compare two LCS graphs."),
        },
        {NULL}  // sentinel
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

    PyObject* relations = PyDict_New();
    if (relations == NULL)
    {
        return NULL;
    } // if
    for (size_t i = 0; i < sizeof(GVA_RELATION_LABELS) / sizeof(*GVA_RELATION_LABELS); ++i)
    {
        PyDict_SetItemString(relations, GVA_RELATION_LABELS[i], PyUnicode_FromString(GVA_RELATION_LABELS[i]));
    } // for

    PyObject* module_enum = PyImport_ImportModule("enum");
    if (module_enum == NULL)
    {
        return NULL;
    } // if

    PyObject* enum_class = PyObject_CallMethod(module_enum, "Enum", "sO", "Relation", relations);
    if (enum_class == NULL)
    {
        return NULL;
    } // if

    Py_DECREF(module_enum);
    Py_DECREF(relations);

    if (PyModule_AddObject(module, "Relation", enum_class) < 0)
    {
        return NULL;
    } // if

    return module;
} // PyInit_algebra_ext

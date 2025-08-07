#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stddef.h>     // NULL, size_t

#include "lcsgraph.h"    // LCSgraph, LCSgraph_*
#include "variant.h"     // Variant, Variant_*


static PyModuleDef algebra_ext_module =
{
    PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_doc = PyDoc_STR("Algebra C extension"),
    .m_size = -1,
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

    return module;
} // PyInit_algebra_ext

#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Py*

#include <stdbool.h>    // bool
#include <stddef.h>     // NULL, offsetof, size_t
#include <string.h>     // memcpy

#include "../include/compare.h"     // gva_compare_graphs
#include "../include/extractor.h"   // gva_canonical
#include "../include/lcs_graph.h"   // GVA_LCS_Graph, GVA_Node, gva_lcs_graph_*, gva_edges
#include "../include/std_alloc.h"   // gva_std_allocator
#include "../include/string.h"      // gva_string_destroy
#include "../include/types.h"       // GVA_NULL, gva_uint
#include "../include/variant.h"     // GVA_VARIANT_*, GVA_Variant, gva_variant_*
#include "../src/array.h"           // ARRAY_DESTROY, array_length


// DEBUG: PySys_FormatStderr


typedef struct
{
    PyObject_VAR_HEAD
    GVA_LCS_Graph graph;
} LCSgraph_Object;


typedef struct
{
    PyObject_VAR_HEAD
    size_t           node;
    gva_uint         edge;
    LCSgraph_Object* graph;
} LCSgraph_Edge_Iterator;


typedef struct
{
    PyObject_VAR_HEAD
    Py_ssize_t start;
    Py_ssize_t end;
    Py_ssize_t len;
    char*      sequence;
} Variant_Object;


typedef struct
{
    PyTypeObject* LCSgraph_Type;
    PyTypeObject* LCSgraph_Edge_Iterator_Type;
    PyTypeObject* Variant_Type;
    PyObject*     Relations;
} algebra_ext_state;


static inline algebra_ext_state*
get_algebra_ext_state(PyObject* module)
{
    return (algebra_ext_state*) PyModule_GetState(module);
} // get_algebra_ext_state


static inline PyObject*
LCSgraph_new(PyTypeObject* subtype, PyObject* Py_UNUSED(args), PyObject* Py_UNUSED(kwargs))
{
    return subtype->tp_alloc(subtype, 0);
} // LCSgraph_new


static inline void
LCSgraph_dealloc(LCSgraph_Object* self)
{
    // FIXME: obeserved ownership
    gva_string_destroy(gva_std_allocator, self->graph.observed);
    gva_lcs_graph_destroy(gva_std_allocator, self->graph);
    Py_TYPE(self)->tp_free((PyObject*) self);
} // LCSgraph_dealloc


static PyObject*
LCSgraph_from_variants(PyObject* cls, PyObject* args, PyObject* kwargs)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule((PyTypeObject*) cls));
    if (state == NULL)
    {
        return NULL;
    } // if

    char const* reference = NULL;
    Py_ssize_t len_ref = 0;
    PyObject* variant_list = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#O!|from_variants",
                                     (char*[]) {"reference", "variants", NULL},
                                     &reference, &len_ref, &PyList_Type, &variant_list))
    {
        return NULL;
    } // if

    size_t const n = PyList_GET_SIZE(variant_list);
    GVA_Variant* variants = gva_std_allocator.allocate(gva_std_allocator.context, NULL, 0, n * sizeof(*variants));
    if (variants == NULL)
    {
        return PyErr_NoMemory();
    } // if

    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(variant_list); ++i)
    {
        Variant_Object* item = NULL;
        if (!PyArg_Parse(PyList_GET_ITEM(variant_list, i), "O!", state->Variant_Type, &item))
        {
            variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, n * sizeof(*variants), 0);
            return NULL;
        } // if
        variants[i] = (GVA_Variant) {item->start, item->end, {item->len, item->sequence}};
    } // for

    LCSgraph_Object* const self = (LCSgraph_Object*) PyObject_CallNoArgs(cls);
    if (self == NULL)
    {
        variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, n * sizeof(*variants), 0);
        return NULL;
    } // if

    self->graph = gva_lcs_graph_from_variants(gva_std_allocator, len_ref, reference, n, variants);

    variants = gva_std_allocator.allocate(gva_std_allocator.context, variants, n * sizeof(*variants), 0);
    return (PyObject*) self;
} // LCSgraph_from_variants


static PyObject*
LCSgraph_canonical(LCSgraph_Object* self)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if

    GVA_Variant* variants = gva_canonical(gva_std_allocator, self->graph);

    PyObject* const result = PyList_New(array_length(variants));
    if (result == NULL)
    {
        variants = ARRAY_DESTROY(gva_std_allocator, variants);
        return NULL;
    } // if

    for (size_t i = 0; i < array_length(variants); ++i)
    {
        PyObject* const str = PyUnicode_FromStringAndSize(variants[i].sequence.str, variants[i].sequence.len);
        if (str == NULL)
        {
            Py_DECREF(result);
            variants = ARRAY_DESTROY(gva_std_allocator, variants);
            return NULL;
        } // if
        PyObject* const args = Py_BuildValue("nnO", variants[i].start, variants[i].end, str);
        if (args == NULL)
        {
            Py_DECREF(str);
            Py_DECREF(result);
            variants = ARRAY_DESTROY(gva_std_allocator, variants);
            return NULL;
        } // if
        PyObject* const item = PyObject_CallObject((PyObject*) state->Variant_Type, args);
        Py_DECREF(args);
        if (item == NULL)
        {
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
LCSgraph_local_supremal(LCSgraph_Object* self)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if

    PyObject* const result = PyList_New(array_length(self->graph.local_supremal) - 1);
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

        PyObject* const str = PyUnicode_FromStringAndSize(variant.sequence.str, variant.sequence.len);
        if (str == NULL)
        {
            Py_DECREF(result);
            return NULL;
        } // if
        PyObject* const args = Py_BuildValue("nnO", variant.start, variant.end, str);
        if (args == NULL)
        {
            Py_DECREF(str);
            Py_DECREF(result);
            return NULL;
        } // if
        PyObject* const item = PyObject_CallObject((PyObject*) state->Variant_Type, args);
        Py_DECREF(args);
        if (item == NULL)
        {
            Py_DECREF(result);
            return NULL;
        } // if
        PyList_SET_ITEM(result, i, item);
    } // for
    return result;
} // LCSgraph_local_supremal


static PyObject*
LCSgraph_supremal(LCSgraph_Object* self)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if

    PyObject* const str = PyUnicode_FromStringAndSize(self->graph.supremal.sequence.str, self->graph.supremal.sequence.len);
    if (str == NULL)
    {
        return NULL;
    } // if
    PyObject* const args = Py_BuildValue("nnO", self->graph.supremal.start, self->graph.supremal.end, str);
    if (args == NULL)
    {
        Py_DECREF(str);
        return NULL;
    } // if
    PyObject* const result = PyObject_CallObject((PyObject*) state->Variant_Type, args);
    Py_DECREF(args);
    return result;
} // LCSgraph_supremal


static PyObject*
LCSgraph_edges(LCSgraph_Object* self)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if
    return PyObject_CallOneArg((PyObject*) state->LCSgraph_Edge_Iterator_Type, (PyObject*) self);
} // LCSgraph_edges


static PyObject*
LCSgraph_Edge_Iterator_new(PyTypeObject* subtype, PyObject* args, PyObject* kwargs)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(subtype));
    if (state == NULL)
    {
        return NULL;
    } // if

    LCSgraph_Object* graph = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|:LCSgraph_Edge_Iterator",
                                     (char*[]) {"LCSgraph", NULL},
                                     state->LCSgraph_Type, &graph))
    {
        return NULL;
    } // if

    LCSgraph_Edge_Iterator* const self = (LCSgraph_Edge_Iterator*) subtype->tp_alloc(subtype, 0);
    if (self == NULL)
    {
        return NULL;
    } // if

    Py_INCREF(graph);
    self->graph = graph;
    self->node = 0;
    self->edge = graph->graph.nodes[self->node].edges;
    return (PyObject*) self;
} // LCSgraph_Edge_Iterator_new


inline static void
LCSgraph_Edge_Iterator_dealloc(LCSgraph_Edge_Iterator* self)
{
    Py_XDECREF(self->graph);
    Py_TYPE(self)->tp_free((PyObject*) self);
} // LCSgraph_Edge_Iterator_dealloc


static PyObject*
LCSgraph_Edge_Iterator_next(LCSgraph_Edge_Iterator* self)
{
    if (self->edge == GVA_NULL)
    {
        self->node += 1;
        if (self->node >= array_length(self->graph->graph.nodes))
        {
            Py_CLEAR(self->graph);
            return NULL;
        } // if
        self->edge = self->graph->graph.nodes[self->node].edges;
    } // if

    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if

    GVA_Node const head = self->graph->graph.nodes[self->node];
    GVA_Node const tail = self->graph->graph.nodes[self->graph->graph.edges[self->edge].tail];
    GVA_Variant variant;
    gva_uint const count = gva_edges(self->graph->graph.observed.str,
        head, tail,
        self->node == self->graph->graph.source, self->graph->graph.nodes[self->graph->graph.edges[self->edge].tail].edges == GVA_NULL,
        &variant);

    self->edge = self->graph->graph.edges[self->edge].next;

    PyObject* const str = PyUnicode_FromStringAndSize(variant.sequence.str, variant.sequence.len);
    if (str == NULL)
    {
        return NULL;
    } // if
    PyObject* const args = Py_BuildValue("nnO", variant.start, variant.end, str);
    if (args == NULL)
    {
        Py_DECREF(str);
        return NULL;
    } // if
    PyObject* const var = PyObject_CallObject((PyObject*) state->Variant_Type, args);
    Py_DECREF(args);
    if (var == NULL)
    {
        Py_DECREF(str);
        return NULL;
    } // if

    return Py_BuildValue("{s: (i, i, i), s: (i, i, i), s: O, s: i}",
        "head", head.row, head.col, head.length,
        "tail", tail.row, tail.col, tail.length,
        "variant", var,
        "count", count
    );
} // LCSgraph_Edge_Iterator_next


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

    Variant_Object* const self = (Variant_Object*) subtype->tp_alloc(subtype, 0);
    if (self == NULL)
    {
        return NULL;
    } // if

    self->sequence = PyMem_Malloc(len + 1);
    if (self->sequence == NULL)
    {
        Py_DECREF(self);
        return PyErr_NoMemory();
    } // if

    self->start = start;
    self->end = end;
    self->len = len;
    memcpy(self->sequence, sequence, len);
    self->sequence[len] = '\0';
    return (PyObject*) self;
} // Variant_new


static inline void
Variant_dealloc(Variant_Object* self)
{
    PyMem_Free(self->sequence);
    Py_TYPE(self)->tp_free((PyObject*) self);
} // Variant_dealloc


static inline PyObject*
Variant_repr(Variant_Object* self)
{
    return PyUnicode_FromFormat(GVA_VARIANT_FMT, GVA_VARIANT_PRINT(((GVA_Variant) {self->start, self->end, {self->len, self->sequence}})));
} // Variant_repr


static PyObject*
Variant_richcompare(Variant_Object* self, Variant_Object* other, int op)
{
    algebra_ext_state* const state = get_algebra_ext_state(PyType_GetModule(Py_TYPE(self)));
    if (state == NULL)
    {
        return NULL;
    } // if

    if(Py_IS_TYPE(other, state->Variant_Type) == 0)
    {
        Py_RETURN_NOTIMPLEMENTED;
    } // if

    bool const equal = self == other ||
                       gva_variant_eq((GVA_Variant) {self->start, self->end, {self->len, self->sequence}},
                                      (GVA_Variant) {other->start, other->end, {other->len, other->sequence}});
    switch (op)
    {
        //case Py_LT:
        //case Py_LE:
        case Py_EQ:
            if (equal)
            {
                Py_RETURN_TRUE;
            } // if
            Py_RETURN_FALSE;
        case Py_NE:
            if (equal)
            {
                Py_RETURN_FALSE;
            } // if
            Py_RETURN_TRUE;
        //case Py_GT:
        //case Py_GE:
    } // switch
    Py_RETURN_NOTIMPLEMENTED;
} // Variant_richcompare


static PyType_Spec LCSgraph_Spec =
{
    .name = "algebra_ext.LCSgraph",
    .basicsize = sizeof(LCSgraph_Object),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
    .slots = (PyType_Slot[])
    {
        {Py_tp_doc, PyDoc_STR("LCS graph class for storing all minimal alignments.")},
        {Py_tp_new, LCSgraph_new},
        {Py_tp_dealloc, LCSgraph_dealloc},
        {Py_tp_methods, (PyMethodDef[])
            {
                {
                    "from_variants",
                    (PyCFunction) LCSgraph_from_variants,
                    METH_CLASS | METH_VARARGS | METH_KEYWORDS,
                    PyDoc_STR("from_variants(reference, variants)\n\n    Iteratively find the supremal LCS graph for an allele by repeatedly widening a range of influence."),
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
                {
                    "edges",
                    (PyCFunction) LCSgraph_edges,
                    METH_NOARGS,
                    PyDoc_STR("Edges iterator."),
                },
                {NULL, NULL, 0, NULL}  // sentinel
            },
        },
        {0, NULL}  // sentinel
    },
};


static PyType_Spec LCSgraph_Edge_Iterator_Spec =
{
    .name = "algebra_ext.LCSgraphEdgeIterator",
    .basicsize = sizeof(LCSgraph_Edge_Iterator),
    .flags = Py_TPFLAGS_DEFAULT,
    .slots = (PyType_Slot[])
    {
        {Py_tp_doc, PyDoc_STR("LCS graph edge iterator.")},
        {Py_tp_new, LCSgraph_Edge_Iterator_new},
        {Py_tp_dealloc, LCSgraph_Edge_Iterator_dealloc},
        {Py_tp_iter, PyObject_SelfIter},
        {Py_tp_iternext, LCSgraph_Edge_Iterator_next},
        {0, NULL}  // sentinel
    },
};


static PyType_Spec Variant_Spec =
{
    .name = "algebra_ext.Variant",
    .basicsize = sizeof(Variant_Object),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
    .slots = (PyType_Slot[])
    {
        {Py_tp_doc, PyDoc_STR("Variant(start, end, inserted)\n\n    Variant class for deletions/insertions.")},
        {Py_tp_new, Variant_new},
        {Py_tp_dealloc, Variant_dealloc},
        {Py_tp_repr, Variant_repr},
        {Py_tp_richcompare, Variant_richcompare},
        {Py_tp_members, (PyMemberDef[])
            {
                {"end", Py_T_PYSSIZET, offsetof(Variant_Object, end), Py_READONLY, PyDoc_STR("The end position (not included) of the deleted part.")},
                {"sequence", Py_T_STRING, offsetof(Variant_Object, sequence), Py_READONLY, PyDoc_STR("The inserted sequence.")},
                {"start", Py_T_PYSSIZET, offsetof(Variant_Object, start), Py_READONLY, PyDoc_STR("The start position (included) of the deleted part (zero-based).")},
                {NULL}  // sentinel
            }
        },
        {0, NULL}  // sentinel
    },
};


static int
algebra_ext_clear(PyObject* module)
{
    algebra_ext_state* const state = get_algebra_ext_state(module);
    Py_CLEAR(state->LCSgraph_Type);
    Py_CLEAR(state->LCSgraph_Edge_Iterator_Type);
    Py_CLEAR(state->Variant_Type);
    Py_CLEAR(state->Relations);
    return 0;
} // algebra_ext_clear


static inline void
algebra_ext_free(void* module)
{
    algebra_ext_clear(module);
} // algebra_ext_free


static int
algebra_ext_mod_exec(PyObject* module)
{
    algebra_ext_state* const state = get_algebra_ext_state(module);

    state->LCSgraph_Type = (PyTypeObject*) PyType_FromModuleAndSpec(module, &LCSgraph_Spec, NULL);
    if (state->LCSgraph_Type == NULL)
    {
        Py_DECREF(module);
        return -1;
    } // if

    state->LCSgraph_Edge_Iterator_Type = (PyTypeObject*) PyType_FromModuleAndSpec(module, &LCSgraph_Edge_Iterator_Spec, NULL);
    if (state->LCSgraph_Edge_Iterator_Type == NULL)
    {
        Py_DECREF(module);
        return -1;
    } // if

    state->Variant_Type = (PyTypeObject*) PyType_FromModuleAndSpec(module, &Variant_Spec, NULL);
    if (state->Variant_Type == NULL)
    {
        Py_DECREF(module);
        return -1;
    } // if

    if (PyModule_AddObjectRef(module, "LCSgraph", (PyObject*) state->LCSgraph_Type) < 0)
    {
        Py_DECREF(module);
        return -1;
    } // if

    if (PyModule_AddObjectRef(module, "LCSgraphEdgeIterator", (PyObject*) state->LCSgraph_Edge_Iterator_Type) < 0)
    {
        Py_DECREF(module);
        return -1;
    } // if

    if (PyModule_AddObjectRef(module, "Variant", (PyObject*) state->Variant_Type) < 0)
    {
        Py_DECREF(module);
        return -1;
    } // if

    PyObject* const dict_relations = PyDict_New();
    if (dict_relations == NULL)
    {
        Py_DECREF(module);
        return -1;
    } // if

    for (size_t i = 0; i < sizeof(GVA_RELATION_LABELS) / sizeof(*GVA_RELATION_LABELS); ++i)
    {
        PyObject* const value = PyUnicode_InternFromString(GVA_RELATION_LABELS[i]);
        if (value == NULL)
        {
            Py_DECREF(dict_relations);
            Py_DECREF(module);
            return -1;
        } // if

        if (PyDict_SetItemString(dict_relations, GVA_RELATION_LABELS[i], value) < 0)
        {
            Py_DECREF(value);
            Py_DECREF(dict_relations);
            Py_DECREF(module);
            return -1;
        } // if
        Py_DECREF(value);
    } // for

    PyObject* const module_enum = PyImport_ImportModule("enum");
    if (module_enum == NULL)
    {
        Py_DECREF(dict_relations);
        Py_DECREF(module);
        return -1;
    } // if

    state->Relations = PyObject_CallMethod(module_enum, "Enum", "sO", "Relations", dict_relations);
    Py_DECREF(dict_relations);
    Py_DECREF(module_enum);

    if (state->Relations == NULL)
    {
        Py_DECREF(module);
        return -1;
    } // if

    if (PyModule_AddObjectRef(module, "Relations", (PyObject*) state->Relations) < 0)
    {
        Py_DECREF(module);
        return -1;
    } // if

    return 0;
} // algebra_ext_mod_exec


static PyObject*
algebra_ext_compare(PyObject* self, PyObject* args, PyObject* kwargs)
{
    algebra_ext_state* const state = get_algebra_ext_state(self);
    if (state == NULL)
    {
        return NULL;
    } // if

    char const* reference = NULL;
    Py_ssize_t len_ref = 0;
    LCSgraph_Object* lhs = NULL;
    LCSgraph_Object* rhs = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#O!O!|",
                                     (char*[]) {"reference", "lhs", "rhs", NULL},
                                     &reference, &len_ref, state->LCSgraph_Type, &lhs, state->LCSgraph_Type, &rhs))
    {
        return NULL;
    } // if

    GVA_Relation const relation = gva_compare_graphs(gva_std_allocator, len_ref, reference, lhs->graph, rhs->graph);
    return PyObject_GetAttrString(state->Relations, GVA_RELATION_LABELS[relation]);
} // algebra_ext_compare


static struct PyModuleDef algebra_ext_module =
{
    PyModuleDef_HEAD_INIT,
    .m_name = "algebra_ext",
    .m_size = sizeof(algebra_ext_state),
    .m_doc = PyDoc_STR("Genetic Variant Algebra C extension."),
    .m_methods = (PyMethodDef[])
    {
        {
            "compare",
            (PyCFunction) algebra_ext_compare,
            METH_VARARGS | METH_KEYWORDS,
            PyDoc_STR("compare(reference, lhs, rhs)\n\n    Determine the relation between two LCS graphs."),
        },
        {NULL, NULL, 0, NULL}  // sentinel
    },
    .m_slots = (PyModuleDef_Slot[])
    {
        {Py_mod_exec, algebra_ext_mod_exec},
        {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
        {0, NULL}  // sentinel
    },
    .m_clear = algebra_ext_clear,
    .m_free = algebra_ext_free,
};


PyMODINIT_FUNC
PyInit_algebra_ext(void)
{
    return PyModuleDef_Init(&algebra_ext_module);
} // PyInit_algebra_ext

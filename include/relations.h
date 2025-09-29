#ifndef GVA_RELATIONS_H
#define GVA_RELATIONS_H


typedef enum
{
    GVA_DISJOINT = 0,
    GVA_OVERLAP,
    GVA_CONTAINS,
    GVA_IS_CONTAINED,
    GVA_EQUIVALENT,
} GVA_Relation;


static char const* const GVA_RELATION_LABELS[] =
{
    [GVA_DISJOINT]     = "disjoint",
    [GVA_OVERLAP]      = "overlap",
    [GVA_CONTAINS]     = "contains",
    [GVA_IS_CONTAINED] = "is_contained",
    [GVA_EQUIVALENT]   = "equivalent",
};


#endif // GVA_RELATIONS_H
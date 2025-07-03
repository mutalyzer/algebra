#ifndef GVA_FASTSTABBER_H
#define GVA_FASTSTABBER_H


typedef struct Node
{
    uint32_t start;
    uint32_t end;

    struct Node* restrict parent;
    struct Node* restrict leftsibling;
    struct Node* restrict rightchild;
} Node;


typedef struct Stack
{
    Node* restrict node;
    struct Stack* restrict next;
} Stack;


struct Bed_entry
{
    char reference[32];
    size_t count;
    size_t node_count;
    int start;
    int end;
}; // Bed_entry



Node*
stack_top(Stack const* const restrict stack);

Stack*
stack_push(Stack* const restrict stack, Node* const restrict node);

Stack*
stack_pop(Stack* const restrict stack);

Node*
create_node(uint32_t const start, uint32_t const end);

Node*
add_node(Node* const restrict last_node, Node* const restrict node);

size_t
draw_edge(Node const* const restrict node, void* const restrict arg);

size_t
draw_node(Node const* const restrict node, void* const restrict arg);

size_t
traverse(Node* const restrict node,
         size_t (*func)(Node const* const restrict, void* const restrict),
         void* const restrict arg);

size_t
write_dot(FILE* const restrict fp, Node* const restrict root);

void
destroy(Node* const restrict node);

void
read_bed(FILE* const restrict fp,
         uint32_t const max,
         Node** tree,
         Node* start_table[]);

size_t
trav(Node* const node,
     uint32_t const start,
     uint32_t const end,
     Stack* restrict* const restrict output);

size_t
stab(Node* const node,
     uint32_t const start,
     uint32_t const end,
     Stack* restrict* const restrict output);

void
query_bed(FILE* const restrict fp,
          Node* const start_table[restrict]);


#endif  // GVA_FASTSTABBER_H

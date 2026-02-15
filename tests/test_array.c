#include <assert.h>     // assert
#include <stddef.h>     // NULL, size_t
#include <stdlib.h>     // EXIT_SUCCESS

#include "../include/std_alloc.h"   // gva_std_allocator

#include "../src/array.h"   // ARRAY_*, array_*
#include "../src/array.c"   // ARRAY_*, array_*


static void
test_array_init_null(void)
{
    int* a = NULL;
    assert(array_length(a) == 0);

    a = ARRAY_DESTROY(gva_std_allocator, a);
    assert(a == NULL);

    a = array_init(gva_std_allocator, 0, sizeof(*a));
    assert(a == NULL);
} // test_array_init


static void
test_array_append_null(void)
{
    int* a = NULL;

    assert(ARRAY_APPEND(gva_std_allocator, a, 0) == 1);
    assert(array_header(a)->capacity == 1);

    assert(ARRAY_APPEND(gva_std_allocator, a, 1) == 2);
    assert(array_header(a)->capacity == 2);

    assert(ARRAY_APPEND(gva_std_allocator, a, 2) == 3);
    assert(array_header(a)->capacity == 4);

    assert(ARRAY_APPEND(gva_std_allocator, a, 3) == 4);
    assert(array_header(a)->capacity == 4);

    assert(array_length(a) == 4);
    assert(a[0] == 0);
    assert(a[1] == 1);
    assert(a[2] == 2);
    assert(a[3] == 3);

    a = ARRAY_DESTROY(gva_std_allocator, a);
} // test_array_append_null


static void
test_array_ensure(void)
{
    int* a = array_init(gva_std_allocator, 2, sizeof(*a));
    assert(array_header(a)->capacity == 2);

    a = array_ensure(gva_std_allocator, a, sizeof(*a), 10);
    assert(array_header(a)->capacity == 16);

    a = ARRAY_DESTROY(gva_std_allocator, a);
} // test_array_ensure


static void
test_array_usage(void)
{
    struct S
    {
        int x;
        int y;
    }* a = NULL;

    for (size_t i = 0; i < 1000; ++i)
    {
        ARRAY_APPEND(gva_std_allocator, a, ((struct S) {i, i}));
    } // for

    a = ARRAY_DESTROY(gva_std_allocator, a);
} // test_array_usage


int
main(int argc, char* argv[static argc])
{
    (void) argv;

    test_array_init_null();
    test_array_append_null();
    test_array_ensure();
    test_array_usage();

    return EXIT_SUCCESS;
} // main
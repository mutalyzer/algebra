#include <assert.h>     // assert
#include <stddef.h>     // NULL
#include <stdlib.h>     // EXIT_SUCCESS

#include "../include/std_alloc.h"   // gva_std_allocate


static void
test_std_allocate(void)
{
    int* a = gva_std_allocate(NULL, NULL, 0, 0);
    assert(a == NULL);

    a = gva_std_allocate(NULL, a, 0, 1000);
    assert(a != NULL);

    a = gva_std_allocate(NULL, a, 1000, 0);
} // test_std_allocate


int
main(int argc, char* argv[static argc])
{
    (void) argv;

    test_std_allocate();

    return EXIT_SUCCESS;
} // main
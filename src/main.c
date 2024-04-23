#include <stddef.h>     // NULL, size_t
#include <stdio.h>      // printf
#include <stdlib.h>     // EXIT_*

#include "../include/alloc.h"       // VA_Allocator
#include "../include/array.h"       // va_array_*
#include "../include/string.h"      // VA_String, va_string_*
#include "../include/variant.h"     // VA_Variant
#include "../include/std_alloc.h"   // std_allocator


static VA_String
patch(VA_String const reference, size_t const len, VA_Variant variants[len], VA_String const external)
{
    VA_String observed = {0, NULL};

    size_t start = 0;
    for (size_t i = 0; i < len; ++i)
    {
        observed = va_string_concat(&std_allocator, observed, (VA_String) {variants[i].start - start, reference.data + start});
        observed = va_string_concat(&std_allocator, observed, (VA_String) {variants[i].obs_end - variants[i].obs_start, external.data + variants[i].obs_start});
        start = variants[i].end;
    } // for

    if (start < reference.length)
    {
        observed = va_string_concat(&std_allocator, observed, (VA_String) {reference.length - start, reference.data + start});
    } // if
    return observed;
} // patch


int
main(int argc, char* argv[argc + 1])
{
    (void) argv;

    VA_String const observed = patch((VA_String) {10, "ATTACCATTA"}, 1, (VA_Variant[]) {{3, 4, 0, 1}}, (VA_String) {1, "G"});

    printf("%.*s\n", (int) observed.length, observed.data);

    std_allocator.alloc(NULL, observed.data, 0, 0);

    return EXIT_SUCCESS;
} // main

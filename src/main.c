#include <stdio.h>      // stderr, fprintf, printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // strlen

#include "../include/parser.h"  // va_parse_*


int
main(int argc, char* argv[argc + 1])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage %s expression\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    printf("%zu\n", va_parse_spdi(strlen(argv[1]), argv[1]));

    return EXIT_SUCCESS;
} // main

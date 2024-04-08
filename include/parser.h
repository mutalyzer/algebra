#ifndef PARSER_H
#define PARSER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t

#include "variant.h"    // VA_Variant


size_t
va_parse_hgvs(size_t const len, char const expression[static len],
              size_t const len_var, VA_Variant variants[static len_var]);


size_t
va_parse_spdi(size_t const len, char const expression[static len],
              VA_Variant variants[static 1]);


#ifdef __cplusplus
} // extern "C"
#endif

#endif

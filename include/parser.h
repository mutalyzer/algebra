#ifndef PARSER_H
#define PARSER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t


size_t
va_parse_hgvs(size_t const len, char const expression[static len]);


size_t
va_parse_spdi(size_t const len, char const expression[static len]);


#ifdef __cplusplus
} // extern "C"
#endif

#endif

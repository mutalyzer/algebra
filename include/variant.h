#ifndef VARIANT_H
#define VARIANT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t


typedef struct VA_Variant
{
    uint32_t start;
    uint32_t end;
    uint32_t obs_start;
    uint32_t obs_end;
} VA_Variant;


size_t
va_variant_len(VA_Variant const variant);


#ifdef __cplusplus
} // extern "C"
#endif

#endif

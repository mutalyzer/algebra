#ifndef VARIANT_H
#define VARIANT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <stdint.h>     // uint32_t


typedef struct
{
    uint32_t start;
    uint32_t end;
    uint32_t obs_start;
    uint32_t obs_end;
} VA_Variant;


bool
va_variant_eq(VA_Variant const lhs, VA_Variant const rhs);


size_t
va_variant_len(VA_Variant const variant);


size_t
va_patch(size_t const len_ref, char const reference[static len_ref],
         size_t const len_obs, char const observed[static len_obs],
         size_t const len_var, VA_Variant const variants[len_var],
         size_t const len, char string[static len]);


#ifdef __cplusplus
} // extern "C"
#endif

#endif

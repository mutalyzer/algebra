#include <stddef.h>     // NULL, size_t

#include "wasm.h"   // WASM_EXPORT
#include "../../include/edit.h"    // varalg_edit


WASM_EXPORT("edit")
size_t
edit(size_t const lhs, size_t const rhs)
{
    console_log(13, "Hello, World!");
    console_log(42, NULL);

    return varalg_edit(lhs, rhs);
} // edit

CC       = clang

CFLAGS   = -std=c99 -Wall -Wextra -Wpedantic -O3 \
           --target=wasm32 -nostdlib \
           -Wl,--no-entry, -Wl,--allow-undefined, -Wl,--import-memory, \
           -Wl,--export=__heap_base

algebra.wasm: wrapper.c wasm.h Makefile
	$(CC) $(CFLAGS) -o $@ ../src/edit.c $<

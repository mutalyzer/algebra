export function wasmInstance(path, initial_pages=2) {
    const utf8decoder = new TextDecoder();
    let heap_blocks;

    const imports = { env : {
        memory: new WebAssembly.Memory({
            initial: initial_pages,
        }),

        console_log: function(len, msg) {
            if (len === 0) {
                console.log(msg);
                return;
            } // if
            console.log(utf8decoder.decode(new Uint8Array(imports.env.memory.buffer, msg, len)));
        },

        free: function(address) {
            console.groupCollapsed(`free(0x${address.toString(16)})`);
            if (address === null || address === 0) {
                console.log("freeing null");
                console.groupEnd();
                return;
            } // if
            const len = heap_blocks.length;
            for (let idx = 0; idx < len; ++idx) {
                let heap_block = heap_blocks[idx];
                if (heap_block.address === address) {
                    console.log(`found heap block ${idx}`);
                    if (heap_block.free) {
                        console.warn(`double free at 0x${address.toString(16)}`);
                        console.groupEnd();
                        return;
                    } // if
                    heap_block.free = true;
                    if (idx > 0 && heap_blocks[idx - 1].free) {
                        console.log(`coalesce with heap block ${idx - 1}`);
                        heap_blocks[idx - 1].size += heap_block.size;
                        heap_blocks.splice(idx, 1);
                        idx -= 1;
                        heap_block = heap_blocks[idx];
                    } // if
                    if (idx < len - 1 && heap_blocks[idx + 1].free) {
                        console.log(`coalesce with heap block ${idx + 1}`);
                        heap_blocks[idx + 1].address = heap_block.address;
                        heap_blocks[idx + 1].size += heap_block.size;
                        heap_blocks.splice(idx, 1);
                    } // if
                    console.log("heap_blocks:", JSON.parse(JSON.stringify(heap_blocks)));
                    console.groupEnd();
                    return;
                } // if
            } // for
            console.warn(`illegal free at 0x${address.toString(16)}`);
            console.groupEnd();
        },

        malloc: function(size) {
            console.groupCollapsed(`malloc(${size})`);
            if (size <= 0) {
                console.warn("invalid size");
                console.groupEnd();
                return null;
            } // if
            const word_size = imports.env.word_size;
            if (size % word_size > 0) {
                size += word_size - size % word_size;
            } // if
            console.log(`for ${word_size}-byte aligned ${size} bytes`);
            const len = heap_blocks.length;
            for (let idx = 0; idx < len; ++idx) {
                const heap_block = heap_blocks[idx];
                const diff = heap_block.size - size;
                if (heap_block.free && diff >= 0) {
                    console.log(`found suitable heap block ${idx}`);
                    heap_block.free = false;
                    const fragment_size = 0;
                    if (diff > fragment_size) {
                        console.log(`splitting ${size} bytes from heap block ${idx}`);
                        heap_block.size = size;
                        heap_blocks.splice(idx + 1, 0, {
                            address: heap_block.address + size,
                            size: diff,
                            free: true,
                        });
                    } // if
                    console.log("heap blocks:", JSON.parse(JSON.stringify(heap_blocks)));
                    console.log(`address: 0x${heap_block.address.toString(16)}`);
                    console.groupEnd();
                    return heap_block.address;
                } // if
            } // for
            const page_size = 1 << 16;
            const pages = Math.floor(size / page_size + 1);
            console.log(`grow memory with ${pages} 16KiB pages`);
            const memory = imports.env.memory;
            try {
                memory.grow(pages);
            } // try
            catch(error) {
                console.warn("out of memory");
                console.groupEnd();
                return null;
            } // catch
            const buffer = memory.buffer;
            const buffer_size = buffer.byteLength;
            memory.words = new Uint32Array(buffer, 0, buffer_size / env.word_size);
            memory.bytes = new Uint8Array(buffer, 0, buffer_size);
            console.log(`total memory: ${memory.buffer.byteLength / (1024 * 1024)} MiB`);
            const heap_block = heap_blocks[heap_blocks.length - 1];
            heap_blocks.push({
                address: heap_block.address + heap_block.size,
                size: pages * page_size,
                free: true,
            });
            console.log("heap blocks:", JSON.parse(JSON.stringify(heap_blocks)));
            console.groupEnd();
            return imports.env.malloc(size);
        },
    }}; // imports

    return WebAssembly.instantiateStreaming(fetch(path), imports).then(obj => {
        const env = imports.env;
        const memory = env.memory;
        const buffer = memory.buffer;
        const buffer_size = buffer.byteLength;
        const exports = obj.instance.exports;
        const word_size = env.word_size = exports.__word_size();
        const heap_base = exports.__heap_base;
        console.assert(word_size === 4);  // For now we assume wasm32
        console.assert(heap_base % word_size === 0);
        heap_blocks = [{
            address: heap_base.value,
            size: buffer_size - heap_base,
            free: true,
        }];
        memory.words = new Uint32Array(buffer, 0, buffer_size / word_size);
        memory.bytes = new Uint8Array(buffer, 0, buffer_size);
        return {
            env: env,
            exports: exports,
        };
    });
} // wasmInstance

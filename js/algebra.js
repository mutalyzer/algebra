import { wasmInstance } from "./wasm.js";


export const algebra = {};


await wasmInstance("./algebra.wasm").then(wasm => {
    algebra.edit = function(lhs, rhs) {
        return wasm.exports.edit(lhs, rhs ?? 42);
    },
});






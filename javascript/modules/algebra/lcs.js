import { wasmInstance } from "../wasm.js";


export let edit;


await wasmInstance("./algebra.wasm").then(wasm => {
    edit = (lhs, rhs) => wasm.exports.edit(lhs, rhs ?? 42);
});

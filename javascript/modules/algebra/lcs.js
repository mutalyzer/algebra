import { wasmInstance } from "../wasm.js";


let edit;


await wasmInstance("./algebra.wasm").then(wasm => {
    edit = (lhs, rhs) => wasm.exports.edit(lhs, rhs ?? 42);
});


export { edit };

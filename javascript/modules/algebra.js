import { wasmInstance } from "./wasm.js";
import * as lcs from "./algebra/lcs.js";


export const algebra = await wasmInstance("./algebra.wasm").then(wasm => {
    return {
        lcs,
        hello: (lhs, rhs) => lhs * (rhs ?? 42),
    };
});

from setuptools import Extension, setup


setup(
    ext_modules = [
        Extension("algebra_ext",
                  extra_compile_args=[
                    "-Wextra",
                    "-Wno-cast-function-type",
                    "-std=c99",
                  ],
                  sources=[
                    "algebra_extmodule.c",
                    "../src/align.c",
                    "../src/array.c",
                    "../src/bitset.c",
                    "../src/compare.c",
                    "../src/edit.c",
                    "../src/extractor.c",
                    "../src/lcs_graph.c",
                    "../src/string.c",
                    "../src/variant.c",
                  ]),
    ],
)
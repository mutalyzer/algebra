from setuptools import Extension, setup


setup(
    ext_modules = [
        Extension("algebra_ext",
                  extra_compile_args=[
                    "-Wextra",
                    "-Wpedantic",
                    "-Wno-cast-function-type",
                    "-std=c99",
                    "-O3",
                    "-fanalyzer",
                  ],
                  sources=[
                    "ext/algebra_ext.c",
                    "ext/lcsgraph.c",
                    "ext/variant.c",
                    "../src/align.c",
                    "../src/array.c",
                    "../src/extractor.c",
                    "../src/lcs_graph.c",
                    "../src/string.c",
                    "../src/variant.c",
                  ]),
    ],
)

from setuptools import Extension, setup


setup(
    ext_modules = [
        Extension("algebra_ext",
                  extra_compile_args=[
                    "-Wextra",
                    "-Wpedantic",
                    "-std=c99",
                  ],
                  sources=[
                    "ext/algebra_ext.c",
                    "ext/lcsgraph.c",
                    "ext/variant.c",
                    "../src/align.c",
                    "../src/array.c",
                    "../src/extract.c",
                    "../src/lcs_graph.c",
                    "../src/string.c",
                    "../src/variant.c",
                  ]),
    ],
)

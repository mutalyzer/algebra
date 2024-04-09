from setuptools import Extension, setup


setup(
    ext_modules = [
        Extension("algebra_ext",
                  extra_compile_args=[
                    "-Wextra",
                    "-Wpedantic",
                    "-std=c11",
                  ],
                  sources=[
                    "ext/algebra_ext.c",
                    "src/parser.c",
                    "src/variant.c",
                  ]),
    ],
)

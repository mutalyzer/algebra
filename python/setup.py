from setuptools import Extension, setup


compile_args = [
    "-Wextra",
    "-Wpedantic",
    "-std=c99",
]


setup(
    ext_modules = [
        Extension("algebra_ext", extra_compile_args=compile_args,
                  sources=[
                    "ext/algebra_ext.c"
                  ]),
        Extension("lcs_ext", extra_compile_args=compile_args,
                  sources=[
                    "ext/lcs_ext.c",
                    "src/edit.c"
                  ]),
    ],
)

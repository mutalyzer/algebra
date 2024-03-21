from setuptools import Extension, setup


setup(
    ext_modules = [
        Extension("algebra_ext",
                  sources=[
                    "ext/algebra_ext.c"
                  ],
                 ),
        Extension("lcs_ext",
                  sources=[
                    "ext/lcs_ext.c", "src/edit.c"
                  ],
                 ),
    ],
)

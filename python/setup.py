from setuptools import setup, Extension


setup(
    name = "algebra",
    ext_modules = [
        Extension("algebra",
            sources = [
                "src/edit.c",
                "ext/wrapper.c",
            ],
            extra_compile_args = [
                "-Wextra",
                "-Wpedantic",
                "-std=c11",
            ],
        )
    ],
)

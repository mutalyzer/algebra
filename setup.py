from setuptools import setup, find_packages

setup(
    name="mutalyzer-algebra",
    version="0.1.0",
    description="A library and CLI for A Boolean Algebra for Genetic Variants.",
    url="https://github.com/mutalyzer/algebra",
    author="Mark Santcroos",
    author_email="m.a.santcroos@lumc.nl",
    license="MIT",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "algebraist=algebra.cli:main",
        ],
    },
    install_requires="normalizer @ git+ssh://git@github.com/mutalyzer/normalizer.git",
)

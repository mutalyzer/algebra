from setuptools import find_packages, setup


setup(
    name="mutalyzer-algebra",
    version="1.0.2",
    license="MIT",
    author='Jonathan Vis, Mark Santcroos',
    author_email='j.k.vis@lumc.nl, m.a.santcroos@lumc.nl',
    packages=find_packages(),
    extras_require={
        "dev": [
            "coverage",
            "flake8",
            "pytest",
        ]
    },
    entry_points={
        "console_scripts": [
            "algebra=algebra.__main__:main",
        ],
    },
)

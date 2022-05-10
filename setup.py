from setuptools import find_packages, setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="mutalyzer-algebra",
    version="1.0.3",
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
    long_description=long_description,
    long_description_content_type='text/markdown'
)

from setuptools import find_packages, setup


setup(
    name="algebra",
    version="2.0.1",
    license="MIT",
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

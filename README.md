mutalyzer-algebra
=================
A Boolean Algebra for Genetic Variants  

A set of Boolean relations: equivalence; containment, i.e., either a
variant is fully contained in another or a variant fully contains another;
overlap, i.e., two variants have (at least) one common element; and
disjoint, i.e., no common elements that allows for a comprehensive
classification of the relation for every pair of variants by taking all
minimal Longest Common Subsequence (LCS) alignments into account.

[Jonathan K. Vis, Mark A. Santcroos, Walter A. Kosters and Jeroen F.J. Laros.
"A Boolean Algebra for Genetic Variants." (2021)](https://arxiv.org/abs/2112.14494)

Installation
------------

Use pip to install from the Python Package Index (PyPI).

```bash
python -m pip install mutalyzer-algebra
```

Or directly from GitHub for development (after cloning in an active
virtual environment).

```bash
python -m pip install --upgrade --editable .[dev]
```

Testing
-------

Run the tests.

```bash
python -m coverage run -m pytest
```

Usage
-----

Use the command-line interface.

```bash
algebra --reference "AAAAA" compare --lhs-hgvs "1_2insTA" --rhs-hgvs "2_3insT"
```

Or as a Python package.

```python
from algebra import compare
from algebra.variants import Parser

reference = "AAAAA"
lhs = Parser("1_2insTA").hgvs()
rhs = Parser("2_3insT").hgvs()

# returns: Relation.DISJOINT
compare(reference, lhs, rhs)
```

See Also
--------

A web interface with integration with [Mutalyzer](https://github.com/mutalyzer): [Mutalyzer Algebra](https://v3.mutalyzer.nl/algebra)

[Mutalyzer Algebra on PyPI](https://pypi.org/project/mutalyzer-algebra/)

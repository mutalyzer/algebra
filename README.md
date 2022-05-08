algebra
=======
A Boolean Algebra for Genetic Variants  

A set of Boolean relations that allows for a comprehensive classification
of the relation for every pair of variants by taking all minimal
Longest Common Subsequence (LCS) alignments into account.

[Jonathan K. Vis, Mark Santcroos, Walter A. Kosters and Jeroen F.J. Laros. "A Boolean Algebra for Genetic Variants." (2021)](https://arxiv.org/abs/2112.14494)

Installation
------------

Use pip to install from the Python Package Index (PyPI).

```bash
python -m pip install mutalyzer-algebra
```

Or directly from github for development (after cloning in an active
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
from algebra.relations import compare
from algebra.variants.parser import Parser
from algebra.variants.variant import patch

reference = "AAAAA"
lhs = patch(reference, Parser("1_2insTA").hgvs())
rhs = patch(reference, Parser("2_3insT").hgvs())

# returns: Relation.DISJOINT
compare(reference, lhs, rhs)
```

See Also
--------

A web interface with integration with [Mutalyzer](https://github.com/mutalyzer).

[Mutalyzer Algebra](https://v3.mutalyzer.nl/algebra)

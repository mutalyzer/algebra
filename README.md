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
"A Boolean Algebra for Genetic Variants." In: *Bioinformatics* (2023).](https://doi.org/10.1093/bioinformatics/btad001)

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
from algebra.variants import parse_hgvs


reference = "AAAAA"
lhs = parse_hgvs("1_2insTA")
rhs = parse_hgvs("2_3insT")

# returns: Relation.DISJOINT
compare(reference, lhs, rhs)


reference = "CATATATC"
lhs = parse_hgvs("2_7AT[4]")  # observed: CATATATATC
rhs = parse_hgvs("5_6insT")   # observed: CATATTATC

# returns: Relation.CONTAINS
compare(reference, lhs, rhs)
```

Extracting variants from sequences.

```python
from algebra.extractor import extract_sequence, to_hgvs


reference = "CATATATC"
observed = "CATATATATC"

canonical, _ = extract_sequence(reference, observed)
# returns: 2_7AT[4]
to_hgvs(canonical, reference)
```

Variant normalization.

```python
from algebra.extractor import extract, to_hgvs
from algebra.variants import parse_hgvs


reference = "CATATATC"
variant = parse_hgvs("6_7dupAT")

canonical, _ = extract(reference, variant)
# returns: 2_7AT[4]
to_hgvs(canonical, reference)
```

See Also
--------

A web interface with integration with [Mutalyzer](https://github.com/mutalyzer): [Mutalyzer Algebra](https://mutalyzer.nl/algebra)

[Mutalyzer Algebra on PyPI](https://pypi.org/project/mutalyzer-algebra/)

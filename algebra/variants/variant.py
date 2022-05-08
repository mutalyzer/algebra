"""Variant class and related functions.

Variants are represented as deletion/insertions. The deleted part is
given as a zero-based half-open interval with regard to some reference
sequence (not stored together with a variant). The insertion is a
sequence (string) of inserted symbols.

See Also
--------
algebra.variants.parser : To construct variants from HGVS and SPDI
representations.
"""


from itertools import combinations


class Variant:
    """Variant class for deletion/insertions."""

    def __init__(self, start, end, sequence=""):
        """Create a variant.

        Parameters
        ----------
        start : int
            The start position (included) of the deleted part (zero-based).
        end : int
            The end position (not included) of the deleted part.
        sequence : str, optional
            The inserted sequence.

        Raises
        ------
        TypeError
            If the parameters are not of the correct type.
        ValueError
            If the interval [`start`, `end`) is invalid, e.g.,
            `start` > `end`.
        """

        if not isinstance(start, int):
            raise TypeError("start must be an integer")
        if not isinstance(end, int):
            raise TypeError("end must be an integer")
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")

        if start < 0:
            raise ValueError("start must be greater or equal to 0")
        if start > end:
            raise ValueError("start must not be after end")

        self.start = start
        self.end = end
        self.sequence = sequence

    def __bool__(self):
        return self.end - self.start > 0 or len(self.sequence) > 0

    def __eq__(self, other):
        return (self.start == other.start and self.end == other.end and
                self.sequence == other.sequence)

    def __hash__(self):
        return hash((self.start, self.end, self.sequence))

    def __len__(self):
        return self.end - self.start + len(self.sequence)

    def __lt__(self, other):
        if ((self.start < other.start and self.end > other.start) or
                (other.start < self.start and other.end > self.start) or
                (self.start == other.start and self.end == other.end and
                    len(self.sequence) > 0)):
            # deletions and insertions are handled asymmetrically
            raise ValueError("variants overlap")

        return self.start < other.start or self.end < other.end

    def __repr__(self):
        return f"<{self.start},{self.end}/{self.sequence}>"

    def atomics(self):
        """Generate all atomic representations.

        Generates all alternative representations for a deletion/insertion
        by using separate deletions and insertions.

        Yields
        ------
        list
            A list of variants representing the original variant.
        """

        for combo in combinations(range(len(self)), len(self.sequence)):
            variants = []
            c = 0
            pos = self.start
            variant = Variant(pos, pos)
            for i in range(len(self.sequence)):
                if combo[i] > c:
                    if variant:
                        variants.append(variant)
                    for j in range(pos, pos + combo[i] - c):
                        variants.append(Variant(j, j + 1))
                    pos += combo[i] - c
                    c = combo[i]
                    variant = Variant(pos, pos, self.sequence[i])
                else:
                    variant.sequence += self.sequence[i]
                c += 1

            if variant:
                variants.append(variant)
            for i in range(pos, self.end):
                variants.append(Variant(i, i + 1))

            yield variants

    def is_disjoint(self, other):
        """Check if two variants are disjoint, i.e., no common deletion or
        insertion."""
        if (self.start < other.end and other.start < self.end and
                self.start < self.end and other.start < other.end):
            return False
        return (self.start > other.end or other.start > self.end or
                set(self.sequence).isdisjoint(set(other.sequence)))

    def to_hgvs(self, reference=None, only_substitutions=True):
        """The variant representation in HGVS [1]_.

        Parameters
        ----------
        reference : str or None, optional
            The reference sequence for the variant. Substitutions (3A>T)
            include the deleted nucleotide in HGVS. If omitted deleted
            symbols will not be filled in.
        only_substitutions : bool, optional
            By default only includes deleted nucleotides for HGVS
            substitutions. If set to `False` all deleted symbols will be
            filled in.

        References
        ----------
        [1] https://varnomen.hgvs.org/.
        """
        if self.end - self.start == 0:
            if len(self.sequence) == 0:
                return "="
            return f"{self.start}_{self.start + 1}ins{self.sequence}"

        deleted = ""
        substitution = ""
        if reference is not None:
            if not only_substitutions:
                deleted = reference[self.start:self.end]
            substitution = reference[self.start:self.end]

        if self.end - self.start == 1:
            if len(self.sequence) == 0:
                return f"{self.start + 1}del{deleted}"
            if len(self.sequence) == 1:
                return f"{self.start + 1}{substitution}>{self.sequence}"
            return f"{self.start + 1}del{deleted}ins{self.sequence}"

        if len(self.sequence) == 0:
            return f"{self.start + 1}_{self.end}del{deleted}"

        return f"{self.start + 1}_{self.end}del{deleted}ins{self.sequence}"

    def to_spdi(self, reference):
        """The variant representation in SPDI [1]_.

        References
        ----------
        [1] Holmes JB, Moyer E, Phan L, Maglott D, Kattman B.
        SPDI: data model for variants and applications at NCBI.
        Bioinformatics. 2020 Mar 1;36(6):1902-1907.
        """
        return (f"{reference}:{self.start}:{self.end - self.start}:"
                f"{self.sequence}")


def patch(reference, variants):
    """Apply a list of variants to a reference sequence to obtain an
    observed sequence.

    Parameters
    ----------
    reference : str
        The reference sequence to which the `variants` are applied.
    variants : list
        A list of variants.

    Raises
    ------
    ValueError
        If variants are overlapping.

    Returns
    -------
    str
        The observed sequence.
    """

    def slices(reference, variants):
        start = 0

        for variant in sorted(variants):
            yield reference[start:variant.start] + variant.sequence
            start = variant.end

        if len(reference[start:]) > 0:
            yield reference[start:]

    return "".join(slices(reference, variants))


def to_hgvs(variants, reference=None, only_substitutions=True, sequence_prefix=True, sort=True):
    """An allele representation of a list of variants in HGVS [1]_.

    Parameters
    ----------
    variants : list
        A list of variants.
    reference : str or None, optional
        The reference sequence.

    Other Parameters
    ----------------
    only_substitutions : bool, optional
        Deleted symbols are filled in only for HGVS substitutions by default.
    sequence_prefix : bool, optional
        Prefix the allele description with the reference sequence
        (not RefSeq ID).
    sort : bool, optional
        The `variants` are sorted by default. Disable if already sorted.

    Returns
    -------
    str
        The allele representation in HGVS.

    References
    ----------
    [1] https://varnomen.hgvs.org/recommendations/DNA/variant/alleles/.
    """

    prefix = ""
    if reference is not None and sequence_prefix:
        prefix = f"{reference}:g."

    if len(variants) == 0:
        return f"{prefix}="

    if len(variants) == 1:
        return f"{prefix}{variants[0].to_hgvs(reference, only_substitutions)}"

    return f"{prefix}[{';'.join([variant.to_hgvs(reference, only_substitutions) for variant in (sorted(variants) if sort else variants)])}]"

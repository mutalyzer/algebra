"""A hand-crafted recursive descent parser for (simple) genomic variants.

Parses (simple) genomic variants in HGVS [1]_ or SPDI [2]_ format.

See Also
--------
algebra.variants.variant : The variant class.

References
----------
[1] https://varnomen.hgvs.org/.
[2] J.B. Holmes, E. Moyer, L. Phan, D. Maglott and B. Kattman. "SPDI:
data model for variants and applications at NCBI".
In: Bioinformatics 36.6 (2019), pp. 1902-1907.
"""


from .variant import DNA_NUCLEOTIDES, Variant, reverse_complement


class Parser:
    """Parser class."""

    def __init__(self, expression):
        """Create a parser for an `expression`.

        Raises
        ------
        TypeError
            If `expression` is not a string.
        """

        if not isinstance(expression, str):
            raise TypeError("expression must be a string")

        self.pos = 0
        self.expression = expression

    def _match(self, word):
        if not self.expression[self.pos:self.pos + len(word)] == word:
            return False

        self.pos += len(word)
        return True

    def _match_digit(self):
        if self.pos >= len(self.expression):
            raise ValueError("unexpected end of expression")
        if not self.expression[self.pos].isdigit():
            raise ValueError(f"expected digit at {self.pos}")

        self.pos += 1
        return int(self.expression[self.pos - 1])

    def _match_nucleotide(self):
        if self.pos >= len(self.expression):
            raise ValueError("unexpected end of expression")
        if not self.expression[self.pos] in DNA_NUCLEOTIDES:
            raise ValueError(f"expected nucleotide at {self.pos}")

        self.pos += 1
        return self.expression[self.pos - 1]

    def _match_number(self):
        number = self._match_digit()
        while True:
            try:
                number = number * 10 + self._match_digit()
            except ValueError:
                return number

    def _match_sequence(self):
        start = self.pos
        self._match_nucleotide()
        while True:
            try:
                self._match_nucleotide()
            except ValueError:
                return self.expression[start:self.pos]

    def _match_variant(self, reference):
        start = self._match_number() - 1
        end = None
        if self._match("_"):
            end = self._match_number()
            if start >= end:
                raise ValueError(f"invalid range at {self.pos}")

        if self._match("dup"):
            if end is None:
                end = start + 1
            duplicated = None
            try:
                duplicated = self._match_sequence()
            except ValueError:
                pass

            if reference is None:
                raise NotImplementedError("duplications without the reference context are not supported")
            if duplicated is not None:
                if len(duplicated) != end - start:
                    raise ValueError(f"inconsistent duplicated length at {self.pos}")

            if start < 0 or len(reference) < end:
                raise ValueError("invalid range in reference")
            if duplicated is not None and reference[start:end] != duplicated:
                raise ValueError(f"'{duplicated}' not found in reference at {start}")

            return Variant(end, end, reference[start:end])

        if self._match("inv"):
            if end is None:
                end = start + 1
            if reference is None:
                raise NotImplementedError("inversions without the reference context are not supported")
            if len(reference) < end:
                raise ValueError("invalid range in reference")

            return Variant(start, end, reverse_complement(reference[start:end]))

        if self._match("del"):
            if end is None:
                end = start + 1
            deleted = None
            try:
                deleted = self._match_sequence()
            except ValueError:
                pass

            if deleted is not None:
                if len(deleted) != end - start:
                    raise ValueError(f"inconsistent deleted length at {self.pos}")
                if reference is not None and reference[start:end] != deleted:
                    raise ValueError(f"'{deleted}' not found in reference at {start}")

            if self._match("ins"):
                return Variant(start, end, self._match_sequence())
            return Variant(start, end)

        if self._match("ins"):
            if end is None or end - start != 2:
                raise ValueError(f"invalid inserted range at {self.pos}")
            return Variant(start + 1, start + 1, self._match_sequence())

        if end is not None:
            # new HGVS-style repeat notation
            repeat_unit = self._match_sequence()
            if not self._match("["):
                raise ValueError(f"expected '[' at {self.pos}")
            repeat_number = self._match_number()
            if not self._match("]"):
                raise ValueError(f"expected ']' at {self.pos}")
            return Variant(start, end, repeat_number * repeat_unit)

        sequence = None
        try:
            sequence = self._match_sequence()
        except ValueError:
            pass
        if (reference is not None and sequence is not None and
                reference[start:start + len(sequence)] != sequence):
            raise ValueError(f"'{sequence}' not found in reference at {start}")

        if self._match("="):
            return Variant(start, start)
        if sequence is not None and self._match("["):
            # dbSNP HGVS repeat notation
            if reference is None:
                raise NotImplementedError("dbSNP repeats without reference context are not supported")
            repeat_number = self._match_number()
            if not self._match("]"):
                raise ValueError(f"expected ']' at {self.pos}")
            found = 0
            while reference[start + found * len(sequence):start + (found + 1) * len(sequence)] == sequence:
                found += 1

            return Variant(start, start + found * len(sequence), sequence * repeat_number)

        if not self._match(">"):
            raise ValueError(f"expected '>' at {self.pos}")
        return Variant(start, start + 1, self._match_sequence())

    def hgvs(self, reference=None, skip_reference=False):
        """Parse an expression as HGVS.

        Only simple (deletions, insertions, substitutions,
        deletion/insertions and repeats) are supported in a genomic
        coordinate system.

        Parameters
        ----------
        reference : str or None, optional
            The reference sequence useful for handling inversions and/or
            duplications.
        skip_reference : bool, optional
            Skips a reference sequence identifier and coordinate system.

        Raises
        ------
        ValueError
            If a syntax (or a simple sematic) error occurs.
        NotImplementedError
            If a syntax construct is not supported.

        Returns
        -------
        list
            A sorted list of variants (allele).
        """

        if skip_reference:
            while self.pos < len(self.expression) and not self._match(":"):
                self.pos += 1
            if not self._match("g."):
                raise ValueError(f"expected 'g.' at {self.pos}")

        if self._match("="):
            if not self.pos == len(self.expression):
                raise ValueError(f"expected end of expression at {self.pos}")
            return []

        variants = []
        if self._match("["):
            variant = self._match_variant(reference)
            if variant:
                variants.append(variant)
            while self._match(";"):
                variant = self._match_variant(reference)
                if variant:
                    variants.append(variant)
            if not self._match("]"):
                raise ValueError(f"expected ']' at {self.pos}")
        else:
            variant = self._match_variant(reference)
            if variant:
                variants.append(variant)

        if not self.pos == len(self.expression):
            raise ValueError(f"expected end of expression at {self.pos}")
        self.pos = 0
        return sorted(variants)

    def spdi(self):
        """Parse an expression as SPDI.

        Raises
        ------
        ValueError
            If a syntax error occurs.

        Returns
        -------
        `Variant`
            The variant.
        """

        # ignore reference sequence identifier
        while self.pos < len(self.expression) and not self._match(":"):
            self.pos += 1

        position = self._match_number()
        if not self._match(":"):
            raise ValueError(f"expected ':' at {self.pos}")

        length = 0
        try:
            length = self._match_number()
        except ValueError:
            try:
                length = len(self._match_sequence())
            except ValueError:
                pass
        if not self._match(":"):
            raise ValueError(f"expected ':' at {self.pos}")

        inserted = ""
        try:
            inserted = self._match_sequence()
        except ValueError:
            pass
        if not self.pos == len(self.expression):
            raise ValueError(f"expected end of expression at {self.pos}")
        self.pos = 0
        return [Variant(position, position + length, inserted)]

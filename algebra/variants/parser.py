"""A hand-crafted recursive descent parser for (simple) genomic variants.

Parses (simple) genomic variants in HGVS [1]_ or SPDI [2]_ format.

See Also
--------
algebra.variants.variant : The `Variant` class.

References
----------
[1] https://varnomen.hgvs.org/.
[2] J.B. Holmes, E. Moyer, L. Phan, D. Maglott and B. Kattman. "SPDI:
data model for variants and applications at NCBI".
In: Bioinformatics 36.6 (2019), pp. 1902-1907.
"""


from .variant import DNA_NUCLEOTIDES, Variant, reverse_complement


def parse_hgvs(expression, reference=None):
    """Parse an expression as HGVS.

    Only simple (deletions, insertions, substitutions,
    deletion/insertions and repeats) are supported in a genomic
    coordinate system.

    Other Parameters
    ----------------
    reference : str or None, optional
        The reference sequence useful for handling inversions,
        duplications and/or repeats.

    Raises
    ------
    TypeError
        If `expression` is not a string.
    ValueError
        If a syntax (or a simple semantic) error occurs.
    NotImplementedError
        If a syntax construct is not supported.

    Returns
    -------
    list
        A sorted list of variants (allele).
    """

    def match(word):
        nonlocal pos
        if pos > len(expression) - len(word):
            raise ValueError("unexpected end of expression")
        if expression[pos:pos + len(word)] != word:
            raise ValueError(f"expected '{word}' at {pos + 1}")
        pos += len(word)
        return word

    def match_plus(predicate, label=None):
        nonlocal pos
        if pos >= len(expression):
            raise ValueError("unexpected end of expression")
        if not predicate(expression[pos]):
            raise ValueError(f"expected {label if label is not None else predicate} at {pos + 1}")
        start = pos
        pos += 1
        while pos < len(expression) and predicate(expression[pos]):
            pos += 1
        return expression[start:pos]

    def match_optional(word):
        try:
            return match(word) == word
        except ValueError:
            return False

    def match_number():
        return int(match_plus(lambda ch: ch.isdigit(), "digit"))

    def match_sequence():
        return match_plus(lambda ch: ch in DNA_NUCLEOTIDES, "nucleotide")

    def match_location():
        start = match_number()
        end = match_number() if match_optional("_") else start
        return start - 1, end

    def match_variant(reference):
        start, end = match_location()
        ctx_pos = pos

        if match_optional("dup"):
            try:
                sequence = match_sequence()
            except ValueError:
                if reference is None:
                    raise NotImplementedError(f"duplication without reference context at {ctx_pos + 1}") from None
                if end > len(reference):
                    raise ValueError("invalid range in reference") from None
                sequence = reference[start:end]
            else:
                if len(sequence) != end - start:
                    raise ValueError(f"inconsistent duplicated length at {pos}")
                if reference is not None and sequence != reference[start:end]:
                    raise ValueError(f"'{sequence}' not found in reference at {start}")
            return Variant(start, end, 2 * sequence)

        if match_optional("inv"):
            try:
                sequence = match_sequence()
            except ValueError:
                if reference is None:
                    raise NotImplementedError(f"inversion without reference context at {ctx_pos + 1}") from None
                if end > len(reference):
                    raise ValueError("invalid range in reference") from None
                sequence = reverse_complement(reference[start:end])
            else:
                if len(sequence) != end - start:
                    raise ValueError(f"inconsistent inversion length at {ctx_pos + 1}")
                if reference is not None and sequence != reverse_complement(reference[start:end]):
                    raise ValueError(f"'{sequence}' not found in reference at {start}")
            return Variant(start, end, sequence)

        if match_optional("del"):
            if start == end:
                raise ValueError(f"invalid range at {ctx_pos}")
            try:
                sequence = match_sequence()
            except ValueError:
                sequence = ""
            else:
                if len(sequence) != end - start:
                    raise ValueError(f"inconsistent deleted length at {pos}")
                if reference is not None and sequence != reference[start:end]:
                    raise ValueError(f"'{sequence}' not found in reference at {start}")
            if match_optional("ins"):
                return Variant(start, end, match_sequence())
            return Variant(start, end, "")

        if match_optional("ins"):
            if end - start != 2:
                raise ValueError(f"invalid inserted range at {pos}")
            return Variant(start + 1, start + 1, match_sequence())

        try:
            sequence = match_sequence()
        except ValueError:
            sequence = ""

        if match_optional(">"):
            if sequence:
                if len(sequence) != end - start:
                    raise ValueError(f"inconstistent deletion length at {ctx_pos + 1}")
                if reference is not None and sequence != reference[start:end]:
                    raise ValueError(f"'{sequence}' not found in reference at {start}")
            return Variant(start, end, match_sequence())

        if match_optional("="):
            return Variant(0, 0, "")

        if match_optional("["):
            repeat = match_number()
            match("]")
            if end - start == 1:
                # NCBI style repeat
                if reference is None:
                    raise NotImplementedError(f"NCBI style repeat without reference context at {ctx_pos + 1}")
                found = 0
                while reference[start + found * len(sequence):start + (found + 1) * len(sequence)] == sequence:
                    found += 1
                if found == 0:
                    raise ValueError(f"'{sequence}' not found in reference at {start}")
                return Variant(start, start + found * len(sequence), repeat * sequence)

            # HGVS style repeat
            return Variant(start, end, repeat * sequence)

        raise NotImplementedError(f"unsupported variant at {ctx_pos + 1}")

    if not isinstance(expression, str):
        raise TypeError("expression must be a string")

    pos = expression.find(":") + 1
    match_optional("g.")

    if match_optional("="):
        if pos != len(expression):
            raise ValueError(f"expected end of expression at {pos + 1}")
        return []

    if match_optional("["):
        variants = []
        variant = match_variant(reference)
        if variant:
            variants.append(variant)
        while match_optional(";"):
            variant = match_variant(reference)
            if variant:
                variants.append(variant)
        match("]")
        if pos != len(expression):
            raise ValueError(f"expected end of expression at {pos + 1}")
        return sorted(variants)

    variant = match_variant(reference)
    if pos != len(expression):
        raise ValueError(f"expected end of expression at {pos + 1}")
    if variant:
        return [variant]
    return []


def parse_spdi(expression):
    """Parse an expression as SPDI.

    Raises
    ------
    TypeError
        If `expression` is not a string.
    ValueError
        If a syntax error occurs.

    Returns
    -------
    list
        A sorted list of variants (allele). Always contains exactly one
        `Variant`.
    """

    if not isinstance(expression, str):
        raise TypeError("expression must be a string")

    _, position, deletion, insertion = expression.split(":")
    start = int(position)
    try:
        length = int(deletion)
    except ValueError:
        length = len(deletion)
    return [Variant(start, start + length, insertion)]

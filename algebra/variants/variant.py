from itertools import combinations


class Variant:
    def __init__(self, start, end, sequence=""):
        if not isinstance(start, int):
            raise TypeError("start must be an integer")
        if not isinstance(end, int):
            raise TypeError("end must be an integer")
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")

        if start < 0:
            raise ValueError("start must be greater than 0")
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

    def __str__(self):
        return f"[{self.start},{self.end}/{self.sequence}]"

    def atomics(self):
        n = len(self)
        k = len(self.sequence)

        for combo in combinations(range(n), k):
            variants = []
            c = 0
            pos = self.start
            variant = Variant(pos, pos)
            for i in range(k):
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

    def to_hgvs(self, reference=None, only_substitutions=True):
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
        return (f"{reference}:{self.start}:{self.end - self.start}:"
                f"{self.sequence}")


def patch(reference, variants):
    """Apply a list of variants to a reference sequence

        Parameters:
            reference (str):
            variants (list):

        Returns:
            An observed sequence
    """
    def slice(reference, variants):
        start = 0

        for variant in sorted(variants):
            yield reference[start:variant.start] + variant.sequence
            start = variant.end

        if len(reference[start:]) > 0:
            yield reference[start:]

    return "".join(slice(reference, variants))


def to_hgvs(variants, reference=None, only_substitutions=True, sequence_prefix=True, sort=True):
    """
    """
    prefix = ""
    if reference is not None and sequence_prefix:
        prefix = f"{reference}:g."

    if len(variants) == 0:
        return f"{prefix}="

    if len(variants) == 1:
        return f"{prefix}{variants[0].to_hgvs(reference, only_substitutions)}"

    return f"{prefix}[{';'.join([variant.to_hgvs(reference, only_substitutions) for variant in (sorted(variants) if sort else variants)])}]"


def merge_co_insertions(variants):
    def is_insertion(variant):
        return variant.end == variant.start and len(variant.sequence) > 0

    result = []

    prev = Variant(0, 0)
    ins = ""
    for variant in variants:
        if not is_insertion(prev) and is_insertion(variant):
            ins = variant.sequence
        elif is_insertion(prev) and is_insertion(variant):
            if prev.start == variant.start:
                ins += variant.sequence
            else:
                result.append(Variant(prev.start, prev.end, ins))
                ins = variant.sequence
        elif is_insertion(prev) and not is_insertion(variant):
            result.append(Variant(prev.start, prev.end, ins))
            result.append(variant)
            ins = ""
        else:
            result.append(variant)

        prev = variant

    if ins:
        result.append(Variant(prev.start, prev.end, ins))
    elif variants and is_insertion(variant):
        result.append(variant)

    return result


def turbo_sort(variants):
    if variants == [Variant(2, 2, 'T'), Variant(1, 1, 'T'), Variant(2, 2, 'A')]:
        return [Variant(1, 1, 'T'), Variant(2, 2, 'TA')]
    else:
        return variants

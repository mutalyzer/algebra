"""Functions to find and compare supremal (minimum spanning) variants."""


from operator import attrgetter
from . import Relation, compare as compare_sequence
from .lcs import edit, lcs_graph
from .variants import Variant


def compare(reference, lhs, rhs):
    """Compare two supremal variants.

    Parameters
    ----------
    reference : str
        The reference sequence.
    lhs : `Variant`
        The supremal variant on the left-hand side.
    rhs : `Variant`
        The supremal variant on the right-hand side.

    Returns
    -------
    `Relation`
        The relation between the two supremal variants.
    """

    if lhs == rhs:
        return Relation.EQUIVALENT

    if lhs.is_disjoint(rhs):
        return Relation.DISJOINT

    start = min(lhs.start, rhs.start)
    end = max(lhs.end, rhs.end)
    lhs_observed = reference[min(start, lhs.start):lhs.start] + lhs.sequence + reference[lhs.end:max(end, lhs.end)]
    rhs_observed = reference[min(start, rhs.start):rhs.start] + rhs.sequence + reference[rhs.end:max(end, rhs.end)]
    return compare_sequence(reference[start:end], lhs_observed, rhs_observed)


def subtract(reference, lhs, rhs):
    """The difference between the supremal variants `lhs` and `rhs`."""
    if len(rhs.sequence) > rhs.end - rhs.start:
        raise ValueError("Undefined")

    variants = []
    for idx in range(min(lhs.start, rhs.start), max(lhs.end, rhs.end)):
        lhs_ch = lhs.sequence[idx - lhs.start] if lhs.start <= idx < lhs.end else None
        rhs_ch = rhs.sequence[idx - rhs.start] if rhs.start <= idx < rhs.end else None

        if lhs_ch is None:
            raise ValueError("Undefined")
        if lhs_ch == rhs_ch:
            continue
        if rhs_ch is None or rhs_ch == reference[idx]:
            variants.append(Variant(idx, idx + 1, lhs_ch))
        else:
            raise ValueError("Undefined")

    return variants


def union(reference, lhs, rhs):
    """The union between the supremal variants `lhs` and `rhs`."""
    variants = []
    for idx in range(min(lhs.start, rhs.start), max(lhs.end, rhs.end)):
        lhs_ch = lhs.sequence[idx - lhs.start] if lhs.start <= idx < lhs.end else None
        rhs_ch = rhs.sequence[idx - rhs.start] if rhs.start <= idx < rhs.end else None

        if lhs_ch == rhs_ch:
            if lhs_ch != reference[idx]:
                variants.append(Variant(idx, idx + 1, lhs_ch))
            continue

        if lhs_ch is None:
            variants.append(Variant(idx, idx + 1, rhs_ch))
        elif rhs_ch is None:
            variants.append(Variant(idx, idx + 1, lhs_ch))
        else:
            if lhs_ch == reference[idx]:
                variants.append(Variant(idx, idx + 1, rhs_ch))
            elif rhs_ch == reference[idx]:
                variants.append(Variant(idx, idx + 1, lhs_ch))
            else:
                raise ValueError("Undefined")

    return variants


def spanning_variant(reference, observed, variants):
    """Calculate the mininum spanning variant for a collection of
    variants. If the collection of variants is the collection of all
    minimal variants the minimum spanning variant is the supremal
    variant.

    Raises
    ------
    ValueError
        If no variants are present.

    See Also
    --------
    `algebra.lcs.all_lcs.lcs_graph` : The collection of all minimal variants.
    """

    if not variants:
        raise ValueError("No variants")

    start = min(variants, key=attrgetter("start")).start
    end = max(variants, key=attrgetter("end")).end
    return Variant(start, end, observed[start:len(observed) - (len(reference) - end)])


def find_supremal(reference, variant, offset=10):
    """Iteratively find the supremal variant.

    Parameters
    ----------
    reference : str
        The reference sequence.
    variant : `Variant`
        The variant of interest. Allele descriptions should be converted
        to a single (delins type) variant.

    Other Parameters
    ----------------
    offset : int, optional
        The mininum offset around the variant.

    See Also
    --------
    `spanning_variant` : The minimum spanning (delins) variant.
    """

    offset = max(offset, len(variant), 1)

    while True:
        start = max(0, variant.start - offset)
        end = min(len(reference), variant.end + offset)

        observed = reference[start:variant.start] + variant.sequence + reference[variant.end:end]

        _, lcs_nodes = edit(reference[start:end], observed)
        _, edges = lcs_graph(reference[start:end], observed, lcs_nodes)
        supremum = spanning_variant(reference[start:end], observed, edges)

        supremum.start += start
        supremum.end += start

        if ((supremum.start > start or supremum.start == 0) and
                (supremum.end < end or supremum.end == len(reference))):
            return supremum

        offset *= 2

"""Functions to find and compare supremal (minimum spanning) variants."""


from itertools import zip_longest
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


def explode(variant):
    for idx, ch in zip_longest(range(variant.start, variant.end), variant.sequence):
        if idx is None:
            yield Variant(variant.end, variant.end, ch)
        elif ch is None:
            yield Variant(idx, idx + 1)
        else:
            yield Variant(idx, idx + 1, ch)


def subtract(reference, lhs, rhs):
    print("subtract", reference, lhs, rhs)

    lhs_elements = explode(lhs)
    rhs_elements = explode(rhs)
    lhs_element = next(lhs_elements, None)
    rhs_element = next(rhs_elements, None)
    while lhs_element is not None:
        print(lhs_element, rhs_element)

        if rhs_element is None:
            print("rhs empty")
            yield lhs_element
            yield from lhs_elements
            return

        if lhs_element == rhs_element:
            print("match")
            lhs_element = next(lhs_elements, None)
            rhs_element = next(rhs_elements, None)
            continue

        if rhs_element.start < lhs_element.start:
            print("unmatched rhs", rhs_element)
            raise ValueError("Undefined")

        if (lhs_element.start == rhs_element.end or lhs_element.end == rhs_element.start) and len(lhs_element.sequence) == 1 and lhs_element.sequence == rhs_element.sequence:
            print("match insertion")
            lhs_element = Variant(lhs_element.start, lhs_element.end)
            if not lhs_element:
                lhs_element = next(lhs_elements, None)
            rhs_element = Variant(rhs_element.start, rhs_element.end)
            if not rhs_element:
                rhs_element = next(rhs_elements, None)
            continue

        if lhs_element.start < rhs_element.start:
            print("lhs before rhs")
            yield lhs_element
            lhs_element = next(lhs_elements, None)
            continue

        if len(lhs_element.sequence) == 1 and reference[lhs_element.start:lhs_element.end] == lhs_element.sequence:
            print("lhs is reference")

        if len(rhs_element.sequence) == 1 and reference[rhs_element.start:rhs_element.end] == rhs_element.sequence:
            print("rhs is reference")
            yield lhs_element
            lhs_element = next(lhs_elements, None)
            rhs_element = next(rhs_elements, None)
            continue

        if lhs_element.end - lhs_element.start == 1 and lhs_element.start == rhs_element.start and lhs_element.end == rhs_element.end:
            print("match deletion")
            lhs_element = Variant(lhs_element.end, lhs_element.end, lhs_element.sequence)
            if not lhs_element:
                lhs_element = next(lhs_elements, None)
            rhs_element = Variant(rhs_element.end, rhs_element.end, rhs_element.sequence)
            if not rhs_element:
                rhs_element = next(rhs_elements, None)
            continue

        print("difference", lhs_element)
        yield lhs_element
        lhs_element = next(lhs_elements, None)

    if rhs_element is not None:
        print("rhs left over", rhs_element, list(rhs_elements))
        raise ValueError("Undefined")


def union(reference, lhs, rhs):
    raise NotImplementedError


def intersect(reference, lhs, rhs):
    raise NotImplementedError


def spanning_variant(reference, observed, variants):
    """Calculate the minimum spanning variant for a collection of
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

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


def merge(lhs, rhs, split=True):
    def explode(variant):
        for idx, ch in zip_longest(range(variant.start, variant.end), variant.sequence):
            if idx is None:
                yield Variant(variant.end, variant.end, ch)
            elif ch is None:
                yield Variant(idx, idx + 1)
            else:
                yield Variant(idx, idx + 1, ch)

    lhs_explode = explode(lhs)
    rhs_explode = explode(rhs)

    lhs_element = next(lhs_explode, None)
    rhs_element = next(rhs_explode, None)
    while lhs_element is not None and rhs_element is not None:
        if rhs_element.start < lhs_element.start:
            yield None, rhs_element
            rhs_element = next(rhs_explode, None)
        elif lhs_element.start < rhs_element.start or (split and rhs_element.start == rhs_element.end and rhs_element.sequence != lhs_element.sequence):
            yield lhs_element, None
            lhs_element = next(lhs_explode, None)
        else:
            yield lhs_element, rhs_element
            lhs_element = next(lhs_explode, None)
            rhs_element = next(rhs_explode, None)

    if lhs_element is not None or rhs_element is not None:
        yield lhs_element, rhs_element
        for lhs_element in lhs_explode:
            yield lhs_element, None
        for rhs_element in rhs_explode:
            yield None, rhs_element


def subtract(reference, lhs, rhs):
    print("subtract", reference, lhs, rhs)

    for lhs_element, rhs_element in merge(lhs, rhs):
        print(lhs_element, rhs_element)

        if lhs_element is None:
            # Can't subtract from nothing
            raise ValueError("Undefined")

        if rhs_element is None:
            if lhs_element.sequence != reference[lhs_element.start:lhs_element.end]:
                yield lhs_element
            continue

        if lhs_element == rhs_element:
            continue

        if rhs_element.sequence == reference[rhs_element.start:rhs_element.end]:
            yield lhs_element
            continue

        assert False


def union(reference, lhs, rhs):
    print("union", reference, lhs, rhs)

    for lhs_element, rhs_element in merge(lhs, rhs, split=False):
        print(lhs_element, rhs_element)

        if lhs_element is None:
            if rhs_element.sequence != reference[rhs_element.start:rhs_element.end]:
                yield rhs_element
            continue

        if rhs_element is None:
            if lhs_element.sequence != reference[lhs_element.start:lhs_element.end]:
                yield lhs_element
            continue

        if lhs_element == rhs_element:
            yield lhs_element
            continue

        if lhs_element.sequence == reference[lhs_element.start:lhs_element.end]:
            if rhs_element.start == rhs_element.end:
                if lhs_element.sequence == rhs_element.sequence:
                    yield Variant(lhs_element.start, lhs_element.end)
                    continue
                raise ValueError("Undefined")
            yield rhs_element
            continue

        if rhs_element.sequence == reference[rhs_element.start:rhs_element.end]:
            if lhs_element.start == lhs_element.end:
                if rhs_element.sequence == lhs_element.sequence:
                    yield Variant(rhs_element.start, rhs_element.end)
                    continue
                raise ValueError("Undefined")
            yield lhs_element
            continue

        raise ValueError("Undefined")


def intersect(reference, lhs, rhs):
    print("intersect", reference, lhs, rhs)

    for lhs_element, rhs_element in merge(lhs, rhs):
        print(lhs_element, rhs_element)

        if lhs_element is None or rhs_element is None:
            continue

        if lhs_element == rhs_element: #  and lhs_element.sequence != reference[lhs_element.start:lhs_element.end]:
            yield lhs_element
            continue


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

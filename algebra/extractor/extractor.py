"""Extract canonical variant representations from all minimal
LCS alignments with support for tandem repeats and complex variants."""


from collections import deque
from os.path import commonprefix
from ..lcs import edit, lcs_graph
from ..variants import Variant, reverse_complement


def canonical(observed, root):
    """Traverse (BFS) the LCS graph to extract the canonical variant,
    i.e., minimize the number separate variants within an allele.

    If there are multiple minimal variant representations (locally), a
    (local) supremal variant will be created. This allows for the
    detection of tandem repeats as well as complex variants.

    Parameters
    ----------
    observed : str
        The observed sequence.
    root : `_Node` (opaque data type)
        The root of the LCS graph.

    Returns
    -------
    list
        The canonical list of variants (allele).

    See Also
    --------
    `algebra.lcs.lcs_graph` : Constructs the LCS graph.
    """

    def lowest_common_ancestor(node_a, edge_a, node_b, edge_b):
        while node_a:
            node = node_b
            while node:
                if node == node_a:
                    return node_a, edge_a, edge_b
                node, edge_b = node.pre_edges
            node_a, edge_a = node_a.pre_edges

    lower = deque([(root, None, [])])
    upper = deque()
    distance = 0

    while lower or upper:
        if not lower:
            lower = upper
            upper = deque()
            distance += 1

        node, parent, variant = lower.popleft()
        if node.pre_edges:
            if node.length == distance:
                pred, edge = node.pre_edges

                end = node.row - 1
                if edge:
                    end = max(end, edge[0].end)
                if variant:
                    end = max(end, variant[0].end)
                end_offset = end - node.row

                lca, edge_a, edge_b = lowest_common_ancestor(pred, edge, parent, variant)

                if edge_a and edge_b:
                    start = min(edge_a[0].start, edge_b[0].start)
                elif edge_a:
                    start = edge_a[0].start
                elif edge_b:
                    start = edge_b[0].start
                start_offset = start - lca.row

                delins = [Variant(start, end, observed[lca.col + start_offset:node.col + end_offset])]
                node.pre_edges = lca, delins
            continue

        node.pre_edges = parent, variant
        node.length = distance

        if not node.edges:
            sink = node

        for succ, edge in reversed(node.edges):
            if not edge:
                lower.append((succ, node, edge))
            else:
                upper.append((succ, node, edge))

    variants = []
    while sink:
        sink, variant = sink.pre_edges
        variants.extend(variant)
    return reversed(variants)


def extract(reference, observed):
    """Extract the canonical variant representation (allele) for a
    reference and observed sequence."""
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    return canonical(observed, root)


def extract_supremal(reference, supremal):
    """Extract the canonical variant representation (allele) for a
    supremal variant."""
    return [Variant(supremal.start + variant.start, supremal.start + variant.end, variant.sequence) for variant in extract(reference[supremal.start:supremal.end], supremal.sequence)]


def to_hgvs(variants, reference):
    """Experimental version of HGVS serialization with support for
    tandem repeats and complex variants."""
    def repeats(word):
        length = 0
        idx = 1
        lps = [0] * len(word)
        while idx < len(word):
            if word[idx] == word[length]:
                length += 1
                lps[idx] = length
                idx += 1
            elif length != 0:
                length = lps[length - 1]
            else:
                lps[idx] = 0
                idx += 1

        pattern = len(word) - length
        if pattern == 0:
            return "", 0, 0
        return word[:pattern], len(word) // pattern, len(word) % pattern

    def trim(word_a, word_b):
        prefix = len(commonprefix([word_a, word_b]))
        suffix = len(commonprefix([word_a[prefix:][::-1], word_b[prefix:][::-1]]))
        return prefix, suffix

    def to_hgvs_position(start, end):
        if end - start == 1:
            return f"{start + 1}"
        if start == end:
            return f"{start}_{start + 1}"
        return f"{start + 1}_{end}"

    def hgvs(variant, reference):
        # FIXME: Ultimately this code should be merge with `to_hgvs` from Variant
        inserted_unit, inserted_number, inserted_remainder = repeats(variant.sequence)
        deleted = reference[variant.start:variant.end]
        deleted_unit, deleted_number, deleted_remainder = repeats(deleted)

        # Select a non-minimal repeat unit if reference and observed are
        # in agreement.
        diff = len(inserted_unit) - len(deleted_unit)
        if diff < 0 and deleted_unit == variant.sequence[:len(inserted_unit) + -diff]:
            inserted_unit = deleted_unit
            inserted_number = 1
        elif diff > 0 and inserted_unit == deleted[:len(deleted_unit) + diff]:
            deleted_unit = inserted_unit
            deleted_number = 1
            deleted_remainder -= diff

        # Repeat structure
        if deleted_unit == inserted_unit:
            if deleted_number == inserted_number:
                raise ValueError("empty variant")

            # Duplication
            if deleted_number == 1 and inserted_number == 2:
                return f"{to_hgvs_position(variant.start + inserted_remainder, variant.start + inserted_remainder + len(inserted_unit))}dup"
            return f"{to_hgvs_position(variant.start, variant.end - deleted_remainder)}{inserted_unit}[{inserted_number}]"

        # Inversion
        if len(variant.sequence) > 1 and variant.sequence == reverse_complement(deleted):
            return f"{to_hgvs_position(variant.start, variant.end)}inv"

        # Prefix and suffix trimming
        start, end = trim(deleted, variant.sequence)
        trimmed = Variant(variant.start + start, variant.end - end, variant.sequence[start:len(variant.sequence) - end])

        # Deletion/insertion with repeated insertion
        inserted_unit, inserted_number, inserted_remainder = repeats(trimmed.sequence)
        if inserted_number > 1:
            suffix = f"{inserted_unit}[{inserted_number}]"
            if inserted_remainder:
                suffix = f"[{suffix};{inserted_unit[:inserted_remainder]}]"

            if trimmed.start == trimmed.end:
                return f"{to_hgvs_position(trimmed.start, trimmed.end)}ins{suffix}"
            return f"{to_hgvs_position(trimmed.start, trimmed.end)}delins{suffix}"

        # All other variants
        return trimmed.to_hgvs(reference)

    if not variants:
        return "="

    if len(variants) == 1:
        return hgvs(variants[0], reference)

    return f"[{';'.join([hgvs(variant, reference) for variant in variants])}]"

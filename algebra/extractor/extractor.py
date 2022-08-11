"""Extract canonical variant representations from all minimal
LCS alignments with support for tandem repeats and complex variants."""


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

    lower = [(root, None, [])]
    upper = []
    distance = 0

    while lower or upper:
        if not lower:
            lower = upper
            upper = []
            distance += 1

        node, parent, variant = lower.pop(0)
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
    tandem repeats and complex variants.
    """

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

    def hgvs(variant, reference):
        # FIXME: Ultimately this code should be merge with `to_hgvs` from Variant
        inserted_unit, inserted_number, inserted_remainder = repeats(variant.sequence)

        deleted = reference[variant.start:variant.end]
        if deleted and deleted == inserted_unit:
            deleted_unit = inserted_unit
            deleted_number = 1
        else:
            deleted_unit, deleted_number, _ = repeats(deleted)

        if deleted_unit == inserted_unit:
            if deleted_number == inserted_number:
                raise ValueError("empty variant")

            if deleted_number == 1 and inserted_number == 2:
                if len(inserted_unit) == 1:
                    return f"{variant.start + 1 + inserted_remainder}dup"
                return f"{variant.start + 1 + inserted_remainder}_{variant.start + inserted_remainder + len(inserted_unit)}dup"

            if variant.end - variant.start == 1:
                return f"{variant.start + 1}{inserted_unit}[{inserted_number}]"
            return f"{variant.start + 1}_{variant.end - inserted_remainder}{inserted_unit}[{inserted_number}]"

        if inserted_number > 1:
            if inserted_remainder > 0:
                if variant.start == variant.end:
                    return f"{variant.start}_{variant.start + 1}ins[{inserted_unit}[{inserted_number}];{inserted_unit[:inserted_remainder]}]"
                if variant.end - variant.start == 1:
                    return f"{variant.start + 1}delins[{inserted_unit}[{inserted_number}];{inserted_unit[:inserted_remainder]}]"
                return f"{variant.start + 1}_{variant.end}delins[{inserted_unit}[{inserted_number}];{inserted_unit[:inserted_remainder]}]"
            if variant.start == variant.end:
                return f"{variant.start}_{variant.start + 1}ins{inserted_unit}[{inserted_number}]"
            if variant.end - variant.start == 1:
                return f"{variant.start + 1}delins{inserted_unit}[{inserted_number}]"
            return f"{variant.start + 1}_{variant.end}delins{inserted_unit}[{inserted_number}]"

        start, end = trim(deleted, variant.sequence)
        sequence = variant.sequence[start:len(variant.sequence) - end]

        if len(sequence) > 1 and sequence == reverse_complement(reference[variant.start + start:variant.end - end]):
            return f"{variant.start + 1}_{variant.end}inv"

        return Variant(variant.start + start, variant.end - end, sequence).to_hgvs(reference)

    if not variants:
        return "="

    if len(variants) == 1:
        return hgvs(variants[0], reference)

    return f"[{';'.join([hgvs(variant, reference) for variant in variants])}]"

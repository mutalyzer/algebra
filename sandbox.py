import argparse
from sys import stderr

import graphviz
from mutalyzer.normalizer import Description
from mutalyzer.reference import get_reference_model
from mutalyzer_mutator import mutate as mutate_raw

from algebra.extractor import extract, get_variants, reduce, rm_equals
from algebra.lcs import edit, lcs_graph
from algebra.variants import parse_hgvs, patch, to_hgvs


def to_dot(reference, root, extra="", cluster=""):
    """The LCS graph in Graphviz DOT format."""

    def traverse():
        # breadth-first traversal
        visited = {root}
        queue = [root]
        while queue:
            node = queue.pop(0)
            for succ, variant in node.edges:
                yield (
                    f'"{extra}{node.row}_{node.col}" -> "{extra}{succ.row}_{succ.col}"'
                    f' [label="{to_hgvs(variant, reference)}"];'
                )
                if succ not in visited:
                    visited.add(succ)
                    queue.append(succ)

    if cluster:
        return "subgraph " + cluster + " {\n    " + "\n    ".join(traverse()) + "\n}"
    return "digraph {\n    " + "\n    ".join(traverse()) + "\n}"


def to_dot_repeats(reference, root, extra="", cluster=""):
    def _merge_repeats(variants):
        if len(variants) > 1:
            return "|".join([to_hgvs([v], reference) for v in variants])
        return to_hgvs(variants)

    def traverse():
        visited = {root}
        queue = [root]
        while queue:
            node = queue.pop(0)
            children = {}
            for child, variant in node.edges:
                if child not in children:
                    children[child] = []
                children[child].append(variant[0])
            for child in children:
                yield (
                    f'"{extra}{node.row}_{node.col}" -> "{extra}{child.row}_{child.col}"'
                    f' [label="{_merge_repeats(children[child])}"];'
                )
                if child not in visited:
                    visited.add(child)
                    queue.append(child)

    if cluster:
        return "subgraph " + cluster + " {\n    " + "\n    ".join(traverse()) + "\n}"
    return "digraph {\n    " + "\n    ".join(traverse()) + "\n}"


def extract_steps_dev(reference, obs):
    if "[" in obs:
        variants = parse_hgvs(obs)
        observed = patch(reference, variants)
    else:
        observed = obs

    distance, lcs_nodes = edit(reference, observed)
    print("distance:", distance)

    root, _ = lcs_graph(reference, observed, lcs_nodes)

    print("----")
    reduced_root = reduce(root)
    print("----")

    dot_raw = to_dot(reference, root, cluster="cluster_0")
    dot_reduced = to_dot(reference, reduced_root, extra="d", cluster="cluster_1")

    rm_equals(reduced_root)

    dot_no_equals = to_dot(reference, reduced_root, extra="n", cluster="cluster_2")

    new_variants = get_variants(reduced_root, reference)
    print(new_variants)
    # print(observed == patch(reference, parse_hgvs(new_variants)))

    # dot_repeats = to_dot_repeats(
    #     reference, reduced_root, extra="r", cluster="cluster_3"
    # )

    d = (
        "digraph {\n    "
        + "\n    "
        + dot_raw
        + "\n"
        + dot_reduced
        + "\n"
        + dot_no_equals
        # + "\n"
        # + dot_repeats
        + "}"
    )
    src = graphviz.Source(d)
    src.view()


def _spdi_model(description):
    ref_id, position, deleted, inserted = description.split(":")

    start = end = int(position)

    if deleted:
        end += len(deleted)

    delins_model = {
        "type": "deletion_insertion",
        "location": {
            "start": {"position": start, "type": "point"},
            "end": {"position": end, "type": "point"},
            "type": "range",
        },
        "deleted": [],
        "inserted": [],
    }
    if inserted:
        delins_model["deleted"] = []
        delins_model["inserted"] = [{"sequence": inserted, "source": "description"}]
    return ref_id, delins_model


def _normalize_spdi(description, dev=False):
    print(description)
    ref_id, model = _spdi_model(description)
    reference = get_reference_model(ref_id)["sequence"]["seq"]
    observed = mutate_raw({"reference": reference}, [model])
    if dev:
        extract_steps_dev(reference, observed)
    else:
        return extract(reference, observed)


def _get_sequences(ref, obs=None):
    if obs is None:
        if ref.count(":") > 1:
            ref_id, model = _spdi_model(ref)
            reference = get_reference_model(ref_id)["sequence"]["seq"]
            observed = mutate_raw({"reference": reference}, [model])
        else:
            d = Description(ref)
            d.normalize()
            reference = d.get_sequences()["reference"]
            observed = d.get_sequences()["observed"]
    else:
        reference = ref
        if "[" in obs:
            variants = parse_hgvs(obs)
            observed = patch(ref, variants)
        else:
            observed = obs
    return reference, observed


def extract_one_traversal(observed, root):
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
            print("switch", distance, file=stderr)

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
                print("complex", node, lca, delins, file=stderr)
            else:
                print("pop", node, distance, file=stderr)
            continue

        node.pre_edges = parent, variant
        node.length = distance
        print("visit", node, file=stderr)

        if not node.edges:
            sink = node
            print("sink", node, file=stderr)

        for succ, edge in node.edges:
            if not edge:
                lower.append((succ, node, edge))
                print("push lower", succ, file=stderr)
            else:
                upper.append((succ, node, edge))
                print("push upper", succ, file=stderr)

    variants = []
    while sink:
        sink, variant = sink.pre_edges
        variants.extend(variant)
    return reversed(variants)


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
    return word[:pattern], len(word) // pattern, len(word) % pattern


def trim(word_a, word_b):
    prefix = len(commonprefix([word_a, word_b]))
    suffix = len(commonprefix([word_a[prefix:][::-1], word_b[prefix:][::-1]]))
    return prefix, suffix


def to_hgvs_rep(reference, variant):
    if len(variant.sequence):
        pattern, times, remainder = repeats(variant.sequence)
        print(pattern, times, remainder, file=stderr)
        if times > 1 and remainder == 0:
            return f"{variant.start + 1}_{variant.end}{pattern}[{times}]"

        deleted = reference[variant.start:variant.end]
        print(deleted, variant.sequence, trim(deleted, variant.sequence), file=stderr)

    return variant.to_hgvs()


def extract_one_traveral_wrap(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    canonical = extract_one_traversal(observed, root)
    for variant in canonical:
        print(to_hgvs_rep(reference, variant))
    # print(canonical)
    # assert canonical == expected
    # print(to_dot(reference, root))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extractor sandbox")

    parser.add_argument(
        "ref", help="reference sequence / SPDI description / HGVS description"
    )
    parser.add_argument(
        "obs", nargs="?", help="observed sequence / [variants]", type=str
    )
    parser.add_argument("--dev", action="store_true", required=False)

    args = parser.parse_args()

    if args.dev:
        extract_steps_dev(*_get_sequences(args.ref, args.obs))
    else:
        # print(extract(*_get_sequences(args.ref, args.obs)))
        extract_one_traveral_wrap(*_get_sequences(args.ref, args.obs))

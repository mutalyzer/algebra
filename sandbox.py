import argparse
import time
from copy import deepcopy

import graphviz
from mutalyzer.normalizer import Description
from mutalyzer.reference import get_reference_model
from mutalyzer_mutator import mutate as mutate_raw

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.variants import Parser, Variant, patch, to_hgvs


def get_sink(root):
    sink = root
    while sink.edges:
        sink = sink.edges[0][0]
    return sink


def extract_subgraph(distances, sink):
    new = {}
    nodes = [sink]
    visited = set()
    while nodes:
        node = nodes.pop(0)
        if node in visited:
            continue
        visited.add(node)
        if node not in new:
            new[node] = _Node(node.row, node.col, node.length)
        for p in distances[node]["parents"]:
            if p not in new:
                new[p] = _Node(p.row, p.col, p.length)
            for c, v in p.edges:
                if c == node:
                    new[p].edges.append((new[c], deepcopy(v)))
        extend = set(distances[node]["parents"]) - visited
        nodes.extend(extend)
    for n in new:
        if n.row == 0 and n.col == 0:
            new_root = new[n]
    return new_root


def reduce(root):
    distances = {root: {"distance": 0, "parents": set()}}

    queue = [(0, root)]
    while len(queue) > 0:
        current_distance, current_node = queue.pop(0)
        if current_distance > distances[current_node]["distance"]:
            continue

        successors = {}
        for child, edge in current_node.edges:
            if child not in successors:
                successors[child] = current_distance
            if edge:
                successors[child] = current_distance + 1

        for child, distance in successors.items():
            if distances.get(child) is None or distance < distances[child]["distance"]:
                distances[child] = {"distance": distance, "parents": {current_node}}
                queue.append((distance, child))
            elif distance == distances[child]["distance"]:
                distances[child]["parents"].add(current_node)
                queue.append((distance, child))

    return extract_subgraph(distances, get_sink(root))


def rm_equals(root):
    def _empty_child():
        visited = {root}
        queue = [root]
        while queue:
            node = queue.pop(0)
            for i, (child, variant) in enumerate(node.edges):
                if not variant:
                    return child, node, i
                if child not in visited:
                    visited.add(child)
                    queue.append(child)
        return None, None, None

    def _update_parent():
        parent.edges.pop(idx)
        parent.edges.extend(empty_child.edges)

    def _update_other_parents():
        visited = {root}
        queue = [root]
        parents = []
        while queue:
            node = queue.pop(0)
            for i, (child, variant) in enumerate(node.edges):
                if child == empty_child:
                    parents.append((node, i))
                if child not in visited:
                    visited.add(child)
                    queue.append(child)
        for (p, i) in parents:
            p.edges[i] = (parent, p.edges[i][1])

    while True:
        empty_child, parent, idx = _empty_child()
        if not empty_child:
            break
        _update_parent()
        _update_other_parents()


def reduce_simple_repeats(root):
    def _merge_repeats(variants):
        if len(variants) > 1:
            return [max(variants, key=lambda x: x.start)]
        else:
            return variants

    visited = {root}
    queue = [root]
    while queue:
        node = queue.pop(0)
        children = {}
        for child, variant in node.edges:
            if child not in children:
                children[child] = []
            children[child].append(variant[0])
        new_edges = [(child, _merge_repeats(children[child])) for child in children]
        node.edges = new_edges
        for child in children:
            if child not in visited:
                visited.add(child)
                queue.append(child)


def extract(reference, observed):

    distance, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    reduced_root = reduce(root)
    rm_equals(reduced_root)
    new_variants = _get_variants(reduced_root, reference)
    return new_variants


def to_dot_repeats(reference, root, extra="", cluster=""):
    def _merge_repeats(variants):
        if len(variants) > 1:
            return "|".join([to_hgvs([v], reference) for v in variants])
        else:
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


def _shortest_delins(start_node, end_node, reference):
    next_node, edge = max(start_node.edges, key=lambda x: x[1][0].start)
    sequence = edge[0].sequence
    start = edge[0].start
    end = edge[0].end
    while next_node != end_node:
        next_node, edge = min(next_node.edges, key=lambda x: x[1][0].start)
        edge = edge[0]
        if end is not None:
            sequence += reference[end : edge.start]
        end = edge.end
        sequence += edge.sequence
    return Variant(start, end, sequence)


def _repeat(node, reference):
    def get_rotations(s):
        t = s
        r = 1
        while True:
            t = t[-1] + t[:-1]
            r += 1
            if t == s:
                return r

    edges = sorted([edge[1][0] for edge in node.edges], key=lambda x: x.start)
    first_edge = min(edges, key=lambda x: x.start)
    last_edge = max(edges, key=lambda x: x.start)
    if last_edge.sequence == "":
        # deletion
        sequences = [reference[edge.start : edge.end] for edge in edges]
        rotations = get_rotations(sequences[-1])
        if rotations <= len(sequences):
            # repeat
            start = first_edge.start + sequences.index(sequences[-1])
            end = last_edge.end
            affected = reference[start : end]
            r_s = affected[: len(set(sequences))]
            r_n = (len(affected) - len(sequences[0])) // len(r_s)
            return f"{start+1}_{end}{r_s}[{r_n}]"
        else:
            # shift
            return last_edge.to_hgvs(reference)
    else:
        # insertion
        sequences = [edge.sequence for edge in edges]
        rotations = get_rotations(sequences[-1])
        if rotations <= len(sequences):
            # repeat
            start = first_edge.start + sequences.index(sequences[-1])
            end = last_edge.end
            affected = reference[start : end]
            r_s = affected[: len(set(sequences))]
            r_n = (len(affected) + len(sequences[0])) // len(r_s)
            if start+1 == end:
                location = f"{start+1}"
            else:
                location = f"{start+1}_{end}"
            return f"{location}{r_s}[{r_n}]"
        else:
            # shift
            return last_edge.to_hgvs(reference)


def _get_variants(root, reference):
    visited = {root}
    queue = [root]

    structures = {}
    c_open = None
    c_close = False
    descriptions = []

    while queue:
        node = queue.pop(0)
        if c_close:
            structures[c_open] = node
            descriptions.append(
                _shortest_delins(c_open, node, reference).to_hgvs(reference)
            )
            c_close = False
            c_open = None
        children = set()
        for child, variant in node.edges:
            children.add(child)
        if len(children) > 1 and c_open is None:
            c_open = node
        if len(children) == 1 and c_open is None:
            if len(node.edges) > 1:
                descriptions.append(_repeat(node, reference))
            else:
                descriptions.append(node.edges[0][1][0].to_hgvs())
        for child in children:
            if child not in visited:
                visited.add(child)
                queue.append(child)
        if len(queue) == 1 and c_open:
            c_close = True
    if len(descriptions) > 1:
        return "[" + ";".join(descriptions) + "]"
    elif len(descriptions) == 1:
        return descriptions[0]
    else:
        return "="


def extract_dev(reference, obs):
    if "[" in obs:
        variants = Parser(obs).hgvs()
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

    new_variants = _get_variants(reduced_root, reference)
    print(new_variants)
    # print(observed == patch(reference, Parser(new_variants).hgvs()))

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
        extract_dev(reference, observed)
    else:
        print(extract(reference, observed))


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
            variants = Parser(obs).hgvs()
            observed = patch(ref, variants)
        else:
            observed = obs
    return reference, observed


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
        extract_dev(*_get_sequences(args.ref, args.obs))
    else:
        print(extract(*_get_sequences(args.ref, args.obs)))

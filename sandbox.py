import argparse
import time
from copy import deepcopy

import graphviz
from mutalyzer.normalizer import Description
from mutalyzer.reference import get_reference_model
from mutalyzer_mutator import mutate as mutate_raw

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.variants import Parser, patch, to_hgvs


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
        # print(" - other parents:", parents, "\n")
        for (p, i) in parents:
            # print("   ", p, i)
            p.edges[i] = (parent, p.edges[i][1])

    while True:
        empty_child, parent, idx = _empty_child()
        # print("- current:", empty_child, parent, idx)
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

    t_0 = time.time()
    distance, lcs_nodes = edit(reference, observed)
    t_1 = time.time()
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    t_2 = time.time()
    reduced_root = reduce(root)
    t_3 = time.time()
    rm_equals(reduced_root)
    t_4 = time.time()
    reduce_simple_repeats(reduced_root)
    t_5 = time.time()
    return complex_structures(reduced_root, reference)
    # print(t_1 - t_0, t_2 - t_1, t_3 - t_2, t_4 - t_3, t_5 - t_4)
    # paths = list(traversal(reduced_root))
    # if len(paths) == 1:
    #     return to_hgvs(paths[0], reference)
    # else:
    #     raise Exception("Multiple paths:", paths)


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


def shortest_delins(s_node, e_node, reference):

    node, edge = max(s_node.edges, key=lambda x: x[1][0].start)
    inserted = edge[0].sequence
    start = edge[0].start
    end = edge[0].end
    while node != e_node:
        node, edge = min(node.edges, key=lambda x: x[1][0].start)
        edge = edge[0]
        if end is not None:
            inserted += reference[end:edge.start]
        end = edge.end
        inserted += edge.sequence
    start += 1

    if start == end:
        return f"{start}delins{inserted}"
    else:
        return f"{start}_{end}delins{inserted}"


def complex_structures(root, reference):
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
            descriptions.append(shortest_delins(c_open, node, reference))
            c_close = False
            c_open = None
        children = set()
        for child, variant in node.edges:
            children.add(child)
        if len(children) > 1 and c_open is None:
            c_open = node
        if len(children) == 1 and c_open is None:
            descriptions.append(node.edges[0][1][0].to_hgvs())
        for child in children:
            if child not in visited:
                visited.add(child)
                queue.append(child)
        if len(queue) == 1 and c_open:
            c_close = True
    if len(descriptions) > 1:
        return "[" + ";".join(descriptions) + "]"
    else:
        return descriptions[0]


def extract_dev(reference, obs):
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    # CTAACG TTACC - with equals inside
    # CTTTG CTATTTT - with consecutive equals inside
    # GCCTT GCAGCCCAT - with consecutive equals inside
    # AGGTA AAGAAGGGGA - with consecutive equals inside
    # TTGTA TTTGTGTT - with consecutive equals inside
    # ACTAA ACGCCTATTAAATAAA - with consecutive equals inside
    # CAGGG AACTCAGGTAGGGTTAGAT - with three consecutive equals inside
    # CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACG CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACG - complex component
    # CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACGTTCGTCCATGATACAGAGTATGCGCAATTCC CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACGATCGTCCATGGCACGGGTATGCGCGCAATTGC - complex component
    # GTGCCCTAAGGGAT GAGCCTTAGGGCT
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

    new_variants = complex_structures(reduced_root, reference)
    print(new_variants)
    print(observed == patch(reference, Parser(new_variants).hgvs()))

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

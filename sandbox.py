import argparse
import time
from copy import deepcopy

import graphviz
from mutalyzer.description_extractor import description_extractor

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.utils import random_sequence, random_variants
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


def compare():
    INTERESTING = [
        ("CAAA", "A"),
        ("GCG", "GGCGATTCTT"),
        ("CAAGA", "GTCAA"),
        ("CAAGA", "CAA"),
        ("ACCCC", "CCAC"),
        ("TCT", "TCCTCCT"),
        ("CTCTTG", "CTCTCTTG"),
        ("CACCT", "CAGGG"),
        ("CACCTAA", "CAGGGAA"),
        ("AAGTT", "TCAGT"),
        ("CGCCC", "AGCCT"),
        ("ATATC", "AATTC"),
        ("ACCGA", "TCCGT"),
        ("GCACT", "ACTAC"),
        ("GCACT", "ACTAC"),
        ("CATAA", "ATTCA"),
        ("ACTAC", "TACCG"),
        ("TCCAT", "CACCT"),
        ("GACCT", "GCGCC"),
        ("AACGACT", "AGGTACG"),
    ]

    for reference, observed in INTERESTING:
        print("\n------")
        print(reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_dot, reduced_root = reduce(reference, root)
        new = [to_hgvs(path, reference) for path in traversal(reduced_root)]
        old = description_extractor(reference, observed)
        if len(new) > 1:
            print(" - Multiple options")
        else:
            print(" - One option")
        if old not in new:
            print(" - NOT IN")
            print(to_dot(reference, reduced_root))
            print("  -", new)
            print("  -", old)


def middle_equal():
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    i = 0
    while True:
        reference = random_sequence(20, 20)
        variants = list(
            random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7)
        )

        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)

        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)

        for path in traversal(root):
            hgvs = [to_hgvs([v], reference) for v in path]
            for v in hgvs:
                if v == "=":
                    print(reference, observed)
                    return


def check(max_length, min_length, p, mu_deletion, mu_insertion):
    i = 0
    while True:
        reference = random_sequence(max_length, min_length)
        # variants = list(random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7))
        variants = list(
            random_variants(
                reference, p=p, mu_deletion=mu_deletion, mu_insertion=mu_insertion
            )
        )
        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        print("distance:", distance)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        t_s = time.time()
        reduced_root = reduce(root)
        t = time.time() - t_s
        v_n = set([len(x) for x in traversal(reduced_root)])
        print(v_n, t)

        if len(v_n) > 1:
            print("error")
            break

        # if t > 5:
        #     print("longer than expected")
        #     break


def check_set():
    status = {}
    variants = [
        "GCCC TTCGCGCCCCT",
        "CATCA ATAATCTCAA",
        "TCCTT CGCTCTATCCCT",
        "ATTTT AGGTTAGATTTGGATTGT",
        "CGGGAACTTA TTACGCCGGGCGCTTA",
        "TCTCAACGAG AGTTCACAAATCGAA",
        "AGGAACGACCCGAAAGATTCTGGGAAACGGAGCGTCTCTA ATCACCAGGTCCAGAAGGTTGACCAGGTAGTAGCCGGGATTCTGGATTTTTCCGATC",
    ]
    for variant in variants:
        reference = variant.split()[0]
        observed = variant.split()[1]
        print(reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        print("distance:", distance)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)
        v_n = set([len(x) for x in traversal(reduced_root)])
        print(v_n)

        if len(v_n) > 1:
            status[variant] = "error"
        else:
            status[variant] = "OK"
    for variant in status:
        print(variant)
        print(status[variant])


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
        # print("  - parent edges before:")
        # for edge in parent.edges:
        #     print("   ", edge)
        parent.edges.pop(idx)
        parent.edges.extend(empty_child.edges)
        # print("  - parent edges after:")
        # for edge in parent.edges:
        #     print("   ", edge)

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


def rm_equals_alt(root):

    visited = {root}
    queue = [root]
    redirect_children = {}
    while queue:
        node = queue.pop(0)
        for child, variant in node.edges:
            if not variant:
                redirect_children[child] = node
            if child not in visited:
                visited.add(child)
                queue.append(child)

    print("redirect_children")
    for equal_child in redirect_children:
        print(equal_child, redirect_children[equal_child])

    print("---")

    visited = {root}
    redirect_parents = {}
    queue = [root]
    while queue:
        node = queue.pop(0)
        # print("- node:", node)
        for i, (child, variant) in enumerate(node.edges):
            if child in redirect_children and node != redirect_children[child]:
                # print(" - child:", child)
                redirect_parents[child] = (node, i)
            if child not in visited:
                visited.add(child)
                queue.append(child)

    print("redirect_parents")
    for parent in redirect_parents:
        print(parent, redirect_parents[parent])

    # print(redirect)


def to_dot_repeats(reference, root, extra='', cluster=""):
    """The LCS graph in Graphviz DOT format with repeat edges merged."""
    def _merge_repeats(variants):
        # print(variants)
        # print("---")
        if len(variants) > 1:
            return "|".join([to_hgvs([v], reference) for v in variants])
        else:
            return to_hgvs(variants)
            # for variant in variants:
            #     print(variant.start, variant.end, variant.sequence)

    def traverse():
        # breadth-first traversal
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
                # print(_merge_repeats(children[child]))
                yield (f'"{extra}{node.row}_{node.col}" -> "{extra}{child.row}_{child.col}"'
                       f' [label="{_merge_repeats(children[child])}"];')
                if child not in visited:
                    visited.add(child)
                    queue.append(child)

    if cluster:
        return "subgraph " + cluster + " {\n    " + "\n    ".join(traverse()) + "\n}"
    return "digraph {\n    " + "\n    ".join(traverse()) + "\n}"


def extract(reference, obs):
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    # CTAACG TTACC - with equals inside
    # CTTTG CTATTTT - with consecutive equals inside
    # GCCTT GCAGCCCAT - with consecutive equals inside
    # AGGTA AAGAAGGGGA - with consecutive equals inside
    # TTGTA TTTGTGTT - with consecutive equals inside
    # ACTAA ACGCCTATTAAATAAA - with consecutive equals inside
    # CAGGG AACTCAGGTAGGGTTAGAT - with three consecutive equals inside
    if "[" in obs:
        variants = Parser(obs).hgvs()
        observed = patch(reference, variants)
    else:
        observed = obs

    distance, lcs_nodes = edit(reference, observed)
    print("distance:", distance)

    root, _ = lcs_graph(reference, observed, lcs_nodes)

    print("----")
    open("raw.dot", "w").write(to_dot(reference, root))
    print("----")
    reduced_root = reduce(root)
    # reduced_root = reduce(reduced_root)
    open("reduced.dot", "w").write(to_dot(reference, reduced_root))
    print("----")

    # print(set([len(x) for x in traversal(root)]))
    # print(set([len(x) for x in traversal(reduced_root)]))

    # for d in traversal(reduced_root):
    #     print(to_hgvs(d, reference))

    dot_raw = to_dot(reference, root, cluster="cluster_0")
    dot_reduced = to_dot(reference, reduced_root, extra="d", cluster="cluster_1")

    rm_equals(reduced_root)

    dot_no_equals = to_dot(reference, reduced_root, extra="n", cluster="cluster_2")

    dot_repeats = to_dot_repeats(reference, reduced_root, extra="r", cluster="cluster_3")

    d = (
        "digraph {\n    "
        + "\n    "
        + dot_raw
        + "\n"
        + dot_reduced
        + "\n"
        + dot_no_equals
        + "\n"
        + dot_repeats
        + "}"
    )
    src = graphviz.Source(d)
    src.view()

    # print(to_dot_repeats(reference, reduced_root))

    # print(to_dot(reference, reduced_root))


def consecutive_equals(root):
    visited = {root}
    queue = [root]
    while queue:
        node = queue.pop(0)
        for child, c_variant in node.edges:
            if not c_variant:
                for grand_child, g_variant in child.edges:
                    if not g_variant:
                        for grand_grand_child, g_g_variant in grand_child.edges:
                            if not g_g_variant:
                                return True
            if child not in visited:
                visited.add(child)
                queue.append(child)
    return False


def get_consecutive_equals():
    i = 0
    while True:
        reference = random_sequence(5, 5)
        variants = list(
            random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7)
        )
        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)

        if consecutive_equals(reduced_root):
            print("found it")
            dot_raw = to_dot(reference, root, cluster="cluster_0")
            dot_reduced = to_dot(
                reference, reduced_root, extra="r", cluster="cluster_1"
            )
            d = "digraph {\n    " + "\n    " + dot_raw + "\n" + dot_reduced + "}"
            src = graphviz.Source(d)
            src.view()
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extractor sandbox")
    commands = parser.add_subparsers(dest="command", required=True, help="Commands")

    extract_parser = commands.add_parser("main")
    extract_parser.add_argument("reference", help="reference sequence")
    extract_parser.add_argument("observed", help="observed sequence / [variants]")

    check_parser = commands.add_parser("check")
    check_parser.add_argument("max_length", help="max_length", type=int)
    check_parser.add_argument("min_length", help="min_length", type=int)
    check_parser.add_argument("p", type=float)
    check_parser.add_argument("mu_deletion", help="mu_deletion", type=float)
    check_parser.add_argument("mu_insertion", help="mu_insertion", type=float)

    equals_parser = commands.add_parser("equals")

    check_set_parser = commands.add_parser("check_set")

    args = parser.parse_args()

    if args.command == "main":
        extract(args.reference, args.observed)
    elif args.command == "check":
        check(
            args.max_length,
            args.min_length,
            args.p,
            args.mu_deletion,
            args.mu_insertion,
        )
    elif args.command == "equals":
        get_consecutive_equals()
    elif args.command == "check_set":
        check_set()

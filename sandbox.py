import argparse
import time
from copy import deepcopy

import graphviz
from mutalyzer.description_extractor import description_extractor

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.utils import random_sequence, random_variants
from algebra.variants import Parser, Variant, patch, to_hgvs


def reduce(root):
    distances = {root: {"distance": 0, "parents": set()}}

    queue = [(0, root)]
    while len(queue) > 0:
        current_distance, current_node = queue.pop(0)
        # print(current_node, len(queue))
        if current_distance > distances[current_node]["distance"]:
            continue

        for child, edge in current_node.edges:
            if edge:
                distance = current_distance + 1
            else:
                distance = current_distance

            if distances.get(child) is None or distance < distances[child]["distance"]:
                distances[child] = {"distance": distance, "parents": {current_node}}
                # print("- update distance: ", child, distance)
                queue.append((distance, child))
            elif distance == distances[child]["distance"]:
                distances[child]["parents"].add(current_node)
                # print("- add distance: ", child, distance)
                queue.append((distance, child))

    # print("---")
    sink = root
    while sink.edges:
        sink = sink.edges[0][0]
    # print("---")
    # print(sink)
    # for d in distances:
    #     print(d, distances[d])

    new = {}
    nodes = [sink]
    visited = set()
    while nodes:
        node = nodes.pop(0)
        if node in visited:
            continue
        visited.add(node)
        # print("current:", node)
        if node not in new:
            # print(" added node:", node)
            new[node] = _Node(node.row, node.col, node.length)
        for p in distances[node]["parents"]:
            if p not in new:
                new[p] = _Node(p.row, p.col, p.length)
            #     print(" added node:", p)
            # print("  ", p.edges)
            for c, v in p.edges:
                if c == node:
                    new[p].edges.append((new[c], deepcopy(v)))
                    # print("  - edge:", p, c)
        extend = set(distances[node]["parents"]) - visited
        nodes.extend(extend)
    #     print("extend with:", extend)
    #     print(nodes)
    #     print()
    # print(new)
    for n in new:
        if n.row == 0 and n.col == 0:
            new_root = new[n]
    # print(new_root)

    return new_root


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

        if t > 1:
            print("longer than expected")
            break


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
    visited = {root}
    queue = [root]
    redirect = {}
    while queue:
        node = queue.pop(0)
        print("visit:", node)
        empty_children = [i for i, (child, variant) in enumerate(node.edges) if not variant]
        # node.edges = [(succ, variant) in node.edges for succ, variant in node.edges if not variant]
        print(empty_children)
        for succ, variant in node.edges:
            if not variant:
                print("empty", node, succ)
                # redirect[succ] = node
                # node.edges.extend(succ.edges)
            if succ not in visited:
                visited.add(succ)
                queue.append(succ)


def main(reference, obs):
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    # CTAACG TTACC - with equals inside
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
    dot_reduced = to_dot(reference, reduced_root, extra="r", cluster="cluster_1")
    d = "digraph {\n    " + "\n    " + dot_raw + "\n" + dot_reduced + "}"
    src = graphviz.Source(d)
    # src.view()
    rm_equals(reduced_root)
    print(to_dot(reference, reduced_root))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extractor sandbox")
    commands = parser.add_subparsers(dest="command", required=True, help="Commands")

    main_parser = commands.add_parser("main")
    main_parser.add_argument("reference", help="reference sequence")
    main_parser.add_argument("observed", help="observed sequence / [variants]")

    check_parser = commands.add_parser("check")
    check_parser.add_argument("max_length", help="max_length", type=int)
    check_parser.add_argument("min_length", help="min_length", type=int)
    check_parser.add_argument("p", type=float)
    check_parser.add_argument("mu_deletion", help="mu_deletion", type=float)
    check_parser.add_argument("mu_insertion", help="mu_insertion", type=float)

    check_set_parser = commands.add_parser("check_set")

    args = parser.parse_args()

    if args.command == "main":
        main(args.reference, args.observed)
    elif args.command == "check":
        check(
            args.max_length,
            args.min_length,
            args.p,
            args.mu_deletion,
            args.mu_insertion,
        )
    elif args.command == "check_set":
        check_set()

import sys
from copy import deepcopy

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.utils import random_sequence, random_variants
from algebra.variants import Parser, Variant, patch, to_hgvs
from mutalyzer.description_extractor import description_extractor


def reduce(root):
    distances = {root: {"distance": 0, "parents": set()}}

    queue = [(0, root)]
    visited = {}
    while len(queue) > 0:
        current_distance, current_node = queue.pop(0)
        # print(current_node, len(queue))

        if current_node in visited:
            continue

        if current_distance > distances[current_node]["distance"]:
            continue

        for child, edge in current_node.edges:
            if edge:
                distance = current_distance + 1
            else:
                distance = current_distance

            if distances.get(child) is None or distance < distances[child]["distance"]:
                distances[child] = {"distance": distance, "parents": {current_node}}
                queue.append((distance, child))
            elif distance == distances[child]["distance"]:
                distances[child]["parents"].add(current_node)
                queue.append((distance, child))
    # print("---")
    sink = root
    while sink.edges:
        sink = sink.edges[0][0]
    # print("---")
    # print(sink)
    # for d in distances:
    #     print(d, distances[d])
        # new_nodes_set.update(distances[d]["parents"])

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


def check():
    i = 0
    while True:
        reference = random_sequence(30, 30)
        variants = list(random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7))
        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)
        v_n = set([len(x) for x in traversal(reduced_root)])
        print(v_n)

        if len(v_n) > 1:
            print("error")
            break


def main():
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    if len(sys.argv) < 3:
        reference = random_sequence(100, 100)
        variants = list(
                    random_variants(reference, p=0.02, mu_deletion=1.5,
                                     mu_insertion=1.7)
        )
        observed = patch(reference, variants)
    else:
        reference = sys.argv[1]
        if "[" in sys.argv[2]:
            variants = Parser(sys.argv[2]).hgvs()
            observed = patch(reference, variants)
        else:
            observed = sys.argv[2]

    distance, lcs_nodes = edit(reference, observed)

    root, _ = lcs_graph(reference, observed, lcs_nodes)

    print("----")
    open("raw.dot", "w").write(to_dot(reference, root))
    print("----")
    reduced_root = reduce(root)
    open("reduced.dot", "w").write(to_dot(reference, reduced_root))
    print("----")
    # print(set([len(x) for x in traversal(root)]))
    print(set([len(x) for x in traversal(reduced_root)]))


if __name__ == "__main__":
    main()
    # check()


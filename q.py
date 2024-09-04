import json
import random
import subprocess
import sys

from algebra import LCSgraph
from algebra.extractor import local_supremal
from algebra.extractor.extractor import canonical
from algebra.utils import random_sequence, to_dot

from algebra.variants import Variant


from itertools import chain

from functools import reduce

import networkx as nx


def post_dominators(G, source):
    def intersect(u, v):
        while u != v:
            while nodes[u]["dfn"] < nodes[v]["dfn"]:
                u = nodes[u]["ipdom"]
            while nodes[u]["dfn"] > nodes[v]["dfn"]:
                v = nodes[v]["ipdom"]
        return u

    nodes = {}
    for idx, u in enumerate(nx.dfs_postorder_nodes(G, source)):
        if not G.succ[u]:
            nodes[u] = {"ipdom": u, "dfn": len(G.nodes()) - idx - 1}
        else:
            nodes[u] = {"ipdom": reduce(intersect, list(G.succ[u])), "dfn": len(G.nodes()) - idx - 1}

    return nodes



def main():
    debug = False

    if len(sys.argv) > 1:
        reference = sys.argv[1]
        observed = sys.argv[2]
        debug = True
    else:
        reference = random_sequence(32, 8)
        observed = random_sequence(32, 8)

    print(reference, observed)

    stdout = subprocess.run(["./a.out", reference, observed], check=True, stdout=subprocess.PIPE).stdout.decode("utf-8")
    try:
        cgraph = json.loads(stdout)
    except json.decoder.JSONDecodeError as exc:
        print(f'\"{stdout}\"')
        raise exc

    print(cgraph)
    # print(cgraph["nodes"])
    # print(cgraph["edges"])

    DG = nx.MultiDiGraph()
    DG.add_nodes_from([(k, cgraph["nodes"][k]) for k in cgraph["nodes"]])
    DG.add_edges_from([(e['head'], e['tail'], {'variant': e['variant'], 'weight': 0 if e['variant'] == 'lambda' else 1}) for e in cgraph["edges"]])

    # print(DG.nodes)
    # print(DG.nodes["s0"])
    #
    # print(DG.edges)
    # print(DG.edges[('s1', 's2')])

    # print(list(nx.all_shortest_paths(DG, "s7", "s0")))
    all_paths = nx.all_shortest_paths(DG, cgraph["source"], "s0", weight="weight")
    # print("all paths", list(all_paths))
    # targets = []
    targets2 = []
    for x in set(chain(*all_paths)):
        # targets.append(LCSgraph.Node(DG.nodes[x]["row"], DG.nodes[x]["col"], DG.nodes[x]["length"]))
        targets2.append(x)

    # print("targets to keep", targets)

    to_remove2 = []
    for n in DG.nodes:
        if n not in targets2:
            to_remove2.append(n)
    for n in to_remove2:
        DG.remove_node(n)

    graph = LCSgraph.from_sequence(reference, observed)
    # print("local supremal", local_supremal(reference, graph))

    # print(list(graph.nodes()))
    # x = graph.nodes()
    # for y in x:
    #     print("edges", y.edges)
    #
    #     to_remove = []
    #     for z in y.edges:
    #         print("z", z)
    #         if z[0] not in targets:
    #             print("remove")
    #             to_remove.append(z)
    #         print("z", z)
    #     for q in to_remove:
    #         y.edges.remove(q)
    #     print("edges after", y.edges)


    # print("remaining nodes", list(DG.nodes()))
    # print("remaining edges", list(DG.edges()))

    # print(nx.dominance_frontiers(DG, cgraph["source"]))
    # print(nx.dominance_frontiers(DG.reverse(), "s0"))
    # print(nx.immediate_dominators(DG, cgraph["source"]))
    # print(nx.immediate_dominators(DG.reverse(), "s0"))
    i = nx.immediate_dominators(DG.reverse(), "s0")
    print(i)

    lhs = cgraph["source"]
    variant = []
    while lhs != i[lhs]:
        rhs = i[lhs]
        print(lhs, DG.nodes[lhs], rhs, DG.nodes[rhs])
        if (DG.nodes[lhs]["row"] + DG.nodes[lhs]["length"] == DG.nodes[rhs]["row"] + DG.nodes[rhs]["length"] and
                DG.nodes[lhs]["col"] + DG.nodes[lhs]["length"] == DG.nodes[rhs]["col"] + DG.nodes[rhs]["length"]):
            print("lambda?")
        else:
            print(DG.edges(lhs, data="variant"))
            del_start = min([int(e[2].split(":")[0]) for e in filter(lambda e: e[2] != "lambda", DG.edges(lhs, data="variant"))])
            del_end = max([int(e[2].split(":")[1].split("/")[0]) for e in filter(lambda e: e[2] != "lambda", DG.in_edges(rhs, data="variant"))])
            # print("edge min", min([int(e[2].split(":")[0]) for e in chain(DG.edges(lhs, data="variant"), DG.edges(rhs, data="variant"))]))

            # del_end = max([int(e[2].split(":")[1].split("/")[0]) for e in filter(lambda e: e[2] != "lambda", chain(DG.edges(lhs, data="variant"), DG.edges(rhs, data="variant")))])

            # print("rhs edge min", min([int(e[2].split(":")[0]) for e in DG.edges(rhs, data="variant")]))
            # print("rhs edge min", min([int(e[2].split(":")[0].split("/")[0]) for e in DG.edges(rhs, data="variant")]))
            # del_start = DG.nodes[lhs]["row"] + DG.nodes[lhs]["length"]
            # del_start = min(DG.nodes[lhs]["row"] + DG.nodes[lhs]["length"], DG.nodes[rhs]["row"])
            # del_end = DG.nodes[rhs]["row"]
            # del_end = max(DG.nodes[lhs]["row"] + DG.nodes[lhs]["length"], DG.nodes[rhs]["row"])
            # ins_start = DG.nodes[lhs]["col"] + DG.nodes[lhs]["length"]
            # ins_start = DG.nodes[lhs]["col"]
            ins_start = DG.nodes[lhs]["col"] + del_start - DG.nodes[lhs]["row"]
            ins_end = DG.nodes[rhs]["col"] + del_end - DG.nodes[rhs]["row"]
            # if lhs == "s2" and rhs == "s1":
            #     # observed[lca.col + start - lca.row - shift:node.col + end - node.row - shift])
            #     #           21        28       27
            #     # ins_start = 22
            #     # ins_end = 24
            #     print(DG.nodes[lhs]["col"], del_start, DG.nodes[lhs]["row"])
            #     ins_start = DG.nodes[lhs]["col"] + del_start - DG.nodes[lhs]["row"]
            #     #                 21                 28        27
            #     ins_end = 24
            print("del", del_start, del_end)
            print("ins", ins_start, ins_end)
            v = Variant(min(del_start, del_end), max(del_start, del_end), observed[min(ins_start, ins_end): max(ins_start, ins_end)])
            print(v)
            variant.append(v)
        lhs = rhs
    print(variant)

    # G = nx.MultiDiGraph([(5, 2), (5, 4), (2, 3), (2, 0), (2, 1), (4, 2), (3, 0), (3, 0), (1, 0)])
    # for k, v in post_dominators(G, 5).items():
    # for k, v in post_dominators(DG, cgraph["source"]).items():
    #     print(k, v)

    print(canonical(graph))

    assert variant == canonical(graph)


if __name__ == "__main__":
    main()

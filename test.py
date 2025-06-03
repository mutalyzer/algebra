import json
# import random
import signal
import subprocess
import sys

from itertools import chain
#import networkx as nx

from algebra import LCSgraph
from algebra.extractor import local_supremal
from algebra.extractor.extractor import canonical
from algebra.utils import random_sequence, to_dot
from algebra.variants import Variant


def handler(signum, frame):
    raise KeyboardInterrupt


def nx_extract(reference, observed):
    # Create graph in C
    stdout = subprocess.run(["./a.out", reference, observed, "true"], check=True, stdout=subprocess.PIPE).stdout.decode("utf-8")
    try:
        cgraph = json.loads(stdout)
    except json.decoder.JSONDecodeError as exc:
        print(f'\"{stdout}\"')
        raise exc

    # Set source and guess sink
    source = cgraph["source"]
    i = iter(cgraph["nodes"])
    while True:
        try:
            sink = next(i)
        except StopIteration:
            break
        if sink != source:
            break

    # Convert C string variant into Variant object
    def str2var(desc):
        if desc == "lambda":
            return None

        start, remainder = desc.split(":")
        end, inserted = remainder.split("/")
        return Variant(int(start), int(end), inserted)

    # Convert C graph to NetworkX
    mdg = nx.MultiDiGraph()
    mdg.add_nodes_from([(k, cgraph["nodes"][k]) for k in cgraph["nodes"]])
    mdg.add_edges_from([(e['head'], e['tail'], {'variant': str2var(e['variant']), 'weight': 0 if str(e['variant']) == 'lambda' else 1}) for e in cgraph["edges"]])

    # Find all shortest paths and create a new subgraph with only those nodes
    all_paths = nx.all_shortest_paths(mdg, source, sink, weight="weight")
    mdg = mdg.subgraph(set(chain(*all_paths)))

    # Calculate the dominators on the reversed graph
    postdoms = nx.immediate_dominators(mdg.reverse(), sink)

    # Extract local supremals from reduced graph
    local_supremal = []
    lhs = source
    while lhs != postdoms[lhs]:
        rhs = postdoms[lhs]
        if not (mdg.nodes[lhs]["row"] + mdg.nodes[lhs]["length"] == mdg.nodes[rhs]["row"] + mdg.nodes[rhs]["length"] and
                mdg.nodes[lhs]["col"] + mdg.nodes[lhs]["length"] == mdg.nodes[rhs]["col"] + mdg.nodes[rhs]["length"]):
            del_start = min([e[2].start for e in filter(lambda e: e[2], mdg.out_edges(lhs, data="variant"))])
            del_end = max([e[2].end for e in filter(lambda e: e[2], mdg.in_edges(rhs, data="variant"))])
            ins_start = mdg.nodes[lhs]["col"] + del_start - mdg.nodes[lhs]["row"]
            ins_end = mdg.nodes[rhs]["col"] + del_end - mdg.nodes[rhs]["row"]
            v = Variant(del_start, del_end, observed[ins_start: ins_end])
            local_supremal.append(v)
        lhs = rhs

    # Validate against original Python implementation
    for lhs, rhs in zip(reversed(cgraph["canonical"]), local_supremal):
        assert lhs == str(rhs), f'canonical elements differ: {lhs} vs {rhs}'


def check(reference, observed, debug=False, timeout=2):
    #graph = LCSgraph.from_sequence(reference, observed)
    graph = LCSgraph(reference, observed)
    valgrind = []
    if debug:
        print("\n".join(to_dot(reference, graph, labels=False, hgvs=False)))
        print(canonical(graph))
        valgrind = ["valgrind", "--leak-check=full", "--error-exitcode=1"]

    stdout = subprocess.run([*valgrind, "./a.out", reference, observed], check=True, stderr=subprocess.DEVNULL, stdout=subprocess.PIPE, timeout=timeout).stdout.decode("utf-8")

    try:
        cgraph = json.loads(stdout)
    except json.decoder.JSONDecodeError as exc:
        print(f'\"{stdout}\"')
        raise exc

    #assert str(graph.supremal) == cgraph["supremal"], f'supremals differ: {cgraph["supremal"]} vs {graph.supremal}'

    nodes = list(graph.nodes())
    edges = list(graph.bfs_traversal())

    assert len(cgraph["nodes"]) == len(nodes), f'# cgraph nodes {len(cgraph["nodes"])} vs # graph nodes {len(nodes)}'
    assert len(cgraph["edges"]) == len(edges), f'# cgraph edges {len(cgraph["edges"])} vs # graph edges {len(edges)}'

    for cnode in cgraph["nodes"].values():
        try:
            cnode["node"] = next(node for node in nodes if node.row == cnode["row"] and node.col == cnode["col"])
        except StopIteration:
            raise ValueError(f'{cnode} not in graph nodes') from None

    for cedge in cgraph["edges"]:
        head = cgraph["nodes"][cedge["head"]]["node"]
        tail = cgraph["nodes"][cedge["tail"]]["node"]
        try:
            next((t, v) for t, v in head.edges if t == tail and str(v) == cedge["variant"])
        except StopIteration:
            raise ValueError(f'{cedge} not in edges of {head}') from None

    #assert len(cgraph["local_supremal"]) == len(local_supremal(reference, graph)), f'local supremal length'

    #for lhs, rhs in zip(cgraph["local_supremal"], local_supremal(reference, graph)):
    #    assert lhs == str(rhs), f'local supremal elements differ: {lhs} vs {rhs}'

    #if debug:
    #    print(list(reversed(cgraph["canonical"])), canonical(graph))


    #print()
    #print(cgraph["canonical"])
    #print(canonical(graph))
    #assert len(cgraph["canonical"]) == len(canonical(graph)), f'canonical length'

    #for lhs, rhs in zip(reversed(cgraph["canonical"]), canonical(graph)):
    #    assert lhs == str(rhs), f'canonical elements differ: {lhs} vs {rhs}'


def main():
    if len(sys.argv) == 3:
        print(sys.argv[1], sys.argv[2])
        return check(sys.argv[1], sys.argv[2], debug=False)

    signal.signal(signal.SIGINT, handler)

    failed = []
    n = 0
    while True:
        prefix = ""
        suffix = ""
        #if random.random() > 0.5:
        #    prefix = random_sequence(1000)
        #    suffix = random_sequence(1000)
        reference = prefix + random_sequence(142) + suffix
        observed = prefix + random_sequence(142) + suffix

        if n % 1_000 == 0:
            print(n)

        try:
            check(reference, observed)
            # nx_extract(reference, observed)
            n += 1
        except KeyboardInterrupt:
            break
        except Exception as exc:
            failed.append({"reference": reference, "observed": observed, "exception": exc})
            break

    print()
    print(n, failed)


if __name__ == "__main__":
    main()

import json
import subprocess
import sys

from itertools import chain
import networkx as nx

from algebra import LCSgraph
from algebra.extractor.extractor import canonical
from algebra.utils import random_sequence
from algebra.variants import Variant


def main():

    if len(sys.argv) > 1:
        reference = sys.argv[1]
        observed = sys.argv[2]
    else:
        reference = random_sequence(32, 0)
        observed = random_sequence(32, 0)
    print(reference, observed)

    # Create graph in C
    stdout = subprocess.run(["./a.out", reference, observed], check=True, stdout=subprocess.PIPE).stdout.decode("utf-8")
    try:
        cgraph = json.loads(stdout)
    except json.decoder.JSONDecodeError as exc:
        print(f'\"{stdout}\"')
        raise exc

    # Set source and guess sink
    source = cgraph["source"]
    i = iter(cgraph["nodes"])
    while True:
        sink = next(i)
        if sink != source:
            break

    # Convert C graph to NetworkX
    mdg = nx.MultiDiGraph()
    mdg.add_nodes_from([(k, cgraph["nodes"][k]) for k in cgraph["nodes"]])
    mdg.add_edges_from([(e['head'], e['tail'], {'variant': e['variant'], 'weight': 0 if e['variant'] == 'lambda' else 1}) for e in cgraph["edges"]])

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
            del_start = min([int(e[2].split(":")[0]) for e in filter(lambda e: e[2] != "lambda", mdg.edges(lhs, data="variant"))])
            del_end = max([int(e[2].split(":")[1].split("/")[0]) for e in filter(lambda e: e[2] != "lambda", mdg.in_edges(rhs, data="variant"))])
            ins_start = mdg.nodes[lhs]["col"] + del_start - mdg.nodes[lhs]["row"]
            ins_end = mdg.nodes[rhs]["col"] + del_end - mdg.nodes[rhs]["row"]
            v = Variant(del_start, del_end, observed[ins_start: ins_end])
            local_supremal.append(v)
        lhs = rhs
    print(local_supremal)

    # Validate against original Python implementation
    graph = LCSgraph.from_sequence(reference, observed)
    assert local_supremal == canonical(graph)


if __name__ == "__main__":
    main()

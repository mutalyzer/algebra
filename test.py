import json
import random
import signal
import subprocess
import sys

from algebra import LCSgraph
from algebra.extractor import local_supremal
from algebra.extractor.extractor import canonical
from algebra.utils import random_sequence, to_dot


def handler(signum, frame):
    raise KeyboardInterrupt


def check(reference, observed, debug=False, timeout=1):
    graph = LCSgraph.from_sequence(reference, observed)
    valgrind = []
    if debug:
        print("\n".join(to_dot(reference, graph, labels=False, hgvs=False)))
        print(canonical(graph))
        valgrind = ["valgrind", "--leak-check=full", "--error-exitcode=1"]

    stdout = subprocess.run([*valgrind, "./a.out", reference, observed], check=True, stdout=subprocess.PIPE, timeout=timeout).stdout.decode("utf-8")

    try:
        cgraph = json.loads(stdout)
    except json.decoder.JSONDecodeError as exc:
        print(f'\"{stdout}\"')
        raise exc

    assert str(graph.supremal) == cgraph["supremal"], f'supremals differ: {cgraph["supremal"]} vs {graph.supremal}'

    nodes = list(graph.nodes())
    edges = list(graph.bfs_traversal())

    assert len(cgraph["nodes"]) == len(nodes), f'# cgraph nodes {len(cgraph["nodes"])} vs # graph nodes {len(nodes)}'
    assert len(cgraph["edges"]) == len(edges), f'# cgraph edges {len(cgraph["edges"])} vs # graph edges {len(edges)}'

    for cnode in cgraph["nodes"].values():
        try:
            cnode["node"] = next(node for node in nodes if node == LCSgraph.Node(cnode["row"], cnode["col"], cnode["length"]))
        except StopIteration:
            raise ValueError(f'{cnode} not in graph nodes') from None

    for cedge in cgraph["edges"]:
        head = cgraph["nodes"][cedge["head"]]["node"]
        tail = cgraph["nodes"][cedge["tail"]]["node"]
        try:
            next((t, v) for t, v in head.edges if t == tail and str(v) == cedge["variant"])
        except StopIteration:
            raise ValueError(f'{cedge} not in edges of {head}') from None

    assert len(cgraph["local_supremal"]) == len(local_supremal(reference, graph)), f'local supremal length'

    for lhs, rhs in zip(cgraph["local_supremal"], local_supremal(reference, graph)):
        assert lhs == str(rhs), f'local supremal elements differ: {lhs} vs {rhs}'

    if debug:
        print(list(reversed(cgraph["canonical"])), canonical(graph))

    assert len(cgraph["canonical"]) == len(canonical(graph)), f'canonical length'

    for lhs, rhs in zip(reversed(cgraph["canonical"]), canonical(graph)):
        assert lhs == str(rhs), f'canonical elements differ: {lhs} vs {rhs}'


def main():
    if len(sys.argv) == 3:
        return check(sys.argv[1], sys.argv[2], debug=True)

    signal.signal(signal.SIGINT, handler)

    failed = []
    n = 0
    while True:
        prefix = ""
        suffix = ""
        #if random.random() > 0.5:
        #    prefix = random_sequence(1000)
        #    suffix = random_sequence(1000)
        reference = prefix + random_sequence(30) + suffix
        observed = prefix + random_sequence(30) + suffix

        try:
            print(reference, observed)
            check(reference, observed)
            n += 1
        except KeyboardInterrupt:
            break;
        except Exception as exc:
            failed.append({"reference": reference, "observed": observed, "exception": exc})
            break

    print()
    print(n, failed)


if __name__ == "__main__":
    main()

from collections import deque

from ..lcs.all_lcs_dfa import edit
from ..lcs.all_lcs_dfa import lcs_graph_dfa as lcs_graph
from ..variants import Variant


def canonical(observed, root):
    def lowest_common_ancestor(n_a, e_a, n_b, e_b):
        while n_a:
            n = n_b
            while n:
                if n == n_a:
                    return n_a, e_a, e_b
                n, e_b = visited[n]
            n_a, e_a = visited[n_a]

    queue = deque([(root, None, [])])

    visited = {}

    distances = {root: 0}

    multilevel_delins = set()

    while queue:

        node, parent, edge = queue.popleft()

        if not node.edges:
            sink = node

        if node in visited:
            other_parent, other_edge = visited[node]

            if (
                other_parent == parent
                or node in multilevel_delins
                or (parent.row - parent.col != other_parent.row - other_parent.col)
            ):
                lca, edge_a, edge_b = lowest_common_ancestor(other_parent, other_edge, parent, edge)

                start = min(edge_a.start, edge_b.start)
                start_offset = start - lca.row

                end = max(node.row - 1, edge.end, other_edge.end)
                end_offset = end - node.row

                delins = Variant(start, end, observed[lca.col + start_offset : node.col + end_offset])

                visited[node] = [lca, delins]
                if parent.row - parent.col != other_parent.row - other_parent.col:
                    multilevel_delins.add(node)

            else:
                lca, edge_a, edge_b = lowest_common_ancestor(*visited[other_parent], *visited[parent])

                start = min(edge_a.start, edge_b.start)
                start_offset = start - lca.row

                end = max(parent.row - 1, visited[parent][1].end, visited[other_parent][1].end)
                end_offset = end - parent.row

                delins = Variant(start, end, observed[lca.col + start_offset : parent.col + end_offset])

                visited[other_parent] = [lca, delins]
                visited[parent] = [lca, delins]

            continue

        else:
            visited[node] = [parent, edge]

        for succ_node, succ_edge in reversed(node.edges):
            if (
                succ_node not in distances
                or distances[succ_node] == distances[node] + 1
            ):
                distances[succ_node] = distances[node] + 1
                queue.append((succ_node, node, succ_edge))

    variants = []
    while sink:
        sink, variant = visited[sink]
        if variant:
            variants.append(variant)

    return reversed(variants)


def extract(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)

    return canonical(observed, root)

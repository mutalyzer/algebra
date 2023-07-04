from collections import deque

from ..lcs.all_lcs import edit
from ..lcs.all_lcs_dfa import lcs_graph_dfa as lcs_graph
from ..variants import Variant


def canonical(observed, root):
    def lowest_common_ancestor(n_a, e_a, n_b, e_b):
        while n_a:
            n = n_b
            while n:
                if n == n_a:
                    return n_a, e_a, e_b
                n, e_b, _ = visited[n]
            n_a, e_a, _ = visited[n_a]

    queue = deque([(root, None, None, 0)])

    visited = {}

    while queue:

        node, parent, edge, distance = queue.popleft()

        if not node.edges:
            sink = node

        if node in visited:
            other_parent, other_edge, other_distance = visited[node]

            if distance > other_distance:
                continue

            lca, edge_a, edge_b = lowest_common_ancestor(other_parent, other_edge, parent, edge)

            start = min(edge_a.start, edge_b.start)
            start_offset = start - lca.row

            if other_parent != parent and parent.row == other_parent.row and parent.col == other_parent.col:
                end = max(parent.row - 1, visited[parent][1].end, visited[other_parent][1].end)
                end_offset = end - parent.row

                delins = Variant(start, end, observed[lca.col + start_offset : parent.col + end_offset])

                visited[other_parent] = (lca, delins, distance - 1)
                # does not seem to be required to: visited[parent] = (lca, delins, distance - 1)

            else:
                end = max(node.row - 1, edge.end, other_edge.end)
                end_offset = end - node.row

                delins = Variant(start, end, observed[lca.col + start_offset : node.col + end_offset])

                visited[node] = (lca, delins, distance)


        else:
            visited[node] = (parent, edge, distance)

            for succ_node, succ_edge in node.edges:
                queue.append((succ_node, node, succ_edge, distance + 1))

    variants = []
    while sink:
        sink, variant, _ = visited[sink]
        if variant:
            variants.append(variant)

    return reversed(variants)


def extract(reference, observed):
    _, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)

    return canonical(observed, root)

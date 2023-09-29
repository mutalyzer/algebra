from .. import Variant


def find_dominators(node, outgoing, incoming=0, depth=0, visited=None):
    """Generate the (post) dominator nodes.

    Parameters
    ----------
    node : `_Node` (opaque data type)
        The current working node (initially the root of the LCS graph).
    outgoing: int
        The minimum start position of the node variants.
    incoming: int
        The maximum end position of the variants pointing towards the node.

    Other Parameters
    ----------------
    depth: int
        The node BFS traversal level.
    visited: dict
        Already visited nodes.

    Yields
    ------
    tuple
        (dominator node, maximum incoming position, minimum outgoing position)
    """

    if not visited:
        visited = {}

    if node in visited:
        visited[node]["depth"] = min(visited[node]["depth"], depth)
        visited[node]["incoming"] = max(visited[node]["incoming"], incoming)
        return

    visited[node] = {
        "pdom": {node},
        "depth": depth,
        "incoming": incoming,
        "outgoing": outgoing,
    }

    pdom = set()
    for child, edge in node.edges:
        yield from find_dominators(child, outgoing, edge.end, depth + 1, visited)
        if not pdom:
            pdom = visited[child]["pdom"]
        pdom = pdom.intersection(visited[child]["pdom"])
        visited[node]["outgoing"] = min(visited[node]["outgoing"], edge.start)
    visited[node]["pdom"] = pdom.union(visited[node]["pdom"])

    if not depth:
        for dominator in sorted(visited[node]["pdom"], key=lambda x: visited[x]["depth"]):
            yield dominator, visited[dominator]["incoming"], visited[dominator]["outgoing"]

def local_suprema(reference, observed, root):
    """Generate the local supremal variants.

    Parameters
    ----------
    reference : str
        The reference sequence.
    observed : str
        The observed sequence.
    root : `_Node` (opaque data type)
        The root of the LCS graph.

    Yields
    ------
    Variant
        The local supremal variant.
    """

    post_dominators = iter(find_dominators(root, len(reference)))
    parent = next(post_dominators)

    for child in post_dominators:
        yield Variant(
            parent[2],
            child[1],
            observed[parent[0].col + parent[2] - parent[0].row:child[0].col + child[1] - child[0].row]
        )
        parent = child

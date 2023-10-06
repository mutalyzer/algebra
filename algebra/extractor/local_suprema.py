from .. import Variant


def find_dominators(root, start, end=0, node = None, visited=None):
    """Generate the (post) dominator nodes.

    Parameters
    ----------
    root : `_Node` (opaque data type)
        The root of the LCS graph.
    start: int
        The minimum start position of the node variants.
    end: int
        The maximum end position of the variants entering the node.
    node : `_Node` (opaque data type)
        The current node.

    Other Parameters
    ----------------
    visited: dict
        Already visited nodes.

    Yields
    ------
    tuple
        (dominator node, maximum start position (in), minimum end position (out))
    """

    if node is None:
        node = root

    if not visited:
        visited = {}

    if node in visited:
        visited[node]["start"] = max(visited[node]["start"], end)
        return

    visited[node] = {
        "pdom": {node},
        "start": end,
        "end": start,
    }

    pdom = set()
    for child, edge in node.edges:
        yield from find_dominators(root, start, edge.end, child, visited)
        if not pdom:
            pdom = visited[child]["pdom"]
        pdom = pdom.intersection(visited[child]["pdom"])
        visited[node]["end"] = min(visited[node]["end"], edge.start)
    visited[node]["pdom"] = pdom.union(visited[node]["pdom"])

    if node == root:
        for dominator in sorted(visited[node]["pdom"], key=lambda x: x.row):
            yield dominator, visited[dominator]["start"], visited[dominator]["end"]


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
            observed[parent[0].col + parent[2] - parent[0].row:
                     child[0].col + child[1] - child[0].row]
        )
        parent = child

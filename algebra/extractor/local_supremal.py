from operator import attrgetter

from .. import Variant


def local_supremal(reference, observed, root):
    """Extract the local supremal representation.

    Parameters
    ----------
    reference : str
        The reference sequence.
    observed : str
        The observed sequence.
    root : `_Node` (opaque data type)
        The root of the LCS graph.

    Returns
    -------
    list
        The list of variants (allele).

    See Also
    --------
    `algebra.lcs.lcs_graph` : Constructs the LCS graph.
    """

    def post_dominators(node, start, visited):
        if node in visited:
            visited[node]["start"] = max(visited[node]["start"], start)
            return visited

        visited[node] = {
            "post": {node},
            "start": start,
            "end": len(reference),
        }

        post = set()
        for child, variant in node.edges:
            post_dominators(child, variant.end, visited)
            if not post:
                post = visited[child]["post"]
            post = post.intersection(visited[child]["post"])

            visited[node]["end"] = min(visited[node]["end"], variant.start)

        visited[node]["post"] = post.union(visited[node]["post"])
        return visited

    visited = post_dominators(root, 0, {})

    variants = []
    pred = None
    for node in sorted(visited[root]["post"], key=attrgetter("row")):
        if pred:
            start = visited[pred]["end"]
            end = visited[node]["start"]
            variants.append(
                Variant(start, end, observed[start + pred.col - pred.row:end + node.col - node.row])
            )
        pred = node
    return variants

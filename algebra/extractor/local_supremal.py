"""Extract a local supremal variant representation. Where local supremal
is defined as the variants between matching symbols in all minimal
alignments."""


from .. import Variant


def local_supremal(reference, observed, graph):
    """Extract the local supremal representation.

    Parameters
    ----------
    reference : str
        The reference sequence.
    observed : str
        The observed sequence.
    graph : `_Node` (opaque data type)
        The LCS graph.

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

    shift = graph.row
    visited = post_dominators(graph, 0, {})
    variants = []
    parent = None
    for node in sorted(visited[graph]["post"]):
        if parent:
            start = visited[parent]["end"]
            end = visited[node]["start"]
            variants.append(
                Variant(start, end,
                        observed[parent.col + start - parent.row - shift:node.col + end - node.row - shift]))
        parent = node

    return variants

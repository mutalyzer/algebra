"""Extract canonical variants."""


from algebra_ext import LCSgraph, canonical


def extract(reference, variants):
    graph = LCSgraph.from_variants(reference, variants)
    return canonical(graph), graph

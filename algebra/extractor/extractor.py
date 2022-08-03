from copy import deepcopy

from algebra.lcs import edit, lcs_graph
from algebra.lcs.all_lcs import _Node
from algebra.variants import Variant


def _get_sink(root):
    sink = root
    while sink.edges:
        sink = sink.edges[0][0]
    return sink


def _extract_subgraph(distances, sink):
    new = {}
    nodes = [sink]
    visited = set()
    while nodes:
        node = nodes.pop(0)
        if node in visited:
            continue
        visited.add(node)
        if node not in new:
            new[node] = _Node(node.row, node.col, node.length)
        for p in distances[node]["parents"]:
            if p not in new:
                new[p] = _Node(p.row, p.col, p.length)
            for c, v in p.edges:
                if c == node:
                    new[p].edges.append((new[c], deepcopy(v)))
        extend = set(distances[node]["parents"]) - visited
        nodes.extend(extend)
    for n in new:
        if n.row == 0 and n.col == 0:
            new_root = new[n]
    return new_root


def reduce(root):
    """
    Reduce the graph to the paths that provide a minimal number of variants.
    """
    distances = {root: {"distance": 0, "parents": set()}}

    queue = [(0, root)]
    while len(queue) > 0:
        current_distance, current_node = queue.pop(0)
        if current_distance > distances[current_node]["distance"]:
            continue

        successors = {}
        for child, edge in current_node.edges:
            if child not in successors:
                successors[child] = current_distance
            if edge:
                successors[child] = current_distance + 1

        for child, distance in successors.items():
            if distances.get(child) is None or distance < distances[child]["distance"]:
                distances[child] = {"distance": distance, "parents": {current_node}}
                queue.append((distance, child))
            elif distance == distances[child]["distance"]:
                distances[child]["parents"].add(current_node)
                queue.append((distance, child))

    return _extract_subgraph(distances, _get_sink(root))


def rm_equals(root):
    """
    Remove the equal edges.
    """

    def _empty_child():
        visited = {root}
        queue = [root]
        while queue:
            node = queue.pop(0)
            for i, (child, variant) in enumerate(node.edges):
                if not variant:
                    return child, node, i
                if child not in visited:
                    visited.add(child)
                    queue.append(child)
        return None, None, None

    def _update_parent():
        parent.edges.pop(idx)
        parent.edges.extend(empty_child.edges)

    def _update_other_parents():
        visited = {root}
        queue = [root]
        parents = []
        while queue:
            node = queue.pop(0)
            for i, (child, variant) in enumerate(node.edges):
                if child == empty_child:
                    parents.append((node, i))
                if child not in visited:
                    visited.add(child)
                    queue.append(child)
        for (p, i) in parents:
            p.edges[i] = (parent, p.edges[i][1])

    while True:
        empty_child, parent, idx = _empty_child()
        if not empty_child:
            break
        _update_parent()
        _update_other_parents()


def _shortest_delins(start_node, end_node, reference):
    """
    Get the deletion insertion variant between the start and the end nodes.
    """
    next_node, edge = max(start_node.edges, key=lambda x: x[1][0].start)
    sequence = edge[0].sequence
    start = edge[0].start
    end = edge[0].end
    while next_node != end_node:
        next_node, edge = min(next_node.edges, key=lambda x: x[1][0].start)
        edge = edge[0]
        if end is not None:
            sequence += reference[end : edge.start]
        end = edge.end
        sequence += edge.sequence
    return Variant(start, end, sequence)


def _repeat(node, reference):
    """
    Determine if the multiple edges of a node form a full or partial repeat.
    """

    def get_rotations(s):
        t = s
        r = 1
        while True:
            t = t[-1] + t[:-1]
            r += 1
            if t == s:
                return r

    edges = sorted([edge[1][0] for edge in node.edges], key=lambda x: x.start)
    first_edge = min(edges, key=lambda x: x.start)
    last_edge = max(edges, key=lambda x: x.start)
    if last_edge.sequence == "":
        # deletion
        sequences = [reference[edge.start : edge.end] for edge in edges]
        rotations = get_rotations(sequences[-1])
        if rotations <= len(sequences):
            # repeat
            start = first_edge.start + sequences.index(sequences[-1])
            end = last_edge.end
            affected = reference[start:end]
            r_s = affected[: len(set(sequences))]
            r_n = (len(affected) - len(sequences[0])) // len(r_s)
            return f"{start+1}_{end}{r_s}[{r_n}]"
        else:
            # shift
            return last_edge.to_hgvs(reference)
    else:
        # insertion
        sequences = [edge.sequence for edge in edges]
        rotations = get_rotations(sequences[-1])
        if rotations <= len(sequences):
            # repeat
            start = first_edge.start + sequences.index(sequences[-1])
            end = last_edge.end
            affected = reference[start:end]
            r_s = affected[: len(set(sequences))]
            r_n = (len(affected) + len(sequences[0])) // len(r_s)
            if start + 1 == end:
                location = f"{start+1}"
            else:
                location = f"{start+1}_{end}"
            return f"{location}{r_s}[{r_n}]"
        else:
            # shift
            return last_edge.to_hgvs(reference)


def get_variants(root, reference):
    """
    Traverse the graph to extract the 'canonical/normalized' variants.

    1. A node has only one child:
       a. with one edge -> a simple the variant.
       b. with multpile edges -> a form of a repeat variant.
    2. A node has multiple children -> a complex structure (deletion insertion).
    """
    visited = {root}
    queue = [root]

    structures = {}
    c_open = None
    c_close = False
    descriptions = []

    while queue:
        node = queue.pop(0)
        if c_close:
            # complex structure (2)
            structures[c_open] = node
            descriptions.append(
                _shortest_delins(c_open, node, reference).to_hgvs(reference)
            )
            c_close = False
            c_open = None
        children = set()
        for child, variant in node.edges:
            children.add(child)
        if len(children) > 1 and c_open is None:
            # start complex structure (2)
            c_open = node
        if len(children) == 1 and c_open is None:
            if len(node.edges) > 1:
                # repeat (1.b.)
                descriptions.append(_repeat(node, reference))
            else:
                # simple variant(1.a.)
                descriptions.append(node.edges[0][1][0].to_hgvs(reference))
        for child in children:
            if child not in visited:
                visited.add(child)
                queue.append(child)
        if len(queue) == 1 and c_open:
            # end complex structure when next node is popped (2)
            c_close = True
    if len(descriptions) > 1:
        return "[" + ";".join(descriptions) + "]"
    elif len(descriptions) == 1:
        return descriptions[0]
    else:
        return "="


def extract(reference, observed):

    distance, lcs_nodes = edit(reference, observed)
    root, _ = lcs_graph(reference, observed, lcs_nodes)

    reduced_root = reduce(root)
    rm_equals(reduced_root)

    return get_variants(reduced_root, reference)

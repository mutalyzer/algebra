from collections import deque


def incoming_as_dot_number(root):
    node_count = 0
    root.incoming = 0
    visited = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for succ, variant in node.edges:
            if succ not in visited:
                node_count += 1
                visited.add(succ)
                succ.incoming = node_count
                queue.append(succ)



def add_pre_edges(root):
    root.pre_edges = []
    visited = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        if not node.edges:
            sink = node
        for succ, variant in node.edges:
            if succ not in visited:
                succ.pre_edges = []
                visited.add(succ)
                queue.append(succ)
            succ.pre_edges.append((node, variant))
    return sink


def limit_minimal_accordion(root, lim_left, lim_right):
    incoming_as_dot_number(root)
    sink = add_pre_edges(root)

    level_left = {root}
    visited_left = {}

    level_right = {sink}
    visited_right = {}

    next_level_left = set()
    next_level_right = set()

    for node in level_left:
        for succ, edge in node.edges:
            if (
                (edge.start == edge.end and edge.start >= lim_left)
                or (edge.end > edge.start > lim_left)
            ) and succ not in visited_left:
                visited_left[succ] = (node, edge)
                next_level_left.add(succ)
    if next_level_left.intersection(level_right):
        return

    for node in level_right:
        for succ, edge in node.pre_edges:
            if lim_right > edge.end and succ not in visited_right:
                visited_right[succ] = (node, edge)
                next_level_right.add(succ)
    if next_level_left.intersection(next_level_right):
        next_level_left.intersection(next_level_right)

    while next_level_left or next_level_right:
        level_left = next_level_left
        next_level_left = set()
        for node in level_left:
            for succ, edge in node.edges:
                if succ in visited_right:
                    variants = [edge]
                    tracking_node = node
                    while tracking_node in visited_left:
                        variants.insert(0, visited_left[tracking_node][1])
                        tracking_node = visited_left[tracking_node][0]
                    tracking_node = succ
                    while tracking_node in visited_right:
                        variants.append(visited_right[tracking_node][1])
                        tracking_node = visited_right[tracking_node][0]
                    return variants
                elif succ not in visited_left:
                    visited_left[succ] = (node, edge)
                    next_level_left.add(succ)

        level_right = next_level_right
        next_level_right = set()

        for node in level_right:
            for succ, edge in node.pre_edges:
                if succ in visited_left:
                    variants = [edge]
                    tracking_node = node
                    while tracking_node in visited_left:
                        variants.insert(0, visited_left[tracking_node][1])
                        tracking_node = visited_left[tracking_node][0]
                    tracking_node = succ
                    while tracking_node in visited_right:
                        variants.append(visited_right[tracking_node][1])
                        tracking_node = visited_right[tracking_node][0]
                    return variants
                elif succ not in visited_right:
                    visited_right[succ] = (node, edge)
                    next_level_right.add(succ)


def limit_minimal_dfs(root, lim_left, lim_right, first=True):
    """Generate all the paths (or just the first) within the reference limits.

    Parameters
    ----------
    root : `_Node` (opaque data type)
        The root of the LCS graph.
    lim_left : int
        The left most reference position not to be affected by the variants.
    lim_right : int
        The right most reference position not to be affected by the variants.

    Other Parameters
    ----------------
    first: bool
        The node BFS traversal level.


    Yields
    ------
    Variant
        The local supremal variant.
    """

    def traverse(node, path):
        if not node.edges:
            yield path

        for succ, edge in node.edges:
            if edge.start > lim_left and lim_right >= edge.end:
                yield from traverse(succ, path + [edge])

    for limited in traverse(root, []):
        yield limited
        if first:
            return

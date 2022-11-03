from copy import deepcopy


def to_hgvs(word, max_cover, inv, pmrs, overlap):
    def _update_nodes(node, pattern):
        print("  - add", node, pattern, "from", current_node, nodes[current_node])
        paths = deepcopy(nodes[current_node])
        for path in paths:
            if path[-1] == "":
                path[-1] = pattern
            elif "]" not in path[-1] and "]" not in pattern:
                path[-1] = pattern + path[-1]
            else:
                path.append(pattern)
        if node in nodes:
            nodes[node].extend(paths)
        else:
            nodes[node] = paths
        print("  - after", node, nodes[node])

    def _next(node):
        if node == -1:
            return
        for primitive in inv[node]:
            print("- repeat:", primitive)
            unit_length = pmrs[primitive][1]
            possible_next = current_node - unit_length * 2
            occurrences = 2
            print(" - possible_next:", possible_next)
            if possible_next >= -1:
                pattern = word[current_node - unit_length + 1: current_node + 1]
                _update_nodes(possible_next, f"{pattern}[{occurrences}]")
            else:
                while possible_next >= pmrs[primitive][0] and possible_next > 0:
                    if max_cover[possible_next] + occurrences * unit_length == max_cover[current_node]:
                        pattern = word[current_node - unit_length + 1 : current_node + 1]
                        _update_nodes(possible_next, f"{pattern}[{occurrences}]")
                    occurrences += 2
                    possible_next -= unit_length
                print("\n\n - possible_next:", possible_next)
        if current_node == 0:
            previous_max_cover = 0
        else:
            previous_max_cover = max_cover[current_node - 1]
        if max_cover[current_node] == previous_max_cover:
            _update_nodes(current_node - 1, word[current_node])

    print("----")
    nodes = {len(max_cover) - 1: [[""]]}

    while True:
        current_node = max(nodes.keys())
        print("\ncurrent_node:", current_node)
        print(nodes)
        _next(current_node)
        if len(nodes) == 1:
            print("\n\n---")
            paths = [";".join(reversed(path)) for path in nodes[current_node]]
            for path in paths:
                print(path)
            print("done")
            return paths
        else:
            nodes.pop(current_node)

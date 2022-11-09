from copy import deepcopy
import sys
from algebra.extractor.cover import find_pmrs, overlapping, inv_array, print_tables, brute_cover, cover




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


def to_hgvs_2(word, max_cover, inv, pmrs, overlap):
    def repeats(idx, sol=[], indent=""):
        print(indent, "idx:", idx, sol)
        print(indent, "----------")
        if idx == -1:
            yield sol
        else:
            for p in inv[idx]:
                print(indent, "pattern:", p, pmrs[p])
                if overlap[p] > 0:
                    end = overlap[p] + 1
                else:
                    end = pmrs[p][0] + 1
                for i in range(pmrs[p][0], end):
                    print(indent, "i:", i, "; idx:", idx, "; idx - i + 1:", idx - i + 1)
                    if (idx - i + 1) % pmrs[p][1] == 0:
                        if i == 0 or (i > 0 and idx - i + 1 + max_cover[i-1] == max_cover[idx]):
                            print(indent, "extend with:", (p, i, idx))
                            sol.append((p, i, idx))
                            yield from repeats(i-1, sol, indent + "  ")
                            sol.pop(-1)
                            print(indent, "back to idx:", idx, sol)

    print("\n\n====\n")
    for r in repeats(len(word)-1, []):
        print("\n==>", r, "\n")


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} word", file=sys.stderr)
        return

    word = sys.argv[1]
    n = len(word)
    print(n, word)

    pmrs = sorted(find_pmrs(word))
    overlap = overlapping(pmrs)
    for idx, pmr in enumerate(pmrs):
        print(f"{pmr},  # {idx:2}: {word[pmr[0]:pmr[0] + pmr[1]]} : {overlap[idx]}")

    inv = inv_array(n, pmrs)
    max_cover = cover(word, pmrs)
    print_tables(n, word, inv, max_cover)
    print(brute_cover(word, pmrs))
    to_hgvs_2(word, max_cover, inv, pmrs, overlap)


if __name__ == '__main__':
    main()

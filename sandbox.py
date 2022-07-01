import argparse
import json
import time
from copy import deepcopy
from multiprocessing import Pool

import graphviz
from mutalyzer.description_extractor import description_extractor
from mutalyzer.mutator import mutate
from mutalyzer.reference import get_reference_model
from mutalyzer.spdi_converter import spdi_converter, spdi_to_hgvs
from mutalyzer_hgvs_parser import to_model
from mutalyzer_mutator import mutate as mutate_raw

from algebra.lcs import edit, lcs_graph, to_dot
from algebra.lcs.all_lcs import _Node, traversal
from algebra.utils import random_sequence, random_variants
from algebra.variants import Parser, patch, to_hgvs


def get_sink(root):
    sink = root
    while sink.edges:
        sink = sink.edges[0][0]
    return sink


def extract_subgraph(distances, sink):
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

    return extract_subgraph(distances, get_sink(root))


def reduce_repeats(root):
    """Preferring the nodes with the most repeats."""
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
                successors[child] += 1

        for child, distance in successors.items():
            if distances.get(child) is None or distance > distances[child]["distance"]:
                distances[child] = {"distance": distance, "parents": {current_node}}
                queue.append((distance, child))
            elif distance == distances[child]["distance"]:
                distances[child]["parents"].add(current_node)
                queue.append((distance, child))

    return extract_subgraph(distances, get_sink(root))


def reduce_simple_repeats(root):
    def _merge_repeats(variants):
        # print(variants)
        # print("---")
        if len(variants) > 1:
            # print("more", [max(variants, key=lambda x: x.start)])
            return [max(variants, key=lambda x: x.start)]
        else:
            # print("one", variants)
            return variants
            # for variant in variants:
            #     print(variant.start, variant.end, variant.sequence)

    visited = {root}
    queue = [root]
    while queue:
        node = queue.pop(0)
        children = {}
        for child, variant in node.edges:
            if child not in children:
                children[child] = []
            children[child].append(variant[0])
        new_edges = [(child, _merge_repeats(children[child])) for child in children]
        node.edges = new_edges
        for child in children:
            if child not in visited:
                visited.add(child)
                queue.append(child)


def compare():
    def _convert(description):
        print(to_model(description, start_rule="variants"))

    INTERESTING = [
        "CAAA A",
        "GCG GGCGATTCTT",
        "CAAGA GTCAA",
        "CAAGA CAA",
        "ACCCC CCAC",
        "TCT TCCTCCT",
        "CTCTTG CTCTCTTG",
        "CACCT CAGGG",
        "CACCTAA CAGGGAA",
        "AAGTT TCAGT",
        "CGCCC AGCCT",
        "ATATC AATTC",
        "ACCGA TCCGT",
        "GCACT ACTAC",
        "GCACT ACTAC",
        "CATAA ATTCA",
        "ACTAC TACCG",
        "TCCAT CACCT",
        "GACCT GCGCC",
        "AACGACT AGGTACG",
    ]

    for t in INTERESTING:
        reference = t.split()[0]
        observed = t.split()[1]
        print("\n------")
        print(reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)
        new = [to_hgvs(path, reference) for path in traversal(reduced_root)]
        old = description_extractor(reference, observed)
        _convert(old)
        if len(new) > 1:
            print(" - Multiple options")
        else:
            print(" - One option:", new)
        if old not in new:
            print(" - NOT IN")
            print(to_dot(reference, reduced_root))
            print("  - new:", new)
            print("  - old:", old)


def middle_equal():
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    i = 0
    while True:
        reference = random_sequence(20, 20)
        variants = list(
            random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7)
        )

        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)

        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)

        for path in traversal(root):
            hgvs = [to_hgvs([v], reference) for v in path]
            for v in hgvs:
                if v == "=":
                    print(reference, observed)
                    return


def check(max_length, min_length, p, mu_deletion, mu_insertion):
    i = 0
    while True:
        reference = random_sequence(max_length, min_length)
        # variants = list(random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7))
        variants = list(
            random_variants(
                reference, p=p, mu_deletion=mu_deletion, mu_insertion=mu_insertion
            )
        )
        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        print("distance:", distance)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        t_s = time.time()
        reduced_root = reduce(root)
        t = time.time() - t_s
        v_n = set([len(x) for x in traversal(reduced_root)])
        print(v_n, t)

        if len(v_n) > 1:
            print("error")
            break

        # if t > 5:
        #     print("longer than expected")
        #     break


def check_set():
    status = {}
    variants = [
        "GCCC TTCGCGCCCCT",
        "CATCA ATAATCTCAA",
        "TCCTT CGCTCTATCCCT",
        "ATTTT AGGTTAGATTTGGATTGT",
        "CGGGAACTTA TTACGCCGGGCGCTTA",
        "TCTCAACGAG AGTTCACAAATCGAA",
        "AGGAACGACCCGAAAGATTCTGGGAAACGGAGCGTCTCTA ATCACCAGGTCCAGAAGGTTGACCAGGTAGTAGCCGGGATTCTGGATTTTTCCGATC",
    ]
    for variant in variants:
        reference = variant.split()[0]
        observed = variant.split()[1]
        print(reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        print("distance:", distance)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)
        v_n = set([len(x) for x in traversal(reduced_root)])
        print(v_n)

        if len(v_n) > 1:
            status[variant] = "error"
        else:
            status[variant] = "OK"
    for variant in status:
        print(variant)
        print(status[variant])


def rm_equals(root):
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
        # print("  - parent edges before:")
        # for edge in parent.edges:
        #     print("   ", edge)
        parent.edges.pop(idx)
        parent.edges.extend(empty_child.edges)
        # print("  - parent edges after:")
        # for edge in parent.edges:
        #     print("   ", edge)

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
        # print(" - other parents:", parents, "\n")
        for (p, i) in parents:
            # print("   ", p, i)
            p.edges[i] = (parent, p.edges[i][1])

    while True:
        empty_child, parent, idx = _empty_child()
        # print("- current:", empty_child, parent, idx)
        if not empty_child:
            break
        _update_parent()
        _update_other_parents()


def rm_equals_alt(root):

    visited = {root}
    queue = [root]
    redirect_children = {}
    while queue:
        node = queue.pop(0)
        for child, variant in node.edges:
            if not variant:
                redirect_children[child] = node
            if child not in visited:
                visited.add(child)
                queue.append(child)

    print("redirect_children")
    for equal_child in redirect_children:
        print(equal_child, redirect_children[equal_child])

    print("---")

    visited = {root}
    redirect_parents = {}
    queue = [root]
    while queue:
        node = queue.pop(0)
        # print("- node:", node)
        for i, (child, variant) in enumerate(node.edges):
            if child in redirect_children and node != redirect_children[child]:
                # print(" - child:", child)
                redirect_parents[child] = (node, i)
            if child not in visited:
                visited.add(child)
                queue.append(child)

    print("redirect_parents")
    for parent in redirect_parents:
        print(parent, redirect_parents[parent])

    # print(redirect)


def to_dot_repeats(reference, root, extra="", cluster=""):
    """The LCS graph in Graphviz DOT format with repeat edges merged."""

    def _merge_repeats(variants):
        # print(variants)
        # print("---")
        if len(variants) > 1:
            return "|".join([to_hgvs([v], reference) for v in variants])
        else:
            return to_hgvs(variants)
            # for variant in variants:
            #     print(variant.start, variant.end, variant.sequence)

    def traverse():
        # breadth-first traversal
        visited = {root}
        queue = [root]
        while queue:
            node = queue.pop(0)
            children = {}
            for child, variant in node.edges:
                if child not in children:
                    children[child] = []
                children[child].append(variant[0])
            for child in children:
                # print(_merge_repeats(children[child]))
                yield (
                    f'"{extra}{node.row}_{node.col}" -> "{extra}{child.row}_{child.col}"'
                    f' [label="{_merge_repeats(children[child])}"];'
                )
                if child not in visited:
                    visited.add(child)
                    queue.append(child)

    if cluster:
        return "subgraph " + cluster + " {\n    " + "\n    ".join(traverse()) + "\n}"
    return "digraph {\n    " + "\n    ".join(traverse()) + "\n}"


def extract_main(reference, obs):
    # CATATAGT "[4_5del;7delinsGA]"
    # CTAA TTA - with equals inside
    # CTAACG TTACC - with equals inside
    # CTTTG CTATTTT - with consecutive equals inside
    # GCCTT GCAGCCCAT - with consecutive equals inside
    # AGGTA AAGAAGGGGA - with consecutive equals inside
    # TTGTA TTTGTGTT - with consecutive equals inside
    # ACTAA ACGCCTATTAAATAAA - with consecutive equals inside
    # CAGGG AACTCAGGTAGGGTTAGAT - with three consecutive equals inside
    # CGACTGACGTTACCGAAGTTTTTTGTACAGTCGACTGACG CGACTGACATTACCGAAGTTTTTTTGTACAGGGTTCTGACG - complex component
    if "[" in obs:
        variants = Parser(obs).hgvs()
        observed = patch(reference, variants)
    else:
        observed = obs

    distance, lcs_nodes = edit(reference, observed)
    print("distance:", distance)

    root, _ = lcs_graph(reference, observed, lcs_nodes)

    print("----")
    open("raw.dot", "w").write(to_dot(reference, root))
    print("----")
    reduced_root = reduce(root)
    open("reduced.dot", "w").write(to_dot(reference, reduced_root))
    print("----")

    # print(set([len(x) for x in traversal(root)]))
    # print(set([len(x) for x in traversal(reduced_root)]))

    # for d in traversal(reduced_root):
    #     print(to_hgvs(d, reference))

    dot_raw = to_dot(reference, root, cluster="cluster_0")
    dot_reduced = to_dot(reference, reduced_root, extra="d", cluster="cluster_1")

    rm_equals(reduced_root)

    dot_no_equals = to_dot(reference, reduced_root, extra="n", cluster="cluster_2")

    dot_repeats = to_dot_repeats(
        reference, reduced_root, extra="r", cluster="cluster_3"
    )

    root_reduced_repeats = reduce_repeats(reduced_root)

    dot_reduced_repeats = to_dot_repeats(
        reference, root_reduced_repeats, extra="rr", cluster="cluster_4"
    )

    d = (
        "digraph {\n    "
        + "\n    "
        + dot_raw
        + "\n"
        + dot_reduced
        + "\n"
        + dot_no_equals
        + "\n"
        + dot_repeats
        + "\n"
        + dot_reduced_repeats
        + "}"
    )
    src = graphviz.Source(d)
    src.view()

    # print(to_dot_repeats(reference, reduced_root))

    # print(to_dot(reference, reduced_root))


def consecutive_equals(root):
    visited = {root}
    queue = [root]
    while queue:
        node = queue.pop(0)
        for child, c_variant in node.edges:
            if not c_variant:
                for grand_child, g_variant in child.edges:
                    if not g_variant:
                        for grand_grand_child, g_g_variant in grand_child.edges:
                            if not g_g_variant:
                                return True
            if child not in visited:
                visited.add(child)
                queue.append(child)
    return False


def get_consecutive_equals():
    i = 0
    while True:
        reference = random_sequence(5, 5)
        variants = list(
            random_variants(reference, p=0.32, mu_deletion=1.5, mu_insertion=1.7)
        )
        observed = patch(reference, variants)
        i += 1
        print(i, reference, observed)
        distance, lcs_nodes = edit(reference, observed)
        root, _ = lcs_graph(reference, observed, lcs_nodes)
        reduced_root = reduce(root)

        if consecutive_equals(reduced_root):
            print("found it")
            dot_raw = to_dot(reference, root, cluster="cluster_0")
            dot_reduced = to_dot(
                reference, reduced_root, extra="r", cluster="cluster_1"
            )
            d = "digraph {\n    " + "\n    " + dot_raw + "\n" + dot_reduced + "}"
            src = graphviz.Source(d)
            src.view()
            break


def extract(reference, observed):

    t_0 = time.time()
    distance, lcs_nodes = edit(reference, observed)
    t_1 = time.time()
    root, _ = lcs_graph(reference, observed, lcs_nodes)
    t_2 = time.time()
    reduced_root = reduce(root)
    t_3 = time.time()
    rm_equals(reduced_root)
    t_4 = time.time()
    reduce_simple_repeats(reduced_root)
    t_5 = time.time()
    # print(t_1 - t_0, t_2 - t_1, t_3 - t_2, t_4 - t_3, t_5 - t_4)
    return list(traversal(reduced_root))


def _compare(work_data):
    i, spdi_variant, reference, model = work_data
    output = f"{i} | {spdi_variant} | "
    # print(f"- {i}: {spdi_variant}")

    observed_2 = mutate_raw({"reference": reference}, [model])
    mutalyzer_output = spdi_converter(spdi_variant)

    descriptions = extract(reference, observed_2)
    if len(descriptions) > 1:
        output += "Error | multiple descriptions. | " + str(descriptions) + " | "
        output += mutalyzer_output["normalized_description"].split(".")[2]
        # print("Error: multiple descriptions.")
        # print(descriptions)
        # extract_main(reference, observed_2)
    elif (
        to_hgvs(descriptions[0], reference)
        != mutalyzer_output["normalized_description"].split(".")[2]
    ):
        output += "Error | different description when compared to Mutalyzer | "
        output += to_hgvs(descriptions[0], reference) + " | "
        output += mutalyzer_output["normalized_description"].split(".")[2]

        # print("Error: different description when compared to Mutalyzer")
        # print(to_hgvs(descriptions[0], reference))
        # print(mutalyzer_output["normalized_description"].split(".")[2])
    else:
        output += "OK"
    print(output)


def _spdi_model(description):
    ref_id, position, deleted, inserted = description.split(":")

    start = end = int(position)

    if deleted:
        end += len(deleted)

    delins_model = {
        "type": "deletion_insertion",
        "location": {
            "start": {"position": start, "type": "point"},
            "end": {"position": end, "type": "point"},
            "type": "range",
        },
        "deleted": [],
        "inserted": [],
    }
    if inserted:
        delins_model["deleted"] = []
        delins_model["inserted"] = [{"sequence": inserted, "source": "description"}]
    return ref_id, delins_model


def _normalize_spdi(description):
    print(description)
    ref_id, model = _spdi_model(description)
    reference = get_reference_model(ref_id)["sequence"]["seq"]
    observed = mutate_raw({"reference": reference}, [model])
    extract_main(reference, observed)


def compare_spdi(file_path, file_results=None):
    def get_results():
        if file_results:
            with open(file_results) as f:
                for line in f:
                    values = line.split("|")
                    if len(values) < 3:
                        continue
                    i, d, status = line.split("|")[:3]
                    i = i.strip()
                    d = d.strip()
                    status = d.strip()
                    results[d] = status

    def get_models():
        with open(file_path) as f:
            for line in f:
                variants.append(line.strip())

        for spdi_variant in variants:
            _convert(spdi_variant)

    def _convert(v):
        ref_id, delins_model = _spdi_model(v)
        variant_models[v] = delins_model
        reference_ids.add(ref_id)

    results = {}
    get_results()
    # print(len(results))

    variant_models = {}
    reference_ids = set()

    variants = []

    get_models()
    reference = get_reference_model(list(reference_ids)[0])["sequence"]["seq"]

    work = [
        (i, v, reference, variant_models[v])
        for i, v in enumerate(variants)
        if v not in results
    ]

    # print(len(work))
    for w in work:
        _compare(w)
    # p = Pool(1)
    # p.map(_compare, work)


def compare_spdi_summary(file_results):
    reference = get_reference_model("NG_016465.4")["sequence"]["seq"]

    ok = 0
    mutalyzer = 0
    multiple = 0
    descriptions_multiple = []
    timeout = 0
    descriptions_timeout = []
    other = 0
    total = 0

    duplications = 0
    duplications_as_insertions = 0
    duplications_same = 0
    repeats = 0
    mut_other = 0
    descriptions_mut_other = []

    with open(file_results) as f:
        for line in f:
            total += 1
            parts = line.split(" | ")
            if "OK" in line:
                ok += 1
            elif "Mutalyzer" in line:
                mutalyzer += 1
                i, d, status, message, new, old = parts
                old = old.strip()
                if "dup" in old:
                    duplications += 1
                    if "ins" in new:
                        duplications_as_insertions += 1
                        inserted = new.split("ins")[1]
                        start = int(new.split("_")[0])
                        if len(inserted) == 1 and reference[start - 1] == inserted:
                            if f"{start}dup" == old:
                                duplications_same += 1
                        elif reference[start - len(inserted) : start] == inserted:
                            if f"{start-len(inserted)+1}_{start}dup" == old:
                                duplications_same += 1
                elif "[" in old and ";" not in old:
                    repeats += 1
                else:
                    # print(line)
                    mut_other += 1
                    descriptions_mut_other.append((parts[1], new, old))

            elif "multiple" in line:
                multiple += 1
                descriptions_multiple.append((parts[1], parts[5].strip()))
            elif "timeout" in line:
                descriptions_timeout.append(parts[1])
                timeout += 1
            else:
                other += 1

    print("\nOverview")
    print("----------")
    l_s = 30
    m_s = 10
    print(f"{'Same as Mutalyzer':{l_s}} {ok:{m_s}} {ok/total*100:5.3}%")
    print(
        f"{'Different than Mutalyzer':{l_s}} {mutalyzer:>{m_s}} {mutalyzer/total*100:5.3}%"
    )
    print(f"{'Multiple paths':{l_s}} {multiple:>{m_s}}")
    print(f"{'Long processing times':{l_s}} {timeout:>{m_s}}")
    print(f"{'Other':{l_s}} {other:>{m_s}} ")
    print(f"{'Total':{l_s}} {total:>{m_s}}")

    print("\nDifferent than Mutalyzer")
    print("----------")

    l_s = 20
    m_s = 10
    print(
        f"{'duplications':{l_s}} {duplications:{m_s}} {duplications/mutalyzer*100:5.3}%"
    )
    print(f"{' duplications same':{l_s}} {duplications_same:{m_s}}")
    print(f"{'repeats':{l_s}} {repeats:{m_s}} {repeats/mutalyzer*100:5.3}%")
    print(f"{'mut_other':{l_s}} {mut_other:{m_s}} {mut_other/mutalyzer*100:5.3}%")

    print("\nDifferent than Mutalyzer - other")
    print("----------")
    for description in descriptions_mut_other:
        print(description)

    print("\nMultiple paths")
    print("----------")
    for description in descriptions_multiple:
        print(description)

    print("\nLong processing times")
    print("----------")
    for description in descriptions_timeout:
        print(description)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extractor sandbox")
    commands = parser.add_subparsers(dest="command", required=True, help="Commands")

    extract_parser = commands.add_parser("main")
    extract_parser.add_argument("reference", help="reference sequence")
    extract_parser.add_argument("observed", help="observed sequence / [variants]")

    extract_parser = commands.add_parser("extract")
    extract_parser.add_argument("reference", help="reference sequence")
    extract_parser.add_argument("observed", help="observed sequence / [variants]")

    check_parser = commands.add_parser("check")
    check_parser.add_argument("max_length", help="max_length", type=int)
    check_parser.add_argument("min_length", help="min_length", type=int)
    check_parser.add_argument("p", type=float)
    check_parser.add_argument("mu_deletion", help="mu_deletion", type=float)
    check_parser.add_argument("mu_insertion", help="mu_insertion", type=float)

    equals_parser = commands.add_parser("equals")

    check_set_parser = commands.add_parser("check_set")

    compare_parser = commands.add_parser("compare")

    compare_spdi_parser = commands.add_parser("compare_spdi")
    compare_spdi_parser.add_argument("--file", help="file path")
    compare_spdi_parser.add_argument("--results", help="already results file path")

    compare_spdi_summary_parser = commands.add_parser("compare_spdi_summary")
    compare_spdi_summary_parser.add_argument("results", help="results file path")

    spdi_parser = commands.add_parser("spdi")
    spdi_parser.add_argument("description", help="SPDI description")

    args = parser.parse_args()

    if args.command == "main":
        extract_main(args.reference, args.observed)
    elif args.command == "extract":
        print(extract(args.reference, args.observed))
    elif args.command == "check":
        check(
            args.max_length,
            args.min_length,
            args.p,
            args.mu_deletion,
            args.mu_insertion,
        )
    elif args.command == "equals":
        get_consecutive_equals()
    elif args.command == "check_set":
        check_set()
    elif args.command == "compare":
        compare()
    elif args.command == "compare_spdi":
        compare_spdi(args.file, args.results)
    elif args.command == "compare_spdi_summary":
        compare_spdi_summary(args.results)
    elif args.command == "spdi":
        _normalize_spdi(args.description)

from algebra.algebra import edit, build


def influence_interval(ref, obs):

    if ref == obs:
        raise Exception("Equal to the reference")

    d, g = edit(ref, obs)
    ops = build(g, ref, obs)
    min_pos = min(ops, key=lambda t: t[0])[0]
    max_pos = max(ops, key=lambda t: t[0])[0]

    return min_pos, max_pos

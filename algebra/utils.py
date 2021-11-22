from normalizer.description import Description
from normalizer.reference import get_reference_model
from algebra.algebra import heur_func
from algebra.simple import edit


def ref2seq(ref):
    r = get_reference_model(ref)
    return r["sequence"]["seq"]


def hgvs2ref(hgvs):
    d = Description(hgvs, stop_on_error=True)
    d.normalize()
    return d.references[d.input_model["reference"]["id"]]["sequence"]["seq"]


def hgvs2obs(hgvs):
    d = Description(hgvs, stop_on_error=True)
    d.normalize()
    return d.references["observed"]["sequence"]["seq"]


def overlay_heuristic(reference, observed, matrix):
    for i in range(len(reference) + 1):
        for j in range(len(observed) + 1):
            matrix[i][j] += heur_func(reference, observed, i, j)

    return matrix


def print_matrix_tex(matrix, reference, observed, max_dist):

    print(f'\\begin{{array}}{{rc|{"c" * (len(observed) + 1)}}}')

    print('\\renewcommand*{\\arraystretch}{1.2}')

    tmp = ' & '.join(f'\\phantom{{\\circled{{{str(i)}}}}}' for i in range(len(observed) + 1))
    print(f' & & {tmp} \\\\')

    tmp = ' & '.join(f'{str(i)}' for i in range(len(observed) + 1))
    print(f' & & {tmp} \\\\')

    tmp = ' & '.join(f'\\texttt{{{ch}}}' for ch in observed)
    print(f' & & \\texttt{{.}} & {tmp} \\\\')

    print(f'\\hline')

    for row in range(len(reference) + 1):
        if row == 0:
            print(f'{str(row)} & \\texttt{{.}} & ', end='')
        else:
            print(f'{str(row)} & \\texttt{{{reference[row - 1]}}} & ', end='')

        tmp = []
        for col in range(len(observed) + 1):
            if matrix[row][col] <= max_dist or max_dist == 0:
                if row > 0 and col > 0 and reference[row - 1] == observed[col - 1]:  # Match
                    tmp.append(f'\\circled{{{matrix[row][col]}}}')
                else:
                    tmp.append(f'{matrix[row][col]}')
            else:
                tmp.append('')
        print(' & '.join(tmp), end=' \\\\\n')
    print(f'\\end{{array}}')


def print_tex(reference, observation, heur=False, max_dist=0):
    matrix = edit(reference, observation)

    if not heur:
        print_matrix_tex(matrix, reference, observation, max_dist)
    else:
        heur_matrix = overlay_heuristic(reference, observation, matrix)
        print_matrix_tex(heur_matrix, reference, observation, max_dist)

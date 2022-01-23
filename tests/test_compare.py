from algebra import algebra


def pytest_generate_tests(metafunc):
    values = []
    ids = []
    with open('tests/tests.txt') as file:
        for line in file.readlines():
            args = line.strip().split(', ')
            values.append(args)
            ids.append('-'.join([str(arg) for arg in args]))

    names = ['reference', 'lhs', 'rhs', 'expected']
    metafunc.parametrize(argnames=names, argvalues=values, ids=ids)


def test_compare(reference, lhs, rhs, expected):
    relation, _, _ = algebra.compare(reference, lhs, rhs)
    assert relation == expected

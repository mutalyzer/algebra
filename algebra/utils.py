import random


def random_sequence(max_length, min_length=0, alphabet="ACGT", weights=None):
    return "".join(random.choices(alphabet, weights=weights, k=random.randint(min_length, max_length)))

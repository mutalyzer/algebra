import random


def random_sequence(max_length, min_length=0, alphabet="ACGT"):
    return "".join(random.choice(alphabet) for _ in range(random.randint(min_length, max_length)))

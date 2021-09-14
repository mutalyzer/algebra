from normalizer.description import Description
from normalizer.reference import get_reference_model


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

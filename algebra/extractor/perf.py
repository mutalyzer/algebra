from sys import stderr
from time import perf_counter_ns


class Perf:
    _tics = []

    def __new__(self):
        return getattr(self, "_instance", super().__new__(self))

    def __enter__(self):
        self.tic()

    def __exit__(self, *_):
        self.toc()

    def tic(self):
        self._tics.append(perf_counter_ns())

    def toc(self):
        if self._tics:
            ns = perf_counter_ns() - self._tics.pop()
            print(f"Elapsed time {ns / 1_000_000_000:.3f} seconds.", file=stderr)
            return ns


tic = Perf().tic
toc = Perf().toc

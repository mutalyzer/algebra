from .variant import Variant


class Parser:
    def __init__(self, expression):
        if not isinstance(expression, str):
            raise TypeError("expression must be a string")

        self.pos = 0
        self.expression = expression

    def _match(self, word):
        if not self.expression[self.pos:self.pos + len(word)] == word:
            return False

        self.pos += len(word)
        return True

    def _match_digit(self):
        if self.pos >= len(self.expression):
            raise ValueError("unexpected end of expression")
        if not self.expression[self.pos].isdigit():
            raise ValueError(f"expected digit at {self.pos}")

        self.pos += 1
        return int(self.expression[self.pos - 1])

    def _match_nucleotide(self):
        if self.pos >= len(self.expression):
            raise ValueError("unexpected end of expression")
        if not self.expression[self.pos] in "ACGT":
            raise ValueError(f"expected nucleotide at {self.pos}")

        self.pos += 1
        return self.expression[self.pos - 1]

    def _match_number(self):
        number = self._match_digit()
        while True:
            try:
                number = number * 10 + self._match_digit()
            except ValueError:
                return number

    def _match_sequence(self):
        start = self.pos
        self._match_nucleotide()
        while True:
            try:
                self._match_nucleotide()
            except ValueError:
                return self.expression[start:self.pos]

    def _match_variant(self):
        start = self._match_number()
        if start == 0:
            raise ValueError(f"invalid position at {self.pos}")
        start -= 1
        end = None
        if self._match("_"):
            end = self._match_number()
            if start >= end - 1:
                raise ValueError(f"invalid range at {self.pos}")

        if self._match("del"):
            if end is None:
                end = start + 1
            deleted = None
            try:
                deleted = self._match_sequence()
            except ValueError:
                pass

            if deleted is not None and len(deleted) != end - start:
                raise ValueError(f"inconsistent deleted length at {self.pos}")
            if self._match("ins"):
                return Variant(start, end, self._match_sequence())
            return Variant(start, end)

        if self._match("ins"):
            if end is None or end - start != 2:
                raise ValueError(f"invalid inserted range at {self.pos}")
            return Variant(start + 1, start + 1, self._match_sequence())

        if end is not None:
            raise ValueError(f"invalid variant at {self.pos}")

        try:
            self._match_nucleotide()
        except ValueError:
            pass
        if not self._match(">"):
            raise ValueError(f"expected '>' at {self.pos}")
        return Variant(start, start + 1, self._match_nucleotide())

    def hgvs(self):
        if self._match("="):
            if not self.pos == len(self.expression):
                raise ValueError(f"expected end of expression at {self.pos}")
            return []

        if self._match("["):
            variants = [self._match_variant()]
            while self._match(";"):
                variants.append(self._match_variant())
            if not self._match("]"):
                raise ValueError(f"expected ']' at {self.pos}")
        else:
            variants = [self._match_variant()]

        if not self.pos == len(self.expression):
            raise ValueError(f"expected end of expression at {self.pos}")
        return variants

    def spdi(self):
        self._match_sequence()
        if not self._match(":"):
            raise ValueError(f"expected ':' at {self.pos}")

        position = self._match_number()
        if not self._match(":"):
            raise ValueError(f"expected ':' at {self.pos}")

        length = 0
        try:
            length = self._match_number()
        except ValueError:
            try:
                length = len(self._match_sequence())
            except ValueError:
                pass
        if not self._match(":"):
            raise ValueError(f"expected ':' at {self.pos}")

        inserted = ""
        try:
            inserted = self._match_sequence()
        except ValueError:
            pass
        if not self.pos == len(self.expression):
            raise ValueError(f"expected end of expression at {self.pos}")
        return [Variant(position, position + length, inserted)]

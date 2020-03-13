from io import StringIO
import sys
from typing import Any, Callable, Dict, List, Tuple

from .sparse_table import SparseTable, Pos

SupportTable = SparseTable[float]
"""
Typedef to represent a sparse support table implemented as a
dictionary (for now).
"""


def load_support_table(input_file: StringIO) -> SupportTable:
    """
    Deserialize a support table as serialized by the Perl prototype
    implementation. The format of the first column ("key") is "<row>.<column>".
    The second column ("support_value") is simply a floating point value.

        key	support_value
        0.0	0.0104846536699391
        0.1	0.0104846536699391
        0.2	0.0104846536699391
        ...

    This function assumes that the table was serialized with column headers.

    >>> from io import StringIO
    >>> s = load_support_table(StringIO("key\tsupport_value\\n1.2\t0.1"))
    >>> s[Pos(1, 2)]
    0.1
    >>> s = load_support_table(StringIO("1.2\t0.1"))
    Traceback (most recent call last):
      ...
    ValueError: invalid headers: ['1.2', '0.1']
    >>> s = load_support_table(StringIO("key support_value\\n1.2 0.1"))
    >>> s[Pos(1, 2)]
    0.1
    """
    # Verify the headers
    headers = next(input_file).split()
    if headers != ["key", "support_value"]:
        raise ValueError(f"invalid headers: {headers}")

    support_table: SupportTable = SupportTable()

    # Populate the actual table values
    for line in input_file:
        # Skip empty lines
        if len(line) == 0:
            continue

        key, value = line.split()
        row, col = key.split(".")
        value_pos = Pos(int(row), int(col))
        support_table[value_pos] = float(value)

    return support_table


def save_support_table(
    support_table: SupportTable,
    output: Callable[[str], Any] = sys.stdout.write,
) -> None:
    """
    Serialize the support table in a manner that is compatible with the
    serialization implemented in the Perl prototype.

    Use `output` to capture output for testing, otherwise the default value
    should be sufficient.

    >>> import io
    >>> s = SupportTable()
    >>> s[Pos(0, 0)] = 1.0
    >>> s[Pos(0, 1)] = 2.0
    >>> r = io.StringIO()
    >>> save_support_table(s, r.write)
    >>> print(r.getvalue())
    key support_value
    0.0 1.0
    0.1 2.0
    <BLANKLINE>
    """
    output("key support_value")
    for value_pos in support_table:
        support_value = support_table[value_pos]
        output(f"\n{value_pos[0]}.{value_pos[1]} {support_value}")
    output("\n")


def fill_support_table() -> SupportTable:
    pass

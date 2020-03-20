import math
import sys
from timeit import default_timer as timer
from typing import Any, Callable, TextIO, List, Optional, Tuple

from .constants import NAN_STRING, DEFAULT_CHANGE_PROB, DEFAULT_CHUNK_SIZE
from .sparse_table import DictSparseTable, Dim, Pos, SparseTable
from .support_table import SupportTable

ProbTable = SparseTable[float]


OriginTable = SparseTable[int]


def load_prob_origin_tables(
    input_file: TextIO,
) -> Tuple[ProbTable, OriginTable]:
    """
    Deserialize a probability table as serialized by the Perl prototype
    and Python implementations. The format of the first column ("key") is
    "<row>.<column>" (as with the support table). The second column
    ("prob_value") is simply a floating point probability value. The third
    column is the row in the previous column from which the value was
    derived, called "origin".

        key	prob_value origin
        0.0	0.0104846536699391 0
        0.1	0.0104846536699391 0
        0.2	0.0104846536699391 1
        ...

    This function assumes that the table was serialized with column headers.

    >>> from io import StringIO
    >>> p, o = load_prob_origin_tables(StringIO("key\tprob_value\torigin\\n1.2\t0.1\t0"))
    >>> Pos(1, 2) in p
    True
    >>> p[Pos(1, 2)]
    0.1
    >>> Pos(1, 2) in o
    True
    >>> o[Pos(1, 2)]
    0
    >>> load_prob_origin_tables(StringIO("1.2\t0.1\t0"))
    Traceback (most recent call last):
      ...
    ValueError: invalid headers: ['1.2', '0.1', '0']
    >>> p, o = load_prob_origin_tables(StringIO("key prob_value origin\\n1.2 0.1 0"))
    >>> Pos(1, 2) in p
    True
    >>> p[Pos(1, 2)]
    0.1
    >>> Pos(1, 2) in o
    True
    >>> o[Pos(1, 2)]
    0
    """
    headers = next(input_file).split()
    if headers != ["key", "prob_value", "origin"]:
        raise ValueError(f"invalid headers: {headers}")

    prob_table: ProbTable = DictSparseTable()
    origin_table: OriginTable = DictSparseTable()

    for line in input_file:
        # Ignore blank lines
        if len(line) == 0:
            continue

        key, prob, origin = line.split()
        row, col = key.split(".")
        value_pos = Pos(int(row), int(col))
        prob_table[value_pos] = float(prob)

        if origin == NAN_STRING:
            parsed_origin = -1
        else:
            parsed_origin = int(origin)
        origin_table[value_pos] = parsed_origin

    return prob_table, origin_table


def save_prob_origin_tables(
    prob_table: ProbTable,
    origin_table: OriginTable,
    output: Callable[[str], Any] = sys.stdout.write,
) -> None:
    """
    Serialize the probability and origin tables, as returned by
    `fill_prob_table`, in a manner that is compatible with the serialization
    implemented in the Perl prototype.

    Use `output` to capture output for testing, otherwise the default value
    should be sufficient.

    >>> import io
    >>> p = DictSparseTable()
    >>> p[Pos(0, 0)] = 1.0
    >>> p[Pos(0, 1)] = 2.0
    >>> o = DictSparseTable()
    >>> o[Pos(0, 0)] = 1
    >>> o[Pos(0, 1)] = 2
    >>> r = io.StringIO()
    >>> save_prob_origin_tables(p, o, r.write)
    >>> print(r.getvalue())
    key prob_value origin
    0.0 1.0 1
    0.1 2.0 2
    <BLANKLINE>
    """
    output("key prob_value origin")
    for value_pos in prob_table.keys():
        prob_value = prob_table[value_pos]
        origin_value = origin_table[value_pos]

        if origin_value == -1:
            origin_string = NAN_STRING
        else:
            origin_string = str(origin_value)

        output(f"\n{value_pos[0]}.{value_pos[1]} {prob_value} {origin_string}")
    output("\n")


def fill_prob_origin_tables(
    support_matrix: SupportTable,
    columns: Optional[List[int]] = None,
    row_count: Optional[int] = None,
    benchmark: bool = False,
    change_prob: float = DEFAULT_CHANGE_PROB,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> Tuple[ProbTable, OriginTable]:
    """
    Fills in the probability score matrix from the support matrix. Also fills
    the origin matrix for convenience. Returns a tuple that contains the
    probability matrix and the origin matrix, in that order..

    support_matrix -- a column-major nested list (matrix) of support scores
    columns        -- the list of sequence column indices, treated as
                      left-most edges of the chunk window, that should be
                      included
    row_count      -- the total number of rows that will be processed, will
                      be computed from the other arguments if omitted
    benchmark      -- whether or not to collect and output benchmarking data
                      to the logger
    change_prob    -- probability of changing sequence
    chunk_size     -- the width of the "window" analyzed at one time

    We omit filling the first column because we want its probabilities to be
    zero anyway.

    The basic algorithm is described below. All calculations happen in log
    space.

    look at all i's in j-1
        mult by confidence in current cell
        if comes from same i, mult by higher prob
        else - mult by lower prob /(numseqs-1) -> so sum of all probs == 1
     return max
    """
    same_prob = 1.0 - change_prob

    logged_change_prob = math.log(change_prob)
    logged_same_prob = math.log(same_prob)

    if benchmark:
        start = timer()

    prob_table: ProbTable = DictSparseTable()
    origin_table: OriginTable = OriginTable

    # TODO: We may be able to get row_count directly from support_table now?
    if columns is None or row_count is None:
        max_col = 0
        max_row = 0
        for (row, col) in support_matrix:
            if col > max_col:
                max_col = col
            if row > max_row:
                max_row = row

        if columns is None:
            columns = list(range(max_col + 1))

        if row_count is None:
            row_count = max_row + 1

    for row_index in range(row_count):
        origin_table[Pos(row_index, 0)] = -1
        prob_table[Pos(row_index, 0)] = 0.0

    # colIndex is j in the Perl version
    col_index = 1

    # Skip the first column because we already filled it above.
    # We need to run off the end of the last column by `chunk_size` since we
    # slide the window along one column at a time. So we augment the selected
    # columns with the necessary additional columns.
    # colIndexIndex is col in the Perl version
    for col_index_index in range(1, len(columns) + chunk_size - 1):
        in_provided_columns = col_index_index < len(columns)

        # Matches behavior on lines 819-823
        if in_provided_columns:
            col_index = columns[col_index_index]
        else:
            col_index += 1

        for row_index in range(row_count):
            max_score = -math.inf
            max_score_index = 0

            # Mapped from line 831
            support_value = support_matrix.get(Pos(row_index, col_index), None)
            if support_value is not None and support_value != 0:
                supportLog = math.log(support_value)

                for innerRowIndex in range(row_count):
                    # Omitting check on line 841 because we expanded columns
                    if in_provided_columns:
                        probValue = prob_table.get(
                            (innerRowIndex, columns[col_index_index - 1]), None
                        )
                    else:
                        probValue = prob_table.get(
                            (innerRowIndex, col_index - 1), None
                        )

                    if probValue is not None:
                        score = supportLog + probValue

                        if row_index == innerRowIndex:
                            # Same sequence (higher probability)
                            score += logged_same_prob
                        else:
                            # Different sequence (lower probability)
                            score += logged_change_prob

                        # Update our current max score
                        if score > max_score:
                            max_score = score
                            max_score_index = innerRowIndex

                prob_table[(row_index, col_index)] = max_score
                origin_table[(row_index, col_index)] = max_score_index

    if benchmark:
        end = timer()
        from logging import getLogger
        logger = getLogger(__name__)
        logger.info(f"benchmark: {end - start}s")

    return (prob_table, origin_table)

from typing import Dict, List, NamedTuple, Tuple

RowColumn = Tuple[int, int]

ConfidenceMatrix = Dict[RowColumn, float]

ConsensusMatrix = Dict[RowColumn, int]

StrandMatrix = Dict[RowColumn, str]

SupportMatrix = Dict[RowColumn, float]


class ConsensusMatrixContainer(NamedTuple):
    """
    A container for a consensus matrix along with its context.
    """

    active_columns: List[int]
    """
    A list of column indices in the matrix that are not empty.
    This lets us avoid looping through unnecessary columns.
    """

    active_rows: Dict[int, List[int]]
    """
    Mapping from column index (starting at zero) to a list of row
    indices in that column that are "active" (have data in them).
    This lets us avoid looping through unnecessary rows.
    """

    matrix: ConsensusMatrix
    """
    TODO: Clarify what the values are, perhaps?
    
    Hash implementation of a sparse 2D matrix used along with DP
    matrices. The key is a 2-tuple of row and column indices
    (respectively) and the values are the alignment positions in the
    consensus subfamily sequence.
    """


class CollapsedMatrices(NamedTuple):
    row_num_update: int
    """
    TODO: What, exactly, is this?
    updates number of rows in matrices
    """

    consensus_matrix: ConsensusMatrix
    # consensus_matrix_collapse: collapsed version

    strand_matrix: StrandMatrix
    # strand_matrix_collapse: Hash implementation of sparse 2D matrix used along with DP matrices.
    # Tuple[str, int] as key. String is the subfamily name of the row, int is the column in matrix.
    # This is a collapsed matrix with no redundant subfamilies as rows. Each cell in matrix is
    # the strand of the consensus sequence that aligned at the corresponding column position of
    # the sequence.

    support_matrix: SupportMatrix
    # support_matrix_collapse: Collapsed version of Support matrix. There may be duplicate rows of the
    # same subfamily, collapsing the matrices puts the duplicates into the same row. Hash
    # implementation with tuple[str, int] as key. String is the subfamily name of the row, int
    # is the column.

    subfamilies: List[str]
    # subfams_collapse: Collapsed version of the array Subfams, any duplicate subfams are consolidated
    # into one.

    active_rows: Dict[int, List[int]]
    # active_cells_collapse: Not all rows in each column hold values. Dictionary that holds column
    # number as the key, and an array of which rows hold values for that column.

    subfamily_indices: Dict[str, int]
    #     subfams_collapse_temp: maps subfam names to it's new row number

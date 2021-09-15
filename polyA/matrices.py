from typing import Dict, List, NamedTuple, Tuple

RowColumn = Tuple[int, int]

AlignMatrix = Dict[RowColumn, float]

ConfidenceMatrix = Dict[RowColumn, float]

ConsensusMatrix = Dict[RowColumn, int]

NodeConfidenceMatrix = Dict[Tuple[str, int], float]

StrandMatrix = Dict[RowColumn, str]

SupportMatrix = Dict[RowColumn, float]

SubfamCol = Tuple[str, int]

RowConsensus = Tuple[int, int]

SubfamAlignmentsMatrix = Dict[SubfamCol, RowConsensus]


class ConsensusMatrixContainer(NamedTuple):
    """
    A container for a consensus matrix along with its context.
    """

    active_rows: Dict[int, List[int]]
    """
    Mapping from column index to a list of row
    indices in that column that are "active" (have data in them).
    This lets us avoid looping through unnecessary rows.
    """

    matrix: ConsensusMatrix
    """    
    Hash implementation of a sparse 2D matrix used along with DP
    matrices. The key is a 2-tuple of row and column indices
    (respectively) and the values are the alignment positions in the
    consensus subfamily sequence.
    """


class CollapsedMatrices(NamedTuple):
    row_num_update: int
    """
    all matrices used in DP are "collapsed" meaning if there are duplicate rows with the same
    subfamily they are collapsed into the same row. We choose which original row to use in the
    collapsed version by doing a mini DP to find the most probable.

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

    subfam_alignments_matrix: SubfamAlignmentsMatrix
    # subfam_alignment_matrix: Hash implementation of sparse 2D matrix.
    # Tuple[str, int] as key. String is the subfamily name of the row, int is the sequence position.
    # This is a collapsed matrix with no redundant subfamilies as rows. Each cell in the matrix is
    # the original subfam/chrom alignment row and the consensus sequence pos that aligned at the
    # corresponding column position of the sequence.

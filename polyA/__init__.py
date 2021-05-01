"""
PolyA is a sequence annotation adjudicator.
"""

__version__ = "0.1.0"

from .alignment import Alignment, get_skip_state
from .calc_repeat_scores import calculate_repeat_scores
from .calculate_score import (
    calculate_score,
    calculate_complexity_adjusted_score,
)
from .collapse_matrices import collapse_matrices
from .confidence_cm import confidence_cm, confidence_only
from .constants import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_SHARD_GAP,
    PROB_SKIP_TR,
    PROB_SKIP,
    CHANGE_PROB,
    NAN_STRING,
    SAME_PROB_LOG,
    SKIP_ALIGN_SCORE,
    START_ID,
)
from .edges import edges
from .extract_nodes import extract_nodes
from .fill_align_matrix import fill_align_matrix
from .fill_confidence_matrix import (
    fill_confidence_matrix,
)  # , trailing_edges_info
from .fill_consensus_position_matrix import fill_consensus_position_matrix
from .fill_node_confidence import fill_node_confidence
from .fill_path_graph import fill_path_graph
from .fill_probability_matrix import fill_probability_matrix
from .fill_support_matrix import fill_support_matrix
from .get_path import get_path
from .lambda_provider import (
    ConstantLambdaProvider,
    EaselLambdaProvider,
    LambdaProvider,
)
from .load_alignments import load_alignments, shard_overlapping_alignments
from .matrices import (
    ConsensusMatrix,
    CollapsedMatrices,
    ConsensusMatrixContainer,
    ConfidenceMatrix,
    SupportMatrix,
    RowColumn,
    StrandMatrix,
)
from .output import Output
from .pad_sequences import pad_sequences
from .printers import (
    print_matrix_hash,
    print_matrix_support,
    print_results,
    print_results_chrom,
    print_results_sequence,
    print_results_soda,
)
from .prior_counts import read_prior_counts
from .substitution_matrix import (
    load_substitution_matrices,
    SubMatrix,
    SubMatrixCollection,
)
from .ultra_provider import (
    UltraProvider,
    UltraOutput,
    ApplicationUltraProvider,
    TandemRepeat,
)
from ._runners import run_full, run_confidence

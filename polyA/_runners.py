from logging import Logger
import math
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple

from .alignment import Alignment
from .calc_repeat_scores import calculate_repeat_scores
from .collapse_matrices import collapse_matrices
from .confidence_cm import confidence_only
from .constants import SAME_PROB_LOG, SKIP_ALIGN_SCORE
from .edges import edges
from .extract_nodes import extract_nodes
from .fill_align_matrix import fill_align_matrix
from .fill_confidence_matrix import fill_confidence_matrix
from .fill_consensus_position_matrix import fill_consensus_position_matrix
from .fill_node_confidence import fill_node_confidence
from .fill_path_graph import fill_path_graph
from .fill_probability_matrix import fill_probability_matrix
from .fill_support_matrix import fill_support_matrix
from .get_path import get_path
from .pad_sequences import pad_sequences
from .print_helpers import find_consensus_lengths
from .printers import Printer
from .substitution_matrix import SubMatrix, SubMatrixCollection
from .ultra_provider import TandemRepeat


def run_confidence(
    alignments: List[Alignment],
    lambs: List[float],
) -> None:
    # command line option to just output confidence values for
    # single annotation instead of do whole algorithm
    """
    Calculate confidence values for single annotations instead of
    running the entire algorithm.

    :param alignments: list of alignments to run on
    :param lambs: the values of lambda to use for each alignment (from Easel)
    """
    if len(alignments) != len(lambs):
        raise RuntimeError(
            "number of alignments must match number of lambda values"
        )

    subfams = []
    subfam_rows = []
    scores = []

    for i, a in enumerate(alignments):
        subfams.append(a.subfamily)
        subfam_rows.append(i)
        scores.append(a.score)

    confidence_list = confidence_only(scores, lambs)

    # ignore this because the types work, but mypy doesn't know that
    confidence_list, subfams_copy = zip(*sorted(zip(confidence_list, subfams)))  # type: ignore

    stdout.write(f"query_label\tconfidence\n")
    for i in range(len(subfams_copy) - 1, 0, -1):
        if confidence_list[i]:
            stdout.write(f"{subfams_copy[i]}\t{confidence_list[i]}\n")


def _validate_target(target: Alignment) -> None:
    if target.chrom_length == 1:
        # artifical seq chr0:0000-0000
        Logger(__name__).warning(
            """No chromosome position information found
            (this is OK for artificial sequences but --viz and --heatmap will fail)"""
        )
    elif target.chrom_length <= 25:
        raise RuntimeError("Target sequence length must be >25")
    elif target.chrom_length < 1000:
        Logger(__name__).warning(
            """Target region is <1000 in length,
            did you mean to use such a short region?"""
        )


def _handle_single_alignment(
    target: Alignment,
    printer: Printer,
) -> Tuple[int, int]:
    """
    If there is only one subfam in the alignment file, no need
    to run anything because we know that subfam is the annotation.
    """
    from uuid import uuid4

    node_id = uuid4().hex
    if printer.use_matrix_position:
        printer.print_results_simple(
            1,
            target.stop - target.start + 1,
            node_id,
            target.subfamily,
        )
    elif printer.use_sequence_position:
        printer.print_results_simple(
            target.start,
            target.stop,
            node_id,
            target.subfamily,
        )
    else:
        printer.print_results_simple(
            target.start + target.chrom_start,
            target.stop + target.chrom_start,
            node_id,
            target.subfamily,
        )
    return target.start, target.stop


def _change_probs(
    seq_count: int, trans_penalty: int
) -> Tuple[float, float, float]:
    change_prob = 10 ** (-1 * trans_penalty)
    change_prob_log = math.log(change_prob / (seq_count - 1))

    # jumping in and then out of the skip state counts as 1 jump
    change_prob_skip = change_prob_log / 2

    # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump
    same_prob_skip = change_prob_log / 30

    return change_prob_log, change_prob_skip, same_prob_skip


def run_full(
    alignments: List[Alignment],
    tandem_repeats: List[TandemRepeat],
    chunk_size: int,
    trans_penalty: int,
    sub_matrix_scores: SubMatrixCollection,
    subfam_counts: Dict[str, float],
    shard_start: int,
    shard_stop: int,
    prev_start: int,
    prev_stop: int,
    printer: Printer,
) -> Tuple[int, int]:
    # chunk start and stop are positions in seq
    seq_count = len(alignments)

    target = alignments[1]
    _validate_target(target)

    if len(tandem_repeats) == 0 and seq_count == 2:
        # only one alignment other than the skip state
        return _handle_single_alignment(
            target,
            printer,
        )

    change_prob_log, change_prob_skip, same_prob_skip = _change_probs(
        seq_count, trans_penalty
    )

    # We use a whole bunch of parallel lists to store alignment data. This may
    # change in the future, but for now it is an artifact of an earlier version
    # of the code and it has some benefits in that we can more easily iterate
    # over a particular field. Each of these lists holds data for exactly one
    # alignment (one row in the matrices we build later).
    alignment_subfamilies: List[str] = []
    alignment_chromosomes: List[str] = []
    alignment_scores: List[int] = []
    alignment_strands: List[str] = []
    alignment_start_positions: List[int] = []
    alignment_stop_positions: List[int] = []
    alignment_consensus_start_positions: List[int] = []
    alignment_consensus_stop_positions: List[int] = []
    alignment_subfamily_sequences: List[str] = []
    alignment_chromosome_sequences: List[str] = []
    alignment_flanks: List[int] = []
    alignment_substitution_matrices: List[SubMatrix] = []
    alignment_gap_inits: List[float] = []
    alignment_gap_exts: List[float] = []
    alignment_lambdas: List[float] = []

    for a in alignments:
        alignment_subfamilies.append(a.subfamily)
        alignment_chromosomes.append(a.chrom_name)
        alignment_scores.append(a.score)
        alignment_strands.append(a.strand)
        alignment_start_positions.append(a.start)
        alignment_stop_positions.append(a.stop)
        alignment_consensus_start_positions.append(a.consensus_start)
        alignment_consensus_stop_positions.append(a.consensus_stop)
        alignment_subfamily_sequences.append(a.subfamily_sequence)
        alignment_chromosome_sequences.append(a.sequence)
        alignment_flanks.append(a.flank)
        alignment_gap_inits.append(a.gap_init)
        alignment_gap_exts.append(a.gap_ext)

        # fixme (george): is this a list of dicts?
        # fixme (george): change it to a dict of dicts?
        if a.sub_matrix_name in sub_matrix_scores:
            alignment_substitution_matrices.append(
                sub_matrix_scores[a.sub_matrix_name]
            )
            alignment_lambdas.append(sub_matrix_scores[a.sub_matrix_name].lamb)
        else:
            alignment_substitution_matrices.append(SubMatrix("skip", 0.0))
            alignment_lambdas.append(0.0)

    tr_consensus_matrix: List[str] = []
    if tandem_repeats:
        for tr in tandem_repeats:
            tr_consensus_matrix.append(tr.consensus)
            # add TR starts and stops
            # check they do not exceed chunk boundary
            if tr.start < shard_start:
                alignment_start_positions.append(shard_start)
            else:
                alignment_start_positions.append(tr.start)
            if tr.stop > shard_stop:
                alignment_stop_positions.append(shard_stop)
            else:
                alignment_stop_positions.append(tr.stop)

    # The implicit start and stop positions for the skip state are the index
    # prior to the minimum start and the index after the maximum stop. In other
    # words, the skip state has N+2 columns, where N is the number of columns
    # spanned by the alignments.
    alignment_start_positions[0] = min(alignment_start_positions[1:]) - 1
    alignment_stop_positions[0] = max(alignment_stop_positions[1:]) + 1

    # save original alignments before padding for SODA viz
    original_subfamily_sequences = list(alignment_subfamily_sequences)
    original_chromosome_sequences = list(alignment_chromosome_sequences)

    pad_sequences(
        chunk_size,
        alignment_subfamily_sequences,
        alignment_chromosome_sequences,
    )

    start_all, stop_all = edges(
        alignment_start_positions, alignment_stop_positions
    )

    # Plus one because the range is inclusive, plus two because
    # we have the implicit first and last positions that only
    # have data for the skip state.
    column_count = stop_all - start_all + 1 + 2

    # number of rows in matrices
    row_count: int = seq_count

    align_matrix = fill_align_matrix(
        alignment_lambdas,
        column_count,
        start_all,
        chunk_size,
        alignment_gap_inits,
        alignment_gap_exts,
        SKIP_ALIGN_SCORE,
        alignment_subfamily_sequences,
        alignment_chromosome_sequences,
        alignment_start_positions,
        alignment_stop_positions,
        [sm.scores for sm in alignment_substitution_matrices],
        [sm.background_freqs for sm in alignment_substitution_matrices],
    )

    (active_cells, consensus_matrix) = fill_consensus_position_matrix(
        row_count,
        column_count,
        start_all,
        alignment_subfamily_sequences,
        alignment_chromosome_sequences,
        alignment_start_positions,
        alignment_stop_positions,
        alignment_consensus_start_positions,
        alignment_strands,
    )

    repeat_scores: Dict[int, float] = {}
    if len(tandem_repeats) > 0:
        repeat_scores = calculate_repeat_scores(
            tandem_repeats,
            chunk_size,
            start_all,
            row_count,
            active_cells,
            align_matrix,
            consensus_matrix,
            shard_start,
            shard_stop,
        )

        # add TRs to subfams
        for _ in tandem_repeats:
            alignment_subfamilies.append("Tandem Repeat")
            alignment_strands.append("+")
            row_count += 1

    confidence_matrix = fill_confidence_matrix(
        column_count,
        subfam_counts,
        alignment_subfamilies,
        active_cells,
        align_matrix,
    )

    support_matrix = fill_support_matrix(
        chunk_size,
        alignment_start_positions,
        alignment_stop_positions,
        confidence_matrix,
    )

    non_empty_columns = [c for c in range(column_count)]

    collapsed_matrices = collapse_matrices(
        row_count,
        start_all,
        non_empty_columns,
        alignment_subfamilies,
        alignment_strands,
        alignment_start_positions,
        alignment_stop_positions,
        active_cells,
        support_matrix,
        consensus_matrix,
    )

    support_matrix_collapse = collapsed_matrices.support_matrix
    subfams_collapse = collapsed_matrices.subfamilies
    active_cells_collapse = collapsed_matrices.active_rows
    consensus_matrix_collapse = collapsed_matrices.consensus_matrix
    strand_matrix_collapse = collapsed_matrices.strand_matrix
    subfams_collapse_index = collapsed_matrices.subfamily_indices
    row_count = collapsed_matrices.row_num_update
    subfam_alignments_collapse = collapsed_matrices.subfam_alignments_matrix

    align_matrix.clear()
    confidence_matrix.clear()
    support_matrix.clear()
    consensus_matrix.clear()

    if tandem_repeats:
        # give different TRs consensus positions that don't allow them to be stitched
        tr_consensus_pos = 1000000
        prev_tr_col = 0
        for tr_col in repeat_scores:
            if prev_tr_col != tr_col - 1:
                tr_consensus_pos -= 500
            consensus_matrix_collapse[
                row_count - 1, tr_col + 1
            ] = tr_consensus_pos
            prev_tr_col = tr_col

    (
        prob_matrix_last_column,
        origin_matrix,
        same_subfam_change_matrix,
    ) = fill_probability_matrix(
        same_prob_skip,
        SAME_PROB_LOG,
        change_prob_log,
        change_prob_skip,
        non_empty_columns,
        collapsed_matrices,
    )

    # node_ids for each nucleotide will be assigned during DP backtrace
    node_ids = [""] * column_count

    (changes_position, changes) = get_path(
        non_empty_columns,
        node_ids,
        subfams_collapse,
        prob_matrix_last_column,
        active_cells_collapse,
        origin_matrix,
        same_subfam_change_matrix,
    )

    # keep the original annotation for reporting results
    changes_orig = changes.copy()
    changes_position_orig = changes_position.copy()
    non_empty_columns_orig = non_empty_columns.copy()

    node_confidence: Dict[Tuple[str, int], float] = {}
    node_confidence_orig: Dict[Tuple[str, int], float] = {}
    path_graph: List[int] = []

    # Finding inserted elements and stitching original elements
    # Steps-
    # 1. Find confidence for nodes
    # 2. Create path graph - Find alternative paths through graph and add those edges
    # 3. Extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
    # 5. Annotate again with removed nodes
    #   ** stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left
    prev_num_nodes: int = 0
    count: int = 0
    while True:
        count += 1
        node_count = len(changes)

        node_confidence.clear()  # reuse old node_confidence matrix
        node_confidence = fill_node_confidence(
            node_count,
            start_all,
            alignment_gap_inits,
            alignment_gap_exts,
            non_empty_columns,
            alignment_start_positions,
            alignment_stop_positions,
            changes_position,
            alignment_subfamilies,
            alignment_subfamily_sequences,
            alignment_chromosome_sequences,
            subfam_counts,
            alignment_substitution_matrices,
            repeat_scores,
            len(tandem_repeats),
        )

        # store original node confidence for reporting results
        if count == 1:
            node_confidence_orig = node_confidence.copy()

        # breakout of loop if there are 2 or less nodes left
        if node_count <= 2 or node_count == prev_num_nodes:
            break

        path_graph.clear()  # reuse old path_graph

        path_graph = fill_path_graph(
            node_count,
            non_empty_columns,
            changes,
            changes_position,
            consensus_matrix_collapse,
            strand_matrix_collapse,
            node_confidence,
            subfams_collapse_index,
        )

        # test to see if there are nodes in the graph that have more than one incoming or outgoing edge,
        # if so - keep looping, if not - break out of the loop
        # if they are all 0, break out of the loop
        test: bool = False
        for j in range(node_count):
            for i in range(j - 1):
                if path_graph[i * node_count + j] == 1:
                    test = True

        if not test:
            break

        extract_nodes(
            node_count,
            non_empty_columns,
            changes_position,
            path_graph,
        )

        # run DP calculations again with nodes corresponding to inserted elements removed
        # ignores removed nodes because they are no longer in NonEmptyColumns
        (
            prob_matrix_last_column,
            origin_matrix,
            same_subfam_change_matrix,
        ) = fill_probability_matrix(
            same_prob_skip,
            SAME_PROB_LOG,
            change_prob_log,
            change_prob_skip,
            non_empty_columns,
            collapsed_matrices,
        )

        (changes_position, changes) = get_path(
            non_empty_columns,
            node_ids,
            subfams_collapse,
            prob_matrix_last_column,
            active_cells_collapse,
            origin_matrix,
            same_subfam_change_matrix,
        )
        prev_num_nodes = node_count

    # handles TRs overlapping shards
    i = 0
    # get first subfam in shard
    tr_overlap_index: int = -1
    while i < len(changes_orig) and changes_orig[i] == "skip":
        i += 1

    if i == len(changes_orig):
        # whole shard is skip state
        last_subfam_start = -1
        last_subfam_stop = -1
    else:
        first_subfam_start = (
            non_empty_columns_orig[changes_position_orig[i]] + start_all - 1
        )
        # check if first print pos is one after prev stop in seq pos
        if prev_stop + 1 == first_subfam_start:
            tr_overlap_index = i
            tr_start_seq_pos: int = prev_start - start_all + 1
            stdout.write("\033[2K\033[1G")  # remove last printed line

        # get last subfam in shard
        # will find something since whole shard is not a skip state
        i = len(changes_orig) - 1
        while changes_orig[i] == "skip":
            i -= 1
        last_subfam_start = (
            non_empty_columns_orig[changes_position_orig[i]] + start_all - 1
        )
        last_subfam_stop = (
            non_empty_columns_orig[changes_position_orig[i + 1] - 1]
            + start_all
            - 1
        )

    tr_consensus_changes: Dict[int, str] = {}
    if len(tandem_repeats) > 0:
        # label TRs with their consensus sequence
        changes_index: int = 0
        while changes_index < len(changes_orig):
            subfam: str = changes_orig[changes_index]
            if subfam == "Tandem Repeat":
                cur_changes_pos: int = changes_position_orig[changes_index]

                start_seq_pos: int = (
                    non_empty_columns_orig[cur_changes_pos] + start_all - 1
                )

                if changes_index == tr_overlap_index:
                    # change start pos of first subfam to be prev start seq pos
                    non_empty_columns_orig[cur_changes_pos] = tr_start_seq_pos
                    tr_overlap_index = -1

                stop_seq_pos: int = (
                    non_empty_columns_orig[
                        changes_position_orig[changes_index + 1] - 1
                    ]
                    + start_all
                    - 1
                )

                # track consensus label row
                prev_subfam_row: int = subfam_alignments_collapse[
                    subfam, start_seq_pos
                ][0]
                tr_row: int = prev_subfam_row - (
                    len(alignment_subfamilies) - len(tandem_repeats)
                )
                tr_consensus_changes[changes_index] = (
                    "(" + tr_consensus_matrix[tr_row] + ")n#Simple_repeat"
                )

                # check if TR annotation is from multiple TR rows
                for seq_pos in range(start_seq_pos, stop_seq_pos + 1):
                    cur_subfam_row: int = subfam_alignments_collapse[
                        subfam, seq_pos
                    ][0]
                    if prev_subfam_row != cur_subfam_row:
                        # different TR row used
                        changes_index += 1

                        tr_row = cur_subfam_row - (
                            len(alignment_subfamilies) - len(tandem_repeats)
                        )
                        tr_consensus_changes[changes_index] = (
                            "("
                            + tr_consensus_matrix[tr_row]
                            + ")n#Simple_repeat"
                        )

                        # insert for printing
                        changes_orig.insert(changes_index, "Tandem Repeat")
                        node_ids.insert(
                            changes_index,
                            node_ids[non_empty_columns_orig[cur_changes_pos]],
                        )
                        changes_position_orig.insert(
                            changes_index, cur_changes_pos
                        )
                    cur_changes_pos += 1
                    prev_subfam_row = cur_subfam_row
            changes_index += 1

    # prints results
    if printer.use_matrix_position:
        printer.print_results_matrix(
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )
    elif printer.use_sequence_position:
        printer.print_results_sequence(
            start_all,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )
    else:
        printer.print_results_chrom(
            start_all,
            target.chrom_start,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )

    if printer.use_soda_output:
        consensus_lengths = find_consensus_lengths(alignments)
        printer.print_results_soda(
            start_all,
            target.chrom_name,
            target.chrom_start,
            target.chrom_stop,
            alignment_subfamilies,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            consensus_lengths,
            strand_matrix_collapse,
            consensus_matrix_collapse,
            subfams_collapse_index,
            node_confidence_orig,
            node_ids,
            original_subfamily_sequences,
            original_chromosome_sequences,
            alignment_consensus_start_positions,
            alignment_consensus_stop_positions,
            alignment_start_positions,
            alignment_stop_positions,
            subfam_alignments_collapse,
            support_matrix_collapse,
            subfams_collapse,
            column_count,
        )

    return last_subfam_start, last_subfam_stop

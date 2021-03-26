from logging import Logger
from math import log
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple

from polyA import *


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
    confidence_list, subfams_copy = zip(*sorted(zip(confidence_list, subfams)))

    stdout.write(f"query_label\tconfidence\n")
    for i in range(len(subfams_copy) - 1, 0, -1):
        if confidence_list[i]:
            stdout.write(f"{subfams_copy[i]}\t{confidence_list[i]}\n")


def _validate_target(target: Alignment) -> None:
    if target.chrom_length == 0:
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
    print_seq_pos: bool,
    print_matrix_pos: bool,
) -> Tuple[int, int]:
    """
    If there is only one subfam in the alignment file, no need
    to run anything because we know that subfam is the annotation.
    """
    from uuid import uuid4

    id = uuid4().hex
    if print_seq_pos:
        stdout.write(
            f"{target.start}\t{target.stop}\t{id}\t{target.subfamily}\n"
        )
    elif print_matrix_pos:
        stdout.write(
            f"{1}\t{target.stop - target.start + 1}\t{id}\t{target.subfamily}\n"
        )
    else:
        stdout.write(
            f"{target.start + target.chrom_start}\t{target.stop + target.chrom_start}\t{id}\t{target.subfamily}\n"
        )
    return target.start, target.stop


def _change_probs(seq_count: int) -> Tuple[float, float, float]:
    change_prob_log = log(CHANGE_PROB / (seq_count - 1))

    # jumping in and then out of the skip state counts as 1 jump
    change_prob_skip = change_prob_log / 2

    # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump
    same_prob_skip = change_prob_log / 30

    return change_prob_log, change_prob_skip, same_prob_skip


def run_full(
    alignments: List[Alignment],
    tandem_repeats: List[TandemRepeat],
    chunk_size: int,
    outfile_viz: Optional[TextIO],
    outfile_conf: Optional[TextIO],
    print_matrix_pos: bool,
    print_seq_pos: bool,
    sub_matrix_scores: SubMatrixCollection,
    subfam_counts: Dict[str, float],
    shard_start: int,
    shard_stop: int,
    prev_start: int,
    prev_stop: int,
) -> Tuple[int, int]:
    # chunk start and stop are positions in seq
    seq_count = len(alignments)

    target = alignments[1]
    _validate_target(target)

    if len(tandem_repeats) == 0 and seq_count == 2:
        # only one alignment other than the skip state
        last_subfam_start, last_subfam_stop = _handle_single_alignment(
            target, print_seq_pos, print_matrix_pos
        )
        return last_subfam_start, last_subfam_stop

    change_prob_log, change_prob_skip, same_prob_skip = _change_probs(seq_count)

    # fixme - rename, don't call 'matrix'
    subfamily_matrix: List[str] = []
    chromosome_matrix: List[str] = []
    scores_matrix: List[int] = []
    strands_matrix: List[str] = []
    starts_matrix: List[int] = []
    stops_matrix: List[int] = []
    consensus_starts_matrix: List[int] = []
    consensus_stops_matrix: List[int] = []
    subfamily_sequences_matrix: List[str] = []
    chromosome_sequences_matrix: List[str] = []
    flanks_matrix: List[int] = []
    substitution_matrices: List[SubMatrix] = []
    gap_inits: List[float] = []
    gap_exts: List[float] = []
    lambda_values: List[float] = []

    for alignment in alignments:
        subfamily_matrix.append(alignment.subfamily)
        chromosome_matrix.append(alignment.chrom_name)
        scores_matrix.append(alignment.score)
        strands_matrix.append(alignment.strand)
        starts_matrix.append(alignment.start)
        stops_matrix.append(alignment.stop)
        consensus_starts_matrix.append(alignment.consensus_start)
        consensus_stops_matrix.append(alignment.consensus_stop)
        subfamily_sequences_matrix.append(alignment.subfamily_sequence)
        chromosome_sequences_matrix.append(alignment.sequence)
        flanks_matrix.append(alignment.flank)
        gap_inits.append(alignment.gap_init)
        gap_exts.append(alignment.gap_ext)

        # fixme (george): is this a list of dicts?
        # fixme (george): change it to a dict of dicts?
        if alignment.sub_matrix_name in sub_matrix_scores:
            substitution_matrices.append(
                sub_matrix_scores[alignment.sub_matrix_name]
            )
            lambda_values.append(
                sub_matrix_scores[alignment.sub_matrix_name].lamb
            )
        else:
            substitution_matrices.append(SubMatrix("skip", 0.0))
            lambda_values.append(0.0)

    # precomputes consensus seq length for PrintResultsViz()
    consensus_lengths_matrix: Dict[str, int] = {}
    if outfile_viz:
        for i in range(1, len(flanks_matrix)):
            if strands_matrix[i] == "+":
                consensus_lengths_matrix[subfamily_matrix[i]] = (
                    consensus_stops_matrix[i] + flanks_matrix[i]
                )
            else:
                consensus_lengths_matrix[subfamily_matrix[i]] = (
                    consensus_starts_matrix[i] + flanks_matrix[i]
                )

    tr_consensus_matrix: List[str] = []
    if len(tandem_repeats) > 0:
        for tr in tandem_repeats:
            tr_consensus_matrix.append(tr.consensus)
            # add TR starts and stops
            # check they do not exceed chunk boundary
            if tr.start < shard_start:
                starts_matrix.append(shard_start)
            else:
                starts_matrix.append(tr.start)
            if tr.stop > shard_stop:
                stops_matrix.append(shard_stop)
            else:
                stops_matrix.append(tr.stop)

    # save original alignments before padding for SODA viz
    subfam_alignments = list(subfamily_sequences_matrix)
    chrom_alignments = list(chromosome_sequences_matrix)

    start_all, stop_all = pad_sequences(
        chunk_size,
        starts_matrix,
        stops_matrix,
        subfamily_sequences_matrix,
        chromosome_sequences_matrix,
    )

    # number of rows in matrices
    rows: int = seq_count

    cols, align_matrix = fill_align_matrix(
        lambda_values,
        start_all,
        chunk_size,
        gap_inits,
        gap_exts,
        SKIP_ALIGN_SCORE,
        subfamily_sequences_matrix,
        chromosome_sequences_matrix,
        starts_matrix,
        [sm.scores for sm in substitution_matrices],
    )

    (
        non_empty_columns,
        active_cells,
        consensus_matrix,
    ) = fill_consensus_position_matrix(
        cols,
        rows,
        start_all,
        subfamily_sequences_matrix,
        chromosome_sequences_matrix,
        starts_matrix,
        stops_matrix,
        consensus_starts_matrix,
        strands_matrix,
    )

    # skip state has no active rows
    active_cells[0] = [0]

    # add skip state pad at end
    align_matrix[0, cols] = SKIP_ALIGN_SCORE
    non_empty_columns.append(cols)
    non_empty_columns.append(cols)
    active_cells[cols] = [0]
    active_cells[cols] = [0]
    cols += 1

    repeat_scores: Optional[Dict[int, float]] = None

    if len(tandem_repeats) > 0:
        repeat_scores = calculate_repeat_scores(
            tandem_repeats,
            chunk_size,
            start_all,
            rows,
            active_cells,
            align_matrix,
            consensus_matrix,
            shard_start,
            shard_stop,
        )

        # add skip states for TR cols
        for tr_col in repeat_scores:
            align_matrix[0, tr_col] = SKIP_ALIGN_SCORE

        # add TRs to subfams
        for _ in tandem_repeats:
            subfamily_matrix.append("Tandem Repeat")
            strands_matrix.append("+")
            rows += 1

        confidence_matrix = fill_confidence_matrix_tr(
            non_empty_columns,
            subfam_counts,
            subfamily_matrix,
            active_cells,
            active_cells,
            repeat_scores,
            align_matrix,
        )

        # check if TR columns were added after last alignment
        # TR cols before alignments were accounted for in PadSeqs
        max_col_index: int = max(repeat_scores)
        if max_col_index + 1 > cols:
            cols = max_col_index + 1

        # add TR cols to NonEmptyColumns
        for tr_col in repeat_scores:
            col_set = set(non_empty_columns)  # only add new columns
            col_set.add(tr_col)
            non_empty_columns = list(col_set)
        non_empty_columns.sort()

    else:
        confidence_matrix = fill_confidence_matrix(
            non_empty_columns,
            subfam_counts,
            subfamily_matrix,
            active_cells,
            align_matrix,
        )

    # add skip state to consensus matrix
    # wait till after incase NonEmptyColumns is updated by TR stuff
    for j in non_empty_columns:
        consensus_matrix[0, j] = 0

    support_matrix = fill_support_matrix(
        rows,
        chunk_size,
        start_all,
        non_empty_columns,
        starts_matrix,
        stops_matrix,
        subfamily_matrix,
        confidence_matrix,
        consensus_matrix,
    )

    collapsed_matrices = collapse_matrices(
        rows,
        start_all,
        non_empty_columns,
        subfamily_matrix,
        strands_matrix,
        starts_matrix,
        stops_matrix,
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
    rows = collapsed_matrices.row_num_update
    subfam_alignments_collapse = collapsed_matrices.subfam_alignments_matrix

    align_matrix.clear()
    confidence_matrix.clear()
    support_matrix.clear()
    consensus_matrix.clear()

    if repeat_scores is not None:
        # give different TRs consensus positions that don't allow them to be stitched
        tr_consensus_pos = 1000000
        prev_tr_col = 0
        for tr_col in repeat_scores:
            if prev_tr_col != tr_col - 1:
                tr_consensus_pos -= 500
            consensus_matrix_collapse[rows - 1, tr_col + 1] = tr_consensus_pos
            prev_tr_col = tr_col

    # save original cols for printing heatmap
    cols_orig = cols
    (
        ProbMatrixLastColumn,
        OriginMatrix,
        SameSubfamChangeMatrix,
    ) = fill_probability_matrix(
        same_prob_skip,
        SAME_PROB_LOG,
        change_prob_log,
        change_prob_skip,
        non_empty_columns,
        collapsed_matrices,
    )

    # node_ids for each nucleotide will be assigned during DP backtrace
    node_ids = [""] * cols

    changes_orig: List[str] = []
    changes_position_orig: List[int] = []
    non_empty_columns_orig: List[int] = []

    (changes_position, changes) = get_path(
        non_empty_columns,
        node_ids,
        changes_orig,
        changes_position_orig,
        non_empty_columns_orig,
        subfams_collapse,
        ProbMatrixLastColumn,
        active_cells_collapse,
        OriginMatrix,
        SameSubfamChangeMatrix,
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
            gap_inits,
            gap_exts,
            non_empty_columns,
            starts_matrix,
            stops_matrix,
            changes_position,
            subfamily_matrix,
            subfamily_sequences_matrix,
            chromosome_sequences_matrix,
            subfam_counts,
            substitution_matrices,
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
        j: int = 0
        while j < node_count:
            i: int = 0
            while i < j - 1:
                if path_graph[i * node_count + j] == 1:
                    test = True
                i += 1
            j += 1

        if not test:
            break

        cols = extract_nodes(
            cols, node_count, non_empty_columns, changes_position, path_graph
        )

        # run DP calculations again with nodes corresponding to inserted elements removed
        # ignores removed nodes because they are no longer in NonEmptyColumns
        (
            ProbMatrixLastColumn,
            OriginMatrix,
            SameSubfamChangeMatrix,
        ) = fill_probability_matrix(
            same_prob_skip,
            SAME_PROB_LOG,
            change_prob_log,
            change_prob_skip,
            non_empty_columns,
            collapsed_matrices,
        )

        changes.clear()
        changes_position.clear()

        (changes_position, changes) = get_path(
            non_empty_columns,
            node_ids,
            changes,
            changes_position,
            non_empty_columns,
            subfams_collapse,
            ProbMatrixLastColumn,
            active_cells_collapse,
            OriginMatrix,
            SameSubfamChangeMatrix,
        )
        prev_num_nodes = node_count

    # handles TRs overlapping shards
    i = 0
    # get first subfam in shard
    TR_overlap_index: int = -1
    while i < len(changes_orig) and changes_orig[i] == "skip":
        i += 1

    if i == len(changes_orig):
        # whole shard is skip state
        last_subfam_start: int = -1
        last_subfam_stop: int = -1
    else:
        first_subfam_start = (
            non_empty_columns_orig[changes_position_orig[i]] + start_all - 1
        )
        # check if first print pos is one after prev stop in seq pos
        if prev_stop + 1 == first_subfam_start:
            TR_overlap_index = i
            TR_start_seq_pos: int = prev_start - start_all + 1
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

                if changes_index == TR_overlap_index:
                    # change start pos of first subfam to be prev start seq pos
                    non_empty_columns_orig[cur_changes_pos] = TR_start_seq_pos
                    TR_overlap_index == -1
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
                    len(subfamily_matrix) - len(tandem_repeats)
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
                            len(subfamily_matrix) - len(tandem_repeats)
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
    if print_matrix_pos:
        print_results(
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )
    elif print_seq_pos:
        print_results_sequence(
            start_all,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )
    else:
        print_results_chrom(
            start_all,
            target.chrom_start,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            node_ids,
        )

    if outfile_viz:
        print_results_soda(
            start_all,
            outfile_viz,
            outfile_conf,
            target.chrom_name,
            target.chrom_start,
            target.chrom_stop,
            subfamily_matrix,
            changes_orig,
            tr_consensus_changes,
            changes_position_orig,
            non_empty_columns_orig,
            consensus_lengths_matrix,
            strand_matrix_collapse,
            consensus_matrix_collapse,
            subfams_collapse_index,
            node_confidence_orig,
            node_ids,
            subfam_alignments,
            chrom_alignments,
            subfam_alignments_collapse,
            support_matrix_collapse,
            subfams_collapse,
            cols_orig,
        )

    return last_subfam_start, last_subfam_stop

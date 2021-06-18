from typing import Dict, List, Tuple

from .confidence_cm import confidence_cm
from .performance import timeit
from .substitution_matrix import SubMatrix
from .calculate_score import (
    calculate_score,
    calculate_complexity_adjusted_score,
)
from .sum_repeat_scores import sub_repeat_scores


@timeit()
def fill_node_confidence(
    nodes: int,
    start_all: int,
    gap_inits: List[float],
    gap_exts: List[float],
    columns: List[int],
    starts: List[int],
    stops: List[int],
    changes_position: List[int],
    subfams: List[str],
    subfam_seqs: List[str],
    chrom_seqs: List[str],
    subfam_countss: Dict[str, float],
    sub_matrices: List[SubMatrix],
    repeat_scores: Dict[int, float],
    tr_count: int,
) -> Dict[Tuple[str, int], float]:
    """
    for a particular annotated node, takes all competing alignments and calculates
    the confidence for that node

    first fills matrix with node alignment scores, then reuses matrix for confidence scores
    Using changes_position indentifies boundaries for all nodes, and computes confidence
    values for each node.

    input:
    all input needed for CalcScore() and ConfidenceCM()
    nodes: number of nodes
    changes_pos: node boundaries
    columns: all non empty columns in matrices

    output:
    node_confidence: Hash implementation of sparse 2D matrix that holds confidence values
    for whole nodes. Used during stitching process. Tuple[str, int] is key that maps a subfamily
    and node number to a confidence score.

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> strts = [0, 0, 0]
    >>> stps = [10, 10, 10]
    >>> change_pos = [1, 3, 7, 10]
    >>> names = ["skip", "n1", "n2"]
    >>> s_seqs = ['', 'AAA-TTTT-T', 'TTTTTTTTTT']
    >>> c_seqs = ['', 'TTTTTTTTTT', 'TTTTTTTTTT']
    >>> counts = {"skip": .33, "n1": .33, "n2": .33}
    >>> sub_mat = SubMatrix("", 0.1227)
    >>> sub_mat.scores = {"AA": 10, "AT": -1, "TA": -1, "TT": 10, "..":0}
    >>> sub_mat.background_freqs = None
    >>> sub_mats = [sub_mat] * 3
    >>> rep_scores = {}
    >>> node_conf = fill_node_confidence(3, 0, [0, -25, -25], [0, -5, -5], non_cols, strts, stps, change_pos, names, s_seqs, c_seqs, counts, sub_mats, rep_scores, 0)
    >>> node_conf
    {('n1', 0): 0.19999999999999998, ('n2', 0): 0.7999999999999999, ('n1', 1): 0.058823529411764705, ('n2', 1): 0.9411764705882353, ('n1', 2): 0.1111111111111111, ('n2', 2): 0.8888888888888888}
    >>> s_seqs = ['', 'AAA-T--TT-', 'TTTTTTTTTT']
    >>> c_seqs = ['', 'TTTTTTTTTT', 'TTTTTTTTTT']
    >>> node_conf2 = fill_node_confidence(3, 0, [0, -25, -25], [0, -5, -5], non_cols, strts, stps, change_pos, names, s_seqs, c_seqs, counts, sub_mats, rep_scores, 0)
    >>> node_conf2
    {('n1', 0): 0.19999999999999998, ('n2', 0): 0.7999999999999999, ('n1', 1): 0.001949317738791423, ('n2', 1): 0.9980506822612085, ('n1', 2): 0.1111111111111111, ('n2', 2): 0.8888888888888888}
    """
    node_confidence: Dict[Tuple[str, int], float] = {}
    node_confidence_temp: Dict[Tuple[int, int], float] = {}
    active: Dict[int, List[int]] = {}

    # no need to look at the skip nodes
    node_non_empty: List[int] = []

    # matrix colunms doesn't always equal sequence position because of gaps
    # for each matrix column, compute the gap offset for the sequence
    # sequence position = matrix position + offset
    gap_offset: List[List[int]] = [[] for _ in range(len(subfams))]
    for chrom_index in range(len(chrom_seqs)):
        offset = 0
        for seq_index in range(len(chrom_seqs[chrom_index])):
            if chrom_seqs[chrom_index][seq_index] == "-":
                offset += 1
            else:
                gap_offset[chrom_index].append(offset)

    for node_index in range(nodes):
        begin_node = columns[changes_position[node_index]]
        end_node = columns[changes_position[node_index + 1] - 1]
        range_in_columns = (
            changes_position[node_index + 1] - changes_position[node_index]
        )
        if node_index == nodes - 1:
            # for the last node, include the end column
            end_node += 1

        for subfam_index in range(1, len(subfams) - tr_count):
            subfam_start = starts[subfam_index] - start_all + 1
            subfam_stop = stops[subfam_index] - start_all + 1

            sub_matrix = sub_matrices[subfam_index]
            lamb = sub_matrix.lamb
            gap_init = gap_inits[subfam_index]
            gap_ext = gap_exts[subfam_index]

            if subfam_start <= end_node and subfam_stop >= begin_node:
                # subfam in node, calculate alignment score
                last_prev_subfam = ""
                last_prev_chrom = ""
                alignment_index_start = begin_node - subfam_start

                last_index = 0

                if (
                    0
                    <= alignment_index_start - 1
                    <= len(gap_offset[subfam_index])
                ):
                    chrom_offset = gap_offset[subfam_index][
                        alignment_index_start - 1
                    ]
                    last_prev_subfam = subfam_seqs[subfam_index][
                        alignment_index_start - 1 + chrom_offset
                    ]
                    last_prev_chrom = chrom_seqs[subfam_index][
                        alignment_index_start - 1 + chrom_offset
                    ]

                for j in range(
                    changes_position[node_index],
                    changes_position[node_index] + range_in_columns,
                ):
                    i = columns[j] - begin_node
                    if (
                        0
                        <= alignment_index_start + i
                        < len(gap_offset[subfam_index])
                    ):
                        chrom_offset = gap_offset[subfam_index][
                            alignment_index_start + i
                        ]

                        first_index = alignment_index_start + i + chrom_offset
                        break

                for j in range(
                    changes_position[node_index] + range_in_columns - 1,
                    changes_position[node_index],
                    -1,
                ):
                    i = columns[j] - begin_node
                    if (
                        0
                        <= alignment_index_start + i
                        < len(gap_offset[subfam_index])
                    ):
                        chrom_offset = gap_offset[subfam_index][
                            alignment_index_start + i
                        ]
                        last_index = alignment_index_start + i + chrom_offset
                        break

                # add padding first in case no part of seq is used
                chrom_seq = (
                    "." + chrom_seqs[subfam_index][first_index : last_index + 1]
                )
                subfam_seq = (
                    "."
                    + subfam_seqs[subfam_index][first_index : last_index + 1]
                )

                char_complexity_adjustments = (
                    calculate_complexity_adjusted_score(
                        sub_matrix.background_freqs, subfam_seq, chrom_seq, lamb
                    )
                )

                align_score = lamb * calculate_score(
                    gap_ext,
                    gap_init,
                    subfam_seq,
                    chrom_seq,
                    last_prev_subfam,
                    last_prev_chrom,
                    sub_matrix.scores,
                    char_complexity_adjustments,
                )

                node_confidence_temp[subfam_index, node_index] = align_score

                if node_index in active:
                    active[node_index].append(subfam_index)
                else:
                    active[node_index] = [subfam_index]
                    node_non_empty.append(node_index)

        # TRs
        for subfam_index in range(len(subfams) - tr_count, len(subfams)):
            rep_sum_score: float = 0.0
            tr_start = starts[subfam_index] - start_all + 1
            tr_stop = stops[subfam_index] - start_all + 1

            if not (tr_start > end_node or tr_stop < begin_node):
                rep_sum_score = sub_repeat_scores(
                    begin_node, end_node, repeat_scores
                )

            node_confidence_temp[subfam_index, node_index] = rep_sum_score

            if node_index in active:
                active[node_index].append(subfam_index)
            else:
                active[node_index] = [subfam_index]
                node_non_empty.append(node_index)

    # reuse same matrix and compute confidence scores for the nodes
    for node_index in node_non_empty:
        temp: List[float] = []
        for row_index in active[node_index]:
            temp.append(node_confidence_temp[row_index, node_index])
        confidence_temp: List[float] = confidence_cm(
            temp, subfam_countss, subfams, active[node_index], tr_count, True
        )
        for row_index2 in range(len(confidence_temp)):
            node_confidence_temp[
                active[node_index][row_index2], node_index
            ] = confidence_temp[row_index2]

    # collapse node_confidence down same way supportmatrix is collapsed - all seqs of
    # the same subfam are put in the same row
    for node_index in node_non_empty:
        for row_index in active[node_index]:
            if (subfams[row_index], node_index) in node_confidence:
                node_confidence[
                    subfams[row_index], node_index
                ] += node_confidence_temp[row_index, node_index]
            else:
                node_confidence[
                    subfams[row_index], node_index
                ] = node_confidence_temp[row_index, node_index]

    node_confidence_temp.clear()

    return node_confidence

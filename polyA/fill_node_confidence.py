from typing import Dict, List, Tuple

from polyA.calculate_score import calculate_score
from polyA.confidence_cm import confidence_cm
from sum_repeat_scores import sum_repeat_scores


def fill_node_confidence(
    nodes: int,
    start_all: int,
    gap_init: int,
    gap_ext: int,
    lamb: float,
    infilee: str,
    columns: List[int],
    starts: List[int],
    stops: List[int],
    changes_position: List[int],
    subfams: List[str],
    subfam_seqs: List[str],
    chrom_seqs: List[str],
    subfam_countss: Dict[str, float],
    sub_matrix: Dict[str, int],
    repeat_scores: Dict[int, float],
    tr_count: int,
) -> Dict[Tuple[str, int], float]:
    """
    Using changes_position indentified boundaries for all nodes, and computes confidence
    values for each node. First fills matrix with node alignment scores, then reuses matrix
    for confidence scores.

    input:
    all input needed for CalcScore() and confidence_cm()
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
    >>> change_pos = [0, 3, 7, 10]
    >>> names = ["skip", "n1", "n2"]
    >>> s_seqs = ['', 'AAA-TTTTT-', 'TTTTTTTTTT']
    >>> c_seqs = ['', 'TTTTTTTTTT', 'TTTTTTTTTT']
    >>> counts = {"skip": .33, "n1": .33, "n2": .33}
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> node_conf = fill_node_confidence(3, 0, -25, -5, 0.1227, "infile", non_cols, strts, stps, change_pos, names, s_seqs, c_seqs, counts, sub_mat)
    >>> node_conf
    {('skip', 0): 0.0, ('n1', 0): 0.5, ('n2', 0): 0.5, ('skip', 1): 0.0, ('n1', 1): 0.19999999999999998, ('n2', 1): 0.7999999999999999, ('skip', 2): 0.0, ('n1', 2): 0.19999999999999998, ('n2', 2): 0.7999999999999999}
    """

    node_confidence_temp: List[float] = [
        0.0 for _ in range(len(subfams) * nodes)
    ]
    node_confidence: Dict[Tuple[str, int], float] = {}

    # holds all chrom seqs and the offset needed to get to the correct position in the chrom
    # seq - remember gaps in chrom seq are skipped over for matrix position
    chrom_seq_offset: Dict[int, int] = {}

    # first node
    begin_node0: int = columns[changes_position[0]]

    for subfam_index0 in range(1, len(subfams) - tr_count):

        count: int = 0
        for i in range(
            begin_node0,
            columns[changes_position[1]] - columns[changes_position[0]],
        ):
            if chrom_seqs[subfam_index0][i] == "-":
                count += 1

        chrom_seq_offset[subfam_index0] = count

        end_node0: int = columns[changes_position[1]] + count

        subfam0: str = subfam_seqs[subfam_index0][begin_node0:end_node0]
        chrom0: str = chrom_seqs[subfam_index0][begin_node0:end_node0]

        align_score0: float = 0.0
        # if whole alignment is padding - don't run CalcScore
        if (
            end_node0 >= starts[subfam_index0] - start_all
            and begin_node0 <= stops[subfam_index0] - start_all
        ):
            align_score0 = calculate_score(
                gap_ext, gap_init, subfam0, chrom0, "", "", sub_matrix
            )
        node_confidence_temp[subfam_index0 * nodes + 0] = align_score0

    # first node TRs
    count = 0
    for subfam_index0 in range(len(subfams) - tr_count, len(subfams)):
        rep_sum_score0: float = 0.0
        chrom_seq_offset[subfam_index0] = count  # used for middle nodes

        end_node0 = columns[changes_position[1]] + count

        # get confidence from RepeatScores
        if end_node0 >= starts[subfam_index0] - start_all and begin_node0 <= stops[subfam_index0] - start_all:
            rep_sum_score0 = sum_repeat_scores(begin_node0, end_node0, repeat_scores)
        node_confidence_temp[subfam_index0 * nodes + 0] = rep_sum_score0

    # middle nodes
    for node_index in range(1, nodes - 1):

        for subfam_index in range(1, len(subfams) - tr_count):
            begin_node: int = (
                columns[changes_position[node_index]]
                + chrom_seq_offset[subfam_index]
            )

            count: int = 0
            for i in range(
                begin_node,
                begin_node
                + columns[changes_position[node_index + 1]]
                - columns[changes_position[node_index]],
            ):
                if chrom_seqs[subfam_index][i] == "-":
                    count += 1

            chrom_seq_offset[subfam_index] = (
                chrom_seq_offset[subfam_index] + count
            )

            end_node: int = (
                columns[changes_position[node_index + 1]]
                + chrom_seq_offset[subfam_index]
                + count
            )

            lastprev_subfam: str = subfam_seqs[subfam_index][begin_node - 1]
            lastprev_chrom: str = chrom_seqs[subfam_index][begin_node - 1]
            subfam: str = subfam_seqs[subfam_index][begin_node:end_node]
            chrom: str = chrom_seqs[subfam_index][begin_node:end_node]

            align_score: float = 0.0
            # if whole alignment is padding - don't run CalcScore
            if (
                end_node >= starts[subfam_index] - start_all
                and begin_node <= stops[subfam_index] - start_all
            ):
                align_score = calculate_score(
                    gap_ext,
                    gap_init,
                    subfam,
                    chrom,
                    lastprev_subfam,
                    lastprev_chrom,
                    sub_matrix,
                )

            node_confidence_temp[
                subfam_index * nodes + node_index
            ] = align_score

        # middle nodes TRs
        for subfam_index in range(len(subfams) - tr_count, len(subfams)):
            begin_node: int = columns[changes_position[node_index]] + chrom_seq_offset[subfam_index]
            end_node: int = columns[changes_position[node_index + 1]] + chrom_seq_offset[subfam_index] + count

            rep_sum_score: float = 0.0
            if end_node >= starts[subfam_index] - start_all and begin_node <= stops[subfam_index] - start_all:
                rep_sum_score = sum_repeat_scores(begin_node, end_node, repeat_scores)
            node_confidence_temp[subfam_index * nodes + node_index] = rep_sum_score

    # does last node
    for subfam_index2 in range(1, len(subfams) - tr_count):
        begin_node2: int = (
            columns[changes_position[-2]] + chrom_seq_offset[subfam_index2]
        )

        count: int = 0
        for i in range(
            begin_node2,
            begin_node2
            + columns[changes_position[-1] - 1]
            - columns[changes_position[-2]],
        ):
            if chrom_seqs[subfam_index2][i] == "-":
                count += 1

        chrom_seq_offset[subfam_index2] = (
            chrom_seq_offset[subfam_index2] + count
        )

        end_node2: int = (
            columns[changes_position[-1] - 1]
            + chrom_seq_offset[subfam_index2]
            + count
        )

        lastprev_subfam2: str = subfam_seqs[subfam_index2][begin_node2 - 1]
        lastprev_chrom2: str = chrom_seqs[subfam_index2][begin_node2 - 1]

        subfam2: str = subfam_seqs[subfam_index2][begin_node2 : end_node2 + 1]
        chrom2: str = chrom_seqs[subfam_index2][begin_node2 : end_node2 + 1]

        align_score2: float = 0.0
        # if whole alignment is padding - don't run CalcScore
        if (
            end_node2 >= starts[subfam_index2] - start_all
            and begin_node0 <= stops[subfam_index2] - start_all
        ):
            align_score2 = calculate_score(
                gap_ext,
                gap_init,
                subfam2,
                chrom2,
                lastprev_subfam2,
                lastprev_chrom2,
                sub_matrix,
            )

        node_confidence_temp[subfam_index2 * nodes + nodes - 1] = align_score2

    # last node TRs
    for subfam_index2 in range(len(subfams) - tr_count, len(subfams)):
        begin_node2: int = columns[changes_position[-2]] + chrom_seq_offset[subfam_index2]
        end_node2: int = columns[changes_position[-1] - 1] + chrom_seq_offset[subfam_index2]

        rep_sum_score2: float = 0.0
        if end_node2 >= starts[subfam_index2] - start_all and begin_node0 <= stops[subfam_index2] - start_all:
            rep_sum_score2 = sum_repeat_scores(begin_node2, end_node2 + 1, repeat_scores)
        node_confidence_temp[subfam_index2 * nodes + nodes - 1] = rep_sum_score2

    # reuse same matrix and compute confidence scores for the nodes
    for node_index4 in range(nodes):
        temp: List[float] = []
        for row_index in range(1, len(subfams)):
            temp.append(node_confidence_temp[row_index * nodes + node_index4])

        confidence_temp: List[float] = confidence_cm(
            lamb, infilee, temp, subfam_countss, subfams
        )

        for row_index2 in range(len(confidence_temp)):
            node_confidence_temp[
                (row_index2 + 1) * nodes + node_index4
            ] = confidence_temp[row_index2]

    # collapse node_confidence down same way supportmatrix is collapsed - all seqs of
    # the same subfam are put in the same row
    # not a sparse hash - holds the 0s
    for node_index5 in range(nodes):
        for row_index3 in range(len(subfams)):
            if (subfams[row_index3], node_index5) in node_confidence:
                node_confidence[
                    subfams[row_index3], node_index5
                ] += node_confidence_temp[row_index3 * nodes + node_index5]
            else:
                node_confidence[
                    subfams[row_index3], node_index5
                ] = node_confidence_temp[row_index3 * nodes + node_index5]

    return node_confidence

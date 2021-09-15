from typing import Dict, List, Tuple

from .performance import timeit
from .ultra_provider import TandemRepeat


@timeit()
def calculate_repeat_scores(
    tandem_repeats: List[TandemRepeat],
    chunk_size: int,
    start_all: int,
    row_num: int,
    active_cells: Dict[int, List[int]],
    align_matrix: Dict[Tuple[int, int], float],
    consensus_matrix: Dict[Tuple[int, int], int],
    shard_start: int,
    shard_stop: int,
) -> Dict[int, float]:
    """
    Calculates score (according to ULTRA scoring) for every segment of size
    chunksize for all tandem repeats found in the target sequence.
    Scores are of the surrounding chunksize nucleotides in the sequence.

    At same time, updates AlignMatrix, ActiveCells, and ConsensusMatrix
    to include tandem repeat scores.

    input:
    tandem_repeats: list of tandem repeats from ULTRA output
    chunk_size: size of nucleotide chunks that are scored
    start_all: minimum starting nucleotide position
    row_num: number of rows in matrices
    columns: all columns that are not empty from alignment scores
    active_cells: dictionary that maps column number to a list of active rows
    align_matrix
    consensus_matrix: dictionary that maps a tuple (row, col) to consensus
    position in the subfam/query sequence
    shard_start: start seq pos of current shard
    shard_stop: stop seq pos of current shard

    output:
    repeat_scores: Hash implementation of sparse 1d array. Key is int that maps col
    index in the target sequence to it's tandem repeat score. Used in FillNodeConfidence.

    tr_info - list of dictionaries - [{'Start': num, 'Length': num, 'PositionScoreDelta': '-0.014500:1.27896'}, ...]
    >>> repeats = [TandemRepeat.from_json(m) for m in [{'Start': 5, 'Length': 4, 'PositionScoreDelta': '-0.5:1:1.5:1', 'Consensus': 'AT'}, {'Start': 10, 'Length': 8, 'PositionScoreDelta': '0:0.5:0.5:1.5:1.5:1:0.5:-0.5', 'Consensus': 'GTT'}]]
    >>> align_mat = {}
    >>> active_cols = {}
    >>> rep_scores = calculate_repeat_scores(repeats, 5, 4, 1, active_cols, align_mat, {}, 1, 30)
    >>> rep_scores
    {2: -0.5, 3: 1.0, 4: 1.5, 5: 1.0, 7: 0.0, 8: 0.5, 9: 0.5, 10: 1.5, 11: 1.5, 12: 1.0, 13: 0.5, 14: -0.5}
    >>> align_mat
    {(1, 2): 3.3333333333333335, (1, 3): 3.75, (1, 4): 3.75, (1, 5): 5.833333333333333, (2, 7): 1.6666666666666667, (2, 8): 3.125, (2, 9): 4.0, (2, 10): 5.0, (2, 11): 5.0, (2, 12): 4.0, (2, 13): 3.125, (2, 14): 1.6666666666666667}
    >>> active_cols
    {2: [0, 1], 3: [0, 1], 4: [0, 1], 5: [0, 1], 7: [0, 2], 8: [0, 2], 9: [0, 2], 10: [0, 2], 11: [0, 2], 12: [0, 2], 13: [0, 2], 14: [0, 2]}

    Shard cuts off TR from left side
    >>> repeats = [TandemRepeat.from_json(m) for m in [{'Start': 5, 'Length': 4, 'PositionScoreDelta': '-0.5:1:1.5:1', 'Consensus': 'AT'}, {'Start': 10, 'Length': 8, 'PositionScoreDelta': '0:0.5:0.5:1.5:1.5:1:0.5:-0.5', 'Consensus': 'GTT'}]]
    >>> align_mat = {}
    >>> active_cols = {}
    >>> rep_scores = calculate_repeat_scores(repeats, 5, 4, 1, active_cols, align_mat, {}, 6, 30)
    >>> rep_scores
    {3: 1.0, 4: 1.5, 5: 1.0, 7: 0.0, 8: 0.5, 9: 0.5, 10: 1.5, 11: 1.5, 12: 1.0, 13: 0.5, 14: -0.5}
    >>> align_mat
    {(1, 3): 3.75, (1, 4): 3.75, (1, 5): 5.833333333333333, (2, 7): 1.6666666666666667, (2, 8): 3.125, (2, 9): 4.0, (2, 10): 5.0, (2, 11): 5.0, (2, 12): 4.0, (2, 13): 3.125, (2, 14): 1.6666666666666667}
    >>> active_cols
    {3: [0, 1], 4: [0, 1], 5: [0, 1], 7: [0, 2], 8: [0, 2], 9: [0, 2], 10: [0, 2], 11: [0, 2], 12: [0, 2], 13: [0, 2], 14: [0, 2]}

    Shard cuts off TR from right side in first half of chunk
    >>> repeats = [TandemRepeat.from_json(m) for m in [{'Start': 5, 'Length': 4, 'PositionScoreDelta': '-0.5:1:1.5:1', 'Consensus': 'AT'}, {'Start': 10, 'Length': 8, 'PositionScoreDelta': '0:0.5:0.5:1.5:1.5:1:0.5:-0.5', 'Consensus': 'GTT'}]]
    >>> align_mat = {}
    >>> active_cols = {}
    >>> rep_scores = calculate_repeat_scores(repeats, 5, 4, 1, active_cols, align_mat, {}, 1, 11)
    >>> rep_scores
    {2: -0.5, 3: 1.0, 4: 1.5, 5: 1.0, 7: 0.0, 8: 0.5}
    >>> align_mat
    {(1, 2): 3.3333333333333335, (1, 3): 3.75, (1, 4): 3.75, (1, 5): 5.833333333333333, (2, 7): 1.6666666666666667, (2, 8): 3.125}
    >>> active_cols
    {2: [0, 1], 3: [0, 1], 4: [0, 1], 5: [0, 1], 7: [0, 2], 8: [0, 2]}
    """

    repeat_scores: Dict[
        int, float
    ] = {}  # maps col in target seq to ultra output score
    # for each tandem repeat, get overlapping windowed score
    for i in range(len(tandem_repeats)):
        # get repeat info
        rep = tandem_repeats[i]
        col_index = rep.start - start_all + 1  # can't be less than 1
        length = rep.length
        pos_scores = rep.position_scores

        # calc score for first chunk
        score: float = 0
        j = 0
        k = int((chunk_size - 1) / 2)
        # update repeat_scores for positions 0 to 15 if in boundary
        while j <= k and j < length:
            score += pos_scores[j]
            if shard_start <= rep.start + j <= shard_stop:
                repeat_scores[col_index + j] = pos_scores[j]
            j += 1
        window_size = j

        # check TR start >= chunk start
        # update AlignMatrix and ConsensusMatrix
        if rep.start >= shard_start:
            align_matrix[i + row_num, col_index] = (
                score * chunk_size / window_size
            )
            consensus_matrix[i + row_num, col_index] = rep.start
            # update active_cells
            if col_index in active_cells:
                active_cells[col_index].append(i + row_num)
            else:
                active_cells[col_index] = [0, i + row_num]

        # calc score for the rest of the chunks
        tr_chunk_length = length
        if rep.stop >= shard_stop:
            # TR ends outside of chunk boundary
            tr_chunk_length = shard_stop - rep.start + 1
        for j in range(1, tr_chunk_length):
            # check to remove last score in window
            if j - k - 1 >= 0:
                score -= float(pos_scores[j - k - 1])
                window_size -= 1
            # check to add new score in window
            if j + k < len(pos_scores):
                score += float(pos_scores[j + k])
                window_size += 1
                if shard_start <= rep.start + j + k <= shard_stop:
                    # add new score to repeat_scores
                    repeat_scores[col_index + j + k] = float(pos_scores[j + k])
            # Update matrices if in chunk boundary
            if rep.start + j >= shard_start:
                align_matrix[i + row_num, col_index + j] = (
                    score * chunk_size / window_size
                )
                consensus_matrix[i + row_num, col_index + j] = rep.start + j
                # add active cells
                if col_index + j in active_cells:
                    active_cells[col_index + j].append(i + row_num)
                else:
                    active_cells[col_index + j] = [0, i + row_num]
    return repeat_scores

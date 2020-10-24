from typing import Dict, List, Tuple


def CalcRepeatScores(
    tandem_repeats,
    chunk_size: int,
    start_all: int,
    row_num: int,
    active_cells: Dict[int, List[int]],
    align_matrix: Dict[Tuple[int, int], float],
    consensus_matrix: Dict[Tuple[int, int], int],
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

    output:
    repeat_scores: Hash implementation of sparse 1d array. Key is int that maps col
    index in the target sequence to it's tandem repeat score. Used in FillNodeConfidence.

    tr_info - list of dictionaries - [{'Start': num, 'Length': num, 'PositionScoreDelta': '-0.014500:1.27896'}, {}]
    >>> repeats = [{'Start': 5, 'Length': 4, 'PositionScoreDelta': '-0.5:1:1.5:1'}, {'Start': 10, 'Length': 8, 'PositionScoreDelta': '0:0.5:0.5:1.5:1.5:1:0.5:-0.5'}]
    >>> align_mat = {}
    >>> active_cols = {}
    >>> rep_scores = CalcRepeatScores(repeats, 5, 4, 1, active_cols, align_mat, {})
    >>> align_mat
    {(1, 2): 3.3333333333333335, (1, 3): 3.75, (1, 4): 3.75, (1, 5): 5.833333333333333, (2, 7): 1.6666666666666667, (2, 8): 3.125, (2, 9): 4.0, (2, 10): 5.0, (2, 11): 5.0, (2, 12): 4.0, (2, 13): 3.125, (2, 14): 1.6666666666666667}
    >>> active_cols
    {2: [0, 1], 3: [0, 1], 4: [0, 1], 5: [0, 1], 7: [0, 2], 8: [0, 2], 9: [0, 2], 10: [0, 2], 11: [0, 2], 12: [0, 2], 13: [0, 2], 14: [0, 2]}
    """

    repeat_scores: Dict[
        int, float
    ] = {}  # maps col in target seq to ultra output score
    # for each tandem repeat, get overlapping windowed score
    for i in range(len(tandem_repeats)):
        # get repeat info
        rep = tandem_repeats[i]
        start_rep = rep["Start"] + 1
        col_index = start_rep - (start_all)
        length = rep["Length"]
        pos_scores = rep["PositionScoreDelta"].split(":")

        # calc score for first chunk
        score: float = 0
        j = 0
        k = int((chunk_size - 1) / 2)
        # gets 0 to 15
        while j <= k and j < length:
            score += float(pos_scores[j])
            repeat_scores[col_index + j] = float(pos_scores[j])
            j += 1
        window_size = j

        # update AlignMatrix and ConsensusMatrix
        align_matrix[i + row_num, col_index] = score * chunk_size / window_size
        consensus_matrix[i + row_num, col_index] = start_rep

        # update active_cells
        if col_index in active_cells:
            active_cells[col_index].append(i + row_num)
        else:
            active_cells[col_index] = [0, i + row_num]

        # calc score for the rest of the chunks
        for j in range(1, length):
            # check to remove last score in window
            if j - k - 1 >= 0:
                score -= float(pos_scores[j - k - 1])
                window_size -= 1
            # check to add new score in window
            if j + k < len(pos_scores):
                score += float(pos_scores[j + k])
                window_size += 1
                # add new score to repeat_scores
                repeat_scores[col_index + j + k] = float(pos_scores[j + k])

            align_matrix[i + row_num, col_index + j] = (
                score * chunk_size / window_size
            )
            consensus_matrix[i + row_num, col_index + j] = start_rep + j
            # add active cells
            if col_index + j in active_cells:
                active_cells[col_index + j].append(i + row_num)
            else:
                active_cells[col_index + j] = [0, i + row_num]

    return repeat_scores

from typing import Dict, List, Tuple

from polyA.matrices import ConfidenceMatrix, SupportMatrix, ConsensusMatrix


def fill_support_matrix(
    row_num: int,
    chunk_size: int,
    start_all: int,
    columns: List[int],
    starts: List[int],
    stops: List[int],
    subfams: List[str],
    confidence_matrix: ConfidenceMatrix,
    consensus_matrix: ConsensusMatrix,
) -> SupportMatrix:
    """
    Fills support_matrix using values in confidence_matrix. Average confidence values
    for surrounding chunk_size confidence values - normalized by dividing by number of
    segments.

    score for subfam x at position i is sum of all confidences for subfam x for all segments that
    overlap position i - divided by number of segments

    input:
    row_num: number of rows in matrices
    chunk_size: number of nucleotides that make up a chunk
    start_all: min start position on chromosome/target sequences for whole alignment
    columns: list of non empty columns in matrices
    starts: chrom/target start positions of all alignments
    stops: chrom/target stop positions of all alignments
    confidence_matrix: Hash implementation of sparse 2D matrix that holds confidence values

    output:
    support_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
    tuple[int, int] that maps row, col with the value held in that cell of matrix. Rows are
    subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
    Each cell in matrix is the support score (or average confidence value) for the surrounding
    chunk_size cells in the confidence matrix.

    >>> non_cols = [0,1,2]
    >>> strts = [0, 0]
    >>> stps = [0, 1]
    >>> subs = ['subfam1', 'subfam2']
    >>> conf_mat = {(0, 0): 0.9, (0, 1): 0.5, (0, 2): .5, (1, 0): 0.1, (1, 1): .3, (1, 2): .1}
    >>> cons_mat = {(0, 0): 1, (0, 1): 1, (0, 2): 1, (1, 0): 1, (1, 1): 1, (1, 2): 1}
    >>> fill_support_matrix(2, 3, 0, non_cols, strts, stps, subs, conf_mat, cons_mat)
    {(0, 0): 0.7, (0, 1): 0.6333333333333333, (0, 2): 0.5, (1, 0): 0.1, (1, 1): 0.3, (1, 2): 0.1}
    """

    support_matrix: Dict[Tuple[int, int], float] = {}

    half_chunk: int = int((chunk_size - 1) / 2)

    # skip state
    for col in range(len(columns)):
        col_index: int = columns[col]

        summ: float = 0
        num_segments: int = 0

        for sum_index in range(
            col_index - half_chunk, col_index + half_chunk + 1
        ):
            if (0, sum_index) in consensus_matrix:
                num_segments += 1
                summ += confidence_matrix[0, sum_index]

        support_matrix[0, col_index] = summ / num_segments

    # rest of rows
    for row_index in range(1, row_num):

        start: int = starts[row_index] - start_all + 1
        stop: int = stops[row_index] - start_all + 1

        # if the alignment is small, do it the easy way to avoid out of bound errors
        if (stop - start + 1 < chunk_size) or subfams[
            row_index
        ] == "Tandem Repeat":
            for col in range(len(columns)):
                j = columns[col]

                if (row_index, j) in consensus_matrix:
                    num: int = j
                    summ: float = 0.0
                    numsegments: int = 0
                    while num >= 0 and num >= j:
                        if (row_index, num) in consensus_matrix:
                            summ = summ + confidence_matrix[row_index, num]
                            numsegments += 1
                        num -= 1

                    if numsegments > 0:
                        support_matrix[row_index, j] = summ / numsegments

        # if the alignment is large, do it the fast way
        else:
            # first chunk_size/2 before getting to full chunks
            left_index: int = 0
            for col_index in range(
                start, starts[row_index] + 1 - start_all + half_chunk
            ):
                summ: float = 0.0
                sum_index: int = col_index - left_index
                num_segments: int = 0

                while (
                    sum_index <= col_index + half_chunk
                    and sum_index < columns[-1]
                ):
                    summ += confidence_matrix[row_index, sum_index]
                    sum_index += 1
                    num_segments += 1

                support_matrix[row_index, col_index] = summ / num_segments
                left_index += 1

            # middle part where num segments is chunk_size
            col_index: int = start + half_chunk
            summ: float = 0.0
            sum_index: int = col_index - half_chunk
            while sum_index <= col_index + half_chunk:
                summ += confidence_matrix[row_index, sum_index]
                sum_index += 1

            support_matrix[row_index, col_index] = summ / chunk_size

            # to get next bp average subtract previous confidence, add next confidence
            for col_index in range(
                (1 + start + half_chunk), (stop + 1) - half_chunk
            ):
                last_index: int = col_index - half_chunk - 1
                next_index: int = col_index + half_chunk

                summ = (
                    summ
                    - confidence_matrix[row_index, last_index]
                    + confidence_matrix[row_index, next_index]
                )

                # rounding error when summ gets really small sometimes causes summ to become negative
                # when this happen recompute sum the naiive way
                if summ <= 0:
                    summ = 0
                    sum_index = last_index + 1
                    while sum_index <= next_index:
                        summ += confidence_matrix[row_index, sum_index]
                        sum_index += 1

                support_matrix[row_index, col_index] = summ / chunk_size

            # last chunk_size/2
            right_index: int = half_chunk
            for col_index in range(
                (stop + 1) - half_chunk, stops[row_index] + 1 - start_all + 1
            ):

                summ: float = 0.0
                sum_index: int = col_index - half_chunk
                num_segments: int = 0

                while sum_index < col_index + right_index:
                    summ += confidence_matrix[row_index, sum_index]
                    sum_index += 1
                    num_segments += 1

                support_matrix[row_index, col_index] = summ / num_segments
                right_index -= 1

    return support_matrix

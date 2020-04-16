from typing import Dict, Iterable, Tuple

from .alignments import Alignment
from .constants import DEFAULT_CHUNK_SIZE, DEFAULT_GAP_EXT, DEFAULT_GAP_INIT

AlignScoreMatrix = Dict[Tuple[int, int], float]


def fill_align_score_matrix(
    alignments: Iterable[Alignment], chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> AlignScoreMatrix:
    alignScoreMatrix: AlignScoreMatrix = {}

    for alignmentIndex, alignment in enumerate(alignments):
        seq = alignment.sequence
        seqLen = len(seq)
        subSeq = alignment.subfamily_sequence

        indexI = 0
        indexJ = 0
        offset = chunk_size

        while indexJ + offset < seqLen:
            innerIndex = indexJ
            innerCount = 0

            while innerCount < chunk_size:
                if seq[innerIndex] != "-":
                    innerCount += 1
                innerIndex += 1

            offset = innerIndex - indexJ

            if seq[indexJ] != "-":
                seqSlice = seq[indexJ : indexJ + offset - 1]
                subSeqSlice = subSeq[indexJ : indexJ + offset - 1]

                alignScore: float
                if seqSlice != "." and subSeqSlice[0] != ".":
                    # TODO: Flip these calculate_score arguments around? Does it matter?
                    if indexJ == 0:
                        alignScore = calculate_score(
                            seqSlice, subSeqSlice, "", "",
                        )
                    else:
                        alignScore = calculate_score(
                            seqSlice,
                            subSeqSlice,
                            seq[indexJ - 1],
                            subSeq[indexJ - 1],
                        )

                    alignScoreMatrix[(alignmentIndex, indexI)] = alignScore

                indexI += 1

            indexJ += 1

    return alignScoreMatrix

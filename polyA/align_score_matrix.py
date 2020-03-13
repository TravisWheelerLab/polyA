from typing import Dict, Iterable, Tuple

from .alignments import Alignment
from .constants import DEFAULT_CHUNK_SIZE, DEFAULT_GAP_EXT, DEFAULT_GAP_INIT

AlignScoreMatrix = Dict[Tuple[int, int], float]


# TODO: Make this work or whatever
def calculate_score(
    seq1: str,
    seq2: str,
    seq1Prev: str,
    seq2Prev: str,
    gapExt: int = DEFAULT_GAP_EXT,
    gapInit: int = DEFAULT_GAP_INIT,
) -> float:
    chunkScore = 0.0

    if seq1[0] == "-":
        if seq1Prev == "-":
            # We are in the midst of an existing gap
            chunkScore += gapExt
        else:
            # We have just entered a brand new gap
            chunkScore += gapInit
    elif seq2[0] == "-":
        if seq2Prev == "-":
            # We are in the midst of an existing gap
            chunkScore += gapExt
        else:
            # We have just entered a brand new gap
            chunkScore += gapInit
    elif seq1[0] == "." or seq2[0] == ".":
        # There is nothing at this point in the alignment so the score stays
        # the same.
        pass
    else:
        chunkScore += 0.0  # TODO: Bunch of weird stuff here

    # TODO: The for-loop, but can we merge the inside bits with the above?

    return chunkScore


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

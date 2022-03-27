from typing import Dict

from .alignment import Alignment


def remove_cg_scores(a: Alignment, sub_matrix: Dict[str, int]) -> int:
    # remove CG score from a.score and return
    cg_score = 0
    consensus = a.sequences[1]
    target = a.sequences[0]
    i = 0
    while i < len(consensus) - 1:
        offset = 2
        # check if it contains 2 chars
        consensus_char_window = consensus[i : i + offset]
        while consensus_char_window[-1] == "-":
            offset += 1
            consensus_char_window = consensus[i : i + offset]
        if consensus_char_window.replace("-", "") == "CG":
            target_char_window = target[i : i + offset]
            for char_index in range(len(consensus_char_window)):
                if (
                    consensus_char_window[char_index] != "-"
                    and target_char_window[char_index] != "-"
                ):
                    # remove scoring
                    # seq1=subfam, seq2=chrom
                    # sub_matrix[seq1_char + seq2_char]
                    cg_score += sub_matrix[
                        consensus_char_window[char_index]
                        + target_char_window[char_index]
                    ]
        i += offset - 1
    # return score - cg content scores
    return a.score - cg_score

from typing import Tuple
from subprocess import Popen, PIPE
from collections import Counter
from typing import Dict, List, Optional, TextIO, Tuple, Iterable
from .confidence_cm import confidence_only, confidence_subfam_pairs
from .alignment import Alignment
from os import remove


MERGE_CONF_THRESH = 0.5


def merge_subfams(subfam_A: str, subfam_B: str, subfam_instances_path):
    merged_name = subfam_A + "_" + subfam_B
    subfam_A_instances = subfam_instances_path + subfam_A + ".fa"
    subfam_B_instances = subfam_instances_path + subfam_B + ".fa"
    merged_subfam_instances = subfam_instances_path + merged_name + ".fa"

    # create file with all subfam instances
    outfile = open(merged_subfam_instances, "w")
    with open(subfam_A_instances, "r") as infile:
        outfile.write(infile.read())
    with open(subfam_B_instances, "r") as infile:
        outfile.write(infile.read())
    outfile.close()

    # create MSA from the instances file
    # TODO: eventually use refiner for this
    merged_subfam_msa = subfam_instances_path + merged_name + ".afa"
    with open(merged_subfam_msa, "w") as f_out_merged_msa:
        process = Popen(
            ["mafft", merged_subfam_instances],
            stdout=f_out_merged_msa,
            stderr=PIPE,
        )
        process.communicate()

    hmm_out = subfam_instances_path + "merged_subfam.hmm"
    process = Popen(
        ["hmmbuild", hmm_out, merged_subfam_msa], stdout=PIPE, stderr=PIPE
    )
    process.communicate()
    process = Popen(["hmmemit", "-c", hmm_out], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    merged_consensus_seq = (
        str(stdout).split("consensus")[1].replace("\\n", "\n")
    )
    merged_consensus_seq = merged_consensus_seq[: len(merged_consensus_seq) - 1]
    remove(hmm_out)
    return merged_consensus_seq, merged_name


def test_seq_confidence(
    scores,
    lambs,
    subfams,
    subfam_winners,
    uncertain_subfam_pairs,
):
    # test_seq_confidence updates uncertain_subfam_pairs
    confidence_list = confidence_only(scores, lambs)
    confidence_list, subfams = zip(*sorted(zip(confidence_list, subfams)))

    # check for a clear winner
    if confidence_list[-1] > 0.75:
        subfam_winners[subfams[-1]] += 1
    else:
        uncertainity_pair_thresh = confidence_list[-1] / 3
        # look for uncertain pairs
        for i in range(len(subfams) - 1, 0, -1):
            if confidence_list[i] < uncertainity_pair_thresh:
                break
            for j in range(i - 1, 0, -1):
                if confidence_list[j] < uncertainity_pair_thresh:
                    break
                # otherwise count subfam pair
                # a test seq could have multiple alignments to the same subfam
                if subfams[i] != subfams[j]:
                    sub_pair = [subfams[i], subfams[j]]
                    sub_pair.sort()
                    uncertain_subfam_pairs[tuple(sub_pair)] += 1


def subfam_confidence(
    alignments: List[Alignment],
    lambs: List[float],
    subfam_instances_path: str,
    merge_conf_thresh: float,
):
    subfams: List[str] = []
    scores: List[int] = []
    subfam_lambs: List[float] = []

    prev_test_seq_name: str = ""
    test_seqs: int = 0
    uncertain_subfam_pairs: Dict[Tuple[str, str], int] = Counter()
    subfam_winners: Dict[str, int] = Counter()

    for i, a in enumerate(alignments):
        # determine if chrom_name has changed, then add to list
        test_seq_name = (
            a.chrom_name + ":" + str(a.chrom_start) + "-" + str(a.chrom_stop)
        )
        if test_seq_name != prev_test_seq_name and len(subfams) != 0:
            test_seqs += 1
            # run confidence_only on prev input
            test_seq_confidence(
                scores,
                lambs,
                subfams,
                subfam_winners,
                uncertain_subfam_pairs,
            )

            # create new lists for new test seq
            subfams = []
            scores = []
            subfam_lambs = []
        prev_test_seq_name = test_seq_name
        subfams.append(a.subfamily)
        scores.append(a.score)
        subfam_lambs.append(lambs[i])

    # last test seq
    test_seq_confidence(
        scores,
        lambs,
        subfams,
        subfam_winners,
        uncertain_subfam_pairs,
    )

    subfam_pair_confidence, zero_conf_subfams = confidence_subfam_pairs(
        uncertain_subfam_pairs, subfam_winners
    )

    # sort by subfam i with highest number of uncertain pairs with j
    # ex: AluYb8 uncertain with AluYb9 403 times
    sorted_zero = sorted(
        zero_conf_subfams.items(),
        key=lambda item: sorted(
            item[1].items(), key=lambda item: item[1], reverse=True
        )[0][1],
        reverse=True,
    )

    merged_consensus = ""
    merged_name = ""
    sub_pair: Tuple[str, str] = ("", "")

    # check to merge zero confidence pairs first
    for zero_conf_item in sorted_zero:
        zero_conf_highest_pair = sorted(
            zero_conf_item[1].items(), key=lambda item: item[1], reverse=True
        )[0]
        sub_pair = (zero_conf_item[0], zero_conf_highest_pair[0])
        merged_consensus, merged_name = merge_subfams(
            zero_conf_item[0], zero_conf_highest_pair[0], subfam_instances_path
        )
        return merged_consensus, merged_name, sub_pair

    # merge if confidence value is under some threshold
    sorted_pairs = sorted(
        subfam_pair_confidence.items(), key=lambda item: item[1]
    )

    # (('AluYk2#SINE/Alu', 'AluY#SINE/Alu'), 0.00558659217877095)
    # merge lowest conf pair first
    if len(sorted_pairs) != 0 and sorted_pairs[0][1] < merge_conf_thresh:
        sub_pair = sorted_pairs[0][0]
        merged_consensus, merged_name = merge_subfams(
            sub_pair[0], sub_pair[1], subfam_instances_path
        )
    # return values will be empty if no pairs to merge
    return merged_consensus, merged_name, sub_pair

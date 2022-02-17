from subprocess import Popen, PIPE
from collections import Counter
from typing import Dict, List, Tuple, Any
from .confidence_cm import confidence_only
from .alignment import Alignment
from os import remove


MERGE_CONF_THRESH = 0.5


def test_seq_confidence(
    scores: List[int],
    lambs: List[float],
    subfams: List[str],
    subfam_winners: Dict[str, int],
    uncertain_subfam_pairs: Dict[Tuple[str, str], int],
    winner_group_count: Dict[str, int],
):
    """
    Finds subfamilies that are a potential winner of a test sequence.
    Potential winners have confidence values for the test sequence
    above the max confidence value * k, where k is set to a default
    of 1/3.

    input:
    scores: list of alignment scores for competing annotations of the test seq
    lambs: list of lambda values for score matrix used
    subfams: list of subfam names
    subfam_winners: dictionary that maps a subfam to its clear winner count
    uncertain_subfam_pairs: dictionary that maps a subfam pair to the number
    of times they are both found in a potential winner group
    winner_group_count: dictionary that maps a subfam to its potential winner
    count
    """
    # calc confidence values of subfam alignments for test seq
    confidence_list = confidence_only(scores, lambs)
    confidence_list, subfams = zip(*sorted(zip(confidence_list, subfams)))

    # subfams with conf values >= to max conf / 3
    # are potential winners of this test seq
    potential_winner_thresh = confidence_list[-1] / 3

    # search for subfams above the potential winner threshold
    subfam_winner_group_count = 0
    for i in range(len(subfams) - 1, 0, -1):
        if confidence_list[i] < potential_winner_thresh:
            break
        # track count of potential subfam winners
        subfam_winner_group_count += 1
        winner_group_count[subfams[i]] += 1
        # create pairs of subfams that are both potential winners for the test seq
        # these are considered to be uncertain pairs
        for j in range(i - 1, 0, -1):
            if confidence_list[j] < potential_winner_thresh:
                break
            # single test seq could have 1+ alignments with the same subfam
            if subfams[i] != subfams[j]:
                sub_pair = [subfams[i], subfams[j]]
                sub_pair.sort()
                uncertain_subfam_pairs[tuple(sub_pair)] += 1

    if subfam_winner_group_count == 1:
        # single subfam found above winner thresh is a clear
        # winner of this test seq
        subfam_winners[subfams[-1]] += 1


def confidence_subfam_pairs(
    uncertain_subfam_pair_counts: Dict[Tuple[str, str], int],
    subfam_winner_counts: Dict[str, int],
) -> Tuple[Dict[Tuple[str, str], float], Dict[Any, Any]]:
    """
    Compute confidence values for uncertain subfamily pairs.
    The confidence value is a measure of how certain we are
    that a pair of subfamilies are distinct enough to produce
    reliable annotation. The confidence of an ij-pair is
    calculated by:
        clear_winner_count_i / (uncertain_count_ij + clear_winner_count_i)
    where clear_winner_count_i is the number of times subfam i was
    a clear winner over all test seqs and uncertain_count_ij is the
    number of times subfam i and j were found to be in an uncertain pair.

    input:
    uncertain_subfam_pair_counts: dictionary that maps a subfam pair to the number
    of times they are both found in a potential winner group
    subfam_winner_counts: dictionary that maps a subfam to its clear winner count

    output:
    subfam_pair_confidence: dictionary that maps a subfam pair to
    its confidence value
    zero_conf_subfams: dictionary that maps a subfam to subfams
    it had zero confidence with and the number of occurrences
    """
    subfam_pair_confidence: Dict[Tuple[str, str], float] = {}
    zero_conf_subfams: Dict[str, Any] = {}

    # compute subfam ij- and ji-pair confidence
    for sub_pair, sub_pair_count in uncertain_subfam_pair_counts.items():
        subfam_i = sub_pair[0]
        subfam_j = sub_pair[1]
        # subfam ij-pair
        if subfam_i in subfam_winner_counts.keys():
            # compute subfam pair confidence
            subfam_pair_confidence[(subfam_i, subfam_j)] = subfam_winner_counts[
                subfam_i
            ] / (sub_pair_count + subfam_winner_counts[subfam_i])
        else:
            # subfam i was never a clear winner
            if subfam_i not in zero_conf_subfams.keys():
                zero_conf_subfams[subfam_i] = {}
            zero_conf_subfams[subfam_i][subfam_j] = sub_pair_count
        # subfam ji-pair
        if subfam_j in subfam_winner_counts.keys():
            # compute subfam pair confidence
            subfam_pair_confidence[(subfam_j, subfam_i)] = subfam_winner_counts[
                subfam_j
            ] / (sub_pair_count + subfam_winner_counts[subfam_j])
        else:
            # subfam j was never a clear winner
            if subfam_j not in zero_conf_subfams.keys():
                zero_conf_subfams[subfam_j] = {}
            zero_conf_subfams[subfam_j][subfam_i] = sub_pair_count

    return subfam_pair_confidence, zero_conf_subfams


def merge_subfams(
    subfam_i: str,
    subfam_j: str,
    subfam_instances_path: str,
    subfam_to_merged_num: Dict[str, int],
) -> Tuple[str, str]:
    """
    Merges subfamilies i and j and creates a new
    instance file and consensus sequence for the
    merged subfamily.

    input:
    subfam_i: name of subfam i
    subfam_j: name of subfam j
    subfam_instances_path: path to a dir of subfam instances
    subfam_to_merged_num: dictionary that maps names of prior
    merged subfamilies to a number for instance file look up

    output:
    merged_consensus_seq: consensus sequence of the merged subfam
    merged_subfam_name: name of the merged subfam
    """
    merged_subfam_name = subfam_i + "_" + subfam_j

    # check if subfams i and j are merged subfams
    if subfam_i in subfam_to_merged_num.keys():
        subfam_i = "merged_" + str(subfam_to_merged_num[subfam_i])
    if subfam_j in subfam_to_merged_num.keys():
        subfam_j = "merged_" + str(subfam_to_merged_num[subfam_j])

    # paths to subfam i and j instance files
    subfam_i_instances = subfam_instances_path + subfam_i + ".fa"
    subfam_j_instances = subfam_instances_path + subfam_j + ".fa"

    # path to merged instance file
    merged_subfam_filename = "merged_" + str(len(subfam_to_merged_num) + 1)
    merged_subfam_instances = (
        subfam_instances_path + merged_subfam_filename + ".fa"
    )

    # output fasta file with all subfam i and j instances
    outfile = open(merged_subfam_instances, "w")
    with open(subfam_i_instances, "r") as infile:
        outfile.write(infile.read())
    with open(subfam_j_instances, "r") as infile:
        outfile.write(infile.read())
    outfile.close()

    # output MSA
    # TODO: use refiner in place of mafft when available
    merged_subfam_msa = subfam_instances_path + merged_subfam_filename + ".afa"
    with open(merged_subfam_msa, "w") as f_out_merged_msa:
        process = Popen(
            ["mafft", merged_subfam_instances],
            stdout=f_out_merged_msa,
            stderr=PIPE,
        )
        process.communicate()
    # return msa filename to get the consensus from RepeatModeler
    return merged_subfam_filename, merged_subfam_name


def subfam_confidence(
    alignments: List[Alignment],
    lambs: List[float],
    subfam_instances_path: str,
    merge_stats_path: str,
    subfam_to_merged_num: Dict[str, int],
) -> Tuple[str, str, Tuple[str, str]]:
    """
    Finds and selects a subfamily pair to merge based
    on confidence values from subfamily alignments to
    a set of test sequences.

    input:
    alignments: list of subfam alignments to test seqs
    lambs: list of lambda values for each alignment (from Easel)
    subfam_instances_path: path to a dir of subfam instances
    merge_stats_path: path to a file to output stats from a subfam merge
    subfam_to_merged_num: dictionary that maps names of prior
    merged subfamilies to a number for instance file look up

    output:
    merged_consensus: consensus sequence of the merged subfam
    merged_name: name of the merged subfam
    sub_pair: names of the subfams that were selected for merging
    """
    subfams: List[str] = []
    scores: List[int] = []
    subfam_lambs: List[float] = []

    prev_test_seq_name: str = ""
    test_seqs: int = 0
    subfam_winners: Dict[str, int] = Counter()
    uncertain_subfam_pairs: Dict[Tuple[str, str], int] = Counter()
    winner_group_count: Dict[str, int] = Counter()

    for i, a in enumerate(alignments):
        cur_test_seq_name = (
            a.chrom_name + ":" + str(a.chrom_start) + "-" + str(a.chrom_stop)
        )
        # check if we are on a new test seq
        if cur_test_seq_name != prev_test_seq_name and len(subfams) != 0:
            test_seqs += 1
            # get confidence values for prev test seq
            test_seq_confidence(
                scores,
                lambs,
                subfams,
                subfam_winners,
                uncertain_subfam_pairs,
                winner_group_count,
            )
            # clear inputs for new test seq
            subfams = []
            scores = []
            subfam_lambs = []
        prev_test_seq_name = cur_test_seq_name
        subfams.append(a.subfamily)
        scores.append(a.score)
        subfam_lambs.append(lambs[i])

    # confidence values for last test seq
    test_seq_confidence(
        scores,
        lambs,
        subfams,
        subfam_winners,
        uncertain_subfam_pairs,
        winner_group_count,
    )

    # calc conf values of all uncertain subfam pairs
    subfam_pair_confidence, zero_conf_subfams = confidence_subfam_pairs(
        uncertain_subfam_pairs, subfam_winners
    )

    # sort subfams that were never found to be a clear winner
    # by the highest amount of uncertain pairs with another subfam
    sorted_zero = sorted(
        zero_conf_subfams.items(),
        key=lambda item: sorted(
            item[1].items(), key=lambda item: item[1], reverse=True
        )[0][1],
        reverse=True,
    )

    # return values will stay empty if no subfams should be merged
    merged_msa_file: str = ""
    merged_name: str = ""
    sub_pair: Tuple[str, str] = ("", "")

    # check to merge a subfam pair with zero confidence first
    for zero_conf_item in sorted_zero:
        zero_conf_highest_pair = sorted(
            zero_conf_item[1].items(), key=lambda item: item[1], reverse=True
        )[0]
        sub_pair = (zero_conf_item[0], zero_conf_highest_pair[0])
        merged_msa_file, merged_name = merge_subfams(
            zero_conf_item[0],
            zero_conf_highest_pair[0],
            subfam_instances_path,
            subfam_to_merged_num,
        )
        # output stats from merge
        if merge_stats_path:
            f_stats = open(merge_stats_path, "a")
            f_stats.write(str(sub_pair))
            f_stats.write("\n")
            f_stats.write(
                "winner group count: " + str(winner_group_count[sub_pair[0]])
            )
            f_stats.write("\n")
            f_stats.write(
                "uncertain pair count: " + str(zero_conf_highest_pair[1])
            )
            f_stats.write("\n")
            f_stats.close()
        return merged_msa_file, merged_name, sub_pair

    # sort uncertain subfamily pairs by confidence values
    # sorted_pairs = [((subfam_i,subfam_j),conf),((subfam_i, subfam_j),conf),...]
    sorted_pairs = sorted(
        subfam_pair_confidence.items(), key=lambda item: item[1]
    )

    # check to merge subfam pair with the lowest confidence
    if len(sorted_pairs) != 0 and sorted_pairs[0][1] < MERGE_CONF_THRESH:
        sub_pair = sorted_pairs[0][0]
        merged_msa_file, merged_name = merge_subfams(
            sub_pair[0],
            sub_pair[1],
            subfam_instances_path,
            subfam_to_merged_num,
        )
        clear_winner_counts = ""
        for sub in sub_pair:
            if sub in subfam_winners:
                clear_winner_counts += (
                    str(sub) + ": " + str(subfam_winners[sub]) + "\n"
                )
        # output stats from merge
        if merge_stats_path:
            f_stats = open(merge_stats_path, "a")
            f_stats.write(str(sub_pair) + " " + str(sorted_pairs[0][1]))
            f_stats.write("\n")
            f_stats.write(
                "winner group count: " + str(winner_group_count[sub_pair[0]])
            )
            f_stats.write("\n")
            f_stats.write(
                "uncertain pair count: "
                + str(uncertain_subfam_pairs[tuple(sorted(sub_pair))])
            )
            f_stats.write("\n")
            f_stats.write("clear winner counts")
            f_stats.write("\n")
            f_stats.write(clear_winner_counts)
            f_stats.write("\n")
            f_stats.close()

    return merged_msa_file, merged_name, sub_pair

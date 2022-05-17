from subprocess import Popen, PIPE
from collections import Counter
from typing import Dict, List, Tuple, Any, Iterable, TextIO
from .confidence_cm import confidence_only
from .alignment import Alignment
from .remove_cg_scores import remove_cg_scores
from .substitution_matrix import SubMatrixCollection


def test_seq_confidence(
    scores: List[int],
    lambs: List[float],
    subfams: List[str],
    subfam_winners: Dict[str, int],
    uncertain_subfam_pairs: Dict[Tuple[str, str], int],
    winner_group_count: Dict[str, int],
    winner_group_dist: Dict[str, Dict[int, int]],
    winner_group_thresh: float,
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

    # subfams with conf values >= to max conf * winner_group_thresh
    # are potential winners of this test seq
    potential_winner_thresh = confidence_list[-1] * winner_group_thresh

    # search for subfams above the potential winner threshold
    subfam_winner_group_count = 0
    winner_subs = set()
    for i in range(len(subfams) - 1, 0, -1):
        if confidence_list[i] < potential_winner_thresh:
            break
        # track count of potential subfam winners
        subfam_winner_group_count += 1
        winner_group_count[subfams[i]] += 1
        winner_subs.add(subfams[i])
        # create pairs of subfams that are both potential winners for the test seq
        # these are considered to be uncertain pairs
        for j in range(i - 1, 0, -1):
            if confidence_list[j] < potential_winner_thresh:
                break
            winner_subs.add(subfams[j])
            # single test seq could have 1+ alignments with the same subfam
            if subfams[i] != subfams[j]:
                subfam_pair = [subfams[i], subfams[j]]
                subfam_pair.sort()
                uncertain_subfam_pairs[tuple(subfam_pair)] += 1

    if subfam_winner_group_count == 1:
        # single subfam found above winner thresh is a clear
        # winner of this test seq
        subfam_winners[subfams[-1]] += 1
    # create dist
    # clear winner will have dist of 1
    for sub in winner_subs:
        if sub not in winner_group_dist:
            winner_group_dist[sub] = Counter()
        winner_group_dist[sub][subfam_winner_group_count] += 1


def calc_subfam_pair_independence(
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
    subfam_pair_independence: Dict[Tuple[str, str], float] = {}
    zero_conf_subfams: Dict[str, Any] = {}

    # compute subfam ij- and ji-pair confidence
    for subfam_pair, sub_pair_count in uncertain_subfam_pair_counts.items():
        subfam_i = subfam_pair[0]
        subfam_j = subfam_pair[1]
        # subfam ij-pair
        if subfam_i in subfam_winner_counts.keys():
            # compute subfam pair confidence
            subfam_pair_independence[
                (subfam_i, subfam_j)
            ] = subfam_winner_counts[subfam_i] / (
                sub_pair_count + subfam_winner_counts[subfam_i]
            )
        else:
            # subfam i was never a clear winner
            subfam_pair_independence[(subfam_i, subfam_j)] = 0
            if subfam_i not in zero_conf_subfams.keys():
                zero_conf_subfams[subfam_i] = {}
            zero_conf_subfams[subfam_i][subfam_j] = sub_pair_count
        # subfam ji-pair
        if subfam_j in subfam_winner_counts.keys():
            # compute subfam pair confidence
            subfam_pair_independence[
                (subfam_j, subfam_i)
            ] = subfam_winner_counts[subfam_j] / (
                sub_pair_count + subfam_winner_counts[subfam_j]
            )
        else:
            # subfam j was never a clear winner
            subfam_pair_independence[(subfam_j, subfam_i)] = 0
            if subfam_j not in zero_conf_subfams.keys():
                zero_conf_subfams[subfam_j] = {}
            zero_conf_subfams[subfam_j][subfam_i] = sub_pair_count

    return subfam_pair_independence, zero_conf_subfams


def merge_subfams(
    subfam_i: str,
    subfam_j: str,
    subfam_instances_path: str,
    total_merged_subfams: int,
) -> str:
    """
    Merges subfamilies i and j and creates a new
    instance file and consensus sequence for the
    merged subfamily.

    input:
    subfam_i: name of subfam i
    subfam_j: name of subfam j
    subfam_instances_path: path to a dir of subfam instances
    total_subfams: int

    output:
    merged_consensus_seq: consensus sequence of the merged subfam
    merged_subfam_name: name of the merged subfam
    """
    merged_subfam_name = subfam_i + "," + subfam_j

    # paths to subfam i and j instance files
    subfam_i_instances = subfam_instances_path + subfam_i + ".fa"
    subfam_j_instances = subfam_instances_path + subfam_j + ".fa"

    # path to merged instance file
    merged_subfam_name = "merged_" + str(total_merged_subfams)
    merged_subfam_instances = subfam_instances_path + merged_subfam_name + ".fa"

    # output fasta file with all subfam i and j instances
    outfile = open(merged_subfam_instances, "w")
    with open(subfam_i_instances, "r") as infile:
        outfile.write(infile.read())
    with open(subfam_j_instances, "r") as infile:
        outfile.write(infile.read())
    outfile.close()

    # output MSA
    # TODO: use refiner in place of mafft when available
    merged_subfam_msa = subfam_instances_path + merged_subfam_name + ".afa"
    with open(merged_subfam_msa, "w") as f_out_merged_msa:
        process = Popen(
            ["mafft", merged_subfam_instances],
            stdout=f_out_merged_msa,
            stderr=PIPE,
        )
        process.communicate()
    return merged_subfam_name


def subfam_confidence(
    alignments: Iterable[Alignment],
    lambs: List[float],
    subfam_instances_path: str,
    merge_stats_path: str,
    total_merged_subfams: int,
    winner_group_thresh: float,
    merge_thresh: float,
    ignore_cg_content: bool,
    sub_matrix_scores: SubMatrixCollection,
    all_merged_file: TextIO,
    cur_merged_file: TextIO,
) -> None:
    """
    Finds and selects a subfamily pair to merge based
    on confidence values from subfamily alignments to
    a set of test sequences.

    input:
    alignments: list of subfam alignments to test seqs
    lambs: list of lambda values for each alignment (from Easel)
    subfam_instances_path: path to a dir of subfam instances
    merge_stats_path: path to a file to output stats from a subfam merge
    total_merged_subfams: number of subfams that have been merged
    merged subfamilies to a number for instance file look up
    """
    subfams: List[str] = []
    scores: List[int] = []
    subfam_lambs: List[float] = []

    prev_test_seq_name: str = ""
    test_seqs: int = 0
    subfam_winners: Dict[str, int] = Counter()
    uncertain_subfam_pairs: Dict[Tuple[str, str], int] = Counter()
    winner_group_count: Dict[str, int] = Counter()
    winner_group_dist: Dict[str, Dict[int, int]] = {}

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
                winner_group_dist,
                winner_group_thresh,
            )
            # clear inputs for new test seq
            subfams = []
            scores = []
            subfam_lambs = []
        prev_test_seq_name = cur_test_seq_name
        subfams.append(a.subfamily)
        score = a.score
        if ignore_cg_content:
            score = remove_cg_scores(
                a, sub_matrix_scores[a.sub_matrix_name].scores
            )
        scores.append(score)
        subfam_lambs.append(lambs[i])

    # confidence values for last test seq
    test_seq_confidence(
        scores,
        lambs,
        subfams,
        subfam_winners,
        uncertain_subfam_pairs,
        winner_group_count,
        winner_group_dist,
        winner_group_thresh,
    )

    # calc conf values of all uncertain subfam pairs
    subfam_pair_independence, zero_conf_subfams = calc_subfam_pair_independence(
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
    merged_num: str = ""
    subfam_pair: Tuple[str, str] = ("", "")
    merged_subs = set()

    # check to merge a subfam pair with zero confidence first
    for zero_conf_item in sorted_zero:
        zero_conf_highest_pair = sorted(
            zero_conf_item[1].items(), key=lambda item: item[1], reverse=True
        )[0]
        subfam_pair = (zero_conf_item[0], zero_conf_highest_pair[0])
        subfam_i, subfam_j = subfam_pair[0], subfam_pair[1]
        # FIXME: check subfams haven't been merged in this iteration before merging
        if subfam_i not in merged_subs and subfam_j not in merged_subs:
            merged_subs.add(subfam_i)
            merged_subs.add(subfam_j)
            total_merged_subfams += 1
            merged_num = merge_subfams(
                subfam_i,
                subfam_j,
                subfam_instances_path,
                total_merged_subfams,
            )
            # write to all merged subfams file
            all_merged_file.write(
                # FIXME: used to output final set, could probably clean this up
                str(subfam_i)
                + " "
                + str(subfam_j)
                + " "
                + str(subfam_i)
                + ","
                + str(subfam_j)
            )
            all_merged_file.write("\n")

            # write to cur merged subfams file
            cur_merged_file.write(subfam_i + " " + subfam_j + " " + merged_num)
            cur_merged_file.write("\n")

            # output stats from merge
            if merge_stats_path:
                for sub in subfam_pair:
                    if sub not in subfam_winners:
                        subfam_winners[sub] = 0
                    if sub not in winner_group_count:
                        winner_group_count[sub] = 0
                        winner_group_dist[sub] = {}

                f_stats = open(merge_stats_path, "a")
                f_stats.write(
                    str(subfam_pair)
                    + " "
                    + str(subfam_pair_independence[subfam_pair])
                )
                f_stats.write("\n")
                reverse_pair = (subfam_j, subfam_i)
                f_stats.write(
                    str(reverse_pair)
                    + " "
                    + str(subfam_pair_independence[reverse_pair])
                )
                f_stats.write("\n")
                f_stats.write(
                    "winner group count subfam i: "
                    + str(winner_group_count[subfam_i])
                )
                f_stats.write("\n")
                f_stats.write(
                    "winner group count subfam j: "
                    + str(winner_group_count[subfam_j])
                )
                f_stats.write("\n")
                f_stats.write(
                    "winner group dist subfam i: "
                    + str(sorted(winner_group_dist[subfam_i].items()))
                )
                f_stats.write("\n")
                f_stats.write(
                    "winner group dist subfam j: "
                    + str(sorted(winner_group_dist[subfam_j].items()))
                )
                f_stats.write("\n")
                # same both directions, u_ij
                f_stats.write(
                    "uncertain pair count: "
                    + str(uncertain_subfam_pairs[tuple(sorted(subfam_pair))])
                )
                f_stats.write("\n")
                # w_i, w_j
                f_stats.write(
                    "clear winner count subfam i: "
                    + str(subfam_winners[subfam_i])
                )
                f_stats.write("\n")
                f_stats.write(
                    "clear winner count subfam j: "
                    + str(subfam_winners[subfam_j])
                )
                f_stats.write("\n")
                f_stats.close()

    # sort uncertain subfamily pairs by confidence values
    # sorted_pairs = [((subfam_i,subfam_j),conf),((subfam_i, subfam_j),conf),...]
    sorted_pairs = sorted(
        subfam_pair_independence.items(), key=lambda item: item[1]
    )

    # merge subfam pairs with conf < thresh
    for pair in sorted_pairs:
        conf = pair[1]
        if conf < merge_thresh:
            # could be more pairs under thresh that have not been merged
            subfam_pair = pair[0]
            subfam_i, subfam_j = subfam_pair[0], subfam_pair[1]
            if subfam_i not in merged_subs and subfam_j not in merged_subs:
                merged_subs.add(subfam_i)
                merged_subs.add(subfam_j)
                total_merged_subfams += 1
                merged_num = merge_subfams(
                    subfam_i,
                    subfam_j,
                    subfam_instances_path,
                    total_merged_subfams,
                )
                # write to all merged subfams file
                all_merged_file.write(
                    # FIXME: used to output final set, could probably clean this up
                    str(subfam_i)
                    + " "
                    + str(subfam_j)
                    + " "
                    + str(subfam_i)
                    + ","
                    + str(subfam_j)
                )
                all_merged_file.write("\n")

                # write to cur merged subfams file
                cur_merged_file.write(
                    str(subfam_i) + " " + str(subfam_j) + " " + str(merged_num)
                )
                cur_merged_file.write("\n")

                # output stats from merge
                if merge_stats_path:
                    for sub in subfam_pair:
                        if sub not in subfam_winners:
                            subfam_winners[sub] = 0
                        if sub not in winner_group_count:
                            winner_group_count[sub] = 0
                            winner_group_dist[sub] = {}

                    f_stats = open(merge_stats_path, "a")
                    f_stats.write(
                        str(subfam_pair)
                        + " "
                        + str(subfam_pair_independence[subfam_pair])
                    )
                    f_stats.write("\n")
                    reverse_pair = (subfam_j, subfam_i)
                    f_stats.write(
                        str(reverse_pair)
                        + " "
                        + str(subfam_pair_independence[reverse_pair])
                    )
                    f_stats.write("\n")
                    f_stats.write(
                        "winner group count subfam i: "
                        + str(winner_group_count[subfam_i])
                    )
                    f_stats.write("\n")
                    f_stats.write(
                        "winner group count subfam j: "
                        + str(winner_group_count[subfam_j])
                    )
                    f_stats.write("\n")
                    f_stats.write(
                        "winner group dist subfam i: "
                        + str(sorted(winner_group_dist[subfam_i].items()))
                    )
                    f_stats.write("\n")
                    f_stats.write(
                        "winner group dist subfam j: "
                        + str(sorted(winner_group_dist[subfam_j].items()))
                    )
                    f_stats.write("\n")
                    # same both directions, u_ij
                    f_stats.write(
                        "uncertain pair count: "
                        + str(
                            uncertain_subfam_pairs[tuple(sorted(subfam_pair))]
                        )
                    )
                    f_stats.write("\n")
                    # w_i, w_j
                    f_stats.write(
                        "clear winner count subfam i: "
                        + str(subfam_winners[subfam_i])
                    )
                    f_stats.write("\n")
                    f_stats.write(
                        "clear winner count subfam j: "
                        + str(subfam_winners[subfam_j])
                    )
                    f_stats.write("\n")
                    f_stats.close()
        else:  # no more pairs < thresh
            break

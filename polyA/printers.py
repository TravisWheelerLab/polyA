import json
from math import inf
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple, Union
from uuid import uuid4

from polyA.matrices import SupportMatrix
from .ultra_provider import TandemRepeat


def print_matrix_hash(
    num_col: int,
    num_row: int,
    subfams: List[str],
    matrix: Dict[Tuple[int, int], Union[float, int, str]],
    file: TextIO = stdout,
) -> None:
    """
    Print values contained in a matrix. Tab separated.
    This function exists for debugging purposes only.
    """
    file.write("\t")
    j: int = 0
    while j < num_col:
        file.write(f"{j}\t")
        j += 1
    file.write("\n")

    i: int = 0
    while i < num_row:
        file.write(f"{subfams[i]}\t")
        j: int = 0
        while j < num_col:
            if (i, j) in matrix:
                file.write(f"{matrix[i, j]}")
            else:
                file.write(f"{-inf}")
            file.write("\t")
            j += 1
        file.write("\n")
        i += 1


def print_matrix_support(
    num_col: int,
    start_all: int,
    chrom_start: int,
    matrix: SupportMatrix,
    subfams_collapse: List[str],
    outfile: TextIO = stdout,
) -> None:
    """
    Prints the given support matrix to `file` in a format appropriate
    to produce a heatmap.
    """
    outfile.write("\t")

    start: int = chrom_start + start_all - 1
    j: int = 0
    while j < num_col:
        outfile.write(f"{start}\t")
        start += 1
        j += 1
    outfile.write("\n")

    for k in range(len(subfams_collapse)):
        outfile.write(f"{subfams_collapse[k]}\t")
        j: int = 0
        while j < num_col:
            if (k, j) in matrix:
                outfile.write(str(matrix[k, j]))
            else:
                outfile.write(f"{-inf}")
            outfile.write("\t")
            j += 1
        outfile.write("\n")


def print_results(
    changes_orig: List[str],
    changespos_orig: List[int],
    columns_orig: List[int],
    ids: List[str],
) -> None:
    """
    prints the final results
    prints start and stop in terms of position in matrix

    TODO: Refactor to accept a TextIO instance
    """
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != "skip":
            stdout.write(str(columns_orig[changespos_orig[i]]))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1]))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def print_results_sequence(
    edgestart: int,
    changes_orig: List[str],
    changespos_orig: List[int],
    columns_orig: List[int],
    ids: List[str],
) -> None:
    """
    prints final results
    prints start and stop in terms of input target sequence
    """
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != "skip":
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart - 1))
            stdout.write("\t")
            stdout.write(
                str(columns_orig[changespos_orig[i + 1] - 1] + edgestart - 1)
            )
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def print_results_chrom(
    edgestart: int,
    chrom_start: int,
    changes_orig: List[str],
    changespos_orig: List[int],
    columns_orig: List[int],
    ids: List[str],
) -> None:
    """
    prints final results
    prints start and stop in terms of chromosome/target position
    """

    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != "skip":
            stdout.write(
                str(
                    columns_orig[changespos_orig[i]]
                    + edgestart
                    + chrom_start
                    - 2  # -1 for padding at start of DP, -1 for overlap between edgestart and chrom_start
                )
            )
            stdout.write("\t")
            stdout.write(
                str(
                    columns_orig[changespos_orig[i + 1] - 1]
                    + edgestart
                    + chrom_start
                    - 2  # -1 for padding at start of DP, -1 for overlap between edgestart and chrom_start
                )
            )
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def print_results_soda(
    start_all: int,
    outfile: Optional[TextIO],
    outfile_json: Optional[TextIO],
    chrom: str,
    chrom_start: int,
    chrom_end: int,
    subfams: List[str],
    changes_orig: List[str],
    changes_position_orig: List[int],
    columns_orig: List[int],
    consensus_lengths: Dict[str, int],
    strand_matrix_collapse: Dict[Tuple[int, int], str],
    consensus_matrix_collapse: Dict[Tuple[int, int], int],
    subfams_collapse_index: Dict[str, int],
    node_confidence_orig: Dict[Tuple[str, int], float],
    ids: List[int],
    subfam_alignments: List[str],
    chrom_alignments: List[str],
    subfam_alignments_collapse: Dict[Tuple[str, int], Tuple[int, int]],
    matrix: SupportMatrix,
    subfams_collapse: List[str],
    num_col: int,
) -> None:
    """
    prints the results in the proper format to input into the SODA visualization tool

    format described here:
    https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=joinedRmsk&hgta_table=rmskJoinedCurrent&hgta_doSchema=describe+table+schema
    """
    id: int = 0

    length = len(changes_position_orig) - 1

    used: List[int] = [
        1
    ] * length  # wont print out the results of the same thing twice

    json_dict_id: Dict[str, Dict[str, Dict[str, float]]] = {}

    json_dict = {}
    json_dict["chr"] = chrom
    json_dict["annotations"] = []
    json_dict["heatmap"] = {}

    min_align_start: int = chrom_end
    max_align_end: int = 0
    i = 0
    while i < length:

        sub_id: int = 0
        json_dict_subid: Dict[str, Dict[str, float]] = {}

        if changes_orig[i] != "skip" and used[i]:
            subfam: str = changes_orig[i]
            strand: str = strand_matrix_collapse[
                subfams_collapse_index[subfam],
                columns_orig[changes_position_orig[i]],
            ]

            left_flank: int
            right_flank: int

            if subfam == "Tandem Repeat":
                subfam = "Tandem#Repeat/TR"
                left_flank = 0
                right_flank = 0
            else:
                if strand == "-":
                    left_flank = (
                        consensus_lengths[subfam]
                        - consensus_matrix_collapse[
                            subfams_collapse_index[subfam],
                            columns_orig[changes_position_orig[i]],
                        ]
                    )
                    right_flank = (
                        consensus_matrix_collapse[
                            subfams_collapse_index[subfam],
                            columns_orig[changes_position_orig[i + 1] - 1],
                        ]
                        - 1
                    )
                else:
                    left_flank = (
                        consensus_matrix_collapse[
                            subfams_collapse_index[subfam],
                            columns_orig[changes_position_orig[i]],
                        ]
                        - 1
                    )
                    right_flank = (
                        consensus_lengths[subfam]
                        - consensus_matrix_collapse[
                            subfams_collapse_index[subfam],
                            columns_orig[changes_position_orig[i + 1] - 1],
                        ]
                    )

            align_start: int = chrom_start + (
                columns_orig[changes_position_orig[i]] + start_all
            )
            feature_start: int = align_start - left_flank
            align_stop: int = chrom_start + (
                columns_orig[changes_position_orig[i + 1] - 1] + start_all
            )
            feature_stop: int = align_stop + right_flank

            block_start_matrix: int = columns_orig[changes_position_orig[i]]

            block_count: int = 3
            id += 1

            json_dict_subfam_i: Dict[str, float] = {}

            for subfam_i in range(1, len(subfams)):
                subfamm = subfams[subfam_i]
                if (
                    subfamm,
                    i,
                ) in node_confidence_orig and node_confidence_orig[
                    subfamm, i
                ] > 0.001:
                    json_dict_subfam_i[subfamm] = node_confidence_orig[
                        subfamm, i
                    ]

            json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(
                json_dict_subfam_i.items(), key=lambda x: x[1], reverse=True
            )
            sub_id += 1

            block_start: List[str] = []
            block_size: List[str] = []

            block_start.append("-1")
            block_start.append(str(left_flank + 1))
            block_start.append("-1")

            block_size.append(str(left_flank))
            block_size.append(
                str(
                    columns_orig[changes_position_orig[i + 1] - 1]
                    - columns_orig[changes_position_orig[i]]
                    + 1
                )
            )
            block_size.append(str(right_flank))

            j: int = i + 1
            while j < length:
                if changes_orig[j] != "skip" and subfam != "Tandem#Repeat/TR":

                    if (
                        ids[columns_orig[changes_position_orig[i]]]
                        == ids[columns_orig[changes_position_orig[j]]]
                    ):
                        if strand == "-":
                            right_flank = (
                                consensus_matrix_collapse[
                                    subfams_collapse_index[changes_orig[j]],
                                    columns_orig[
                                        changes_position_orig[j + 1] - 1
                                    ],
                                ]
                                - 1
                            )
                        else:
                            right_flank = (
                                consensus_lengths[subfam]
                                - consensus_matrix_collapse[
                                    subfams_collapse_index[changes_orig[j]],
                                    columns_orig[
                                        changes_position_orig[j + 1] - 1
                                    ],
                                ]
                            )

                        del block_size[-1]
                        block_size.append("0")

                        align_stop: int = chrom_start + (
                            columns_orig[changes_position_orig[j + 1] - 1]
                            + start_all
                        )
                        feature_stop: int = align_stop + right_flank

                        block_start.append(
                            str(
                                columns_orig[changes_position_orig[j]]
                                + 1
                                - block_start_matrix
                                + left_flank
                            )
                        )
                        block_start.append("-1")

                        block_size.append(
                            str(
                                columns_orig[changes_position_orig[j + 1] - 1]
                                - columns_orig[changes_position_orig[j]]
                                + 1
                            )
                        )

                        block_size.append(str(right_flank))

                        block_count += 2

                        used[j] = 0

                        json_dict_subfam_j: Dict[str, float] = {}

                        for subfam_i in range(1, len(subfams)):
                            subfamm = subfams[subfam_i]
                            if (
                                (subfamm, j) in node_confidence_orig
                                and node_confidence_orig[subfamm, j] > 0.001
                            ):
                                json_dict_subfam_j[
                                    subfamm
                                ] = node_confidence_orig[subfamm, j]

                        json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(
                            json_dict_subfam_j.items(),
                            key=lambda x: x[1],
                            reverse=True,
                        )
                        sub_id += 1

                j += 1

            json_dict_id[str(id)] = json_dict_subid

            ucscString = (
                "000 "
                + chrom
                + " "
                + str(feature_start)
                + " "
                + str(feature_stop)
                + " "
                + subfam
                + " 0 "
                + strand
                + " "
                + str(align_start)
                + " "
                + str(align_stop)
                + " 0 "
                + str(block_count)
                + " "
                + (",".join(block_size))
                + " "
                + (",".join(block_start))
                + " "
                + str(id)
            )

            json_annotation = {}
            json_annotation["id"] = id
            json_annotation["blockCount"] = block_count
            json_annotation["ucscString"] = ucscString
            json_annotation["chrStart"] = align_start
            json_annotation["chrEnd"] = align_stop
            block_alignments = []
            if align_start < min_align_start:
                min_align_start = align_start
            if align_stop > max_align_end:
                max_align_end = align_stop

            # get alignments for each block
            if subfam != "Tandem#Repeat/TR":
                # col in seq
                subfam_start_col = align_start - chrom_start - 1
                subfam_stop_col = align_stop - chrom_start - 1
                subfam_rows = [
                    subfam_alignments_collapse[subfam, col]
                    for col in range(subfam_start_col, subfam_stop_col)
                    if (subfam, col) in subfam_alignments_collapse
                ]
                align_changes = [subfam_rows[0]]
                # get changes
                align_length = len(subfam_rows)
                for align_num in range(1, align_length):
                    if (
                        subfam_rows[align_num][0]
                        != subfam_rows[align_num - 1][0]
                    ):
                        align_changes.append(subfam_rows[align_num - 1])
                        align_changes.append(subfam_rows[align_num])
                align_changes.append(subfam_rows[align_length - 1])
                for align_num in range(0, len(align_changes) - 1, 2):
                    block_subfam = align_changes[align_num][
                        0
                    ]  # collapsed subfam row
                    block_sub_alignment = {}
                    block_sub_alignment["chrSeq"] = chrom_alignments[
                        block_subfam
                    ]
                    block_sub_alignment["famSeq"] = subfam_alignments[
                        block_subfam
                    ]
                    block_sub_alignment["alignStart"] = align_changes[
                        align_num
                    ][1]
                    block_sub_alignment["alignEnd"] = align_changes[
                        align_num + 1
                    ][1]
                    # consensus positions skip ahead, ex: [167, 407], [167, 416], ...
                    block_alignments.append(block_sub_alignment)
            json_annotation["alignments"] = block_alignments
            json_dict["annotations"].append(json_annotation)
        used[i] = 0
        i += 1

    json_dict["chrStart"] = min_align_start
    json_dict["chrEnd"] = max_align_end
    # Get heatmap values
    heatmap_dict = {}
    for k in range(len(subfams_collapse)):
        heatmap_vals = []
        j: int = 0
        # values in list
        while j < num_col:
            if (k, j) in matrix:
                heatmap_vals.append(str(matrix[k, j]))
            else:
                heatmap_vals.append("-inf")
            j += 1
        heatmap_dict[subfams_collapse[k]] = heatmap_vals
    json_dict["heatmap"] = heatmap_dict
    # prints  outfile for SODA viz
    outfile.write(json.dumps(json_dict))

    # prints json file with confidence values for each annotation
    outfile_json.write(json.dumps(json_dict_id))

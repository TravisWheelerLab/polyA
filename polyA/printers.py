import json
from math import inf
from sys import stdout
from typing import Dict, List, TextIO, Tuple, Union

from polyA.matrices import SupportMatrix


def print_matrix_hash(
    num_col: int,
    num_row: int,
    subfams: List[str],
    matrix: Dict[Tuple[int, int], Union[float, int, str]],
    file: TextIO = stdout,
) -> None:
    """
    Print values contained in a non-collapsed matrix.
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
    file: TextIO = stdout,
) -> None:
    """
    Prints the given support matrix to `file` in a format appropriate
    to produce a heatmap.
    """
    file.write("\t")

    start: int = chrom_start + start_all
    j: int = 0
    while j < num_col:
        file.write(f"{start}\t")
        start += 1
        j += 1
    file.write("\n")

    for k in range(len(subfams_collapse)):
        file.write(f"{subfams_collapse[k]}\t")
        j: int = 0
        while j < num_col:
            if (k, j) in matrix:
                file.write(str(matrix[k, j]))
            else:
                file.write(f"{-inf}")
            file.write("\t")
            j += 1
        file.write("\n")


def print_results(
    changes_orig: List[str],
    changespos_orig: List[int],
    columns_orig: List[int],
    ids: List[int],
) -> None:
    """
    prints the final results
    prints start and stop in terms of position in matrix
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
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
    ids: List[int],
) -> None:
    """
    prints final results
    prints start and stop in terms of input chrom sequence
    """
    # stdout.write("start\tstop\tID\tname\n")
    # stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != "skip":
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart))
            stdout.write("\t")
            stdout.write(
                str(columns_orig[changespos_orig[i + 1] - 1] + edgestart)
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
    ids: List[int],
) -> None:
    """
    prints final results
    prints start and stop in terms of chromosome/target position
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != "skip":
            stdout.write(
                str(columns_orig[changespos_orig[i]] + edgestart + chrom_start)
            )
            stdout.write("\t")
            stdout.write(
                str(
                    columns_orig[changespos_orig[i + 1] - 1]
                    + edgestart
                    + chrom_start
                )
            )
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def print_results_soda(
    start_all: int,
    outfile: str,
    outfile_json: str,
    chrom: str,
    chrom_start: int,
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

    with open(outfile, "w") as out:

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
                    if node_confidence_orig[subfamm, i] > 0.001:
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
                    if changes_orig[j] != "skip" and subfam != "Tandem Repeat":

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
                                        subfams_collapse_index[subfam],
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
                                    columns_orig[
                                        changes_position_orig[j + 1] - 1
                                    ]
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
                                if node_confidence_orig[subfamm, j] > 0.001:
                                    json_dict_subfam_j[
                                        subfamm
                                    ] = node_confidence_orig[subfamm, j]

                            json_dict_subid[
                                str(id) + "-" + str(sub_id)
                            ] = sorted(
                                json_dict_subfam_j.items(),
                                key=lambda x: x[1],
                                reverse=True,
                            )
                            sub_id += 1

                    j += 1

                json_dict_id[str(id)] = json_dict_subid

                out.write(
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
                out.write("\n")

            used[i] = 0
            i += 1

    # prints json file with confidence values for each annotation
    with open(outfile_json, "w") as out_json:
        out_json.write(json.dumps(json_dict_id))

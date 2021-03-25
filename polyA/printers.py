import json
from math import inf
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple, Union
from uuid import uuid4

from polyA.matrices import SupportMatrix, SubfamAlignmentsMatrix
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
    consensus_starts: List[int],
    consensus_stops: List[int],
    alignments_matrix: SubfamAlignmentsMatrix,
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
    json_dict["heatmap"] = []

    # heatmap
    subfam_ids: Dict[str, str] = {}
    subfam_ids[subfams_collapse[0]] = str(0)
    heatmap_dict = {"name": subfams_collapse[0], "id": 0}
    heatmap_vals = []
    confidence = []
    j: int = 0
    cur_subfam_row: int = 0
    start_offset: int = chrom_start + start_all - 2
    align_start: int = start_offset
    # skip state values
    while j < num_col:
        heatmap_vals.append(str(round(matrix[0, j], 3)))
        j += 1
    confidence.append({"chromStart": align_start, "values": heatmap_vals})
    heatmap_dict["confidence"] = confidence
    heatmap_dict["alignments"] = []
    json_dict["heatmap"].append(heatmap_dict)

    for k in range(1, len(subfams_collapse)):
        if subfams_collapse[k] == "Tandem Repeat":
            subfam_ids["Tandem#Repeat/TR"] = str(k)
        else:
            subfam_ids[subfams_collapse[k]] = str(k)
        heatmap_dict = {"name": subfams_collapse[k], "id": k}
        heatmap_vals = []
        subfam_rows = []
        confidence = []
        alignments = []
        j = 0
        cur_col = 0
        prev_subfam_row = 0
        align_start = start_offset
        while j < num_col:
            if (k, j) in matrix:
                # cur subfam row will be unique
                cur_subfam_row = alignments_matrix[
                    subfams_collapse[k], j + start_all - 1
                ]
                if j == cur_col and cur_subfam_row == prev_subfam_row:
                    # continue to add values
                    heatmap_vals.append(str(round(matrix[k, j], 3)))
                    subfam_rows.append(cur_subfam_row)
                    cur_col += 1
                else:
                    if len(heatmap_vals) > 0:
                        block_sub_alignment = {}
                        if subfams_collapse[k] != "Tandem Repeat":
                            block_sub_alignment["ID"] = prev_subfam_row
                            block_sub_alignment["chrSeq"] = chrom_alignments[
                                prev_subfam_row
                            ]
                            block_sub_alignment["famSeq"] = subfam_alignments[
                                prev_subfam_row
                            ]
                            block_sub_alignment["relativeStart"] = abs(
                                consensus_starts[prev_subfam_row]
                                - subfam_rows[0][1]
                            )
                            block_sub_alignment["relativeEnd"] = abs(
                                consensus_stops[prev_subfam_row]
                                - subfam_rows[-1][1]
                            )
                            block_sub_alignment[
                                "alignStart"
                            ] = consensus_starts[prev_subfam_row]
                            block_sub_alignment["alignEnd"] = consensus_stops[
                                prev_subfam_row
                            ]
                        alignments.append(block_sub_alignment)
                        confidence.append(
                            {"chromStart": align_start, "values": heatmap_vals}
                        )
                    heatmap_vals = [str(round(matrix[k, j], 3))]
                    subfam_rows = [cur_subfam_row]
                    cur_col = j + 1
                    align_start = j + start_offset
                    prev_subfam_row = cur_subfam_row
            j += 1
        if len(heatmap_vals) > 0:
            block_sub_alignment = {}
            if subfams_collapse[k] != "Tandem Repeat":
                block_sub_alignment["ID"] = cur_subfam_row
                block_sub_alignment["chrSeq"] = chrom_alignments[cur_subfam_row]
                block_sub_alignment["famSeq"] = subfam_alignments[
                    cur_subfam_row
                ]
                block_sub_alignment["relativeStart"] = abs(
                    consensus_starts[cur_subfam_row] - subfam_rows[0][1]
                )
                block_sub_alignment["relativeEnd"] = abs(
                    consensus_stops[cur_subfam_row] - subfam_rows[-1][1]
                )
                block_sub_alignment["alignStart"] = consensus_starts[
                    cur_subfam_row
                ]
                block_sub_alignment["alignEnd"] = consensus_stops[
                    cur_subfam_row
                ]
            alignments.append(block_sub_alignment)
            confidence.append(
                {"chromStart": chrom_start, "values": heatmap_vals}
            )
        heatmap_dict["confidence"] = confidence
        heatmap_dict["alignments"] = alignments
        json_dict["heatmap"].append(heatmap_dict)

    # annotations
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
            json_annotation = {}
            json_annotation["bin"] = "0"
            json_annotation["chrom"] = chrom
            json_annotation["chromStart"] = str(feature_start)
            json_annotation["chromEnd"] = str(feature_stop)
            json_annotation["name"] = subfam
            json_annotation["score"] = "0"
            json_annotation["strand"] = strand
            json_annotation["alignStart"] = str(align_start)
            json_annotation["alignEnd"] = str(align_stop)
            json_annotation["reserved"] = "0"
            json_annotation["blockCount"] = str(block_count)
            json_annotation["blockSizes"] = block_size
            json_annotation["blockStarts"] = block_start
            json_annotation["id"] = subfam_ids[subfam]
            json_dict["annotations"].append(json_annotation)

            if align_start < min_align_start:
                min_align_start = align_start
            if align_stop > max_align_end:
                max_align_end = align_stop

        used[i] = 0
        i += 1

    json_dict["chrStart"] = min_align_start
    json_dict["chrEnd"] = max_align_end
    # prints  outfile for SODA viz
    outfile.write(json.dumps(json_dict))

    # prints json file with confidence values for each annotation
    outfile_json.write(json.dumps(json_dict_id))

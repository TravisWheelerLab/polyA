import json
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple, Any

from .matrices import SupportMatrix, SubfamAlignmentsMatrix


class Printer:
    __output_file: TextIO = stdout
    __print_id: bool
    __soda_conf_file: Optional[TextIO]
    __soda_viz_file: Optional[TextIO]
    __use_matrix_position: bool
    __use_sequence_position: bool

    def __init__(
        self,
        output_file: Optional[TextIO] = None,
        print_id: bool = False,
        soda_conf_file: Optional[TextIO] = None,
        soda_viz_file: Optional[TextIO] = None,
        use_matrix_position: bool = False,
        use_sequence_position: bool = False,
    ) -> None:
        if output_file is not None:
            self.__output_file = output_file

        self.__print_id = print_id
        self.__soda_conf_file = soda_conf_file
        self.__soda_viz_file = soda_viz_file
        self.__use_matrix_position = use_matrix_position
        self.__use_sequence_position = use_sequence_position

    @property
    def use_matrix_position(self) -> bool:
        return self.__use_matrix_position

    @property
    def use_sequence_position(self) -> bool:
        return self.__use_sequence_position

    @property
    def use_soda_output(self) -> bool:
        return (
            self.__soda_viz_file is not None
            and self.__soda_conf_file is not None
        )

    def set_soda_files(
        self, viz_file: Optional[TextIO], conf_file: Optional[TextIO]
    ):
        self.__soda_viz_file = viz_file
        self.__soda_conf_file = conf_file

    def print_results_header(self):
        self.__output_file.write("start\tstop\t")

        if self.__print_id:
            self.__output_file.write("id\t")

        self.__output_file.write("name\n")

    def print_results_simple(
        self,
        start: int,
        stop: int,
        node_id: str,
        node_name: str,
    ) -> None:
        """
        A version of the output printer that takes the exact values to output.
        """
        self.__output_file.write(f"{start}\t{stop}\t")

        if self.__print_id:
            self.__output_file.write(f"{node_id}\t")

        self.__output_file.write(f"{node_name}\n")

    def print_results_matrix(
        self,
        changes: List[str],
        tr_changes: Dict[int, str],
        position_changes: List[int],
        columns: List[int],
        ids: List[str],
    ) -> None:
        """
        Print final results. Start and stop are in terms of matrix position.
        """
        for index, change in enumerate(changes):
            if changes[index] == "skip":
                continue

            position_change = position_changes[index]
            position_change_next = position_changes[index + 1]

            start = columns[position_change]
            stop = columns[position_change_next - 1]

            self.__output_file.write(f"{start}\t")
            self.__output_file.write(f"{stop}\t")

            if self.__print_id:
                self.__output_file.write(f"{ids[start]}\t")

            if change == "Tandem Repeat":
                tr_change = tr_changes[index]
                self.__output_file.write(tr_change)
            else:
                self.__output_file.write(change)

            self.__output_file.write("\n")

    def print_results_sequence(
        self,
        start_all: int,
        changes: List[str],
        tr_changes: Dict[int, str],
        position_changes: List[int],
        columns: List[int],
        ids: List[str],
    ) -> None:
        """
        Print final results. Start and stop are in terms of the target sequence.
        """
        for index, change in enumerate(changes):
            if changes[index] == "skip":
                continue

            position_change = position_changes[index]
            position_change_next = position_changes[index + 1]

            start = columns[position_change] + start_all - 1
            stop = columns[position_change_next - 1] + start_all - 1

            self.__output_file.write(f"{start}\t")
            self.__output_file.write(f"{stop}\t")

            if self.__print_id:
                self.__output_file.write(f"{ids[columns[position_change]]}\t")

            if change == "Tandem Repeat":
                tr_change = tr_changes[index]
                self.__output_file.write(tr_change)
            else:
                self.__output_file.write(change)

            self.__output_file.write("\n")

    def print_results_chrom(
        self,
        start_all: int,
        chrom_start: int,
        changes: List[str],
        tr_changes: Dict[int, str],
        position_changes: List[int],
        columns: List[int],
        ids: List[str],
    ) -> None:
        """
        Print final results. Start and stop are in terms of the chromosome /
        target position.
        """

        for index, change in enumerate(changes):
            if change == "skip":
                continue

            position_change = position_changes[index]
            position_change_next = position_changes[index + 1]

            start = columns[position_change] + start_all + chrom_start - 2
            stop = (
                columns[position_change_next - 1] + start_all + chrom_start - 2
            )

            self.__output_file.write(f"{start}\t")
            self.__output_file.write(f"{stop}\t")

            if self.__print_id:
                self.__output_file.write(f"{ids[columns[position_change]]}\t")

            if change == "Tandem Repeat":
                tr_change = tr_changes[index]
                self.__output_file.write(tr_change)
            else:
                self.__output_file.write(change)

            self.__output_file.write("\n")

    def print_results_soda(
        self,
        start_all: int,
        chrom: str,
        chrom_start: int,
        chrom_end: int,
        subfams: List[str],
        changes_orig: List[str],
        tr_consensus_changes: Dict[int, str],
        changes_position_orig: List[int],
        columns_orig: List[int],
        consensus_lengths: Dict[str, int],
        strand_matrix_collapse: Dict[Tuple[int, int], str],
        consensus_matrix_collapse: Dict[Tuple[int, int], int],
        subfams_collapse_index: Dict[str, int],
        node_confidence_orig: Dict[Tuple[str, int], float],
        ids: List[str],
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
        Prints the results in the proper format to input into the SODA visualization
        tool.

        The output format is described here:
        https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=joinedRmsk&hgta_table=rmskJoinedCurrent&hgta_doSchema=describe+table+schema
        """
        if self.__soda_conf_file is None or self.__soda_viz_file is None:
            return

        current_id: int = 0
        length = len(changes_position_orig) - 1

        # wont print out the results of the same thing twice
        used: List[int] = [1] * length

        json_dict_id: Dict[str, Any] = {}

        json_dict: Dict[str, Any] = {
            "chr": chrom,
            "annotations": [],
            "heatmap": [],
        }

        # heatmap
        subfam_ids: Dict[str, int] = {}
        subfam_subids: Dict[str, int] = {}
        subfam_ids[subfams_collapse[0]] = 0
        subfam_subids[subfams_collapse[0]] = 1
        heatmap_dict = {"name": subfams_collapse[0], "id": 0}
        heatmap_vals = []
        confidence = []
        j: int = 0
        cur_subfam_row: int = 0
        align_start: int = chrom_start - 1
        # skip state values
        while j < num_col:
            heatmap_vals.append(round(matrix[0, j], 3))
            j += 1
        confidence.append({"chromStart": align_start, "values": heatmap_vals})
        heatmap_dict["confidence"] = confidence
        heatmap_dict["alignments"] = []
        json_dict["heatmap"].append(heatmap_dict)
        for k in range(1, len(subfams_collapse)):
            subfam_ids[subfams_collapse[k]] = k
            subfam_subids[subfams_collapse[k]] = 1
            heatmap_dict = {"name": subfams_collapse[k], "id": k}
            heatmap_vals = []
            subfam_rows = []
            confidence = []
            alignments = []
            cur_col = 1
            prev_col = 0
            prev_subfam_row = 0
            align_start = chrom_start
            while cur_col < num_col:
                if (k, cur_col) in matrix:
                    # cur subfam row will be unique
                    subfam_row_consensus = alignments_matrix[
                        subfams_collapse[k], cur_col + start_all - 1
                    ]
                    cur_subfam_row = subfam_row_consensus[0]
                    if (
                        cur_col == prev_col
                        and cur_subfam_row == prev_subfam_row
                    ):
                        # continue to add values
                        heatmap_vals.append(round(matrix[k, cur_col], 3))
                        subfam_rows.append(subfam_row_consensus[1])
                        prev_col += 1
                    else:
                        if len(heatmap_vals) > 0:
                            block_sub_alignment: Dict[str, Any] = {}
                            if subfams_collapse[k] != "Tandem Repeat":
                                block_sub_alignment["id"] = prev_subfam_row
                                block_sub_alignment[
                                    "chrSeq"
                                ] = chrom_alignments[prev_subfam_row]
                                block_sub_alignment[
                                    "famSeq"
                                ] = subfam_alignments[prev_subfam_row]
                                block_sub_alignment["relativeStart"] = abs(
                                    consensus_starts[prev_subfam_row]
                                    - subfam_rows[0]
                                )
                                block_sub_alignment["relativeEnd"] = abs(
                                    consensus_stops[prev_subfam_row]
                                    - subfam_rows[-1]
                                )
                                block_sub_alignment[
                                    "alignStart"
                                ] = consensus_starts[prev_subfam_row]
                                block_sub_alignment[
                                    "alignEnd"
                                ] = consensus_stops[prev_subfam_row]
                                alignments.append(block_sub_alignment)
                            confidence.append(
                                {
                                    "chromStart": align_start,
                                    "values": heatmap_vals,
                                }
                            )
                        heatmap_vals = [round(matrix[k, cur_col], 3)]
                        subfam_rows = [subfam_row_consensus[1]]
                        align_start = chrom_start + cur_col + start_all - 2
                        prev_col = cur_col + 1
                        prev_subfam_row = cur_subfam_row
                cur_col += 1
            if len(heatmap_vals) > 0:
                block_sub_alignment = {}
                if subfams_collapse[k] != "Tandem Repeat":
                    block_sub_alignment["id"] = cur_subfam_row
                    block_sub_alignment["chrSeq"] = chrom_alignments[
                        cur_subfam_row
                    ]
                    block_sub_alignment["famSeq"] = subfam_alignments[
                        cur_subfam_row
                    ]
                    block_sub_alignment["relativeStart"] = abs(
                        consensus_starts[cur_subfam_row] - subfam_rows[0]
                    )
                    block_sub_alignment["relativeEnd"] = abs(
                        consensus_stops[cur_subfam_row] - subfam_rows[-1]
                    )
                    block_sub_alignment["alignStart"] = consensus_starts[
                        cur_subfam_row
                    ]
                    block_sub_alignment["alignEnd"] = consensus_stops[
                        cur_subfam_row
                    ]
                    alignments.append(block_sub_alignment)
                confidence.append(
                    {"chromStart": align_start, "values": heatmap_vals}
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

            # Note: PyCharm isn't (as of 2021-06) smart enough to figure
            # out the types here, but they are correct and Mypy gets it right.
            json_dict_subid: Dict[str, List[Tuple[str, float]]] = {}

            if changes_orig[i] != "skip" and used[i]:
                orig_subfam: str = changes_orig[i]
                subfam: str = changes_orig[i]
                strand: str = strand_matrix_collapse[
                    subfams_collapse_index[subfam],
                    columns_orig[changes_position_orig[i]],
                ]

                left_flank: int
                right_flank: int

                if orig_subfam == "Tandem Repeat":
                    subfam = tr_consensus_changes[i]
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

                align_start = (
                    chrom_start
                    + (columns_orig[changes_position_orig[i]] + start_all)
                    - 2
                )
                feature_start: int = align_start - left_flank
                align_stop: int = (
                    chrom_start
                    + (
                        columns_orig[changes_position_orig[i + 1] - 1]
                        + start_all
                    )
                    - 2
                )
                feature_stop: int = align_stop + right_flank

                block_start_matrix: int = columns_orig[changes_position_orig[i]]

                block_count: int = 3
                current_id += 1

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

                j = i + 1
                while j < length:
                    if (
                        changes_orig[j] != "skip"
                        and orig_subfam != "Tandem Repeat"
                    ):

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

                            align_stop = chrom_start + (
                                columns_orig[changes_position_orig[j + 1] - 1]
                                + start_all
                            )
                            feature_stop = align_stop + right_flank

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
                                if (
                                    (subfamm, j) in node_confidence_orig
                                    and node_confidence_orig[subfamm, j] > 0.001
                                ):
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
                json_annotation: Dict[str, Any] = {
                    "bin": "0",
                    "chrom": chrom,
                    "chromStart": str(feature_start),
                    "chromEnd": str(feature_stop),
                    "name": subfam,
                    "score": "0",
                    "strand": strand,
                    "alignStart": str(align_start),
                    "alignEnd": str(align_stop),
                    "reserved": "0",
                    "blockCount": str(block_count),
                    "blockSizes": block_size,
                    "blockStarts": block_start,
                    "id": (
                        str(subfam_ids[orig_subfam])
                        + "-"
                        + str(subfam_subids[orig_subfam])
                    ),
                }
                subfam_subids[orig_subfam] += 1

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
        self.__soda_viz_file.write(json.dumps(json_dict))

        # prints json file with confidence values for each annotation
        self.__soda_conf_file.write(json.dumps(json_dict_id))

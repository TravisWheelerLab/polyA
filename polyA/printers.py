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
        id: int = 0

        length = len(changes_position_orig) - 1

        used: List[int] = [
            1
        ] * length  # wont print out the results of the same thing twice

        json_dict_id: Dict[str, Any] = {}
        json_dict: Dict[str, Any] = {
            "chr": chrom,
            "annotations": [],
            "heatmap": {},
        }

        min_align_start: int = chrom_end
        max_align_end: int = 0
        i = 0
        while i < length:

            sub_id: int = 0
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

                ucsc_string = (
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

                json_annotation: Dict[str, Any] = {
                    "id": id,
                    "blockCount": block_count,
                    "ucscString": ucsc_string,
                    "chrStart": align_start,
                    "chrEnd": align_stop,
                }
                block_alignments = []
                if align_start < min_align_start:
                    min_align_start = align_start
                if align_stop > max_align_end:
                    max_align_end = align_stop

                # get alignments for each block
                if orig_subfam != "Tandem Repeat":
                    # col in seq
                    subfam_start_col = align_start - chrom_start - 1
                    subfam_stop_col = align_stop - chrom_start - 1
                    subfam_rows = [
                        alignments_matrix[subfam, col]
                        for col in range(subfam_start_col, subfam_stop_col)
                        if (subfam, col) in alignments_matrix
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
                    block_sub_alignment: Dict[str, Any]
                    for align_num in range(0, len(align_changes) - 1, 2):
                        block_subfam = align_changes[align_num][
                            0
                        ]  # collapsed subfam row
                        block_sub_alignment = {
                            "chrSeq": chrom_alignments[block_subfam],
                            "famSeq": subfam_alignments[block_subfam],
                            "alignStart": align_changes[align_num][1],
                            "alignEnd": align_changes[align_num + 1][1],
                        }
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
            j = 0
            # values in list
            while j < num_col:
                if (k, j) in matrix:
                    heatmap_vals.append(str(matrix[k, j]))
                else:
                    heatmap_vals.append("-inf")
                j += 1
            heatmap_dict[subfams_collapse[k]] = heatmap_vals
        json_dict["heatmap"] = heatmap_dict

        if self.__soda_viz_file:
            self.__soda_viz_file.write(json.dumps(json_dict))

        if self.__soda_conf_file:
            self.__soda_conf_file.write(json.dumps(json_dict_id))

    def set_soda_files(
        self,
        viz_file: Optional[TextIO] = None,
        conf_file: Optional[TextIO] = None,
    ):
        """
        Set the output files to be used when producing SODA output. Setting a
        file to `None` will close the existing file first.
        """
        if self.__soda_conf_file and self.__soda_conf_file != conf_file:
            self.__soda_conf_file.close()
        self.__soda_conf_file = conf_file

        if self.__soda_viz_file and self.__soda_viz_file != viz_file:
            self.__soda_viz_file.close()
        self.__soda_viz_file = viz_file

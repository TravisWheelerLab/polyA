import json
import importlib.resources
import requests
from sys import stdout
from typing import Dict, List, Optional, TextIO, Tuple, Any

from .matrices import SupportMatrix, SubfamAlignmentsMatrix
from . import soda

polya_soda = importlib.resources.read_text(soda, 'polya-soda.js')

char_complement = {"A": "T", "T": "A", "G": "C", "C": "G"}


def complement(sequence: str) -> str:
    seq_complement: str = ""
    for seq_char in sequence:
        if seq_char in {"A", "T", "C", "G"}:
            seq_complement += char_complement[seq_char]
        else:
            seq_complement += seq_char
    return seq_complement


def calc_relative_start_and_end(
        confidence_start: int,
        alignment_start: int,
        confidence_length: int,
        target_seq: str,
):
    # starts relative to chrom position
    confidence_offset = confidence_start - alignment_start
    relative_start = -1
    relative_end = -1

    # find relative start position
    non_gap_count = 0
    for i, char in enumerate(target_seq):
        if non_gap_count == confidence_offset:
            relative_start = i
            break
        if char != "-":
            non_gap_count += 1

    # find relative end position
    non_gap_count = 0
    for i, char in enumerate(target_seq[relative_start:]):
        if char != "-":
            non_gap_count += 1
        if non_gap_count == confidence_length:
            relative_end = relative_start + i + 1
            break

    return relative_start, relative_end


class Printer:
    __output_file: TextIO = stdout
    __print_id: bool
    __soda_viz_file: Optional[TextIO]
    __use_matrix_position: bool
    __use_sequence_position: bool

    def __init__(
            self,
            output_file: Optional[TextIO] = None,
            print_id: bool = False,
            soda_viz_file: Optional[TextIO] = None,
            use_matrix_position: bool = False,
            use_sequence_position: bool = False,
    ) -> None:
        if output_file is not None:
            self.__output_file = output_file

        self.__print_id = print_id
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
        )

    def set_soda_files(
            self, viz_file: Optional[TextIO]
    ):
        self.__soda_viz_file = viz_file

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
            chrom_starts: List[int],
            chrom_stops: List[int],
            alignments_matrix: SubfamAlignmentsMatrix,
            matrix: SupportMatrix,
            subfams_collapse: List[str],
            num_col: int,
    ) -> None:
        """
        Produce a polya-soda html file with the results.

        """
        if self.__soda_viz_file is None:
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
        # skip state starts one before the first alignment
        conf_chrom_start: int = chrom_start - 2 + start_all
        # skip state values
        while j < num_col:
            heatmap_vals.append(round(matrix[0, j], 3))
            j += 1
        confidence.append({"start": conf_chrom_start, "values": heatmap_vals})
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
            next_expected_col = 0
            prev_subfam_row = 0
            conf_chrom_start = (
                chrom_start  # chrom position for start of conf values
            )
            while cur_col < num_col:
                if (k, cur_col) in matrix:
                    # cur subfam row will be unique
                    # alignments_matrix[subfam, collapse_col + start_all - 1] = (row_index, consensus_pos)
                    # col val between stop_all - start_all + 1 + 2, sequence position
                    # + start_all - 1 -> on chrom
                    # at chrom pos on target, which consensus pos does it correspond to
                    subfam_row_consensus = alignments_matrix[
                        subfams_collapse[k], cur_col + start_all - 1
                    ]
                    cur_subfam_row = subfam_row_consensus[0]
                    if (
                            cur_col == next_expected_col
                            and cur_subfam_row == prev_subfam_row
                    ):
                        # continue to add values
                        heatmap_vals.append(round(matrix[k, cur_col], 3))
                        subfam_rows.append(subfam_row_consensus[1])
                        next_expected_col += 1
                    else:
                        if len(heatmap_vals) > 0:
                            block_sub_alignment: Dict[str, Any] = {}
                            if subfams_collapse[k] != "Tandem Repeat":
                                block_sub_alignment["id"] = prev_subfam_row

                                (
                                    relative_start,
                                    relative_end,
                                ) = calc_relative_start_and_end(
                                    conf_chrom_start,
                                    chrom_starts[prev_subfam_row]
                                    + chrom_start
                                    - 1,
                                    len(heatmap_vals),
                                    chrom_alignments[prev_subfam_row],
                                )
                                block_sub_alignment[
                                    "relativeStart"
                                ] = relative_start
                                block_sub_alignment[
                                    "relativeEnd"
                                ] = relative_end
                                if subfam_rows[0] > subfam_rows[-1]:
                                    # complement and swap align start and stop positions
                                    block_sub_alignment["target"] = complement(
                                        chrom_alignments[prev_subfam_row]
                                    )
                                    block_sub_alignment["query"] = complement(
                                        subfam_alignments[prev_subfam_row]
                                    )
                                    block_sub_alignment[
                                        "alignStart"
                                    ] = consensus_stops[prev_subfam_row]
                                    block_sub_alignment[
                                        "alignEnd"
                                    ] = consensus_starts[prev_subfam_row]
                                else:
                                    block_sub_alignment[
                                        "target"
                                    ] = chrom_alignments[prev_subfam_row]
                                    block_sub_alignment[
                                        "query"
                                    ] = subfam_alignments[prev_subfam_row]
                                    block_sub_alignment[
                                        "alignStart"
                                    ] = consensus_starts[prev_subfam_row]
                                    block_sub_alignment[
                                        "alignEnd"
                                    ] = consensus_stops[prev_subfam_row]
                                block_sub_alignment["start"] = (
                                        chrom_starts[prev_subfam_row]
                                        + chrom_start
                                        - 1
                                )
                                block_sub_alignment["end"] = (
                                        chrom_stops[prev_subfam_row]
                                        + chrom_start
                                        - 1
                                )
                                alignments.append(block_sub_alignment)
                            confidence.append(
                                {
                                    "start": conf_chrom_start,
                                    "values": heatmap_vals,
                                }
                            )
                        heatmap_vals = [round(matrix[k, cur_col], 3)]
                        subfam_rows = [subfam_row_consensus[1]]
                        conf_chrom_start = chrom_start + cur_col + start_all - 2
                        next_expected_col = cur_col + 1
                        prev_subfam_row = cur_subfam_row
                cur_col += 1
            if len(heatmap_vals) > 0:
                block_sub_alignment = {}
                if subfams_collapse[k] != "Tandem Repeat":
                    block_sub_alignment["id"] = cur_subfam_row
                    (
                        relative_start,
                        relative_end,
                    ) = calc_relative_start_and_end(
                        conf_chrom_start,
                        chrom_starts[cur_subfam_row] + chrom_start - 1,
                        len(heatmap_vals),
                        chrom_alignments[cur_subfam_row],
                    )
                    block_sub_alignment["relativeStart"] = relative_start
                    block_sub_alignment["relativeEnd"] = relative_end
                    if subfam_rows[0] > subfam_rows[-1]:
                        # complement and swap align start and end
                        block_sub_alignment["target"] = complement(
                            chrom_alignments[cur_subfam_row]
                        )
                        block_sub_alignment["query"] = complement(
                            subfam_alignments[cur_subfam_row]
                        )
                        block_sub_alignment["alignStart"] = consensus_stops[
                            cur_subfam_row
                        ]
                        block_sub_alignment["alignEnd"] = consensus_starts[
                            cur_subfam_row
                        ]
                    else:
                        block_sub_alignment["target"] = chrom_alignments[
                            cur_subfam_row
                        ]
                        block_sub_alignment["query"] = subfam_alignments[
                            cur_subfam_row
                        ]
                        block_sub_alignment["alignStart"] = consensus_starts[
                            cur_subfam_row
                        ]
                        block_sub_alignment["alignEnd"] = consensus_stops[
                            cur_subfam_row
                        ]
                    block_sub_alignment["start"] = (
                            chrom_starts[cur_subfam_row] + chrom_start - 1
                    )
                    block_sub_alignment["end"] = (
                            chrom_stops[cur_subfam_row] + chrom_start - 1
                    )
                    alignments.append(block_sub_alignment)
                confidence.append(
                    {"start": conf_chrom_start, "values": heatmap_vals}
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

                block_start: List[int] = []
                block_size: List[int] = []

                block_start.append(-1)
                block_start.append(left_flank + 1)
                block_start.append(-1)

                block_size.append(left_flank)
                block_size.append(
                    columns_orig[changes_position_orig[i + 1] - 1]
                    - columns_orig[changes_position_orig[i]]
                    + 1
                )
                block_size.append(right_flank)

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
                            block_size.append(0)

                            align_stop = chrom_start + (
                                    columns_orig[changes_position_orig[j + 1] - 1]
                                    + start_all
                            )
                            feature_stop = align_stop + right_flank

                            block_start.append(
                                columns_orig[changes_position_orig[j]]
                                + 1
                                - block_start_matrix
                                + left_flank
                            )
                            block_start.append(-1)

                            block_size.append(
                                columns_orig[changes_position_orig[j + 1] - 1]
                                - columns_orig[changes_position_orig[j]]
                                + 1
                            )

                            block_size.append(right_flank)

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
                    "bin": 0,
                    "chrom": chrom,
                    "chromStart": feature_start,
                    "chromEnd": feature_stop,
                    "name": subfam,
                    "score": 0,
                    "strand": strand,
                    "alignStart": align_start,
                    "alignEnd": align_stop,
                    "reserved": 0,
                    "blockCount": block_count,
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

        json_dict["start"] = min_align_start
        json_dict["end"] = max_align_end

        genome_url = f"https://sodaviz.org/data/hg38/{chrom}/{min_align_start}/{max_align_end}"
        ucsc_url = f"https://sodaviz.org/data/rmsk/{chrom}/{min_align_start}/{max_align_end}"

        genome_seq_response = requests.get(genome_url)
        ucsc_annotations_response = requests.get(ucsc_url)

        if genome_seq_response.status_code == 200:
            genome_seq = genome_seq_response.text
        else:
            print(f"failed to perform genome get request: {genome_url}")
            print("genome will be missing from the visualization")
            genome_seq = ""

        if ucsc_annotations_response.status_code == 200:
            ucsc_annotations = ucsc_annotations_response.json()
        else:
            print(f"failed to perform ucsc get request: {ucsc_url}")
            print("ucsc annotations will be missing from the visualization")
            ucsc_annotations = []

        json_dict["genomeSeq"] = genome_seq
        json_dict["ucscAnnotations"] = ucsc_annotations

        data_json_string = json.dumps(json_dict)
        bundle = importlib.resources.read_text(soda, 'polya-soda.js')
        title = f"{chrom}: {min_align_start} - {max_align_end}"
        html_blob = f"""<!DOCTYPE html>
         <html>
         <head>
           <title>{title}</title>
         </head>
         <body>
         <h1>{title}</h1>
         <div id="charts"></div>
         <div class="controls">
           <button class="button" id="toggle-confidence">toggle confidence</button>
           <button class="button" id="toggle-alignments">toggle alignments</button>
         </div>
         <script>
           {bundle}
           let data = {data_json_string};
         </script>
         <script>
           let container = new ps.PolyaContainer({{selector: "#charts"}});
           document.getElementById("toggle-confidence").addEventListener("click", () => container.toggleConfidence());
           document.getElementById("toggle-alignments").addEventListener("click", () => container.toggleAlignments());
           container.render(data);
         </script>
         </body>
         </html>"""

        # prints  outfile for SODA viz
        self.__soda_viz_file.write(html_blob)

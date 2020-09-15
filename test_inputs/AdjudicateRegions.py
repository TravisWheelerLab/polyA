from getopt import getopt
from math import log
import re
from sys import argv, stdout, stderr
from typing import Dict, List, Tuple
import os
import json
import subprocess
import tempfile

from polyA.collapse_matrices import collapse_matrices
from polyA.extract_nodes import extract_nodes
from polyA.fill_align_matrix import fill_align_matrix
from polyA.calculate_repeat_scores import calculate_repeat_scores
from polyA.fill_confidence_matrix import fill_confidence_matrix, fill_confidence_matrix_tr
from polyA.fill_consensus_position_matrix import fill_consensus_position_matrix
from polyA.fill_node_confidence import fill_node_confidence
from polyA.fill_path_graph import fill_path_graph
from polyA.fill_probability_matrix import fill_probability_matrix
from polyA.fill_support_matrix import fill_support_matrix
from polyA.get_path import get_path
from polyA.lambda_provider import EaselLambdaProvider
from polyA.load_alignments import load_alignments


# -----------------------------------------------------------------------------------#
#			   MAIN														   			#
# -----------------------------------------------------------------------------------#

from polyA.pad_sequences import pad_sequences
from polyA.printers import (
    print_matrix_support,
    print_results,
    print_results_chrom,
    print_results_sequence,
    print_results_soda,
)

if __name__ == "__main__":
    GapInit: int = -25
    GapExt: int = -5
    Lamb: float = 0.0
    EslPath = ""
    ChunkSize: int = 31
    SameProbLog: float = 0.0
    ChangeProb: float = 10 ** -45
    ChangeProbLog: float = 0.0  # Reassigned later
    ChangeProbSkip: float = 0.0  # Reassigned later
    SameProbSkip: float = 0.0
    SkipAlignScore: float = 0.0

    StartAll: int = 0  # Reassigned later
    StopAll: int = 0  # Reassigned later
    ID: int = 1111

    infile_prior_counts: str = ""
    outfile_viz: str = ""
    outfile_heatmap: str = ""

    # running ultra and esl with raw data
    esl_sfetch_path: str = ""
    ultra_path: str = ""
    seq_file: str = ""
    start_pos: int = 0
    end_pos: int = 0
    chrom_name: str = ""
    # using ultra output
    ultra_output_path: str = ""  # json file

    help: bool = False  # Reassigned later
    prin: bool = False  # Reassigned later
    printMatrixPos: bool = False  # Reassigned later
    printSeqPos: bool = False

    helpMessage: str = f"""
    usage: {argv[0]} alignFile subMatrixFile\n
    ARGUMENTS
        --GapInit[-25]
        --getExt[-5]
        --lambda [will calc from matrix if not included]
        --eslPath [specify path to easel]
        --segmentsize (must be odd) [31]
        --changeprob[1e-45]
        --priorCounts PriorCountsFile
        --ultraPath [specify path to ULTRA executable]
        --seqFile [specify path to genome file]
        --startPos [start of sequence region]
        --endPos [end of sequence region]
        --chromName [name of chromosome] (must match genome file name)
        --ultraOutput [specify path to ULTRA output file]
    
    OPTIONS
        --help - display help message
        --matrixpos - prints subfam changes in matrix position instead of genomic position
        --sequencepos - prints subfam changes in sequence position instead of genomic position
        --viz outfile - prints output format for SODA visualization
        --heatmap outfile - prints probability file for input into heatmap
    """
    raw_opts, args = getopt(argv[1:], "", [
        "GapInit=",
        "GapExt=",
        "skipScore=",
        "lambda=",
        "eslPath=",
        "segmentsize=",
        "changeprob=",
        "priorCounts=",
        "viz=",
        "heatmap=",
        "ultraPath=",
        "eslSfetch=",
        "seqFile=",
        "startPos=",
        "endPos=",
        "chromName=",
        "ultraOutput=",

        "help",
        "matrixpos",
        "seqpos",
    ])
    opts = dict(raw_opts)

    GapInit = int(opts["--GapInit"]) if "--GapInit" in opts else GapInit
    GapExt = int(opts["--GapExt"]) if "--GapExt" in opts else GapExt
    Lamb = float(opts["--lambda"]) if "--lambda" in opts else Lamb
    EslPath = str(opts["--eslPath"]) if "--eslPath" in opts else EslPath
    ChunkSize = (
        int(opts["--segmentsize"]) if "--segmentsize" in opts else ChunkSize
    )
    ChangeProb = (
        float(opts["--changeprob"]) if "--changeprob" in opts else ChangeProb
    )
    infile_prior_counts = (
        str(opts["--priorCounts"])
        if "--priorCounts" in opts
        else infile_prior_counts
    )
    outfile_viz = str(opts["--viz"]) if "--viz" in opts else outfile_viz
    outfile_heatmap = str(opts["--heatmap"]) if "--heatmap" in opts else outfile_heatmap

    # options for using ultra
    esl_sfetch_path = str(opts["--eslSfetch"]) if "--eslSfetch" in opts else esl_sfetch_path
    ultra_path = str(opts["--ultraPath"]) if "--ultraPath" in opts else ultra_path # gets path to exe
    seq_file = str(opts["--seqFile"]) if "--seqFile" in opts else seq_file
    start_pos = str(opts["--startPos"]) if "--startPos" in opts else start_pos
    end_pos = str(opts["--endPos"]) if "--endPos" in opts else end_pos
    chrom_name = str(opts["--chromName"]) if "--chromName" in opts else chrom_name
    ultra_output_path = str(opts["--ultraOutput"]) if "--ultraOutput" in opts else ultra_output_path

    help = "--help" in opts
    printMatrixPos = "--matrixpos" in opts
    printSeqPos = "--seqpos" in opts

    outfile_conf = outfile_viz + ".json"

    if help:
        print(helpMessage)
        exit(0)

    # input is alignment file of hits region and substitution matrix
    infile: str = args[0]
    infile_matrix: str = args[1]

    # Other open was moved down to where we load the alignments file
    with open(infile_matrix) as _infile_matrix:
        in_matrix: List[str] = _infile_matrix.readlines()

    # if command line option included to use prior counts into in confidence calculations
    if infile_prior_counts:
        with open(infile_prior_counts) as _infile_prior_counts:
            in_counts: List[str] = _infile_prior_counts.readlines()

    # if lambda isn't included at command line, run esl_scorematrix to calculate it from scorematrix
    if not Lamb:
        provider = EaselLambdaProvider(EslPath)
        Lamb = provider()

    # reads in the score matrix from file and stores in dict that maps 'char1char2' to the score from the
    # input substitution matrix - ex: 'AA' = 8
    SubMatrix: Dict[str, int] = {}
    line = in_matrix[0]
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    chars = re.split(r"\s+", line)

    count: int = 0
    for line in in_matrix[1:]:
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r"\s+$", "", line)
        subScores = re.split(r"\s+", line)
        for i in range(len(subScores)):
            SubMatrix[chars[count] + chars[i]] = int(subScores[i])
        count += 1
    SubMatrix[".."] = 0

    TR: bool = False
    # TODO: Make running ultra and such an integrated field
    # running esl and ultra
    # should already have a correct fasta file from getting alignments
    if ultra_path and esl_sfetch_path and seq_file and start_pos and end_pos and chrom_name:
        TR = True
        # pre-processing
        subprocess.call([esl_sfetch_path, '--index', seq_file])
        # get smaller region
        file, small_region = tempfile.mkstemp()
        subprocess.call([esl_sfetch_path, '-o', small_region, '-c',
                         start_pos + '..' + end_pos, seq_file, chrom_name])
        # run ULTRA
        pop = subprocess.Popen([ultra_path, '-ss', small_region], stdout=subprocess.PIPE,
                               universal_newlines=True)
        ultra_out, err = pop.communicate()
        ultra_output = json.loads(ultra_out)
        os.remove(small_region)
    elif ultra_path and seq_file:
        TR = True
        # run ULTRA
        pop = subprocess.Popen([ultra_path, '-ss', seq_file], stdout=subprocess.PIPE,
                               universal_newlines=True)
        ultra_out, err = pop.communicate()
        ultra_output = json.loads(ultra_out)
    elif ultra_output_path:
        TR = True
        # my path: /Users/audrey/tr.json
        with open(ultra_output_path) as f:
            ultra_output = json.load(f)

    # Get Tandem Repeats from ULTRA output
    TR_count: int = 0
    if TR:
        TandemRepeats = ultra_output['Repeats']
        if len(TandemRepeats) == 0:
            TR = False
        else:
            TR_count = len(TandemRepeats)

    # maps subfam names to genomic prior_count/total_in_genome from input file
    # used during confidence calculations
    SubfamCounts: Dict[str, float] = {}
    PriorTotal: float = 0
    prob_skip = 0.4  # about 60% of genome is TE derived
    prob_tr = 0.06  # ~6% of genome exepected to be tandem repeat
    if infile_prior_counts:
        SubfamCounts["skip"] = prob_skip
        if TR:
            SubfamCounts["Tandem Repeat"] = prob_tr
            prob_skip += prob_tr

        for line in in_counts[1:]:
            line = re.sub(r"\n", "", line)
            info = re.split(r"\s+", line)
            # FIXME: infile counts do not match subfam names
            count = info[1]
            subfam = info[0]
            SubfamCounts[subfam] = int(count)
            PriorTotal += float(info[1])

        for key in SubfamCounts:
            SubfamCounts[key] = (1 - prob_skip) * SubfamCounts[key] / PriorTotal

    Subfams: List[str] = []
    Chroms: List[str] = []
    Scores: List[int] = []
    Strands: List[str] = []
    Starts: List[int] = []
    Stops: List[int] = []
    ConsensusStarts: List[int] = []
    ConsensusStops: List[int] = []
    SubfamSeqs: List[str] = []
    ChromSeqs: List[str] = []
    Flanks: List[int] = []

    RepeatScores: Dict[int, float] = {}
    AlignMatrix: Dict[Tuple[int, int], float] = {}
    SingleAlignMatrix: Dict[Tuple[int, int], int] = {}
    ConfidenceMatrix: Dict[Tuple[int, int], float] = {}
    SupportMatrix: Dict[Tuple[int, int], float] = {}
    OriginMatrix: Dict[Tuple[int, int], int] = {}
    ConsensusMatrix: Dict[Tuple[int, int], int] = {}
    SameSubfamChangeMatrix: Dict[Tuple[int, int], int] = {}
    ProbMatrixLastColumn: List[float] = []

    NonEmptyColumns: List[int] = []
    ActiveCells: Dict[int, List[int]] = {}

    Changes: List[str] = []
    ChangesPosition: List[int] = []

    IDs: List[int] = []
    ChangesOrig: List[str] = []
    ChangesPositionOrig: List[int] = []
    NonEmptyColumnsOrig: List[int] = []

    # for graph/node part
    NumNodes: int = 0
    NodeConfidence: Dict[Tuple[str, int], float] = {}
    NodeConfidenceOrig: Dict[Tuple[str, int], float] = {}
    PathGraph: List[int] = []
    total: int = 0
    loop: int = 1

    # opens alignment file and stores all Subfams, Scores, Starts, Stops, subfams seqs and Chrom seqs in arrays
    numseqs: int = 0
    with open(infile) as _infile:
        alignments = load_alignments(_infile)
        for alignment in alignments:
            numseqs += 1
            # FIXME: split subfam to match prior counts file for testing - can't always split by #
            Subfams.append(alignment.subfamily)
            Chroms.append(alignment.chrom)
            Scores.append(alignment.score)
            Strands.append(alignment.strand)
            Starts.append(alignment.start)
            Stops.append(alignment.stop)
            ConsensusStarts.append(alignment.consensus_start)
            ConsensusStops.append(alignment.consensus_stop)
            SubfamSeqs.append(alignment.subfamily_sequence)
            ChromSeqs.append(alignment.sequence)
            Flanks.append(alignment.flank)

    # if there is only one subfam in the alignment file, no need to run anything because we know
    # that subfam is the annotation
    # numseqs = 2 because of the skip state
    if numseqs == 2:
        if printMatrixPos:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{0}\t{Stops[1] - Starts[1]}\t1111\t{Subfams[1]}\n")
        else:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{Starts[1]}\t{Stops[1]}\t1111\t{Subfams[1]}\n")
        exit()

    match = re.search(r"(.+)/(\d+)-(\d+)", Chroms[1])  # FIXME: Replace : with / to work with the TR test align files
    Chrom: str = match.groups()[0]
    ChromStart: int = int(match.groups()[1])
    ChromEnd: int = int(match.groups()[2])
    TargetLen: int = ChromEnd - ChromStart

    # bail out if target sq is < 25 nucls
    # warning if less than 1000 nucls
    # warning if no chrom info given - okay for artificial seq inputs
    if TargetLen == 0:
        stderr.write(
            "WARNING - No chromosome position information given.\nThis is okay if running on artificial sequences, but cannot use command line options --viz or --heatmap.\n\n"
        )
    elif TargetLen <= 25:
        stderr.write(
            "ERROR - Target sequence length needs to be > 25 nucleotides.\n\n"
        )
        exit()
    elif TargetLen < 1000:
        stderr.write(
            "WARNING - Did you mean to run this on a target region < 1000 nucleotides?\n\n"
        )

    ChangeProbLog = log(ChangeProb / (numseqs - 1))
    ChangeProbSkip = (
        ChangeProbLog / 2
    )  # jumping in and then out of the skip state counts as 1 jump
    SameProbSkip = (
        ChangeProbLog / 30
    )  # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump

    # precomputes number of rows in matrices
    rows: int = len(Subfams)
    cols: int = 0  # assign cols in FillAlignMatrix

    # precomputes consensus seq length for PrintResultsViz()
    ConsensusLengths: Dict[str, int] = {}
    if outfile_viz:
        for i in range(1, len(Flanks)):
            if Strands[i] == "+":
                ConsensusLengths[Subfams[i]] = ConsensusStops[i] + Flanks[i]
            else:
                ConsensusLengths[Subfams[i]] = ConsensusStarts[i] + Flanks[i]

    if TR:
        # add all TR starts and stops
        for rep in TandemRepeats:
            Starts.append(rep['Start'])  # TR start index
            Stops.append(rep['Start'] + rep['Length'] - 1)  # TR stop index

    StartAll, StopAll = pad_sequences(
        ChunkSize, Starts, Stops, SubfamSeqs, ChromSeqs
    )

    (cols, AlignMatrix) = fill_align_matrix(
        StartAll,
        ChunkSize,
        GapExt,
        GapInit,
        SkipAlignScore,
        SubfamSeqs,
        ChromSeqs,
        Starts,
        SubMatrix,
    )

    (
        NonEmptyColumns,
        ActiveCells,
        ConsensusMatrix,
    ) = fill_consensus_position_matrix(
        cols,
        rows,
        StartAll,
        SubfamSeqs,
        ChromSeqs,
        Starts,
        Stops,
        ConsensusStarts,
        Strands,
    )

    if TR:
        TR_row_index_start = rows
        RepeatScores = calculate_repeat_scores(TandemRepeats, ChunkSize, StartAll, rows, ActiveCells,
                                        AlignMatrix, ConsensusMatrix)
        # add skip states for repeat cols
        for tr_col in RepeatScores:
            AlignMatrix[0, tr_col] = float(SkipAlignScore)
            ConsensusMatrix[0, tr_col] = 0

        for rep in TandemRepeats:
            Subfams.append("Tandem Repeat")
            Strands.append("+")
            rows += 1

        ConfidenceMatrix = fill_confidence_matrix_tr(Lamb, infile_prior_counts, NonEmptyColumns, SubfamCounts, Subfams,
                                                     ActiveCells, RepeatScores, AlignMatrix)
        max_tr_col = max(RepeatScores)
        max_align_col = max(NonEmptyColumns)  # last NonEmptyColumn is not always the max col index
        max_col_index: int = max(max_tr_col, max_align_col)
        # check if TR columns were added after last alignment
        # TR cols before alignments were accounted for in PadSeqs
        if max_col_index + 1 > cols:
            cols = max_col_index + 1

        for tr_col in RepeatScores:
            col_set = set(NonEmptyColumns)  # to only add new columns
            col_set.add(tr_col)
            NonEmptyColumns = list(col_set)

    ConfidenceMatrix = fill_confidence_matrix(
        Lamb,
        infile_prior_counts,
        NonEmptyColumns,
        SubfamCounts,
        Subfams,
        ActiveCells,
        AlignMatrix,
    )

    SupportMatrix = fill_support_matrix(
        rows,
        ChunkSize,
        StartAll,
        NonEmptyColumns,
        Starts,
        Stops,
        ConfidenceMatrix,
    )

    collapsed_matrices = collapse_matrices(
        rows,
        NonEmptyColumns,
        Subfams,
        Strands,
        ActiveCells,
        SupportMatrix,
        ConsensusMatrix,
    )

    SupportMatrixCollapse = collapsed_matrices.support_matrix
    SubfamsCollapse = collapsed_matrices.subfamilies
    ActiveCellsCollapse = collapsed_matrices.active_rows
    ConsensusMatrixCollapse = collapsed_matrices.consensus_matrix
    StrandMatrixCollapse = collapsed_matrices.strand_matrix
    SubfamsCollapseIndex = collapsed_matrices.subfamily_indices

    # if command line option included to output support matrix for heatmap
    if outfile_heatmap:
        with open(outfile_heatmap, "w") as outfile:
            print_matrix_support(
                cols,
                StartAll,
                ChromStart,
                SupportMatrixCollapse,
                SubfamsCollapse,
                file=outfile,
            )

    (
        ProbMatrixLastColumn,
        OriginMatrix,
        SameSubfamChangeMatrix,
    ) = fill_probability_matrix(
        SameProbSkip,
        SameProbLog,
        ChangeProbLog,
        ChangeProbSkip,
        NonEmptyColumns,
        collapsed_matrices,
    )

    # IDs for each nucleotide will be assigned during DP backtrace
    IDs = [0] * cols

    (ID, ChangesPosition, Changes) = get_path(
        ID,
        NonEmptyColumns,
        IDs,
        ChangesOrig,
        ChangesPositionOrig,
        NonEmptyColumnsOrig,
        SubfamsCollapse,
        ProbMatrixLastColumn,
        ActiveCellsCollapse,
        OriginMatrix,
        SameSubfamChangeMatrix,
    )

    # keep the original annotation for reporting results
    ChangesOrig = Changes.copy()
    ChangesPositionOrig = ChangesPosition.copy()
    NonEmptyColumnsOrig = NonEmptyColumns.copy()

    # Finding inserted elements and stitching original elements
    # Steps-
    # 1. Find confidence for nodes
    # 2. Create path graph - Find alternative paths through graph and add those edges
    # 3. Extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
    # 5. Annotate again with removed nodes
    #   ** stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left
    prev_num_nodes: int = 0
    count: int = 0
    while True:
        count += 1
        NumNodes = len(Changes)

        # breakout of loop if there are 2 or less nodes left
        if NumNodes <= 2 or NumNodes == prev_num_nodes:
            break

        NodeConfidence.clear() #reuse old NodeConfidence matrix

        NodeConfidence = fill_node_confidence(NumNodes, StartAll, GapInit, GapExt, Lamb, infile_prior_counts,
                                            NonEmptyColumns, Starts, Stops, ChangesPosition, Subfams, SubfamSeqs,
                                            ChromSeqs, SubfamCounts, SubMatrix, RepeatScores, TR_count)
        # store original node confidence for reporting results
        if count == 1:
            NodeConfidenceOrig = NodeConfidence.copy()

        PathGraph.clear() #reuse old PathGraph
        # Update for TR's - no alternative edges
        PathGraph = fill_path_graph(NumNodes, NonEmptyColumns, Changes, ChangesPosition,
                                  ConsensusMatrixCollapse, StrandMatrixCollapse, NodeConfidence, SubfamsCollapseIndex)

        # test to see if there are nodes in the graph that have more than one incoming or outgoing edge,
        # if so - keep looping, if not - break out of the loop
        # if they are all 0, break out of the loop
        test: bool = False
        j: int = 0
        while j < NumNodes:
            i: int = 0
            while i < j - 1:
                if PathGraph[i * NumNodes + j] == 1:
                    test = True
                i += 1
            j += 1

        if not test:
            break

        cols = extract_nodes(
            cols, NumNodes, NonEmptyColumns, ChangesPosition, PathGraph
        )

        # run DP calculations again with nodes corresponding to inserted elements removed
        # ignores removed nodes because they are no longer in NonEmptyColumns
        (
            ProbMatrixLastColumn,
            OriginMatrix,
            SameSubfamChangeMatrix,
        ) = fill_probability_matrix(
            SameProbSkip,
            SameProbLog,
            ChangeProbLog,
            ChangeProbSkip,
            NonEmptyColumns,
            collapsed_matrices,
        )

        Changes.clear()
        ChangesPosition.clear()

        (ID, ChangesPosition, Changes) = get_path(ID, NonEmptyColumns, IDs, ChangesOrig, ChangesPositionOrig,
                                                 NonEmptyColumnsOrig, SubfamsCollapse, ProbMatrixLastColumn, ActiveCellsCollapse,
                                                 OriginMatrix, SameSubfamChangeMatrix)

    if printMatrixPos:
        print_results(
            ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs
        )
    elif printSeqPos:
        print_results_sequence(
            StartAll, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs
        )
    else:
        print_results_chrom(
            StartAll,
            ChromStart,
            ChangesOrig,
            ChangesPositionOrig,
            NonEmptyColumnsOrig,
            IDs,
        )

    if outfile_viz:
        print_results_soda(
            StartAll,
            outfile_viz,
            outfile_conf,
            Chrom,
            ChromStart,
            Subfams,
            ChangesOrig,
            ChangesPositionOrig,
            NonEmptyColumnsOrig,
            ConsensusLengths,
            StrandMatrixCollapse,
            ConsensusMatrixCollapse,
            SubfamsCollapseIndex,
            NodeConfidenceOrig,
        )

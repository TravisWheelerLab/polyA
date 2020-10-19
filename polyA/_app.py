import json
import re
import subprocess
from getopt import getopt
from math import log
from sys import argv, stderr, stdout

from polyA import (
    pad_sequences,
    print_matrix_support,
    print_results,
    print_results_chrom,
    print_results_sequence,
    print_results_soda,
)
from polyA.calc_repeat_scores import CalcRepeatScores
from polyA.collapse_matrices import collapse_matrices
from polyA.extract_nodes import extract_nodes
from polyA.fill_align_matrix import fill_align_matrix
from polyA.fill_confidence_matrix import *
from polyA.fill_confidence_matrix_tr import FillConfidenceMatrixTR
from polyA.fill_consensus_position_matrix import fill_consensus_position_matrix
from polyA.fill_node_confidence import fill_node_confidence
from polyA.fill_path_graph import fill_path_graph
from polyA.fill_probability_matrix import fill_probability_matrix
from polyA.fill_support_matrix import fill_support_matrix
from polyA.get_path import get_path
from polyA.lambda_provider import EaselLambdaProvider
from polyA.load_alignments import load_alignments


def run():
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
    SkipAlignScore: float = 5.0  # FIXME - not sure what number this should be

    StartAll: int = 0  # Reassigned later
    StopAll: int = 0  # Reassigned later
    ID: int = 1111

    infile_prior_counts: str = ""
    outfile_viz: str = ""
    outfile_heatmap: str = ""

    # running ultra
    ultra_path: str = ""
    seq_file: str = ""
    # given ultra output
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
        --seqFile [specify path to genome region file]
        --ultraOutput [specify path to ULTRA output file from region]
    
    OPTIONS
        --help - display help message
        --matrixpos - prints subfam changes in matrix position instead of genomic position
        --sequencepos - prints subfam changes in sequence position instead of genomic position
        --viz outfile - prints output format for SODA visualization
        --heatmap outfile - prints probability file for input into heatmap
    """

    raw_opts, args = getopt(
        argv[1:],
        "",
        [
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
            "seqFile=",
            "ultraOutput=",
            "help",
            "matrixpos",
            "seqpos",
        ],
    )
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
    outfile_heatmap = (
        str(opts["--heatmap"]) if "--heatmap" in opts else outfile_heatmap
    )

    # options for using ultra
    ultra_path = (
        str(opts["--ultraPath"]) if "--ultraPath" in opts else ultra_path
    )
    seq_file = str(opts["--seqFile"]) if "--seqFile" in opts else seq_file
    ultra_output_path = (
        str(opts["--ultraOutput"])
        if "--ultraOutput" in opts
        else ultra_output_path
    )

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
        provider = EaselLambdaProvider(EslPath, infile_matrix)
        Lamb = provider()

    # reads in the score matrix from file and stores in dict that maps 'char1char2' to the score from the
    # input substitution matrix - ex: 'AA' = 8
    SubMatrix: Dict[str, int] = {}

    # add all ambiguity codes just incase matrix doesnt have them
    nucleotide_codes = "AGCTYRWSKMDVHBXN."
    for code in nucleotide_codes:
        for code2 in nucleotide_codes:
            SubMatrix[code + code2] = 0

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
    # user should already have the correct file of seq region from alignments
    if ultra_path and seq_file:
        TR = True
        # run ULTRA
        pop = subprocess.Popen(
            [ultra_path, "-ss", seq_file],
            stdout=subprocess.PIPE,
            universal_newlines=True,
        )
        ultra_out, err = pop.communicate()
        ultra_output = json.loads(ultra_out)
    elif ultra_output_path:
        TR = True
        # load ULTRA output file
        with open(ultra_output_path) as f:
            ultra_output = json.load(f)

    # get TR info from ULTRA output
    TR_count: int = 0
    if TR:
        TandemRepeats = ultra_output["Repeats"]
        if len(TandemRepeats) == 0:
            # no repeats found
            TR = False
        else:
            TR_count = len(TandemRepeats)

    # maps subfam names to genomic prior_count/total_in_genome from input file
    # used during confidence calculations
    SubfamCounts: Dict[str, float] = {}
    PriorTotal: float = 0
    prob_skip = 0.4  # about 60% of genome is TE derived
    prob_tr = 0.06  # about 6% of genome expected to be a tandem repeat
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

    NonEmptyColumns_trailing: List[int] = []
    ActiveCells_trailing: Dict[int, List[int]] = {}

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

    match = re.search(r"(.+):(\d+)-(\d+)", Chroms[1])
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
            Starts.append(rep["Start"])  # TR start index
            Stops.append(rep["Start"] + rep["Length"] - 1)  # TR stop index

    (StartAll, StopAll) = pad_sequences(
        ChunkSize, Starts, Stops, SubfamSeqs, ChromSeqs
    )

    #add skip state pad at start
    AlignMatrix[0, 0] = SkipAlignScore
    NonEmptyColumns.append(0)
    NonEmptyColumns_trailing.append(0)

    (cols, AlignMatrix) = fill_align_matrix(
        Lamb,
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

    # originally NonEmptyColumns and ActiveCells have trailing edge included
    # redo these later to not include trailing edges
    # TODO: should be able to make this faster
    NonEmptyColumns_trailing, ActiveCells_trailing = trailing_edges_info(
        rows, cols, AlignMatrix
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

    ActiveCells[0] = [0]

    #add skip state pad at end
    AlignMatrix[0, cols] = SkipAlignScore
    NonEmptyColumns.append(cols)
    NonEmptyColumns_trailing.append(cols)
    ActiveCells[cols] = [0]
    ActiveCells_trailing[cols] = [0]
    cols += 1

    if TR:
        RepeatScores = CalcRepeatScores(
            TandemRepeats,
            ChunkSize,
            StartAll,
            rows,
            ActiveCells,
            AlignMatrix,
            ConsensusMatrix,
        )

        # add skip states for TR cols
        for tr_col in RepeatScores:
            AlignMatrix[0, tr_col] = float(SkipAlignScore)

        # add TRs to subfams
        for rep in TandemRepeats:
            Subfams.append("Tandem Repeat")
            Strands.append("+")
            rows += 1

        ConfidenceMatrix = FillConfidenceMatrixTR(
            infile_prior_counts,
            NonEmptyColumns,  # FIXME - does this need to be NonEmptyColumns_trailing?
            SubfamCounts,
            Subfams,
            ActiveCells,
            RepeatScores,
            AlignMatrix,
        )

        # check if TR columns were added after last alignment
        # TR cols before alignments were accounted for in PadSeqs
        max_col_index: int = max(RepeatScores)
        if max_col_index + 1 > cols:
            cols = max_col_index + 1

        # add TR cols to NonEmptyColumns
        for tr_col in RepeatScores:
            col_set = set(NonEmptyColumns)  # only add new columns
            col_set.add(tr_col)
            NonEmptyColumns = list(col_set)
        NonEmptyColumns.sort()

    else:
        ConfidenceMatrix = fill_confidence_matrix(
            infile_prior_counts,
            NonEmptyColumns_trailing,
            SubfamCounts,
            Subfams,
            ActiveCells_trailing,
            AlignMatrix,
        )

    # add skip state to consensus matrix
    # wait till after incase NonEmptyColumns is updated by TR stuff
    for j in NonEmptyColumns:
        ConsensusMatrix[0, j] = 0

    SupportMatrix = fill_support_matrix(
        rows,
        ChunkSize,
        StartAll,
        NonEmptyColumns,
        Starts,
        Stops,
        ConfidenceMatrix,
        ConsensusMatrix,
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
    rows = collapsed_matrices.row_num_update

    if TR:
        # give different TRs consensus positions that don't allow them to be stitched
        tr_consensus_pos = 1000000
        prev_tr_col = 0
        for tr_col in RepeatScores:
            if prev_tr_col != tr_col - 1:
                tr_consensus_pos -= 500
            ConsensusMatrixCollapse[rows - 1, tr_col] = tr_consensus_pos
            prev_tr_col = tr_col

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

        NodeConfidence.clear()  # reuse old NodeConfidence matrix

        NodeConfidence = fill_node_confidence(
            NumNodes,
            StartAll,
            GapInit,
            GapExt,
            Lamb,
            infile_prior_counts,
            NonEmptyColumns,
            Starts,
            Stops,
            ChangesPosition,
            Subfams,
            SubfamSeqs,
            ChromSeqs,
            SubfamCounts,
            SubMatrix,
            RepeatScores,
            TR_count,
        )
        # store original node confidence for reporting results
        if count == 1:
            NodeConfidenceOrig = NodeConfidence.copy()

        PathGraph.clear()  # reuse old PathGraph
        # Update for TR's - no alternative edges
        PathGraph = fill_path_graph(
            NumNodes,
            NonEmptyColumns,
            Changes,
            ChangesPosition,
            ConsensusMatrixCollapse,
            StrandMatrixCollapse,
            NodeConfidence,
            SubfamsCollapseIndex,
        )

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

        (ID, ChangesPosition, Changes) = get_path(
            ID,
            NonEmptyColumns,
            IDs,
            Changes,
            ChangesPosition,
            NonEmptyColumns,
            SubfamsCollapse,
            ProbMatrixLastColumn,
            ActiveCellsCollapse,
            OriginMatrix,
            SameSubfamChangeMatrix,
        )
        prev_num_nodes = NumNodes

        if TR:
            for subfam in Changes:
                assert (
                    subfam != "Tandem Repeat"
                ), "can't add alternative edges to TRs"

    # prints results
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
            IDs,
        )

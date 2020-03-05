from getopt import getopt
from math import inf, log
import re
from sys import argv, stderr, stdout
from typing import Dict, List, Tuple, Union

from polyA.load_alignments import load_alignments

def PrintMatrixHashCollapse(column: int, Hash: Dict) -> None:
    stdout.write("\t")

    j: int = 0
    while j < column:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    j: int = 0
    while j < column:
        if ("skip", j) in Hash:
            stdout.write(Hash["skip", j])
        else:
            stdout.write(f"{-inf:5.3g}")
        stdout.write("\t")
        j += 1
    stdout.write("\n")

    for k in SubFamsCollapse:
        if k != "skip":
            stdout.write(f"{k:<20s}\t")
            j: int = 0
            while j < column:
                if (i, j) in Hash:
                    stdout.write(Hash[i, j])
                else:
                    stdout.write(f"{-inf:5.3g}")
                stdout.write("\t")
                j += 1
            stdout.write("\n")


def PrintMatrixHash(column: int, Hash: Dict) -> None:
    stdout.write("\t")
    j: int = 0
    while j < column:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    i: int = 0
    while i < rows:
        stdout.write(f"{SubFams[i]:<20s}\t")
        j: int = 0
        while j < column:
            if (i, j) in Hash:
                stdout.write(f"{Hash[i, j]}")
            else:
                stdout.write(f"{-inf:5.3g}")
            stdout.write("\t")
            j += 1
        stdout.write("\n")
        i += 1


gapInit: int = -25
gapExt: int = -5
lamb: float = 0.1227  # From command line
chunksize: int = 30
sameProbLog: float = log(1 - (10 ** -45))
changeProb: float = 10 ** -45
changeProbLog: float = 0.0  # Reassigned later
skipAlignScore: int = 30
startall: int = 0  # Reassigned later
stopall: int = 0  # Reassigned later

help: bool = False  # Reassigned later
prin: bool = False  # Reassigned later
printMatrixPos: bool = False  # Reassigned later

helpMessage: str = f"""
usage: {argv[0]} alignFile matrixFile\n
ARGUMENTS
    --gapInit[-25]
    --getExt[-5]
    --lambda [will calc from matrix if not included]
    --segmentsize[30]
    --changeprob[1e-45]

OPTIONS
    --help - display help message
    --printmatrices - output all matrices
    --matrixpos - prints subfam changes in matrix position instead of genomic position
"""

raw_opts, args = getopt(argv[1:], "", [
    "gapInit=",
    "gapExt=",
    "lambda=",
    "segmentsize=",
    "changeprob=",

    "help",
    "printmatrices",
    "matrixPos",
])
opts = dict(raw_opts)

gapInit = int(opts["gapInit"]) if "gapInit" in opts else gapInit
gapExt = int(opts["gapExt"]) if "gapExt" in opts else gapExt
lamb = float(opts["lambda"]) if "lambda" in opts else lamb
chunksize = int(opts["segmentsize"]) if "segmentsize" in opts else chunksize
changeProb = float(opts["changeprob"]) if "changeprob" in opts else changeProb
help = "help" in opts
printt = "printmatrices" in opts
printMatrixPos = "matrixPos" in opts

if help:
    print(helpMessage)
    exit(0)

infile: str = args[0]
infile_matrix: str = args[1]

# Other open was moved down to where we load the alignments file
with open(infile_matrix) as _infile_matrix:
    in_matrix: List[str] = _infile_matrix.readlines()

# We require a --lambda or provide a default so there's no need to run easel

CharPos: Dict[str, int] = {}
line = in_matrix[0]
line = re.sub(r"^\s+", "", line)
chars = re.split(r"\s+", line)
for i in range(len(chars)):
    CharPos[chars[i]] = i
subMatrixCols: int = len(chars)

SubMatrix: Dict[int, int] = {}
count: int = 0
for line in in_matrix[1:]:
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    subScores = re.split(r"\s+", line)
    for i in range(len(subScores)):
        SubMatrix[count*subMatrixCols+i] = int(subScores[i])
    count += 1

SubFams: List[str] = []
Scores: List[int] = []
Strands: List[str] = []
Starts: List[int] = []
Stops: List[int] = []
ConsensusStarts: List[int] = []
ConsensusStops: List[int] = []
SubfamSeqs: List[str] = []
ChromSeqs: List[str] = []

AlignHash: Dict[Tuple[int, int], int] = {}
ConfHash: Dict[Tuple[int, int], float] = {}
SupportHash: Dict[Tuple[int, int], float] = {}
ProbHash: Dict[Tuple[str, int], float] = {}
OriginHash: Dict[Tuple[str, int], str] = {}
ConsensusHash: Dict[Tuple[int, int], int] = {}

subfampath: List[str] = []

RemoveStarts: List[int] = []
RemoveStops: List[int] = []

Changes: List[str] = []
ChangesPos: List[int] = []

# SubFams.append("skip")

# Scores.append(-1)
# Starts.append(-1)
# Strands.append("")
# Stops.append(-1)
# ConsensusStarts.append(-1)
# ConsensusStops.append(-1)
# SubfamSeqs.append("")
# ChromSeqs.append("")

numseqs: int = 0
with open(infile) as _infile:
    alignments = load_alignments(_infile)
    for alignment in alignments:
        numseqs += 1

        SubFams.append(alignment.subfamily)
        Scores.append(alignment.score)
        # TODO: If we need this, add it back in load_alignments
        Strands.append(alignment.strand)
        Starts.append(alignment.start)
        Stops.append(alignment.stop)
        ConsensusStarts.append(alignment.consensus_start)
        ConsensusStops.append(alignment.consensus_stop)
        SubfamSeqs.append(alignment.subfamily_sequence)
        ChromSeqs.append(alignment.sequence)

print(f"stuff = {ChromSeqs[1]}")

changeProbLog = log(changeProb / (numseqs - 1))

IDs: List[int] = []
Active: List[bool] = []
nodeStarts: List[int] = []
nodeStops: List[int] = []
idnum = 1111
for sf in SubFams:
    IDs.append(idnum)
    Active.append(False)
    nodeStarts.append(0)
    nodeStops.append(0)

    idnum += 1

rows: int = len(SubFams)
cols: int = 0

Position = [0] * len(SubfamSeqs)


def CalcScore(seq1: str, seq2: str, lastpreva: str, lastprevb: str):
    chunkscore: int = 0

    if seq1[0] == "-":
        if lastpreva == "-":
            chunkscore += gapExt
        else:
            chunkscore += gapInit
    elif seq2[0] == "-":
        if lastprevb == "-":
            chunkscore += gapExt
        else:
            chunkscore += gapInit
    elif seq1[0] == "." or seq2[0] == ".":
        chunkscore = chunkscore
    else:
        chunkscore += SubMatrix[CharPos[seq1[0]] * subMatrixCols + CharPos[seq2[0]]]

    j: int = 1
    while j < len(seq1):
        if seq1[j] == "-":
            if seq1[j - 1] == "-":
                chunkscore += gapExt
            else:
                chunkscore += gapInit
        elif seq2[j] == "-":
            if seq2[j - 1] == "-":
                chunkscore += gapExt
            else:
                chunkscore += gapInit
        elif seq1[j] == "." or seq2[j] == ".":
            chunkscore = chunkscore
        else:
            if seq1[j] in CharPos and seq2[j] in CharPos:
                chunkscore += SubMatrix[CharPos[seq1[j]] * subMatrixCols + CharPos[seq2[j]]]

        j += 1

    return chunkscore


def FillAlignScoreMatrix(subfams: List[str], chroms: List[str]):
    global cols

    i: int = 1
    while i < len(ChromSeqs):
        subfam1: str = subfams[i]
        chrom1: str = chroms[i]

        j: int = Starts[i] - startall
        index: int = j

        offset: int = chunksize
        alignscore: int = 0

        tempindex: int = j
        tempcount: int = 0

        while tempcount < chunksize:
            if chrom1[tempindex] != "-":
                tempcount += 1
            tempindex += 1

        offset = tempindex - j
        prevoffset = offset

        ChromSlice: str = chrom1[j:j + offset]
        SubfamSlice: str = subfam1[j:j + offset]

        alignscore = CalcScore(SubfamSlice, ChromSlice, "", "")
        AlignHash[i, index] = alignscore

        if cols < index:
            cols = index + 1

        index += 1

        # TODO: Make sure these bounds are right since Python indexing is different
        while j + offset + 1 < len(chrom1):
            tempindex = j
            tempcount = 0

            while tempcount < chunksize:
                if chrom1[tempindex + 1] != "-":
                    tempcount += 1
                tempindex += 1

            offset = tempindex - j

            if chrom1[j + 1] != "-":
                if chrom1[j + 1] != "." and subfam1[j + 1] != ".":
                    if prevoffset != offset:
                        ChromSlice = chrom1[j+1:j+offset+1]
                        SubfamSlice = subfam1[j+1:j + offset + 1]
                        alignscore = CalcScore(SubfamSlice, ChromSlice, subfam1[j], chrom1[j])
                    else:
                        if subfam1[j] == "-":
                            if subfam1[j - 1] == "-":
                                alignscore = alignscore - gapExt
                            else:
                                alignscore = alignscore - gapInit
                        else:
                            alignscore = alignscore - SubMatrix[CharPos[subfam1[j]]*subMatrixCols+CharPos[chrom1[j]]]

                        if subfam1[j + offset] == "-":
                            if subfam1[j + offset - 1] == "-":
                                alignscore = alignscore + gapExt
                            else:
                                alignscore = alignscore + gapInit
                        elif subfam1[j + offset] == "." or chrom1[j + offset] == ".":
                            alignscore = alignscore
                        else:
                            alignscore = alignscore + SubMatrix[CharPos[subfam1[j + offset]] * subMatrixCols + CharPos[chrom1[j + offset]]]

                    if alignscore <= 0:
                        AlignHash[i, index] = 1
                    else:
                        AlignHash[i, index] = alignscore

                if cols < index:
                    cols = index + 1
                index += 1

            j += 1
            prevoffset = offset

        i += 1


def Edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    minStart: int = starts[1]
    maxStop: int = stops[1]

    i: int = 1
    while i < len(starts):
        if starts[i] < minStart:
            minStart = starts[i]
        if stops[i] > maxStop:
            maxStop = stops[i]

        i += 1

    return minStart, maxStop


def padSeqs(start, stop, subfamseq, chromseq):
    global startall, stopall

    startall, stopall = Edges(start, stop)

    i: int = 1
    while i < len(SubfamSeqs):
        leftpad = start[i] - startall
        rightpad = stopall - stop[i]

        chromseq[i] = ("." * leftpad) + f"{chromseq[i]}" + ("." * rightpad)
        subfamseq[i] = ("." * leftpad) + f"{subfamseq[i]}" + ("." * rightpad)

        i += 1


padSeqs(Starts, Stops, SubfamSeqs, ChromSeqs)

FillAlignScoreMatrix(SubfamSeqs, ChromSeqs)

def FillConsensusPosMatrix(consensus: Dict[Tuple[int, int], int],
                           subfams: List[str],
                           chroms: List[str],
                           consensusstart: List[int],
                           consensusstop: List[int]):
    i: int = 1
    while i < rows:
        SubfamArray: str = subfams[i]
        ChromArray: str = chroms[i]

        consensuspos: int = 0
        if Strands[i] == "+":
            consensuspos = consensusstart[i] - 1
            matrixpos: int = 0
            j: int = 0
            while j < len(SubfamArray):
                if SubfamArray[j] != ".":
                    if SubfamArray[j] != "-":
                        consensuspos += 1

                    # TODO: What does the -> in the perl code do?
                    consensus[i, matrixpos] = consensuspos

                if ChromArray[j] != "-":
                    matrixpos += 1
                j += 1
        else:
            consensuspos = consensusstart[i] + 1
            matrixpos: int = 0
            j: int = 0
            while j < len(SubfamArray):
                if SubfamArray[j] != ".":
                    if SubfamArray[j] != "-":
                        consensuspos -= 1
                    consensus[i, matrixpos] = consensuspos

                if ChromArray[j] != "-":
                    matrixpos += 1
                j += 1

        # FIXME: Delete later, just for testing (see perl code, line 985)
        if consensusstop[i] != consensuspos:
            stderr.write("\n\nERROR - consensus seq positions not correct\n\n")

        i += 1


FillConsensusPosMatrix(ConsensusHash, SubfamSeqs, ChromSeqs, ConsensusStarts, ConsensusStops)

Columns: List[int] = []
j: int = 0
while j < cols:
    empty: bool = True
    i: int = 0
    while i < rows:
        if (i, j) in AlignHash:
            empty = False
            i = rows
        else:
            i += 1

    if not empty:
        Columns.append(j)

    AlignHash[0, j] = skipAlignScore

    j += 1

i: int = 0
while i < chunksize - 1:
    Columns.append(j)
    i += 1
    j += 1


def ConfidenceCM(lamb: float, region: List[float]) -> str:
    confidenceString: str = ""

    ScoreTotal: int = 0
    for Score in region:
        if Score > 0:
            convertedScore = Score * lamb
            ScoreTotal += 2 ** convertedScore

    for Score in region:
        if Score > 0:
            convertedScore = Score * lamb
            confidence = ((2 ** convertedScore) / ScoreTotal)

            if confidenceString != "":
                confidenceString += f" {confidence}"
            else:
                confidenceString = f"{confidence}"
        else:
            if confidenceString != "":
                confidenceString += " 0"
            else:
                confidenceString = "0"

    return confidenceString


def FillConfScoreMatrix(alignhash: Dict[Tuple[int, int], int], confhash: Dict[Tuple[int, int], float]):
    i: int = 0
    while i < len(Columns) - chunksize + 1:
        col: int = Columns[i]
        temp: List[int] = []
        row: int = 0
        while row < rows:
            if (row, col) in alignhash:
                temp.append(alignhash[row, col])
            else:
                temp.append(0)
            row += 1

        confidenceTemp = ConfidenceCM(lamb, temp).split(" ")
        row: int = 0
        while row < rows:
            if confidenceTemp[row] != '0':
                confhash[row, col] = float(confidenceTemp[row])
            row += 1
        i += 1

FillConfScoreMatrix(AlignHash, ConfHash)

cols = cols + chunksize - 1


def FillSupportMatrix(supporthash: Dict[Tuple[int, int], float], alignhash: Dict[Tuple[int, int], int], confhash: Dict[Tuple[int, int], float]):
    i: int = 0
    while i < rows:
        tempcol: int = -1
        col: int = 0
        while col < len(Columns) - chunksize + 1:
            j = Columns[col]

            if (i, j) in confhash:
                num: int = j
                summ: float = 0.0
                numsegments: int = 0
                while num >= 0 and num >= j - chunksize + 1:
                    if (i, num) in confhash:
                        summ = summ + confhash[i, num]
                        numsegments += 1
                    num -= 1

                if numsegments > 0:
                    supporthash[i, j] = summ / numsegments
            col += 1

        j: int = cols - chunksize
        while j <= cols:
            num: int = j
            summ: float = 0
            numsegments: int = 0
            while num >= 0 and num >= j - chunksize + 1:
                if (i, num) in confhash:
                    summ = summ + confhash[i, num]
                    numsegments += 1
                num -= 1

            if numsegments > 0:
                supporthash[i, j] = summ / numsegments
            j += 1

        i += 1


FillSupportMatrix(SupportHash, AlignHash, ConfHash)

SupportHashCollapse: Dict[Tuple[str, int], int] = {}
ActiveCellsCollapse: Dict[int, List[str]] = {}
SubFamsCollapse: Dict[str, int] = {}
ConsensusHashCollapse: Dict[Tuple[str, int], int] = {}
StrandHashCollapse: Dict[Tuple[str, int], str] = {}

j = 0
while j < cols:
    activecols: List[str] = []
    ActiveCellsCollapse[j] = activecols

    maxx: float = 0
    maxrow: int = -1
    i: int = 1
    while i < rows:
        if (i, j) in SupportHash and (i, j) in ConsensusHash:
            ConsensusHashCollapse[SubFams[i], j] = ConsensusHash[i, j]
            StrandHashCollapse[SubFams[i], j] = Strands[i]

            if (SubFams[i], j) in SupportHashCollapse:
                SupportHashCollapse[SubFams[i], j] = SupportHashCollapse[SubFams[i], j]

                if SupportHash[i, j] > maxx:
                    maxx = SupportHash[i, j]
                    maxrow = i
            else:
                SupportHashCollapse[SubFams[i], j] = SupportHash[i, j]
                ActiveCellsCollapse[j].append(SubFams[i])

                maxx = SupportHash[i, j]
                maxrow = i

        i += 1

    if maxrow > -1:
        ConsensusHashCollapse[SubFams[maxrow], j] = ConsensusHash[maxrow, j]
        StrandHashCollapse[SubFams[maxrow], j] = Strands[maxrow]

    j += 1

i: int = 0
while i < rows:
    SubFamsCollapse[SubFams[i]] = 0
    i += 1

rows = len(SubFamsCollapse)

# line 381
for k in SubFamsCollapse:
    ProbHash[k, 0] = 0


def FillProbMatrix(probhash: Dict[Tuple[str, int], float], supporthash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str]):
    j: int = 1
    col: int = 1
    while col < len(Columns):
        if col in Columns:
            j = Columns[col]
        else:
            j += 1

        for i in ActiveCellsCollapse[j]:
            max: float = -inf
            maxindex: str = ''
            supportlog: float = log(supporthash[i, j])

            # TODO: Bug around here, trying to populate some nonsense and it hasn't populated the previous nonsense
            # KAITLIN!!!!!

            for row in ActiveCellsCollapse[Columns[col - 1]]:
                score: float = -1
                there: bool = False

                if col - 1 in Columns:
                    score = supportlog + probhash[row, col - 1]
                    there = True
                else:
                    score = supportlog + probhash[row, j - 1]
                    there = True

                if there:
                    if row == i:
                        score = score + sameProbLog
                    else:
                        score = score + changeProbLog

                    if score > max:
                        max = score
                        maxindex = row

            probhash[i, j] = max
            originhash[i, j] = maxindex

        col += 1


FillProbMatrix(ProbHash, SupportHashCollapse, OriginHash)
PrintMatrixHashCollapse(cols, ProbHash)
exit()


def GetPath(probhash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str], subfams: List[str]) -> List[str]:
    maxxx: float = -inf
    maxindex: str = ''
    for i1 in ActiveCellsCollapse[cols - 1]:
        print(i1)
        print(cols - 1)
        print(probhash[i1, cols - 1])
        if maxxx < probhash[i1, cols - 1]:
            maxxx = probhash[i1, cols - 1]
            maxindex = i1

    subfamorder: List[str] = []
    prev: str = originhash[maxindex, cols - 1]
    i: int = cols - 1

    subfamorder.append(maxindex)

    col: int = len(Columns)
    while col > 0:
        if col in Columns:
            i = Columns[col]
        else:
            i -= 1

        if col in Columns:
            if (prev, Columns[col - 1]) in originhash and (prev, i) in originhash:
                subfamorder.append(prev)
                prev = originhash[prev, Columns[col - 1]]
            else:
                subfamorder.append('')
        else:
            if (prev, i - 1) in originhash and (prev, i) in originhash:
                subfamorder.append(prev)
                prev = originhash[prev, i - 1]
            else:
                subfamorder.append('')
        col -= 1

    subfamorder.reverse()
    return subfamorder


subfampath = GetPath(ProbHash, OriginHash, SubFams)


# FIXME: Add parameter types
def GetChanges(changes, changespos):
    prev: str = "skip"
    if subfampath[0] != "":
        prev = subfampath[0]

    i: int = 1
    while i < len(subfampath):
        curr_subfam: str = "skip"
        if subfampath[i] != "":
            curr_subfam = subfampath[i]

            # FIXME: Need the Python equivalent of /g on the regex
            match: Union[re.Match, None] = re.search(r"(.+?)_.+", curr_subfam)
            if match is not None:
                curr_subfam = match.groups()[0]

            if curr_subfam != prev:
                changespos.append(i)
                changes.append(subfampath[i])

            prev = curr_subfam

            i += 1


GetChanges(Changes, ChangesPos)


def PrintChanges(changes, changespos):
    i: int = 0
    while i < len(changes):
        stderr.write(Columns[changespos[i]])
        stderr.write("\t")
        stderr.write(f"{changes[i]}\n")


PrintChanges(Changes, ChangesPos)


def PrintAllMatrices():
    stdout.write("Align Scores\n")
    PrintMatrixHash(cols, AlignHash)
    stdout.write("confidence\n")
    PrintMatrixHash(cols, ConfHash)
    stdout.write("support\n")
    PrintMatrixHash(cols, SupportHash)
    stdout.write("prob\n")
    PrintMatrixHash(cols, ProbHash)
    stdout.write("origin\n")
    PrintMatrixHash(cols, OriginHash)


stderr.write("\n")
if printt:
    PrintAllMatrices()


numnodes: int = 0
NodeConfidenceDict: Dict[Tuple[str, int], float] = {}
pathGraph: List[int] = []
total: int = 0
loop: int = 1

count: int = 0
while (True):
    count += 1
    numnodes = len(Changes)
    if (numnodes <= 2):
        pass # `last;` in perl


    def NodeConfidence(nodeconfidence: Dict[Tuple[str, int], float], subfamseqs: List[str], chromseqs,
                       changespos: List[int]):
        nodeconfidence_temp: List[float] = [0 for _ in range(rows * numnodes)]

        j: int = 1
        while j < len(SubFams):
            b: int = changespos[0] - 1
            e: int = changespos[1]
            subfam: str = subfamseqs[j][b:e]
            chrom: str = chromseqs[j][b:e]
            alignscore: float = CalcScore(subfam, chrom, '', '')
            nodeconfidence_temp[j * numnodes + 0] = alignscore

        i: int = 0
        while i < numnodes - 1:
            j: int = 1
            while j < len(SubFams):
                b: int = changespos[i] - 1
                e: int = changespos[i + 1]
                subfam: str = subfamseqs[j][b:e]
                chrom: str = chromseqs[j][b:e]

                lastpreva: str = subfamseqs[j][changespos[i + 1] - 1]
                lastprevb: str = chromseqs[j][changespos[i + 1] - 1]
                alignscore: float = CalcScore(subfam, chrom, lastpreva, lastprevb)
                nodeconfidence_temp[j * numnodes + i] = alignscore

                j += 1
            i += 1

        j: int = 1
        while j < len(SubFams):
            subfam: str = subfamseqs[j][changespos[-1] - 1:]
            chrom: str = chromseqs[j][changespos[-1] - 1:]
            alignscore: float = CalcScore(subfam, chrom, '', '')
            nodeconfidence_temp[j * numnodes + numnodes - 1] = alignscore
            j += 1

        for j in range(numnodes):
            temp: List[float] = []
            for i in range(len(SubFams)):
                temp.append(nodeconfidence_temp[i * numnodes + j])

            confidenceTemp: List[float] = [float(x) for x in ConfidenceCM(lamb, temp).split(' ')]
            for i in range(len(SubFams)):
                nodeconfidence_temp[i * numnodes + j] = confidenceTemp[i]

        for j in range(numnodes):
            for i in range(len(SubFams)):
                if (SubFams[i], j) in nodeconfidence:
                    nodeconfidence[SubFams[i], j] += nodeconfidence_temp[i * numnodes + j]
                else:
                    nodeconfidence[SubFams[i], j] = nodeconfidence_temp[i * numnodes + j]


    NodeConfidenceDict.clear()
    NodeConfidence(NodeConfidenceDict, SubfamSeqs, ChromSeqs, ChangesPos)


    def FillPathGraph(pathgraph: List[int]):
        for i in range(numnodes * numnodes):
            pathgraph.append(0)

        for i in range(numnodes - 1):
            pathgraph[i * numnodes + i + 1] = 1

        for j in range(numnodes):
            sinkSubfam: str = Changes[j]
            sinkSubfamStart: int = ConsensusHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]
            sinkStrand: str = StrandHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]

            for i in range(j - 1):
                for sourceSubfam in SubFamsCollapse:
                    sourceConf = NodeConfidenceDict[sourceSubfam, i]
                    if (sourceSubfam, Columns[ChangesPos[i + 1]] - 1) in ConsensusHashCollapse:
                        sourceSubfamStop = ConsensusHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]] - 1]
                        sourceStrand = StrandHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]] - 1]

                        if sinkStrand == '+' and sinkStrand == sourceStrand:
                            if sinkSubfam == sourceSubfam and sourceConf >= 0.5:
                                if sourceSubfamStop <= sinkSubfamStart + 50:
                                    pathgraph[i * numnodes + j] = 1
                        elif sinkStrand == '-' and sinkStrand == sourceStrand:
                            if sinkSubfam == sourceSubfam and sourceConf >= 0.5:
                                if sourceSubfamStop >= sinkSubfamStart + 50:
                                    pathgraph[i * numnodes + j] = 1


    pathGraph.clear()
    FillPathGraph(pathGraph)

    test: bool = False
    j: int = 0
    while j < numnodes:
        i: int = 0
        while i < j - 1:
            if (pathGraph[i * numnodes + j] == 1):
                test = True
            i += 1
        j += 1

    if not test:
        pass # `last;` in perl

    RemoveStarts.clear()
    RemoveStops.clear()


    def ExtractNodes(removestarts, removestops, changespos, pathgraph, numnodes: int):
        global cols

        RemoveNodes: List[bool] = [False for _ in range(numnodes)]

        NumEdgesIn: List[int] = [0 for _ in range(numnodes)]
        NumEdgesOut: List[int] = [0 for _ in range(numnodes)]

        for i in range(numnodes):
            for j in range(numnodes):
                NumEdgesIn[j] += pathgraph[i * numnodes + j]
                NumEdgesOut[i] += pathgraph[i * numnodes + j]

        for i in range(numnodes - 1):
            if NumEdgesIn[i] <= 1 and NumEdgesOut[i] <= 1:
                removestarts.append(changespos[i])
                removestops.append(changespos[i + 1])

                RemoveNodes[i] = True

        if NumEdgesIn[numnodes - 1] <= 1 and NumEdgesOut[numnodes - 1] <= 1:
            removestarts.append(changespos[numnodes - 1])
            removestops.append(cols)

            RemoveNodes[numnodes - 1] = True

        i: int = numnodes - 1
        while RemoveNodes[i]:
            cols = changespos[i] - 1
            i -= 1


    ExtractNodes(RemoveStarts, RemoveStops, ChangesPos, pathGraph, numnodes)

    total: int = 0
    i: int = 0
    while i < len(RemoveStops):
        del Columns[RemoveStarts[i] - total:RemoveStops[i] - RemoveStarts[i] + 1]

        total += (RemoveStops[i] - RemoveStarts[i])

        i += 1

    FillProbMatrix(ProbHash, SupportHashCollapse, OriginHash)

    Changes.clear()
    ChangesPos.clear()
    subfampath.clear()

    subfampath = GetPath(ProbHash, OriginHash, SubFams)
    GetChanges(Changes, ChangesPos)
    PrintChanges(Changes, ChangesPos)
    stderr.write("\n")

    # TODO Remove before running
    break

def PrintAllMatrices():
    stdout.write("Align Scores\n")
    PrintMatrixHash(cols, AlignHash)
    stdout.write("confidence\n")
    PrintMatrixHash(cols, ConfHash)
    stdout.write("support\n")
    PrintMatrixHash(cols, SupportHash)
    stdout.write("prob\n")
    PrintMatrixHash(cols, ProbHash)
    stdout.write("origin\n")
    PrintMatrixHash(cols, OriginHash)


def PrintResults():
    stderr.write("\nstart\tstop\tid\tsubfam\n")
    stderr.write("-------------------------------------\n")

    i: int = 1
    while i < len(Active):
        if Active[i]:
            stderr.write(f"{nodeStarts[i]}\t{nodeStops[i]}\t{IDs[i]}\t{SubFams[i]}\n")


def PrintPathGraph():
    i: int = 0
    while i < numnodes:
        stderr.write(f"{Changes[i]} ")
    stderr.write("\n")

    i: int = 0
    while i < numnodes:
        j: int = 0
        while j < numnodes:
            stderr.write(f"{pathGraph[i*numnodes+j]}\t")
        stderr.write(f"{Changes[i]}\n")


# FIXME: Add parameter types
def GetChanges(changes, changespos):
    prev: str = "skip"
    if subfampath[0] != "":
        prev = subfampath[0]

    i: int = 1
    while i < len(subfampath):
        curr_subfam: str = "skip"
        if subfampath[i] != "":
            curr_subfam = subfampath[i]

            # FIXME: Need the Python equivalent of /g on the regex
            match: Union[re.Match, None] = re.search(r"(.+?)_.+", curr_subfam)
            if match is not None:
                curr_subfam = match.groups()[0]

            if curr_subfam != prev:
                changespos.append(i)
                changes.append(subfampath[i])

            prev = curr_subfam

            i += 1


def PrintChanges(changes, changespos):
    i: int = 0
    while i < len(changes):
        stderr.write(Columns[changespos[i]])
        stderr.write("\t")
        stderr.write(f"{changes[i]}\n")


def FillAlignScoreMatrix(subfams: List[str], chroms: List[str]):
    global cols

    i: int = 1
    while i < len(ChromSeqs):
        subfam1: str = subfams[i]
        chrom1:  str = chroms[i]

        j: int = Starts[i] - startall
        index: int = j

        offset: int = chunksize
        alignscore: int = 0

        tempindex: int = j
        tempcount: int = 0

        while tempcount < chunksize:
            if chrom1[tempindex] != "-":
                tempcount += 1
            tempindex += 1

        offset = tempindex - j
        prevoffset = offset

        ChromSlice: str = chrom1[j:j + offset]
        SubfamSlice: str = subfam1[j:j + offset]

        alignscore = CalcScore(SubfamSlice, ChromSlice, "", "")
        AlignHash[i, index] = alignscore

        if cols < index:
            cols = index + 1

        index += 1

        # TODO: Make sure these bounds are right since Python indexing is different
        while j + offset + 1 < len(chrom1):
            tempindex = j
            tempcount = 0

            while tempcount < chunksize:
                if chrom1[tempindex + 1] != "-":
                    tempcount += 1
                tempindex += 1

            offset = tempindex - j

            if chrom1[j + 1] != "-":
                if chrom1[j + 1] != "." and subfam1[j + 1] != ".":
                    if prevoffset != offset:
                        ChromSlice = chrom1[j+1:j+offset+1]
                        SubfamSlice = subfam1[j+1:j + offset + 1]
                        alignscore = CalcScore(SubfamSlice, ChromSlice, subfam1[j], chrom1[j])
                    else:
                        if subfam1[j] == "-":
                            if subfam1[j - 1] == "-":
                                alignscore = alignscore - gapExt
                            else:
                                alignscore = alignscore - gapInit
                        else:
                            alignscore = alignscore - SubMatrix[CharPos[subfam1[j]]*subMatrixCols+CharPos[chrom1[j]]]

                        if subfam1[j + offset] == "-":
                            if subfam1[j + offset - 1] == "-":
                                alignscore = alignscore + gapExt
                            else:
                                alignscore = alignscore + gapInit
                        elif subfam1[j + offset] == "." or chrom1[j + offset] == ".":
                            alignscore = alignscore
                        else:
                            alignscore = alignscore + SubMatrix[CharPos[subfam1[j + offset]] * subMatrixCols + CharPos[chrom1[j + offset]]];

                    if alignscore <= 0:
                        AlignHash[i, index] = 1
                    else:
                        AlignHash[i, index] = alignscore

                if cols < index:
                    cols = index + 1
                index += 1

            j += 1
            prevoffset = offset

        i += 1


def FillConsensusPosMatrix(consensus: Dict[Tuple[int, int], int],
                           subfams: List[str],
                           chroms: List[str],
                           consensusstart: List[int],
                           consensusstop: List[int]):
    i: int = 1
    while i < rows:
        SubfamArray: List[str] = subfams[i].split()
        ChromArray: List[str] = chroms[i].split()

        consensuspos: int = 0
        if Strands[i] == "+":
            consensuspos = consensusstart[i] - 1
            matrixpos: int = 0
            j: int = 0
            while j < len(SubfamArray):
                if SubfamArray[j] != ".":
                    if SubfamArray[j] != "-":
                        consensuspos += 1

                    # TODO: What does the -> in the perl code do?
                    consensus[i, matrixpos] = consensuspos

                if ChromArray[j] != "-":
                    matrixpos += 1
        else:
            consensuspos = consensusstart[i] + 1
            matrixpos: int = 0
            j: int = 0
            while j < len(SubfamArray):
                if SubfamArray[j] != ".":
                    if SubfamArray[j] != "-":
                        consensuspos -= 1
                    consensus[i, matrixpos] = consensuspos

                if ChromArray[j] != "-":
                    matrixpos += 1

        # FIXME: Delete later, just for testing (see perl code, line 985)
        if consensusstop[i] != consensuspos:
            stderr.write("\n\nERROR - consensus seq positions not correct\n\n")


def FillConfScoreMatrix(alignhash: Dict[Tuple[int, int], int], confhash: Dict[Tuple[int, int], float]):
    i: int = 0
    while i < len(Columns) - chunksize + 1:
        col: int = Columns[i]
        temp: List[int] = []
        row: int = 0
        while row < rows:
            if alignhash[row, col]:
                temp.append(alignhash[row, col])
            else:
                temp.append(0)

        confidenceTemp = ConfidenceCM(lamb, temp).split(" ")
        row: int = 0
        while row < rows:
            if confidenceTemp[row] != 0:
                confhash[row, col] = confidenceTemp[row]


def FillProbMatrix(probhash: Dict[Tuple[str, int], float], supporthash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], int]):
    j: int = 1
    col: int = 1
    while col < len(Columns):
        if col in Columns:
            j = Columns[col]
        else:
            j += 1

        for i in ActiveCellsCollapse[j]:
            max: float = -inf
            maxindex: str = ''
            supportlog: float = log(supporthash[i, j])

            for row in ActiveCellsCollapse[Columns[col - 1]]:
                score: float = -1
                there: bool = False

                if col in Columns:
                    score = supportlog + probhash[row, col - 1]
                    there = True
                else:
                    score = supportlog + probhash[row, j - 1]
                    there = True

                if there:
                    if row == i:
                        score = score + sameProbLog
                    else:
                        score = score + changeProbLog

                    if score > max:
                        max = score
                        maxindex = row

            probhash[i, j] = max
            originhash[i, j] = maxindex

        col += 1



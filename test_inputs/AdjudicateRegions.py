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

    stdout.write("skip\t")
    j: int = 0
    while j < column:
        if ("skip", j) in Hash:
            stdout.write(str(Hash["skip", j]))
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
                if (k, j) in Hash:
                    stdout.write(str(Hash[k, j]))
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
skipAlignScore: int = 30 #FIXME - still need to decide what this number is, skip state doesn't work in seqs_fullAlu.align unless skipAlignScore = 120
startall: int = 0  # Reassigned later
stopall: int = 0  # Reassigned later
ID: int = 1111

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
line = re.sub(r"\s+$", "", line)
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

RemoveStarts: List[int] = []
RemoveStops: List[int] = []

Changes: List[str] = []
ChangesPos: List[int] = []

IDs: List[int] = [];
Changes_orig: List[int] = [];
ChangesPos_orig: List[int] = [];
Columns_orig: List[int] = [];

SupportHashCollapse: Dict[Tuple[str, int], int] = {}
ActiveCellsCollapse: Dict[int, List[str]] = {}
SubFamsCollapse: Dict[str, int] = {}
ConsensusHashCollapse: Dict[Tuple[str, int], int] = {}
StrandHashCollapse: Dict[Tuple[str, int], str] = {}



numseqs: int = 0
with open(infile) as _infile:
    alignments = load_alignments(_infile)
    for alignment in alignments:
        numseqs += 1

        SubFams.append(alignment.subfamily)
        Scores.append(alignment.score)
        Strands.append(alignment.strand)
        Starts.append(alignment.start)
        Stops.append(alignment.stop)
        ConsensusStarts.append(alignment.consensus_start)
        ConsensusStops.append(alignment.consensus_stop)
        SubfamSeqs.append(alignment.subfamily_sequence)
        ChromSeqs.append(alignment.sequence)

#if there is only one subfam in the alignment file, no need to run anything because we know
#that subfam is what's there 
#FIXME - add stop of the alignment         
if numseqs == 1:
	stdout.write("1\t")
	stdout.write(SubFams[1])
	stdout.write("\n")
	exit()

changeProbLog = log(changeProb / (numseqs - 1))


rows: int = len(SubFams)
cols: int = 0


def Edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    minStart: int = starts[1]
    maxStop: int = stops[1]

    for i in range(1, len(starts)):
        if starts[i] < minStart:
            minStart = starts[i]
        if stops[i] > maxStop:
            maxStop = stops[i]

    return minStart, maxStop


def padSeqs(start, stop, subfamseq, chromseq):
    global startall, stopall

    startall, stopall = Edges(start, stop)

    for i in range(1, len(SubfamSeqs)):
        leftpad = start[i] - startall
        rightpad = stopall - stop[i]

        chromseq[i] = ("." * leftpad) + f"{chromseq[i]}" + ("." * rightpad)
        subfamseq[i] = ("." * leftpad) + f"{subfamseq[i]}" + ("." * rightpad)



padSeqs(Starts, Stops, SubfamSeqs, ChromSeqs)

	
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
        
    for j in range(1, len(seq1)):
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

    return chunkscore

##FIXME - come back to this one later once I have the results printed - easier for testing 
def FillAlignScoreMatrix(subfams: List[str], chroms: List[str]):
    global cols

    for i in range(1, len(ChromSeqs)):
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
        
        index+=1

        # TODO: Make sure these bounds are right since Python indexing is different
        while j + offset < len(chrom1):
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

                index += 1

            j += 1
            prevoffset = offset
            
        cols = index
        i += 1


FillAlignScoreMatrix(SubfamSeqs, ChromSeqs)

def FillConsensusPosMatrix(consensus: Dict[Tuple[int, int], int],
                           subfams: List[str],
                           chroms: List[str],
                           consensusstart: List[int],
                           consensusstop: List[int]):
                           
    #0s for consensus pos of skip state
    for j in range(cols+chunksize-1):
    	consensus[0, j] = 0
	
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


#instead of an array like in the perl code, this is a dict
#keys are the perl array indicies and values are the matrix columns
#needs to be dict so I can make some of the matrix columns undef - in python if you undef
# from a list it splices out the indices
Columns = []
j: int = 0
for j in range(cols):
	empty = 1;
	for i in range(rows):
		if (i, j) in AlignHash:
			empty = 0
			i = rows
 	 	
	if not empty:
		Columns.append(j)
 
	AlignHash[0, j] = skipAlignScore


#adding the last chunksize-1 columns that are not in the AlignScoreMatrix to Columns
for i in range(chunksize-1):
	j+=1
	Columns.append(j)
	

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
    for i in range(len(Columns)):
        if i >= len(Columns) - chunksize + 1:
            break
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
        for col in range(len(Columns)):
            if col >= len(Columns) - chunksize + 1:
                break
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
# PrintMatrixHash(cols, ConfHash)

#collapses matrices 
for col in range(len(Columns)):
	j: int = Columns[col]
	DupMaxCon: Dict[str, float] = {}
	DupMaxSup: Dict[str, float] = {}
	
	activecols: List[str] = []
	ActiveCellsCollapse[j] = activecols
    
	for i in range(rows):
		if (i,j) in SupportHash and (i,j) in ConsensusHash:
			if (SubFams[i]) in DupMaxCon:
				if SupportHash[i,j] > DupMaxCon[SubFams[i]]:
					DupMaxCon[SubFams[i]] = SupportHash[i,j]
					ConsensusHashCollapse[SubFams[i],j] = ConsensusHash[i,j]
					StrandHashCollapse[SubFams[i], j] = Strands[i]
			else:
				DupMaxCon[SubFams[i]] = SupportHash[i,j]
				ConsensusHashCollapse[SubFams[i],j] = ConsensusHash[i,j]
				StrandHashCollapse[SubFams[i], j] = Strands[i]
		
		if (i,j) in SupportHash:
			if (SubFams[i]) in DupMaxSup:
				if SupportHash[i,j] > DupMaxSup[SubFams[i]]:
					DupMaxSup[SubFams[i]] = SupportHash[i,j]
					SupportHashCollapse[SubFams[i],j] = SupportHash[i,j]
			else:
				DupMaxSup[SubFams[i]] = SupportHash[i,j]
				SupportHashCollapse[SubFams[i],j] = SupportHash[i,j]
				ActiveCellsCollapse[j].append(SubFams[i])

		

for i in range(rows):
    SubFamsCollapse[SubFams[i]] = 0

rows = len(SubFamsCollapse)

# line 381
for k in SubFamsCollapse:
    ProbHash[k, 0] = 0


def FillProbMatrix(probhash: Dict[Tuple[str, int], float], supporthash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str]):
	for j in range(1, len(Columns)):
	 			
		for i in ActiveCellsCollapse[Columns[j]]:
			max: float = -inf
			maxindex: str = ''
			supportlog: float = log(supporthash[i, Columns[j]])
			
			#loop through all the subfams in the previous column
			for row in ActiveCellsCollapse[Columns[j-1]]:
				score: float = supportlog + probhash[row, Columns[j-1]]
				
				if row == i:
# 					score = score + sameProbLog
					# because rows are collapsed, if the consensus seqs are contiguous - treat as if they are not the same row and get the jump penalty
					if (row, Columns[j-1]) in ConsensusHashCollapse and (i, Columns[j]) in ConsensusHashCollapse:
						if ConsensusHashCollapse[row, Columns[j-1]] > ConsensusHashCollapse[i, Columns[j]] + 50:
							score = score + changeProbLog
						else:
							score = score + sameProbLog
					else:
						score = score + sameProbLog
						
				else:
					score = score + changeProbLog
					
				if score > max:
					max = score
					maxindex = row
					
			probhash[i, Columns[j]] = max
			originhash[i, Columns[j]] = maxindex
			

# PrintMatrixHashCollapse(cols, SupportHashCollapse)
# PrintMatrixHashCollapse(cols, ConsensusHashCollapse)
			
#FIXME - these numbers are slightly different than the perl ones... I think it's just
#rounding error?
FillProbMatrix(ProbHash, SupportHashCollapse, OriginHash)

IDs = [0] * cols

def GetPath(probhash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str], subfams: List[str]) -> List[str]:
    global ID
    maxxx: float = -inf
    maxindex: str = ''
    
    for i1 in ActiveCellsCollapse[cols - 1]:
        if maxxx < probhash[i1, cols - 1]:
            maxxx = probhash[i1, cols - 1]
            maxindex = i1
    
    prev: str = originhash[maxindex, cols - 1]
    
    ChangesPos.append(len(Columns))
    
    j: int = len(Columns)-1
    while j > 1:
    	if (prev, Columns[j-1]) in originhash:
    	
    		IDs[Columns[j-1]] = ID
    		
    		if prev != originhash[prev, Columns[j-1]]:
    			ID+=1234
    			ChangesPos.append(j-1)
    			Changes.append(prev)
    			
    		prev = originhash[prev, Columns[j-1]]
    		
    	j-=1
    	
    	
    IDs[Columns[j-1]] = ID	
    ChangesPos.append(j-1)
    Changes.append(prev)
    
    Changes.reverse()
    ChangesPos.reverse()
    
    #changes ID for next round of stitching, so when starts stitching will have unique ID
    ID+=1234


GetPath(ProbHash, OriginHash, SubFams)

#keep the original info for printing at the end 
Changes_orig = Changes.copy()
ChangesPos_orig = ChangesPos.copy()
Columns_orig = Columns.copy()


def PrintChanges(changes, changespos):
    i: int = 0
    while i < len(changes):
        stdout.write(str(Columns[changespos[i]]))
        stdout.write("\t")
        stdout.write(f"{changes[i]}\n")
        i=i+1


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
    
    
def NodeConfidence(nodeconfidence: Dict[Tuple[str, int], float], subfamseqs: List[str], chromseqs, changespos: List[int]):
    nodeconfidence_temp: List[float] = [0 for _ in range(len(SubFams) * numnodes)]
    
    #calculated node confidence for first node - doesn't look back at prev char bc there isn't one
    #starts at 1 bc don't need to calc for skip state
    for j in range(1, len(SubFams)):
        b: int = Columns[changespos[0]]
        e: int = Columns[changespos[1]]
        subfam: str = subfamseqs[j][b:e]
        chrom: str = chromseqs[j][b:e]
        alignscore: float = CalcScore(subfam, chrom, '', '')
        nodeconfidence_temp[j * numnodes + 0] = alignscore	 
	 
	#does rest of nodes - looks back at prev char incase of gap ext
    for i in range(1, numnodes):
    	for j in range(1, len(SubFams)):
            b: int = Columns[changespos[i]]
            e: int = Columns[changespos[i + 1]-1]
            subfam: str = subfamseqs[j][b:e]
            chrom: str = chromseqs[j][b:e]
            lastpreva: str = subfamseqs[j][b-1]#subfamseqs[j][Columns[changespos[i + 1] - 1]]
            lastprevb: str = chromseqs[j][b-1]#chromseqs[j][changespos[i + 1] - 1]
            alignscore: float = CalcScore(subfam, chrom, lastpreva, lastprevb)
            nodeconfidence_temp[j * numnodes + i] = alignscore
    		
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

def PrintPathGraph():
	stdout.write(" ")
	for i in range(numnodes):
		stdout.write(f"{Changes[i]} ")
	stdout.write("\n")
	
	for i in range(numnodes):
		for j in range(numnodes):
			stdout.write(f"{pathGraph[i*numnodes+j]}\t")
		stdout.write(f"{Changes[i]}\n")
	stdout.write("\n")

def PrintNodeConfidence():
	for i in range(numnodes):
		stdout.write(f"{Changes[i]} ")
	stdout.write("\n")

	for subfam in SubFamsCollapse:
		stdout.write(f"{subfam} ")
		for j in range(numnodes):
			if (subfam,j) in NodeConfidenceDict:
				stdout.write(f"{NodeConfidenceDict[subfam,j]} ")
			else:
				stdout.write(f"-inf ")
		stdout.write("\n")



#FIXME - in Seqs_noStitching.align AluJr being stitched, but shouldn't be - think it's an 
# error in consensusHashCollapse the current method for testing isn't working - the 
# consensusHashCollapse is correct, but the switches in probmatrix are happening later 
# than actual so it's looking at a diff consensus seqs postions when testing if source 
# is before sink
##What if we check 50nt around the sinkstart and sourcestop and if any of the pos satisfies 
##the condition stitch them?
def FillPathGraph(pathgraph: List[int]):
    for i in range(numnodes * numnodes):
        pathgraph.append(0)

    for i in range(numnodes - 1):
        pathgraph[i * numnodes + i + 1] = 1

    for j in range(numnodes):
        sinkSubfam: str = str(Changes[j])
        sinkSubfamStart: int = ConsensusHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]
        sinkStrand: str = StrandHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]

        for i in range(j - 1):
            for sourceSubfam in SubFamsCollapse:
            	sourceSubfam = str(sourceSubfam)
            	sourceConf = NodeConfidenceDict[sourceSubfam, i]
            	
            	if (sourceSubfam, Columns[ChangesPos[i + 1]-1]) in ConsensusHashCollapse:
            		sourceSubfamStop = ConsensusHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]-1]]
            		sourceStrand = StrandHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]-1]]
            		
            		if sinkStrand == '+' and sinkStrand == sourceStrand:
            			if (sinkSubfam == sourceSubfam) and (sourceConf >= 0.3):            				
            				if sourceSubfamStop <= sinkSubfamStart + 50:
            					pathgraph[i * numnodes + j] = 1
            					
            		elif sinkStrand == '-' and sinkStrand == sourceStrand:
            			if sinkSubfam == sourceSubfam and sourceConf >= 0.3:
            				if sourceSubfamStop >= sinkSubfamStart + 50:
            					pathgraph[i * numnodes + j] = 1
    
#     PrintPathGraph()
#     PrintNodeConfidence()


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
		removestops.append(len(Columns)-1)
		RemoveNodes[numnodes - 1] = True
		
	#when removing from the end, have to update cols because don't want to do to the end of the matrix anymore 
	i: int = numnodes - 1
	while RemoveNodes[i]:
		cols = Columns[changespos[i] - 1]  #FIXME - do I index cols here or not?
		i -= 1



numnodes: int = 0
NodeConfidenceDict: Dict[Tuple[str, int], float] = {}
pathGraph: List[int] = []
total: int = 0
loop: int = 1


def PrintResults():
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(Changes_orig)):
		if str(Changes_orig[i]) != 'skip':
			stdout.write(str(Columns_orig[ChangesPos_orig[i]]))
			stdout.write("\t")
			stdout.write(str(Columns_orig[ChangesPos_orig[i+1]-1]))
			stdout.write("\t")
			stdout.write(str(IDs[Columns_orig[ChangesPos_orig[i]]]))
			stdout.write("\t")
			stdout.write(str(Changes_orig[i]))
			stdout.write("\n")
	

def PrintResultsSequence():
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(Changes_orig)):
		if str(Changes_orig[i]) != 'skip':
			stdout.write(str(Columns_orig[ChangesPos_orig[i]]+startall))
			stdout.write("\t")
			stdout.write(str(Columns_orig[ChangesPos_orig[i+1]-1]+startall))
			stdout.write("\t")
			stdout.write(str(IDs[Columns_orig[ChangesPos_orig[i]]]))
			stdout.write("\t")
			stdout.write(str(Changes_orig[i]))
			stdout.write("\n")

count: int = 0
while (True):
    count += 1
    numnodes = len(Changes)
    if (numnodes <= 2):
        break # `last;` in perl
	
    NodeConfidenceDict.clear()
    
    NodeConfidence(NodeConfidenceDict, SubfamSeqs, ChromSeqs, ChangesPos)
    
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
        break # `last;` in perl

    RemoveStarts.clear()
    RemoveStops.clear()

    ExtractNodes(RemoveStarts, RemoveStops, ChangesPos, pathGraph, numnodes)
    
    total: int = 0
    for i in range(len(RemoveStops)):
    	del Columns[RemoveStarts[i]-total:RemoveStops[i]-total]
    	total += (RemoveStops[i] - RemoveStarts[i])
        
    FillProbMatrix(ProbHash, SupportHashCollapse, OriginHash)

    Changes.clear()
    ChangesPos.clear()
        
    GetPath(ProbHash, OriginHash, SubFams)


# PrintResultsSequence()
PrintResults()

        


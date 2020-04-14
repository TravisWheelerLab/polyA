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

	#Subfams_collapse is not ordered, so prints the skip state first 
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
chunksize: int = 31
sameProbLog: float = log(1 - (10 ** -45))  #FIXME - this just becomes 0.0... need to be more precise
changeProb: float = 10 ** -45
changeProbLog: float = 0.0  # Reassigned later
changeProbSkip: float = 0.0 # Reassigned later
sameProbSkip: float = 0.0
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

gapInit = int(opts["--gapInit"]) if "--gapInit" in opts else gapInit
gapExt = int(opts["--gapExt"]) if "--gapExt" in opts else gapExt
lamb = float(opts["--lambda"]) if "--lambda" in opts else lamb
chunksize = int(opts["--segmentsize"]) if "--segmentsize" in opts else chunksize
changeProb = float(opts["--changeprob"]) if "--changeprob" in opts else changeProb
help = "--help" in opts
printt = "--printmatrices" in opts
printMatrixPos = "--matrixPos" in opts

if help:
    print(helpMessage)
    exit(0)

#input is alignment file of hits region and substitution matrix 
infile: str = args[0]
infile_matrix: str = args[1]

# Other open was moved down to where we load the alignments file
with open(infile_matrix) as _infile_matrix:
    in_matrix: List[str] = _infile_matrix.readlines()

#FIXME - want to add option to run easel in here
# We require a --lambda or provide a default so there's no need to run easel

#creates a dict that associates character from score matrix file with position in score matrix
#alignment char as key and pos in array as value   
CharPos: Dict[str, int] = {}
line = in_matrix[0]
line = re.sub(r"^\s+", "", line)
line = re.sub(r"\s+$", "", line)
chars = re.split(r"\s+", line)
for i in range(len(chars)):
    CharPos[chars[i]] = i
subMatrixCols: int = len(chars)


#reads in the score matrix from file and stores in 2D array Score matrix 
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


#opens alignment file and stores all Subfams, Scores, Starts, Stops, subfams seqs and Chrom seqs in arrays 
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
#2 because of the skip state 
if numseqs == 2:
	if printMatrixPos:
		stdout.write("start\tstop\tID\tname\n")
		stdout.write("----------------------------------------\n")
		stdout.write(f"{0}\t{Stops[1]-Starts[1]}\t1111\t{SubFams[1]}\n")
	else:
		stdout.write("start\tstop\tID\tname\n")
		stdout.write("----------------------------------------\n")
		stdout.write(f"{Starts[1]}\t{Stops[1]}\t1111\t{SubFams[1]}\n")
	exit()

changeProbLog = log(changeProb / (numseqs - 1))
changeProbSkip = changeProbLog / 2;
sameProbSkip = changeProbLog / 20; # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump

#precomputes global vars rows and cols in matrices 
rows: int = len(SubFams)
cols: int = 0 #assign cols in FillAlignScoreMatrix


#find the min start and max stop for the whole region
def Edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    minStart: int = starts[1]
    maxStop: int = stops[1]

    for i in range(1, len(starts)):
        if starts[i] < minStart:
            minStart = starts[i]
        if stops[i] > maxStop:
            maxStop = stops[i]

    return minStart, maxStop


#pads sequences with '.' and makes sparse matrix - allows regions where seqs do not all have same 
#start and stop positions
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

	#deals with the first character of a segment being a gap character - have to look at last
	#segment to see if this is a gap init or ext
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


#calculates alignment score (according to crossmatch scoring) for every segment and all seqs
#fills align score matrix - pos 0 in the array is the score of the segment that starts at pos 0,
#pos 1 in array is score of segment that starts at pos 1, etc
#computes score for the first segment that does not start with a '.' and from there keeps the 
# base score and adds new char and subtracts old char - if a new gap is introduced, calls CalcScore()
#instead of adding on base score 
##FIXME - come back to this one later once I have the results printed - easier for testing 
def FillAlignScoreMatrix(subfams: List[str], chroms: List[str]):
    global cols
    
    index: int = 0
    
    #chunks can't start on gaps and gaps don't count when getting to the 30 bps
	
    for i in range(1, len(ChromSeqs)):
        subfam1: str = subfams[i]
        chrom1: str = chroms[i]
                
        #grab first chunk of 30 and calculate the raw score 
		
		#starts at the first non '.' char, but offsets it in the matrix based on where
		#the alignments start in the seq - ex: if first alignment in the seq starts at 10,
		#will offset by 10

        j: int = Starts[i] - startall
        index = j+15  #index is the col we are in the align score matrix, $j is the place in @subfam1 and @chrom1

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

		#grabs first chunk - here $ j = pos of first non '.' char
        ChromSlice: str = chrom1[j:j + offset]
        SubfamSlice: str = subfam1[j:j + offset]

		#calculates score for first chunk and puts score in alignhash
        alignscore = CalcScore(SubfamSlice, ChromSlice, "", "")
        AlignHash[i, index] = alignscore # already to scale so don't need to * 31 and / 31
        
        print(index)
        print(SubfamSlice)
        print(ChromSlice)
        
        index+=1
        
        numnucls: int = chunksize  #how many nucls contributed to align score 

        # TODO: Make sure these bounds are right since Python indexing is different
		#move to next chunk by adding next chars score and subtracting prev chars score
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
                    if prevoffset != offset:  #there is a new gap, or a gap was removed from beginning 
                        ChromSlice = chrom1[j+1:j+offset+1]
                        SubfamSlice = subfam1[j+1:j + offset + 1]
                        alignscore = CalcScore(SubfamSlice, ChromSlice, subfam1[j], chrom1[j])
                        
                        tempcount2: int = 0
                        for nuc in ChromSlice:
                        	if nuc != '-' and nuc != '.':
                        		tempcount2 += 1
                        numnucls = tempcount2 #resetting numnucls to 31
                        
                    else:
 						#alignscore from previous segment - prev chars score + next chars score 
						#subtracting prev  chars score - tests if its a gap in the subfam as well

                        if subfam1[j] == "-":
                        	numnucls -= 1
                        	if subfam1[j - 1] == "-":
                        		alignscore = alignscore - gapExt
                        	else:
                        		alignscore = alignscore - gapInit
                        else:
                            alignscore = alignscore - SubMatrix[CharPos[subfam1[j]]*subMatrixCols+CharPos[chrom1[j]]]
                            numnucls -= 1

						#adding next chars score - tests if its a gap in the subfam as well
                        if subfam1[j + offset] == "-":
                        	numnucls += 1
                        	if subfam1[j + offset - 1] == "-":
                        		alignscore = alignscore + gapExt
                        	else:
                        		alignscore = alignscore + gapInit
                        elif subfam1[j + offset - 15] == "." or chrom1[j + offset - 15] == ".":
                            alignscore = -inf
                        elif subfam1[j + offset] == "." or chrom1[j + offset] == ".":
                        	alignscore = alignscore
                        else:
                            alignscore = alignscore + SubMatrix[CharPos[subfam1[j + offset]] * subMatrixCols + CharPos[chrom1[j + offset]]]
                            numnucls += 1
                    
                    if alignscore <= 0:
                        AlignHash[i, index] = 1
                    else:
                        AlignHash[i, index] = int(alignscore / numnucls * chunksize)
                        
                    if alignscore == -inf:
                    	del AlignHash[i, index]

                index += 1

            j += 1
            prevoffset = offset
            
        i += 1
    cols = index


FillAlignScoreMatrix(SubfamSeqs, ChromSeqs)
print(cols)


#FIXME - add first 15 and last 15 into matrix
cols = cols + 30
# PrintMatrixHash(cols, AlignHash)

exit()

#fills parallel array to the Align Matrix that holds the consensus position for each 
# subfam at that position in the alignment
def FillConsensusPosMatrix(consensus: Dict[Tuple[int, int], int],
                           subfams: List[str],
                           chroms: List[str],
                           consensusstart: List[int],
                           consensusstop: List[int]):
                           
    #0s for consensus pos of skip state
    for j in range(cols+chunksize-1):
    	consensus[0, j] = 0
	
	#start at 1 to skip 'skip state'
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
                	#consensus pos only advances when there is not a gap in the subfam seq
                    if SubfamArray[j] != "-":
                        consensuspos += 1

                   #put consensus pos corresponding to pos in matrix in hash
                    consensus[i, matrixpos] = consensuspos

				#matrix position only advances when there is not a gap in the chrom seq
                if ChromArray[j] != "-":
                    matrixpos += 1
                j += 1
        else:  #reverse strand 
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
        i += 1


FillConsensusPosMatrix(ConsensusHash, SubfamSeqs, ChromSeqs, ConsensusStarts, ConsensusStops)


#puts all columns that are not empty into @Columns, so when I loop through hash I can use the 
#vals in @Columns - this will skip over empty columns
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
	
	#assigns skip states an alignment score 
	AlignHash[0, j] = skipAlignScore


#adding the last chunksize-1 columns that are not in the AlignScoreMatrix to Columns
for i in range(chunksize-1):
	j+=1
	Columns.append(j)
	
#send in an array of scores for a segment - output an array of confidence values for the segment
def ConfidenceCM(lamb: float, region: List[float]) -> str:
    confidenceString: str = ""

	#loops through the array once to get the sum of 2^every_hit_score in region 
	#converts the score to account for lambda before summing 
    ScoreTotal: int = 0
    for Score in region:
        if Score > 0:
            convertedScore = Score * lamb
            ScoreTotal += 2 ** convertedScore

	#once region score is summed, loop back through the region array and calculate
	#confidence for each hit
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

#reassign global var col to account for last chunksize-1 cols added to the next matrices 
# cols = cols + chunksize - 1

# Fills support score matrix using values in conf matrix
#score for subfam x at position i is sum of all confidences for subfam x for all segments that 
#overlap position i - divided by number of segments
def FillSupportMatrix(supporthash: Dict[Tuple[int, int], float], alignhash: Dict[Tuple[int, int], int], confhash: Dict[Tuple[int, int], float]):
    
    # i = subfam, j = position
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

		#Columns only goes until the end of align score matrix so it is $chunksize smaller than
		#we need here - do this last part of the loop so fill in the rest of the support matrix	
#         j: int = cols - chunksize
#         while j <= cols:
#             num: int = j
#             summ: float = 0
#             numsegments: int = 0
#             while num >= 0 and num >= j - chunksize + 1:
#                 if (i, num) in confhash:
#                     summ = summ + confhash[i, num]
#                     numsegments += 1
#                 num -= 1
#             if numsegments > 0:
#                 supporthash[i, j] = summ / numsegments
#             j += 1

        i += 1


FillSupportMatrix(SupportHash, AlignHash, ConfHash)


#collapses matrices 
#collapse and combine rows that are the same subfam - just sum their support 
#new support dict has key = subfamname.col 
#also creates a bookkeeping dict that has all the cols as keys and their values 
#are arrays that hold all the active subfams in that col - used so that don't have 
#to loop through all the $i's just to see if a column exists 
for col in range(len(Columns)):
	j: int = Columns[col]
	DupMaxCon: Dict[str, float] = {}
	DupMaxSup: Dict[str, float] = {}
	
	activecols: List[str] = []
	ActiveCellsCollapse[j] = activecols
	
	#sum the support score for rows that are collapsed together
	#find max support score for collapsed rows and use the consensus from that row
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

#update global var rows after collapse 
rows = len(SubFamsCollapse)

#fill first col of probhash with 0s
for k in SubFamsCollapse:
    ProbHash[k, 0] = 0



##FIXME - I think we can change all log values in dp matrix to int and speed up dp
## and all the same and change penalties can become ints 

#fills prob score matrix from support matrix hash
#also fills origin matrix
#skips first col because want it to have a prob of 0
# calculates probability score for a cell in the confidence score matrix
# 	look at all i's in j-1
#	mult by confidence in current cell
# 	if comes from same i, mult by higher prob
# 	else - mult by lower prob /(numseqs-1) -> so sum of all probs == 1
# 	return max
#done in log space
def FillProbMatrix(probhash: Dict[Tuple[str, int], float], supporthash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str]):
	for j in range(1, len(Columns)):
	 			
		for i in ActiveCellsCollapse[Columns[j]]:
			max: float = -inf
			maxindex: str = ''
			supportlog: float = log(supporthash[i, Columns[j]])
			
			#loop through all the subfams in the previous column
			for row in ActiveCellsCollapse[Columns[j-1]]:
				score: float = supportlog + probhash[row, Columns[j-1]]
				prob: int = 0
				
				if row == i:  #staying in same row
					
					if row == 'skip':  #staying in skip
						prob = sameProbSkip
						
					prob = sameProbLog
					
					# print(StrandHashCollapse[row, Columns[j-1]])
# 					print(StrandHashCollapse[i, Columns[j]])

					# because rows are collapsed, if the consensus seqs are contiguous - treat as if they are not the same row and get the jump penalty
					if StrandHashCollapse[row, Columns[j-1]] != StrandHashCollapse[i, Columns[j]]:
						prob = changeProbLog
					else: #if on same strand
						if StrandHashCollapse[row, Columns[j-1]] == '+':
							if ConsensusHashCollapse[row, Columns[j-1]] > ConsensusHashCollapse[i, Columns[j]] + 50:
								prob = changeProbLog
						if StrandHashCollapse[row, Columns[j-1]] == '-':
							if ConsensusHashCollapse[row, Columns[j-1]] + 50 < ConsensusHashCollapse[i, Columns[j]]:
								prob = changeProbLog
						
				else:  #jumping rows
					prob = changeProbLog  
					if row == 'skip' or i == 'skip':  #jumping in or out of skip
						prob = changeProbSkip
					
				score = score + prob
					
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

#using origin matrix, back traces through the 2D array to get the subfam path
#finds where the path switches to a different row and populates @Changes and @ChangesPos
#reverses @Changes and @ChangesPos because it's a backtrace so they are initially backwards
#jumps over removed columns when necessary
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
    
    #already added the last col, but this adds the one before $col so still start at last col
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

#keep the original annotation for reporting results
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
 
#FIXME - if goes into skip state for just one position, this will have an error .. also when
#stitching want to ignore skip states
#fills node confidence matrix 
#first fills matrix with node alignment scores, then reuses matrix for confidence scores    
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
    for i in range(1, numnodes-1):
    	for j in range(1, len(SubFams)):
            b: int = Columns[changespos[i]]
            e: int = Columns[changespos[i + 1]]
            subfam: str = subfamseqs[j][b:e]
            chrom: str = chromseqs[j][b:e]
            lastpreva: str = subfamseqs[j][b-1]#subfamseqs[j][Columns[changespos[i + 1] - 1]]
            lastprevb: str = chromseqs[j][b-1]#chromseqs[j][changespos[i + 1] - 1]
            alignscore: float = CalcScore(subfam, chrom, lastpreva, lastprevb)
            nodeconfidence_temp[j * numnodes + i] = alignscore
    
	#does last node
    for j in range(1, len(SubFams)):
        b: int = Columns[changespos[-2]]
        e: int = Columns[changespos[-1]-1]
        subfam: str = subfamseqs[j][b:e]
        chrom: str = chromseqs[j][b:e]
        lastpreva: str = subfamseqs[j][b-1]
        lastprevb: str = chromseqs[j][b-1]
        alignscore: float = CalcScore(subfam, chrom, lastpreva, lastprevb)
        nodeconfidence_temp[j * numnodes + numnodes-1] = alignscore
 
    #reuse same matrix and compute confidence scores for the nodes	
    for j in range(numnodes):
    	temp: List[float] = []
    	for i in range(len(SubFams)):
    		temp.append(nodeconfidence_temp[i * numnodes + j])
    		
    	confidenceTemp: List[float] = [float(x) for x in ConfidenceCM(lamb, temp).split(' ')]
    	
    	for i in range(len(SubFams)):
    		nodeconfidence_temp[i * numnodes + j] = confidenceTemp[i]

	#collapse nodeconfidence down same way supportmatrix is collapsed - all seqs of 
	#the same subfam are put in the same row
	#not a sparse hash - holds the 0s, but I think this is okay because it won't ever
	#be a very large matrix, and this way we don't have to test if anything exists 	
	# my %nodeconfidence_collapse;
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



def FillPathGraph(pathgraph: List[int]):
    for i in range(numnodes * numnodes):
        pathgraph.append(0)

	# filling beginning path graph with straight line through the nodes
    for i in range(numnodes - 1):
        pathgraph[i * numnodes + i + 1] = 1

    for j in range(numnodes):
        sinkSubfam: str = str(Changes[j])
        sinkSubfamStart: int = ConsensusHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]
        sinkStrand: str = StrandHashCollapse[sinkSubfam, Columns[ChangesPos[j]]]

		#looks at all the preceding nodes, except the one directly before it (source nodes)
        for i in range(j - 1):
        	#look at all the subfams in each node
            for sourceSubfam in SubFamsCollapse:
            	sourceSubfam = str(sourceSubfam)
            	sourceConf = NodeConfidenceDict[sourceSubfam, i]
            	
            	if (sourceSubfam, Columns[ChangesPos[i + 1]-1]) in ConsensusHashCollapse:
            		sourceSubfamStop = ConsensusHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]-1]]
            		sourceStrand = StrandHashCollapse[sourceSubfam, Columns[ChangesPos[i + 1]-1]]
  
  					# adds in edge if the subfam of the sink is at the source node and if it's 
					# confidence >= 1%, and if the source is before the sink in the consensus sequence           		
            		if sinkStrand == '+' and sinkStrand == sourceStrand:
            			if (sinkSubfam == sourceSubfam) and (sourceConf >= 0.3):
            				#FIXME- not sure what this overlap should be .. just allowed 50 for now            				
            				if sourceSubfamStop <= sinkSubfamStart + 50:
            					pathgraph[i * numnodes + j] = 1
            					
            		elif sinkStrand == '-' and sinkStrand == sourceStrand:
            			if sinkSubfam == sourceSubfam and sourceConf >= 0.3:
            				if sourceSubfamStop >= sinkSubfamStart + 50:
            					pathgraph[i * numnodes + j] = 1
    
#     PrintPathGraph()
#     PrintNodeConfidence()


#finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds them to
# @RemoveStarts and @RemoveStops so they can be extracted from the alignment 
def ExtractNodes(removestarts, removestops, changespos, pathgraph, numnodes: int):
	global cols
	
	#boolean for which nodes will be removed
	RemoveNodes: List[bool] = [False for _ in range(numnodes)]
	
	#extracting nodes that only have one incoming and one outgoing edge
	NumEdgesIn: List[int] = [0 for _ in range(numnodes)]
	NumEdgesOut: List[int] = [0 for _ in range(numnodes)]
	
	for i in range(numnodes):
		for j in range(numnodes):
			#$j - incoming
			NumEdgesIn[j] += pathgraph[i * numnodes + j]
			#$j - outgoing
			NumEdgesOut[i] += pathgraph[i * numnodes + j]
			
	for i in range(numnodes - 1):
		if NumEdgesIn[i] <= 1 and NumEdgesOut[i] <= 1:
			removestarts.append(changespos[i])
			removestops.append(changespos[i + 1])
			RemoveNodes[i] = True
	
	#deals with last node, so when $numnodes-1 the last remove stop is the end of the matrix
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


#uses position in matrix
def PrintResults():
	stderr.write("start\tstop\tID\tname\n")
	stderr.write("----------------------------------------\n")
	for i in range(len(Changes_orig)):
		if str(Changes_orig[i]) != 'skip':
			stderr.write(str(Columns_orig[ChangesPos_orig[i]]))
			stderr.write("\t")
			stderr.write(str(Columns_orig[ChangesPos_orig[i+1]-1]))
			stderr.write("\t")
			stderr.write(str(IDs[Columns_orig[ChangesPos_orig[i]]]))
			stderr.write("\t")
			stderr.write(str(Changes_orig[i]))
			stderr.write("\n")

	
#uses position in input sequence
def PrintResultsSequence():
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(Changes_orig)):
		if str(Changes_orig[i]) != 'skip':
			stderr.write(str(Columns_orig[ChangesPos_orig[i]]+startall))
			stderr.write("\t")
			stderr.write(str(Columns_orig[ChangesPos_orig[i+1]-1]+startall))
			stderr.write("\t")
			stderr.write(str(IDs[Columns_orig[ChangesPos_orig[i]]]))
			stderr.write("\t")
			stderr.write(str(Changes_orig[i]))
			stderr.write("\n")
			

#Steps- 
#1.create confidence for nodes
#	will be in a matrix that is #subfams x #nodes
#2.create path graph
#3.find alternative paths through graph and add those edges
#4.extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
#5.annotate again with removed subfams
#   --stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left

count: int = 0
while (True):
    count += 1
    numnodes = len(Changes)
    
    #breakout of loop if there are 2 or less nodes left
    if (numnodes <= 2):
        break 
	
    NodeConfidenceDict.clear()
    
    #initializes and fills node confidence matrix 
    NodeConfidence(NodeConfidenceDict, SubfamSeqs, ChromSeqs, ChangesPos)
    
    pathGraph.clear()
    FillPathGraph(pathGraph)
    
	#test to see if there nodes in the graph that have more that one incoming or outgoing edge,
	#if so keep looping, if not break out of the loop
	#if they are all 0, break out of the loop
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
        break 

    RemoveStarts.clear()
    RemoveStops.clear()

    ExtractNodes(RemoveStarts, RemoveStops, ChangesPos, pathGraph, numnodes)
    
    # removing inserted elements from @Columns so they can be ignored 
    total: int = 0
    for i in range(len(RemoveStops)):
    	del Columns[RemoveStarts[i]-total:RemoveStops[i]-total]
    	# 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
    	total += (RemoveStops[i] - RemoveStarts[i])
     
    # using prob matrix and origin matrix, just skip the cols I'm not interested in and annotate
	# without the removed subfam
	# using old prob matrix and origin matrix        
	# this time ignores inserted subfam because there are values in @RemoveStarts and @RemoveStops
    FillProbMatrix(ProbHash, SupportHashCollapse, OriginHash)

    Changes.clear()
    ChangesPos.clear()
        
    GetPath(ProbHash, OriginHash, SubFams)
    
#     print(Changes)
#     print()


if printMatrixPos:
	PrintResults()
else:
	PrintResultsSequence()

# PrintMatrixHashCollapse(cols, SupportHashCollapse)


        


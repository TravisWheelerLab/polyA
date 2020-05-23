from getopt import getopt
from math import inf, log
import re
from sys import argv, stderr, stdout
from typing import Dict, List, Tuple, Union

from polyA.load_alignments import load_alignments

#-----------------------------------------------------------------------------------#
#			FUNCTIONS																#
#-----------------------------------------------------------------------------------#


#just for debugging so can look at values in matrices
def PrintMatrixHashCollapse(column: int, Hash: Dict, subfamscollapse: Dict[str, int]) -> None:
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

    for k in subfamscollapse:
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


#just for debugging so can look at values in matrices
def PrintMatrixHash(column: int, subfams: List[str], Hash: Dict) -> None:
    stdout.write("\t")
    j: int = 0
    while j < column:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    i: int = 0
    while i < rows:
        stdout.write(f"{subfams[i]:<20s}\t")
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
#pass in SubfamSeqs and ChromSeqs and changes the actual strings int the lists 
def padSeqs(start: List[int], stop: List[int], subfamseq: List[str], chromseq: List[str]) -> Tuple[int, int]:

	edgestart: int = 0
	edgestop: int = 0
	
	(edgestart, edgestop) = Edges(start, stop)
	
	for i in range(1, len(subfamseq)):
		leftpad = start[i] - edgestart
		rightpad = edgestop - stop[i]
		
		chromseq[i] = ("." * leftpad) + f"{chromseq[i]}" + ("." * (rightpad + 15))
		subfamseq[i] = ("." * leftpad) + f"{subfamseq[i]}" + ("." * (rightpad + 15))
		
	return (edgestart, edgestop)
	

#lastpreva and lastprevb will be single chars        
def CalcScore(gapext: int, gapinit: int, submatrixcols: int, seq1: str, seq2: str, lastpreva: str, lastprevb: str, submatrix: Dict[int, int], charpos: Dict[str, int]) -> int:
    chunkscore: int = 0

	#deals with the first character of a segment being a gap character - have to look at last
	#segment to see if this is a gap init or ext
    if seq1[0] == "-":
        if lastpreva == "-":
            chunkscore += gapext
        else:
            chunkscore += gapinit
    elif seq2[0] == "-":
        if lastprevb == "-":
            chunkscore += gapext
        else:
            chunkscore += gapinit
    elif seq1[0] == "." or seq2[0] == ".":
        chunkscore = chunkscore
    else:
        chunkscore += submatrix[charpos[seq1[0]] * submatrixcols + charpos[seq2[0]]]
        
    for j in range(1, len(seq1)):
        if seq1[j] == "-":
            if seq1[j - 1] == "-":
                chunkscore += gapext
            else:
                chunkscore += gapinit
        elif seq2[j] == "-":
            if seq2[j - 1] == "-":
                chunkscore += gapext
            else:
                chunkscore += gapinit
        elif seq1[j] == "." or seq2[j] == ".":
            chunkscore = chunkscore
        else:
            if seq1[j] in CharPos and seq2[j] in charpos:
                chunkscore += submatrix[charpos[seq1[j]] * submatrixcols + charpos[seq2[j]]]

    return chunkscore


#calculates alignment score (according to crossmatch scoring) for every segment and all seqs
#fills align score matrix - pos 0 in the array is the score of the segment that starts at pos 0,
#pos 1 in array is score of segment that starts at pos 1, etc
#computes score for the first segment that does not start with a '.' and from there keeps the 
# base score and adds new char and subtracts old char - if a new gap is introduced, calls CalcScore()
#instead of adding on base score 
##FIXME - come back to this one later once I have the results printed - easier for testing 
def FillAlignScoreMatrix(edgestart: int, chunk_size: int, gapext: int, gapinit: int, skipalignscore: int, submatrixcols: int, subfams: List[str], chroms: List[str], starts: List[int], submatrix: Dict[int, int], charpos: Dict[str, int]) -> Tuple[int, Dict[Tuple[int, int], int]]:

#send in Starts, startall, chunksize, gapExt, gapInit, SubMatrix, subMatrixCols, CharPos

    numcols: int = 0
    index: int = 0
    
    alignmatrix: Dict[Tuple[int, int], int] = {}
    
    #chunks can't start on gaps and gaps don't count when getting to the 30 bps
	
    for i in range(1, len(chroms)):
        subfam1: str = subfams[i]
        chrom1: str = chroms[i]
        		
		#starts at the first non '.' char, but offsets it in the matrix based on where
		#the alignments start in the seq - ex: if first alignment in the seq starts at 10,
		#will offset by 10
        
        #calculates score for the first 16 chunks, chunk 16 is the first one that is 30nt
        #FIXME - this could start with smallest chunk and just add score of next nt to get next, etc
        #FIXME - but right now it takes all the chunks separately and runs them through CalcScore()
        j: int = starts[i] - edgestart
        index = j+15  #index is the col we are in the align score matrix, $j is the place in @subfam1 and @chrom1
        alignscore: int = 0
        tempindex: int = j
        tempcount: int = 0

        for k in range(15, -1, -1):
        
        	offset: int = chunk_size - k
        	alignscore = 0

        	tempindex = j
        	tempcount = 0

        	while tempcount < chunk_size - k:
        		if chrom1[tempindex] != "-":
        			tempcount += 1
        		tempindex += 1
        		
        	offset = tempindex - j
        	prevoffset = offset

			#grabs first chunk - here $ j = pos of first non '.' char
        	ChromSlice: str = chrom1[j:j + offset]
        	SubfamSlice: str = subfam1[j:j + offset]

			#calculates score for first chunk and puts score in alignhash
        	alignscore = CalcScore(gapext, gapinit, submatrixcols, SubfamSlice, ChromSlice, "", "", submatrix, charpos)
        	alignmatrix[i, index-k] = int(alignscore * chunk_size / (chunk_size - k)) # already to scale so don't need to * 31 and / 31
        
        index+=1 
        
        numnucls: int = chunk_size  #how many nucls contributed to align score 

        # TODO: Make sure these bounds are right since Python indexing is different
		#move to next chunk by adding next chars score and subtracting prev chars score
        while j + offset < len(chrom1):
            tempindex = j
            tempcount = 0

            while tempcount < chunk_size:
                if chrom1[tempindex + 1] != "-":
                    tempcount += 1
                tempindex += 1

            offset = tempindex - j

            if chrom1[j + 1] != "-":
                if chrom1[j + 1] != "." and subfam1[j + 1] != ".":
                    if prevoffset != offset:  #there is a new gap, or a gap was removed from beginning 
                        ChromSlice = chrom1[j+1:j+offset+1]
                        SubfamSlice = subfam1[j+1:j + offset + 1]
                        alignscore = CalcScore(gapext, gapinit, submatrixcols, SubfamSlice, ChromSlice, subfam1[j], chrom1[j], submatrix, charpos)
                        
                        tempcount2: int = 0
                        for nuc in ChromSlice:
                        	if nuc != '-' and nuc != '.':
                        		tempcount2 += 1
                        numnucls = tempcount2 #resetting numnucls to 31
                        
                        if numnucls < 16:
                        	alignscore = -inf
                        
                    else:
 						#alignscore from previous segment - prev chars score + next chars score 
						#subtracting prev  chars score - tests if its a gap in the subfam as well

                        if subfam1[j] == "-":
                        	numnucls -= 1
                        	if subfam1[j - 1] == "-":
                        		alignscore = alignscore - gapext
                        	else:
                        		alignscore = alignscore - gapinit
                        else:
                            alignscore = alignscore - submatrix[charpos[subfam1[j]]*submatrixcols+charpos[chrom1[j]]]
                            numnucls -= 1
                            
						#adding next chars score - tests if its a gap in the subfam as well
                        if subfam1[j + offset - 15] == "." or chrom1[j + offset - 15] == ".":
                            alignscore = -inf
                        elif subfam1[j + offset] == "-":
                        	numnucls += 1
                        	if subfam1[j + offset - 1] == "-":
                        		alignscore = alignscore + gapext
                        	else:
                        		alignscore = alignscore + gapinit
                        elif subfam1[j + offset] == "." or chrom1[j + offset] == ".":
                        	alignscore = alignscore
                        else:
                            alignscore = alignscore + SubMatrix[charpos[subfam1[j + offset]] * submatrixcols + charpos[chrom1[j + offset]]]
                            numnucls += 1
                    
                    if alignscore <= 0:
                        alignmatrix[i, index] = 1
                    else:
                        alignmatrix[i, index] = int(alignscore / numnucls * chunk_size)
                        
                    if alignscore == -inf:
                    	del alignmatrix[i, index]
                    	break

                index += 1

            j += 1
            prevoffset = offset
                    
        #max index is assigned to cols
        if numcols < index:
        	numcols = index
        	
    #assigns skip states an alignment score 
    for j in range(numcols):
    	alignmatrix[0, j] = skipalignscore
    	
    return(numcols, alignmatrix)



#fills parallel array to the Align Matrix that holds the consensus position for each 
# subfam at that position in the alignment
def FillConsensusPosMatrix(colnum: int, subfams: List[str], chroms: List[str], consensusstart: List[int], consensusstop: List[int], strands: List[str]) -> Dict[Tuple[int, int], int]:
                           
    consensus: Dict[Tuple[int, int], int] = {}
    
    #0s for consensus pos of skip state
    for j in range(colnum):
    	consensus[0, j] = 0
	
	#start at 1 to skip 'skip state'
    i: int = 1
    while i < rows:
        SubfamArray: str = subfams[i]
        ChromArray: str = chroms[i]

        consensuspos: int = 0
        if strands[i] == "+":
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
        
    return consensus


#puts all columns that are not empty into @Columns, so when I loop through hash I can use the 
#vals in @Columns - this will skip over empty columns
def FillColumns(num_cols: int, num_rows: int, alignmatrix: Dict[Tuple[int, int], int]) -> List[int]:
	columns: List[int] = []
	j: int = 0
	for j in range(num_cols):
		empty = 1;
		for i in range(num_rows):
			if (i, j) in alignmatrix:
				empty = 0
				i = num_rows
 	 	
		if not empty:
			columns.append(j)
		
	return columns


#FIXME - make this return an array of ints instead of a string
#send in an array of scores for a segment - output an array of confidence values for the segment
def ConfidenceCM(lamb: float, region: List[int]) -> str:
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


def FillConfScoreMatrix(rownum: int, lamb: float, columns: List[int], alignhash: Dict[Tuple[int, int], int]) -> Dict[Tuple[int, int], float]:
	confhash: Dict[Tuple[int, int], float] = {}
	
	for i in range(len(columns)):
		if i >= len(columns):
			break
			
		col: int = columns[i]
		temp: List[int] = []
		row: int = 0
		
		while row < rownum:
		
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
	
	return confhash


# Fills support score matrix using values in conf matrix
#score for subfam x at position i is sum of all confidences for subfam x for all segments that 
#overlap position i - divided by number of segments
def FillSupportMatrix(row_num: int, columns: List[int], alignhash: Dict[Tuple[int, int], int], confhash: Dict[Tuple[int, int], float]) -> Dict[Tuple[int, int], float]:
    
    supporthash: Dict[Tuple[int, int], float] = {}

    # i = subfam, j = position
    i: int = 0
    while i < row_num:
        tempcol: int = -1
        for col in range(len(columns)):
            if col >= len(columns):
                break
            j = columns[col]

            if (i, j) in confhash:
                num: int = j
                summ: float = 0.0
                numsegments: int = 0
                while num >= 0 and num >= j:
                    if (i, num) in confhash:
                        summ = summ + confhash[i, num]
                        numsegments += 1
                    num -= 1

                if numsegments > 0:
                    supporthash[i, j] = summ / numsegments
            col += 1
        i += 1
        
    return supporthash


#collapses matrices 
#collapse and combine rows that are the same subfam - just sum their support 
#new support dict has key = subfamname.col 
#also creates a bookkeeping dict that has all the cols as keys and their values 
#are arrays that hold all the active subfams in that col - used so that don't have 
#to loop through all the $i's just to see if a column exists 
def CollapseMatrices(row: int, columns: List[int], subfams: List[str], strands: List[str], supporthash: Dict[Tuple[int, int], float], consensushash: Dict[Tuple[int, int], int]) -> Tuple[int, Dict[Tuple[str, int], int], Dict[Tuple[str, int], str], Dict[Tuple[int, int], float], Dict[str, int], Dict[int, List[str]]]:
	
	row_update: int = 0
	consensushashcollapse: Dict[Tuple[str, int], int] = {}
	strandhashcollapse: Dict[Tuple[str, int], str] = {}
	supporthashcollapse: Dict[Tuple[str, int], float] = {}
	subfamscollapse: Dict[str, int] = {}
	activecellscollapse: Dict[int, List[str]] = {}
	
	#input Columns, SupportHash, ConsensusHash, Strands, 
	#updates rows
	#return ConsensusHashCollapse, StrandHashCollapse, SupportHashCollapse, 
	#SubFamsCollapse, ActiveCellsCollapse
	
	for col in range(len(columns)):
		j: int = columns[col]
		DupMaxCon: Dict[str, float] = {}
		DupMaxSup: Dict[str, float] = {}
	
		activecols: List[str] = []
		activecellscollapse[j] = activecols
	
		#sum the support score for rows that are collapsed together
		#find max support score for collapsed rows and use the consensus from that row
		for i in range(row):
			if (i,j) in supporthash and (i,j) in consensushash:
				if (subfams[i]) in DupMaxCon:
					if supporthash[i,j] > DupMaxCon[subfams[i]]:
						DupMaxCon[subfams[i]] = supporthash[i,j]
						consensushashcollapse[subfams[i],j] = consensushash[i,j]
						strandhashcollapse[subfams[i], j] = strands[i]
				else:
					DupMaxCon[SubFams[i]] = supporthash[i,j]
					consensushashcollapse[subfams[i],j] = consensushash[i,j]
					strandhashcollapse[subfams[i], j] = strands[i]
		
			if (i,j) in supporthash:
				if (subfams[i]) in DupMaxSup:
					if supporthash[i,j] > DupMaxSup[subfams[i]]:
						DupMaxSup[SubFams[i]] = supporthash[i,j]
						supporthashcollapse[subfams[i],j] = supporthash[i,j]
				else:
					DupMaxSup[subfams[i]] = supporthash[i,j]
					supporthashcollapse[subfams[i],j] = supporthash[i,j]
					activecellscollapse[j].append(subfams[i])
					
	for i in range(row):
		subfamscollapse[subfams[i]] = 0

	#update var rows after collapse 
	row_update = len(subfamscollapse)

	return (row_update, consensushashcollapse, strandhashcollapse, supporthashcollapse, subfamscollapse, activecellscollapse)


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
def FillProbMatrix(sameprobskip: float, sameproblog: float, changeproblog: float, changeprobskip: float, columns: List[int], subfamscollapse: Dict[str, int], activecellscollapse: Dict[int, List[str]], supporthash: Dict[Tuple[str, int], float], strandhashcollapse: Dict[Tuple[str, int], str], consensushashcollapse: Dict[Tuple[str, int], int]) -> Tuple[Dict[Tuple[str, int], float], Dict[Tuple[str, int], str]]:
	
	probhash: Dict[Tuple[str, int], float] = {}
	originhash: Dict[Tuple[str, int], str] = {}
	
	#fill first col of probhash with 0s
	for k in subfamscollapse:
		probhash[k, 0] = 0

	for j in range(1, len(columns)):
	 			
		for i in activecellscollapse[columns[j]]:
			max: float = -inf
			maxindex: str = ''
			supportlog: float = log(supporthash[i, columns[j]])
			
			#loop through all the subfams in the previous column
			for row in activecellscollapse[columns[j-1]]:
				score: float = supportlog + probhash[row, columns[j-1]]
				prob: int = 0
				
				if row == i:  #staying in same row
					
					if row == 'skip':  #staying in skip
						prob = sameprobskip
						
					prob = sameproblog

					# because rows are collapsed, if the consensus seqs are contiguous - treat as if they are not the same row and get the jump penalty
					if strandhashcollapse[row, columns[j-1]] != strandhashcollapse[i, columns[j]]:
						prob = changeproblog
					else: #if on same strand
						if strandhashcollapse[row, columns[j-1]] == '+':
							if consensushashcollapse[row, columns[j-1]] > consensushashcollapse[i, columns[j]] + 50:
								prob = changeproblog
						if strandhashcollapse[row, columns[j-1]] == '-':
							if consensushashcollapse[row, columns[j-1]] + 50 < consensushashcollapse[i, columns[j]]:
								prob = changeproblog
						
				else:  #jumping rows
					prob = changeproblog  
					if row == 'skip' or i == 'skip':  #jumping in or out of skip
						prob = changeprobskip
					
				score = score + prob
					
				if score > max:
					max = score
					maxindex = row
					
			probhash[i, columns[j]] = max
			originhash[i, columns[j]] = maxindex
			
	return (probhash, originhash)


#using origin matrix, back traces through the 2D array to get the subfam path
#finds where the path switches to a different row and populates @Changes and @ChangesPos
#reverses @Changes and @ChangesPos because it's a backtrace so they are initially backwards
#jumps over removed columns when necessary
def GetPath(temp_id: int, columns: List[int], ids: List[int], subfams: List[str], activecellscollapse: Dict[int, List[str]], probhash: Dict[Tuple[str, int], float], originhash: Dict[Tuple[str, int], str]) -> Tuple[int, List[int], List[str]]:
    maxxx: float = -inf
    maxindex: str = ''
    
    changespos: List[int] = []
    changes: List[str] = []
    
    for i1 in activecellscollapse[cols - 1]:
        if maxxx < probhash[i1, cols - 1]:
            maxxx = probhash[i1, cols - 1]
            maxindex = i1
    
    prev: str = originhash[maxindex, cols - 1]
    
    changespos.append(len(columns))
    
    #already added the last col, but this adds the one before $col so still start at last col
    j: int = len(columns)-1
    while j > 1:
    	if (prev, columns[j-1]) in originhash:
    	
    		ids[columns[j-1]] = temp_id
    		
    		if prev != originhash[prev, Columns[j-1]]:
    			temp_id+=1234
    			changespos.append(j-1)
    			changes.append(prev)
    			
    		prev = originhash[prev, Columns[j-1]]
    		
    	j-=1
    	
    	
    ids[columns[j-1]] = temp_id
    changespos.append(j-1)
    changes.append(prev)
    
    changes.reverse()
    changespos.reverse()
    
    #changes ID for next round of stitching, so when starts stitching will have unique ID
    temp_id+=1234
    
    return (temp_id, changespos, changes)

#for debugging
def PrintChanges(columns: List[int], changes: List[str], changespos: List[int]) -> None:
    i: int = 0
    while i < len(changes):
        stdout.write(str(columns[changespos[i]]))
        stdout.write("\t")
        stdout.write(f"{changes[i]}\n")
        i=i+1


#FIXME - if goes into skip state for just one position, this will have an error .. also when
#stitching want to ignore skip states
#fills node confidence matrix 
#first fills matrix with node alignment scores, then reuses matrix for confidence scores 
def NodeConfidence(nodes: int, gapext: int, gapinit: int, submatrixcols: int, lamb: float, columns: List[int], subfamseqs: List[str], chromseqs: List[str], changespos: List[int], subfams: List[str], submatrix: Dict[int, int], charpos: Dict[str, int]) -> Dict[Tuple[str, int], float]:   
    nodeconfidence_temp: List[float] = [0 for _ in range(len(subfams) * nodes)]
    
    nodeconfidence: Dict[Tuple[str, int], float] = {}
        
    #calculated node confidence for first node - doesn't look back at prev char bc there isn't one
    #starts at 1 bc don't need to calc for skip state
    for j in range(1, len(subfams)):
        b: int = columns[changespos[0]]
        e: int = columns[changespos[1]]
        subfam: str = subfamseqs[j][b:e]
        chrom: str = chromseqs[j][b:e]
        alignscore: float = CalcScore(gapext, gapinit, submatrixcols, subfam, chrom, '', '', submatrix, charpos)
        nodeconfidence_temp[j * nodes + 0] = alignscore
	 
	#does rest of nodes - looks back at prev char incase of gap ext
    for i in range(1, nodes-1):
    	for j in range(1, len(subfams)):
            b: int = columns[changespos[i]]
            e: int = columns[changespos[i + 1]]
            subfam: str = subfamseqs[j][b:e]
            chrom: str = chromseqs[j][b:e]
            lastpreva: str = subfamseqs[j][b-1]#subfamseqs[j][Columns[changespos[i + 1] - 1]]
            lastprevb: str = chromseqs[j][b-1]#chromseqs[j][changespos[i + 1] - 1]
            alignscore: float = CalcScore(gapext, gapinit, submatrixcols, subfam, chrom, lastpreva, lastprevb, submatrix, charpos)
            nodeconfidence_temp[j * nodes + i] = alignscore
    
	#does last node
    for j in range(1, len(subfams)):
        b: int = columns[changespos[-2]]
        e: int = columns[changespos[-1]-1]
        subfam: str = subfamseqs[j][b:e]
        chrom: str = chromseqs[j][b:e]
        lastpreva: str = subfamseqs[j][b-1]
        lastprevb: str = chromseqs[j][b-1]
        alignscore: float = CalcScore(gapext, gapinit, submatrixcols, subfam, chrom, lastpreva, lastprevb, submatrix, charpos)
        nodeconfidence_temp[j * nodes + nodes-1] = alignscore
 
    #reuse same matrix and compute confidence scores for the nodes	
    for j in range(nodes):
    	temp: List[float] = []
    	for i in range(len(subfams)):
    		temp.append(nodeconfidence_temp[i * nodes + j])
    		
    	confidenceTemp: List[float] = [float(x) for x in ConfidenceCM(lamb, temp).split(' ')]
    	
    	for i in range(len(subfams)):
    		nodeconfidence_temp[i * nodes + j] = confidenceTemp[i]

	#collapse nodeconfidence down same way supportmatrix is collapsed - all seqs of 
	#the same subfam are put in the same row
	#not a sparse hash - holds the 0s, but I think this is okay because it won't ever
	#be a very large matrix, and this way we don't have to test if anything exists 	
	# my %nodeconfidence_collapse;
    for j in range(nodes):
        for i in range(len(subfams)):
            if (subfams[i], j) in nodeconfidence:
                nodeconfidence[subfams[i], j] += nodeconfidence_temp[i * nodes + j]
            else:
            	nodeconfidence[subfams[i], j] = nodeconfidence_temp[i * nodes + j]
            	
    return nodeconfidence


#used for debugging
def PrintPathGraph(nodes: int, changes: List[str], pathgraph: List[int]) -> None:
	stdout.write(" ")
	for i in range(nodes):
		stdout.write(f"{changes[i]} ")
	stdout.write("\n")
	
	for i in range(nodes):
		for j in range(nodes):
			stdout.write(f"{pathgraph[i*nodes+j]}\t")
		stdout.write(f"{changes[i]}\n")
	stdout.write("\n")


#used for debugging
def PrintNodeConfidence(nodes: int, changes: List[str], subfamscollapse: Dict[str, int], nodeconfidencedict: Dict[Tuple[str, int], float]) -> None:
	for i in range(nodes):
		stdout.write(f"{changes[i]} ")
	stdout.write("\n")

	for subfam in subfamscollapse:
		stdout.write(f"{subfam} ")
		for j in range(nodes):
			if (subfam,j) in nodeconfidencedict:
				stdout.write(f"{nodeconfidencedict[subfam,j]} ")
			else:
				stdout.write(f"-inf ")
		stdout.write("\n")


def FillPathGraph(nodes: int, columns: List[int], changes: List[str], changespos: List[int], subfamscollapse: Dict[str, int], consensushashcollapse: Dict[Tuple[str, int], int], strandhashcollapse: Dict[Tuple[str, int], str], nodeconfidencedict: Dict[Tuple[str, int], float]) -> List[int]:
	
    pathgraph: List[int] = []
    
    for i in range(nodes * nodes):
        pathgraph.append(0)

	# filling beginning path graph with straight line through the nodes
    for i in range(nodes - 1):
        pathgraph[i * nodes + i + 1] = 1

    for j in range(nodes):
        sinkSubfam: str = str(changes[j])
        sinkSubfamStart: int = consensushashcollapse[sinkSubfam, columns[changespos[j]]]
        sinkStrand: str = strandhashcollapse[sinkSubfam, columns[changespos[j]]]

		#looks at all the preceding nodes, except the one directly before it (source nodes)
        for i in range(j - 1):
        	#look at all the subfams in each node
            for sourceSubfam in subfamscollapse:
            	sourceSubfam = str(sourceSubfam)
            	sourceConf = nodeconfidencedict[sourceSubfam, i]
            	
            	if (sourceSubfam, columns[changespos[i + 1]-1]) in consensushashcollapse:
            		sourceSubfamStop = consensushashcollapse[sourceSubfam, columns[changespos[i + 1]-1]]
            		sourceStrand = strandhashcollapse[sourceSubfam, columns[changespos[i + 1]-1]]
  
  					# adds in edge if the subfam of the sink is at the source node and if it's 
					# confidence >= 1%, and if the source is before the sink in the consensus sequence           		
            		if sinkStrand == '+' and sinkStrand == sourceStrand:
            			if (sinkSubfam == sourceSubfam) and (sourceConf >= 0.3):
            				#FIXME- not sure what this overlap should be .. just allowed 50 for now            				
            				if sourceSubfamStop <= sinkSubfamStart + 50:
            					pathgraph[i * nodes + j] = 1
            					
            		elif sinkStrand == '-' and sinkStrand == sourceStrand:
            			if sinkSubfam == sourceSubfam and sourceConf >= 0.3:
            				if sourceSubfamStop >= sinkSubfamStart + 50:
            					pathgraph[i * nodes + j] = 1
            					
    return pathgraph


#finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds them to
# @RemoveStarts and @RemoveStops so they can be extracted from the alignment 
#pass in columns and updates it to ignore removed nodes
# updates cols if removes node from end 
def ExtractNodes(col: int, nodes: int, columns: List[int], changespos: List[int], pathgraph: List[int]) -> int:
	
	removestarts: List[int] = []
	removestops: List[int] = []

	#boolean for which nodes will be removed
	RemoveNodes: List[bool] = [False for _ in range(nodes)]
	
	#extracting nodes that only have one incoming and one outgoing edge
	NumEdgesIn: List[int] = [0 for _ in range(nodes)]
	NumEdgesOut: List[int] = [0 for _ in range(nodes)]
	
	for i in range(nodes):
		for j in range(nodes):
			#$j - incoming
			NumEdgesIn[j] += pathgraph[i * nodes + j]
			#$j - outgoing
			NumEdgesOut[i] += pathgraph[i * nodes + j]
			
	for i in range(nodes - 1):
		if NumEdgesIn[i] <= 1 and NumEdgesOut[i] <= 1:
			removestarts.append(changespos[i])
			removestops.append(changespos[i + 1])
			RemoveNodes[i] = True
	
	#deals with last node, so when $numnodes-1 the last remove stop is the end of the matrix
	if NumEdgesIn[nodes - 1] <= 1 and NumEdgesOut[nodes - 1] <= 1:
		removestarts.append(changespos[nodes - 1])
		removestops.append(len(columns)-1)
		RemoveNodes[nodes - 1] = True
		
	#when removing from the end, have to update cols because don't want to do to the end of the matrix anymore 
	i: int = nodes - 1
	while RemoveNodes[i]:
		col = columns[changespos[i] - 1]  #FIXME - do I index cols here or not?
		i -= 1
		
	# removing inserted elements from @Columns so they can be ignored 
	total: int = 0
	for i in range(len(removestops)):
		del columns[removestarts[i]-total:removestops[i]-total]
		# 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
		total += (removestops[i] - removestarts[i])

	return col


#uses position in matrix
def PrintResults(changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
	stderr.write("start\tstop\tID\tname\n")
	stderr.write("----------------------------------------\n")
	for i in range(len(changes_orig)):
		if str(changes_orig[i]) != 'skip':
			stderr.write(str(columns_orig[changespos_orig[i]]))
			stderr.write("\t")
			stderr.write(str(columns_orig[changespos_orig[i+1]-1]))
			stderr.write("\t")
			stderr.write(str(ids[columns_orig[changespos_orig[i]]]))
			stderr.write("\t")
			stderr.write(str(changes_orig[i]))
			stderr.write("\n")

	
#uses position in input sequence
def PrintResultsSequence(edgestart: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(changes_orig)):
		if str(changes_orig[i]) != 'skip':
			stderr.write(str(columns_orig[changespos_orig[i]]+edgestart))
			stderr.write("\t")
			stderr.write(str(columns_orig[changespos_orig[i+1]-1]+edgestart))
			stderr.write("\t")
			stderr.write(str(ids[columns_orig[changespos_orig[i]]]))
			stderr.write("\t")
			stderr.write(str(changes_orig[i]))
			stderr.write("\n")


#---------------------------------------------------------------------------------------#
#GLOBALS - yes I know these are bad, will be passed into functions to make port easier	#
#---------------------------------------------------------------------------------------#

gapInit: int = -25
gapExt: int = -5
Lamb: float = 0.1227  # From command line
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
    "matrixPos",
])
opts = dict(raw_opts)

gapInit = int(opts["--gapInit"]) if "--gapInit" in opts else gapInit
gapExt = int(opts["--gapExt"]) if "--gapExt" in opts else gapExt
Lamb = float(opts["--lambda"]) if "--lambda" in opts else Lamb
chunksize = int(opts["--segmentsize"]) if "--segmentsize" in opts else chunksize
changeProb = float(opts["--changeprob"]) if "--changeprob" in opts else changeProb
help = "--help" in opts
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

Columns: List[int] = [];

# RemoveStarts: List[int] = []
# RemoveStops: List[int] = []

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

#for graph/node part
numnodes: int = 0
NodeConfidenceDict: Dict[Tuple[str, int], float] = {}
pathGraph: List[int] = []
total: int = 0
loop: int = 1



#-----------------------------------------------------------------------------------#
#			MAIN																	#
#-----------------------------------------------------------------------------------#

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
sameProbSkip = changeProbLog / 10; # 10% of the jump penalty, staying in skip state for 20nt "counts" as one jump

#precomputes global vars rows and cols in matrices 
rows: int = len(SubFams)
cols: int = 0 #assign cols in FillAlignScoreMatrix

(startall, Stopall) = padSeqs(Starts, Stops, SubfamSeqs, ChromSeqs)
                    
(cols, AlignHash) = FillAlignScoreMatrix(startall, chunksize, gapExt, gapInit, skipAlignScore, subMatrixCols, SubfamSeqs, ChromSeqs, Starts, SubMatrix, CharPos)

ConsensusHash = FillConsensusPosMatrix(cols, SubfamSeqs, ChromSeqs, ConsensusStarts, ConsensusStops, Strands)

Columns = FillColumns(cols, rows, AlignHash)

ConfHash = FillConfScoreMatrix(rows, Lamb, Columns, AlignHash)

SupportHash = FillSupportMatrix(rows, Columns, AlignHash, ConfHash)

(rows, ConsensusHashCollapse, StrandHashCollapse, SupportHashCollapse, SubFamsCollapse, ActiveCellsCollapse) = CollapseMatrices(rows, Columns, SubFams, Strands, SupportHash, ConsensusHash)

(ProbHash, OriginHash) = FillProbMatrix(sameProbSkip, sameProbLog, changeProbLog, changeProbSkip, Columns, SubFamsCollapse, ActiveCellsCollapse, SupportHashCollapse, StrandHashCollapse, ConsensusHashCollapse)

IDs = [0] * cols

(ID, ChangesPos, Changes) = GetPath(ID, Columns, IDs, SubFams, ActiveCellsCollapse, ProbHash, OriginHash)

#keep the original annotation for reporting results
Changes_orig = Changes.copy()
ChangesPos_orig = ChangesPos.copy()
Columns_orig = Columns.copy()			

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
    NodeConfidenceDict = NodeConfidence(numnodes, gapExt, gapInit, subMatrixCols, Lamb, Columns, SubfamSeqs, ChromSeqs, ChangesPos, SubFams, SubMatrix, CharPos)
        
    pathGraph.clear()
    pathGraph = FillPathGraph(numnodes, Columns, Changes, ChangesPos, SubFamsCollapse, ConsensusHashCollapse, StrandHashCollapse, NodeConfidenceDict)
    
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

    cols = ExtractNodes(cols, numnodes, Columns, ChangesPos, pathGraph)
         
    # using prob matrix and origin matrix, just skip the cols I'm not interested in and annotate
	# without the removed subfam
	# using old prob matrix and origin matrix        
	# this time ignores inserted subfam because there are values in @RemoveStarts and @RemoveStops
    (ProbHash, OriginHash) = FillProbMatrix(sameProbSkip, sameProbLog, changeProbLog, changeProbSkip, Columns, SubFamsCollapse, ActiveCellsCollapse, SupportHashCollapse, StrandHashCollapse, ConsensusHashCollapse)

    Changes.clear()
    ChangesPos.clear()
        
    (ID, ChangesPos, Changes) = GetPath(ID, Columns, IDs, SubFams, ActiveCellsCollapse, ProbHash, OriginHash)


if printMatrixPos:
	PrintResults(Changes_orig, ChangesPos_orig, Columns_orig, IDs)
else:
	PrintResultsSequence(startall, Changes_orig, ChangesPos_orig, Columns_orig, IDs)


        


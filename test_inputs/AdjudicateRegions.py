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
def PrintMatrixHashCollapse(num_col: int, matrix: Dict, subfams_collapse: Dict[str, int]) -> None:
    stdout.write("\t")

    j: int = 0
    while j < num_col:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

	#Subfams_collapse is not ordered, so prints the skip state first 
    stdout.write("skip\t")
    j: int = 0
    while j < num_col:
        if ("skip", j) in matrix:
            stdout.write(str(matrix["skip", j]))
        else:
            stdout.write(f"{-inf:5.3g}")
        stdout.write("\t")
        j += 1
    stdout.write("\n")

    for k in subfams_collapse:
        if k != "skip":
            stdout.write(f"{k:<20s}\t")
            j: int = 0
            while j < num_col:
                if (k, j) in matrix:
                    stdout.write(str(matrix[k, j]))
                else:
                    stdout.write(f"{-inf:5.3g}")
                stdout.write("\t")
                j += 1
            stdout.write("\n")


#just for debugging so can look at values in matrices
def PrintMatrixHash(num_col: int, num_row: int, subfams: List[str], matrix: Dict) -> None:
    stdout.write("\t")
    j: int = 0
    while j < num_col:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    i: int = 0
    while i < num_row:
        stdout.write(f"{subfams[i]:<20s}\t")
        j: int = 0
        while j < num_col:
            if (i, j) in matrix:
                stdout.write(f"{matrix[i, j]}")
            else:
                stdout.write(f"{-inf:5.3g}")
            stdout.write("\t")
            j += 1
        stdout.write("\n")
        i += 1


#find the min start and max stop for the whole region
def Edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    min_start: int = starts[1]
    max_stop: int = stops[1]

    for i in range(1, len(starts)):
        if starts[i] < min_start:
            min_start = starts[i]
        if stops[i] > max_stop:
            max_stop = stops[i]

    return min_start, max_stop


#pads sequences with '.' and makes sparse matrix - allows regions where seqs do not all have same 
#start and stop positions
#pass in SubfamSeqs and ChromSeqs and changes the actual strings int the lists 
def PadSeqs(start: List[int], stop: List[int], subfam_seqs: List[str], chrom_seqs: List[str]) -> Tuple[int, int]:

	edge_start: int = 0
	edge_stop: int = 0
	
	(edge_start, edge_stop) = Edges(start, stop)
	
	for i in range(1, len(subfam_seqs)):
		left_pad: int = start[i] - edge_start
		right_pad: int = edge_stop - stop[i]
		
		chrom_seqs[i] = ("." * left_pad) + f"{chrom_seqs[i]}" + ("." * (right_pad + 15))
		subfam_seqs[i] = ("." * left_pad) + f"{subfam_seqs[i]}" + ("." * (right_pad + 15))
		
	return (edge_start, edge_stop)
	

#lastpreva and lastprevb will be single chars        
def CalcScore(gap_ext: int, gap_init: int, submatrix_cols: int, seq1: str, seq2: str, prev_char_seq1: str, prev_char_seq2: str, sub_matrix: Dict[int, int], char_pos: Dict[str, int]) -> int:
    chunk_score: int = 0

	#deals with the first character of a segment being a gap character - have to look at last
	#segment to see if this is a gap init or ext
    if seq1[0] == "-":
        if prev_char_seq1 == "-":
            chunk_score += gap_ext
        else:
            chunk_score += gap_init
    elif seq2[0] == "-":
        if prev_char_seq2 == "-":
            chunk_score += gap_ext
        else:
            chunk_score += gap_init
    elif seq1[0] == "." or seq2[0] == ".":
        chunk_score = chunk_score
    else:
        chunk_score += sub_matrix[char_pos[seq1[0]] * submatrix_cols + char_pos[seq2[0]]]
        
    for j in range(1, len(seq1)):
        if seq1[j] == "-":
            if seq1[j - 1] == "-":
                chunk_score += gap_ext
            else:
                chunk_score += gap_init
        elif seq2[j] == "-":
            if seq2[j - 1] == "-":
                chunk_score += gap_ext
            else:
                chunk_score += gap_init
        elif seq1[j] == "." or seq2[j] == ".":
            chunk_score = chunk_score
        else:
            if seq1[j] in char_pos and seq2[j] in char_pos:
                chunk_score += sub_matrix[char_pos[seq1[j]] * submatrix_cols + char_pos[seq2[j]]]

    return chunk_score


#calculates alignment score (according to crossmatch scoring) for every segment and all seqs
#fills align score matrix - pos 0 in the array is the score of the segment that starts at pos 0,
#pos 1 in array is score of segment that starts at pos 1, etc
#computes score for the first segment that does not start with a '.' and from there keeps the 
# base score and adds new char and subtracts old char - if a new gap is introduced, calls CalcScore()
#instead of adding on base score 
##FIXME - come back to this one later once I have the results printed - easier for testing 
def FillAlignMatrix(edge_start: int, chunk_size: int, gap_ext: int, gap_init: int, skip_align_score: int, sub_matrix_cols: int, subfams: List[str], chroms: List[str], starts: List[int], sub_matrix: Dict[int, int], char_pos: Dict[str, int]) -> Tuple[int, Dict[Tuple[int, int], int]]:

    num_cols: int = 0
    col_index: int = 0
    
    align_matrix: Dict[Tuple[int, int], int] = {}
    
    #chunks can't start on gaps and gaps don't count when getting to the 30 bps
	
    for i in range(1, len(chroms)):
        subfam_seq: str = subfams[i]
        chrom_seq: str = chroms[i]
        		
		#starts at the first non '.' char, but offsets it in the matrix based on where
		#the alignments start in the seq - ex: if first alignment in the seq starts at 10,
		#will offset by 10
        
        #calculates score for the first 16 chunks, chunk 16 is the first one that is 30nt
        #FIXME - this could start with smallest chunk and seq_indexust add score of next nt to get next, etc
        #FIXME - but right now it takes all the chunks separately and runs them through CalcScore()
        seq_index: int = starts[i] - edge_start
        col_index = seq_index+int((ChunkSize-1)/2)  #col_index is the col we are in the align score matrix, $seq_index is the place in @subfam_seq and @chrom_seq
        align_score: int = 0
        temp_index: int = seq_index
        temp_count: int = 0

        for k in range(int((ChunkSize-1)/2), -1, -1):
        
        	offset: int = chunk_size - k
        	align_score = 0

        	temp_index = seq_index
        	temp_count = 0

        	while temp_count < chunk_size - k:
        		if chrom_seq[temp_index] != "-":
        			temp_count += 1
        		temp_index += 1
        		
        	offset = temp_index - seq_index
        	prev_offset = offset

			#grabs first chunk - here $ seq_index = pos of first non '.' char
        	chrom_slice: str = chrom_seq[seq_index:seq_index + offset]
        	subfam_slice: str = subfam_seq[seq_index:seq_index + offset]

			#calculates score for first chunk and puts score in alignhash
        	align_score = CalcScore(gap_ext, gap_init, sub_matrix_cols, subfam_slice, chrom_slice, "", "", sub_matrix, char_pos)
        	align_matrix[i, col_index-k] = int(align_score * chunk_size / (chunk_size - k)) # already to scale so don't need to * 31 and / 31
        
        col_index+=1 
        
        num_nucls: int = chunk_size  #how many nucls contributed to align score 

        # TODO: Make sure these bounds are right since Python indexing is different
		#move to next chunk by adding next chars score and subtracting prev chars score
        while seq_index + offset < len(chrom_seq):
            temp_index = seq_index
            temp_count = 0

            while temp_count < chunk_size:
                if chrom_seq[temp_index + 1] != "-":
                    temp_count += 1
                temp_index += 1

            offset = temp_index - seq_index

            if chrom_seq[seq_index + 1] != "-":
                if chrom_seq[seq_index + 1] != "." and subfam_seq[seq_index + 1] != ".":
                    if prev_offset != offset:  #there is a new gap, or a gap was removed from beginning 
                        chrom_slice = chrom_seq[seq_index+1:seq_index+offset+1]
                        subfam_slice = subfam_seq[seq_index+1:seq_index + offset + 1]
                        align_score = CalcScore(gap_ext, gap_init, sub_matrix_cols, subfam_slice, chrom_slice, subfam_seq[seq_index], chrom_seq[seq_index], sub_matrix, char_pos)
                        
                        temp_count2: int = 0
                        for nuc in chrom_slice:
                        	if nuc != '-' and nuc != '.':
                        		temp_count2 += 1
                        num_nucls = temp_count2 #resetting num_nucls to 31
                        
                        if num_nucls <= int((ChunkSize-1)/2):
                        	align_score = -inf
                        
                    else:
 						#align_score from previous segment - prev chars score + next chars score 
						#subtracting prev  chars score - tests if its a gap in the subfam as well

                        if subfam_seq[seq_index] == "-":
                        	num_nucls -= 1
                        	if subfam_seq[seq_index - 1] == "-":
                        		align_score = align_score - gap_ext
                        	else:
                        		align_score = align_score - gap_init
                        else:
                            align_score = align_score - sub_matrix[char_pos[subfam_seq[seq_index]]*sub_matrix_cols+char_pos[chrom_seq[seq_index]]]
                            num_nucls -= 1
                            
						#adding next chars score - tests if its a gap in the subfam as well
                        if subfam_seq[seq_index + offset - int((ChunkSize-1)/2)] == "." or chrom_seq[seq_index + offset - int((ChunkSize-1)/2)] == ".":
                            align_score = -inf
                        elif subfam_seq[seq_index + offset] == "-":
                        	num_nucls += 1
                        	if subfam_seq[seq_index + offset - 1] == "-":
                        		align_score = align_score + gap_ext
                        	else:
                        		align_score = align_score + gap_init
                        elif subfam_seq[seq_index + offset] == "." or chrom_seq[seq_index + offset] == ".":
                        	align_score = align_score
                        else:
                            align_score = align_score + sub_matrix[char_pos[subfam_seq[seq_index + offset]] * sub_matrix_cols + char_pos[chrom_seq[seq_index + offset]]]
                            num_nucls += 1
                    
                    if align_score <= 0:
                        align_matrix[i, col_index] = 1
                    else:
                        align_matrix[i, col_index] = int(align_score / num_nucls * chunk_size)
                        
                    if align_score == -inf:
                    	del align_matrix[i, col_index]
                    	break

                col_index += 1

            seq_index += 1
            prev_offset = offset
                    
        #max col_index is assigned to cols
        if num_cols < col_index:
        	num_cols = col_index
        	
    #assigns skip states an alignment score 
    for j in range(num_cols):
    	align_matrix[0, j] = skip_align_score
    	
    return(num_cols, align_matrix)



#fills parallel array to the Align Matrix that holds the consensus position for each 
# subfam at that position in the alignment
def FillConsensusPositionMatrix(col_num: int, subfams: List[str], chroms: List[str], consensus_starts: List[int], consensus_stops: List[int], strands: List[str]) -> Dict[Tuple[int, int], int]:
                           
    consensus_matrix: Dict[Tuple[int, int], int] = {}
    
    #0s for consensus pos of skip state
    for j in range(col_num):
    	consensus_matrix[0, j] = 0
	
	#start at 1 to skip 'skip state'
    for row_index in range(1, rows):

        consensus_pos: int = 0
        if strands[row_index] == "+":
            consensus_pos = consensus_starts[row_index] - 1
            col_index: int = 0
            
            for seq_index in range(len(subfams[row_index])):
                if subfams[row_index][seq_index] != ".":
                	#consensus pos only advances when there is not a gap in the subfam seq
                    if subfams[row_index][seq_index] != "-":
                        consensus_pos += 1

                   #put consensus pos corresponding to pos in matrix in hash
                    consensus_matrix[row_index, col_index] = consensus_pos

				#matrix position only advances when there is not a gap in the chrom seq
                if chroms[row_index][seq_index] != "-":
                    col_index += 1
                    
        else:  #reverse strand 
            consensus_pos2 = consensus_starts[row_index] + 1
            col_index2: int = 0
            
            for seq_index2 in range(len(subfams[row_index])):
                if subfams[row_index][seq_index2] != ".":
                    if subfams[row_index][seq_index2] != "-":
                        consensus_pos2 -= 1
                    consensus_matrix[row_index, col_index2] = consensus_pos2

                if chroms[row_index][seq_index2] != "-":
                    col_index2 += 1
        
    return consensus_matrix


#puts all columns that are not empty into @NonEmptyColumns, so when I loop through hash I can use the 
#vals in @NonEmptyColumns - this will skip over empty columns
def FillColumns(num_cols: int, num_rows: int, align_matrix: Dict[Tuple[int, int], int]) -> List[int]:
	columns: List[int] = []
	j: int = 0
	for j in range(num_cols):
		empty = 1;
		for i in range(1, num_rows):
			if (i, j) in align_matrix:
				empty = 0
				i = num_rows
 	 	
		if not empty:
			columns.append(j)
		
	return columns


#FIXME - make this return an array of ints instead of a string
#send in an array of scores for a segment - output an array of confidence values for the segment
def ConfidenceCM(lambdaa: float, region: List[int]) -> str:
    confidence_string: str = ""

	#loops through the array once to get the sum of 2^every_hit_score in region 
	#converts the score to account for lambda before summing 
    score_total: int = 0
    for score in region:
        if score > 0:
            converted_score = score * lambdaa
            score_total += 2 ** converted_score

	#once region score is summed, loop back through the region array and calculate
	#confidence for each hit
    for score in region:
        if score > 0:
            converted_score = score * lambdaa
            confidence = ((2 ** converted_score) / score_total)

            if confidence_string != "":
                confidence_string += f" {confidence}"
            else:
                confidence_string = f"{confidence}"
        else:
            if confidence_string != "":
                confidence_string += " 0"
            else:
                confidence_string = "0"

    return confidence_string


def FillConfidenceMatrix(row_num: int, lamb: float, columns: List[int], align_matrix: Dict[Tuple[int, int], int]) -> Dict[Tuple[int, int], float]:
	confidence_matrix: Dict[Tuple[int, int], float] = {}
	
	for i in range(len(columns)):
		if i >= len(columns):
			break
			
		col_index: int = columns[i]
		temp_region: List[int] = []
		
		for row_index in range(row_num):
			if (row_index, col_index) in align_matrix:
				temp_region.append(align_matrix[row_index, col_index])
			else:
				temp_region.append(0)
			
		temp_confidence: List[str] = ConfidenceCM(lamb, temp_region).split(" ")

		for row_index2 in range(rows):
			if temp_confidence[row_index2] != '0':
				confidence_matrix[row_index2, col_index] = float(temp_confidence[row_index2])
	
	return confidence_matrix


# Fills support score matrix using values in conf matrix
#score for subfam x at position i is sum of all confidences for subfam x for all segments that 
#overlap position i - divided by number of segments
def FillSupportMatrix(row_num: int, columns: List[int], align_matrix: Dict[Tuple[int, int], int], confidence_matrix: Dict[Tuple[int, int], float]) -> Dict[Tuple[int, int], float]:
    
    support_matrix: Dict[Tuple[int, int], float] = {}

    for row_index in range(row_num):
        tempcol: int = -1
        for col in range(len(columns)):
            if col >= len(columns):
                break
            col_index: int = columns[col]

            if (row_index, col_index) in confidence_matrix:
                num: int = col_index
                summ: float = 0.0
                num_segments: int = 0
                while num >= 0 and num >= col_index:
                    if (row_index, num) in confidence_matrix:
                        summ = summ + confidence_matrix[row_index, num]
                        num_segments += 1
                    num -= 1

                if num_segments > 0:
                    support_matrix[row_index, col_index] = summ / num_segments
        
    return support_matrix


#collapses matrices 
#collapse and combine rows that are the same subfam - just sum their support 
#new support dict has key = subfamname.col 
#also creates a bookkeeping dict that has all the cols as keys and their values 
#are arrays that hold all the active subfams in that col - used so that don't have 
#to loop through all the $i's just to see if a column exists 
def CollapseMatrices(row_num: int, columns: List[int], subfams: List[str], strands: List[str], support_matrix: Dict[Tuple[int, int], float], consensus_matrix: Dict[Tuple[int, int], int]) -> Tuple[int, Dict[Tuple[str, int], int], Dict[Tuple[str, int], str], Dict[Tuple[int, int], float], Dict[str, int], Dict[int, List[str]]]:
	
	row_num_update: int = 0
	consensus_matrix_collapse: Dict[Tuple[str, int], int] = {}
	strand_matrix_collapse: Dict[Tuple[str, int], str] = {}
	support_matrix_collapse: Dict[Tuple[str, int], float] = {}
	subfams_collapse: Dict[str, int] = {}
	active_cells_collapse: Dict[int, List[str]] = {}
	
	for col in range(len(columns)):
		col_index: int = columns[col]
		dup_max_consensus: Dict[str, float] = {}
		dup_max_support: Dict[str, float] = {}
	
		active_cols: List[str] = []
		active_cells_collapse[col_index] = active_cols
	
		#sum the support score for row_nums that are collapsed together
		#find max support score for collapsed row_nums and use the consensus from that row_num
		for row_index in range(row_num):
			if (row_index, col_index) in support_matrix and (row_index, col_index) in consensus_matrix:
				if (subfams[row_index]) in dup_max_consensus:
					if support_matrix[row_index, col_index] > dup_max_consensus[subfams[row_index]]:
						dup_max_consensus[subfams[row_index]] = support_matrix[row_index, col_index]
						consensus_matrix_collapse[subfams[row_index], col_index] = consensus_matrix[row_index, col_index]
						strand_matrix_collapse[subfams[row_index], col_index] = strands[row_index]
				else:
					dup_max_consensus[Subfams[row_index]] = support_matrix[row_index, col_index]
					consensus_matrix_collapse[subfams[row_index],col_index] = consensus_matrix[row_index, col_index]
					strand_matrix_collapse[subfams[row_index], col_index] = strands[row_index]
		
			if (row_index, col_index) in support_matrix:
				if (subfams[row_index]) in dup_max_support:
					if support_matrix[row_index, col_index] > dup_max_support[subfams[row_index]]:
						dup_max_support[Subfams[row_index]] = support_matrix[row_index, col_index]
						support_matrix_collapse[subfams[row_index],col_index] = support_matrix[row_index,col_index]
				else:
					dup_max_support[subfams[row_index]] = support_matrix[row_index,col_index]
					support_matrix_collapse[subfams[row_index],col_index] = support_matrix[row_index,col_index]
					active_cells_collapse[col_index].append(subfams[row_index])
					
	for i in range(row_num):
		subfams_collapse[subfams[i]] = 0
		
	#update var row_nums after collapse 
	row_num_update = len(subfams_collapse)

	return (row_num_update, consensus_matrix_collapse, strand_matrix_collapse, support_matrix_collapse, subfams_collapse, active_cells_collapse)


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
def FillProbabilityMatrix(same_prob_skip: float, same_prob: float, change_prob: float, change_prob_skip: float, columns: List[int], subfams_collapse: Dict[str, int], active_cells_collapse: Dict[int, List[str]], support_matrix_collapse: Dict[Tuple[str, int], float], strand_matrix_collapse: Dict[Tuple[str, int], str], consensus_matrix_collapse: Dict[Tuple[str, int], int]) -> Tuple[Dict[Tuple[str, int], float], Dict[Tuple[str, int], str]]:
## all probabilities are in log	space
## row indices are subfam names
	prob_matrix: Dict[Tuple[str, int], float] = {}
	origin_matrix: Dict[Tuple[str, int], str] = {}
	
	#fill first col of prob_matrix with 0s
	for k in subfams_collapse:
		prob_matrix[k, NonEmptyColumns[0]] = 0
		
	for columns_index in range(1, len(columns)):
	 			
		for row_index in active_cells_collapse[columns[columns_index]]:
			max: float = -inf
			max_index: str = ''
			support_log: float = log(support_matrix_collapse[row_index, columns[columns_index]])
			
			#row = prev_row_index
			#loop through all the subfams in the previous column
			for prev_row_index in active_cells_collapse[columns[columns_index-1]]:
				score: float = support_log + prob_matrix[prev_row_index, columns[columns_index-1]]
				prob: int = 0
				
				if prev_row_index == row_index:  #staying in same row
					
					if prev_row_index == 'skip':  #staying in skip
						prob = same_prob_skip
						
					prob = same_prob

					# because rows are collapsed, if the consensus seqs are contiguous - treat as if they are not the same row and get the jump penalty
					if strand_matrix_collapse[prev_row_index, columns[columns_index-1]] != strand_matrix_collapse[row_index, columns[columns_index]]:
						prob = change_prob
					else: #if on same strand
						if strand_matrix_collapse[prev_row_index, columns[columns_index-1]] == '+':
							if consensus_matrix_collapse[prev_row_index, columns[columns_index-1]] > consensus_matrix_collapse[row_index, columns[columns_index]] + 50:
								prob = change_prob
						if strand_matrix_collapse[prev_row_index, columns[columns_index-1]] == '-':
							if consensus_matrix_collapse[prev_row_index, columns[columns_index-1]] + 50 < consensus_matrix_collapse[row_index, columns[columns_index]]:
								prob = change_prob
						
				else:  #jumping rows
					prob = change_prob  
					if prev_row_index == 'skip' or row_index == 'skip':  #jumping in or out of skip
						prob = change_prob_skip
					
				score = score + prob
					
				if score > max:
					max = score
					max_index = prev_row_index
					
			prob_matrix[row_index, columns[columns_index]] = max
			origin_matrix[row_index, columns[columns_index]] = max_index
			
	return (prob_matrix, origin_matrix)


#using origin matrix, back traces through the 2D array to get the subfam path
#finds where the path switches to a different row and populates @Changes and @ChangesPosition
#reverses @Changes and @ChangesPosition because it's a backtrace so they are initially backwards
#jumps over removed columns when necessary
def GetPath(num_col: int, temp_id: int, columns: List[int], ids: List[int], subfams: List[str], active_cells_collapse: Dict[int, List[str]], prob_matrix: Dict[Tuple[str, int], float], origin_matrix: Dict[Tuple[str, int], str]) -> Tuple[int, List[int], List[str]]:
    maxxx: float = -inf
    max_row_index: str = ''
    
    changes_position: List[int] = []
    changes: List[str] = []
    
    for i in active_cells_collapse[num_col - 1]:
        if maxxx < prob_matrix[i, num_col - 1]:
            maxxx = prob_matrix[i, num_col - 1]
            max_row_index = i
    
    prev_row_index: str = origin_matrix[max_row_index, num_col - 1]
    
    changes_position.append(len(columns))
    
    #already added the last col, but this adds the one before $col so still start at last col
    for columns_index in range(len(columns)-1, 0, -1):
    	if (prev_row_index, columns[columns_index-1]) in origin_matrix:
    	
    		ids[columns[columns_index-1]] = temp_id
    		
    		if prev_row_index != origin_matrix[prev_row_index, columns[columns_index-1]]:
    			temp_id+=1234
    			changes_position.append(columns_index-1)
    			changes.append(prev_row_index)
    			
    		prev_row_index = origin_matrix[prev_row_index, columns[columns_index-1]]
    		
    ids[columns[0]] = temp_id
    changes_position.append(0)
    changes.append(prev_row_index)
    
    changes.reverse()
    changes_position.reverse()
    
    #changes ID for next round of stitching, so when starts stitching will have unique ID
    temp_id+=1234
    
    return (temp_id, changes_position, changes)
    

#for debugging
def PrintChanges(columns: List[int], changes: List[str], changes_position: List[int]) -> None:
    i: int = 0
    while i < len(changes):
        stdout.write(str(columns[changes_position[i]]))
        stdout.write("\t")
        stdout.write(f"{changes[i]}\n")
        i=i+1


#FIXME - if goes into skip state for just one position, this will have an error .. also when
#stitching want to ignore skip states
#fills node confidence matrix 
#first fills matrix with node alignment scores, then reuses matrix for confidence scores 
def FillNodeConfidence(nodes: int, gap_ext: int, gap_init: int, sub_matrix_cols: int, lamb: float, columns: List[int], subfam_seqs: List[str], chrom_seqs: List[str], changes_position: List[int], subfams: List[str], sub_matrix: Dict[int, int], char_pos: Dict[str, int]) -> Dict[Tuple[str, int], float]:   
    node_confidence_temp: List[float] = [0 for _ in range(len(subfams) * nodes)]
    
    node_confidence: Dict[Tuple[str, int], float] = {}
        
    #calculated node confidence for first node - doesn't look back at prev char bc there isn't one
    #starts at 1 bc don't need to calc for skip state
    for subfam_index in range(1, len(subfams)):
        begin_node: int = columns[changes_position[0]]
        end_node: int = columns[changes_position[1]]
        subfam: str = subfam_seqs[subfam_index][begin_node:end_node]
        chrom: str = chrom_seqs[subfam_index][begin_node:end_node]
        align_score: float = CalcScore(gap_ext, gap_init, sub_matrix_cols, subfam, chrom, '', '', sub_matrix, char_pos)
        node_confidence_temp[subfam_index * nodes + 0] = align_score
	 
	#does rest of nodes - looks back at prev char incase of gap ext
    for node_index2 in range(1, nodes-1):
    	for subfam_index2 in range(1, len(subfams)):
            begin_node2: int = columns[changes_position[node_index2]]
            end_node2: int = columns[changes_position[node_index2 + 1]]
            subfam: str = subfam_seqs[subfam_index2][begin_node2:end_node2]
            chrom: str = chrom_seqs[subfam_index2][begin_node2:end_node2]
            last_prev_subfam2: str = subfam_seqs[subfam_index2][begin_node2-1]#subfam_seqs[j][NonEmptyColumns[changes_position[i + 1] - 1]]
            last_prev_chrom2: str = chrom_seqs[subfam_index2][begin_node2-1]#chrom_seqs[j][changes_position[i + 1] - 1]
            align_score: float = CalcScore(gap_ext, gap_init, sub_matrix_cols, subfam, chrom, last_prev_subfam2, last_prev_chrom2, sub_matrix, char_pos)
            node_confidence_temp[subfam_index2 * nodes + node_index2] = align_score
    
	#does last node
    for node_index3 in range(1, len(subfams)):
        begin_node3: int = columns[changes_position[-2]]
        end_node3: int = columns[changes_position[-1]-1]
        subfam: str = subfam_seqs[node_index3][begin_node3:end_node3]
        chrom: str = chrom_seqs[node_index3][begin_node3:end_node3]
        last_prev_subfam3: str = subfam_seqs[node_index3][begin_node3-1]
        last_prev_chrom3: str = chrom_seqs[node_index3][begin_node3-1]
        align_score: float = CalcScore(gap_ext, gap_init, sub_matrix_cols, subfam, chrom, last_prev_subfam3, last_prev_chrom3, sub_matrix, char_pos)
        node_confidence_temp[node_index3 * nodes + nodes-1] = align_score
 
    #reuse same matrix and compute confidence scores for the nodes	
    for node_index4 in range(nodes):
    	temp: List[float] = []
    	for row_index in range(len(subfams)):
    		temp.append(node_confidence_temp[row_index * nodes + node_index4])
    		
    	confidence_temp: List[float] = [float(x) for x in ConfidenceCM(lamb, temp).split(' ')]
    	
    	for row_index2 in range(len(subfams)):
    		node_confidence_temp[row_index2 * nodes + node_index4] = confidence_temp[row_index2]

	#collapse node_confidence down same way supportmatrix is collapsed - all seqs of 
	#the same subfam are put in the same row
	#not a sparse hash - holds the 0s, but I think this is okay because it won't ever
	#be a very large matrix, and this way we don't have to test if anything exists 	

    for node_index5 in range(nodes):
        for row_index3 in range(len(subfams)):
            if (subfams[row_index3], node_index5) in node_confidence:
                node_confidence[subfams[row_index3], node_index5] += node_confidence_temp[row_index3 * nodes + node_index5]
            else:
            	node_confidence[subfams[row_index3], node_index5] = node_confidence_temp[row_index3 * nodes + node_index5]
            	
    return node_confidence


#used for debugging
def PrintPathGraph(nodes: int, changes: List[str], path_graph: List[int]) -> None:
	stdout.write(" ")
	for i in range(nodes):
		stdout.write(f"{changes[i]} ")
	stdout.write("\n")
	
	for i in range(nodes):
		for j in range(nodes):
			stdout.write(f"{path_graph[i*nodes+j]}\t")
		stdout.write(f"{changes[i]}\n")
	stdout.write("\n")


#used for debugging
def PrintNodeConfidence(nodes: int, changes: List[str], subfams_collapse: Dict[str, int], node_confidence: Dict[Tuple[str, int], float]) -> None:
	for i in range(nodes):
		stdout.write(f"{changes[i]} ")
	stdout.write("\n")

	for subfam in subfams_collapse:
		stdout.write(f"{subfam} ")
		for j in range(nodes):
			if (subfam,j) in node_confidence:
				stdout.write(f"{node_confidence[subfam,j]} ")
			else:
				stdout.write(f"-inf ")
		stdout.write("\n")


def FillPathGraph(nodes: int, columns: List[int], changes: List[str], changes_position: List[int], subfams_collapse: Dict[str, int], consensus_matrix_collapse: Dict[Tuple[str, int], int], strand_matrix_collapse: Dict[Tuple[str, int], str], node_confidence: Dict[Tuple[str, int], float]) -> List[int]:
	
    path_graph: List[int] = []
    
    for i in range(nodes * nodes):
        path_graph.append(0)

	# filling beginning path graph with straight line through the nodes
    for i in range(nodes - 1):
        path_graph[i * nodes + i + 1] = 1
        
        
    for sink_node_index in range(nodes):
        sink_subfam: str = str(changes[sink_node_index])
        sink_subfam_start: int = consensus_matrix_collapse[sink_subfam, columns[changes_position[sink_node_index]]]
        sink_strand: str = strand_matrix_collapse[sink_subfam, columns[changes_position[sink_node_index]]]

		#looks at all the preceding nodes, except the one directly before it (source nodes)
        for source_node_index in range(sink_node_index - 1):
        	#look at all the subfams in each node
            for source_subfam in subfams_collapse:
            	source_subfam = str(source_subfam)
            	sourceConf = node_confidence[source_subfam, source_node_index]
            	
            	if (source_subfam, columns[changes_position[source_node_index + 1]-1]) in consensus_matrix_collapse:
            		source_subfam_stop = consensus_matrix_collapse[source_subfam, columns[changes_position[source_node_index + 1]-1]]
            		source_strand = strand_matrix_collapse[source_subfam, columns[changes_position[source_node_index + 1]-1]]
  
  					# adds in edge if the subfam of the sink is at the source node and if it's 
					# confidence >= 1%, and if the source is before the sink in the consensus sequence           		
            		if sink_strand == '+' and sink_strand == source_strand:
            			if (sink_subfam == source_subfam) and (sourceConf >= 0.3):
            				#FIXME- not sure what this overlap should be .. just allowed 50 for now            				
            				if source_subfam_stop <= sink_subfam_start + 50:
            					path_graph[source_node_index * nodes + sink_node_index] = 1
            					
            		elif sink_strand == '-' and sink_strand == source_strand:
            			if sink_subfam == source_subfam and sourceConf >= 0.3:
            				if source_subfam_stop >= sink_subfam_start + 50:
            					path_graph[source_node_index * nodes + sink_node_index] = 1
            					
    return path_graph


#finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds them to
# @RemoveStarts and @RemoveStops so they can be extracted from the alignment 
#pass in columns and updates it to ignore removed nodes
# updates cols if removes node from end 
def ExtractNodes(num_col: int, nodes: int, num_columns: List[int], changes_position: List[int], path_graph: List[int]) -> int:
	
	remove_starts: List[int] = []
	remove_stops: List[int] = []
	
	updated_num_col: int = num_col

	#boolean for which nodes will be removed
	remove_nodes: List[bool] = [False for _ in range(nodes)]
	
	#extracting nodes that only have one incoming and one outgoing edge
	num_edges_in: List[int] = [0 for _ in range(nodes)]
	num_edges_out: List[int] = [0 for _ in range(nodes)]
	
	for row in range(nodes):
		for col in range(nodes):
			num_edges_in[col] += path_graph[row * nodes + col]
			num_edges_out[row] += path_graph[row * nodes + col]
			
	for node in range(nodes - 1):
		if num_edges_in[node] <= 1 and num_edges_out[node] <= 1:
			remove_starts.append(changes_position[node])
			remove_stops.append(changes_position[node + 1])
			remove_nodes[node] = True
	
	#deals with last node, so when $NumNodes-1 the last remove stop is the end of the matrix
	if num_edges_in[nodes - 1] <= 1 and num_edges_out[nodes - 1] <= 1:
		remove_starts.append(changes_position[nodes - 1])
		remove_stops.append(len(num_columns)-1)
		remove_nodes[nodes - 1] = True
		
	#when removing from the end, have to update num_cols because don't want to do to the end of the matrix anymore 
	col_index: int = nodes - 1
	while remove_nodes[col_index]:
		updated_num_col = num_columns[changes_position[col_index] - 1]  #FIXME - do I index num_cols here or not?
		col_index -= 1
		
	# removing inserted elements from @NonEmptyColumns so they can be ignored 
	total: int = 0
	for i in range(len(remove_stops)):
		del num_columns[remove_starts[i]-total:remove_stops[i]-total]
		# 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
		total += (remove_stops[i] - remove_starts[i])

	return updated_num_col


#uses position in matrix
def PrintResults(changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(changes_orig)):
		if str(changes_orig[i]) != 'skip':
			stdout.write(str(columns_orig[changespos_orig[i]]))
			stdout.write("\t")
			stdout.write(str(columns_orig[changespos_orig[i+1]-1]))
			stdout.write("\t")
			stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
			stdout.write("\t")
			stdout.write(str(changes_orig[i]))
			stdout.write("\n")

	
#uses position in input sequence
def PrintResultsSequence(edgestart: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
	stdout.write("start\tstop\tID\tname\n")
	stdout.write("----------------------------------------\n")
	for i in range(len(changes_orig)):
		if str(changes_orig[i]) != 'skip':
			stdout.write(str(columns_orig[changespos_orig[i]]+edgestart))
			stdout.write("\t")
			stdout.write(str(columns_orig[changespos_orig[i+1]-1]+edgestart))
			stdout.write("\t")
			stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
			stdout.write("\t")
			stdout.write(str(changes_orig[i]))
			stdout.write("\n")


#---------------------------------------------------------------------------------------#
#GLOBALS - yes I know these are bad, will be passed into functions to make port easier	#
#---------------------------------------------------------------------------------------#

GapInit: int = -25
GapExt: int = -5
Lamb: float = 0.1227  # From command line
ChunkSize: int = 31
SameProbLog: float = log(1 - (10 ** -45))  #FIXME - this just becomes 0.0... need to be more precise
ChangeProb: float = 10 ** -45
ChangeProbLog: float = 0.0  # Reassigned later
ChangeProbSkip: float = 0.0 # Reassigned later
SameProbSkip: float = 0.0
SkipAlignScore: int = 30 #FIXME - still need to decide what this number is, skip state doesn't work in seqs_fullAlu.align unless SkipAlignScore = 120
StartAll: int = 0  # Reassigned later
StopAll: int = 0  # Reassigned later
ID: int = 1111

help: bool = False  # Reassigned later
prin: bool = False  # Reassigned later
printMatrixPos: bool = False  # Reassigned later

helpMessage: str = f"""
usage: {argv[0]} alignFile matrixFile\n
ARGUMENTS
    --GapInit[-25]
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
    "GapInit=",
    "GapExt=",
    "lambda=",
    "segmentsize=",
    "changeprob=",

    "help",
    "matrixPos",
])
opts = dict(raw_opts)

GapInit = int(opts["--GapInit"]) if "--GapInit" in opts else GapInit
GapExt = int(opts["--GapExt"]) if "--GapExt" in opts else GapExt
Lamb = float(opts["--lambda"]) if "--lambda" in opts else Lamb
ChunkSize = int(opts["--segmentsize"]) if "--segmentsize" in opts else ChunkSize
ChangeProb = float(opts["--changeprob"]) if "--changeprob" in opts else ChangeProb
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
SubMatrixCols: int = len(chars)


#reads in the score matrix from file and stores in 2D array Score matrix 
SubMatrix: Dict[int, int] = {}
count: int = 0
for line in in_matrix[1:]:
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    subScores = re.split(r"\s+", line)
    for i in range(len(subScores)):
        SubMatrix[count*SubMatrixCols+i] = int(subScores[i])
    count += 1

Subfams: List[str] = []
Scores: List[int] = []
Strands: List[str] = []
Starts: List[int] = []
Stops: List[int] = []
ConsensusStarts: List[int] = []
ConsensusStops: List[int] = []
SubfamSeqs: List[str] = []
ChromSeqs: List[str] = []

AlignMatrix: Dict[Tuple[int, int], int] = {}
ConfidenceMatrix: Dict[Tuple[int, int], float] = {}
SupportMatrix: Dict[Tuple[int, int], float] = {}
ProbMatrix: Dict[Tuple[str, int], float] = {}
OriginMatrix: Dict[Tuple[str, int], str] = {}
ConsensusMatrix: Dict[Tuple[int, int], int] = {}

NonEmptyColumns: List[int] = [];

Changes: List[str] = []
ChangesPosition: List[int] = []

IDs: List[int] = [];
ChangesOrig: List[int] = [];
ChangesPositionOrig: List[int] = [];
NonEmptyColumnsOrig: List[int] = [];

SupportMatrixCollapse: Dict[Tuple[str, int], int] = {}
ActiveCellsCollapse: Dict[int, List[str]] = {}
SubfamsCollapse: Dict[str, int] = {}
ConsensusMatrixCollapse: Dict[Tuple[str, int], int] = {}
StrandMatrixCollapse: Dict[Tuple[str, int], str] = {}

#for graph/node part
NumNodes: int = 0
NodeConfidence: Dict[Tuple[str, int], float] = {}
PathGraph: List[int] = []
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

        Subfams.append(alignment.subfamily)
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
		stdout.write(f"{0}\t{Stops[1]-Starts[1]}\t1111\t{Subfams[1]}\n")
	else:
		stdout.write("start\tstop\tID\tname\n")
		stdout.write("----------------------------------------\n")
		stdout.write(f"{Starts[1]}\t{Stops[1]}\t1111\t{Subfams[1]}\n")
	exit()

ChangeProbLog = log(ChangeProb / (numseqs - 1))
ChangeProbSkip = ChangeProbLog / 2;
SameProbSkip = ChangeProbLog / 10; # 10% of the jump penalty, staying in skip state for 20nt "counts" as one jump

#precomputes global vars rows and cols in matrices 
rows: int = len(Subfams)
cols: int = 0 #assign cols in FillAlignMatrix

(StartAll, Stopall) = PadSeqs(Starts, Stops, SubfamSeqs, ChromSeqs)
                    
(cols, AlignMatrix) = FillAlignMatrix(StartAll, ChunkSize, GapExt, GapInit, SkipAlignScore, SubMatrixCols, SubfamSeqs, ChromSeqs, Starts, SubMatrix, CharPos)

ConsensusMatrix = FillConsensusPositionMatrix(cols, SubfamSeqs, ChromSeqs, ConsensusStarts, ConsensusStops, Strands)

NonEmptyColumns = FillColumns(cols, rows, AlignMatrix)

ConfidenceMatrix = FillConfidenceMatrix(rows, Lamb, NonEmptyColumns, AlignMatrix)

SupportMatrix = FillSupportMatrix(rows, NonEmptyColumns, AlignMatrix, ConfidenceMatrix)

(rows, ConsensusMatrixCollapse, StrandMatrixCollapse, SupportMatrixCollapse, SubfamsCollapse, ActiveCellsCollapse) = CollapseMatrices(rows, NonEmptyColumns, Subfams, Strands, SupportMatrix, ConsensusMatrix)


(ProbMatrix, OriginMatrix) = FillProbabilityMatrix(SameProbSkip, SameProbLog, ChangeProbLog, ChangeProbSkip, NonEmptyColumns, SubfamsCollapse, ActiveCellsCollapse, SupportMatrixCollapse, StrandMatrixCollapse, ConsensusMatrixCollapse)

IDs = [0] * cols

(ID, ChangesPosition, Changes) = GetPath(cols, ID, NonEmptyColumns, IDs, Subfams, ActiveCellsCollapse, ProbMatrix, OriginMatrix)

#keep the original annotation for reporting results
ChangesOrig = Changes.copy()
ChangesPositionOrig = ChangesPosition.copy()
NonEmptyColumnsOrig = NonEmptyColumns.copy()			

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
    NumNodes = len(Changes)
    
    #breakout of loop if there are 2 or less nodes left
    if (NumNodes <= 2):
        break 
	
    NodeConfidence.clear()
    
    #initializes and fills node confidence matrix
    NodeConfidence = FillNodeConfidence(NumNodes, GapExt, GapInit, SubMatrixCols, Lamb, NonEmptyColumns, SubfamSeqs, ChromSeqs, ChangesPosition, Subfams, SubMatrix, CharPos)
        
    PathGraph.clear()
    PathGraph = FillPathGraph(NumNodes, NonEmptyColumns, Changes, ChangesPosition, SubfamsCollapse, ConsensusMatrixCollapse, StrandMatrixCollapse, NodeConfidence)
    
	#test to see if there nodes in the graph that have more that one incoming or outgoing edge,
	#if so keep looping, if not break out of the loop
	#if they are all 0, break out of the loop
    test: bool = False
    j: int = 0
    while j < NumNodes:
        i: int = 0
        while i < j - 1:
            if (PathGraph[i * NumNodes + j] == 1):
                test = True
            i += 1
        j += 1

    if not test:
        break

    cols = ExtractNodes(cols, NumNodes, NonEmptyColumns, ChangesPosition, PathGraph)
         
    # using prob matrix and origin matrix, just skip the cols I'm not interested in and annotate
	# without the removed subfam
	# using old prob matrix and origin matrix        
	# this time ignores inserted subfam because there are values in @RemoveStarts and @RemoveStops
    (ProbMatrix, OriginMatrix) = FillProbabilityMatrix(SameProbSkip, SameProbLog, ChangeProbLog, ChangeProbSkip, NonEmptyColumns, SubfamsCollapse, ActiveCellsCollapse, SupportMatrixCollapse, StrandMatrixCollapse, ConsensusMatrixCollapse)

    Changes.clear()
    ChangesPosition.clear()
        
    (ID, ChangesPosition, Changes) = GetPath(cols, ID, NonEmptyColumns, IDs, Subfams, ActiveCellsCollapse, ProbMatrix, OriginMatrix)


if printMatrixPos:
	PrintResults(ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)
else:
	PrintResultsSequence(StartAll, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)


        


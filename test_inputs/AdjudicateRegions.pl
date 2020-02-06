#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use Getopt::Long;
use POSIX ();
use Time::HiRes;
    
#------------------------------------------------------------------------------------
#test client -
# 	Usage: 
# 		usage: $0 alignFile matrix
# 		
# 	Input:
# 		-align file of region
# 		-Score matrix file
# 	
# 	-output:
# 		positon in alignment, subfamily switch
#------------------------------------------------------------------------------------


# my $nextafter = POSIX::nextafter(0.0, .01);
# POSIX::nexttoward(1.5, 1.6); - supposed to be more accurate
# can't get these to work without updating to perl v 30
# ignoring this for now because no point in fixing this in perl when going to rewrite in
# a different language 


my $gapInit = -25;
my $gapExt = -5;
my $lambda;
my $chunksize = 30;
my $sameProbLog = log(1-(10**-45)); #FIXME - this just becomes 0, because sameProb becomes 1... 
my $changeProb = (10**-45);
my $changeProbLog;
my $skipAlignScore = 30; #FIXME - have to figure out what the best score for this is
my $startall;
my $stopall;


my $help;
my $print;
my $printMatrixPos;


my $helpMessage = 
"\nusage: $0 alignFile matrixFile\n
ARGUMENTS
	--gapInit [-25]
	--getExt [-5]
	--lambda [will calc from matrix if not included]
	--segmentsize [30]
	--changeprob [1e-45]
	
OPTIONS
	--help - display help message
	--printmatrices - output all matrices
	--matrixpos - prints subfam changes in matrix position instead of genomic position
";
 
GetOptions (
	'gapInit=i' => \$gapInit,
    'gapExt=i' => \$gapExt,
    'lambda=f' => \$lambda,
    'segmentsize=i' => \$chunksize,
    'changeprob=f' => \$changeProb,
           
    'help' => \$help,
    'printmatrices' => \$print,
    'matrixPos' => \$printMatrixPos,
);
           
if($help){
	die $helpMessage;
}


#input is alignment file of hits region and substitution matrix 
my $infile = shift;
my $infile_matrix = shift;

open(my $in,  "<",  $infile)
        or die "Could not open $infile";
open(my $in_matrix,  "<",  $infile_matrix)
        or die "Could not open $infile_matrix";
                
        
#calculates lambda using esl-scorematrix if not included at command line         
if(!$lambda){	
	my $esloutput  = `/Users/kaitlincarey/git/hmmer/easel/esl_scorematrix --dna $infile_matrix`;
		
	if($esloutput =~ /.+?\nLambda\s+=\s([0-9.]+)\n[\s\S]+/g){
		$lambda = $1;
	}
}


#creates a dict that associates character from score matrix file with position in score matrix
#alignment char as key and pos in array as value   
my %CharPos;        
# <$in_matrix>;
my $line = <$in_matrix>;
$line =~ s/^\s+//;
my @chars = split(/\s+/, $line);
for(my $i = 0; $i < @chars; $i++){
	$CharPos{$chars[$i]} = $i;
}
my $subMatrixCols = scalar @chars;


#reads in the score matrix from file and stores in 2D array Score matrix 
my @SubMatrix;
my $count = 0;
while($line = <$in_matrix>){
	$line =~ s/^\s+//;
	my @subScores = split(/\s+/, $line);
	for(my $col = 0; $col < @subScores; $col++){
		$SubMatrix[$count*$subMatrixCols+$col] = $subScores[$col];
	}
	$count++;
}       
close $in_matrix;

   
my @Subfams; 
my @Scores;
my @Strands;
my @Starts;
my @Stops;
my @ConsensusStarts;
my @ConsensusStops;
my @SubfamSeqs;
my @ChromSeqs;

my %AlignHash;
my %ConfHash;
my %SupportHash;
my %ProbHash;
my %OriginHash;
my %ConsensusHash;

my @subfampath;
my @Position;

my @RemoveStarts;
my @RemoveStops;

my @Changes;
my @ChangesPos;


push @Subfams, 'skip';

push @Scores, '';
push @Starts, '';
push @Strands, '';
push @Stops, '';
push @ConsensusStarts, '';
push @ConsensusStops, '';
push @SubfamSeqs, '';
push @ChromSeqs, '';


$/ = "Align: ";
  
#opens alignment file and stores all Subfams, Scores, Starts, Stops, subfams seqs and Chrom seqs in arrays 
my $numseqs = 0;  
while(my $region = <$in>){
	$region =~ s/Align: //;
 	while($region =~ /(.+?)\t.+?\t(\d+)\t([+-C])\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t*\n>.+?\n([\s\S]+)\n>.+?\n([\s\S]+)/g){
		
		$numseqs++;
		my $subfam = $1;
		my $score = $2;
		my $strand = $3;
		my $start = $4;
		my $stop = $5;
		my $consensusStart = $6;
		my $consensusStop = $7;
						
		my $chromSeq = $8;
		my $subfamSeq = $9;
		$chromSeq =~ s/[\n]//g;
		$subfamSeq =~ s/[\n]//g;
						
 		push @Subfams, $subfam;
 		push @Scores, $score;
 		push @Strands, $strand;
 		push @Starts, $start;
 		push @Stops, $stop;
 		push @ConsensusStarts, $consensusStart;
 		push @ConsensusStops, $consensusStop;
		push @SubfamSeqs, $subfamSeq;
		push @ChromSeqs, $chromSeq;	
	}		   
}
close $in;


$changeProbLog = log($changeProb/($numseqs-1));

#arrays for which subfams are active, meaning they have been identified as the subfam for
#part of the seq - and the IDs and start and stop positions in the alignment of the 
# active subfams
my @IDs;
my @Active;
my @nodeStarts;
my @nodeStops;
my $idnum = 1111; #makes the IDs 4 digit nums
for(my $i = 0; $i < @Subfams; $i++){
	push @IDs, $idnum;
	push @Active, 0;
	push @nodeStarts, 0;
	push @nodeStops, 0;
	
	$idnum++;
}


#precomputing rows and cols in align score matrix
my $rows = @Subfams;  
my $cols = 0;  #assign cols in FillAlignScoreMatrix

padSeqs(\@Starts, \@Stops, \@SubfamSeqs, \@ChromSeqs);

@Position = (0) x @SubfamSeqs;
FillAlignScoreMatrix(\@SubfamSeqs, \@ChromSeqs); 

FillConsensusPosMatrix(\%ConsensusHash, \@SubfamSeqs, \@ChromSeqs, \@ConsensusStarts, \@ConsensusStops);

 
#BUG HERE - FIXME - if you hard code in here and fill those gaps, the issue is resolved
# for Seqs_doubleNested.align
# $AlignHash{'129.969'} = $AlignHash{'129.970'};
# for Seqs_doubleNesting2.align
# $AlignHash{'129.1498'} = $AlignHash{'129.1497'};


#puts all columns that are not empty into @Columns, so when I loop through hash I can use the 
#vals in @Columns - this will skip over empty columns
my @Columns;
my $j;
for($j = 0; $j < $cols; $j++){

	my $empty = 1;
	for(my $i = 0; $i < $rows; $i++){
		if(exists $AlignHash{$i.'.'.$j}){
			$empty = 0;
			$i = $rows;
		}
	}
 	 	
	if(!$empty){
		push @Columns, $j;
	}
	
	#FIXME - instead of adding a skip state score, just give the skip state a raw align score and 
	#let the algorithm calculate it's skip score 
	$AlignHash{'0.'.$j} = $skipAlignScore;
}

#adding the last $chunksize-1 columns that are not in the AlignScoreMatrix to @Columns
for(my $i = 0; $i < $chunksize-1; $i++){
	push @Columns, $j;
	$j++;
}


#printing columns array to take as input into python version of FillProbMatrix
# for(my $i = 0; $i < @Columns; $i++){
# 	print "$Columns[$i]\n";
# }


FillConfScoreMatrix(\%AlignHash, \%ConfHash); 

$cols = $cols + $chunksize-1;

FillSupportMatrix(\%SupportHash, \%AlignHash, \%ConfHash);


#collapse and combine rows that are the same subfam - just sum their support 
#new support dict has key = subfamname.col 
#also creates a bookkeeping dict that has all the cols as keys and their values 
#are arrays that hold all the active subfams in that col - used so that don't have 
#to loop through all the $i's just to see if a column exists 
my %SupportHash_collapse;
my %ActiveCells_collapse;
my %Subfams_collapse;
my %ConsensusHash_collapse;
my %StrandHash_collapse;

for(my $j = 0; $j < $cols; $j++){

	my @activecols;
	$ActiveCells_collapse{$j} = [@activecols];
	
	my $max = 0;
	my $maxrow;
	for(my $i = 1; $i < $rows; $i++){
		if(exists $ConsensusHash{$i.'.'.$j}){	
			
			$ConsensusHash_collapse{$Subfams[$i].'.'.$j} = $ConsensusHash{$i.'.'.$j};
			$StrandHash_collapse{$Subfams[$i].'.'.$j} = $Strands[$i];
			
			#sum the support score for rows that are collapsed together
			#find max support score for collapsed rows and use the consensus from that row
			if(exists $SupportHash_collapse{$Subfams[$i].'.'.$j}){
				$SupportHash_collapse{$Subfams[$i].'.'.$j} = $SupportHash_collapse{$Subfams[$i].'.'.$j} + $SupportHash{$i.'.'.$j};
			
				#if it already exists check if this support score is larger than the max
				#if it us update the max values for the collapsed array 
				if($SupportHash{$i.'.'.$j} > $max){
					
					#just want to see if the support scores we are choosing from are
					#close or one clearly stands out 
 					#print "$max\t$SupportHash{$i.'.'.$j}\n";
					
					$max = $SupportHash{$i.'.'.$j};
					$maxrow = $i;
				}
				
			}else{
				$SupportHash_collapse{$Subfams[$i].'.'.$j} = $SupportHash{$i.'.'.$j};
				push @{$ActiveCells_collapse{$j}}, $Subfams[$i];
				
				#if it doesn't already exist this is the max that will go into the collapsed array
				$max = $SupportHash{$i.'.'.$j};
				$maxrow = $i;
			}
			
			
		}
	}
	
	#if $maxrow doesn't exist there were no active rows in that col 
	if($maxrow){
		#for collapsed rows inserted updated consensus pos for whichever one had the max support 
		$ConsensusHash_collapse{$Subfams[$maxrow].'.'.$j} = $ConsensusHash{$maxrow.'.'.$j};
		$StrandHash_collapse{$Subfams[$maxrow].'.'.$j} = $Strands[$maxrow];
	}
	
}

for(my $i = 0; $i < $rows; $i++){
	$Subfams_collapse{$Subfams[$i]} = 0;
}

# PrintMatrixHash($cols, %SupportHash);
# PrintMatrixHash($cols, %ConsensusHash);


#update number of rows after collapse 
$rows = scalar keys %Subfams_collapse;

# PrintMatrixHashCollapse($cols, %ConsensusHash_collapse);
# PrintMatrixHashCollapse($cols, %StrandHash_collapse);


	
# 	
# # prints out the book keeping hash that holds columns as keys and an array of all 
# #the active subfams for that row as values 	
# print STDERR "\n";
# foreach my $key (keys %ActiveCells_collapse){
# 	print STDERR "$key\t";
# 	my @value_array = @{$ActiveCells_collapse{$key}};
# 	
# 	for(my $i = 0; $i < @value_array; $i++){
# 		print STDERR "$value_array[$i] ";
# 	}
# 	print STDERR "\n";
# }


 
#fill first col of probhash with 0s
# for(my $i = 0; $i < $rows; $i++){
for my $i (keys %Subfams_collapse){
	$ProbHash{$i.'.0'} = 0;
}
 
FillProbMatrix(\%ProbHash, \%SupportHash_collapse, \%OriginHash);

# PrintMatrixHashCollapse($cols, %OriginHash);


@subfampath = GetPath(\%ProbHash, \%OriginHash, \@Subfams);
GetChanges(\@Changes, \@ChangesPos);
# GetBounds(\@Changes, \@ChangesPos);  #FIXME - this no longer works with the new hashes
PrintChanges(\@Changes, \@ChangesPos);
print STDERR "\n";

if($print){
	PrintAllMatrices();
}





#Steps- 
#1.create confidence for nodes
#	will be in a matrix that is #subfams x #nodes
#2.create path graph
#3.find alternative paths through graph and add those edges
#4.extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
#5.annotate again with removed subfams
#   --stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left

my $numnodes;
my %NodeConfidence;
my @pathGraph;
my $total;
my $loop = 1;


$count = 0;
while(1){
	$count++;
	
	#FIXME - take this out.. just does loop once for debugging	
	# if($count == 2){
# 		last;
# 	}

	$numnodes = @Changes;
	
		
	#breakout of loop if there are 2 or less nodes left
	if($numnodes <= 2){
		last;
	}
	
	
	undef %NodeConfidence;
	#initializes and fills node confidence matrix 
	NodeConfidence(\%NodeConfidence, \@SubfamSeqs, \@ChromSeqs, \@ChangesPos);

	undef @pathGraph;
	FillPathGraph(\@pathGraph);
		
	#test to see if there nodes in the graph that have more that one incoming or outgoing edge,
	#if so keep looping, if not break out of the loop
	my $test = 0;
	for(my $j = 0; $j < $numnodes; $j++){ 
		for(my $i = 0; $i < $j-1; $i++){
			#if they are all 0, break out of the loop
			if($pathGraph[$i*$numnodes+$j] == 1){
				$test = 1;
			}
		}
	}
	
	if($test == 0){
# 		print STDERR "breaking out\n";
		last;
	}
	
	undef @RemoveStarts;
	undef @RemoveStops;
	
	ExtractNodes(\@RemoveStarts, \@RemoveStops, \@ChangesPos, \@pathGraph, $numnodes);

	# removing inserted Alus from @Columns so they can be ignored 
	$total = 0;
	for(my $i = 0; $i < @RemoveStops; $i++){
		splice @Columns, $RemoveStarts[$i]-$total, ($RemoveStops[$i] - $RemoveStarts[$i]);
	
	# 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
		$total+=($RemoveStops[$i] - $RemoveStarts[$i]); 
	}
	
	# using prob matrix and origin matrix, just skip the cols I'm not interested in and annotate
	# without the removed subfam (Alu in this case)
	# using old prob matrix and origin matrix        
	# this time ignores Alus because there are values in @RemoveStarts and @RemoveStops
	FillProbMatrix(\%ProbHash, \%SupportHash_collapse, \%OriginHash);

	undef @Changes;
	undef @ChangesPos;
	undef @subfampath;
	
	@subfampath = GetPath(\%ProbHash, \%OriginHash, \@Subfams);
	GetChanges(\@Changes, \@ChangesPos);
	PrintChanges(\@Changes, \@ChangesPos);
	print STDERR "\n";	
		
}

# PrintResults();



#----------------------------------------------------------------------------------#
# 					SUBROUTINES				  			 						   #
#----------------------------------------------------------------------------------#


#initializes 2D array with specified number of rows and cols
sub InitializeMatrix{
    my ($Matrix, $row, $col) = (@_);
     #initialize matrix with all 0s  
    
    for(my $i = 0; $i < $row*$col; $i++){
    	@$Matrix[$i] = 0;
    }
    
}


#prints 2D array with subfams as row labels 
sub PrintMatrix{
    my ($Matrix) = (@_);
    
    print "\t";
    for(my $j = 0; $j <= $cols; $j++){
    	print "$j\t";
    }
    print "\n";
     
    for(my $i = 0; $i < $rows; $i++){
        printf "%-20s\t", $Subfams[$i];
        for(my $j = 0; $j < $cols; $j++){
#             printf "%5.3g", @$Matrix[$i*$cols+$j];
            print @$Matrix[$i*$cols+$j];
            print "\t";
        }
        print "\n";
    }
}


sub PrintMatrixHash{
	my ($column, %Hash) = (@_);
		 
	print "\t";
    for(my $j = 0; $j < $column; $j++){
    	print "$j\t";
    }
    print "\n";
	 
	for(my $i = 0; $i < $rows; $i++){
		printf "%-20s\t", $Subfams[$i];
		for(my $j = 0; $j < $column; $j++){
			if(exists $Hash{$i.'.'.$j}){
# 				printf "%5.3g", $Hash{$i.'.'.$j};
				print $Hash{$i.'.'.$j};
			}else{
				printf "%5.3g", '-inf';
			}
			print "\t";
		}
		print "\n";
	}
}


#print collapsed matrix - subroutine PrintMatrixHash() doesn't work because it uses
# $i and $j to test if the cell exists - with this need to use subfam name to check
sub PrintMatrixHashCollapse{
	my ($column, %Hash) = (@_);
	
	print "\t";
	for(my $j = 0; $j < $column; $j++){
    	print "$j\t";
	}
	print "\n";
	
	#%Subfams_collapse is not ordered, so prints the skip state first 
	printf "%-20s\t", "skip";
	for(my $j = 0; $j < $column; $j++){
		if(exists $Hash{'skip.'.$j}){
			print $Hash{'skip.'.$j};
		}else{
			printf "%5.3g", '-inf';
		}
		print "\t";
	}
	print "\n";
	
	
	for my $i (keys %Subfams_collapse){
		if($i ne 'skip'){
			printf "%-20s\t", $i;
			for(my $j = 0; $j < $column; $j++){
				if(exists $Hash{$i.'.'.$j}){
					print $Hash{$i.'.'.$j};
				}else{
					printf "%5.3g", '-inf';
				}
				print "\t";
			}
			print "\n";
		}
	}
}



#prints all the matrices used in the dynamic programming pass
sub PrintAllMatrices{
	print "Align Scores\n";
	PrintMatrixHash($cols, %AlignHash);
	print "confidence\n";
	PrintMatrixHash($cols, %ConfHash);
	print "support\n";
	PrintMatrixHash($cols, %SupportHash);
	print "prob\n";
	PrintMatrixHash($cols, %ProbHash);
	print "origin\n";
	PrintMatrixHash($cols, %OriginHash);
}


#prints out all the subfams that are indentified, where they start and stop in the seq, 
#and thier IDs - subfams that are stitched together have same ID
sub PrintResults{

	print STDERR "\nstart\tstop\tid\tsubfam\n";
	print STDERR "-------------------------------------\n";
	
	#don't need to print out the start and stop of the skip state, so start at $i = 1
	for(my $i = 1 ; $i < @Active; $i++){
		if($Active[$i]){
			print STDERR "$nodeStarts[$i]\t$nodeStops[$i]\t$IDs[$i]\t$Subfams[$i]\n";
		}
	}
}


#prints the 2D matrix that represents the path graph
sub PrintPathGraph(){

	for(my $i = 0 ; $i < $numnodes; $i++){
		print STDERR "$Changes[$i] ";
	}
	print STDERR "\n";
	
	for(my $i = 0; $i < $numnodes; $i++){
		for(my $j = 0; $j < $numnodes; $j++){
			print STDERR "$pathGraph[$i*$numnodes+$j]\t";
		}
		print STDERR " $Changes[$i]\n";
	}	
}


#stores matrix pos changes in an array along with the subfam its changing to 
#if there is one subfam that is in different rows in the matrices this will account for that
#by looking at the strings for the subfam names
sub GetChanges{	
	my ($changes, $changespos) = (@_);
		
	my $prev = 'skip';
	if($subfampath[0] ne ''){
# 		$prev = $Subfams[$subfampath[0]];
		$prev = $subfampath[0];
	}
	
	for(my $i = 1; $i < @subfampath; $i++){
		
		my $curr_subfam = 'skip';
		
		if($subfampath[$i] ne ""){
		
# 			$curr_subfam = $Subfams[$subfampath[$i]];
			$curr_subfam = $subfampath[$i];
		
			#if the subfam is an L1 in the form L1*_* just use the L1 name and not the end part
			#this is used to stitch the L1s back together 
			if($curr_subfam =~ /(.+?)_.+/g){
				$curr_subfam = $1;
			}
			
			if($curr_subfam ne $prev){
							
				push @$changespos, $i;
				push @$changes, $subfampath[$i];
			}
			
			$prev = $curr_subfam;

		}			
	}	
}




#FIXME - need to check on multiple examples and make sure it's reporting the correct start
# and end bounds 
#Finds starts/stop positions in the sequence of the different subfams/nodes
sub GetBounds{	
	my ($changes, $changespos) = (@_);
		
	my $prev = 'skip';
	if($subfampath[0] ne ''){
		$prev = $subfampath[0];
	}
	
	my $i;
	for($i = 1; $i < @subfampath; $i++){
		
		my $curr_subfam = 'skip';
		if($subfampath[$i] ne ''){
			$curr_subfam = $subfampath[$i];
		}
		
		if($curr_subfam =~ /(.+?)_.+/g){
			$curr_subfam = $1;
		}
		
		if($curr_subfam ne $prev){
			
			if($curr_subfam ne 'skip'){
				$nodeStarts[$subfampath[$i]] = $Columns[$i];
			}
			
			if($subfampath[$i-1] ne ''){
				$nodeStops[$subfampath[$i-1]] = $Columns[$i-1];
			}
		}
		
		$prev = $curr_subfam;
	}
	
	#FIXME - not sure if this gives the correct end value for the last node 
	$nodeStops[$subfampath[$i-1]] = $Columns[-1];
}


#prints matrix positions and subfams where the path changes to a different subfam
sub PrintChanges{
	my ($changes, $changespos) = (@_);
	
	for(my $i = 0; $i < @$changes; $i++){
		#changed the position to it's position in @columns, this make it so all the positions 
		#are correct even when full length subfams get spliced out
		print STDERR $Columns[@$changespos[$i]];
		print STDERR "\t";
		print STDERR "@$changes[$i]\n";
		
	}	
}


#find the min start and max stop for the whole region
sub Edges{
	my ($starts, $stops) = (@_);
	
	my $minStart = 9**9**9;
	my $maxStop = 0;
	for(my $i = 1; $i < @$starts; $i++){
		if(@$starts[$i] < $minStart){
			$minStart = @$starts[$i];
		}
		if(@$stops[$i] > $maxStop){
			$maxStop = @$stops[$i];
		}
	}
	return ($minStart, $maxStop);	
}


#pads sequences with '.' and makes sparse matrix - allows regions where seqs do not all have same 
#start and stop positions
sub padSeqs{
	my ($start, $stop, $subfamseq, $chromseq) = (@_);
	
	($startall, $stopall) = Edges($start, $stop);
					
	for(my $i = 1; $i < @SubfamSeqs; $i++){
		my $leftpad = @$start[$i] - $startall;
		my $rightpad = $stopall - @$stop[$i];	
	
		@$chromseq[$i] = ("." x $leftpad).@$chromseq[$i].("." x $rightpad);
		@$subfamseq[$i] = ("." x $leftpad).@$subfamseq[$i].("." x $rightpad);
		
	}
}



#calculates alignment score (according to crossmatch scoring) for every segment and all seqs
#fills align score matrix - pos 0 in the array is the score of the segment that starts at pos 0,
#pos 1 in array is score of segment that starts at pos 1, etc
#computes score for the first segment that does not start with a '.' and from there keeps the 
# base score and adds new char and subtracts old char - if a new gap is introduced, calls CalcScore()
#instead of adding on base score 
sub FillAlignScoreMatrix{
	my ($subfams, $chroms) = (@_);
	
	for(my $i = 1; $i < @ChromSeqs; $i++){
		
		my @subfam1 = split(//, @$subfams[$i]);
		my @chrom1 = split(//, @$chroms[$i]);
			
		#grab first chunk of 30 and calculate the raw score 
		
		#starts at the first non '.' char, but offsets it in the matrix based on where
		#the alignments start in the seq - ex: if first alignment in the seq starts at 10,
		#will offset by 10
		my $j = $Starts[$i] - $startall;
		my $index = $j; #index is the col we are in the align score matrix, $j is the place in @subfam1 and @chrom1

		my $offset = $chunksize;
		my $alignscore;
			
		my $tempindex = $j;
		my $tempcount = 0;
			
		while($tempcount < $chunksize){
			
			if($chrom1[$tempindex] ne '-'){
				$tempcount++;
			}
			
			$tempindex++;
		}
 				
		$offset = $tempindex - $j;
		my $prevoffset = $offset;
		
		#grabs first chunk - here $ j = pos of first non '.' char
		my @ChromSlice = @chrom1[$j .. $j+$offset-1];
		my @SubfamSlice = @subfam1[$j .. $j+$offset-1];
			
		#calculates score for first chunk and puts score in alignhash
		$alignscore = CalcScore(\@SubfamSlice, \@ChromSlice, '', '');
		$AlignHash{$i.'.'.$index} = $alignscore;
		
		if($cols < $index){
			$cols = $index+1;
		}
		$index++;
		
		
		#move to next chunk by adding next chars score and subtracting prev chars score
		
		while($j+$offset+1 < @chrom1){
				
			$tempindex = $j;
			$tempcount = 0;
			
			while($tempcount < $chunksize){
			
				if($chrom1[$tempindex+1] ne '-'){
					$tempcount++;
				}
				$tempindex++;
			}
				
			$offset = $tempindex - $j;
				
			if($chrom1[$j+1] ne '-'){
			
				if($chrom1[$j+1] ne '.' and $subfam1[$j+1] ne '.'){
				
					if($prevoffset != $offset){ #there is a new gap, or a gap was removed from beginning 
						my @ChromSlice = @chrom1[$j+1 .. $j+$offset];
						my @SubfamSlice = @subfam1[$j+1 .. $j+$offset];
						$alignscore = CalcScore(\@SubfamSlice, \@ChromSlice, $subfam1[$j], $chrom1[$j]);
						
					}else{
						#alignscore from previous segment - prev chars score + next chars score 
						
						#subtracting prev  chars score - tests if its a gap in the subfam as well
						if($subfam1[$j] eq '-'){
							if($subfam1[$j-1] eq '-'){
								$alignscore = $alignscore - $gapExt;
							}else{
								$alignscore = $alignscore - $gapInit;
							}
						}else{
							$alignscore = $alignscore - $SubMatrix[$CharPos{$subfam1[$j]}*$subMatrixCols+$CharPos{$chrom1[$j]}];
						}
					
					
						#adding next chars score - tests if its a gap in the subfam as well
						if($subfam1[$j+$offset] eq '-'){
							if($subfam1[$j+$offset-1] eq '-'){
								$alignscore = $alignscore + $gapExt;
							}else{
								$alignscore = $alignscore + $gapInit;
							}
						}elsif($subfam1[$j+$offset] eq '.' or $chrom1[$j+$offset] eq '.'){
								$alignscore = $alignscore;
						}else{
							$alignscore = $alignscore + $SubMatrix[$CharPos{$subfam1[$j+$offset]}*$subMatrixCols+$CharPos{$chrom1[$j+$offset]}];
						}
										
					}
					
					if($alignscore <= 0){
						$AlignHash{$i.'.'.$index} = 1;
					}else{
 						$AlignHash{$i.'.'.$index} = $alignscore;
 					} 					
		
				}
 				
				#finds max index value - which is assigned to num cols in align score matrix
				if($cols < $index){
					$cols = $index+1;
				}
			 	$index++;
			}
			$j++;
			$prevoffset = $offset;
		}	
		
	}
	
			
	#keep gaps in mind 
	#chunks can't start on gaps and gaps don't count when getting to the 30 bps 
}



#fills parallel array to the Align Matrix that holds the consensus position for each 
# subfam at that position in the alignment
sub FillConsensusPosMatrix{
	my ($consensus, $subfams, $chroms, $consensusstart, $consensusstop) = (@_);
	
	#start at 1 to skip 'skip state'
	for(my $i = 1; $i < $rows; $i++){	
		my @SubfamArray = split(//, @$subfams[$i]);
		my @ChromArray = split(//, @$chroms[$i]);
		
		my $consensuspos;
		if($Strands[$i] eq '+'){
			
			$consensuspos = @$consensusstart[$i]-1;
			my $matrixpos = 0;
			my $j = 0;
			for(my $j = 0; $j < @SubfamArray; $j++){
		
				if($SubfamArray[$j] ne '.'){
		
					#consensus pos only advances when there is not a gap in the subfam seq
					if($SubfamArray[$j] ne '-'){
						$consensuspos++;
					}
			
					#put consensus pos corresponding to pos in matrix in hash
					$consensus->{$i.".".$matrixpos} = $consensuspos;
				}
				#matrix position only advances when there is not a gap in the chrom seq
				if($ChromArray[$j] ne '-'){
					$matrixpos++;
				}
			}		
		
		}else{ #reverse strand 
		
			$consensuspos = @$consensusstart[$i]+1;
			my $matrixpos = 0;
			my $j = 0;
			for(my $j = 0; $j < @SubfamArray; $j++){
		
				if($SubfamArray[$j] ne '.'){
		
					#consensus pos only advances when there is not a gap in the subfam seq
					if($SubfamArray[$j] ne '-'){
						$consensuspos--;
					}
			
					#put consensus pos corresponding to pos in matrix in hash
					$consensus->{$i.".".$matrixpos} = $consensuspos;
				}

				#matrix position only advances when there is not a gap in the chrom seq
				if($ChromArray[$j] ne '-'){
					$matrixpos++;
				}
			}		
		}
		
		#FIXME - delete this later, just use it for testing 
		if(@$consensusstop[$i] != $consensuspos){
			print STDERR "\n\nERROR - consensus seq positions not correct\n\n";
		}
	
	}
}


#send in 2 seqs (char array) - return alignment score 
sub CalcScore{
	my ($seq1, $seq2, $lastpreva, $lastprevb) = (@_);
		
	my $chunkscore = 0;
	
	#deals with the first character of a segment being a gap character - have to look at last
	#segment to see if this is a gap init or ext
	if(@$seq1[0] eq '-'){
		if($lastpreva eq '-'){
			$chunkscore = $chunkscore + $gapExt;
		}else{
			$chunkscore = $chunkscore + $gapInit;
		}
	}elsif(@$seq2[0] eq '-'){
		if($lastprevb eq '-'){
			$chunkscore = $chunkscore + $gapExt;
		}else{
			$chunkscore = $chunkscore + $gapInit;
		}
	}elsif(@$seq1[0] eq '.' or @$seq2[0] eq '.'){
			$chunkscore = $chunkscore;
	}else{
		$chunkscore = $chunkscore + $SubMatrix[$CharPos{@$seq1[0]}*$subMatrixCols+$CharPos{@$seq2[0]}];
	}
		
	for(my $j = 1; $j < @$seq1; $j++){
		if(@$seq1[$j] eq '-'){
			if(@$seq1[$j-1] eq '-'){
				$chunkscore = $chunkscore + $gapExt;
			}else{
				$chunkscore = $chunkscore + $gapInit;
			}
		}elsif(@$seq2[$j] eq '-'){
			if(@$seq2[$j-1] eq '-'){
				$chunkscore = $chunkscore + $gapExt;
			}else{
				$chunkscore = $chunkscore + $gapInit;
			}
		}elsif(@$seq1[$j] eq '.' or @$seq2[$j] eq '.'){
			$chunkscore = $chunkscore;
		}else{
			if((exists $CharPos{@$seq1[$j]}) and (exists $CharPos{@$seq2[$j]})){
				$chunkscore = $chunkscore + $SubMatrix[$CharPos{@$seq1[$j]}*$subMatrixCols+$CharPos{@$seq2[$j]}];
				
			}
		}
	}
	
	# if($chunkscore <= 0){
# 		return 1;
# 	}
	
	return $chunkscore;	
}


#calculates alignment score with the crossmatch adjustment
#if score of a segment is less than 0, return 0 as score
sub CalcScoreAdjusted{
	my ($seq1, $seq2, $lastpreva, $lastprevb) = (@_);
		
	my $rawscore = CalcScore($seq1, $seq2, $lastpreva, $lastprevb);
	
	my $score = int($rawscore + (tSum($seq1, $seq2)/$lambda) + 0.999);
	
	if($score <= 0){
		return 1;
	}
	
	return $score;
}


#query sequence counts - query is chrom
#needed for score adjustment
sub GetCounts{
	my ($seq1, $seq2) = (@_);
	
	#actg
	my @counts = (0,0,0,0);
	
	for(my $k = 0; $k < @$seq2; $k++){
		if(@$seq1[$k] ne '-'){
			if(@$seq2[$k] eq 'A' or @$seq2[$k] eq 'a'){
				$counts[0]++;
			}elsif(@$seq2[$k] eq 'C' or @$seq2[$k] eq 'c'){
				$counts[1]++;
			}elsif(@$seq2[$k] eq 'T' or @$seq2[$k] eq 't'){
				$counts[2]++;
			}elsif(@$seq2[$k] eq 'G' or @$seq2[$k] eq 'g'){
				$counts[3]++;
			}
		}
	}
	
	#ensures no log of 0 errors, if a segment doesn't have a certain letter instead of a 0 
	#in the counts array I put a 1
	for(my $k = 0; $k < @counts; $k++){
		if($counts[$k] == 0){
			$counts[$k] = 1;
		}	
	}
	return @counts;
}


#calculates t_factor using query sequence composition counts 
#needed for score adjustment
sub tFactor{
	my ($seq1, $seq2) = (@_);
	
	my @querycounts = GetCounts($seq1, $seq2);
	
	my $t_factor = 0;
	
	foreach my $count (@querycounts){
		if($count == 0){
			$t_factor = $t_factor;
		}else{
			$t_factor = $t_factor + ($count * log($count));
		}	
	}
	
	my $sum = sum @querycounts;
	$t_factor = $t_factor - (($sum) * log($sum));
	
	return "$t_factor\n";
}


#calculates t_sum
#needed for score adjustment
sub tSum{
	my ($seq1, $seq2) = (@_);
	
	#actg
	#hard coded these in here for now
	my @backgroundFreqs = (0.295, 0.205, 0.295, 0.205); #from matrix

	my @querycounts = GetCounts($seq1, $seq2);
	
	my $t_sum = 0;
	
	for(my $k = 0; $k < @backgroundFreqs; $k++){
		$t_sum = $t_sum + ($querycounts[$k] * log ($backgroundFreqs[$k]));	
	}
	
	return $t_sum - tFactor($seq1, $seq2);
}


#FIXME - for now this just does score from CM confidence - add hmm later, and maybe RM
#send in an array of scores for a segment - output an array of confidence values for the segment
sub ConfidenceCM{

	my ($lambda, $region) = (@_);
	my $confidenceString = '';
	
	#loops through the array once to get the sum of 2^every_hit_score in region 
	#converts the score to account for lambda before summing 
	my $ScoreTotal = 0;
	foreach my $Score (@$region){
	
		if($Score > 0){
			my $convertedScore = $Score * $lambda;
			$ScoreTotal = $ScoreTotal + (2 ** $convertedScore);
		}
	}
	
	#once region score is summed, loop back through the region array and calculate
	#confidence for each hit
	foreach my $Score (@$region){
		if($Score > 0){
			my $convertedScore = $Score * $lambda;
			my $confidence = ((2**$convertedScore) / $ScoreTotal);	
			if($confidenceString ne ""){
				$confidenceString = $confidenceString." ".$confidence;
			}else{
				$confidenceString = $confidence;
			}
		}else{
			if($confidenceString ne ""){
				$confidenceString = $confidenceString." ".'0';
			}else{
				$confidenceString = '0';
			}
		}
	}
	return $confidenceString;
}


#FIXME - make it so don't have to loop through all the i's and j's to fill this - must be better way
#fill confidence score matrix based on align score matrix hash
sub FillConfScoreMatrix{
	my ($alignhash, $confhash) = (@_);
	
	
	for(my $i = 0; $i < @Columns-($chunksize-1); $i++){
		my $col = $Columns[$i];
		my @temp;
		for(my $row = 0; $row < $rows; $row++){
			if($alignhash->{$row.'.'.$col}){
				push @temp, $alignhash->{$row.'.'.$col};
			}else{
				push @temp, 0;
			}
		}
	
		my @confidenceTemp = split(/ /, ConfidenceCM($lambda, \@temp));
		for(my $row = 0; $row < $rows; $row++){
			if($confidenceTemp[$row] != 0){
				$confhash->{$row.'.'.$col} = $confidenceTemp[$row];
			}
		}	
	}
}


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
sub FillProbMatrix{
	my ($probhash, $supporthash, $originhash) = (@_);
	
	#prints input into the python version of FillProbMatrix
	# print "key support_value\n";
# 	for(my $i = 0; $i < $rows; $i++){
# 		for(my $j = 0; $j < $cols; $j++){
# 			if(exists $supporthash->{$i.'.'.$j}){
# 				print $i.'.'.$j;
# 				print " ";
# 				print $supporthash->{$i.'.'.$j};
# 				print "\n";
# 			}
# 		}
# 	}
# 	my $time  = Time::HiRes::gettimeofday();
		
	my $j = 1; 
	for(my $col = 1; $col < @Columns; $col++){
	
		#if a column is undefined, increments j instead of trying to index $Columns
		#use $col when indexing $columns, if not - use $j (columns already indexed)
		if(exists $Columns[$col]){
			$j = $Columns[$col];
		}else{
			$j++;
		}
		
		
		foreach my $i (@{$ActiveCells_collapse{$j}}){
			
			my $max = -inf;
			my $maxindex;
	
				my $supportlog = log($supporthash->{$i.'.'.$j});
				 
				#loops through all the active rows in the previous column
				foreach my $row (@{$ActiveCells_collapse{$Columns[$col-1]}}){
				
					my $score;
					my $there = 0; #to see if the score was calculated
					
					if(exists $Columns[$col]){
						$score = $supportlog + $probhash->{$row.'.'.($Columns[$col-1])};
						$there = 1;
					}else{
						$score = $supportlog + $probhash->{$row.'.'.($j-1)};
						$there = 1;
					}
				
 					if($there){
						if($row eq $i){
							#come from same seq - higher prob
							$score = $score +  $sameProbLog;
						}else{
							$score = $score + $changeProbLog;
						}
											
						#find the max score, it's index goes into the origin matrix 
						if($score > $max){
							$max = $score;
							$maxindex = $row;
						}
					}
					
				}
							
				#fills the 2D loc array with the origin location 
				$probhash->{$i.'.'.$j} = $max;
				$originhash->{$i.'.'.$j} = $maxindex;
			
		}
		
	}
	#prints the output correct output put FillProbMatrix
	# print "key\tprob_value\torigin\n";
# 	for(my $i = 0; $i < $rows; $i++){
# 		for(my $j = 0; $j < $cols; $j++){
# 			if(exists $supporthash->{$i.'.'.$j} and $supporthash->{$i.'.'.$j} != 0){
# 				print $i.'.'.$j;
# 				print "\t";
# 				print $probhash->{$i.'.'.$j};
# 				print "\t";
# 				
# 				if($j != 0){
# 					print $originhash->{$i.'.'.$j};
# 				}else{
# 					print 'NaN';
# 				}
# 				print "\n";
# 			}
# 		}
# 	}
# 	
# 	print STDERR "prob scores - in function : ";
# 	printf STDERR ("%.2fs\n", abs($time - Time::HiRes::gettimeofday()));
	
}


# Fills support score matrix using values in conf matrix
#score for subfam x at position i is sum of all confidences for subfam x for all segments that 
#overlap position i - divided by number of segments
sub FillSupportMatrix{
	my ($supporthash, $alignhash, $confhash) = (@_);
		
	# i = subfam, j = position
	for(my $i = 0; $i < $rows; $i++){
		my $tempcol;
		for(my $col = 0; $col < @Columns-($chunksize-1); $col++){
			my $j = $Columns[$col];
			
			if(exists $confhash->{$i.'.'.$j}){
				my $num = $j;
				my $sum = 0;
				my $numsegments = 0;
				while($num >= 0 and $num >= $j-$chunksize+1){
					if(exists $confhash->{$i.'.'.$num}){
						$sum = $sum + $confhash->{$i.'.'.$num};
					
						$numsegments++;
					
					}
					$num--;
				}
						
				if($numsegments){
					$supporthash->{$i.'.'.$j} = $sum/$numsegments;
				}
			}
		}	
		
		#@Columns only goes until the end of align score matrix so it is $chunksize smaller than
		#we need here - do this last part of the loop so fill in the rest of the support matrix	
		for(my $j = $cols-$chunksize; $j <= $cols; $j++){
			
			my $num = $j;
			my $sum = 0;
			my $numsegments = 0;
			while($num >= 0 and $num >= $j-$chunksize+1){
				if(exists $confhash->{$i.'.'.$num}){
					$sum = $sum + $confhash->{$i.'.'.$num};
					
					$numsegments++;
					
				}
				$num--;
			}
						
			if($numsegments){
				$supporthash->{$i.'.'.$j} = $sum/$numsegments;
			}
		}
	}
}


#using origin matrix, back traces through the 2D array to get the subfam path
#reverses it and returns array of the path 
#jumps over removed columns when necessary
sub GetPath{
	my ($prob, $origin, $subfams) = (@_);
	
	my $max = -inf;	
	my $maxindex;
# 	for(my $i = 0; $i < $rows; $i++){
	foreach my $i (@{$ActiveCells_collapse{$cols-1}}){
		if($max < $prob->{$i.'.'.($cols-1)}){
			$max = $prob->{$i.'.'.($cols-1)};
			$maxindex = $i;
		}	
	}
		
	my @subfamorder;
	
	my $prev = $origin->{$maxindex.'.'.($cols-1)};
	
		
	my $i = $cols -1;	# $i = column loop is on
	push @subfamorder, $maxindex;
	
	for(my $col = @Columns; $col > 0; $col--){
	
		#deals with last $chunksize-1 cols that don't exist in @Columns because they are not in align score matrix
		if(exists $Columns[$col]){
			$i = $Columns[$col];
		}else{
			$i--;
		}

		#deals with last $chunksize-1 cols that don't exist in @Columns because they are not in align score matrix
		if(exists $Columns[$col]){
			if((exists $origin->{$prev.'.'.($Columns[$col-1])}) and (exists $origin->{$prev.'.'.($i)})){
# 				push @subfamorder, @$subfams[$prev];
				push @subfamorder, $prev;
				$prev = $origin->{$prev.'.'.($Columns[$col-1])};
			}else{
				push @subfamorder, "";
			}
		}else{
			if((exists $origin->{$prev.'.'.($i-1)}) and (exists $origin->{$prev.'.'.($i)})){
# 				push @subfamorder, @$subfams[$prev];
				push @subfamorder, $prev;
				$prev = $origin->{$prev.'.'.($i-1)};
			}else{
				push @subfamorder, "";
			}
		}
	}	
	return reverse @subfamorder;
}

#fills node confidence matrix 
#first fills matrix with node alignment scores, then reuses matrix for confidence scores 
sub NodeConfidence{
	my ($nodeconfidence, $subfamseqs, $chromseqs, $changespos) = (@_);
	
	my @nodeconfidence_temp;
	#initialize node confidence
	for(my $i = 0; $i < $rows*$numnodes; $i++){
		$nodeconfidence_temp[$i] = 0; 
	}

	#compute alignment scores for node region
	
	#does first node in the sequence - does not look back to previous char incase of gap ext
	for(my $j = 1; $j < @Subfams; $j++){
		my @subfam = split(//, substr @$subfamseqs[$j], @$changespos[0]-1, @$changespos[1]-@$changespos[0]);
		my @chrom = split(//, substr @$chromseqs[$j], @$changespos[0]-1, @$changespos[1]-@$changespos[0]);
		my $alignscore = CalcScore(\@subfam, \@chrom, '', '');
		$nodeconfidence_temp[$j*$numnodes+0] = $alignscore;
	}

	#does rest of nodes - looks back at prev char incase of gap ext		
	for(my $i = 0; $i < $numnodes-1; $i++){
		for(my $j = 1; $j < @Subfams; $j++){
			my @subfam = split(//, substr @$subfamseqs[$j], @$changespos[$i]-1, @$changespos[$i+1]-@$changespos[$i]);
			my @chrom = split(//, substr @$chromseqs[$j], @$changespos[$i]-1, @$changespos[$i+1]-@$changespos[$i]);
			my $alignscore = CalcScore(\@subfam, \@chrom, (substr $SubfamSeqs[$j], @$changespos[$i+1]-1, 1), (substr @$chromseqs[$j], @$changespos[$i+1]-1, 1));
			$nodeconfidence_temp[$j*$numnodes+$i] = $alignscore;
		}
	}
	
	#does last node
	for(my $j = 1; $j < @Subfams; $j++){
		my @subfam = split(//, (substr @$subfamseqs[$j], @$changespos[-1]-1));
		my @chrom = split(//, (substr @$chromseqs[$j], @$changespos[-1]-1));
		my $alignscore = CalcScore(\@subfam, \@chrom, '', '');
		$nodeconfidence_temp[$j*$numnodes+$numnodes-1] = $alignscore;
	}

	#reuse same matrix and compute confidence scores for the nodes
	for(my $j = 0; $j < $numnodes; $j++){
		my @temp;
		for(my $i = 0; $i < @Subfams; $i++){  ##FIXME - should this start at 1 (because of skip state)?
			push @temp, $nodeconfidence_temp[$i*$numnodes+$j];
		}
	
		my @confidenceTemp = split(/ /, ConfidenceCM($lambda, \@temp));
		for(my $i = 0; $i < @Subfams; $i++){   ##FIXME - should this start at 1 (because of skip state)?
			$nodeconfidence_temp[$i*$numnodes+$j] = $confidenceTemp[$i];
		}
	}
	
	
	#collapse nodeconfidence down same say supportmatrix is collapsed - all seqs of 
	#the same subfam are put in the same row
	#not a sparse hash - holds the 0s, but I think this is okay because it won't ever
	#be a very large matrix, and this way we don't have to test if anything exists 	
# 	my %nodeconfidence_collapse;
	for(my $j = 0; $j < $numnodes; $j++){
		for(my $i = 0; $i < @Subfams; $i++){
			if(exists $nodeconfidence->{$Subfams[$i].'.'.$j}){
				$nodeconfidence->{$Subfams[$i].'.'.$j} = $nodeconfidence->{$Subfams[$i].'.'.$j} + $nodeconfidence_temp[$i*$numnodes+$j];
			}else{
				$nodeconfidence->{$Subfams[$i].'.'.$j} = $nodeconfidence_temp[$i*$numnodes+$j];			
			}	
		}
	}
	
}


##FIXME - currently ignores L1 special cases, and also ignores checking to see if the 
# source is before the sink in the consensus seq
sub FillPathGraph{
	my ($pathgraph) = (@_);

	#initializing path graph
	for(my $i = 0; $i < $numnodes*$numnodes; $i++){
		push @$pathgraph, 0;
	}

	# filling beginning path graph with straight line through the nodes
	for(my $i = 0; $i < $numnodes-1; $i++){
		@$pathgraph[$i*$numnodes+$i+1] = 1;
	}	
		
# 	FIXME - have to loop through all the subfam names at each node to see if they match - might be a 
# 	better way to match them later?
	for(my $j = 0; $j < $numnodes; $j++){
		my $sinkSubfam = $Changes[$j];
		
		my $sinkSubfamStart = $ConsensusHash_collapse{$sinkSubfam.'.'.$Columns[$ChangesPos[$j]]};
		my $sinkStrand = $StrandHash_collapse{$sinkSubfam.'.'.$Columns[$ChangesPos[$j]]};

# 		print STDERR "sink subfam: $sinkSubfam\n";
# 		print STDERR "sink subfam start: $sinkSubfamStart\n";		
						
		#looks at all the preceding nodes, except the one directly before it (source nodes)
		for(my $i = 0; $i < $j-1; $i++){
					
			#look at all the subfams in each node 
			foreach my $sourceSubfam (keys %Subfams_collapse){
				
				my $sourceConf = $NodeConfidence{$sourceSubfam.'.'.$i};
				
				if( exists $ConsensusHash_collapse{$sourceSubfam.'.'.($Columns[$ChangesPos[$i+1]]-1)}){
					my $sourceSubfamStop = $ConsensusHash_collapse{$sourceSubfam.'.'.($Columns[$ChangesPos[$i+1]]-1)};
					my $sourceStrand = $StrandHash_collapse{$sourceSubfam.'.'.($Columns[$ChangesPos[$i+1]]-1)};

# 					print STDERR "$sourceSubfam\t";
# 					print STDERR "$Columns[$ChangesPos[$i+1]]-1\t";
# 					print STDERR "$sourceSubfamStop\n";
					
# 					adds in edge if the subfam of the sink is at the source node and if it's 
# 					confidence >= 1%, and if the source is before the sink in the consensus sequence 
					if($sinkStrand eq '+' and $sinkStrand eq $sourceStrand){
						if($sinkSubfam eq $sourceSubfam and $sourceConf >= .01){
 							#FIXME- not sure what this overlap should be .. just allowed 50 for now
 							if($sourceSubfamStop <= $sinkSubfamStart+50){
								@$pathgraph[$i*$numnodes+$j] = 1;	
							}
						}
					}elsif($sinkStrand eq '-' and $sinkStrand eq $sourceStrand){
						if($sinkSubfam eq $sourceSubfam and $sourceConf >= .01){
 							#FIXME- not sure what this overlap should be .. just allowed 50 for now
 							if($sourceSubfamStop >= $sinkSubfamStart+50){
								@$pathgraph[$i*$numnodes+$j] = 1;	
							}
						}
					}
					
				}	
				
			}
		}
# 		print STDERR "\n\n";
	}
	
# 	PrintPathGraph();
}


#finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds them to
# @RemoveStarts and @RemoveStops so they can be extracted from the alignment 
sub ExtractNodes{
	my ($removestarts, $removestops, $changespos, $pathgraph, $numnodes) = (@_);

	#boolean for which nodes will be removed
	my @RemoveNodes = (0)x$numnodes;
	
	#extracting nodes that only have one incoming and one outgoing edge
	my @NumEdgesIn;
	my @NumEdgesOut;
	for(my $i = 0; $i < $numnodes; $i++){
		push @NumEdgesIn, 0;
		push @NumEdgesOut, 0;
	}

	for(my $i = 0; $i < $numnodes; $i++){
		for(my $j = 0; $j < $numnodes; $j++){
			#$j - incoming
			$NumEdgesIn[$j] += @$pathgraph[$i*$numnodes+$j];
		
			#$i - outgoing
			$NumEdgesOut[$i] += @$pathgraph[$i*$numnodes+$j];

		}
	}

	#$numnodes-1 because if $i = $numnodes-1 @$changespos[$i+1] is undefined because it
	#isn't a change position, it is the end of the matrix
	for(my $i = 0; $i < $numnodes-1; $i++){
		if($NumEdgesIn[$i] <= 1 and $NumEdgesOut[$i] <= 1){
			push @$removestarts, @$changespos[$i];
			push @$removestops, @$changespos[$i+1];
			
			$RemoveNodes[$i] = 1;
		}
	}	
	
	#deals with last node, so when $numnodes-1 the last remove stop is the end of the matrix
	if($NumEdgesIn[$numnodes-1] <= 1 and $NumEdgesOut[$numnodes-1] <= 1){
		push @$removestarts, @$changespos[$numnodes-1];
		push @$removestops, $cols;
		
		$RemoveNodes[$numnodes-1] = 1;
	}
	
	#if remove the end nodes, need to ignore the last columns in the matrix that
	#correspond with those nodes - this makes it so when back tracking through the
	#matrix we start at the corrext column instead of starting at the end 
	my $i = $numnodes-1;
	while($RemoveNodes[$i]){
		$cols = @$changespos[$i] - 1;
		$i--;
	}	
}






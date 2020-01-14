# This is a piece of the AdjudicateRegions.pl file that loads the substitution
# matrix so that we can compare it against the result of the Python code.

my $infile_matrix = shift;
open(my $in_matrix,  "<",  $infile_matrix)
        or die "Could not open $infile_matrix";
        
# CharacterPositions
my %CharPos;        
my $line = <$in_matrix>;
$line =~ s/^\s+//;
my @chars = split(/\s+/, $line);
for(my $i = 0; $i < @chars; $i++){
	$CharPos{$chars[$i]} = $i;
}
my $subMatrixCols = scalar @chars;

# SubstitutionMatrix
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

use Data::Dumper;
print Dumper(\%CharPos);
print "\n";
print "@SubMatrix\n";

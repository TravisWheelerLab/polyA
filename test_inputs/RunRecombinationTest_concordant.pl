#!/usr/bin/perl
use strict;
use warnings;

my %pairs;

my @files = <"/Users/kaitlincarey/Desktop/NOTEBOOK/2020/DiscordantPaper/Aug3-RecombinationPlot/ALIGN_CONCORDANT/*.format">;

my $count = 0;
my $recomb_count = 0;
foreach my $file (@files) {
    $count = $count + 1;

    print "$count\n";

    print "$file\n";

    my $SD_num;
    if($file =~ /ALIGN_CONCORDANT\/(\d+)_.+format/g){
#    if($file =~ /TEST\/(\d+)_.+format/g){
        $SD_num = $1;
    }

    my $output = `poetry run python3 AdjudicateRegions.py --seqpos --lambda .1227 $file 25p41g_edited.matrix`;

    my @output_array = split(/\n/, $output);
    my $recomb = scalar @output_array;

    if(exists $pairs{$SD_num}){
        $pairs{$SD_num} = $pairs{$SD_num} + $recomb;
#        print "$pairs{$SD_num}\n";
    }else{
        $pairs{$SD_num} = $recomb
    }

}


foreach my $key_num (keys %pairs){
    if($pairs{$key_num} > 2){
        $recomb_count++;
    }
}


print "recombinations: $recomb_count\n";
print "total pairs: $count\n";


#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $infile_align = shift;
my $outfile = shift;

open(my $in_align,  "<",  $infile_align)
        or die "Could not open $infile_align";
open(my $out,  ">",  $outfile)
        or die "Could not open $outfile";      

$/ = "Transitions / transversions";

my $prev_chrom = '';

my $count = 0;
while(my $region = <$in_align>){
    $region =~ s/Transitions \/ transversions//;
#     $count++;
        
    if($region =~ / \=*.+?\n(\s*?\d+\s+[0-9]+\.[0-9]+\s+[0-9.]+\s+[0-9.]+\s+.+?)\n\n([\s\S]+)/g){
 
		$count++;
 		
		my $info = $1;
		my $alignment = $2;	
		$info =~ s/^\s+//;
				
		my @infoline = split(/\s+/, $info);
				
		my $score = $infoline[0];
# 		my $chrom = $infoline[4];  
		my $chrom = 'chr0:0000-0000';  #doesn't have chrom info because it's an artificial seq
		my $start = $infoline[5];
		my $stop = $infoline[6];
		my $subfam;
		my $consensusStart;
		my $consensusStop;
		my $strand;
		my $flank;
		

		#if 6 contains (), then start/stop = 7/8
		if($infoline[6] =~ /\(\d+\)/){
			$start = $infoline[6];
			$stop = $infoline[7];
		}
		
		
		if(scalar @infoline == 12){ #not compliment
			
			$strand = '+';
			$subfam = $infoline[8];
			$consensusStart = $infoline[9];
			$consensusStop = $infoline[10];
			$flank = $infoline[11];
			
			#if 10 contains (), then start/stop = 7/8
			if($infoline[9] =~ /\(\d+\)/){
				$flank = $infoline[9];
				$start = $infoline[10];
				$stop = $infoline[11];
			}
			
		}elsif(scalar @infoline == 13){ #compliment
		
			$strand = '-';
			$subfam = $infoline[9];
			$consensusStart = $infoline[10];
			$consensusStop = $infoline[11];
			$flank = $infoline[12];
			
			#if 11 contains (), then start/stop = 7/8
			if($infoline[10] =~ /\(\d+\)/){
				$flank = $infoline[10];
				$consensusStart = $infoline[11];
				$consensusStop = $infoline[12];
			}
					
		}else{
			print STDERR "\nERROR... regex not working\n\n";
		}

		if($flank =~ /\((\d+)\)/g){
			$flank = $1;
		}	
				
		print $out "Align: $subfam\t$chrom\t$score\t$strand\t$start\t$stop\t$consensusStart\t$consensusStop\t$flank\n"; 						
 			
 		my @Alignment = split(/\n/, $alignment);
 		
 			 				
 		my $chromSeq = "";
 		my $subfamSeq = "";
 				
 		for(my $i = 0; $i < @Alignment; $i = $i + 4){
 			if($Alignment[$i] =~ /.+?\s\d+\s(.+?)\s\d+/){
 				$chromSeq = $chromSeq.$1;
 			}
 			if($Alignment[$i+2] =~ /.+?\s\d+\s(.+?)\s\d+/){
 				$subfamSeq = $subfamSeq.$1;
 			}
 		}
 			
 					
 		my @ChromSeq = split(//, $chromSeq);
 		my @SubfamSeq = split(//, $subfamSeq);
 				
 		print $out ">$chrom\n";
 		my $j = 0;
 		for($j = 0; $j < @ChromSeq; $j++){
 			print $out "$ChromSeq[$j]";
 			if(($j+1) % 60 == 0){
 				print $out "\n";
 			}
 		}
 		
 		if($j % 60 != 0){
 			print $out "\n";
 		}
 		
 				
 		print $out ">$subfam\n";
 		$j = 0;
 		for($j = 0; $j < @SubfamSeq; $j++){
 			print $out "$SubfamSeq[$j]";
 			if(($j+1) % 60 == 0){
 				print $out  "\n";
 			}
 		}
 		if($j % 60 != 0){
 			print $out "\n";
 		}
	}
}
close $in_align;

# print "\n$count\n";



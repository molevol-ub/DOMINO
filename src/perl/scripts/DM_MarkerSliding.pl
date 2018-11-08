#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
	require Parallel::ForkManager;
	require File::Path; use File::Path qw(remove_tree);
}
#################################################################
my $domino_version ="DOMINO v1.1 ## Revised 07-11-2018";
my $scripts_path = $FindBin::Bin."/../";
my $domino_Scripts = $scripts_path."scripts";
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
#################################################################

## Arguments
my $path = $ARGV[0];
my $id = $ARGV[1]; 
my $seq = $ARGV[2]; 
my $output_file = $ARGV[3];
my $dna_seq = $seq; my $seqlen = length($seq);
#################################################################

#################################################################
## Get general parameters and files
my $hash_parameters = DOMINO::get_parameters($path."/", "markers");
#print Dumper $hash_parameters;
#################################################################

open (OUT, ">$output_file");

my $window_var_CONS = $$hash_parameters{'marker'}{'window_var_CONS'}[0];
my $variable_divergence= $$hash_parameters{'marker'}{'variable_divergence'}[0];

## Loop
for (my $i=0; $i < $seqlen; $i += $$hash_parameters{'marker'}{'SLIDING_user'}[0]) {
	# AAATATGACTATCGATCGATCATGCTACGATCGATCGATCGTACTACTACGACTGATCGATCGATCGACGACTGAC
	# 		P1		P2|P3							P4|P5						P6
	my $coord_P1 = $i; #print "P1: $coord_P1\n";
	for (my $h=$$hash_parameters{'marker'}{'window_size_CONS_min'}[0]; $h <= $$hash_parameters{'marker'}{'window_size_CONS_max'}[0]; $h += $$hash_parameters{'marker'}{'C-SI_inc'}[0]) {
		my $coord_P2 = $i + $h; 			#print "P2: $coord_P2\n";
		my $coord_P3 = $coord_P2 + 1; 		#print "P3: $coord_P3\n";
		if ($coord_P3 > $seqlen) {next;} 

		for (my $j = $$hash_parameters{'marker'}{'window_size_VARS_min'}[0]; $j <= $$hash_parameters{'marker'}{'window_size_VARS_max'}[0]; $j += $$hash_parameters{'marker'}{'V-SI_inc'}[0]) {
			my $coord_P4 = $coord_P3 + $j; #print "P4: $coord_P4\n";
			my $coord_P5 = $coord_P4 + 1; #print "P5: $coord_P5\n";
			my $coord_P6 = $coord_P5 + $h; #print "P6: $coord_P6\n";

			if ($coord_P4 > $seqlen){next;} 
			if ($coord_P5 > $seqlen){next;} 
			if ($coord_P6 > $seqlen){next;}				

			#print "\n******\n";
			#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\nPositions:\n";
			#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
			#print "******\n\n";

			## Get substrings
			my $length_string_P1_P2 = $coord_P2 - $coord_P1;
			my $string_P1_P2 = substr($dna_seq, $coord_P1, $length_string_P1_P2);

			my $length_string_P3_P4 = $coord_P4 - $coord_P3;
			my $string_P3_P4 = substr($dna_seq, $coord_P3, $length_string_P3_P4);

			my $length_string_P5_P6 = $coord_P6 - $coord_P5;
			my $string_P5_P6 = substr($dna_seq, $coord_P5, $length_string_P5_P6);
			
			## Count the variable positions in each 
			my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//; ## Conserved 
			my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//; ## Variable
			my $count_string_P5_P6 = $string_P5_P6  =~ tr/1//; ## Conserved
			my $count_N_P1_P2 = $string_P1_P2 =~ tr/N//; ## Conserved 
			my $count_N_P5_P6 = $string_P5_P6  =~ tr/N//; ## Conserved
			
			## More than variations allowed in conserved
			if ($count_string_P1_P2 > $window_var_CONS) { next; } 
			if ($count_string_P5_P6 > $window_var_CONS) { next; }
			if ($count_N_P1_P2 > $window_var_CONS) { next; } 
			if ($count_N_P5_P6 > $window_var_CONS) { next; } 

			#Get marker profile					
			my $total_length = $coord_P6 - $coord_P1;
			my $marker_profile = $string_P1_P2.$string_P3_P4.$string_P5_P6;
			
			### Check variation
			my $expected_var_sites;	my $flag_error=0;
			if ($variable_divergence) {
				# If a minimun divergence, get the expected variable sites for the length
				$expected_var_sites = int($variable_divergence * $total_length);
				unless ($count_string_P3_P4 >= $expected_var_sites) { $flag_error++; }			
			} else {
				# User provides a number of variable positions or a range
				if ($expected_var_sites > $$hash_parameters{'marker'}{'variable_positions_user_min'}[0]) {
					unless ($expected_var_sites < $$hash_parameters{'marker'}{'variable_positions_user_max'}[0]) {
						$flag_error++;
					}} else { $flag_error++; }
			}
			
			## Missing % of bases missing
			my $missing_count = $marker_profile =~ tr/N//;
			my $missing_count2 = $marker_profile =~ tr/-//;
			$missing_count += $missing_count2;
			my $missing_count_percent = ($missing_count/$total_length)*100;
			my $percent_total_length = $total_length * $$hash_parameters{'marker'}{'missing_allowed'}[0];  ## Default 0.05

			if ($flag_error > 0) {	next; } #if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; }
			if ($missing_count_percent < $percent_total_length) {
				print OUT "$id\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\n";

				#if ($debugger) {
					#print "\n\n***********************************************\n";
					#print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
					#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\n\nPositions:\n";
					#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
					#print "\nSubsets\n";
					#print "Conserved (P1-P2): ".length($string_P1_P2)."\n"; print "$string_P1_P2\n";
					#print "Variable (P3-P4): ".length($string_P3_P4)."\n";  print $string_P3_P4."\n";
					#print "Conserved (P5-P6): ".length($string_P5_P6)."\n"; print $string_P5_P6."\n";
					#print "\nVariations:\n";
					#print "Count_string_P1_P2: $count_string_P1_P2\n";
					#print "Count_string_P3_P4: $count_string_P3_P4\n";
					#print "Count_string_P5_P6: $count_string_P5_P6\n";
					#print "\nMarker:\n";
					#print "Total Length:$total_length\tVariable Positions:$count_string_P3_P4\tExpected:$expected_var_sites\n";
					#print "Missing allowed: $percent_total_length %\tMissing:$missing_count_percent %\n";
					#print "***********************************************\n";
				#}
} else { next;	}}}}
close (OUT);
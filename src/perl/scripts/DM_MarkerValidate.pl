#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
}
#################################################################
## Overlaps and maximizes domino markers obtained
my $file = $ARGV[0]; 				## markers_shared 
my $dir = $ARGV[1];
my $DOMINO_markers_file = $ARGV[2];
my $ref_taxa_all = $ARGV[3];
my $path = $ARGV[4];
#################################################################
my $domino_version ="DOMINO v1.1 ## Revised 07-11-2018";
my $hash_parameters = DOMINO::get_parameters($path."/");
my %new_domino_files = %{ DOMINO::get_DOMINO_files($path."/") };
##print Dumper \%new_domino_files; exit();
#################################################################

#################################################################
## Check each group of overlapping markers
open (MARKERS, "$DOMINO_markers_file") or die "Could not open file $DOMINO_markers_file [DM_MarkerScan: check_DOMINO_marker]";
my $j = 1; my %hash_array;
while (1) {
	my $chunk;
	if (!eof(MARKERS)) { 
		$/ = "\/\/"; ## Telling perl where a new line starts
		$chunk = <MARKERS>; 
		push(@{$hash_array{$j}}, $chunk);
	} 
	$j++; 
	last if eof(MARKERS);
}
$/ = "\n";	
#print Dumper \%hash_array;
#################################################################

my $marker_number = 0;
foreach my $group (sort {$a <=> $b} keys %hash_array) {
	## Split chunk push into array
	my @array_markers_tmp = @{ $hash_array{$group} };
	my @array_markers;
	for (my $i=0; $i < scalar @array_markers_tmp; $i++) {
		#print "CHUNK!\n"; print $array_markers_tmp[$i]."\n"; print "********\n";
		my @array = split("\n", $array_markers_tmp[$i]);
		for (my $j=0; $j < scalar @array; $j++) {
			chomp($array[$j]);
			if ($array[$j] =~ ".*:.*") {
				push(@array_markers, $array[$j]);
	}}}
	## For each group of markers, evaluate the features
	for (my $i=0; $i < scalar @array_markers; $i++) {
		my @contig_name = split("##", $array_markers[$i]);
		my @array = split(":", $contig_name[1]);
		my $coord_P1 = $array[0]; my $coord_P2 = $array[1];
		my $coord_P3 = $array[2]; my $coord_P4 = $array[3];
		my $coord_P5 = $array[4]; my $coord_P6 = $array[5];
		my $taxa = $array[6];
		# print $contig_name[0]."\t".join(",",@array)."\n"; exit();
		
		## Retrieve msa for this marker			
		my %hash; my %fasta_msa_sub;
		if ($$hash_parameters{"marker"}{"option"}[0] eq "msa_alignment") {
			my $fasta_msa_sub_ref = DOMINO::readFASTA_hash($ref_taxa_all);
			%fasta_msa_sub = %{ $fasta_msa_sub_ref };
		} else {
			my @desired_taxa = split(",", $taxa);
			for (my $j=0; $j < scalar @desired_taxa; $j++) {
				next if ($ref_taxa_all eq $desired_taxa[$j]);
				my $fasta_each_taxa = $new_domino_files{$desired_taxa[$j]}{"PROFILE::Ref:".$ref_taxa_all}[0]."/".$contig_name[0]."_sequence.fasta";
				if (-f $fasta_each_taxa) {
					$fasta_msa_sub{$desired_taxa[$j]} = $fasta_each_taxa;
				} else { 
					$fasta_msa_sub{$desired_taxa[$j]} = "missing"; 
			}}
			my $reference_file_contig = $new_domino_files{$ref_taxa_all}{'REF_DIR'}[0]."/".$contig_name[0].".fasta";
			$fasta_msa_sub{$ref_taxa_all} = $reference_file_contig;
		}
		foreach my $keys (sort keys %fasta_msa_sub ) {			
			if ($fasta_msa_sub{$keys} =~ /missing/) {next;}
			if ($fasta_msa_sub{$keys} =~ /.*fasta/) {
				open (FILE, $fasta_msa_sub{$keys});
				my ($titleline, $sequence);
				$/ = ">"; ## Telling perl where a new line starts
				while (<FILE>) {		
					next if /^#/ || /^\s*$/;
					chomp;
					($titleline, $sequence) = split(/\n/,$_,2);
					next unless ($sequence && $titleline);
					chomp $sequence;
					$sequence =~ s/\s+//g; $sequence =~ s/\r//g; $titleline =~ s/\r//g;
				}
				close(FILE); $/ = "\n";
				my ($seq_id, $seq) = DOMINO::fetching_range_seqs($keys, $coord_P1, $coord_P6, $sequence);
				$hash{$keys} = $seq;
			} else {
				my ($seq_id, $seq) = DOMINO::fetching_range_seqs($keys, $coord_P1, $coord_P6, $fasta_msa_sub{$keys});
				$hash{$keys} = $seq;
		}}
		my $seq_name;
		
		## Check pairwise MSA
		my $valueReturned = DOMINO::check_marker_pairwise(\%hash, $$hash_parameters{'marker'}{'MCT'}[0], $$hash_parameters{'marker'}{'variable_positions_user_min'}[0], $$hash_parameters{'marker'}{'variable_positions_user_max'}[0], $$hash_parameters{'marker'}{'variable_divergence'}[0], $$hash_parameters{'marker'}{'polymorphism_user'}[0]);
		if ($valueReturned == 1) { ## if it is variable for each pairwise comparison
			## Get variable positions for the whole marker
			my $array_ref_returned = DOMINO::check_marker_ALL(\%hash, "Ref", $$hash_parameters{'marker'}{'missing_allowed'}[0], $$hash_parameters{'marker'}{'polymorphism_user'}[0], $$hash_parameters{'marker'}{'dnaSP_flag'}[0]);
			if ($array_ref_returned eq 'NO') { 
			} else {
				$marker_number++;
				my $msa_file = $dir."/".$contig_name[0]."_marker_".$marker_number.".fasta";
				open (MSA_OUT, ">$msa_file");
				foreach my $seqs (sort keys %hash) { print MSA_OUT ">".$seqs."\n".$hash{$seqs}."\n"; }
				close (MSA_OUT);
				my @arrayReturned = @$array_ref_returned; 	
				## 0: taxa join(::) ## 1: variable sites  ## 2: length of the marker ## 3: profile string ## 4: effective length

				my $msa_file_array = $dir."/".$contig_name[0]."_marker_".$marker_number."_ARRAY.txt";
				open (ARRAY, ">$msa_file_array"); print ARRAY ">".$contig_name[0]."_marker_".$marker_number."_ARRAY\n".$arrayReturned[3]."\n"; close (ARRAY);

				## Print into a file						
				my $percentage = ($arrayReturned[1]/$arrayReturned[4])*100;
				my $variation_sprintf = sprintf ("%.3f", $percentage);
				my $marker2write = $contig_name[0]."_#$marker_number";
				my @array2print = ($marker2write, $array[0].":".$array[1], $array[2].":".$array[3], $array[4].":".$array[5], $arrayReturned[2], $variation_sprintf, $array[6]); 
				my $string = join("\t", @array2print);
				open (OUT, ">>$file"); print OUT $string."\n"; close(OUT); last;		
}}}} close (MARKERS);
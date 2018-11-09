#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
}
my $domino_version ="DOMINO v1.1 ## Revised 07-11-2018";
my $scripts_path = $FindBin::Bin."/../";
my $samtools_path = $scripts_path."samtools-1.3.1/samtools";

##########################################################################################
##  This function generates a PILEUP and filters it										##
##	Jose Fco. Sanchez Herrero, 08/06/2015 jfsanchezherrero@ub.edu						##
##########################################################################################
my $sorted_bam = $ARGV[0]; 
my $contig_file = $ARGV[1];
my $reference_id = $ARGV[2];
my $taxa = $ARGV[3];
my $returned_outfolder = $ARGV[4];
my $path = $ARGV[5];
##########################################################################################

if (!@ARGV) {print "Missing arguments\n"; exit();}

## get general parameters
my $hash_parameters = DOMINO::get_parameters($path."/", "mapping");
#print Dumper \$hash_parameters;

my @temp_name = split ("\.sorted.bam", $sorted_bam);
my ($ID, @sam);
my $input_pileup = $temp_name[0].".profile";
my $pileup_command = $samtools_path." mpileup -f ".$contig_file." -o ".$input_pileup." ".$sorted_bam." 2> ".$$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0];
#&debugger_print("SAMTOOLS command: ".$pileup_command);
my $sytem_command_pileup = system ($pileup_command);
if ($sytem_command_pileup != 0) {
	DOMINO::printError("Exiting the script. Some error happened when calling SAMtools for generating the PILEUP for the file $contig_file...\n", $$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0]); DOMINO::dieNicely();
}
unless (-d $returned_outfolder) { mkdir $returned_outfolder, 0755; } 
#&debugger_print("Changing dir to $returned_outfolder");
##########################################################################################

## Get reference fasta information
my $reference_hash_fasta = DOMINO::readFASTA_hash($contig_file);
## Get DNA code
my %ambiguity_DNA_codes = %{ DOMINO::ambiguity_DNA_codes() };
##########################################################################################

## get ids
my $tmp = $input_pileup."_tmp_dir";
unless (-d $tmp) { mkdir $tmp, 0755 };
my $uniq_ids = $tmp."/uniq_ids.txt";
my $command = "awk \'{print \$1}\' $input_pileup \| sort \| uniq > $uniq_ids"; #print $command."\n";
system($command);
#print "\t\t- Filtering the PILEUPs generated for $reference_id vs. $taxa\n";
open ("IN", "<$uniq_ids");
while (<IN>) {
	my $line = $_; chomp $line;
	my $output = $tmp."/".$line.".profile";
	system("grep $line $input_pileup > $output");
	&read_pileup($output, $line);
}
close (IN);

sub read_pileup {

	my $file = $_[0];
	my $contig = $_[1];
	my @array_positions; my @fasta_positions;
	my $array_positions_ref = &initilize_contig($contig, $reference_hash_fasta);
	@array_positions = @$array_positions_ref;
	@fasta_positions = @array_positions;			

	open (PILEUP,"<$file"); while (<PILEUP>){

		print $_."\n";

		my $line = $_; chomp $line;
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;
		next if $line=~ m/\[REDUCE RESULT\]/;
		$line =~ s/\s/\t/g;
		my @pileup = split(/\t/,$line); ##	HKU61D301.clean.MID4_Nem_098_c8679	161	t	3	,,.	FA=
		my $contig = $pileup[0];
	
		my $pos_base; my $num_pos_base = $pileup[1]; my $num_pos_array = $num_pos_base -1;
	
		## Get array for each position
		if ($pileup[3] != 0) {
			my $read_base = $pileup[4]; my $ref_base = $pileup[2];
			my @base_record = split(//, $read_base);
			my (%posibilities, %polymorphism);
			my @base_parse;
			my $base_counter=0;
			for (my $i = 0; $i < scalar @base_record; $i++) {  
				if ($base_record[$i] =~ m/\^/) { $i++; next; } ## Starting to map a new read
				$base_record[$i]= uc($base_record[$i]);
	
				if ($base_record[$i] =~ m/A|G|C|T|a|g|c|t/) {
					$polymorphism{$base_record[$i]}++;
					$base_counter++;
					push (@base_parse, $base_record[$i]);
				} elsif ($base_record[$i] eq "." || $base_record[$i] eq ",") {
					$polymorphism{$ref_base}++;
					$base_counter++;
					push (@base_parse, $ref_base);
				} elsif ($base_record[$i] =~ m/(\+|\-)/) {# INDEL
					##	Example of INDEL: 
					##	HKU61D301.clean.MID4_Nem_098_c8679	280	g	2	,.+1T	I? 
					##  Contig_47_MID1	4	T	1	a+47cggatcgatcaaagtaagatatcatacttggaaggcaacatgcacgt	1	1	Position: 1
					my $length_base_record = length($read_base);
					my $indel;
					my $out_of_range = 0;
					my $indel_2 = $i+2;
					my $indel_3 = $i+3;
					if ($indel_2 >= $length_base_record) { $out_of_range = 1; }
					if ($indel_3 >= $length_base_record) { $out_of_range = 1; }
					if ($out_of_range == 0) {
						if ($base_record[$i+2] =~ /\d+/) { ## +/-12ATGAGACGATCC
							$indel = $base_record[$i+1].$base_record[$i+2];
							$i++; $i++;
						} elsif ($base_record[$i+3] =~ /\d+/) { ## +/-121ATGAGACGATCC...ATGAGACGATCC
							$indel = $base_record[$i+1].$base_record[$i+2].$base_record[$i+2];
							$i++; $i++;	$i++;
						} else { ## +/-1A
							$indel = $base_record[$i+1];
							$i++;
					}} else { ## +/-1A
						$indel = $base_record[$i+1];
						$i++;
					}
					$i = $i+$indel;
					next;
				} else { next; }
			}
			# Debug			
			#print Dumper @base_record; #print Dumper @base_parse; print Dumper %polymorphism;
			
			if ($ref_base eq "N") {
				$array_positions[$num_pos_array] = 'N'; ## Not informative enough  
				my @array_keys = keys %polymorphism;
				my $position;
				if (scalar @array_keys == 1) { ## a unique base is mapping
					$position = $array_keys[0];
				} else { ## get ambiguous code
					$position = DOMINO::get_amb_code(\%polymorphism);
				}
				$fasta_positions[$num_pos_array] = $position; 
				# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
				next;	
			}
			
			if ($base_counter == 1) { ## There is only base mapping...
				unless ($$hash_parameters{'mapping'}{'map_contig_files'}[0] || $$hash_parameters{'mapping'}{'LowCoverageData'}[0]) {
					$array_positions[$num_pos_array] = 'N'; ## Not informative enough
					$fasta_positions[$num_pos_array] = 'N'; ## Not informative enough
					# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
					next;
			}} ## Let it be informative if DOMINO simulations or mapping contig files				
			
			my @array_keys = keys %polymorphism;
			if (scalar @array_keys == 1) { ## a unique base is mapping
				if ($ambiguity_DNA_codes{$ref_base}) {
					my $flag = 0; my @bases;
					for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$ref_base}}; $j++) {
						foreach my $keys (sort keys %polymorphism) {
							if ($keys eq $ambiguity_DNA_codes{$ref_base}[$j]) { $flag = 1; }
					}}
					if ($flag == 1) { 
						## Contig1	4	M	10	cccccccccc	## type 4
						$array_positions[$num_pos_array] = '0'; 
						$fasta_positions[$num_pos_array] = $ref_base;
					} else { 
						## Contig1	5	R	10	cccccccccc	## type 5
						$array_positions[$num_pos_array] = '1';
						$fasta_positions[$num_pos_array] = $array_keys[0];											
				}} else {
					if ($base_parse[0] eq $ref_base) {  
						##Contig1	2	T	10	..........	## type 2
						$array_positions[$num_pos_array] = '0';
						$fasta_positions[$num_pos_array] = $ref_base;						
					} else {  
						## Contig1	3	T	10	cccccccccc	## type 3
						$array_positions[$num_pos_array] = '1';
						$fasta_positions[$num_pos_array] = $base_parse[0];
			}}} elsif (scalar @array_keys == 2) { ## Maybe true polymorphism or NGS error
				my ($value_return, $base_return) = &check_array($ref_base, \%polymorphism);
				$array_positions[$num_pos_array] = $value_return; 
				$fasta_positions[$num_pos_array] = $base_return;
			} elsif (scalar @array_keys == 3) { ## Something odd...
				my ($last_key, $small_key, $highest_value, $smallest_value);
				my $flag=0;				
				my @array = sort values %polymorphism;
				for (my $i=1; $i < scalar @array; $i++) {
					if ($array[0] == $array[$i]) { $flag++; }
				}
				if ($flag != 0) { ## There are two with the same frequence! DISCARD!
					$array_positions[$num_pos_array] = 'N';
					$fasta_positions[$num_pos_array] = 'N';
				} else {
					## Get lowest value, discard it and check the rest.
					foreach my $keys (sort keys %polymorphism) {
						if ($polymorphism{$keys} eq $array[0]) {
							delete $polymorphism{$keys}; last;
					}}
					my ($value_return, $base_return) = &check_array($ref_base, \%polymorphism);
					$array_positions[$num_pos_array] = $value_return; 
					$fasta_positions[$num_pos_array] = $base_return;
			}} elsif (scalar @array_keys > 3) { 
				$array_positions[$num_pos_array] = 'N';
				$fasta_positions[$num_pos_array] = 'N';
		}} else { ## No base is mapping this reference
			##Contig1	1	A	0	## type 1
			$array_positions[$num_pos_array] = 'N';
			$fasta_positions[$num_pos_array] = 'N';
		}
		# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
		
	}
	close(PILEUP);

	&print_coordinates(\@array_positions, \$contig, $reference_id, $returned_outfolder);
	&print_fasta_coordinates(\@fasta_positions, \$contig, $reference_id, $returned_outfolder); ## Print array into file $previous_contig
	sleep(2);
}

sub initilize_contig {		
	my $name_contig = $_[0];
	my $reference_hash_fasta = $_[1];
	
	## Initialize new contig
	my @array_positions_sub;
	if (${$reference_hash_fasta}{$name_contig}) {
		my $tmp_size = length(${$reference_hash_fasta}{$name_contig});
		@array_positions_sub = ("-") x $tmp_size;
	}
	return \@array_positions_sub;		
}

sub print_coordinates {
	## Print array into file $previous_contig
	my $coord_array_ref = $_[0];
	my $contig_name = $_[1];
	my $ref_id = $_[2];
	my $dir_tmp = $_[3];

	my @coord_array = @$coord_array_ref;
	my $seq_contig = join "", @coord_array;
	my $var_sites = $seq_contig =~ tr/1/1/;

	if ($var_sites != 0) {
		my $array_file = $dir_tmp."/".$$contig_name."_ARRAY.txt";
		open (FH, ">$array_file"); print FH ">$$contig_name\n$seq_contig\n"; close(FH);	
	}
}

sub print_fasta_coordinates {
	## Print array into file $previous_contig
	my $coord_array_ref = $_[0];
	my $contig_name = $_[1];
	my $ref_id = $_[2];
	my $dir_tmp = $_[3];
	
	my @coord_array_sub = @$coord_array_ref;
	my $seq_contig = join "", @coord_array_sub;
	my $array_file = $dir_tmp."/".$$contig_name."_sequence.fasta";
	open (FH, ">$array_file"); print FH ">$$contig_name\n$seq_contig\n"; close(FH);	
}

sub check_array {
	my $ref_base = $_[0];
	my $ref_poly_hash = $_[1];

	my %polymorphism = %$ref_poly_hash;
	my $total=0;
	my $polymorphism_user;
	if ($$hash_parameters{'mapping'}{'polymorphism'}) { $polymorphism_user++; }
	my ($last_key, $highest_value, $smallest_value);
	foreach my $keys (sort keys %polymorphism) {
		$total += $polymorphism{$keys};
		if (!$highest_value) {
			$highest_value = $polymorphism{$keys};
			$smallest_value = $polymorphism{$keys};
			$last_key = $keys;							
		} elsif ($polymorphism{$keys} > $highest_value) {
			$highest_value = $polymorphism{$keys};
			$last_key = $keys;							
		} elsif ($polymorphism{$keys} < $smallest_value) {
			$smallest_value = $polymorphism{$keys};
	}}
	
	print $ref_base."\t".$total."\t".$smallest_value."\t".$highest_value."\n";
	print Dumper $ref_poly_hash;	

	if ($total > 170) { ## will be out of range
	 	if ($$hash_parameters{'mapping'}{'noDiscard'}[0]) { # no discard contigs by coverage
			## make it equivalent to max:169
			my $new_highest_value = int(($highest_value/$total)*100);
			print "Total = $total\t$highest_value\t$new_highest_value\t169\n";
			$highest_value = $new_highest_value;
	 		$total = 169;
	 	} else { return ("N","N"); }
	}
	## Check wether there are more than 8 positions mapping
	## and if not if there at least for the smallest two
	## bases.
	
	my $prob = DOMINO::binomial($highest_value, $total, 0.5);
	my $significance_level = 1 - 0.95;
	if ($ambiguity_DNA_codes{$ref_base}) {
		## There is an ambiguous code in the reference
		my %ref_base_match;
		my $flag = 0;
		for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$ref_base}}; $j++) {
			foreach my $keys (sort keys %polymorphism) {
				if ($keys eq $ambiguity_DNA_codes{$ref_base}[$j]) {
					$flag++; $ref_base_match{$ambiguity_DNA_codes{$ref_base}[$j]}++;
		}}}
		if ($flag == 2) { 
			## Only if both possibilities are the same as the ambiguous code
			if ($prob < $significance_level) { 
				return ('0', $last_key); ##Contig1	10	R	10	AAAAAAAAAG
			} else { 
				## Contig1	9	R	10	AAAAAAGGGG
				## Contig1	9	R	10	AAGG
				## Contig1	9	R	10	AG
				if ($polymorphism_user) {
					return ('1', $ref_base);
				} else {
					return ('0', $ref_base);
		}}} else {
			if ($flag == 1) { ## Only one type of read is within the ambiguous code
				if ($prob < $significance_level) { 
					if ($ref_base_match{$last_key}) { 
						return ('0', $last_key); ## Contig1	11	R	10	AAAAAAAAAC
					} else { 
						return ('1', $last_key); ## Contig1	13	R	10	TTTTTTTTTTA
				}} else { 
					## Contig1	12	R	10	AAAAACCCCC
					## Contig1	12	R	10	AACC
					## Contig1	12	R	10	AC
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					if ($polymorphism_user) {							
						return ('1', $pos);
					} else {
						return ('0', $pos);
			}}} else { ## None of the reads mapping are similar to the ambiguous reference
				if ($prob < $significance_level) { 	
					return ('1', $last_key); ## Contig1	13	R	10	TTTTTTTTTTC
				} else { 
					## Contig1	13	R	10	TTTTTTCCCCC
					## Contig1	13	R	10	TTCC
					## Contig1	13	R	10	TC
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					return ('1', $pos);
	}}}} else { 
		## There is no ambiguity code in the reference base
		my $flag = 0;
		if ($polymorphism{$ref_base}) { ## if any mapping reads are the same as reference
			if ($total >= 8) {
				if ($prob < $significance_level) {
					if ($last_key eq $ref_base) { 
						return ('0', $last_key); ## Contig1	6	T	10	.........a
					} else { 
						return ('1', $last_key); ## Contig1	7	T	10	ccccccccc.
				}} else { 
					## Contig1	8	T	10	....cc..cc
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					if ($polymorphism_user) {
						return ('1', $pos);
					} else {
						return ('0', $pos);
			}}} else { ## less than 8 reads mapping here
				if ($smallest_value >= 2) {
					## Contig1	8	T	5	...cc
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					if ($polymorphism_user) {
						return ('1', $pos);
					} else {
						return ('0', $pos);
				}} else {
					if ($highest_value == $smallest_value) {
						return ('N', 'N'); ## Contig1	8	T	2	.c
					} elsif ($last_key eq $ref_base) { 
						return ('0', $ref_base); ## Contig1	8	T	4	...c
					} else {
						return ('1', $last_key); ## Contig1	8	T	4	.c
		}}}} else {
			## None of the reads are similar to reference
			if ($total >= 8) { 
				if ($prob < $significance_level) { 	
					return ('1', $last_key); ## Contig1	15	T	10	ccccccccccca
				} else { 
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					return ('1', $pos); ## Contig1	15	T	10	ccccccgggggg
			}} else {
				if ($smallest_value >= 2) {
					my $pos = DOMINO::get_amb_code(\%polymorphism);
					return ('1', $pos); ## Contig1	8	T	5	GGGcc
				} else {
					if ($highest_value == $smallest_value) {
						return ('1', 'N'); ## Contig1	8	T	2	Gc
					} else {
						return ('1', $last_key); ## Contig1	8	T	3	GGc
}}}}}}

__END__
################################################################################################################################################
PILEUP Explanation

	Each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities. 
	At the read base column, 
	+ a dot stands for a match to the reference base on the forward strand, 
	+ a comma for a match on the reverse strand, 
	+ `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. 
	+ A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. 
		The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:
		seq2 156 A 11  .$......+2AG.+2AG.+2AGGG    <975;:<<<<<

		Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. Here is an exmaple of a 4bp deletions from the reference, supported by two reads:
		seq3 200 A 20 ,,,,,..,.-4CACC.-4CACC....,.,,.^~. ==<<<<<<<<<<<::<;2<<
	
	+ A symbol `^' marks the start of a read segment which is a contiguous subsequence on the read separated by `N/S/H' CIGAR operations. 
	+ The ASCII of the character following `^' minus 33 gives the mapping quality. 
	+ A symbol `$' marks the end of a read segment. 
	
	Start and end markers of a read are largely inspired by Phil Green's CALF format. These markers make it possible to reconstruct the read sequences from pileup.
	SAMtools can optionally append mapping qualities to each line of the output. This makes the output much larger, but is necessary when a subset of sites are selected.
	
	Contig_25_MID3  489     C       19      ,.,.,.,...t,.,.,,..     AC9F9F;FIH86I;351II
	Contig_25_MID3  490     A       18      ,.,.,.,...,,.,,,..      AC7F9F:FIH76I851II
	Contig_25_MID3  491     G       18      cCcCcCc.CCccCccc.C      <C;F;F<FIF1:I5<6II
	Contig_25_MID3  492     G       18      ,.,.,.,...,,.,,,..      4E5E7F:FIF.9I3:9II
	Contig_25_MID3  493     A       19      ,.,.,.,...,.,,,..^$.^$. 4B5E7F9FIF9H3:6II55
	Contig_25_MID3  494     A       20      ,.,.,.,...,,.,,,....    4B6E7F9FIF.9H371II55
	Contig_25_MID3  495     C       20      tTtTtTtTGTttTtttTTTT    9B9F;F;FIF19H771II55
	Contig_25_MID3  496     A       20      ,.,.,.,...,,.,,,....    9A;F9B=FIF16F540II88
	
	Contig1	1	A	0	## type 1

	Contig1	2	T	10	..........	## type 2

	Contig1	3	T	10	cccccccccc	## type 3

	Contig1	4	M	10	cccccccccc	## type 4

	Contig1	5	R	10	cccccccccc	## type 5

	Contig1	6	T	10	.........a	## type 6

	Contig1	7	T	10	ccccccccc.	## type 7

	Contig1	8	T	10	....cc..cc	## type 8

	Contig1	9	T	10	ccccccccccca ## type 9

	Contig1	10	T	10	ccccccgggggg ## type 10

	Contig1	11	R	10	AAAAAAGGGG	## type 11

	Contig1	12	R	10	AAAAAAAAAG	## type 12

	Contig1	13	R	10	AAAAAAAAAC	## type 13

	Contig1	14	R	10	AAAAACCCCC	## type 14

	Contig1	15	R	10	TTTTTTTTTTA	## type 15

	Contig1	16	R	10	TTTTTTTTTTC	## type 16

	Contig1	17	R	10	TTTTTTCCCCC	## type 17

	Contig1	18	T	10	ACAAGGGGGG	## type 18

	Contig1	19	T	10	ACCCAAGGGGGG	## type 19

	Contig1	20	R	10	ACAAAAGGGGGG	## type 20

	Contig1	21	R	10	ACCCAAGGGGGG	## type 21

	Contig1	22	R	10	ACCCAAGGGGGGTTTTTT	## type 22

	Contig1	23	R	10	GCCCCCCCCCCCTTTTTT	## type 23

#####################################################################################################################

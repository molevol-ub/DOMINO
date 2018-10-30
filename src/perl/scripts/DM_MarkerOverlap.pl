#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
}

## Overlaps and maximizes domino markers obtained
my $file = $ARGV[0]; 				## markers_shared 
my $mergeArray_file = $ARGV[1];
my $file2return = $ARGV[2];
my $path = $ARGV[3];

my $domino_version ="DOMINO v1.1 ## Revised 30-10-2018";
my $hash_parameters = DOMINO::get_parameters($path."/");

#&debugger_print("Checking file $file [DM_MarkerScan: check_overlapping_markers]");

my $contig_id; my %tmp_hash; my $marker_counter_tmp = 0; my @sequences;
open (FILE, $file); while (<FILE>) {
	my $line = $_; chomp $line; 
	$line =~ s/ /\t/;
	my @array_lines = split ("\t", $line);		
	$contig_id = $array_lines[0];
	my @a = split(":", $array_lines[1]); ## conserved region1 coord
	my @b = split(":", $array_lines[2]); ## variable region coord
	my @c = split(":", $array_lines[3]); ## conserved region2 coord
	my @coordinates = ($a[0],$a[1],$b[0],$b[1],$c[0],$c[1]);
	my $string = join(",", @coordinates);
	my $taxa;
	if ($array_lines[4]) { $taxa = $array_lines[4];
	} else { $taxa = $$hash_parameters{'marker'}{'taxa_string'}[0];  }		
	push (@{ $tmp_hash{ $contig_id }{ $taxa }{ $a[0] } }, $string); ## Keep record of taxa and coordinates
} 
close(FILE);

# Debug 
#print Dumper \%tmp_hash;
	
my %coord_seen;
foreach my $contig (sort keys %tmp_hash) {
	foreach my $taxa (sort keys %{ $tmp_hash{$contig} }) {		
		foreach my $marker (keys %{ $tmp_hash{$contig}{$taxa} }) {			
			if ($coord_seen{$contig}{$taxa}{$marker}) {next;}
			my $bool = 1;
			my ($counter, $bad_counter) = 0;
			while ($bool) {
				$counter += $$hash_parameters{'marker'}{'SLIDING_user'}[0];
				my $new_coord = $marker + $counter;
				if ($coord_seen{$contig}{$taxa}{$new_coord}) {next;}
				if ($tmp_hash{$contig}{$taxa}{$new_coord}) {
					push (@{ $tmp_hash{$contig}{$taxa}{$marker} }, @{ $tmp_hash{$contig}{$taxa}{$new_coord} });
					$coord_seen{$contig}{$taxa}{$new_coord}++;
					$tmp_hash{$contig}{$taxa}{$new_coord} = 1;
				} else {
					$bad_counter++;
					if ($bad_counter > 3) { ## We would consider the same marker if overlapping 3 SLIDING_user!!
						($bool,$counter,$bad_counter) = 0;
}}}}}}
my %tmp_coord;
foreach my $contig (sort keys %tmp_hash) {
	foreach my $taxa (keys %{ $tmp_hash{$contig} }) {
		foreach my $marker (keys %{ $tmp_hash{$contig}{$taxa} }) {
			if ($tmp_hash{$contig}{$taxa}{$marker} == 1) {next;}
			my @sort_uniq = do { my %seen; grep { !$seen{$_}++ } @{ $tmp_hash{$contig}{$taxa}{$marker} } };
			push (@{ $tmp_coord{$contig}{$taxa}{$marker} }, @sort_uniq);
}}}
undef %coord_seen; undef %tmp_hash; ## release RAM

## Set range values
my $range = $$hash_parameters{'marker'}{'window_size_VARS_max'}[0] - $$hash_parameters{'marker'}{'window_size_VARS_min'}[0]; my @length;
if ($range >= 500) { 	 
	@length = (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500); 
} elsif ($range < 500) { 
	@length = (50, 100, 200, 300, 400, 500);
}

# Debug print Dumper \%tmp_coord
my $hash_ref = DOMINO::readFASTA_hash($mergeArray_file); 
open (OUT, ">$file2return"); 
foreach my $contig (sort keys %tmp_coord) {
	foreach my $taxa (keys %{ $tmp_coord{$contig} }) {
		foreach my $keys_markers (keys %{ $tmp_coord{$contig}{$taxa} }) {		
			my @array_coordinates = @{ $tmp_coord{$contig}{$taxa}{$keys_markers} };
			my (@good_ones, %hash2print, %marker_seen); 
			for (my $i=0; $i < scalar @array_coordinates; $i++) {
				my $coordinate = $array_coordinates[$i];
				my @array_coord = split (",", $coordinate);
				for (my $j = (scalar @array_coordinates) - 1; $j >= 0; $j--) {
					my $coordinate2 = $array_coordinates[$j];
					my @array_coord2 = split (",", $coordinate2);
					my @array2check = ($array_coord[0], $array_coord2[1], $array_coord2[2], $array_coord2[3], $array_coord2[4], $array_coord2[5]);
					my $string = join(":", @array2check).":".$taxa;
					#print $$hash_ref{$contig}."\n";
					if ($marker_seen{$string}) {next;}
					my $result = &check_given_marker(\@array2check, $$hash_ref{$contig});						
				
					if ($result ne 1) {
						my $id;
						for (my $j = 0; $j < scalar @length; $j++) {
							if ($result <= $length[$j] ) { 		$id = "$length[$j - 1] - $length[$j]"; last;
							} elsif ($result > $length[-1]) {	$id = "bigger"; last;
						}}
						push (@{ $hash2print{$array_coord[0]}{$id}}, $string);
						# print $keys_markers."\t".$array_coord[0]."\t".$array_coord2[5]."\t".$result."\t".$id."\t".$string."\n"; 
						$marker_seen{$string}++;
			}}}
			foreach my $keys (keys %hash2print) {
				foreach my $lent (sort keys %{ $hash2print{$keys} }) {
					my @array = @{ $hash2print{$keys}{$lent} };					
					for (my $i=0; $i < scalar @array; $i++) { 
						print OUT $contig."##".$array[$i]."\n";
					} print OUT "//\n";
}}}}}
close(OUT); 

DOMINO::print_success_Step("marker_overlap");


sub check_given_marker {
	
	my $array_Coord = $_[0];
	my $dna_seq = $_[1];
	my $total_length_sub = length($dna_seq);
	
	my @coordinates = @{ $array_Coord };
	my $coord_P1 = $coordinates[0]; my $coord_P2 = $coordinates[1];
	my $coord_P3 = $coordinates[2]; my $coord_P4 = $coordinates[3];
	my $coord_P5 = $coordinates[4]; my $coord_P6 = $coordinates[5];
	
	# Conserved
	my $length_string_P1_P2 = $coord_P2 - $coord_P1;
	my $string_P1_P2 = substr($dna_seq, $coord_P1, $length_string_P1_P2);
	my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//; ## Conserved 
	
	if ($length_string_P1_P2 < $$hash_parameters{'marker'}{'window_size_CONS_min'}[0]) {return 1;}
	if ($length_string_P1_P2 > $$hash_parameters{'marker'}{'window_size_CONS_max'}[0]) { 
		#print Dumper $array_Coord; print "ERROR length_string_P1_P2! $length_string_P1_P2 > $window_size_CONS_max\n"; return 1;
		if ($$hash_parameters{'marker'}{'window_size_CONS_max'}[0] != $$hash_parameters{'marker'}{'window_size_CONS_min'}[0]) { return 1; } ## If user specifies a range..
	}
	if ($count_string_P1_P2 > $$hash_parameters{'marker'}{'window_var_CONS'}[0]) { 
		#print Dumper $array_Coord; print "ERROR count_string_P1_P2! $count_string_P1_P2 > $window_var_CONS\n";
		return 1;
	} 
		
	# Variable
	my $length_string_P3_P4 = $coord_P4 - $coord_P3;
	if ($length_string_P3_P4 > $$hash_parameters{'marker'}{'window_size_VARS_max'}[0]) { 
		#print Dumper $array_Coord; print "ERROR! length_string_P3_P4 $length_string_P3_P4 > $window_size_VARS_max\n"; 
		return 1;
	}
	my $string_P3_P4 = substr($dna_seq, $coord_P3, $length_string_P3_P4);
	my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//; ## Variable
	
	if ($$hash_parameters{'marker'}{'variable_divergence'}[0]) {
		# If a minimun divergence, get the expected variable sites for the length
		my $expected_var_sites = int($$hash_parameters{'marker'}{'variable_divergence'}[0] * $total_length_sub);
		unless ($count_string_P3_P4 >= $expected_var_sites) { return 1; }			
	}

	# Conserved
	my $length_string_P5_P6 = $coord_P6 - $coord_P5;
	my $string_P5_P6 = substr($dna_seq, $coord_P5, $length_string_P5_P6);
	if (!$string_P5_P6) {
		#print Dumper $array_Coord; #print "Length: $length_string_P5_P6\n"; #print $dna_seq;
		return 1;
	}
	my $count_string_P5_P6 = $string_P5_P6  =~ tr/1//; ## Conserved
	if ($length_string_P5_P6 < $$hash_parameters{'marker'}{'window_size_CONS_min'}[0]) { return 1;}
	if ($length_string_P5_P6 > $$hash_parameters{'marker'}{'window_size_CONS_max'}[0]) { 
		#print Dumper $array_Coord; print "ERROR! length_string_P5_P6 $length_string_P5_P6 > $window_size_CONS_max\n"; return 1;
		if ($$hash_parameters{'marker'}{'window_size_CONS_max'}[0] != $$hash_parameters{'marker'}{'window_size_CONS_min'}[0]) { return 1; } ## If user specifies a range..
	}
	if ($count_string_P5_P6 > $$hash_parameters{'marker'}{'window_var_CONS'}[0]) {
		#print Dumper $array_Coord; print "ERROR! count_string_P5_P6 $count_string_P5_P6 > $window_var_CONS\n";
		return 1;
	}
	my $total_length = $coord_P6 - $coord_P1;
	return $total_length;
}


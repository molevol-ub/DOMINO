#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
	require File::Copy;
	require Parallel::ForkManager;
	require File::Path; use File::Path qw(remove_tree);
}
#################################################################
my $domino_version ="DOMINO v1.1 ## Revised 07-11-2018";
my $scripts_path = $FindBin::Bin."/../";
my $domino_Scripts = $FindBin::Bin;
my $BLAST = $scripts_path."NCBI_BLAST/";

#################################################################

## Arguments
my $path = $ARGV[0];
my $step_time = $ARGV[1];
#################################################################

#################################################################
## Get general parameters and files
my $hash_parameters = DOMINO::get_parameters($path."/", "markers");
my %domino_cluster_files = %{ DOMINO::get_DOMINO_files($path."/","markers") };
	#print Dumper $hash_parameters;
	#print Dumper \%domino_cluster_files;

my $marker_dirname = $$hash_parameters{'marker'}{'folder'}[0];
my $num_proc_user = $$hash_parameters{'marker'}{'cpu'}[0];

my %new_domino_files;
my %new_dowmino_param;
#################################################################

#############################################################
## Clusterize markers using different taxa as reference ##
#############################################################
my $blast_dir = $marker_dirname."/clustering";
mkdir $blast_dir, 0755; chdir $blast_dir; #&debugger_print("Changing dir to $blast_dir");

# read dir
my $files_dir_ref = DOMINO::readDir($marker_dirname);
my @markers_folders = @$files_dir_ref;

# output files
my $all_coordinates_file = $blast_dir."/all_coordinates.fasta"; open (ALL_coordinates, ">$all_coordinates_file");
my $all_contigs_file = $blast_dir."/all_contigs.fasta"; open (ALL_CONTIGS, ">$all_contigs_file"); 
print "+ Merging different DOMINO markers according to the taxa of reference...\n";

# get sequences
my @CoordMarker;
foreach my $keys (sort keys %domino_cluster_files) {
	if (!$domino_cluster_files{$keys}{'coordinates'}) {next;}
	if ($domino_cluster_files{$keys}{'taxa'}) {
		push (@CoordMarker, $domino_cluster_files{$keys}{'coordinates'}[0]);
		my $coordinate_file = $domino_cluster_files{$keys}{'markers'}[0];
		if (-e -r -s $coordinate_file) {
			my $hash = DOMINO::readFASTA_hash($coordinate_file);
			foreach my $seq (keys %{$hash}) {
				print ALL_coordinates ">".$seq."_taxa_".$keys."\n".$$hash{$seq}."\n";
		}}
		my $contig_file = $domino_cluster_files{$keys}{'CONTIGS'}[0];
		if (-e -r -s $contig_file) {
			my $size_file = DOMINO::get_size($contig_file);
			open (FH, $contig_file);
			my $chunk; read(FH, $chunk, $size_file); 
			print ALL_CONTIGS $chunk;
}}} close(ALL_CONTIGS); close(ALL_coordinates);

## Use BLAST for clustering sequences
print "+ Generate a BLAST database...\n"; 
my ($blast_DB, $blast_DB_message) = DOMINO::makeblastdb($all_coordinates_file, $BLAST, $$hash_parameters{"mapping"}{"mapping_markers_errors_details"}[0]);
#&debugger_print($blast_DB_message);
if ($blast_DB eq "1") {
	my $msg="Early termination of the DOMINO Marker Scan...\n\nPlease note that DOMINO could not find any markers for the parameters provided. Please Re-Run DOMINO using other parameters\n\n\n"; 
	print $msg; DOMINO::printError($msg, $$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0]); exit();
} 

## Parallelize BLAST
print "+ BLAST search now...\n";
my $size_file_BLAST = DOMINO::get_size($all_coordinates_file);
my $parts2split_BLAST = int($size_file_BLAST/$num_proc_user);
my $fasta_files_split_BLAST = DOMINO::fasta_file_splitter($all_coordinates_file, $parts2split_BLAST, "fasta", $blast_dir);

## SEND THREAD 
my $pm_BLAST = new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
for (my $j = 0; $j < scalar @{ $fasta_files_split_BLAST }; $j++) {
	my $pid = $pm_BLAST->start($j) and next;
	my $total = scalar @{ $fasta_files_split_BLAST };
	my $int = $j+1;
	print "\t- Sending BLASTn command for clustering [$int/$total]...\n";
	my $blast_search_tmp = "blast_search.txt_tmp".$j; 
	my ($blastn, $blastn_message) = DOMINO::blastn($all_coordinates_file, $blast_DB, $blast_search_tmp, $BLAST);
	#&debugger_print($blastn_message);
	$pm_BLAST->finish(); # finish for each contig
} 
$pm_BLAST->wait_all_children; #each marker
my $blast_search = "blast_search.txt";
system("cat blast_search.txt_tmp* >> $blast_search; rm blast_search.txt_tmp*");

## Filter BLAST results
print "+ Filtering BLAST search now...\n";
my $contig_length_Ref = DOMINO::readFASTA_hash($all_coordinates_file);
my (%markers2keep, @markers_seen, %clusterized_contigs_keep);
my $first_hit = 0;
my $aln_overlapped = $$hash_parameters{'marker'}{'window_size_VARS_max'}[0] + $$hash_parameters{'marker'}{'window_size_CONS_max'}[0];
print "+ Clustering markers with > $aln_overlapped bp overlapped...\n";
open (BLAST, $blast_search); while (<BLAST>) {
	my $line = $_;
	chomp $line;
	my @array = split("\t", $line);
	my $query = $array[0]; my $subject = $array[1];
	if ($query eq $subject) { next;} # same sequence
	my $query_string; my $subject_taxa;
	my $subject_string; my $taxa_query;
	if ($query =~ /(.*)\_taxa\_(.*)/) { $query_string = $1;   $taxa_query = $2; }
	if ($subject =~ /(.*)\_taxa\_(.*)/){$subject_string = $1; $subject_taxa = $2; }
	if ($subject_taxa eq $taxa_query) { next; } ## if same taxa
	if ($array[10] < 1e-20 && $array[3] > $aln_overlapped) {    ## how to clusterize...
		if ($query =~ /(.*)\_coord\_(.*)\_taxa\_(.*)/) { $clusterized_contigs_keep{$1}++;}
		if ($first_hit == 0) {
			$first_hit++;
			push( @{$markers2keep{$query_string} }, $subject_string); 
			push (@markers_seen, $subject_string);
		} else {
			my $flag_this = 0;
			foreach my $keys (sort keys %markers2keep) {
				if (grep /.*$subject_string.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
				if (grep /.*$query_string.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
			}
		if ($flag_this == 0) { 	push(@{$markers2keep{$query_string}}, $subject_string);  push (@markers_seen, $subject_string); }
}}} close (BLAST);
foreach my $keys (sort keys %$contig_length_Ref) { 
	if ($keys =~ /((.*)\_coord\_(.*))\_taxa\_(.*)/) { 
		unless (grep /$1/, @markers_seen) { 
			$markers2keep{$1}++; $clusterized_contigs_keep{$2}++;
}}}

## Printing definitely Results
my $definitely_results_dirname = $marker_dirname."/DOMINO_markers_Results";
mkdir $definitely_results_dirname, 0755; chdir $definitely_results_dirname; #&debugger_print("Changing dir to $definitely_results_dirname");
print "+ Getting the information ready for the visualization of results...\n";
my $contig_def_results_sequences = "DM_contigs.fasta";
open (CONTIGS_DEF, ">$contig_def_results_sequences");
my $hash_contigs_ref = DOMINO::readFASTA_hash($all_contigs_file);
my %hash_contigs = %$hash_contigs_ref;
foreach my $keys (sort keys %hash_contigs) {
if ($clusterized_contigs_keep{$keys}) { print CONTIGS_DEF ">".$keys."\n".$hash_contigs{$keys}."\n";	}}
close(CONTIGS_DEF);

my %hash4markers2keep;
print "+ Filtering coordinates files...\n\n";
for (my $h = 0; $h < scalar @CoordMarker; $h++) {
	open (tmp_COOR, $CoordMarker[$h]);
	while (<tmp_COOR>) {
		my $line = $_;
		chomp $line;
		next if $line=~ m/^\s*$/o;
		if ($line =~ /.*Vari/) {next;}
		my @array_split = split("\t",$line);
		my @contig_split = split("\_#",$array_split[0]);
		if ($clusterized_contigs_keep{$contig_split[0]}) {
			my @array_start = split(":", $array_split[1]);
			my @array_end = split(":", $array_split[3]);
			my $coord_string = $array_start[0]."_".$array_end[1];
			my $coord2check = $contig_split[0]."_coord_".$coord_string;
			if ($markers2keep{$coord2check}) {
				$hash4markers2keep{$contig_split[0]}{$coord_string} = $line;
}}} close(tmp_COOR); }

my $coordinates_def_results = $definitely_results_dirname."/DM_markers-coordinates.txt";
open (COOR, ">$coordinates_def_results");
print COOR "Contig\t\tConserved_Region\tVariable_Region\tConserved_Region\tMapping_Taxa\Length\tDivergence\n";
my %rename;
foreach my $contigs (sort keys %hash4markers2keep) {
	my $counter = 0;
	foreach my $markers (sort keys %{ $hash4markers2keep{$contigs} }) {
		$counter++;
		my $string2change = $hash4markers2keep{$contigs}{$markers};
		my $stringchanged;
		my @array = split("\t", $string2change);
		my $marker2search = $array[0];
		$marker2search =~ s/_#/_marker_/;
		my $tmp = $contigs."_#".$counter;
		$array[0] = $tmp;
		$rename{$marker2search} = $tmp;
		for (my $i=0; $i < scalar @array; $i++) { $stringchanged .= $array[$i]."\t"; }			
		print COOR $stringchanged."\n";
} print COOR "\n"; }
close(COOR); &time_log(); print "\n";

## Get MSA markers
my $markers_msa_folder = $definitely_results_dirname."/MSA_markers"; mkdir $markers_msa_folder, 0755;
foreach my $keys (sort keys %domino_cluster_files) {
	if (!$domino_cluster_files{$keys}{'coordinates'}) {next;}
	if ($domino_cluster_files{$keys}{'taxa'}) {
		my $MSA_markers_each_Taxa = $domino_cluster_files{$keys}{'MSA_markers'}[0];
		my $array_Ref = DOMINO::readDir($MSA_markers_each_Taxa);
		for (my $i=0; $i < scalar @{ $array_Ref }; $i++) {
			my @name = split(".fasta", $$array_Ref[$i]);
			if ( $rename{$name[0]} ) {
				File::Copy::copy($MSA_markers_each_Taxa."/".$$array_Ref[$i], $markers_msa_folder."/".$rename{$name[0]}.".fasta");	
}}}}

## Print excel for clusterized results
print "+ Generating an Excel file for DOMINO markers coordinates...\n";

my $domino_Scripts_excel = $domino_Scripts."/DM_PrintExcel.pl";
my $command = "perl $domino_Scripts_excel $path $coordinates_def_results $definitely_results_dirname";
print "\n[ System Call: ".$command." ]\n";
system($command);
################################################################################################################################################


###########################
####### SUBROUTINES #######
###########################

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}
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
my $domino_version ="DOMINO v1.1 ## Revised 30-10-2018";
my $scripts_path = $FindBin::Bin."/../";
my $domino_Scripts = $scripts_path."scripts";
#################################################################

## Arguments
my $path = $ARGV[0];
my $step_time = $ARGV[1];
#################################################################

#################################################################
## Get general parameters and files
my $hash_parameters = DOMINO::get_parameters($path."/");
my %domino_align_files = %{ DOMINO::get_DOMINO_files($path."/") };
#print Dumper $hash_parameters;
#print Dumper $domino_align_files_Ref;

my $align_dirname = $$hash_parameters{'mapping'}{'folder'}[0];
my $num_proc_user = $$hash_parameters{'mapping'}{'cpu'}[0];
my $msa_dirname = $align_dirname."/MSA_files";
my $profile_dir = $marker_dirname."/PROFILES";	mkdir $profile_dir, 0755;
my $msa_dir_tmp = $marker_dirname."/MSA_markers_tmp"; mkdir $msa_dir_tmp, 0755;
my $dir_Dump_file = $marker_dirname."/DUMP_files"; mkdir $dir_Dump_file, 0755;
my $msa_dir = $marker_dirname."/MSA_markers"; mkdir $msa_dir, 0755;

my %new_domino_files;
#################################################################

my $array_files_fasta_msa_ref = DOMINO::readDir($domino_files{'MSA_files'}{'folder'}[0]);
my @array_files_fasta_msa = @$array_files_fasta_msa_ref;
my $total_files = scalar @array_files_fasta_msa;
print "+ Checking files in folder generated...\n";
my $counter = 0;	
my $pm_MARKER_MSA_files =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
$pm_MARKER_MSA_files->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
$pm_MARKER_MSA_files->run_on_start( sub { my ($pid,$ident)=@_; } );
for (my $i=0; $i < scalar @array_files_fasta_msa; $i++) {
	if ($array_files_fasta_msa[$i] eq "." || $array_files_fasta_msa[$i] eq ".." || $array_files_fasta_msa[$i] eq ".DS_Store"  ) {next;}
	$counter++; 
	if ($total_files > 100) {
		my $perc = sprintf( "%.3f", ( $counter/$total_files )*100 );
		print "\t- Checking each sequence: [ $perc % ]...\r";
	} else { print "\t- Checking each sequence: [$counter/$total_files]...\r";}	

	my $pid = $pm_MARKER_MSA_files->start($i) and next;
	my %domino_files_msa;
	my $file_path = $domino_files{'MSA_files'}{'folder'}[0]."/".$array_files_fasta_msa[$i];
	unless (-f -e -r -s $file_path) {
		DOMINO::printError("File $file_path is not readable or empty. It would be skipped....\nPlease discarded from the folder...\n");
	} else {
		#print "\nChecking now: $array_files_fasta_msa[$i]\n";
		my ($region_id, $string2print_ref, $hash_ref_msa);
		if ($stacks_file) {
			#### STACKS
			$hash_ref_msa = DOMINO::readFASTA_hash($file_path);
			my $filename = $array_files_fasta_msa[$i];	
			my %hash = %$hash_ref_msa; my ( %hash2fill, %hash2return);
			foreach my $keys (sort keys %hash) {
			if ($keys =~ /(CLocus\_\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
				$region_id = $1; my $sample = $2;
				my $sequence = $hash{$keys};
				push (@{ $hash2fill{$sample} }, $keys.":::".$sequence);
			}}
			foreach my $keys (sort keys %hash2fill) {
				my $alleles = scalar @{ $hash2fill{$keys} };
				if ($alleles == 1) {
					my @array = split(":::", $hash2fill{$keys}[0]);
					$hash2return{$keys} = $array[1];					
				} elsif ($alleles == 2){					
					my @array = split(":::", $hash2fill{$keys}[0]);
					my @allele1 = split("", $array[1]);
					my @array2 = split(":::", $hash2fill{$keys}[1]);
					my @allele2 = split("", $array2[1]);
					my @string;
					for (my $i=0; $i < scalar @allele1; $i++) {
						my @tmp; push(@tmp, $allele1[$i]); push(@tmp, $allele2[$i]);
						my @tmp_uniq_sort = do { my %seen; grep { !$seen{$_}++ } @tmp };
						if (scalar @tmp_uniq_sort == 1) {
							if ($tmp_uniq_sort[0] eq 'N') {  		push (@string, 'N');
							} elsif ($tmp_uniq_sort[0] eq '-') { 	push (@string, '-');
							} else { 								push (@string, $allele1[$i]);
						}} else {
							my $hash;
							$$hash{$allele1[$i]}++; $$hash{$allele2[$i]}++;
							push (@string, &get_amb_code($hash));
					}}
					my $tmp = join ("", @string);
					$hash2return{$keys} = $tmp;
					#print $array[1]."\n"; print $array2[1]."\n"; print $tmp."\n\n";
				} elsif ($alleles > 2) {
					DOMINO::printError("ERROR: $region_id: [ $keys ] contains more than 2 alleles\nDOMINO would not die here but this locus would be skipped....\n");
					$pm_MARKER_MSA_files->finish;
			}}
			#print "\n\n"; print Dumper \%hash2return;

			undef %hash2fill;
			unless ($dnaSP_flag) {
				my $valueReturned = &check_marker_pairwise(\%hash2return);
				if ($valueReturned != 1) { $pm_MARKER_MSA_files->finish; }
			}
			$string2print_ref = &check_marker_ALL(\%hash2return, "Ref"); # Check there is a minimun variation
			if ($string2print_ref eq 'NO' ) { $pm_MARKER_MSA_files->finish;}
		} else {
			#### OTHER MSAs
			my @path_file = split("\.fasta", $array_files_fasta_msa[$i]);
			$region_id = $path_file[0];
			$hash_ref_msa = DOMINO::readFASTA_hash($file_path);				
			my $taxa4marker = scalar keys %{ $hash_ref_msa };
			if ($taxa4marker < $minimum_number_taxa_covered) { $pm_MARKER_MSA_files->finish; }
			unless ($dnaSP_flag) {
				my $valueReturned = DOMINO::check_marker_pairwise($hash_ref_msa, $$hash_parameters{'marker'}{'MCT'}[0], $$hash_parameters{'marker'}{'variable_positions_user_min'}[0], $$hash_parameters{'marker'}{'variable_positions_user_max'}[0], $$hash_parameters{'marker'}{'variable_divergence'}[0]);
				if ($valueReturned != 1) { $pm_MARKER_MSA_files->finish; }
			}
			$string2print_ref = &check_marker_ALL($file_path); # Check there is a minimun variation
			if ($string2print_ref eq 'NO' ) { $pm_MARKER_MSA_files->finish;}
		}
		
		#print Dumper $string2print_ref;
		## taxa_included	var_sites	length	profile	effective_length
		##	0					1		2		3		4			
		
		my ($taxa, $var_sites, $length_string, $string_profile, $effective_length) = @{ $string2print_ref };
		my @array_taxa_split = split(",", $taxa);
		## Control steps
		unless ($variable_divergence) { 
			if ($var_sites > $$hash_parameters{'marker'}{'variable_positions_user_max'}[0]) {$pm_MARKER_MSA_files->finish;} 
			if ($var_sites < $$hash_parameters{'marker'}{'variable_positions_user_min'}[0]) {$pm_MARKER_MSA_files->finish;} 
		}
		
		if ($select_markers) {
			#### SELECTION MODE
			
			## Check marker
			my $variation_perc = ($var_sites/$effective_length)*100;
			my $h = sprintf ("%.3f", $variation_perc);
			my $string = $region_id."\t".$taxa."\t".$var_sites."\t".$effective_length."\t".$h;		
			push (@{ $domino_files_msa{$region_id}{'string'} }, $string);
										
			## Print profile
			my $profile_dir_file = $profile_dir."/".$region_id."_profile.txt";
			open (PRF, ">$profile_dir_file"); print PRF ">".$region_id."\n".$string_profile."\n"; close (PRF);
			#push (@{ $domino_files_msa{$region_id}{'profile'} }, $profile_dir_file);
		
			#print $file_path."\n"; print Dumper \%domino_files_msa;
			# Print format: same as input and *mmfas 				
			my $msa_fasta = $msa_dir."/".$region_id.".fasta";
			open (OUT_MSA, ">$msa_fasta");
			foreach my $keys (sort keys %{ $hash_ref_msa }) { print OUT_MSA ">".$keys."\n".$$hash_ref_msa{$keys}."\n"; }
			close (OUT_MSA);
			#push (@{ $domino_files_msa{$region_id}{'markers'} }, $msa_fasta);

			my $dump_folder_files = $dir_Dump_file."/dump_markers_".$region_id.".txt";
			# Dump into file # print Dumper \%domino_files_msa;
			DOMINO::printDump(\%domino_files_msa, $dump_folder_files);	
		} elsif ($identify_markers) { ## Identify markers in MSA alignments
			#### DISCOVERY MODE

			my $profile_dir_file = $profile_dir."/".$region_id.".txt";
			open (OUT, ">$profile_dir_file"); print OUT ">".$region_id."\n".$string_profile."\n"; close(OUT);
			push (@{ $domino_files_msa{$region_id}{'profile'} }, $profile_dir_file);
	

			my $mergeProfile = $profile_dir."/".$region_id."_merged_ARRAY.txt";
			my $string = $window_size_VARS_range;$string =~ s/\:\:/-/; my $string2 = $window_size_CONS_range; $string2 =~ s/\:\:/-/;	
			my $mergeCoord;

			if ($variable_divergence) { $mergeCoord = $profile_dir."/".$region_id."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
			} else { 					$mergeCoord = $profile_dir."/".$region_id."-VPmin_".$domino_params{'marker'}{'variable_positions_user_min'}[0]."-VPmax_".$domino_params{'marker'}{'variable_positions_user_max'}[0]."-CL_".$string2."-CD_".$domino_params{'marker'}{'window_var_CONS'}[0]."-VL_".$string.".tab";
			}
			push (@{ $domino_files_msa{$region_id}{'mergeProfile'} }, $mergeProfile);
			push (@{ $domino_files_msa{$region_id}{'mergeCoord'} }, $mergeCoord);
			open (OUT_COORD, ">$mergeCoord");

			## Identify markers in MSA alignments
			my $infoReturned = DOMINO::sliding_window_conserve_variable(\$region_id, \$string_profile); 
			#print $$fileReturned."\n";
			if (!$infoReturned) { $pm_MARKER_MSA_files->finish;
			} else { 
				my @array = @$infoReturned;
				for (my $j=0; $j < scalar @array; $j++) {
					print OUT_COORD $array[$j]."\n";
			}}
			close (OUT_COORD);
			
			my $file_markers_collapse = DOMINO::check_overlapping_markers($mergeCoord, \$profile_dir_file);
			push (@{ $domino_files_msa{$region_id}{'markers_Merge'} }, $file_markers_collapse);
			
			# Retrieve fasta sequences...
			my $output_file = $msa_dir_tmp."/".$region_id."_markers_retrieved.txt";				
			my $array_Ref = DOMINO::check_DOMINO_marker($output_file, $msa_dir_tmp, $file_markers_collapse, $file_path);
			unless (scalar @$array_Ref == 0) { 
				push (@{ $domino_files_msa{$region_id}{'markers_files'} }, @$array_Ref);
				my $dump_folder_files = $dir_Dump_file."/dump_markers_".$region_id.".txt";
				push (@{ $domino_files_msa{$region_id}{'markers'} }, $output_file);
				# Dump into file # print Dumper \%domino_files_msa;
				DOMINO::printDump(\%domino_files_msa, $dump_folder_files);	
}}} $pm_MARKER_MSA_files->finish();}
$pm_MARKER_MSA_files->wait_all_children;
print "\n\n";
print "**********************************************\n";
print "**** All checking processes have finished ****\n";
print "**********************************************\n\n";	
print "\n"; &time_log(); print "\n"; chdir $marker_dirname;

#################################################################################
##	Once the coordinates are found, print different files with the information ##
#################################################################################	
### open Output and Error file
my $output_file_coord = "DM_markers-summary.txt"; open (OUT_coord,">$output_file_coord")or die "Cannot write the file $output_file_coord [DM_MarkerScan: Write OUT_coord]";
print OUT_coord "Region\t\tTaxa_included\tVariable_Positions\tEffective_length\tVariation(%)\n";
my $out_file = "DM_markers";
if ($pyRAD_file) { $out_file .= ".loci"; } elsif ($stacks_file) {$out_file .= ".fa";} else {$out_file .= ".fasta";}
open (OUT, ">$out_file");
print "+ Printing selected markers in $output_file_coord and $out_file...\n";

if ($identify_markers) {
	my %hashRetrieve;
	my $array_files = DOMINO::readDir($dir_Dump_file);
	my @dump_files = @{ $array_files };
	for (my $j=0; $j < scalar @dump_files; $j++) {
		if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
		open (DUMP_IN, "$dir_Dump_file/$dump_files[$j]");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			#&debugger_print($line);
			push (@{ $hashRetrieve{$array[0]}{$array[1]}}, $array[2]);
	} close (DUMP_IN); }
	foreach my $regions (sort keys %hashRetrieve) {
		if ($hashRetrieve{$regions}{'markers'}) {
			my @array_coord = @{$hashRetrieve{$regions}{'markers'}};
			for (my $i=0; $i < scalar @array_coord; $i++) {
				open (FILE, $array_coord[$i]); while (<FILE>) { print OUT_coord $_; } close(FILE);					
		}}
		if ($hashRetrieve{$regions}{'markers_files'}) {
			my @array_markers = @{$hashRetrieve{$regions}{'markers_files'}};
			for (my $j=0; $j < scalar @array_markers; $j++) {
				File::Copy::move($array_markers[$j], $msa_dir);
}}}} else {		
	my %hashRetrieve;
	my $array_files_dump = DOMINO::readDir($dir_Dump_file);
	my @dump_files = @{ $array_files_dump };
	for (my $j=0; $j < scalar @dump_files; $j++) {
		if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
		open (DUMP_IN, "$dir_Dump_file/$dump_files[$j]"); while (<DUMP_IN>) { print OUT_coord $_; } close (DUMP_IN); 
	}
	my $array_files_msa = DOMINO::readDir($msa_dir);
	my @msa_files = @{ $array_files_msa };
	for (my $j=0; $j < scalar @msa_files; $j++) {			
		if ($msa_files[$j] eq '.' || $msa_files[$j] eq '..' || $msa_files[$j] eq '.DS_Store') { next;}			
		if ($pyRAD_file) { ## print Loci
			my ($msa_hash_fasta_ref) = DOMINO::readFASTA_hash("$msa_dir/$msa_files[$j]"); ## Obtain reference of a hash
			foreach my $keys (sort keys %{$msa_hash_fasta_ref}) { print OUT ">".$keys."\t".$$msa_hash_fasta_ref{$keys}."\n"; } print OUT "//\n";			
		} elsif ($stacks_file) { ## print stacks file
			open (MSA, "$msa_dir/$msa_files[$j]"); while (<MSA>) { print OUT $_; } close (MSA);
		} else { open (MSA, "$msa_dir/$msa_files[$j]"); while (<MSA>) { print OUT $_; } close (MSA); print OUT "//\n";			
}}} close(OUT_coord); close (OUT);

## USE THE SUBROUTINE print_Excel and control if radseq_like_data
print "+ Done...\n+ Retrieving informative locus has been done...\n+ Generating an Excel file for DOMINO markers identified...\n";
my $excelBook = DOMINO::print_Excel(\$output_file_coord, \$marker_dirname);

unless ($avoidDelete_tmp_files) { 
	###########################
	## Delete temporary file ##
	###########################
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Deleting temporary files ", "#"); DOMINO::printHeader("", "#");
	remove_tree($msa_dir_tmp);
	remove_tree($dir_Dump_file);
}	
File::Copy::move($param_Detail_file_markers, $marker_dirname."/");
if (-z $mapping_markers_errors_details) { remove_tree($mapping_markers_errors_details); 
} else { File::Copy::move($mapping_markers_errors_details, $marker_dirname."/"); }	## Finish and exit
DOMINO::finish_time_stamp($start_time); print "\n\n Job done succesfully, exiting the script\n\n\n"; exit();

###########################
####### SUBROUTINES #######
###########################

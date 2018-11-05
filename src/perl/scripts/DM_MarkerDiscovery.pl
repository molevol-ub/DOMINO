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
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
#################################################################

## Arguments
my $path = $ARGV[0];
my $step_time = $ARGV[1];
my $ref_taxa = $ARGV[2];
my $marker_dir = $ARGV[3];
#################################################################

#################################################################
## Get general parameters and files
my $hash_parameters = DOMINO::get_parameters($path."/");
my %domino_marker_files = %{ DOMINO::get_DOMINO_files($path."/") };
	#print Dumper $hash_parameters;
	#print Dumper \%domino_marker_files;

my $align_dirname = $$hash_parameters{'mapping'}{'folder'}[0];
my $num_proc_user = $$hash_parameters{'mapping'}{'cpu'}[0];
my $msa_dirname = $align_dirname."/MSA_files";

my %new_domino_files;
my %new_dowmino_param;
#################################################################

mkdir $marker_dir, 0755; chdir $marker_dir; #&debugger_print("Changing dir to $marker_dir");

### Get variables
my $totalContigs2use4markers = $$hash_parameters{"marker"}{"totalContigs2use4markers"}[0];
my $variable_divergence = $$hash_parameters{"marker"}{"variable_divergence"}[0];
my $minimum_number_taxa_covered = $$hash_parameters{"marker"}{"MCT"}[0];
my $option = $$hash_parameters{"marker"}{"option"}[0];
my $subset_offset_user = $$hash_parameters{"marker"}{"subset_offset_user"}[0];

#######################################
###		 Copy necessary files		### 
#######################################
print "+ Retrieve necessary files\n";
if ($option eq "msa_alignment") { 
	push ( @{ $new_domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/merged.profile_ARRAY.txt");
} else {
	my @uniq_sort_taxa = do { my %seen; grep { !$seen{$_}++ } @{ $domino_marker_files{'taxa'}{'user_Taxa'} } };
	my $name = join("_", @uniq_sort_taxa);
	push ( @{ $new_domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/$name.profile_ARRAY.txt");
}

##################################
## Get the the reference fasta  ##
##################################
print "+ Checking the file specified as reference fasta... $domino_marker_files{$ref_taxa}{'contigs'}[0]...OK\n";
print "+ Generating symbolic links for necessary files\n"; 
my $ref_Fasta = $ref_taxa.".fasta";
system("ln -s $domino_marker_files{$ref_taxa}{'contigs'}[0] $ref_Fasta");

############################################################
###	Generate an array information file for the reference ### 
############################################################
print "+ Reading the reference fasta file...\n";
my $fasta_seqs = DOMINO::retrieve_info(\@{ $domino_marker_files{$ref_taxa}{"hash_reference_file"} }, 1);

## Print each contig into a file
print "+ Splitting file into multiple files to speed the computation...\n";
print "+ Using parallel threads ($num_proc_user CPUs)...\n";			
my $reference_dir = $marker_dir."/REF_DIR"; mkdir $reference_dir, 0755;
push (@{ $new_domino_files{$ref_taxa}{'REF_DIR'}}, $reference_dir);
my $size_file = DOMINO::get_size($ref_Fasta);
my $parts2split = int($size_file/$num_proc_user);
my $fasta_files_split = DOMINO::fasta_file_splitter($ref_Fasta, $parts2split, "fasta", $reference_dir);
	#&debugger_print("Total Size: $size_file\nCharacters to split: $parts2split"); #&debugger_print("Ref", $fasta_files_split);		

## Get size for each contig
my $pm_SPLIT_FASTA =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
my $total_fasta_files_split = scalar @{ $fasta_files_split };
my $count_fasta_files_split=0;
for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
	$count_fasta_files_split++;
	print "\t- Splitting file: [$count_fasta_files_split/$total_fasta_files_split]\n";
	my $pid = $pm_SPLIT_FASTA->start($i) and next;
	open(FILE, $$fasta_files_split[$i]) || die "Could not open the $$fasta_files_split[$i]... [DM_MarkerScan: fasta_files_split size]\n";
	my $size_file = $$fasta_files_split[$i]."_size";
	open (SIZE, ">$size_file");
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
		my ($titleline, $sequence) = split(/\n/,$_,2);
		next unless ($sequence && $titleline);
		chomp $sequence; $sequence =~ s/\n//g; $titleline =~ s/\s/\t/g;
		print SIZE length($sequence)."\t".$titleline."\n";
	} close(FILE); $/ = "\n"; close (SIZE);
	$pm_SPLIT_FASTA->finish($i); # pass an exit code to finish
}
$pm_SPLIT_FASTA->wait_all_children;

print "\n+ Concatenate contig size for each splitted file...\n";
my $ref_Fasta_size = $ref_taxa.".size";
my $tmp_Fasta_size = $ref_taxa.".size_tmp";
for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
	my $size_file = $$fasta_files_split[$i]."_size";
	system ("cat $size_file >> $tmp_Fasta_size"); system ("rm $size_file");
}
print "\n+ Sorting contig size...\n"; system ("sort -nr $tmp_Fasta_size >> $ref_Fasta_size"); system ("rm $tmp_Fasta_size");

## get ordered size of contigs
my $total_contigs=0; my @size_contigs;
open (FILE, $marker_dir."/".$ref_Fasta_size); while (<FILE>) { $total_contigs++; } close (FILE); ## total contigs provided	

## check all or some
if ($totalContigs2use4markers == -1) { ## use all contigs
	$totalContigs2use4markers = $total_contigs;		 
} elsif ( $totalContigs2use4markers > $total_contigs) {
	$totalContigs2use4markers = $total_contigs;		 
}
#&debugger_print("totalContigs2use4markers\ttotal_contigs?\n".$totalContigs2use4markers."\t".$total_contigs);

if ($total_contigs > $totalContigs2use4markers) { 
	print "\n\nATTENTION: There are too many contigs. DOMINO would only check for markers in the largest 50.000 ones\n";
	print "ATTENTION: Provide option -all for using the whole set\n\n";
	
	print "+ Only retrieve up to $totalContigs2use4markers...\n";
	my $counter_stop=0;
	open (FILE, $marker_dir."/".$ref_Fasta_size);
	my $Fasta_ids2retrieve = $marker_dir."/".$ref_taxa.".ids2retrieve";
	open (IDS, ">$Fasta_ids2retrieve");
	while (<FILE>) {
		chomp; my @array = split("\t", $_);
		push (@size_contigs, $array[1]);
		$counter_stop++;
		print IDS $array[1]."\n";
		if ($counter_stop == $totalContigs2use4markers) { last;}
	}
	close (FILE); close (IDS);
	DOMINO::mothur_retrieve_FASTA_seqs($ref_Fasta, $marker_dir, $Fasta_ids2retrieve, $mothur_path); ## file generated:
	my $ref_Fasta_ids2retrieve = $marker_dir."/".$ref_taxa.".pick.fasta";
	
	## Get size for each contig
	my $size_file_retrieved = DOMINO::get_size($ref_Fasta_ids2retrieve);
	my $parts2split_retrieved = int($size_file_retrieved/$num_proc_user);
	my $fasta_files_split_retrieved = DOMINO::fasta_file_splitter($ref_Fasta_ids2retrieve, $parts2split_retrieved, "fasta", $marker_dir);

	my $pm_SPLIT_FASTA =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	my $total_fasta_files_split_retrieved = scalar @{ $fasta_files_split_retrieved };
	my $count_fasta_files_split_retrieved=0;
	for (my $i=0; $i < scalar @{ $fasta_files_split_retrieved }; $i++) {
		$count_fasta_files_split_retrieved++;
		print "\t- Splitting file: [$count_fasta_files_split_retrieved/$total_fasta_files_split_retrieved]\n";
		my $pid = $pm_SPLIT_FASTA->start($i) and next;
		open(FILE, $$fasta_files_split_retrieved[$i]) || die "Could not open the $$fasta_files_split_retrieved[$i]... [DM_MarkerScan: fasta_files_split_retrieved]\n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			my $file = $reference_dir."/".$titleline.".fasta";
			open (OUT, ">$file"); print OUT $titleline."\n".$sequence."\n"; close(OUT);
			chomp $sequence; $sequence =~ s/\n//g; $titleline =~ s/\s/\t/g;
		} $/ = "\n"; close (SIZE);
		remove_tree($$fasta_files_split_retrieved[$i]);
		$pm_SPLIT_FASTA->finish($i); # pass an exit code to finish
	}
	$pm_SPLIT_FASTA->wait_all_children;
	
} else { ## use all contigs

	$totalContigs2use4markers = $total_contigs;
	
	## push into array ids
	open (FILE, $marker_dir."/".$ref_Fasta_size);
	while (<FILE>) { chomp; my @array = split("\t", $_); push (@size_contigs, $array[1]); } close (FILE);
	
	## If all contigs to be used, generate individual fasta for each
	my $pm_SPLIT_FASTA =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	my $total_fasta_files_split = scalar @{ $fasta_files_split };
	my $count_fasta_files_split=0;
	for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
		$count_fasta_files_split++;
		print "\t- Splitting file: [$count_fasta_files_split/$total_fasta_files_split]\n";
		my $pid = $pm_SPLIT_FASTA->start($i) and next;
		open(FILE, $$fasta_files_split[$i]) || die "Could not open the $$fasta_files_split[$i]...[DM_MarkerScan: pm_SPLIT_FASTA]\n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence; $sequence =~ s/\n//g; $titleline =~ s/\s/\t/g;
			my $file = $reference_dir."/".$titleline.".fasta";
			open (OUT, ">$file"); print OUT $titleline."\n".$sequence."\n"; close(OUT);
		} close(FILE); $/ = "\n";
		$pm_SPLIT_FASTA->finish($i); # pass an exit code to finish
	} 
	$pm_SPLIT_FASTA->wait_all_children;	
}

## remove files	
my $new_count_fasta_files_split=0;
for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
	$new_count_fasta_files_split++;
	remove_tree($$fasta_files_split[$i]);
	print "\t- Remove file: [$new_count_fasta_files_split/$total_fasta_files_split]\r";
} print "\n";

## Number of contigs in each subset 
my $subset_offset = $subset_offset_user;

#########################################
## Merge PILEUP information arrays     ##
#########################################
print "\n"; DOMINO::printHeader(" Fetching information from all the PROFILEs generated ", "#");
my %pileup_files;
my (@clean_filtered_sam_files, @reference_bam_files, @pileup_Arrays);
foreach my $reads (sort keys %domino_marker_files) {
	unless ($domino_marker_files{$reads}{'taxa'}) {next; }
	if ($domino_marker_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa}) { push ( @reference_bam_files, @{ $domino_marker_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa}});}		
	if ($domino_marker_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa}) { push ( @clean_filtered_sam_files, @{ $domino_marker_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa}});}
	next if ($reads eq $ref_taxa);
}
my $PILEUP_merged_folder_abs_path = $marker_dir."/PROFILE_merge_species";
mkdir $PILEUP_merged_folder_abs_path, 0755; chdir $PILEUP_merged_folder_abs_path; 
#&debugger_print("Changing dir to $PILEUP_merged_folder_abs_path");
push (@{ $new_domino_files{$ref_taxa}{'PROFILE_FOLDER'}}, $PILEUP_merged_folder_abs_path);

print "+ Checking profiles of variation for each contig and merging information...\n";
print "+ Using a sliding window approach...\n"; print "+ Using parallel threads ($num_proc_user CPUs)...\n";			

## Check each markers using threads
my $dir_Dump_file = $PILEUP_merged_folder_abs_path."/DUMP_files"; mkdir $dir_Dump_file, 0755;
my $dir2print_markers = $PILEUP_merged_folder_abs_path."/MSA_fasta_tmp"; mkdir $dir2print_markers, 0755;
	
## Check for markers: USING THREADS, ONE FOR EACH BLOCK OF CONTIGS
my $pm_MARKER_PILEUP =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
my $SETS;
if ($totalContigs2use4markers % 2) { 	$SETS = int($totalContigs2use4markers/$subset_offset) + 1;			## Even number
} else { 								$SETS = int($totalContigs2use4markers/$subset_offset); }	## Odd number 	
print "+ Dataset would be splitted for speeding computation into $SETS subsets...\n";

my $counter=0; 
#&debugger_print("TOTAL contigs:\n".$totalContigs2use4markers);

## Generate subsets of given amount of contigs to avoid collapsing system with so many files
for (my $set=1; $set <= $SETS; $set++) {
	my @subset_array; my $tmp=1;
	for ($counter=$counter; $counter < $totalContigs2use4markers; $counter++) {
		#&debugger_print("Counter\tTotal\n".$counter."\t".$totalContigs2use4markers);
		push(@subset_array, $size_contigs[$counter]);
		#&debugger_print("SETS:$SETS\tset:$set\tCounter:".$counter."\ttmp: $tmp\tContig:  $size_contigs[$counter]");
		if ($tmp eq $subset_offset) { $counter++; last; }
		$tmp++;
	}
	#&debugger_print("Ref", \@subset_array);
	
	## SEND THREAD 
	my $pid = $pm_MARKER_PILEUP->start($set) and next;
	
	my (%pileup_files_threads, %contigs_pileup_fasta);
	foreach my $reads (sort keys %domino_marker_files) {
		next if ($reads eq $ref_taxa);
		next if ($reads eq 'taxa'); next if ($reads eq 'main'); next if ($reads eq 'original');
		next if ($reads eq 'assembly'); next if ($reads eq 'clean_data'); next if ($reads eq 'mapping');
		next if ($reads eq 'markers'); 

		for (my $j=0; $j < scalar @subset_array; $j++) {
			my $file = $domino_marker_files{$reads}{"PROFILE::Ref:".$ref_taxa}[0]."/".$subset_array[$j]."_ARRAY.txt";
		
			unless (-e -r -s $file) { next; } # skip
			my $tmp_hash_reference = DOMINO::readFASTA_hash($file);
			my %tmp_fasta = %{$tmp_hash_reference};
			foreach my $seqs (sort keys %tmp_fasta) {
				push (@{ $contigs_pileup_fasta{$subset_array[$j]} }, $tmp_fasta{$seqs});
	}}}		

	my $mergeProfile = $PILEUP_merged_folder_abs_path."/SET_$set"."_merged_ARRAY.txt";
	my $string = $$hash_parameters{'marker'}{"window_size_VARS_range"}[0];$string =~ s/\:\:/-/; 
	my $string2 = $$hash_parameters{'marker'}{"window_size_CONS_range"}[0]; $string2 =~ s/\:\:/-/;	
	my $mergeCoord;
	if ($variable_divergence) { $mergeCoord = $PILEUP_merged_folder_abs_path."/SET_$set"."-VD_".$variable_divergence."-CL_".$string2."-CD_".$$hash_parameters{'marker'}{'window_var_CONS'}[0]."-VL_".$string.".tab";
	} else { 					$mergeCoord = $PILEUP_merged_folder_abs_path."/SET_$set"."-VPmin_".$$hash_parameters{'marker'}{'variable_positions_user_min'}[0]."-VPmax_".$$hash_parameters{'marker'}{'variable_positions_user_max'}[0]."-CL_".$string2."-CD_".$$hash_parameters{'marker'}{'window_var_CONS'}[0]."-VL_".$string.".tab";
	}
	push (@{ $pileup_files_threads{"SET_$set"}{'mergeProfile'} }, $mergeProfile);
	push (@{ $pileup_files_threads{"SET_$set"}{'mergeCoord'} }, $mergeCoord);
	my $markers_shared_file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.tmp.txt";

	## NAME for output merged files
	my ($output_merged_file, $error_merged_file, $file);
	if ($variable_divergence) {	$file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared-VD_".$variable_divergence."-CL_".$string2."-CD_".$$hash_parameters{'marker'}{'window_var_CONS'}[0]."-VL_".$string.".tab";
	} else {  					$file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared-VPmin_".$$hash_parameters{'marker'}{'variable_positions_user_min'}[0]."-VPmax_".$$hash_parameters{'marker'}{'variable_positions_user_max'}[0]."-CL_".$string2."-CD_".$$hash_parameters{'marker'}{'window_var_CONS'}[0]."-VL_".$string.".tab";
	}		
	$output_merged_file = $file.".out";  $error_merged_file = $file.".err";
	push (@{ $pileup_files_threads{"SET_$set"}{'eachTaxaCoord'} }, $output_merged_file);

	## Merging variable and conserved information into a unique array, profile and generate coordinates
	my $missing_allowed_species = $minimum_number_taxa_covered;
	my $SLIDING_file = $file."_sliding_window_MERGED.txt";
	open (OUT_COORD, ">$mergeCoord"); open (SHARED, ">$markers_shared_file");
	foreach my $seqs (keys %contigs_pileup_fasta) {
		my $size = $$fasta_seqs{$seqs};
		my $tmp_string;			
		for (my $i = 0; $i < scalar $size; $i++) {
			my (@tmp, $pb);
			for (my $j = 0; $j < scalar @{ $contigs_pileup_fasta{$seqs} }; $j++){
				$pb = substr($contigs_pileup_fasta{$seqs}[$j], $i, 1);
				push (@tmp, $pb);					
			}
			my $flag = 0; my $missing_count = 0; my $missing_flag = 'N';
			my $scalar_keys = scalar @tmp;
			for (my $keys = 0; $keys < scalar @tmp; $keys++) {
				if ($tmp[$keys] eq 1) { # One specie has a variable position
					$flag = 1; 
				} elsif ($tmp[$keys] eq 'N') { ## this position is missing for that specie
					$missing_count++;
				}}	
			my $tmp = $$hash_parameters{'mapping'}{'number_sp'}[0] - $missing_count;
			if ($tmp < $missing_allowed_species) {
				$tmp_string .= $missing_flag; 
			} else { $tmp_string .= $flag; }
			undef @tmp;
		}
		
		my $var_sites = $tmp_string =~ tr/1/1/; ## count variable sites
		my $cons_sites = $tmp_string =~ tr/0/0/; ## count conserved sites
		if ($var_sites != 0 && $cons_sites != 0) { 
			open (FH, ">>$mergeProfile"); print FH ">$seqs\n$tmp_string\n"; close(FH);
		} else { next; } 
		my $infoReturned = DOMINO::sliding_window_conserve_variable(\$seqs, \$tmp_string);
		if ($infoReturned) { 
			open (TMP_COORD, ">$SLIDING_file");
			for (my $j=0; $j < scalar @$infoReturned; $j++) {
				print OUT_COORD $$infoReturned[$j]."\n";
				print TMP_COORD $$infoReturned[$j]."\n";
				print SHARED $$infoReturned[$j]."\t".$ref_taxa."\n";
			}
			close (TMP_COORD);
		} else { delete $contigs_pileup_fasta{$seqs}; next; ## if empty next
		}
		
		delete $contigs_pileup_fasta{$seqs}; undef $infoReturned;
		
		######################################################################
		## Check the coordinates foreach taxa against the merge statistics  ##
		######################################################################
		foreach my $taxa (sort keys %domino_marker_files) {
			if ($taxa eq 'taxa') {next;}
			unless ($domino_marker_files{$taxa}{'taxa'}) { next; }
			if ($taxa eq $ref_taxa) {next;}	
			#print "Ref_taxa: $ref_taxa Taxa: $taxa\n";
			my $pileup_each_taxa = $domino_marker_files{$taxa}{"PROFILE::Ref:".$ref_taxa}[0]."/".$seqs."_ARRAY.txt";
			## For each taxa confirm profile
			if (-f $pileup_each_taxa) { &get_coordinates_each_taxa(\$pileup_each_taxa, $SLIDING_file, $taxa, \$output_merged_file, \$error_merged_file);
	}}}
	close (OUT_COORD); close (SHARED);

	undef %contigs_pileup_fasta; # no longer necessary
	unless (-e -r -s $mergeCoord) { #if empty file
		undef %pileup_files_threads; $pm_MARKER_PILEUP->finish(); 
	} 
	unless (-e -r -s $mergeProfile) { #if empty file
		undef %pileup_files_threads; $pm_MARKER_PILEUP->finish(); 
	} 

	##########################################
	## Get Coordinates of Molecular Markers ##
	##########################################		
	my $shared_markers_others = $pileup_files_threads{"SET_$set"}{'eachTaxaCoord'}[0];
	my $markers_shared_file_sort = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.tmp_sort.txt";
	system("cat $shared_markers_others >> $markers_shared_file ; sort $markers_shared_file > $markers_shared_file_sort");
	my %coord_contig;
	open (MERGE_COORD,"<$markers_shared_file_sort") or die "Cannot open file $markers_shared_file_sort [DM_MarkerScan: Check Merge coordinates]";
	while(<MERGE_COORD>){
		chomp;
		my $line = $_;
		next if ($line =~ m/^Contig\tCons1_ini:.*/o);
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;
		my @array = split("\t",$line);
		push (@{ $coord_contig{$array[0]}{$array[1].";".$array[2].";".$array[3]}}, $array[4]);
	}
	close (MERGE_COORD);	
	if (!%coord_contig) { 
		undef %pileup_files_threads; 
		$pm_MARKER_PILEUP->finish(); ## finish thread 
	}
	
	## Print in tmp file for sorting and obtaining unique
	chdir $PILEUP_merged_folder_abs_path;
	my $markers_shared = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.txt";
	open(TMP, ">$markers_shared");
	for my $scaffold (keys %coord_contig) {    		
		foreach my $marker (keys %{ $coord_contig{$scaffold} }) {
			my @taxa = @{ $coord_contig{$scaffold}{$marker} };
			## Check we find markers for all the taxa desired
			if (scalar @taxa < $minimum_number_taxa_covered) { next; }   		
			## Write DOMINO Markers Coordinates in tmp txt file
			my @string = split(";", $marker);
			my $string = join("\t", @string);
			my @sort_taxa = sort(@taxa);
			my $string_taxa = join(",", @sort_taxa); #taxa providing variation
			print TMP "$scaffold\t$string\t$string_taxa\n"; 
		}
	} close(TMP);
	undef %coord_contig;
	unless (-e -r -s $markers_shared) { undef %pileup_files_threads; undef %contigs_pileup_fasta; $pm_MARKER_PILEUP->finish();  }

	## Collapse markers
	my $file_markers_collapse = DOMINO::check_overlapping_markers($markers_shared, $pileup_files_threads{"SET_$set"}{'mergeProfile'}[0]);

	# Retrieve fasta sequences...
	my $output_file = $PILEUP_merged_folder_abs_path."/SET_".$set."_markers_retrieved.txt";
	my $markers_print_ref = DOMINO::check_DOMINO_marker($output_file, $dir2print_markers, $file_markers_collapse, $ref_taxa);

	unless (scalar @$markers_print_ref == 0) { 
		push (@{ $pileup_files_threads{"SET_$set"}{'markers'} }, $output_file);
		push (@{ $pileup_files_threads{"SET_$set"}{'markers_files'} }, @$markers_print_ref);
		my $dump_folder_files = $dir_Dump_file."/dump_markers_SET_".$set.".txt";
		DOMINO::printDump(\%pileup_files_threads, $dump_folder_files);	
		undef %pileup_files_threads; undef %contigs_pileup_fasta;
	}
	$pm_MARKER_PILEUP->finish(); # finish for each contig
} 
$pm_MARKER_PILEUP->wait_all_children; #each marker
print "\n\n";
print "******************************************************\n";
print "**** All parallel parsing processes have finished ****\n";
print "******************************************************\n\n";
&time_log(); print "\n";

## Retrieve info of all markers identified...
my $dump_files = DOMINO::readDir($dir_Dump_file);
for (my $i=0; $i < scalar @$dump_files; $i++) {
	if ($$dump_files[$i] eq '.' || $$dump_files[$i] eq '..' || $$dump_files[$i] eq '.DS_Store') { next;}
	$$dump_files[$i] = $dir_Dump_file."/".$$dump_files[$i];
}
my %markers2retrieve;
for (my $j=0; $j < scalar @{ $dump_files }; $j++) {
	if ($$dump_files[$j] eq '.' || $$dump_files[$j] eq '..' || $$dump_files[$j] eq '.DS_Store') { next;}
	open (DUMP_IN, "$$dump_files[$j]");
	while (<DUMP_IN>) {
		my $line = $_; chomp $line; my @array = split("\t", $line);
		push (@{ $markers2retrieve{$array[0]}{$array[1]}}, $array[2]);
} close (DUMP_IN); }
# print Dumper \%markers2retrieve; #

#################################################################
## Get the information ready for the user to visualize contigs ##
#################################################################
print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Getting the information ready to present ", "#"); DOMINO::printHeader("", "#"); 
my $markers_msa_folder = $marker_dir."/MSA_markers";  mkdir $markers_msa_folder, 0755; chdir $marker_dir;

## Open output file for markers coordinates
my $output_file;	
if ($option eq "genome") { $output_file = $marker_dir."/DM_markers-coordinates.txt";		
} else { $output_file = $marker_dir."/DM_markers-coordinates_ref_".$ref_taxa.".txt";	}
open (OUT_markers, ">$output_file"); 

## Open output file for contigs
my $output_file_putative_contigs = $marker_dir."/DM_contigs.fasta";
open (OUT, ">$output_file_putative_contigs"); 

## Open output file for markers sequences
my $file_coordinates = $marker_dir."/DM_sequence-markers.fasta"; 
open (OUT_coord, ">$file_coordinates");

push(@{ $new_domino_files{$ref_taxa}{'MSA_markers'}}, $markers_msa_folder);
push(@{ $new_domino_files{$ref_taxa}{'markers'}}, $file_coordinates);
push(@{ $new_domino_files{$ref_taxa}{'CONTIGS'}}, $output_file_putative_contigs);
push(@{ $new_domino_files{$ref_taxa}{'coordinates'}}, $output_file);

print "+ Copying reference fasta contigs...\n+ Printing reference sequences in $output_file_putative_contigs...\n";
foreach my $subset (sort keys %markers2retrieve) {		
	## Move msa markers
	my @array_markers = @{ $markers2retrieve{$subset}{'markers_files'} };
	for (my $i=0; $i < scalar @array_markers; $i++) {
		File::Copy::move($array_markers[$i], $markers_msa_folder);	
	}
	
	my @ALL_contigs;
	# Get all contigs involved
	my $marker_contigs = $markers2retrieve{$subset}{'markers'}[0];
	open (IN_marker, "<$marker_contigs"); 
	while (<IN_marker>) { 
		print OUT_markers $_; ## print coordinates
		# Print sequence coordinates
		my $line = $_; chomp $line;
		my @array = split("\t", $line); 
		my @contig_name = split("_#", $array[0]);
		push (@ALL_contigs, $contig_name[0]);
	} close (IN_marker);
	
	# sort uniq
	my @uniq_contigs = do { my %seen; grep { !$seen{$_}++ } @ALL_contigs };
		
	## Printing contigs
	my %hash_contigs;
	for (my $c = 0; $c < scalar @uniq_contigs; $c++) {
		my $in_file = $reference_dir."/".$uniq_contigs[$c].".fasta";
		if (-e -r -s $in_file) { 
			my $hash_contigs = DOMINO::readFASTA_hash($in_file);
			$hash_contigs{$uniq_contigs[$c]} = $$hash_contigs{$uniq_contigs[$c]};
			print OUT ">".$uniq_contigs[$c]."\n".$$hash_contigs{$uniq_contigs[$c]}."\n";			
	}}
	## Get markers	
	open (IN_marker, "<$marker_contigs"); 
	while (<IN_marker>) { 
		# Print sequence coordinates
		my $line = $_; chomp $line;
		my @array = split("\t", $line); my @array1 = split(":", $array[1]); my @array3 = split(":", $array[3]);
		my @contig_name = split("_#", $array[0]);
		my ($seq_id, $seq, $array) = &fetching_range_seqs($contig_name[0], $array1[0], $array3[1], $hash_contigs{$contig_name[0]});
		print OUT_coord $seq_id."\n".uc($seq)."\n";
	} close (IN_marker);

} close(OUT); close(OUT_markers); close (OUT_coord);
print "+ Marker development for $ref_taxa is finished here...\n\n"; &time_log(); print "\n";

### dump info
my @array = ( $align_dirname."/DOMINO_dump_information.txt" );
%new_domino_files = %{ DOMINO::retrieve_info(\@array, \%new_domino_files)}; 
if (-r -e -s $array[0]) { remove_tree($array[0]); 
DOMINO::printDump(\%new_domino_files, $array[0]);} 

DOMINO::print_success_Step("marker_discovery");


###########################
####### SUBROUTINES #######
###########################

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}

sub get_coordinates_each_taxa {
	
	##########################################################################################
	##	 																					##
	##  This function checks the variation of each taxa against the merge PILEUP			##
	##  just to make sure not a unique specie is inflating the variation 					##
	## 		        																		##
	##	Cristina Frias 05/06/2014 	cristinafriaslopez@ub.edu								##
	## 		        																		##
	##########################################################################################
	
	my $file_MID_array = $_[0];
	my $merge_file_coord = $_[1];
	my $sp_id = $_[2];
	my $output_file = $_[3];
	my $error_file = $_[4];

	### declare vars script
	my ($contig_name, $seq_MID);
	my $coverage_var_pct = 0.75;

	### open FILE MID
	#print "\tChecking $sp_id coordinates...\n";	
	my $hash_ref_coordinates = DOMINO::readFASTA_hash($$file_MID_array);
	my %hash_seq = %{$hash_ref_coordinates};

	##### open output
	open (OUT, ">>$$output_file") or DOMINO::printError("Could not open output file $$output_file"); 
	open (ERR, ">>$$error_file") or DOMINO::printError("Could not open error file $$error_file");

	### open coord and check coord of file 1
	open (COORD,"<$merge_file_coord") or DOMINO::printError("Could not open merge coordenates $$merge_file_coord");
	while(<COORD>){
		chomp;
		my $line = $_;
		next if ($line =~ m/^Contig\tCons1_ini.*/o);
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;

		###################################################################
		###		 Get the coordinates of the different regions 			### 
		###################################################################

		my @regions= split(/\t/,$line);
		my $contig_MID = $regions[0];
		my $cons1 = $regions[1];
		my $var = $regions[2];
		my $cons2 = $regions[3];

		### split data of regions
		my ($cons1_START, $cons1_END) = split(/:/,$cons1);
		my ($VAR_START, $VAR_END) = split(/:/,$var);
		my ($cons2_START, $cons2_END) = split(/:/,$cons2);
		
		my ($coord_P1, $coord_P2) = split(/:/, $cons1);
		my ($coord_P3, $coord_P4) = split(/:/, $var);
		my ($coord_P5, $coord_P6) = split(/:/, $cons2);

		## Check if exists
		if (exists ($hash_seq{$contig_MID}) ){ $seq_MID = $hash_seq{$contig_MID};
		} else { print ERR "$contig_MID is not mapping for file $file_MID_array\n"; next; }

		###################################################################
		###			 Get the sequence of the different regions 			### 
		###################################################################
	
		### CONS1 ###
		my $size_P1_P2 = $coord_P2 - $coord_P1;
		my $string_P1_P2 = substr ($seq_MID, $coord_P1, $size_P1_P2);
		my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//;

		### VAR ###
		my $size_P3_P4 = $coord_P4 - $coord_P3;
		my $string_P3_P4 = substr ($seq_MID, $coord_P3, $size_P3_P4);
		my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//;

		### CONS2 ###
		my $size_P5_P6 = $coord_P6 - $coord_P5;
		my $string_P5_P6 = substr ($seq_MID, $coord_P5, $size_P1_P2);
		my $count_string_P5_P6 = $string_P5_P6 =~ tr/1//;

		###################################################################
		### 		Check the composition of the regions				###
		###################################################################
		
		## More than variations allowed in conserved
		if ($count_string_P1_P2 > $$hash_parameters{"marker"}{"window_var_CONS"}[0]) { next; } 
		if ($count_string_P5_P6 > $$hash_parameters{"marker"}{"window_var_CONS"}[0]) { next; }

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
		my $percent_total_length = $total_length * $$hash_parameters{"marker"}{"missing_allowed"}[0];  ## Default 0.05

		#if ($debugger) {
		#	print "\n\n***********************************************\n";
		#	print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
		#	print "Contig: $contig_MID\n\nPositions:\n";
		#	print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
		#	print "\nSubsets\n";
		#	print "Conserved (P1-P2): ".length($string_P1_P2)."\n"; print "$string_P1_P2\n";
		#	print "Variable (P3-P4): ".length($string_P3_P4)."\n";  print $string_P3_P4."\n";
		#	print "Conserved (P5-P6): ".length($string_P5_P6)."\n"; print $string_P5_P6."\n";
		#	print "\nVariations:\n";
		#	print "Count_string_P1_P2: $count_string_P1_P2\n";
		#	print "Count_string_P3_P4: $count_string_P3_P4\n";
		#	print "Count_string_P5_P6: $count_string_P5_P6\n";
		#	print "\nMarker:\n";
		#	print "Total Length:$total_length\tVariable Positions:$count_string_P3_P4\tExpected:$expected_var_sites\n";
		#	print "Missing allowed: $percent_total_length %\tMissing:$missing_count_percent %\n";
		#	print "***********************************************\n";
		#}

		if ($flag_error > 0) {	
			print ERR "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\t$sp_id\n"; next; 
		} else {
			### Print coordinates if the meet requirements
			if ($missing_count_percent < $percent_total_length) {
				print OUT "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\t$sp_id\n";
	}}} close(COORD); close(OUT); close(ERR);
}

sub fetching_range_seqs {
	
	my $contig_name = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $sequence = $_[3];
	
	my $length = $end - $start;	
	my $sub_seq = substr($sequence, $start, $length);
	my $seq_id = ">".$contig_name."_coord_".$start."_".$end;
	return($seq_id, $sub_seq);	
}
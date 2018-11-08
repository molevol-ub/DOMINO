#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
	require List::Uniq; use List::Uniq qw(uniq);
}

my $path = $ARGV[0];
my $step_time = $ARGV[1];

# Get memory
if(!defined($path)) { print "ERROR: No input files are provided: [$0]\n"; exit; } ## give error message for DOMINO debug
my $domino_version ="DOMINO v1.1 ## Revised 24-10-2018";

my $hash_parameters = DOMINO::get_parameters($path,"assembly");
my $domino_mira_files_Ref = DOMINO::get_DOMINO_files($path, "assembly");
my %domino_mira_files = %{$domino_mira_files_Ref};

my $assembly_directory = $$hash_parameters{"assembly"}{"dir"}[0];
my $dirname = $$hash_parameters{"assembly"}{"folder"}[0];
my $noOfProcesses=$$hash_parameters{"assembly"}{"CPU"}[0];
my $error_log = $$hash_parameters{'assembly'}{'errorLog'}[0];

my $scripts_path = $FindBin::Bin."/../";
my $MIRA_exec = $scripts_path."mira_v4.0/bin/mira";
my $CAP3_exec = $scripts_path."cap3/bin/cap3";
my $BLAST = $scripts_path."NCBI_BLAST/";
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";

foreach my $keys (sort keys %domino_mira_files) {
	next if $keys eq "original"; next if $keys eq "main";
	my $MIRA_manifest_file;
	my @species_files = @{ $domino_mira_files{$keys}{'original'} };
	
	if (scalar @species_files == 2) {
		$MIRA_manifest_file = &Generate_manifest_file($keys, 'Illumina', $species_files[0], $species_files[1]);	
	} elsif ($$hash_parameters{"assembly"}{"file_type"}[0] == 5) {		
		$MIRA_manifest_file = &Generate_manifest_file($keys, 'Illumina', $species_files[0]);	
	} elsif ($$hash_parameters{"assembly"}{"file_type"}[0] == 3) {
		$MIRA_manifest_file = &Generate_manifest_file($keys, '454', $species_files[0]);	
	}
	
	my $abs_path_manifest_file = $assembly_directory."/".$MIRA_manifest_file;
	push (@{ $domino_mira_files{$keys}{'manifest'}}, $abs_path_manifest_file);
	#&debugger_print("MIRA manifest: $abs_path_manifest_file\n");
	print "- Calling MIRA now for $keys assembly...\n- It might take a while...\n";
	my $mira_exe = $MIRA_exec." -t $noOfProcesses ".$abs_path_manifest_file;
	$mira_exe .= " > $error_log"; ## show on screen or maybe print to a file

	print "- Sending MIRA command for $keys...\n"; #&debugger_print("MIRA command: $mira_exe\n\n"); 
	my $system_call = system($mira_exe);
	if ($system_call != 0) { DOMINO::printError("Something happened when calling MIRA for assembly reads..."); DOMINO::dieNicely(); }
	print "\n"; &time_log();
	push(@{ $domino_mira_files{$keys}{'contigs_fasta'}}, $domino_mira_files{$keys}{'DIR'}[0]."/".$keys."_d_results/".$keys."_out.unpadded.fasta");
	push(@{ $domino_mira_files{$keys}{'contigs_qual'}}, $domino_mira_files{$keys}{'DIR'}[0]."/".$keys."_d_results/".$keys."_out.unpadded.fasta.qual");
	push(@{ $domino_mira_files{$keys}{'readTagList'}}, $domino_mira_files{$keys}{'DIR'}[0]."/".$keys."_d_info/".$keys."_info_readtaglist.txt");
	push(@{ $domino_mira_files{$keys}{'contigsMIRA'}}, $domino_mira_files{$keys}{'DIR'}[0]."/assembly_id-".$keys.".contigs-MIRA.fasta");

	#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mira_files); print "\n";

	## Obtain reads/contigs identified as repeats
	## Extract read repeats
	
	chdir $domino_mira_files{$keys}{'DIR'}[0];
	my $clean_folder_path = DOMINO::get_earliest("clean_data", $path);
	unless (-d $clean_folder_path) {mkdir $clean_folder_path, 0755;}
	my $preclean_folder = $clean_folder_path."/preRepeatScanner_DOMINO_clean_data";
	unless (-d $preclean_folder) { mkdir $preclean_folder, 0755; }
	
	my (@array_repetitive_reads, @array_repetitive_contigs);
	open(IN, $domino_mira_files{$keys}{'readTagList'}[0]) or DOMINO::printError( "Couldn't open read tag list for $keys taxa for reading and discarding repeats");
	while(<IN>) {
		# do whatever needs to be done
		chomp;
		my $line=$_;
		next if ($line=~/^#/); next if ($line=~/^\s$/);
#file ...readtaglist.txt
##
# conName       cFromPadded     cToPadded       cFromUnpadded   cToUnpadded     type    rName   rFromPadded     rToPadded       rFromUnpadded   rToUnpadded     comment
#
##clean_reads.id-sp3_c1   546     453     546     453     HAF3    Seq652_sp3      52      145     52      145     
##clean_reads.id-sp3_c1   452     452     452     452     HAF2    Seq652_sp3      146     146     146     146     
##clean_reads.id-sp3_c1   451     396     451     396     HAF6    Seq653_sp3      147     202     147     202     

		my @fields = split(/\s+/, $line);
		my ($contig, $flag, $read, @array_reads);
		$contig = $fields[0]; $flag = $fields[5]; $read = $fields[6];
		if (($flag eq 'HAF6') || ($flag eq 'HAF7') || ($flag eq 'MNRr')){ ## Discard reads identified as highly repetitive
			push (@array_repetitive_reads, $read);
			push (@array_repetitive_contigs, $contig);
	}}
	close(IN);

	## Print contigs to discard
	my $accns_file_contigs =  "ids2remove_contigs_".$keys."tab"; open(ACCNS_C, ">$accns_file_contigs");
	my @array_repetitive_contigs_sort = sort(@array_repetitive_contigs); my @array_repetitive_contigs_sort_uniq = uniq(@array_repetitive_contigs_sort);
	for (my $i=0; $i < scalar @array_repetitive_contigs_sort_uniq; $i++) { print ACCNS_C $array_repetitive_contigs_sort_uniq[$i]."\n"; } close(ACCNS_C);
	
	## Print reads to discard
	my $accns_file_reads = "ids2remove_reads_".$keys.".tab"; open(ACCNS_R, ">$accns_file_reads"); 
	my @array_repetitive_reads_sort = sort(@array_repetitive_reads); my @array_repetitive_reads_sort_uniq = uniq(@array_repetitive_reads_sort);
	for (my $j=0; $j < scalar @array_repetitive_reads_sort_uniq; $j++) { print ACCNS_R $array_repetitive_reads_sort_uniq[$j]."\n"; } close(ACCNS_R);

	## Discard reads and contigs using mothur
	if (-z $accns_file_reads) { print "\n+ No reads discarded as there were no repeats identified during MIRA assembly...\n";		
	} else {
		for (my $j=0; $j < scalar @species_files; $j++) {
			print "\n+ Calling mothur executable for discarding reads in file $species_files[$j]...\n\n";
			DOMINO::mothur_remove_seqs($domino_mira_files{$keys}{'original'}[$j], 'YES', $domino_mira_files{$keys}{'DIR'}[0], $accns_file_reads, $mothur_path);
			my $folder_files = DOMINO::readDir($domino_mira_files{$keys}{'DIR'}[0]);
			my @file = grep /.*pick\.fastq/, @{$folder_files};
			File::Copy::move($domino_mira_files{$keys}{'original'}[$j], $preclean_folder);
			File::Copy::move($file[0], $domino_mira_files{$keys}{'original'}[$j]);
	}}
	my $contigs2cluster;
	if (-z $accns_file_contigs) {
		print "+ No contigs discarded as there were no repeats identified during MIRA assembly...\n";
		$contigs2cluster = $domino_mira_files{$keys}{'contigs_fasta'}[0];		
	} else {
		print "\n+ Calling mothur executable for discarding contigs in file ".$domino_mira_files{$keys}{'contigs_fasta'}[0]."...\n";
		DOMINO::mothur_remove_seqs($domino_mira_files{$keys}{'contigs_fasta'}[0], 'NO', $domino_mira_files{$keys}{'DIR'}[0], $accns_file_contigs, $mothur_path);
		my $folder_files = DOMINO::readDir($domino_mira_files{$keys}{'DIR'}[0]);
		my @file = grep /.*pick\.fasta/, @{$folder_files};
		$contigs2cluster = $file[0];
	}	
	
	## Use BLAST for clustering sequences
	print "\n\n+ Generate a BLAST database for $contigs2cluster...\n"; 
	my ($db_generated, $db_generated_message) = DOMINO::makeblastdb($contigs2cluster, $BLAST, $error_log);
	#&debugger_print($db_generated_message);
	my $blast_search = "blast_search.txt";
	print "+ BLAST search now...\n"; 
	my ($blastn, $blastn_message) = DOMINO::blastn($contigs2cluster, $db_generated, $blast_search, $BLAST);
	#&debugger_print($blastn);
	if ($blastn != 0) { DOMINO::printError("BLASTN failed...\n"); exit(); } 
	my $contig_length_Ref = DOMINO::readFASTA_hash($contigs2cluster);
	my $perc_aln_desired = 0.85; my $iden_desired = 85;
	
	## Filter BLAST results
	print "+ Filtering BLAST search now...\n";
	my (%contigs_keep, @contigs_seen);
	my $first_hit = 0;
	open (BLAST, $blast_search); while (<BLAST>) {
		my $line = $_;
		chomp $line;
		my ($query, $subject, $perc_iden, $aln, $mismatch, $gap_open, $query_start, $query_end, $subject_start, $subject_end, $eval, $score ) = split("\t", $line);
		my $length_sub = length($$contig_length_Ref{$subject});
		my $perc_aln = $length_sub*$perc_aln_desired;	
		if ($query eq $subject) { next;  }
		if ($eval < 1e-50 && $aln >= $perc_aln && $perc_iden >= $iden_desired) { 
			if ($first_hit == 0) {
				$first_hit++;
				push(@{$contigs_keep{$query}}, $subject); push (@contigs_seen, $subject);
			} else {
				my $flag = 0;
				foreach my $keys (keys %contigs_keep) {
					if (grep /.*$subject.*/, @{$contigs_keep{$keys}}) { $flag = 1; last; }
					if (grep /.*$query.*/, @{$contigs_keep{$keys}}) { $flag = 1; last; }
				}
				if ($flag == 0) { 
					push(@{$contigs_keep{$query}}, $subject); push (@contigs_seen, $subject);	
	}}}}
	close (BLAST);
	
	foreach my $keys (sort keys %$contig_length_Ref) { unless (grep /$keys/, @contigs_seen) { push (@{ $contigs_keep{$keys}}, $keys); } }
	my $clean_contigs = $domino_mira_files{$keys}{'contigsMIRA'}[0];
	open (CLN, ">$clean_contigs");
	foreach my $keys (sort keys %contigs_keep) {
		my $largest = length($$contig_length_Ref{$keys});
		my $largest_id = $keys;	
		for (my $i=0; $i < scalar @{$contigs_keep{$keys}}; $i++) {
			my $seq = $contigs_keep{$keys}[$i];
			my $len = length( $$contig_length_Ref{$seq} );
			if ($len > $largest) { $largest = $len; $largest_id = $seq; }
		} 
		print CLN ">$largest_id\n$$contig_length_Ref{$largest_id}"; 
	}
	close(CLN);	
	print "+ Done...\n\n";
	chdir $assembly_directory;
	DOMINO::printHeader(" Assembly finished for $keys ","#"); &time_log(); print "\n\n";
	
	my $mira_dump_hash = $domino_mira_files{$keys}{'DIR'}[0]."/dumper_files_threads.txt";
	DOMINO::printDump(\%domino_mira_files, $mira_dump_hash);
}
print "\n\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" MIRA Assembly Step finished ","#");  DOMINO::printHeader("","#"); print "\n\n";
#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mira_files); print "\n";

if ($$hash_parameters{'assembly'}{'CAP3'}) { print "\n\nGetting ready for scaffolding step using CAP3...\n"; mkdir $$hash_parameters{'assembly'}{'CAP3_directory'}[0], 0755; }
## Generates folders and generates link for files of each taxa into them
foreach my $taxa (keys %domino_mira_files) {	
	next if $taxa eq "original"; next if $taxa eq "main";

	my $contigsMIRA_Fasta_file = $domino_mira_files{$taxa}{'contigsMIRA'}[0];
	my $FINAL_fasta_file = $$hash_parameters{"assembly"}{"folder"}[0]."/assembly_id-".$taxa.".contigs.fasta";			
	push( @{ $domino_mira_files{$taxa}{'FINAL'}}, $FINAL_fasta_file);
	
	if ($$hash_parameters{'assembly'}{'CAP3'}) { # If cap3 is used for scaffolding move files and generate folders
		#&debugger_print("DOMINO files"); &debugger_print("Ref", \%domino_mira_files);
		chdir $$hash_parameters{'assembly'}{'CAP3_directory'}[0];
		my $qual_tmp_file = $domino_mira_files{$taxa}{'contigs_qual'}[0];
		my $contigsMIRA_qual_file = $domino_mira_files{$taxa}{'contigsMIRA'}[0].".qual";
		push (@{$domino_mira_files{$taxa}{'contigsMIRA_qual'}}, $contigsMIRA_qual_file);
		print "Fetching qual values for $contigsMIRA_Fasta_file\n";
		my $hash_identifiers_ref = DOMINO::readFASTA_IDSfile($contigsMIRA_Fasta_file);
		my %hash_identifiers = %{$hash_identifiers_ref};
		$/ = ">"; ## Telling perl where a new line starts
		open (QUAL, $qual_tmp_file);
		open (OUT, ">$contigsMIRA_qual_file");
		while (<QUAL>) {
			next if /^#/ || /^\s*$/;
			chomp;
			my ($seq_id, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $seq_id);
			if ($hash_identifiers{$seq_id}) { 
				print OUT ">".$seq_id."\n".$sequence;
				#print ">".$seq_id."\n";
				}
		} $/ = "\n";
		close (QUAL); close (OUT);
																
		## Generate directories
		my $MID_directory = $$hash_parameters{'assembly'}{'CAP3_directory'}[0]."/".$taxa;
		unless (-d $MID_directory) { mkdir $MID_directory, 0755; }                                               
		#&debugger_print("ln -s $contigsMIRA_Fasta_file $MID_directory/"); 
		system("ln -s $contigsMIRA_Fasta_file $MID_directory/"); 
		#&debugger_print("ln -s $contigsMIRA_qual_file $MID_directory/"); 
		system("ln -s $contigsMIRA_qual_file $MID_directory/"); 
		push(@{$domino_mira_files{$taxa}{'CAP3_dir'}}, $MID_directory);
		chdir $MID_directory;
			
		###########################################################################
		###	cap3 Scaffolding of the contigs MIRA generated for each taxa	###
		###########################################################################
		#&debugger_print("DOMINO files"); &debugger_print("Ref", \%domino_mira_files);
		print "\n\n"; DOMINO::printHeader("","#"); 
		DOMINO::printHeader(" CAP3 Assembly Step Started ","#"); DOMINO::printHeader("","#"); print "\n\n";
		my @tmp_array = split("/", $contigsMIRA_Fasta_file);
		my $fasta_name_contigsMIRA = $tmp_array[-1];
		my $command_CAP3 = $CAP3_exec." ".$fasta_name_contigsMIRA." -o ".$$hash_parameters{'assembly'}{'overCAP3'}[0]." -p ".$$hash_parameters{'assembly'}{'simCAP3'}[0]." 2> $error_log";
		#&debugger_print("CAP3 command: $command_CAP3\n"); 
		my $system_CAP3 = system($command_CAP3);
		if ($system_CAP3 != 0) { DOMINO::printError("Some error happened when calling CAP3 for assembly reads..."); DOMINO::dieNicely(); }
		
		###########################################
		### Get contigs and singlets assembled 	###
		###########################################
		my $singlets_file_CAP3 = $fasta_name_contigsMIRA.".cap.singlets"; push (@{$domino_mira_files{$taxa}{'singletsCAP3'}}, $singlets_file_CAP3);
		my $contigs_file_CAP3 = $fasta_name_contigsMIRA.".cap.contigs"; push (@{$domino_mira_files{$taxa}{'contigsCAP3'}}, $contigs_file_CAP3);

		my $tmp_fasta = "tmp.fasta"; open (OUT_fasta, ">$tmp_fasta"); 
		open (CONTIGS, $contigs_file_CAP3); while (<CONTIGS>) { print OUT_fasta $_; } close (CONTIGS);
		open (SINGLETS, $singlets_file_CAP3); while (<SINGLETS>) { print OUT_fasta $_; } close (SINGLETS);close (OUT_fasta);

		#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mira_files); print "\n";
		&change_seq_names($tmp_fasta, $FINAL_fasta_file, $taxa);
		File::Copy::move($FINAL_fasta_file, $dirname);
		print "\n"; &time_log();			
	} else { 
		&change_seq_names($contigsMIRA_Fasta_file, $FINAL_fasta_file, $taxa);
	}
	
	my $dump_hash = $domino_mira_files{$taxa}{'DIR'}[0]."/dumper_files_threads.txt";
	DOMINO::printDump(\%domino_mira_files, $dump_hash);
}

###############################
## Generating some statistics #
###############################
chdir $dirname;
my $array_Ref = DOMINO::readDir($dirname);
my @dirname_files = @$array_Ref;
foreach my $taxa (keys %domino_mira_files) {
	next if $taxa eq "original"; next if $taxa eq "main";
	print "\n+ Generating some statistics for Assembly file: $domino_mira_files{$taxa}{'FINAL'}[0]\n";
	if ($$hash_parameters{'assembly'}{'CAP3'}) {
		print "+ Statistics for MIRA assembly:\n";
		my $stats_file_MIRA = DOMINO::Contig_Stats($domino_mira_files{$taxa}{'contigsMIRA'}[0]);
		push (@{ $domino_mira_files{$taxa}{'MIRA_stats'} }, $stats_file_MIRA);
		print "+ Statistics for CAP3 scaffolding:\n";
	}
	my $stats_file = DOMINO::Contig_Stats($domino_mira_files{$taxa}{'FINAL'}[0]);
	print $stats_file."\n";
	push (@{ $domino_mira_files{$taxa}{'FINAL_stats'} }, $stats_file);	
	print "\n\n";
}

DOMINO::print_success_Step("MIRA");

sub Generate_manifest_file {
	
	##########################################################################################
	##  This function generates a manifest file for Illumina or 454 for MIRA assembler		##
	##	Jose Fco. Sanchez Herrero, 10/06/2014 jfsanchezherrero@ub.edu						##
	##########################################################################################

	my $name_of_project = $_[0]; my $technology = $_[1];
	my $fastq_file = $_[2];	my $fastq_file2 = $_[3];
	
	DOMINO::printHeader(" Generate MIRA manifest file for $name_of_project ","%"); 
	my $manifest_file = "manifest_$name_of_project.txt";
	print "- Generating the manifest file $manifest_file...\n";
	
	open (OUT, ">$manifest_file");
	print OUT "# Manifest file for $technology Project $name_of_project to use MIRA\n";
	print OUT "# First part: defining some basic things\n";
	print OUT "project = $name_of_project\n";
	print OUT "job = genome, denovo, accurate\n";
	
	my $parameter = "parameters = --hirep_something -NW:cnfs=warn\\\n";  
	$parameter .= "\tCOMMON_SETTINGS -GE:not=$noOfProcesses:amm=on:kpmf=20 -OUT:ors=no:orc=no:otc=no:orw=no:rtd=yes -CL:ascdc=no\\\n";

	my $mrs = $$hash_parameters{'assembly'}{'mrs'}[0];	
	if ($technology eq "454") { 
		$parameter .= "\t454_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	} else { 
		$parameter .= "\tSOLEXA_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	}

	if ($fastq_file2) { ## PE reads
		print OUT "readgroup = ".$name_of_project."_reads\nautopairing\n";
		print OUT "data = $fastq_file $fastq_file2\n";
	} else {
		print OUT $parameter."\nreadgroup = ".$name_of_project."_reads\ndata = $fastq_file\n";
	}

	if ($technology eq "454") { 
		print OUT "technology = 454\n";
	} else { 
		print OUT "technology = solexa\n"; }
	close (OUT);
	
	return $manifest_file;	
}

sub change_seq_names {
	
	my $file = $_[0];
	my $name_file = $_[1];
	my $identifier = $_[2];
	
	$/ = ">"; ## Telling perl where a new line starts
	open (FILE, $file) or die "[Can not open file \"$file\": DOMINO_Scripts: $0 line $.\n";
	open (OUT, ">", "$name_file");
	my $count;
	while (<FILE>) {
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($seq_id, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $seq_id);
    	$count++;
    	$seq_id = ">Contig_".$count."_".$identifier;
    	
    	my @array_seq = split ("\n", $sequence);
    	my $string_seq;
    	for (my $i = 0; $i < scalar @array_seq; $i++) { $string_seq .= $array_seq[$i]; }
		print OUT $seq_id."\n".uc($string_seq)."\n";    	
	}
	close (FILE); close (OUT);	$/ = "\n";
}

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}
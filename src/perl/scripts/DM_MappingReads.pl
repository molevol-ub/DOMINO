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
my $samtools_path = $scripts_path."samtools-1.3.1/samtools";
my $bowtie_path = $scripts_path."bowtie2-2.2.9/";
my $domino_Scripts = $scripts_path."scripts";

#################################################################

## Arguments
my $path = $ARGV[0];
my $step_time = $ARGV[1];
#################################################################

#################################################################
## Get general parameters and files
my $hash_parameters = DOMINO::get_parameters($path."/", "mapping");
my %domino_mapping_files = %{ DOMINO::get_DOMINO_files($path."/", "mapping") };
#print Dumper $hash_parameters; #print Dumper \%domino_mapping_files; #exit();

my $align_dirname = $$hash_parameters{'mapping'}{'folder'}[0];
my $num_proc_user = $$hash_parameters{'mapping'}{'cpu'}[0];
my $noDiscard = $$hash_parameters{'mapping'}{'noDiscard'}[0];
my %new_domino_files;
#################################################################

#######################################################
### Get Pre-assemble taxa read contigs of each taxa ### 
#######################################################
print "\n"; DOMINO::printHeader(" Get FASTQ files of the contigs generated ", "%"); print "\n";
#&debugger_print("Change dir to: ".$align_dirname);

## Mapping of the reads, all taxa used as reference
foreach my $reference_identifier (sort keys %domino_mapping_files) {
	unless ($domino_mapping_files{$reference_identifier}{'contigs'}) { next; }
	chdir $align_dirname; #&debugger_print("Change dir to: ".$align_dirname);

	my (@sam_files, @clean_sam_files, @sorted_bam); 
		
	## Get the name and identifier of each reference fasta used	
	my @temp_contigs_name = split ("/", $domino_mapping_files{$reference_identifier}{'contigs'}[0]);
	my $contigs_fasta = $temp_contigs_name[$#temp_contigs_name];
	my $ref_Fasta = $reference_identifier.".fasta";
	
	###############################
	###	 Copy necessary files	### 
	###############################
	print "+ Generating symbolic links for necessary files\n"; 
	my $dir = $align_dirname."/".$reference_identifier; mkdir $dir, 0755; chdir $dir;
	system("ln -s $domino_mapping_files{$reference_identifier}{'contigs'}[0] $ref_Fasta");
	push (@{ $domino_mapping_files{$reference_identifier}{'dir'} }, $dir);
	## Generate a directory for each one
	print "+ Using as reference: $contigs_fasta\tID: $reference_identifier...OK\n";
	print "+ Generating a new directory $dir....OK\n\n";
	#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mapping_files); 

	###################################
	###		 Index Contig file		### 
	###################################
	DOMINO::printHeader(" Indexing Contig File for mapping Reference ", "%");
	# Index contig reference file using Bowtie
	my $reference_tag = "reference_".$reference_identifier;
	print "- Reference: $ref_Fasta...\n";
	my $bowtie_index_call = $bowtie_path."bowtie2-build --threads $num_proc_user -f ".$domino_mapping_files{$reference_identifier}{'contigs'}[0]." ".$reference_tag;   
	#&debugger_print("BOWTIE2 command: ".$bowtie_index_call);
	my $index_result = system ($bowtie_index_call);
	if ($index_result != 0) {
		DOMINO::printError("Exiting the script. Some error happened when calling bowtie for indexing the file...\n", $$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0]); 
		DOMINO::print_fail_Step("index_genome_$reference_tag");
		DOMINO::dieNicely();
	} else { DOMINO::print_success_Step("index_genome_$reference_tag"); #$domino_success_steps{$reference_tag}{"index_genome"}++;
	}
	print "\n"; &time_log();	print "\n";
	push (@{ $domino_mapping_files{$reference_identifier}{"index_genome"} }, $dir."/".$reference_tag); ## get symbolic link if it is already indexed

	###########################
	###	Align taxa Reads	### 
	###########################
	print "\n";	DOMINO::printHeader(" Aligning Reads Individually ", "%"); print "\n";
	my @clean_fastq_files;
	my @tmp_array;
	if ($$hash_parameters{'mapping'}{'user_cleanRead_files'}[0]) { print "+ User clean reads files would be mapped\n";
	} elsif ($$hash_parameters{'mapping'}{'map_contig_files'}[0]) { 			print "+ Contig files would be mapped\n";
	} else { 								print "+ Clean reads files would be mapped\n"; ## Map DOMINO clean reads
	}
	print "+ Obtain information of the reference sequences\n";
	my ($reference_hash_fasta_ref, $message) = DOMINO::readFASTA_hashLength($domino_mapping_files{$reference_identifier}{'contigs'}[0]); ## Obtain reference of a hash
	my $file2dump_seqs = $dir."/contigs_".$reference_identifier."_length.txt";
	push (@{ $domino_mapping_files{$reference_identifier}{"hash_reference_file"} }, $file2dump_seqs);
	DOMINO::printDump($reference_hash_fasta_ref,$file2dump_seqs,1);
	
	## Set DOMINO to use as many CPUs as provided
	## maximize process if > 10
	my ($split_CPU, $subprocesses);
	my $number_sp = $$hash_parameters{'mapping'}{'number_sp'}[0];
	if ($num_proc_user > 10) { 	
		$split_CPU = int($num_proc_user/$number_sp);  
		$subprocesses = $number_sp; 
		if ($split_CPU == 0) { $split_CPU = 1;}
	} else {
		$split_CPU = $num_proc_user; $subprocesses = 1; 			
	}			
	
	## Get files for later dump
	foreach my $reads_here (sort keys %domino_mapping_files) {
		unless ($domino_mapping_files{$reads_here}{'reads'}) { next; }
		my $file2dump = $dir."/dump_file_mapping_".$reads_here."_Ref_".$reference_identifier.".txt";
		#print "File to dump: ".$file2dump."\n";
		push (@{ $domino_mapping_files{$reads_here}{"DUMP_Mapping::Ref:".$reference_identifier} }, $file2dump);
	}
	
	my $pm_read_Reference =  new Parallel::ForkManager($subprocesses); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	my $total_subprocesses = $number_sp; my $count_subprocesses=0;
	print "+ Mapping process would be divided into $number_sp parts using up to $split_CPU/$num_proc_user CPUs provided\n";
	foreach my $reads (sort keys %domino_mapping_files) {
		unless ($domino_mapping_files{$reads}{'reads'}) { next; }
		$count_subprocesses++;
		print "\t- Mapping reads ($reads) vs reference ($reference_identifier): [$count_subprocesses/$total_subprocesses]\n";
		my $out_log_file = $dir."/reads_".$reads."-reference_".$reference_identifier."_logfile.txt";
		open (LOG, ">$out_log_file");
		chdir $domino_mapping_files{$reference_identifier}{'dir'}[0];
		print LOG "\nChange dir to: ".$domino_mapping_files{$reference_identifier}{'dir'}[0]."\n";
		my $pid = $pm_read_Reference->start() and next;
		my %domino_mapping_files_split_mapping;
		push(@{ $domino_mapping_files_split_mapping{$reads}{'LOG'}}, $out_log_file);
		
		## Mapping Parameters
		my $R_group_id = '--rg-id '.$reads;
		my $R_group_name = ' --rg '.$reads;
		my $threads = ' -p '.$split_CPU;
		my $mismatches = ' -N 1 --np 0'; ## Do not add penalty if read/ref got an ambiguous base
		my $read_gap_open = ' --rdg '.$$hash_parameters{'mapping'}{'rdgopen'}[0].','.$$hash_parameters{'mapping'}{'rdgexten'}[0];
		my $ref_gap_open = ' --rfg '.$$hash_parameters{'mapping'}{'rfgopen'}[0].','.$$hash_parameters{'mapping'}{'rfgexten'}[0];
		my $mismatch_penalty = ' --mp '.$$hash_parameters{'mapping'}{'mis_penalty'}[0];
		my $mapping_file = $domino_mapping_files{$reads}{'reads'}[0];
		my $botwie_system = $bowtie_path."bowtie2";
		if ($$hash_parameters{'mapping'}{'bowtie_local'}[0]) { $botwie_system .= " --local"; }
		my $sam_name = $dir."/".$reference_tag."-taxa_".$reads.".sam";
		
		print LOG "+ Aligning reads for $reads against $reference_identifier as reference\n";
		if ($$hash_parameters{'mapping'}{'type_input'}[0] eq 'pair_end') {
			my $second_Read_file = $domino_mapping_files{$reads}{'reads'}[1];
			print LOG "+ Reads 1: $mapping_file\n+ Reads 2: $second_Read_file\n";
			$botwie_system .= " -x ".$reference_tag." -q -1 $mapping_file -2 $second_Read_file -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
		
		} elsif ($$hash_parameters{'mapping'}{'type_input'}[0] eq 'single_end') { ## Illumin single end, 454
			print LOG "+ Reads 1: $mapping_file...\n";
			if ($$hash_parameters{'mapping'}{'map_contig_files'}[0]) { ## Mapping contigs
				$botwie_system .= " -x ".$reference_tag." -f -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
			} else { 
				$botwie_system .= " -x ".$reference_tag." -q -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
		}}
		print LOG "+ SAM file: $sam_name\n+ Mapping now...\n\n"; #&debugger_print("BOWTIE2 command: ".$botwie_system); 
		
		### Map Reads
		my $error_bowtie = $dir."/reads_".$reads."-reference_".$reference_identifier."_mapping_logfile.txt";
		print LOG "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nMapping Statistics:\n\n"; 
		$botwie_system .= " 2> ".$error_bowtie;
		my $system_bowtie_call = system ($botwie_system);
		if ($system_bowtie_call != 0) {
			DOMINO::printError("Exiting the script. Some error happened when calling bowtie for mapping the file $mapping_file...\n", $$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0]); 
			DOMINO::print_fail_Step("mapping_$reads");
			DOMINO::dieNicely(); 
		} else { DOMINO::print_success_Step("mapping_$reads"); #$domino_success_steps{$reference_tag}{"map_$reads"}++
		}
								
		push (@{$domino_mapping_files_split_mapping{$reads}{"SAM::Ref:".$reference_identifier}}, $sam_name);
		open (IN, $error_bowtie); while (<IN>) {print LOG $_; } close(IN);
		print LOG "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
		print LOG "\n\n+ Mapping finished for taxa $reads against $reference_identifier\n";		
		
		#################################
		## Get the the reference fasta ##
		#################################
		print LOG "\n+ Checking mapping reads in ".$sam_name."...\n";
		
		## Generate a sam for each contig
		my @temp = split ("\.sam", $sam_name);
		chdir $dir; my $dir_tmp = $temp[0]."_SPLIT"; mkdir $dir_tmp, 0755; chdir $dir_tmp;
		print LOG "\nChange dir to: ".$dir_tmp."\n";
		#print Dumper $reference_hash_fasta_ref;
		my @number_contigs = sort (keys %$reference_hash_fasta_ref);
		my $scalar = scalar @number_contigs;
		print LOG "+ This SAM file contains $scalar referense sequences...\n";
		system("ln -s $sam_name"); 
		my @temp_name = split ("/", $sam_name); my $sam_base_name = $temp_name[-1];
		
		my $sorted_bam_file = DOMINO::generate_bam($sam_name, $split_CPU);
		DOMINO::generate_index_bam($sorted_bam_file);

		if ($scalar == 1) {
			push (@{ $domino_mapping_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $temp_name[-1]);
		} else {
			print LOG "+ Using parallel threads ($split_CPU CPUs)...\n";
			print LOG "+ Splitting SAM file into several parts to speed the computation...\n"; 
			my $parts = int($scalar/$split_CPU); ## Split SAM into as many CPUs provided
			my @commands; my $iteration = 0; 
			while (1) {
				if (@number_contigs) {
					my @array1 = splice(@number_contigs, 0, $parts + 1);
					my $string = join(" ", @array1);
					my @temp_1 = split ("\.sorted.bam", $sorted_bam_file);
					my @temp_2 = split ("/", $temp_1[0]);
					my $sam_file_part = $dir_tmp."/".$temp_2[-1]."_part-".$iteration.".sam";
					#print "Contigs: ".$string."\n";
					my $command = $samtools_path." view -@ $split_CPU -Sh -o $sam_file_part $sorted_bam_file $string";
					push (@{ $domino_mapping_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $sam_file_part);
					push (@commands, $command); $iteration++; 
			} else { last; }}
			
			my $pm_SAM_split =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			my $total_commands = scalar @commands;
			my $count_commands=0;
			for (my $a=0; $a < scalar @commands; $a++) {
				$count_commands++;
				print LOG "\t- Splitting SAM $sam_base_name: [$count_commands/$total_commands]\n";
				my $pid = $pm_SAM_split->start() and next; 
				#&debugger_print("SAMTOOLS command: $commands[$a]");	
				my $system_call = system ($commands[$a]);
				$pm_SAM_split->finish(); # pass an exit code to finish
			}
			$pm_SAM_split->wait_all_children;
			#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mapping_files_split_mapping); 		
		}
		
		### Remove multimapping reads	### 
		##  DOMINO checks the SAM files generated and discards bad reads mapping, unmapping reads or multimapping reads. ##
		my @array_files_split = @{ $domino_mapping_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier}};				
		print LOG "\n+ Cleaning reads now in parallel threads ($split_CPU CPUs)...\n";
		## Get files for later dump
		for (my $i=0; $i < scalar @array_files_split; $i++) {
			my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part_".$i.".txt";
			#print "File to dump: ".$file2dump."\n";
			push (@{ $domino_mapping_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier} }, $file2dump);
		}		
		
		my $pm_SAM_parts =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
		my $scalar_array_files_split = scalar @array_files_split; my $count_split=0;
		for (my $i=0; $i < scalar @array_files_split; $i++) {
			$count_split++;
			print LOG "\t- Checking each splitted file for $sam_base_name: [$count_split/$scalar_array_files_split]\n";	
			my $pid = $pm_SAM_parts->start() and next;
			my %domino_mapping_files_SAM_parts; my $discard_reads = 0; my $good_reads = 0; my $total_reads = 0;
			open (SAM, "<$array_files_split[$i]");
			my @temp = split ("\.sam", $array_files_split[$i]);
			my $output_sam = $temp[0]."_clean.sam"; 
			push (@{ $domino_mapping_files_SAM_parts{$reads}{"CLEAN_SAM_Parts::Ref:".$reference_identifier} }, $output_sam);
			open (SAM_OUT, ">$output_sam");	
			while (<SAM>) {
				chomp; my $line = $_;
				if ($line =~ /^\@.*/) { print SAM_OUT $line."\n"; next; }		
				my @sam = split ("\t", $line); 
				$total_reads++;				
				if ($line =~ /.*SA:Z:.*/) {   ## Discard Multimapping reads
					$discard_reads++; 
				} elsif ($sam[2] eq '*') { ## Discard unmapped reads
					## When calling Bowtie2 we use --no-unal not to show unaligned reads
					## but still we would check it, also for user input BAM files			
					$discard_reads++;
				} else { ## Check if read is efficiently mapping 
					## Even Bowtie2 is mapping end to end and allowing no clipping
					## we would check if any hard or soft clipping is introduced
					## (inherited from a previous version using BWA)
					
					my $cigar=$sam[5]; my $xcent = 0; my $xcentmax = 0; my $pos = 0; my $NOpos = 0;
					my $int= $$hash_parameters{'marker'}{'cigar_pct'}[0];
					my $cigar_pct_convert = $int/100;
					
					#&debugger_print("\nLINE: ".$line."\nCIGAR: ".$cigar."\nMAX_CIGAR: ".$cigar_pct_convert."\n");
					while ($cigar !~ /^$/){         		
						if ($cigar =~ /^([0-9]+[MIDSH])/){
							my $cigar_part = $1;								
							#&debugger_print("CIGAR Part: $cigar_part\n");								
							if ($cigar_part =~ /(\d+)M/){ $pos += $1;
							} elsif ($cigar_part =~ /(\d+)I/){ $NOpos += $1;
							} elsif ($cigar_part =~ /(\d+)D/){ $NOpos += $1;
							} elsif ($cigar_part =~ /(\d+)S/){ $NOpos+= $1;
							} elsif ($cigar_part =~ /(\d+)H/){ $NOpos+= $1;
							} else { # die "Unexpected cigar: $cigar\n";
							}# close if
							$cigar =~ s/$cigar_part//;
					}}
					my $total = $pos + $NOpos; 				 ## total positions
					$xcentmax = $total * $cigar_pct_convert; ## Max bad CIGAR
					$xcent = $NOpos; ## Bad CIGAR
					if($xcent <= $xcentmax){
						print SAM_OUT $line."\n"; $good_reads++;
					} else {
						#&debugger_print("XCENT > XCENTMAX\nDISCARD READ!!\n"); 
						if ($$hash_parameters{'mapping'}{'bowtie_local'}[0]) { print SAM_OUT $line."\n"; next; } 
						$discard_reads++;
			}}} close(SAM); close(SAM_OUT);	
			#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mapping_files_SAM_parts); 
		
			###################################################################
			## Generate sorted bam files in order to be able to get coverage ##
			###################################################################
			my $clean_sorted_bam = DOMINO::generate_bam($array_files_split[$i], $split_CPU);
			push (@{ $domino_mapping_files_SAM_parts{$reads}{"clean_BAM_Parts::Ref:".$reference_identifier} }, $clean_sorted_bam);

			## Generate Coverage statistics for the alignment file
			my @tmp_bam_name = split ("\.sorted.bam", $clean_sorted_bam);
			my $coverage_file = $tmp_bam_name[0]."_coverage_stats.txt";
			push (@{ $domino_mapping_files_SAM_parts{$reads}{"coverage_Parts::Ref:".$reference_identifier} }, $coverage_file);

			my $coverage_samtools_command = $samtools_path." depth ".$clean_sorted_bam." > ".$coverage_file;
			my $system_coverage_call = system ($coverage_samtools_command); 
			#&debugger_print("SAMTOOLS command: $coverage_samtools_command");
			if ($system_coverage_call != 0) {
				DOMINO::printError("Exiting the script. Some error happened when calling SAMtools for obtaining coverage of file $sorted_bam[$i]...\n", $$hash_parameters{'mapping'}{'mapping_markers_errors_details'}[0]); DOMINO::dieNicely();
			}
			
			unless (-z $coverage_file) {
				#print "\t- Filtering Coverage Stats...\n\n";
				my ($sum_coverage_each, $total_positions_each);
				my %max_cov;
				open (COVERAGE, "<$coverage_file");
				while (<COVERAGE>) {
					my $line = $_;
					chomp $line;
					my @array = split (/\s+/,$line);
					if ($array[2]) {
						my $coverage_position = $array[2];
						$sum_coverage_each += $coverage_position;
						$total_positions_each++;
						my $contig = $array[0];
						if (defined($max_cov{$contig})) {
							if ($coverage_position > $max_cov{$contig}) { $max_cov{$contig} = $coverage_position; }
						} else { $max_cov{$contig} = $coverage_position; }
				}}
				close(COVERAGE);
				
				## Print max coverage to a file
				my $max_cov_file = $coverage_file."_max_coverage";
				push (@{$domino_mapping_files_SAM_parts{$reads}{"max_coverage_Parts::Ref:".$reference_identifier} }, $max_cov_file);
				open (OUT_COV, ">$max_cov_file");
				foreach my $contig (sort keys %max_cov) {
					print OUT_COV $contig.":".$max_cov{$contig}."\n";
				} close (OUT_COV);
				undef %max_cov;
				
				## Push useful info
				my $coverage_each = $sum_coverage_each.":".$total_positions_each;
				push (@{ $domino_mapping_files_SAM_parts{$reads}{"mean_coverage_Parts::Ref:".$reference_identifier} }, $coverage_each);
				push (@{$domino_mapping_files_SAM_parts{$reads}{"discard_reads::Ref:".$reference_identifier} }, $discard_reads);
				push (@{$domino_mapping_files_SAM_parts{$reads}{"good_reads::Ref:".$reference_identifier} }, $good_reads);
				push (@{$domino_mapping_files_SAM_parts{$reads}{"total_reads::Ref:".$reference_identifier} }, $total_reads);

				# Dump info into file
				DOMINO::printDump(\%domino_mapping_files_SAM_parts, $domino_mapping_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier}[$i]);
			}
			$pm_SAM_parts->finish(); # pass an exit code to finish
		}
		$pm_SAM_parts->wait_all_children;
		my @dump_files1 = @{ $domino_mapping_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier} };
		my $hash_dump1 = DOMINO::retrieve_info(\@dump_files1, \%domino_mapping_files_split_mapping);
		%domino_mapping_files_split_mapping = %{$hash_dump1};

		## Get total reads: discard, good and total
		my $discard_reads = 0; my $good_reads = 0; my $total_reads = 0; my @array;
		push (@array, $domino_mapping_files_split_mapping{$reads}{"discard_reads::Ref:".$reference_identifier});
		push (@array, $domino_mapping_files_split_mapping{$reads}{"good_reads::Ref:".$reference_identifier});
		push (@array, $domino_mapping_files_split_mapping{$reads}{"total_reads::Ref:".$reference_identifier});
		for (my $i=0; $i < scalar @array; $i++) {
			my @array_2 = @{$array[$i]};
			for (my $j=0; $j < scalar @array_2; $j++) {
				if ($i==0) { $discard_reads += $array_2[$j];
				} elsif ($i==1) { $good_reads += $array_2[$j];
				} else { $total_reads += $array_2[$j];
		}}}		
		
		##############################################################################################	
		## Calculate the probability of being under a poisson distribution with the mean of our data
		##############################################################################################
		my @stats = @{ $domino_mapping_files_split_mapping{$reads}{"mean_coverage_Parts::Ref:".$reference_identifier}};
		my ($sum_coverage, $total_positions);
		for (my $i=0; $i < scalar @stats; $i++) {
			my @stats_each = split(":", $stats[$i]);
			$sum_coverage += $stats_each[0];
			$total_positions += $stats_each[1];
		}
		my $mean_coverage = $sum_coverage/$total_positions;
		my $mean = sprintf ("%.3f", $mean_coverage);

		my @max_cov_files = @{ $domino_mapping_files_split_mapping{$reads}{"max_coverage_Parts::Ref:".$reference_identifier} };
		my $num_contigs; my %discard_contigs; my $contigs_discarded;
		for (my $h=0; $h < scalar @max_cov_files; $h++) {
			## Read Coverage file
			open (COV_READ, $max_cov_files[$h]);
			while (<COV_READ>) {
				chomp;
				my $line = $_;
				my @array = split(":", $line); 
				$num_contigs++; # $array[0] contig id in case we need					
				my $prob_poisson;
				if ($array[1] > 169) { ## Factorial would be out of range!
					$prob_poisson = 0;	
				} else { $prob_poisson = DOMINO::Poisson_distribution($array[1], $mean_coverage); }		
				if ($prob_poisson < $$hash_parameters{'mapping'}{'level_significance_coverage_distribution'}[0]) { # Discard
					$contigs_discarded++; $discard_contigs{$array[0]}++;
				} elsif ($prob_poisson eq 'nan') { # Discard
					$contigs_discarded++; $discard_contigs{$array[0]}++;
				} elsif ($array[1] == 0) {  # Discard
					$contigs_discarded++; $discard_contigs{$array[0]}++;
		}}}

		my $percentage_discard_Reads; my $h; my $perc_contig; my $h_cont;
		if ($discard_reads) { $percentage_discard_Reads = (1 - ($discard_reads/$total_reads))*100; $h = sprintf ("%.3f", $percentage_discard_Reads);} else { $h=100; $discard_reads=0;}
		if ($contigs_discarded) { $perc_contig = (1 - ($contigs_discarded/$num_contigs))*100; $h_cont = sprintf ("%.3f", $perc_contig); } else { $contigs_discarded = 0; $h_cont = 100; }
		
		## Adjust SAM files
		print LOG "\n+ Adjusting the SAM/BAM files using parallel threads ($split_CPU CPUs)\n+ Splitted files would be used...\n";
		my @parts_clean_sam = @{ $domino_mapping_files_split_mapping{$reads}{"CLEAN_SAM_Parts::Ref:".$reference_identifier}};

		## Get files for later dump
		for (my $i=0; $i < scalar @parts_clean_sam; $i++) {
			my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part2_".$i.".txt";
			#print "File to dump: ".$file2dump."\n";
			push (@{ $domino_mapping_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} }, $file2dump);
		}	
		
		my $pm_SAM_PILEUP =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
		my $total_parts_clean_sam = scalar @parts_clean_sam;
		my $count_totalParts=0;
		for (my $j=0; $j < scalar @parts_clean_sam; $j++) {
			$count_totalParts++;
			print LOG "\t- Adjusting SAM/BAM part file for $sam_base_name: [$count_totalParts/$total_parts_clean_sam]\n";
			my @basename = split("/", $parts_clean_sam[$j]);
			my @name = split(".sam", $basename[-1]);
			my $pid = $pm_SAM_PILEUP->start($name[0]) and next;
			my %domino_mapping_files_SAM_PILEUP;
			my @tmp_sam = split("\.clean.sam", $parts_clean_sam[$j]);
			my $sam_filter = $tmp_sam[0]."_filtered.sam";
			open (SAM_OUT, ">$sam_filter"); open (SAM, "<$parts_clean_sam[$j]");
			while (<SAM>) {
				chomp; my $line = $_;
				if ($line =~ /^@.*/ ) { print SAM_OUT $line."\n";	
				} else {
					my @array = split (/\s+/,$line);
					if ($noDiscard) {
						print SAM_OUT $line."\n";	
					} else {
						if (!$discard_contigs{$array[2]}) { 
							print SAM_OUT $line."\n";	
						}
					}
				}
			}
			close(SAM_OUT); close(SAM); undef %discard_contigs;
			#print LOG "- File checked: Contigs and Reads discarded...\n";
			push (@{ $domino_mapping_files_SAM_PILEUP{$reads}{"FILTERED_SAM_Parts::Ref:".$reference_identifier} }, $sam_filter);
			my $bam_filtered_returned = DOMINO::generate_bam($sam_filter, 1);
			push (@{ $domino_mapping_files_SAM_PILEUP{$reads}{"FILTERED_BAM_Parts::Ref:".$reference_identifier} }, $bam_filtered_returned);
			unless ($reads eq $reference_identifier) { ## DO NOT GENERATE FILTER PROFILE FOR REFERENCE
				print LOG "- Generate a PILEUP file for $sam_filter...\n";
				my $dir_returned = &generate_filter_PILEUP($bam_filtered_returned, $domino_mapping_files{$reference_identifier}{'contigs'}[0], $reference_identifier, $reads);
				print LOG "Finish PILEUP for $bam_filtered_returned\n";
				push (@{ $domino_mapping_files_SAM_PILEUP{$reads}{"PROFILE::Ref:".$reference_identifier} }, $dir_returned);
			}
			# Dump info into file
			DOMINO::printDump(\%domino_mapping_files_SAM_PILEUP, $domino_mapping_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier}[$j]);
			$pm_SAM_PILEUP->finish($name[0]); # pass an exit code to finish
		}
		$pm_SAM_PILEUP->wait_all_children;
		
		my @dump_files = @{ $domino_mapping_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} };
		my $hash_dump = DOMINO::retrieve_info(\@dump_files, \%domino_mapping_files_split_mapping);
		%domino_mapping_files_split_mapping = %{$hash_dump};

		print LOG "\n\n\n";			
		my $stats_file = $dir."/mapping_statistics_Ref_".$reference_identifier."_Reads_".$reads.".txt";
		push (@{ $domino_mapping_files_split_mapping{$reference_identifier}{'stats'} }, $stats_file);
		print LOG "================ Filtering Statistics ======================\n"; 
		open (STATS, ">$stats_file");			
		print STATS "File: $sam_name\n"; 																print LOG "File: $sam_name\n";
		print STATS "Reference: ".$reference_identifier."\n"; 											print LOG "Reference: ".$reference_identifier."\n";
		print STATS "Reads: ".$reads."\n";																print LOG "Reads: ".$reads."\n";
		print STATS "==========================================================\n";						print LOG "==========================================================\n";
		print STATS "Total Contigs: $num_contigs\n";													print LOG "Total Contigs: $num_contigs\n";
		print STATS "Total Reads mapping: $total_reads\n";												print LOG "Total Reads mapping: $total_reads\n";
		print STATS "Coverage Mean: x".$mean."\n";														print LOG "Coverage Mean: x".$mean."\n";
		print STATS "\n******* CONTIGS *******\n";														print LOG "\n******* CONTIGS *******\n";
		print STATS "Contigs discarded: ".$contigs_discarded."\n";										print LOG "Contigs discarded: ".$contigs_discarded."\n";
		print STATS "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n"; 	print LOG "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n";
		print STATS "\n******* READS *******\n";														print LOG "\n******* READS *******\n";
		print STATS "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";	print LOG "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";
		print STATS "Reads remaining: $good_reads\t( ".$h." %)\n";										print LOG "Reads remaining: $good_reads\t( ".$h." %)\n";
		print LOG "============================================================\n"; 
		print STATS "============================================================\n"; 
		print LOG "\n\n\n"; close (STATS);			

		DOMINO::printDump(\%domino_mapping_files_split_mapping, $domino_mapping_files{$reads}{"DUMP_Mapping::Ref:".$reference_identifier}[0]);
		$pm_read_Reference->finish(); # pass an exit code to finish
	} ## foreach reads
	$pm_read_Reference->wait_all_children;			
	
	
	## Get files for later dump
	foreach my $reads_here (sort keys %domino_mapping_files) {
		unless ($domino_mapping_files{$reads_here}{'reads'}) { next; }
		my @dump_files = @{ $domino_mapping_files{$reads_here}{"DUMP_Mapping::Ref:".$reference_identifier} };
		%new_domino_files = %{ DOMINO::retrieve_info(\@dump_files, \%domino_mapping_files) }; 
		#print Dumper \%new_domino_files;
	}
	#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_mapping_files);
	print "\n\n"; DOMINO::printHeader("", "+");DOMINO::printHeader(" Mapping finished for Reference $reference_identifier ", "+");DOMINO::printHeader("", "+"); &time_log(); print "\n";		
	DOMINO::print_success_Step("mapping_ref_".$reference_identifier);

} # foreach reference

### dump info
my %hash = %{ DOMINO::get_uniq_hash(\%new_domino_files) };
my $dump_file = $align_dirname."/DOMINO_dump_information.txt";
if (-r -e -s $dump_file) { remove_tree($dump_file); DOMINO::printDump(\%hash, $dump_file);}
DOMINO::print_success_Step("mapping");

###########################
####### SUBROUTINES #######
###########################

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}

sub generate_filter_PILEUP {
	my $sorted_bam = $_[0]; my $contig_file = $_[1]; 
	my $reference_id = $_[2]; my $taxa = $_[3];
	
	my $dir_path = $domino_mapping_files{$reference_id}{'dir'}[0];
	my $tmp = $dir_path."/ARRAY_files_".$taxa."_PROFILE"; 

	my $domino_Scripts_GeneratePileup = $domino_Scripts."/DM_GeneratePileup.pl";
	my $command = "perl $domino_Scripts_GeneratePileup $sorted_bam $contig_file $reference_id $taxa $tmp $path";  
	#print $command."\n"; 
	system($command);	
	chdir $dir_path; #&debugger_print("Changing dir to $dir_path");
	return($tmp);	
}

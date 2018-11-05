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
my $domino_align_files_Ref = DOMINO::get_DOMINO_files($path."/");
my %domino_align_files = %{$domino_align_files_Ref};

#print Dumper $hash_parameters;
#print Dumper $domino_align_files_Ref;

my $align_dirname = $$hash_parameters{'mapping'}{'folder'}[0];
my $num_proc_user = $$hash_parameters{'mapping'}{'cpu'}[0];
my $msa_dirname = $align_dirname."/MSA_files";

my %new_domino_files;
my %new_dowmino_param;
#################################################################

#####################################
### Check MSA: file/folder/RADseq ###
#####################################
mkdir $msa_dirname, 0755; chdir $align_dirname;
push (@{ $domino_align_files{'MSA_files'}{'folder'} }, $msa_dirname); 
#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_align_files);
if ($$hash_parameters{"marker"}{"RADseq"}) {
	my @file_name = split("/", $domino_align_files{'RADseq'}{'file'}[0]);
	my $file_path = $align_dirname."/".$file_name[-1];
	system("ln -s $domino_align_files{'RADseq'}{'file'}[0] $file_path");
	#&debugger_print("RADseq data file provided...");
	print "\n\n"; DOMINO::printHeader(" Checking RADseq file provided ", "%");
	
	## Split files
	print "+ Splitting file into multiple files to speed the computation...\n";
	my $size = DOMINO::get_size($file_path);
	my $chars = int($size/$num_proc_user);
	#	&debugger_print("Total Size: $size\nCharacters to split: $chars");
	#	&debugger_print("File: $file_path");		

	my $files_ref;
	if ($$hash_parameters{'marker'}{'pyRAD'}) { $files_ref = &loci_file_splitter($file_path, $chars,'loci');
	} else { $files_ref = DOMINO::fasta_file_splitter($file_path, $chars, 'fa'); }		
	push (@{ $domino_align_files{'RADseq'}{'parts'}}, @$files_ref);		
	#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_align_files);

	for (my $i=0; $i < scalar @$files_ref; $i++) {
		my $file2dump = $align_dirname."/dump_file_split_Part_".$i.".txt";
		push (@{ $domino_align_files{'all'}{'dump_file_split'} }, $file2dump);
	}

	## Implement threads
	print "\n+ For each loci a MSA file would be generated...\n+ Parsing splitted files...\n+ Using parallel threads ($num_proc_user CPUs)...\n";

	##############	
	## pyRAD
	##############
	if ($$hash_parameters{'marker'}{'pyRAD'}) {		
	#	&debugger_print("pyRAD file provided...");
		### RADSEQ like data ####
		## Parse pyRAD loci file provided
		my $pm_pyRAD =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
		$pm_pyRAD->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code\n"; } );
		$pm_pyRAD->run_on_start( sub { my ($pid,$ident)=@_; print "\t- pyRAD analysis for file $ident and PID=$pid started\n"; } );
		for (my $i=0; $i < scalar @$files_ref; $i++) {
			my @basename = split("/", $$files_ref[$i]);
			my @name = split(".loci", $basename[-1]);
	
			my $pid = $pm_pyRAD->start($name[-1]) and next; 
			my %domino_align_files_pyRAD_split;
			my $counter = 1; my %hash;
			open(FILE, $$files_ref[$i]) || die "Could not open the $$files_ref[$i] ... [DM_MarkerScan: pyRAD]\n";
			while (<FILE>) {		
				next if /^#/ || /^\s*$/;
				my $line = $_;
				chomp $line;
				#&debugger_print($line); &debugger_print("Counter: $counter");
				if ($line =~ /\/\//) { 
					my $file = $msa_dirname."/".$name[0]."_loci_".$counter.".fasta";
					foreach my $keys (sort keys %hash) { 
						if ($domino_align_files{$keys}{'taxa'} || $domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
							open (OUT, ">>$file"); print OUT ">".$keys."\n".$hash{$keys}."\n"; close(OUT);
					}} $counter++; undef %hash; next;
				}
				$line =~ s/\s+/\t/g; 
				$line =~ s/\%/>/g; ## Sometimes there is this symbol or at least in my test set
				my @array = split("\t", $line);
				$array[0] =~ s/\>//;
				$hash{$array[0]} = $array[1];
				if ($domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
					$domino_align_files_pyRAD_split{$array[0]}{'taxa'}[0]++;
				}
			} close (FILE); undef %hash;
			#print Dumper \%domino_align_files_pyRAD_split;
			DOMINO::printDump(\%domino_align_files_pyRAD_split, $domino_align_files{'all'}{'dump_file_split'}[$i]);
			$pm_pyRAD->finish($name[-1]); # pass an exit code to finish
		}
		$pm_pyRAD->wait_all_children; 
		print "\n\n";
		print "***************************************************\n";
		print "**** All pyRAD parsing processes have finished ****\n";
		print "***************************************************\n\n";		
	
	##############	
	## STACKS
	##############	
	} elsif ($$hash_parameters{'marker'}{'STACKS'}) {
		## Parse STACKS file provided
		#&debugger_print("STACKS file provided...");
		my $pm_STACKS =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
		$pm_STACKS->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code \n"; } );
		$pm_STACKS->run_on_start( sub { my ($pid,$ident)=@_; print "\t- STACKS analysis for file $ident and PID=$pid started\n"; } );
		for (my $i=0; $i < scalar @$files_ref; $i++) {
			my @basename = split("/", $$files_ref[$i]);
			my @name = split(".fa", $basename[-1]);
			my $pid = $pm_STACKS->start($name[-1]) and next; 
			
			my %domino_align_files_STACKS_split;
			my (%hash, $new_id); my $first = 0; my $previous;
			$/ = ">"; ## Telling perl where a new line starts
			open(FILE, $$files_ref[$i]) || die "Could not open the $$files_ref[$i] ... [DM_MarkerScan: STACKS]\n";
			while (<FILE>) {		
				next if /^#/ || /^\s*$/; chomp;
				my ($titleline, $sequence) = split(/\n/,$_,2);
				next unless ($sequence && $titleline);
				chomp $sequence;
				$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
				$titleline =~ s/\r//g;
				#&debugger_print($titleline."\t".$sequence);
				if ($titleline =~ /(CLocus\_\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
					my $CLocus = $1; my $sample = $2;
					if ($domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') { $domino_align_files_STACKS_split{$sample}{'taxa'}[0]++; }
					if ($first == 0) { 
						push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence); $first++;  
						$previous = $CLocus; next;
					}							
					if ($hash{$CLocus}) {
						push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence); $previous = $CLocus;
					} else {
						#&debugger_print("Ref", \%hash);
						my %Clocus_hash = %{ $hash{$previous} };
						#print "CLocus hash:\n"; print Dumper \%Clocus_hash;								
						## Parse
						my $out_file = $msa_dirname."/".$previous.".fasta";
						open (OUT, ">>$out_file");
						foreach my $keys (sort keys %Clocus_hash) {
							if ($domino_align_files{$keys}{'taxa'} || $domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
								my @array = @{$Clocus_hash{$keys}};
								for (my $i=0; $i < scalar @array; $i++) {
									my @split = split(":::", $array[$i]);
									print OUT ">".$split[0]."\n".$split[1]."\n";
						}}} close (OUT);
						push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence);
						$previous = $CLocus;
			}}} close(FILE); $/ = "\n";
			
			## Parse last
			my %Clocus_hash = %{ $hash{$previous} };
			my $out_file = $msa_dirname."/".$previous.".fasta";
			open (OUT, ">>$out_file");
			foreach my $keys (sort keys %Clocus_hash) {
				if ($domino_align_files{$keys}{'taxa'} || $domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
					my @array = @{$Clocus_hash{$keys}};
					for (my $i=0; $i < scalar @array; $i++) {
						my @split = split(":::", $array[$i]);
						print OUT ">".$split[0]."\n".$split[1]."\n";
			}}} close (OUT);					
			DOMINO::printDump(\%domino_align_files_STACKS_split, $domino_align_files{'all'}{'dump_file_split'}[$i]);
			$pm_STACKS->finish($name[-1]); # pass an exit code to finish
		}
		$pm_STACKS->wait_all_children; 
		print "\n\n";
		print "****************************************************\n";
		print "**** All STACKS parsing processes have finished ****\n";
		print "****************************************************\n\n";		
	} 
	print "\n"; &time_log(); print "\n";

} else {

	### MSA file or folder provided
	#mkdir $msa_dirname, 0755;
	if ($$hash_parameters{'MSA'}{'file'}) {
		## Check the alignment format provided...
		print "\n\n"; DOMINO::printHeader(" Checking Alignment file provided ", "%");
	#	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_align_files);
		my $species_alignment_ref = &read_phylip_aln($domino_align_files{'MSA'}{'file'}[0]);
		if ($domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
			undef $domino_align_files{'taxa'}{'user_Taxa'};		
			for (my $i=0; $i < scalar @$species_alignment_ref; $i++) {
				push (@{ $domino_align_files{$$species_alignment_ref[$i]}{'taxa'}}, 1);
				push (@{ $domino_align_files{'taxa'}{'user_Taxa'} }, $$species_alignment_ref[$i]);
		}}
		print "\n"; &time_log(); print "\n";		
	} elsif ($$hash_parameters{'MSA_folder'}{'folder'}) {				
		chdir $msa_dirname;
		print "\n\n"; DOMINO::printHeader(" Checking Alignment folder provided ", "%");
		my @array_files = @{ $domino_align_files{'MSA_folder'}{'files'} };
		my $tmp = $align_dirname."/tmp"; mkdir $tmp, 0755;
		for (my $i=0; $i < scalar @array_files; $i++) {
			my $file2dump = $tmp."/dump_file_split_Part_".$i.".txt";
			#print "File to dump: ".$file2dump."\n";
			push (@{ $domino_align_files{'all'}{'dump_file_split'} }, $file2dump);
		}
		#&debugger_print("MSA folder provided...");
		print "+ Using parallel threads ($num_proc_user CPUs)...\n";			
		my $pm_MSA_folder =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
		$pm_MSA_folder->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code \n"; } );
		$pm_MSA_folder->run_on_start( sub { my ($pid,$ident)=@_; print "\t- Alignment analysis for file $ident and PID=$pid started \n"; } );
		for (my $i = 0; $i < scalar @array_files; $i++) {
			my @species_alignment;
			my $pid = $pm_MSA_folder->start($array_files[$i]) and next; 
			my $file_path = $domino_align_files{'MSA_folder'}{'folder'}[0]."/".$array_files[$i];
			my %alignment;			
			my %domino_align_files_MSA_folder;
			if ($array_files[$i] =~ /(.*)\.fasta/) {
				my $name = $1;
				open(FILE, $file_path) || die "Could not open the $file_path... [DM_MarkerScan: MSA folder provided]\n";
				$/ = ">"; ## Telling perl where a new line starts
				while (<FILE>) {		
					next if /^#/ || /^\s*$/;
					chomp;
					my ($titleline, $sequence) = split(/\n/,$_,2);
					next unless ($sequence && $titleline);
					chomp $sequence; chomp $titleline;
					$titleline =~ s/\r//g;						
					if ($domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
						unless (grep /$titleline/, @{ $domino_align_files_MSA_folder{'taxa'}{'user_Taxa'} }) { 
							$domino_align_files_MSA_folder{$titleline}{'taxa'}[0]++;
							push (@{ $domino_align_files_MSA_folder{'taxa'}{'user_Taxa'} }, $titleline);
						}
						$alignment{$titleline} = $sequence;
					} elsif ($domino_align_files{$titleline}{'taxa'}) { $alignment{$titleline} = $sequence; }
				} close(FILE); $/ = "\n";
				
				if ($domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') { 
					system("ln -s $file_path");
				} else {
					my $out_file = $msa_dirname."/parsed_".$name.".fasta";
					open (OUT, ">$out_file");
					foreach my $seqs (sort keys %alignment) {
						print OUT ">".$seqs."\n".$alignment{$seqs}."\n";
					} close (OUT);
			}}

			DOMINO::printDump(\%domino_align_files_MSA_folder, $domino_align_files{'all'}{'dump_file_split'}[$i]);
			$pm_MSA_folder->finish($i); # pass an exit code to finish
		}
		$pm_MSA_folder->wait_all_children; 
		print "\n\n";
		print "*********************************************\n";
		print "**** All parsing processes have finished ****\n";
		print "*********************************************\n\n";
		my $files = scalar @array_files;
		print "\n\n+ Parsing of the $files files has been done...\n"; print "\n"; &time_log(); print "\n";
	}
}

if ($domino_align_files{'all'}{'dump_file_split'}) {
	my @dump_files = @{ $domino_align_files{'all'}{'dump_file_split'} };
	for (my $j=0; $j < scalar @dump_files; $j++) {
		if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
		open (DUMP_IN, "$dump_files[$j]");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			push (@{ $domino_align_files{$array[0]}{$array[1]}}, $array[2]);
} close (DUMP_IN); }}
	
unless ($$hash_parameters{'marker'}{'taxa_string'}) {
	print "\n\n"; DOMINO::printHeader("","#"); print "NOTE:\n\n";
	print "\t+ No taxa names were provided so DOMINO have parsed and checked for the names...\n";	
	print "\t+ DOMINO would use all the taxa available...\n\n";
	my @taxa_sp;
	foreach my $keys (sort keys %domino_align_files) { 
		if ($domino_align_files{$keys}{'taxa'}) {
			DOMINO::printDetails("\tName: $keys\n", $$hash_parameters{'mapping'}{'mapping_parameters'}, $$hash_parameters{'mapping'}{'mapping_markers_errors_details'};
			$number_sp++;
			push (@taxa_sp, $keys);
	}}
	print "\n\n"; DOMINO::printHeader("","#"); print "\n\n";
	if (!$minimum_number_taxa_covered) { 
		$minimum_number_taxa_covered = 0;  ## Force to be any taxa
	}
	push (@{ $new_dowmino_param{'marker'}{'number_taxa'}}, $number_sp);
	push (@{ $new_dowmino_param{'marker'}{'MCT'}}, $minimum_number_taxa_covered);
	$MID_taxa_names = join(",", @taxa_sp);
}
#&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_align_files);	

### dump info
my $dump_file = $align_dirname."/DOMINO_dump_information.txt";
if (-r -e -s $dump_file) { remove_tree($dump_file); DOMINO::printDump(\%new_domino_files, $dump_file);}

### dump parameters
my $dump_param = $align_dirname."/DOMINO_dump_information.txt";
if (-r -e -s $dump_param) { remove_tree($dump_file); DOMINO::printDump(\%new_dowmino_param, $dump_param);}

DOMINO::print_success_Step("parseAlign");


###########################
####### SUBROUTINES #######
###########################

sub loci_file_splitter {
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	my $file = $_[0];
	my $block = $_[1];
	my $ext = $_[2]; # fasta, fastq, loci, fa

	open (FH, "<$file") or die "Could not open file $file [DOMINO.pm:loci_file_splitter]";
	print "\t- Blocks of $block characters would be generated...\n";
	my $j = 0; my @files;
	while (1) {
		my $chunk;
	   	my @tmp = split (".".$ext, $file);
		my $file_name = $tmp[0];
		my $block_file = $file_name."_part-".$j."_tmp.".$ext;
		print "\t- Printing file $block_file\n";
		push (@files, $block_file);
		open(OUT, ">$block_file") or die "Could not open destination file [DOMINO.pm:loci_file_splitter]";
		if (!eof(FH)) { 
			read(FH, $chunk,$block);  
			#if ($j > 0) { $chunk = $chunk; }
			print OUT $chunk;
		} ## Print the amount of chars	
		if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken	
		if (!eof(FH)) { 
			$/ = "\/\/"; ## Telling perl where a new line starts
			$chunk = <FH>; 
			chop $chunk;
			print OUT $chunk;
			print OUT "/\n"; 
			$/ = "\n";
			$chunk = <FH>; 
		} ## print the sequence if it is broken
		$j++; close(OUT); last if eof(FH);
	}
	close(FH); return (\@files);	
}

sub read_phylip_aln {
	
	######################################################################
	##  This function process a MSA into a hash of sequences and it 	##
	## 	calls another function to generate the variation profiles.		##
	##	Jose F. Sanchez 02/2/2016										##
	######################################################################
	
	my $file = $_[0];
	my $region_provided = $_[1];
	
	open (INPUT, $file);
	my ($length, $species, %aln, $counter, $region);
	while (<INPUT>) {
		next if /^#/;
		next if /^\s+$/;
		my $line = $_;
		chomp $line;
		if ($line =~ /\s*(\d+)\s+(\d+)/) {
			$species = $1;
			$length = $2;
			$counter++;
			$region = "msa_".$counter;	
			next;
		}			
		if ($line =~ /(\S+)\s+(.*)/) {			
			$aln{$region}{0}{"name"} = $1;
			my $seq = $2;
			$seq =~ s/\s//g;
			$aln{$region}{0}{"seq"} = $seq;
			my $left_species = $species - 1;
			for (my $i=0; $i < $left_species; $i++) {
				my $next_line= <INPUT>;
				chomp $next_line;
				if ($next_line =~ /(\S+)\s+(.*)/) {
					my $id = $i+1;
					$aln{$region}{$id}{"name"} = $1;
					my $seq = $2;
					$seq =~ s/\s//g;
					$aln{$region}{$id}{"seq"} = $seq;
		}}} else {
			my @array = split("\t", $line);
			my $sequence = $aln{$region}{0}{"seq"};
			$sequence .= $array[0];
			$sequence =~ s/\n//;
			$aln{$region}{0}{"seq"} = $sequence;
			my $left_species = $species - 1;
			for (my $i=0; $i < $left_species; $i++) {
				my $next_line= <INPUT>;
				chomp $next_line;
				my $id = $i+1;
				my @array = split("\t", $next_line);
				my $sequence = $aln{$region}{$id}{"seq"};
				$sequence .= $array[0];
				$sequence =~ s/\n//;
				$aln{$region}{$id}{"seq"} = $sequence;
	}}}
	close(INPUT);
	my @array_taxa;
	foreach my $keys (sort keys %aln) {
		my %alignment;
		if (!$region_provided) { $region_provided = $keys; }
		my $region_Fasta = $msa_dirname."/".$region_provided.".fasta";
		open (OUT_MSA, ">$region_Fasta");
		foreach my $taxa (sort keys %{ $aln{$keys} }) {
			if ($domino_align_files{$aln{$keys}{$taxa}{"name"}}{'taxa'} || $domino_align_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
				print OUT_MSA ">".$aln{$keys}{$taxa}{"name"}."\n".$aln{$keys}{$taxa}{"seq"}."\n";
			}
			unless (grep /$aln{$keys}{$taxa}{"name"}/, @array_taxa ) {
				push (@array_taxa, $aln{$keys}{$taxa}{"name"});
		}}
		close(OUT_MSA);	
		$region_provided = "";
	}
	return \@array_taxa;
}

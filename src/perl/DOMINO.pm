#!/usr/bin/perl
####################################################################
###	DOMINO: Development of molecular markers in non-model organisms 
####################################################################
# This package provides multiple subroutines for the DOMINO package
## date 24/10/18

package DOMINO;

use FindBin;
use lib $FindBin::Bin."/lib";
require List::MoreUtils;
use List::MoreUtils qw(firstidx);
use Data::Dumper;

#################################################################
my $domino_Scripts = $FindBin::Bin;
my $scripts_path = $FindBin::Bin."/../";
## General binaries variables
my $samtools_path = $scripts_path."samtools-1.3.1/samtools";
my $bowtie_path = $scripts_path."bowtie2-2.2.9/";
my $BLAST = $scripts_path."NCBI_BLAST/";
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
#################################################################

###################
## Call programs ##
###################

sub blastn {
	my $file = $_[0]; my $db = $_[1]; my $results = $_[2]; my $BLAST = $_[3]; 
	my $filter = $BLAST."blastn -query ".$file." -evalue 1e-10 -max_target_seqs 150 -db '".$db."' -out $results -outfmt 6";
	my $message = "BLASTN command: $filter\n"; 
	my $blastn = system($filter);
	return ($blastn, $message);
}

sub makeblastdb {
	my $file = $_[0]; my $BLAST_path = $_[1]; my $error_log = $_[2];
	my $db;
	my $make_blast_db = $BLAST_path."makeblastdb -in ".$file;
	$make_blast_db .= " -dbtype nucl";
	if ($file =~ /(.*)\.fasta/ ) { $db = $1 ; } 
	$make_blast_db .= " -out ".$db." -logfile ".$db.".log 2> $error_log";	
	my $makeblastresult = system($make_blast_db);
	if ($makeblastresult != 0) { print "Generating the database failed when trying to proccess the file... DOMINO would not stop in this step...\n"; return "1"; }
	return ($db, $make_blast_db);
}

sub mothur_remove_seqs {	
	# This subroutine takes as input a FASTQ file AND classifies according to the tags provided, by Roche and trims the seqs
	my $fasta = $_[0]; my $qual = $_[1]; my $directory = $_[2];
	my $ids2remove = $_[3]; my $mothur_path = $_[4];
	my $line;
	if ($qual eq 'YES') { 
		$line = $mothur_path." '#set.dir(output=$directory); remove.seqs(accnos=$ids2remove, fastq=$fasta)'";
	} elsif ($qual eq 'NO') {
		$line = $mothur_path." '#set.dir(output=$directory); remove.seqs(accnos=$ids2remove, fasta=$fasta)'";
	} else { 
		$line = $mothur_path." '#set.dir(output=$directory); remove.seqs(accnos=$ids2remove, fasta=$fasta, qfile=$qual)'";
	}
	print "\n+ Calling mothur executable for discarding reads...\n\n";
	my $system_call = system($line);
}

sub mothur_retrieve_seqs {
	## This sub retrieve a given number of ids and generates a fastq
	my $fasta = $_[0]; my $qual = $_[1]; my $directory = $_[2]; 
	my $ids2retrieve = $_[3]; my $mothur_path = $_[4];
	
	my $line = $mothur_path." '#set.dir(output=$directory); get.seqs(accnos=$ids2retrieve, fasta=$fasta, qfile=$qual)'";
	print "\n+ Calling mothur executable for retrieving reads in file $ids2retrieve...\n\n";
	my $system_call = system($line);
}

sub mothur_retrieve_FASTA_seqs {
	## This sub retrieve a given number of ids and generates a fastq
	my $fasta = $_[0]; my $directory = $_[1]; 
	my $ids2retrieve = $_[2]; my $mothur_path = $_[3];
	
	my $line = $mothur_path." '#set.dir(output=$directory); get.seqs(accnos=$ids2retrieve, fasta=$fasta)'";
	print "\n+ Calling mothur executable for retrieving sequences in file $ids2retrieve...\n\n";
	my $system_call = system($line);
}

sub generate_bam {
	my $sam_file = $_[0];
	my $num_proc_user = $_[1];
	my $avoid = $_[2];
	
	my @temp = split ("\.sam", $sam_file);
	my $name = $temp[0]; my $bam_file = $name.".bam";
	#&debugger_print("- Generating a BAM file for $sam_file\n"); 	
	my $system_samtools_sam2bam = "$samtools_path view --threads $num_proc_user -bS -o $bam_file $sam_file";
	#print $system_samtools_sam2bam."\n"; #&debugger_print("SAMTOOLS command: $system_samtools_sam2bam");	
	my $system_call = system ($system_samtools_sam2bam);
	if ($system_call != 0) {
		if (!$avoid) { &printError("Some error happened when calling SAMTOOLs for SAM -> BAM conversion $sam_file to $bam_file...."); &dieNicely(); }
	}
	my $sorted;
	if ($avoid) { if ($avoid eq 'sam') { $sorted = &generate_sorted_bam($bam_file, $num_proc_user, "sam"); } }
	if (!$avoid) {$sorted = &generate_sorted_bam($bam_file, $num_proc_user);}
	return $sorted;
}

sub generate_sorted_bam {
	my $bam_file = $_[0];
	my $num_proc_user = $_[1];
	my $sam = $_[2];

	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	#&debugger_print("- Sorting the BAM file: $bam_file\n"); 	
	my $sorted;	
	my $system_samtools_sort;
	if ($sam) {
		$sorted = $name.".sorted.sam";	
		$system_samtools_sort = "$samtools_path sort --threads $num_proc_user --output-fmt SAM -o $sorted $bam_file";
	} else {
		$sorted = $name.".sorted.bam";	
		$system_samtools_sort = "$samtools_path sort --threads $num_proc_user -o $sorted $bam_file";
	}	 
	#&debugger_print("SAMTOOLS command: ".$system_samtools_sort);
	my $system_call_2 = system ($system_samtools_sort);
	if ($system_call_2 != 0) {
		&printError("Some error happened when calling SAMTOOLs for sorting BAM file...."); &dieNicely();
	}
	return $sorted;
}

sub generate_sam {
	my $bam_file = $_[0];
	my $num_proc_user = $_[1];

	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	#&debugger_print("- Generating a SAM file for $bam_file\n"); 	
	my $system_samtools_bam2sam = $samtools_path." view --threads $num_proc_user -h ".$bam_file." -o ".$name.".sam";
	#&debugger_print("SAMTOOLS command: $system_samtools_bam2sam");
	my $system_call = system ($system_samtools_bam2sam);
	if ($system_call != 0) {
		&printError("Some error happened when calling SAMTOOLs for BAM -> SAM conversion....");  &dieNicely();
	}
	return $name.".sam";
}

sub generate_index_bam {
	my $bam_file = $_[0];
	#print "\t- Generating an index bam file for $bam_file\n"; 	
	my $index_system = $samtools_path." index ".$bam_file;
	#&debugger_print("SAMTOOLS command: ".$index_system);
	my $index_call = system($index_system);
	if ($index_call != 0) {
		&printError("Some error happened when calling SAMTOOLs for indexing the BAM file [$bam_file]...."); &dieNicely();
	}
}

##################
## Call scripts ##
##################

sub Contig_Stats { 
	
	my $fasta_file = $_[0];
	my @name = split("\.fasta", $fasta_file);
	my $outFile = $name[0]."-statistics.csv";
	my $domino_Scripts_contig = $domino_Scripts."/DM_ContigStats.pl";
	system("perl $domino_Scripts_contig $fasta_file $outFile");
	return $outFile;
}


####################
## Print messages ##
####################

sub printFormat_message {
	print "\n\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\nWhere:\n\txxx: any character or none.Please avoid using dots (.)\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n\n";
}

sub printHeader {
	my $sentence = $_[0]; my $symbol = $_[1];	
	my @length_array = ($symbol) x 97;	
	my @array_sentence = split("",$sentence);
	my $length = scalar @array_sentence;	
	my $start = 49 - ($length/2);	
	for (my $i = 0; $i < scalar @array_sentence; $i++) {
		$length_array[$start+$i] = $array_sentence[$i];
	}	
	my $string = join("", @length_array);
	print $string."\n";
}

sub printDetails {
	my $string = $_[0]; 
	my $param_Detail_file = $_[1];
	my $file2 = $_[2];
	
	my @files;
	if ($file2) {
		@files = ($param_Detail_file, $file2);
	} else {
		push (@files, $param_Detail_file);
	}
	for (my $i=0; $i < scalar @files; $i++) {
		open (PARAM, ">>$files[$i]");
		print PARAM $string;
		close(PARAM);
	}
	print $string;		
}

sub printDump {
	my $ref = $_[0]; 
	my $out_file = $_[1];
	my $simpleHash = $_[2];

	my %tmp_hash = %$ref;
	
	open (DUMP, ">>$out_file"); 
	#&debugger_print("Printing to a file $out_file");
	foreach my $names (sort keys %tmp_hash) {
		if ($simpleHash) { print DUMP $names."\t".$tmp_hash{$names}."\n"; next;}
		my %hash = %{$tmp_hash{$names}};
		foreach my $tags (sort keys %hash) {
			my @array = @{ $hash{$tags} };
			for (my $i=0; $i < scalar @array; $i++) {
				print DUMP $names."\t".$tags."\t".$hash{$tags}[$i]."\n";
	}}} close (DUMP);
}

sub printError {
    my $msg = $_[0];
	print "\n\n";&printHeader(" ERROR ","!!"); print "\n";
    print $msg."\n\nTry \'perl $0 -h|--help or -man\' for more information.\nExit program.\n";
	print "\n\n"; &printHeader("","!!"); &printHeader("","!!"); 
    &printError_log($msg, $error_log);
}

sub printError_log {
	my $message = $_[0]; my $error_log = $_[1];
	open (ERR, ">>$error_log");
	print ERR $message."\n";
	close (ERR);
}

sub printInput_type {

print "\n\nType of input explanation:\nAccording to the type of NGS files provided several options are available but only one option would be provided: -type_input [int]
 1: A single file in Standard Flowgram Format (SFF), 454 file, containing all the reads of the different taxa accordingly tagged
 2: FASTQ files coming from 454 Roche. A single file containing all the reads of the different taxa accordingly tagged
 3: Multiple FASTQ files coming from 454 Roche. Each file contains each taxa reads. 
 4: FASTQ file from Illumina single end, containing all the reads of the different taxa accordingly tagged
 5: Multiple FASTQ file from Illumina single end. Each file contains each taxa reads. 
 6: A single pair of FASTQ files from Illumina paired-end sequencing: Each pair would contain all the reads of the different taxa accordingly tagged. 
 Please tagged left read file with xxx_R1.fastq and right read file as xxx_R2.fastq
 7: Multiple FASTQ files from Illumina paired-end: Each pair of files containing each taxa left and right reads respectively\n\n";
 
 DOMINO::printFormat_message();
 
}

sub finish_time_stamp {

	my $start_time = $_[0];
	my $finish_time = time;
	print "\n\n"; DOMINO::printHeader("","+"); 
	DOMINO::printHeader(" ANALYSIS FINISHED ","+"); 
	DOMINO::printHeader("","+"); 
	print DOMINO::time_stamp();
	my $secs = $finish_time - $start_time; 
	my $hours = int($secs/3600); 
	$secs %= 3600; 	
	my $mins = int($secs/60); 
	$secs %= 60; 
	printf (" Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub debugger_print {
	my $string = $_[0];
	my $ref = $_[1];
	## Print to a file
	if ($string eq "Ref") {
		print "\n******\nDEBUG:\n";
		print Dumper $ref; ## if filehandle OUT: print OUT Dumper $ref;
		print "******\n";
	} else {
		print "\n******\nDEBUG:\n".$string."\n******\n\n";
	}
}

sub dieNicely {
	print DOMINO::time_stamp();	print "\n\nTry perl $0 -man for more information\n\n";
	Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 0 );
}

sub time_stamp { return "[ ".(localtime)." ]"; }

sub time_log {	
	my $given_step_time = $_[0];
	my $current_time = time;

	print DOMINO::time_stamp."\t";
	my $secs = $current_time - $given_step_time; 
	my $hours = int($secs/3600); $secs %= 3600; 
	my $mins = int($secs/60); $secs %= 60; 
	printf ("Step took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
	return \$current_time;
}

sub print_success_Step {
	my $name = $_[0];
	my $new = $name.".success";
	open (OUT, ">$new");
	print OUT "OK\n";
	close OUT;
}

sub print_fail_Step {
	my $name = $_[0];
	my $new = $name.".failed";
	open (OUT, ">$new");
	print OUT "FAIL\n";
	close OUT;
}

####################
## Maths routines ##
####################

sub sum {
	my $array_ref = $_[0];
	my $sum=0;
	for (my $i=0; $i < scalar @{ $array_ref } ; $i++) {
		$sum += $$array_ref[$i];
	}
	return $sum;
}

sub min_max {
	my $array_ref = $_[0];
	my $min=9999999999999999999999999999; my $max=0;
	for (my $i=0; $i < scalar @{ $array_ref }; $i++) {
		if ($$array_ref[$i] < $min) {$min = $$array_ref[$i]; }
		if ($$array_ref[$i] > $max) {$max = $$array_ref[$i]; }
	}
	return ($min, $max)
}

sub Poisson_distribution {
	my ($x, $a) = @_;
	return unless $a >= 0 && $x >= 0 && $x == int($x); 
	return (($a ** $x) * exp(-$a))/&functional_factorial($x);	
}

sub binomial {
	## Returns the probability of an event ocurring $k times in $n attempts where the probability
	## of it occurring in a sample attempt is $p
    my $k = shift(@_); ## Occurrences
    my $n = shift(@_); ## Attempts
    my $p = shift(@_); ## Probability
    my $prob = ($p**$k) * ((1 - $p)**($n - $k)) * &functional_factorial($n) / (&functional_factorial($k) * &functional_factorial($n - $k));
    return $prob;
}

sub functional_factorial {

    my $n = shift(@_);
    my $fact = 1;
    if (($n < 0) or (170 < $n)) { die "Factorial out of range [DM_MarkerScan: functional_factorial]"; }
    for (my $i = 1; $i <= $n; $i++) { $fact *= $i; }
    return $fact;
}

sub convert_ASCII_to_number {

	my $qual_ASCII_string =	$_[0];
	my @qual_array_Ascii = split("", $qual_ASCII_string);
	my $string_qual_nums;
	for (my $q = 0; $q < scalar @qual_array_Ascii; $q++) {
		my $tmp = ord($qual_array_Ascii[$q]) - 33; ## ord: basic perl function: returns number given ASCII
		$string_qual_nums .= $tmp." ";
	}
	return $string_qual_nums;
}

####################
## Get Parameters ##
####################

sub retrieve_info {	
	my $dump_files_ref = $_[0];
	my $hash_Ref = $_[1];
	my %hash;
	unless ($hash_Ref eq 1) { %hash = %{$hash_Ref}; }
	my @dump_files = @{ $dump_files_ref };
	for (my $j=0; $j < scalar @dump_files; $j++) {
		if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
		unless (-e -r -s $dump_files[$j]) { next; }
		open (DUMP_IN, "$dump_files[$j]");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			if ($hash_Ref eq 1) { $hash{$array[0]} = $array[1];
			} else { 
				push (@{ $hash{$array[0]}{$array[1]}}, $array[2]);
		}} close (DUMP_IN);
	} 
	if ($hash_Ref eq 1) { return \%hash; }
	my $hash_ref = &get_uniq_hash(\%hash);
	return $hash_ref;
}

sub get_earliest {
	
	my $option = $_[0];
	my $folder = $_[1];
	#$option == "clean_data, assembly, mapping" 
	
	my $array_files_ref=DOMINO::readDir($folder);
	my @array_files = @$array_files_ref;
	my (%mapping_dirs, $earliest);
	for (my $i=0; $i<scalar @array_files;$i++) {
		if ($array_files[$i] eq "." || $array_files[$i] eq ".." || $array_files[$i] eq ".DS_Store") {next;}
		if ($array_files[$i] =~ /(\d+)\_DM\_$option$/) {
			my $time_log=$1;
			if (!$earliest) { $earliest = $time_log; 
			} else { 
				if ($time_log > $earliest) { 
					$earliest = $time_log
			}}
			$mapping_dirs{$time_log} = $folder."/".$array_files[$i];
	}}
	if (!exists $mapping_dirs{$earliest}) { return 'NO';
	} else { return $mapping_dirs{$earliest}; }	
}

sub readDir {
	my $dir = $_[0];
	opendir(DIR, $dir) or die "ERROR: Can not open folder $dir..."; ## FIX ADD TO ERROR-LOG
	my @dir_files = readdir(DIR);

	## Discard '.', '..' and '.DS_Store'
	my $idx_dot = firstidx { $_ eq '.' } @dir_files;
	if ($idx_dot >= 0) { splice(@dir_files, $idx_dot, 1); }
	my $idx_dot_dot = firstidx { $_ eq '..' } @dir_files;
	if ($idx_dot_dot >= 0) { splice(@dir_files, $idx_dot_dot, 1); }
	my $idx_dot_DSStore = firstidx { $_ eq '.DS_Store' } @dir_files;
	if ($idx_dot_DSStore >= 0) { splice(@dir_files, $idx_dot_DSStore, 1); }	
	my $array_ref = \@dir_files;
	return $array_ref;
}

sub get_parameters {

	my $path = $_[0];
	my %parameters;
	my $array_files_ref=DOMINO::readDir($path);
	my @array_files = @$array_files_ref;
	my (%dirs, $earliest);
	for (my $i=0; $i<scalar @array_files;$i++) {
		if ($array_files[$i] eq "." || $array_files[$i] eq ".." || $array_files[$i] eq ".DS_Store") {next;}
		next if ($array_files[$i] =~ /.*old.*/);
		if ($array_files[$i] =~ /(\d+)\_DM\_(.*)/) {
			unless (-d $path."/".$array_files[$i]) {next; };
			my $time_stamp=$1; my $type = $2;
			$dirs{$type}{$time_stamp} = $path.$array_files[$i];
	}}
	#print Dumper \%dirs;  	
	my @params; my @info; my %initial_files;
	foreach my $dir (sort keys %dirs) {
		my $last;
		foreach my $times (sort {$a<=>$b} keys %{$dirs{$dir}}) {
			$last = $times;	## Only used the last folder for each process			
		}
		if ($dir eq "assembly" || $dir eq "clean_data" || $dir eq "mapping" || $dir eq "markers" ) {
			my $params = $dirs{$dir}{$last}."/DOMINO_dump_param.txt"; 		push (@params, $params);
	}}
	#print Dumper \@params;
	my $hash_param = &retrieve_info(\@params, \%parameters);
	return $hash_param;
}

sub get_DOMINO_files {
	my $path = $_[0];
	my $array_files_ref=DOMINO::readDir($path);
	my @array_files = @$array_files_ref;
	my (%dirs, $earliest);
	for (my $i=0; $i<scalar @array_files;$i++) {
		if ($array_files[$i] eq "." || $array_files[$i] eq ".." || $array_files[$i] eq ".DS_Store") {next;}
		next if ($array_files[$i] =~ /.*old.*/);
		if ($array_files[$i] =~ /(\d+)\_DM\_(.*)/) {
			unless (-d $path.$array_files[$i]) {next;};
			my $time_stamp=$1; my $type = $2;
			$dirs{$type}{$time_stamp} = $path.$array_files[$i];
	}}
	my @params; my @info; my %initial_files;
	foreach my $dir (sort keys %dirs) {
		my $last;
		foreach my $times (sort {$a<=>$b} keys %{$dirs{$dir}}) {
			$last = $times;	## Only used the last folder for each process			
		}
		
		if ($dir eq "assembly" || $dir eq "clean_data" || $dir eq "mapping" || $dir eq "markers" ) {
			my $files = $dirs{$dir}{$last}."/DOMINO_dump_information.txt";	push (@info, $files);
	}}
	my $hash_info = &retrieve_info(\@info, \%initial_files);
	return $hash_info;
}

sub get_uniq_hash {
	my $hash_Ref = $_[0];
	my %new_hash;
	foreach my $keys (keys %{ $hash_Ref }) {
		foreach my $subkeys (keys %{ $$hash_Ref{$keys} }) {
			my @array = @{ $$hash_Ref{$keys}{$subkeys} }; 
			my @uniq_sort = do { my %seen; grep { !$seen{$_}++ } @array};
			push (@{ $new_hash{$keys}{$subkeys} }, @uniq_sort);
	}}
	return \%new_hash;
}

################
## Some Utils ##
################

sub check_ID_length {
	
	my $file = $_[0];
	my $option = $_[1];
	open (F1, "$file") or die "Could not open file $file for reading ID length [DOMINO.pm: check_id_length]";
	
	if ($option eq 'fastq') {
		my $isEOF = 1;
		open (F1, "$file");
		while(<F1>) {
			my @Read = ();
			chomp(my $id = $_);
			last if($id=~ /^\n$/);
			## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
			for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>; }
			chomp(my $QualLine = $Read[2]); chomp(my $SeqLine = $Read[0]);
			my $length_id = length($id);
			if ($length_id < 40) { return (1, "Length < 40 char");
			} else { return (0, "Length > 40 char. Lets rename the files");
		}} close (F1);
	} else {	
		$/ = ">"; ## Telling perl where a new line starts
		while (<F1>) {		
			next if /^#/ || /^\s*$/;
			chomp;
	    	my ($titleline, $sequence) = split(/\n/,$_,2);
	    	next unless ($sequence && $titleline);
			my $length_id = length($titleline);
	    	if ($length_id < 40) { return (1, "Length < 40 char");
	    	} else { return (0, "Length > 40 char. Lets rename the files");
		}} close(F1); $/ = "\n";
}}

sub check_file_format {    
    my $file = $_[0]; my $id; my $count = 3;
    my $fasta = my $fastq = my $qual = 0; my $format = 'unknown';
    open(FILE, $file) or die "Could not open file $file when checking file format... [DOMINO.pm: check_file_format]";
    while (<FILE>) {
        if($count-- == 0) { last;
        } elsif(!$fasta && /^\>\S+\s*/o) { $fasta = 1; $qual = 1;
        } elsif($fasta == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) { $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) { $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) { $id = $1; $fastq = 1;
        } elsif($fastq == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) { $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/o) { $fastq = 3 if($id eq $1 || /^\+\s*$/o);
    }}    
    if($fasta == 2) { $format = 'fasta';
    } elsif($qual == 2) { $format = 'qual';
    } elsif($fastq == 3) { $format = 'fastq'; }
    return $format;
}

sub check_paired_file {

	my $file = $_[0];
	my $option = $_[1];
	my ($tmp_id, $pair, $pair_return, @pair, @types, $type);
	my $count = 0;
	
	if ($option eq "fastq") {
		open (F1, "$file") or die "Could not open file $file for checking paired-end reads [DOMINO.pm: check_paired_file]";
		while(<F1>) {
			my @Read = ();
			chomp(my $id = $_);
			last if($id=~ /^\n$/);
			## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
			for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>; }
			chomp(my $QualLine = $Read[2]);
			chomp(my $SeqLine = $Read[0]);
			if ($id =~ /(.*)\/(\d+)/) { ## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/2
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "/1,/2");
			} elsif ($id =~ /(\S*)\s+(\d+)\:N\:\.*/) {			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918 2:N:0:CCGTCC
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "1:,2:");
			} elsif ($id =~ /(.*)\/(\R|\L)/) { 			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/L
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "L,R");
			}
			if ($count == 5) {last;}
		} close(F1);
		
	} else {
		open (FILE, "$file") or die "Could not open file $file when checking paired-end reads [DOMINO.pm: check_paired_file]";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($id, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $id);
			if ($id =~ /(.*)\/(\d+)/) {
				## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/2
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "/1,/2");
			} elsif ($id =~ /(\S*)\s+(\d+)\:N\:\.*/) {
				## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918 2:N:0:CCGTCC
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "1:,2:");
			} elsif ($id =~ /(.*)\/(\R|\L)/) {
				## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/L
				$tmp_id = $1; $pair = $2; $count++;
				push(@pair, $pair); push(@types, "L,R");
			}
			if ($count == 5) { last; }
		}
		close(FILE);
		$/ = "\n";
	}
	
	for (my $i = 1; $i < scalar @pair; $i++) {
		if ($pair[0] == $pair[$i]) {
			$pair_return = $pair[$i];
		} else { die "The reads provided within the file $file do not belong to the same pair...[DOMINO.pm: check_paired_file]";}
		if ($types[0] eq $types[$i]) {
			$type = $types[$i];
		} else { die "The reads provided within the file $file do not belong to the same pair... [DOMINO.pm: check_paired_file]";}
	}
	return ($pair_return, $type);
}

sub fasta_file_splitter {
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	my $file = $_[0];
	my $block = $_[1];
	my $ext = $_[2]; # fasta, fastq, loci, fa
	my $dir = $_[3];

	open (FH, "<$file") or die "Could not open file $file [DOMINO.pm:fasta_file_splitter]";
	print "\t- Splitting file into blocks of $block characters aprox ...\n";
	my $j = 0; my @files;
	while (1) {
		my $chunk;
	   	my @tmp = split ("\.".$ext, $file);
		my $block_file = $tmp[0]."_part-".$j."_tmp.".$ext;
		push (@files, $block_file);
		open(OUT, ">$block_file") or die "Could not open destination file: $block_file [DOMINO.pm:fasta_file_splitter]";
		if (!eof(FH)) { read(FH, $chunk,$block);  
			if ($j > 0) { $chunk = ">".$chunk; }
			print OUT $chunk;
		} ## Print the amount of chars	
		if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken	
		if (!eof(FH)) { 
			$/ = ">"; ## Telling perl where a new line starts
			$chunk = <FH>; chop $chunk; print OUT $chunk; 
			$/ = "\n";
		} ## print the sequence if it is broken
		$j++; close(OUT); last if eof(FH);
	}
	close(FH);
	return (\@files);
}

sub file_splitter {
	
	my $file = $_[0];
	my $block = $_[1];
	my $ext = $_[2]; # fasta, fastq, loci, fa
	
	my @files;
	
	# Splits a file such a sam or whatever file that could be read for each line
	open (FH, "<$file") or die "Could not open file $file [DOMINO.pm:file_splitter]";
	print "+ Splitting file $file into blocks of $block characters...\n";
	my $j = 0; 
	while (1) {
    	my $chunk;
    	my @tmp = split (".".$ext, $file);
		my $file_name = $tmp[0];
		
	   	my $block_file = $file_name."_part-".$j."_tmp.".$ext;
		print "\t- Printing file: ".$block_file."\n";
    	push (@files, $block_file);
    	open(OUT, ">$block_file") or die "Could not open destination file [DOMINO.pm:file_splitter]";
    	$j++;
    	if (!eof(FH)) { read(FH, $chunk,$block);  print OUT $chunk; } ## Print the amount of chars
    	if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken
    	close(OUT); last if eof(FH);
	}
	close(FH);
	return (\@files);	
}

sub get_seq_sizes {
	
	my $file = $_[0];
	
	my $out = $file."_sizes";
	my $hash_size = &readFASTA_hashLength(\$file);
	open (OUT, ">$out");
	foreach my $keys (keys %{$hash_size}) {
		print OUT $$hash_size{$keys}."\t".$keys."\n";
	}
	close (OUT);
	return $out;	
}

sub get_size { # Multiplatform	
	my $size = -s $_[0]; #To get only size
	return $size;	
}

sub line_splitter {
	# given a separator, splits a given line
	my $separator = $_[0];
	my @array = split($separator, $_[1]); 
	return \@array; 
}

sub qualfa2fq_modified_bwa {

	##########################################################################################
	##  This is a modified version of the perl script provided in the bwa-0.7.8 package		##	
	## 	to convert fasta+qual files into fastq files										##
	##	Jose Fco. Sanchez Herrero, 06/05/2014 jfsanchezherrero@ub.edu						##
	##########################################################################################
	
	my ($fhs, $fhq, $q, $a);
	# fasta = $_[0], qual = $_[1]
	my $taxa_name = $_[2];
	my $name = $_[3];
	my $file_type = $_[4];
	
	my @temp_name = split ("fasta", $_[0]);
	my $fastq_file; my $bool_change_name = 0; my $pair;

	my $codeReturn = DOMINO::check_ID_length($_[0], "fasta");
	my ($pairReturn, $type);
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		($pairReturn, $type) = DOMINO::check_paired_file($_[0], "fasta");
		#print "File: $_[0]\tPair: ".$pairReturn."\t"."Type: ".$type."\n";
		if ($pairReturn eq "L" ) { 
			$codeReturn = 0; ## Force to rewrite seq					
			#&debugger_print("Lets change id $_[0] : pair was $pairReturn");
		} elsif ($pairReturn eq "R") { 
			$codeReturn = 0; ## Force to rewrite seq
			#&debugger_print("Lets change id $_[0] : pair was $pairReturn");
		} # else {
			## There is no need to change file according to pair
	}
	if ($codeReturn == 0) { $bool_change_name = 1; }	
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		($pair, $type) = DOMINO::check_paired_file($_[0], "fasta"); #&debugger_print("Pair: ".$pair."\t"."Type: ".$type);
	}
	if ($name) { $fastq_file = $name; } else { $fastq_file = $temp_name[0]."fastq";}

	## Writing
	open (FASTQ_out, ">$fastq_file");	
	open($fhs, ($_[0])) or die ("Unable to open $_[0] file"); #fasta_file
	open($fhq, ($_[1])) or die ("Unable to open $_[1] file"); #qual_file
	my $int=0;
	print "+ Generating a fastq file $fastq_file now...\n";
	$/ = ">"; <$fhs>; <$fhq>; $/ = "\n";
	while (<$fhs>) {
		my $id = <$fhq>;
  		$/ = ">";
  		$a = <$fhs>; 
  		$q = <$fhq>;
  		chomp($q); chomp($a);
		my @array_seq = split ("\n", $a);
		my $seq;
		for (my $i = 0; $i < scalar @array_seq; $i++) { $seq .= $array_seq[$i]; }
		my @array_qual = split (" ", $q);
		my $qual;
		for (my $j = 0; $j < scalar @array_qual; $j++) {
			if ($array_qual[$j] =~ /\d+/) {
				if ($array_qual[$j] > 80) { $array_qual[$j] = $array_qual[$j] - 30; }
				my $qual_tmp = $array_qual[$j];
				$qual_tmp =~ s/(\d+)/chr($1+33)/eg;	
				$qual .=$qual_tmp;
		}}
		if ($bool_change_name == 1) {
			$int++; my $id_new = "Seq_".$int."_".$taxa_name;
			my $new_pair;
			if ($pair) { 
				if ($pair eq "L" ) { $new_pair = 1;					
				} elsif ($pair eq "R") { $new_pair = 2;					
				} else { $new_pair = $pair; }
				# /1|/2 or 1:N|2:N
				$id_new .= "/".$new_pair;
			}
			if (length($id_new) > 40) {
				$id_new = "Seq_".$int;
				if ($pair) { $id_new .= "/".$new_pair; }
			}
			print FASTQ_out "\@".$id_new."\n".$seq."\n+\n".$qual."\n";
		} else { 
			print FASTQ_out "\@".$id.$seq."\n+\n".$qual."\n";
		}
		$/ = "\n";
	} close($fhs); close($fhq); close(FASTQ_out);
	$/ = "\n"; ## Telling Perl where a new line starts
	return $fastq_file;
}

sub readFASTA_hash {

	my $file = $_[0];
	my %hash;
	open(FILE, $file) || die "Could not open the file $file [DOMINO.pm: readFASTA_hash]\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	chop $sequence;
    	$hash{$titleline} = $sequence;
	}
	close(FILE); $/ = "\n";
	my $hashRef = \%hash;
	return $hashRef;
}

sub readFASTA_hashLength {

	my $file = $_[0];
	my %hash;
	my $counter;
	my $message;
	open(FILE, $file) || die "Could not open the file $file [DOMINO.pm: readFASTA_hashLength] \n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	$titleline =~ s/ /\t/g;
    	my @array_titleline = split("\t",$titleline);    	
		chomp $sequence;
		$sequence =~ s/\s+//g;
		$sequence =~ s/\r//g;
		$titleline =~ s/\r//g;
   		my $size = length($sequence);
	   	$hash{$array_titleline[0]} = $size;
	   	$counter++;
	   	
   		my $length;
   		my $same_counter = 0;
	   	if ($counter == 10) {
	   		my @ids = keys %hash;
			for (my $j=0; $j < scalar @ids; $j++) {
				if ($hash{$ids[0]} == $hash{$ids[$j]}) {
					$length = $hash{$ids[$j]}; $same_counter++;
		}}}
		## There is no need to continue if all the reads are the same
		if ($same_counter > 5) { 
			$message = "Fixed size for reads in $file -- $length";
			my %return = ("UNDEF" => $length); 
			return (\%return, $message); 
		} 
	}
	close(FILE); $/ = "\n"; 
	$message = "Different sequence lengths";
	return (\%hash, $message);
}

sub readFASTQ_IDSfile {
	
	my $file = $_[0];
	my %hash;		
	open (F1, "$file") or &printError("Could not open file $file [DOMINO.pm:readFASTA_IDSfile]") and exit();
	while(<F1>) {
		my @Read = ();
		chomp(my $id = $_);
		last if($id=~ /^\n$/);

		## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
		for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>; }
		chomp(my $QualLine = $Read[2]);
		chomp(my $SeqLine = $Read[0]);
		my $pair; my $seq_id;
		if ($id =~ /(.*)(\/\d+)/) { $seq_id = $1; $pair = $2; }
		$hash{$seq_id}++;
	}
	close(F1);	
	my $hashRefreturn = \%hash;
	return $hashRefreturn;
}

sub readFASTA_IDSfile {

	my $file = $_[0];
	my %hash;
	open(FILE, $file) || die "Could not open the file $file [DOMINO.pm:readFASTA_IDSfile]\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	$titleline =~ s/ /\t/g;
    	my @array_titleline = split("\t",$titleline);    	
		chomp $sequence;
		$sequence =~ s/\s+//g;
		$sequence =~ s/\r//g;
		$titleline =~ s/\r//g;
	   	$hash{$array_titleline[0]}++;
	}
	close(FILE); $/ = "\n"; 
	return \%hash;
}

sub input_option_hash {
	
	my %input_options = (
		1 =>'454_sff', 
		2 =>'454_fastq', 
		3=>'454_multiple_fastq', 
		4 => 'Illumina', 
		5 => 'Illumina_multiple_fastq', 
		6 => 'Illumina_pair_end', 
		7 => 'Illumina_pair_end_multiple_fastq'
	);
	
	return \%input_options;
}

sub ambiguity_DNA_codes {
	my %ambiguity_DNA_codes = (
		"R" => ["A","G"], "Y" => ["C","T"], "K" => ["G","T"], "M" => ["A","C"], "S" => ["G","C"], 
		"W" => ["A","T"], "B" => ["G","C","T"], "D" => ["A","G","T"], "H" => ["A","C","T"], 
		"V" => ["A","C","G"], "N" => ["A","C","G","T"]
	);
	return \%ambiguity_DNA_codes;
}

sub get_amb_code {
	my $ref_hash = $_[0];
	
	my %ambiguity_DNA_codes_reverse;
	my %ambiguity_DNA_codes = %{ &ambiguity_DNA_codes() };
	foreach my $keys (sort keys %ambiguity_DNA_codes) {
		my @array = sort @{$ambiguity_DNA_codes{$keys}};
		my $string = join "",@array;
		$ambiguity_DNA_codes_reverse{$string} = $keys;	
	}

	my @array = keys %{$ref_hash};
	my @array_sorted = sort @array;
	my $string = join "",@array_sorted;
	if ($ambiguity_DNA_codes_reverse{$string}) {
		return $ambiguity_DNA_codes_reverse{$string};
	} else { return "N"; }		
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

sub check_marker_pairwise {

	## Given a MSA file in a hash, check each taxa pairwise
	my $hash_ref = $_[0];
	my $MCT = $_[1];
	my $variable_positions_user_min = $_[2]; 
	my $variable_positions_user_max = $_[3];
	my $variable_divergence = $_[4];
	my $polymorph = $_[5];
	
	my @taxa = keys %$hash_ref;
	my (%seen, %pairwise, %discard);
	for (my $j=0; $j < scalar @taxa; $j++) {
		my $reference = $taxa[$j];
		#if (!$domino_files{$reference}{'taxa'}) {next;}
		foreach my $keys (sort keys %$hash_ref) {
			#if (!$domino_files{$keys}{'taxa'}) {next;}
			if ($seen{$keys}) {next;}
			if ($discard{$keys}) {next;}
			if ($keys eq $reference) {next;}
			
			my $seq2check = $$hash_ref{$keys};			
			my @array_reference = split("",$$hash_ref{$reference});
			my @array2check = split("", $$hash_ref{$keys});
			my @array;
			
			for (my $i=0; $i < scalar @array_reference; $i++) {
				my $reference_nuc = $array_reference[$i];
				my $base2check = $array2check[$i];
				push(@array, &check_reference_bp($reference_nuc, $base2check, ));
			}
			
			my $string = join "", @array;
			my $var_sites_sub = $string =~ tr/1/1/;
			my $con_sites_sub = $string =~ tr/0/0/;
			my $total_sub = $con_sites_sub + $var_sites_sub;			
			
			if ($var_sites_sub == 0) { ## If does not fit the necessary divergence
				$seen{$reference}++; $discard{$keys}++; $pairwise{$reference}{$keys} = "NO!";
				if ($variable_divergence) { 
					#&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t0\t".$variable_divergence);
				} else { 					
					#&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t0\t".$variable_positions_user_range);
				}
				next;
			}
			
			my $percentage_sub = $var_sites_sub/$total_sub;
			if ($variable_divergence) { 
				#&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t".$percentage_sub."\t".$variable_divergence);
			} else {					
				#&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t".$percentage_sub."\t".$variable_positions_user_range);
			}
			
			if ($variable_positions_user_min) {
				if ($var_sites_sub > $variable_positions_user_min) { 
					if ($var_sites_sub < $variable_positions_user_max) { 
						$pairwise{$reference}{$keys} = "YES!";
						#&debugger_print("YES");
					} else { 
						#&debugger_print("NO");
					}
				} else {
					## If does not fit the necessary divergence
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
					$pairwise{$reference}{$keys} = "NO!";
					#&debugger_print("NO");
				}
			} elsif ($variable_divergence) {
				if ($percentage_sub > $variable_divergence) { 
					$pairwise{$reference}{$keys} = "YES!";
					#&debugger_print("YES");
				} else {
					## If does not fit the necessary divergence
					#&debugger_print("NO");
					$pairwise{$reference}{$keys} = "NO!";
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
				}
			}
		}
		$seen{$reference}++; ## avoid checking again the same pair
	}
	my $flag_fitting = 0;
	foreach my $keys (sort keys %pairwise) {
		if ($discard{$keys}) {next;}
		foreach my $k (sort keys %{$pairwise{$keys}}) {
			if ($discard{$k}) {next;}
			if ($pairwise{$keys}{$k} eq 'YES!') {
				$flag_fitting++;
	}}}

	#print Dumper \%pairwise;
	if ($number_sp == 2) {
		if ($flag_fitting == 1) {
			#print "YES!\n"; 
			return '1'; 
		} else {
			#print "NO!\n";
			return '0'; 
	}} else {
		if ($flag_fitting < $MCT) {
			#print "NO!\n";
			return '0'; 
		} else {
			#print "YES!\n"; 
			return '1'; 
	}}	
}

sub check_marker_ALL {
	##########################################################################################
	##	This function checks each region using all taxa provided and generates the 			##	
	##	variation profile for the whole set													##
	## 		        																		##
	##	jfsanchezherrero@ub.edu 09/02/2016													##
	## 		        																		##
	##########################################################################################

	my $file = $_[0];	
	my $ref = $_[1];
	my $missing_allowed = $_[2];
	my $polymorphism_user = $_[3];
	my $dnaSP_flag = $_[4];
	
	#print "check_marker_ALL $file\n";
	my (%hash, $length, @taxa, @length_seqs);
	if ($ref) {
		## get alignment from hash
		foreach my $seqs (sort keys %{ $file }) {
			my @array = split("", $$file{$seqs});
			#if (!$domino_files{$seqs}{'taxa'}) {next;}
			push (@{ $hash{$seqs}}, @array);
			$length = scalar @array;
			push (@length_seqs, $length);
			push (@taxa, $seqs);
		}	
	} else {
		## get alignment from file
		open(FILE, $file) || die "Could not open the file $file ... [DM_MarkerScan: check_marker_ALL] \n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence;
			$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
			$titleline =~ s/\r//g;
			#print $titleline."\n".$sequence."\n";
			my @array = split("", $sequence);
			if (!$domino_files{$titleline}{'taxa'}) {next;}
			push (@{ $hash{$titleline}}, @array);
			$length = scalar @array;
			push (@length_seqs, $length);
			push (@taxa, $titleline);
		}
		close(FILE); $/ = "\n";	
	}

	## Get DNA code
	my %ambiguity_DNA_codes = %{ &ambiguity_DNA_codes() };

	#print Dumper \%hash; ## get into a hash a value [taxa] with an array the marker base by base
	my @tmp_length_uniq = do { my %seen; grep { !$seen{$_}++ } @length_seqs };
	if (scalar @tmp_length_uniq > 1) {
		&printError("There is problem: length of the markers do not match for $file..."); return "";
	} else { $length_seqs[0] = $length; }	
	
	my @profile;
	for (my $i=0; $i < $length; $i++) {
		my $flag_position = 0;
		my @tmp;
		foreach my $seqs (sort keys %hash) { push (@tmp, $hash{$seqs}[$i]); }
		my @tmp_uniq_sort = do { my %seen; grep { !$seen{$_}++ } @tmp };

		if (scalar @tmp_uniq_sort == 1) {
			if ($tmp_uniq_sort[0] eq 'N') { 						push (@profile, 'N');
			} elsif ($tmp_uniq_sort[0] eq '-') { 					push (@profile, '-');
			} elsif ($ambiguity_DNA_codes{$tmp_uniq_sort[0]}) { 	push (@profile, '1');
			} else {  												push (@profile, '0'); 
			}
		} else {
			## We are assuming the calling has been correctly done and
			## the ambiguity codes are due to polymorphism
			my $escape_flag = 0;
			my (@tmp, @amb);
			for (my $j=0; $j < scalar @tmp_uniq_sort; $j++) {
				if ($tmp_uniq_sort[$j] eq '-') { ## Gaps would be codify as -
					$escape_flag++;
				} elsif ($tmp_uniq_sort[$j] eq 'N') { ## Gaps would be codify as -
					$escape_flag++;
				} elsif ($ambiguity_DNA_codes{$tmp_uniq_sort[$j]}) {
					push(@amb, $tmp_uniq_sort[$j]);
				} else { 
					push(@tmp, $tmp_uniq_sort[$j]); 
			}}

			if ($escape_flag) {  push (@profile, '-');			
			} else {
				if (scalar @amb == 0) { ## No ambiguous code
					push (@profile, '1');			
				} elsif (scalar @amb == 1) { ## 1 amb code
					for (my $h=0; $h < scalar @amb; $h++) {
						my $flag_yes = 0;
						for (my $k = 0; $k < scalar @{ $ambiguity_DNA_codes{$amb[$h]} }; $k++) {
							if (grep /$ambiguity_DNA_codes{$amb[$h]}[$k]/, @tmp) { $flag_yes++; }
						}
						if ($flag_yes > 0) {
							if ($polymorphism_user) { 	
								push (@profile, '1'); 	## if polymorphism
							} else { 					
								push (@profile, '0'); ## The ambiguous is the present snps: 	YCT => [ Y > C/T ]
						}} else { 						
							push (@profile, '1');	## The ambiguous is not the present snps  	MGG => [ M > A/C ]
				}}} elsif (scalar @amb > 1) {  			
					push (@profile, '1');   ## Several
	}}}}
	my $string = join ("", @profile); undef %hash; undef @profile;

	my $var_sites = $string =~ tr/1/1/; ## count variable sites
	my $con_sites = $string =~ tr/0/0/; ## count conserved sites
	my $species = join (",", sort @taxa);
	my $count_length = $con_sites + $var_sites;
	if ($var_sites == 0) { 
		#print "NO: $species, $var_sites, $length, $string, $count_length\n";
		#if ($dnaSP_flag) { } else { return 'NO'; } 
		#if we get to provide these markers there is no need to do DOMINO as we will be reporting everything
		return 'NO';
	}
	my $missing = $length - $count_length;
	my $missing_allowed_length = $missing_allowed * $length;
	if ($missing > $missing_allowed_length) { 
		#print "NO: $species, $var_sites, $length, $string, $count_length\n";
		if ($dnaSP_flag) {} else { return 'NO';  }
		#if we get to provide these markers there is no need to do DOMINO as we will be reporting everything
		return 'NO';
	}
	my @array = ($species, $var_sites, $length, $string, $count_length);
	#print Dumper \@array; print "\n";
	return \@array;
}

sub check_reference_bp {

	my $reference_nuc = $_[0];
	my $base2check = $_[1];
	my $polymorphism_user = $_[2];
	
	## Get DNA code
	my %ambiguity_DNA_codes = %{ &ambiguity_DNA_codes() };

	## For this taxa
	if ($base2check eq "-" || $reference_nuc eq "-") { 
		return '-';
	} elsif ($base2check eq "N" || $reference_nuc eq "N") { 
		return 'N';
	} elsif ($reference_nuc ne $base2check) { 
		## Check wether there is an ambiguity code or not
		## and also if user would like some polymorphism
		## and decide whether it is a variable or conserved position
		my $flag = 1;
		if ($ambiguity_DNA_codes{$reference_nuc} || $ambiguity_DNA_codes{$base2check}) { ## one or the other or both
			if ($ambiguity_DNA_codes{$reference_nuc} and $ambiguity_DNA_codes{$base2check}) {
				## Both bases are ambiguous
				for (my $h = 0; $h < scalar @{ $ambiguity_DNA_codes{$base2check}}; $h++) {
					for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$reference_nuc}}; $j++) {
						if ($ambiguity_DNA_codes{$reference_nuc}[$j] eq $ambiguity_DNA_codes{$base2check}[$h]) {
							$flag = 0;
						} else {
							if ($polymorphism_user) {
								$flag = 1;	last;
			}}}}
			} elsif ($ambiguity_DNA_codes{$reference_nuc}) {
				for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$reference_nuc}}; $j++) {
					if ($ambiguity_DNA_codes{$reference_nuc}[$j] eq $base2check) {
						$flag = 0;								
					} else {
						if ($polymorphism_user) {
							$flag = 1;	last;
				}}}
			} elsif ($ambiguity_DNA_codes{$base2check}) {
				for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$base2check}}; $j++) {
					if ($ambiguity_DNA_codes{$base2check}[$j] eq $reference_nuc) {
						$flag = 0;								
					} else {
						if ($polymorphism_user) {
							$flag = 1;	last;
				}}}}
		} else { 
			## If neither the reference or the base to check are ambiguous
			## and they are different, this would be a variable site
			$flag = 1;
		}
		if ($flag == 0) { return '0'; } else { return '1'; }

	} elsif ($reference_nuc eq $base2check) { 					## Both bases are the same, this is a conserved site
		return '0';
	} else {
		return '-';
	}
}

1;
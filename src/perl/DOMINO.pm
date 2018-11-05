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

sub print_Excel {

	my $markers_file_ref = $_[0]; 
	my $path = $_[1];
	
	my $domino_Scripts_excel = $domino_Scripts."/DM_PrintExcel.pl";
	system("perl $domino_Scripts_excel $folder_abs_path $$path $$markers_file_ref");
	my $excel_woorkbook_name = $$path."/DM_markers-summary.xls";
	my $markers_count;
	open (MARKER, "<$$markers_file_ref"); while (<MARKER>) { if ($_ =~ /.*\_\#.*/) {$markers_count++} } close (MARKER);
	
	## Print results and instructions
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" RESULTS ","#"); DOMINO::printHeader("", "#");
	print "\n+ DOMINO has retrieved $markers_count markers\n";
	my $instructions_txt = $marker_dirname."/Instructions.txt";
	my $MSA = $$path."/MSA_markers";
	open (OUT, ">$instructions_txt");
	my $string = "+ Several files and folders has been generated:
\t+ $marker_dirname: contains DOMINO markers detected for each taxa as a reference and the clusterized results.
\t+ $$path: contains the clusterized and definitive results for DOMINO markers.
\t+ $excel_woorkbook_name: contains information of the markers identified and parameters used by DOMINO.
\t+ $MSA folder contains a single file for each marker identified.\n\n";	
	print OUT $string; print $string."\n"; close(OUT);
	return $excel_woorkbook_name;
}

sub sliding_window_conserve_variable {

	my $id = $_[0]; my $seq = $_[1]; 
	my @output_file_info;
	my $dna_seq = $$seq; my $seqlen = length($$seq);
	
	## Loop
	for (my $i=0; $i < $seqlen; $i += $SLIDING_user) {
		# AAATATGACTATCGATCGATCATGCTACGATCGATCGATCGTACTACTACGACTGATCGATCGATCGACGACTGAC
		# 		P1		P2|P3							P4|P5						P6
		my $coord_P1 = $i; #print "P1: $coord_P1\n";
		for (my $h=$domino_params{'marker'}{'window_size_CONS_min'}[0]; $h <= $domino_params{'marker'}{'window_size_CONS_max'}[0]; $h += $CONS_inc) {
			my $coord_P2 = $i + $h; 			#print "P2: $coord_P2\n";
			my $coord_P3 = $coord_P2 + 1; 		#print "P3: $coord_P3\n";
			if ($coord_P3 > $seqlen) {next;} 

			for (my $j = $domino_params{'marker'}{'window_size_VARS_min'}[0]; $j <= $domino_params{'marker'}{'window_size_VARS_max'}[0]; $j += $VAR_inc) {
				my $coord_P4 = $coord_P3 + $j; #print "P4: $coord_P4\n";
				my $coord_P5 = $coord_P4 + 1; #print "P5: $coord_P5\n";
				my $coord_P6 = $coord_P5 + $h; #print "P6: $coord_P6\n";

				if ($coord_P4 > $seqlen){next;} 
				if ($coord_P5 > $seqlen){next;} 
				if ($coord_P6 > $seqlen){next;}				

				#print "\n******\n";
				#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\nPositions:\n";
				#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
				#print "******\n\n";

				## Get substrings
				my $length_string_P1_P2 = $coord_P2 - $coord_P1;
				my $string_P1_P2 = substr($dna_seq, $coord_P1, $length_string_P1_P2);

				my $length_string_P3_P4 = $coord_P4 - $coord_P3;
				my $string_P3_P4 = substr($dna_seq, $coord_P3, $length_string_P3_P4);

				my $length_string_P5_P6 = $coord_P6 - $coord_P5;
				my $string_P5_P6 = substr($dna_seq, $coord_P5, $length_string_P5_P6);
				
				## Count the variable positions in each 
				my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//; ## Conserved 
				my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//; ## Variable
				my $count_string_P5_P6 = $string_P5_P6  =~ tr/1//; ## Conserved
				my $count_N_P1_P2 = $string_P1_P2 =~ tr/N//; ## Conserved 
				my $count_N_P5_P6 = $string_P5_P6  =~ tr/N//; ## Conserved
				
				## More than variations allowed in conserved
				if ($count_string_P1_P2 > $window_var_CONS) { next; } 
				if ($count_string_P5_P6 > $window_var_CONS) { next; }
				if ($count_N_P1_P2 > $window_var_CONS) { next; } 
				if ($count_N_P5_P6 > $window_var_CONS) { next; } 

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
					if ($expected_var_sites > $domino_params{'marker'}{'variable_positions_user_min'}[0]) {
						unless ($expected_var_sites < $domino_params{'marker'}{'variable_positions_user_max'}[0]) {
							$flag_error++;
						}} else { $flag_error++; }
				}
				
				## Missing % of bases missing
				my $missing_count = $marker_profile =~ tr/N//;
				my $missing_count2 = $marker_profile =~ tr/-//;
				$missing_count += $missing_count2;
				my $missing_count_percent = ($missing_count/$total_length)*100;
				my $percent_total_length = $total_length * $missing_allowed;  ## Default 0.05

				if ($flag_error > 0) {	next; } #if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; }
				if ($missing_count_percent < $percent_total_length) {
					push (@output_file_info, "$$id\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6");

					#if ($debugger) {
						#print "\n\n***********************************************\n";
						#print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
						#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\n\nPositions:\n";
						#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
						#print "\nSubsets\n";
						#print "Conserved (P1-P2): ".length($string_P1_P2)."\n"; print "$string_P1_P2\n";
						#print "Variable (P3-P4): ".length($string_P3_P4)."\n";  print $string_P3_P4."\n";
						#print "Conserved (P5-P6): ".length($string_P5_P6)."\n"; print $string_P5_P6."\n";
						#print "\nVariations:\n";
						#print "Count_string_P1_P2: $count_string_P1_P2\n";
						#print "Count_string_P3_P4: $count_string_P3_P4\n";
						#print "Count_string_P5_P6: $count_string_P5_P6\n";
						#print "\nMarker:\n";
						#print "Total Length:$total_length\tVariable Positions:$count_string_P3_P4\tExpected:$expected_var_sites\n";
						#print "Missing allowed: $percent_total_length %\tMissing:$missing_count_percent %\n";
						#print "***********************************************\n";
					#}
	} else { next;	}}}}
	return \@output_file_info;

}

sub check_DOMINO_marker {
	
	my $file = $_[0]; ## OUTPUT
	my $dir = $_[1];	
	my $DOMINO_markers_file = $_[2]; ## markers collapse
	my $ref_taxa_all = $_[3]; 			# if MSA alignment it is a file containing msa

	my @files; 
	
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
			if ($option eq "msa_alignment") {
				my $fasta_msa_sub_ref = DOMINO::readFASTA_hash($ref_taxa_all);
				%fasta_msa_sub = %{ $fasta_msa_sub_ref };
			} else {
				my @desired_taxa = split(",", $taxa);
				for (my $j=0; $j < scalar @desired_taxa; $j++) {
					next if ($ref_taxa_all eq $desired_taxa[$j]);
					my $fasta_each_taxa = $domino_files{$desired_taxa[$j]}{"PROFILE::Ref:".$ref_taxa_all}[0]."/".$contig_name[0]."_sequence.fasta";
					if (-f $fasta_each_taxa) {
						$fasta_msa_sub{$desired_taxa[$j]} = $fasta_each_taxa;
					} else { 
						$fasta_msa_sub{$desired_taxa[$j]} = "missing"; 
				}}
				my $reference_file_contig = $domino_files{$ref_taxa_all}{'REF_DIR'}[0]."/".$contig_name[0].".fasta";
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
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $coord_P1, $coord_P6, $sequence);
					$hash{$keys} = $seq;
				} else {
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $coord_P1, $coord_P6, $fasta_msa_sub{$keys});
					$hash{$keys} = $seq;
			}}
			my $seq_name;
			
			## Check pairwise MSA
			my $valueReturned = &check_marker_pairwise(\%hash);
			if ($valueReturned == 1) { ## if it is variable for each pairwise comparison
				## Get variable positions for the whole marker
				my $array_ref_returned = &check_marker_ALL(\%hash, "Ref");
				if ($array_ref_returned eq 'NO') { 
					remove_tree($msa_file);
				} else {
					$marker_number++;
					my $msa_file = $dir."/".$contig_name[0]."_marker_".$marker_number.".fasta";
					open (MSA_OUT, ">$msa_file");
					foreach my $seqs (sort keys %hash) { print MSA_OUT ">".$seqs."\n".$hash{$seqs}."\n"; }
					close (MSA_OUT);
					push (@files, $msa_file);
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
	return (\@files);
}

sub check_overlapping_markers {

	my $file = $_[0]; 				## markers_shared 
	my $mergeArray_file = $_[1];

	my @array = split(".txt", $file); 
	my $file2return = $array[0]."_overlapped_Markers.txt";

	my $domino_Scripts_MarkerOverlap = $domino_Scripts."/DM_MarkerOverlap.pl";
	my $command = "perl $domino_Scripts_MarkerOverlap $file $mergeArray_file $file2return $folder_abs_path"; #print $command."\n"; 
	system($command);
	
	return $file2return;	
}

1;
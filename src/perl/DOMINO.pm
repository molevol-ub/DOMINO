#!/usr/bin/perl
####################################################################
###	DOMINO: Development of molecular markers in non-model organisms 
####################################################################
# This package provides multiple subroutines for the DOMINO package
## date 7/11/16

package DOMINO;

use FindBin;
use lib $FindBin::Bin."/lib";
require List::MoreUtils;
use List::MoreUtils qw(firstidx);

sub blastn {
	my $file = $_[0]; my $db = $_[1]; my $results = $_[2]; my $BLAST = $_[3]; 
	my $filter = $BLAST."blastn -query ".$file." -evalue 1e-10 -db '".$db."' -out $results -outfmt 6";
	my $message = "BLASTN command: $filter\n"; 
	my $blastn = system($filter);
	return ($blastn, $message);
}

sub check_ID_length {
	
	my $file = $_[0];
	my $option = $_[1];
	open (F1, "$file") or die "Could not open file $file for reading ID length";
	
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
    open(FILE, $file) or die "Could not open file $file when checking file format...";
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
		open (F1, "$file") or die "Could not open file $file for checking paired-end reads";
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
		open (FILE, "$file") or die "Could not open file $file when checking paired-end reads";
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
		} else { die "The reads provided within the file $file do not belong to the same pair...";}
		if ($types[0] eq $types[$i]) {
			$type = $types[$i];
		} else { die "The reads provided within the file $file do not belong to the same pair...";}
	}
	return ($pair_return, $type);
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

sub dieNicely {
	print DOMINO::time_stamp();	print "\n\nTry perl $0 -man for more information\n\n";
	Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 0 );
}

sub fasta_file_splitter {
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	my $file = $_[0];
	my $block = $_[1];
	my $ext = $_[2]; # fasta, fastq, loci, fa
	my $dir = $_[3];

	open (FH, "<$file") or die "Could not open source file. $!";
	print "\t- Splitting file into blocks of $block characters aprox ...\n";
	my $j = 0; my @files;
	while (1) {
		my $chunk;
	   	my @tmp = split ("\.".$ext, $file);
		my $block_file = $tmp[0]."_part-".$j."_tmp.".$ext;
		push (@files, $block_file);
		open(OUT, ">$block_file") or die "Could not open destination file";
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
	open (FH, "<$file") or die "Could not open source file. $!";
	print "+ Splitting file $file into blocks of $block characters...\n";
	my $j = 0; 
	while (1) {
    	my $chunk;
    	my @tmp = split (".".$ext, $file);
		my $file_name = $tmp[0];
		
	   	my $block_file = $file_name."_part-".$j."_tmp.".$ext;
		print "\t- Printing file: ".$block_file."\n";
    	push (@files, $block_file);
    	open(OUT, ">$block_file") or die "Could not open destination file";
    	$j++;
    	if (!eof(FH)) { read(FH, $chunk,$block);  print OUT $chunk; } ## Print the amount of chars
    	if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken
    	close(OUT); last if eof(FH);
	}
	close(FH);
	return (\@files);	
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
			if (!$earliest) { $earliest=$time_log; } else { if ($time_log > $earliest) {$earliest=$time_log}}
			$mapping_dirs{$time_log} = $folder."/".$array_files[$i];
	}}
	if (!exists $mapping_dirs{$earliest}) { return 'NO';
	} else {  return $mapping_dirs{$earliest}; }	
}

sub get_size { # Multiplatform	
	my $size = -s $_[0]; #To get only size
	return $size;	
	# only UNIX
	# Get lines, words and characters
	#my $wc_out = `wc $_[0]`;
	#$wc_out =~ s/\s/-/g;
	#my $array = &splitter("\-", $wc_out);
	#my @word_count;
	#for (my $i=0; $i < scalar @$array; $i++) {
	#	if ($$array[$i] =~ /\S+/) {
	#		push (@word_count, $$array[$i]);			
	#	}
	#} return \@word_count;
}

sub line_splitter {
	# given a separator, splits a given line
	my $separator = $_[0];
	my @array = split($separator, $_[1]); 
	return \@array; 
}

sub loci_file_splitter {
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	my $file = $_[0];
	my $block = $_[1];
	my $ext = $_[2]; # fasta, fastq, loci, fa

	open (FH, "<$file") or die "Could not open source file. $!";
	print "\t- Blocks of $block characters would be generated...\n";
	my $j = 0; my @files;
	while (1) {
		my $chunk;
	   	my @tmp = split (".".$ext, $file);
		my $file_name = $tmp[0];
		my $block_file = $file_name."_part-".$j."_tmp.".$ext;
		print "\t- Printing file $block_file\n";
		push (@files, $block_file);
		open(OUT, ">$block_file") or die "Could not open destination file";
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

sub readFASTA_hash {

	my $file = $_[0];
	my %hash;
	open(FILE, $file) || die "Could not open the $file ...\n";
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
	open(FILE, $file) || die "Could not open the $file ...\n";
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
	open (F1, "$file") or &printError("Could not open file $file") and exit();
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
	open(FILE, $file) || die "Could not open the $file ...\n";
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

sub time_stamp { return "[ ".(localtime)." ]"; }



## TODO
sub seq_counter {
	
	my $file = $_[0];	
	my $option_format = DOMINO::check_file_format($file);
	my ($nSequences, $nLines);
	open (F1, "$file") or &printError("Could not open file $file") and exit(); 
	while (<F1>) { $nLines++; } close(F1);
	if ($option_format eq "fastq") { 
		$nSequences = int ($nLines/4);
	} elsif ($option_format eq "fasta") {
		$nSequences = int ($nLines/2);
	}
	return $nSequences;
}

1;
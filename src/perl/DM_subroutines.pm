#!/usr/bin/perl
package DM_subroutines;
# This package provides multiple subroutines for the DOMINO package 

sub time_stamp {
	my $current_time = time;
	print "[ ".(localtime)." ]\t";
	my $secs = $current_time - $step_time; 
	my $hours = int($secs/3600); 
	$secs %= 3600; 
	my $mins = int($secs/60); 
	$secs %= 60; 
	$step_time = $current_time;
	printf ("Step took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub printError {
    my $msg = $_[0];
	print "\n\n";&print_Header(" ERROR ","!!"); print "\n";
    print $msg."\n\nTry \'perl $0 -h|--help or -man\' for more information.\nExit program.\n";
	print "\n\n"; &print_Header("","!!"); &print_Header("","!!"); 
    &printError_log($msg);
}

sub printError_log {
	my $message = $_[0];
	open (ERR, ">>$error_log");
	print ERR $message."\n";
	close (ERR);
	#print STDERR $message."\n";
}

sub finish_time_stamp {

	my $finish_time = time;
	print "\n\n"; &print_Header("","+"); 
	&print_Header(" ANALYSIS FINISHED ","+"); 
	&print_Header("","+"); 
	print "[".(localtime)." ]\t";
	my $secs = $finish_time - $start_time; 
	my $hours = int($secs/3600); $secs %= 3600; 	
	my $mins = int($secs/60); $secs %= 60; 
	printf ("Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub printFormat_message {
	print "\n\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\nWhere:\n\txxx: any character or none.Please avoid using dots (.)\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n\n";
}

sub print_Header {

	my $sentence = $_[0];
	my $symbol = $_[1];	
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

sub dieNicely {
	print "\n"; &time_stamp();	print "\n";
	Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 0 );
}

sub blastn {

	##########################################################################################
	##	 																					##
	##  This function uses BLASTN to check for the putative contaminants					##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 	jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################

	my $file = $_[0];
	my $db;
	if ($file =~ /(.*)\.fasta/ ) { $db = $1 ; } 
	my $filter = $BLAST."blastn -query ".$file." -evalue 1e-5 -db ".$db." -out blast_search.txt -outfmt 6";
	my $blastn = system($filter);
	return $blastn;
}

sub makeblastdb {

	my $file = $_[0];
	my $db;
	my $make_blast_db;
	$make_blast_db = $BLAST."makeblastdb";
	$make_blast_db .= " -in ".$file;
	$make_blast_db .= " -dbtype nucl";
	if ($file =~ /(.*)\.fasta/ ) { $db = $1 ; } 
	$make_blast_db .= " -out ".$db;
	$make_blast_db .= " -logfile ".$db.".log 2> $error_log";	
	my $makeblastresult = system($make_blast_db);
	if ($makeblastresult != 0) { &printError("Generating the database failed when trying to proccess the file... DOMINO would not stop in this step...\n"); }
	return ($db, $makeblastresult);
}

sub read_dir {
	my $dir = $_[0];
	opendir(DIR, $dir);
	my @dir_files = readdir(DIR);
	my $array_ref = \@dir_files;
	return $array_ref;
}

sub mothur_remove_seqs {	
	# This subroutine takes as input a FASTQ file AND classifies according to the tags provided, by Roche and trims the seqs
	my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
	my $fasta = $_[0];
	my $qual = $_[1];
	my $directory = $_[2];
	my $ids2remove = $_[3];
	
	my $line = $mothur_path." '#set.dir(output=$directory); remove.seqs(accnos=$ids2remove, fasta=$fasta, qfile=$qual)'";
	print "\n+ Calling mothur executable for discarding reads...\n\n";
	my $system_call = system($line);

}

sub mothur_retrieve_seqs {
	## This sub retrieve a given number of ids and generates a fastq
	my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
	my $fasta = $_[0];
	my $qual = $_[1];
	my $directory = $_[2];
	my $ids2retrieve = $_[3];
	
	my $line = $mothur_path." '#set.dir(output=$directory); get.seqs(accnos=$ids2retrieve, fasta=$fasta, qfile=$qual)'";
	print "\n+ Calling mothur executable for retrieving reads in file $ids2retrieve...\n\n";
	my $system_call = system($line);
}

sub line_splitter {
	# given a separator, splits a given line
	my $separator = $_[0];
	my @array = split($separator, $_[1]); 
	return \@array; 
}

sub get_size {
	#my $size = -s $_[0]; To get only size
	# Get lines, words and characters
	my $wc_out = `wc $_[0]`;
	$wc_out =~ s/\s/-/g;
	my $array = &splitter("\-", $wc_out);
	my @word_count;
	for (my $i=0; $i < scalar @$array; $i++) {
		if ($$array[$i] =~ /\S+/) {
			push (@word_count, $$array[$i]);			
		}
	}
	return \@word_count;
}

sub file_splitter {
	# Splits a file such a sam or whatever file that could be read for each line
	open (FH, "<$file") or die "Could not open source file. $!";
	print "\n\nSplitting file into blocks of $block characters...\n";
	my $j = 0; 
	while (1) {
    	my $chunk;
    	print "Processing block $j...\n";
    	my @tmp = split (".txt", $file);
		my $file_name = $tmp[0];
		
	   	my $block_file = $file_name."_part-".$j."_tmp.txt";
    	push (@files, $block_file);
    	open(OUT, ">$block_file") or die "Could not open destination file";
    	$j++;
    	if (!eof(FH)) { read(FH, $chunk,$block);  print OUT $chunk; } ## Print the amount of chars
    	if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken
    	close(OUT); last if eof(FH);
	}
	close(FH);
}

sub fasta_file_splitter {
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	print "\n\nSplitting file into blocks of $block characters...\n";
	my $j = 0; my @files;
	while (1) {
		my $chunk;
		print "Processing block $j...\n";
		my $file_name; my $extension;
		if ($file =~ /(.*)\.(.*)/) {
			$file_name = $1;
			$extension = $2;
		}
		my $block_file = $file_name."_part-".$j."_tmp.".$extension;
		push (@files, $block_file);
		open(OUT, ">$block_file") or die "Could not open destination file";
		if (!eof(FH)) { read(FH, $chunk,$block);  
			if ($j > 0) {
				$chunk = ">".$chunk;
			}
			print OUT $chunk;
		} ## Print the amount of chars
	
		if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken
	
		if (!eof(FH)) { 
			$/ = ">"; ## Telling perl where a new line starts
			$chunk = <FH>; 
			chop $chunk;
			print OUT $chunk; 
			$/ = "\n";
		} ## print the sequence if it is broken
		$j++;
		close(OUT); last if eof(FH);
	}
	close(FH);
}

1;
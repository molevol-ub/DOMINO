#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin."/../lib";
use DOMINO;
use Data::Dumper;

my $fasta = $ARGV[0];
my $output = $ARGV[1];

if(!defined($fasta)) { print "ERROR: No input files are provided\n"; exit; } ## give error message for DOMINO debug

my ($total_Contigs_all_sets, 
	$total_GOOD_Contigs_all_sets, $all_bases, 
	$all_GOOD_bases, @all_contigs, @all_GOOD_contigs, %parts_array);

my @parts = (150, 500, 1000, 5000, 10000);
my @file = split("/", $fasta);
my @name = split("\.fasta", $file[-1]);
my $discarded=0; 
my $gaps = 0; my $gap_length = 0; my $total_len = 0;

open(FILE, $fasta) || die "Could not open the $fasta ...\n"; ## give error message for DOMINO debug
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {		
	next if /^#/ || /^\s*$/;
	chomp;
	my ($titleline, $sequence) = split(/\n/,$_,2);
	next unless ($sequence && $titleline);
    $sequence =~ s/\n//g;
    $sequence = uc($sequence);
	
	## Get total bases A+T+C+G
	my $len = length $sequence;
	$total_len += $len;
	my @bases_count = &baseCount($sequence);
	my $id;
	for (my $j = 0; $j < scalar @parts; $j++) {
		if ($len <= $parts[0]) {
			$id = "0";
			push (@ {$parts_array{$id}{"length"} }, $len);				
			$parts_array{$id}{"TOTAL_COUNT"} += $len;
			$parts_array{$id}{"COUNT"}++;
			$parts_array{$id}{"nucleotides"}{"A"} += $bases_count[0];
			$parts_array{$id}{"nucleotides"}{"T"} += $bases_count[1];
			$parts_array{$id}{"nucleotides"}{"C"} += $bases_count[2];
			$parts_array{$id}{"nucleotides"}{"G"} += $bases_count[3];
			$parts_array{$id}{"N_COUNT"} += $bases_count[4];
			last;
		} elsif ($len <= $parts[$j]) {
			#my $id = "$parts[$j - 1] - $parts[$j]";
			$id = $parts[$j - 1];
			push (@ {$parts_array{$id}{"length"} }, $len);				
			$parts_array{$id}{"TOTAL_COUNT"} += $len;
			$parts_array{$id}{"COUNT"}++;
			$parts_array{$id}{"nucleotides"}{"A"} += $bases_count[0];
			$parts_array{$id}{"nucleotides"}{"T"} += $bases_count[1];
			$parts_array{$id}{"nucleotides"}{"C"} += $bases_count[2];
			$parts_array{$id}{"nucleotides"}{"G"} += $bases_count[3];
			$parts_array{$id}{"N_COUNT"} += $bases_count[4];
			last;
		} elsif ($len > $parts[-1]) { ## Get contigs >50kb (if any)
			## Push rank
			$id = $parts[-1];
			push (@ {$parts_array{$id}{"length"} }, $len);				
			$parts_array{$id}{"TOTAL_COUNT"} += $len;
			$parts_array{$id}{"COUNT"}++;
			
			$parts_array{$id}{"nucleotides"}{"A"} += $bases_count[0];
			$parts_array{$id}{"nucleotides"}{"T"} += $bases_count[1];
			$parts_array{$id}{"nucleotides"}{"C"} += $bases_count[2];
			$parts_array{$id}{"nucleotides"}{"G"} += $bases_count[3];
			$parts_array{$id}{"N_COUNT"} += $bases_count[4];
			last;
		}
	}	
	## Check for NNNNN gaps in the sequences
	while ($sequence =~ /(N{25,})/g) {
		my $gapLen = length($1);
		$gaps++; $gap_length += $gapLen;
	}
}
close(FILE);
$/ = "\n";

my %ALL_nucleotides = ("A" => 0, "C" => 0, "T" => 0, "G" => 0, "N" => 0);
foreach my $keys (sort {$a<=>$b} keys %parts_array) {
	$ALL_nucleotides{"A"} += $parts_array{$keys}{"nucleotides"}{"A"};
	$ALL_nucleotides{"T"} += $parts_array{$keys}{"nucleotides"}{"T"};
	$ALL_nucleotides{"C"} += $parts_array{$keys}{"nucleotides"}{"C"};
	$ALL_nucleotides{"G"} += $parts_array{$keys}{"nucleotides"}{"G"};
	$ALL_nucleotides{"N"} += $parts_array{$keys}{"N_COUNT"};
	$all_bases += $parts_array{$keys}{"TOTAL_COUNT"};
	$total_Contigs_all_sets += $parts_array{$keys}{"COUNT"};
}

##	Set	average   median N95	Count_Genomic_Contigs  		CountContigs_%  pb_this_set	pb_this_set_%	Reads_RNA_Count			 Reads_RNA_Count_%	Unmapped_contigs_Count				Unmapped_contigs_%							
DOMINO::printHeader("","#"); DOMINO::printHeader(" Assembly Statistics ","#");  DOMINO::printHeader("","#"); 
my $csv_file = $output;
print "Assembly Statisitcs for file: $fasta\n";
print "Printing statistics in CSV file: $csv_file\n\n";
print "## General Statistics ##\n";
print "Total sequences: $total_Contigs_all_sets\n";
printf "Total bases: %.3e\n", $all_bases; 
my $all_bases_noNs = $all_bases - $ALL_nucleotides{"N"};
my $total_gaps_length = $gap_length/$all_bases*100;
printf "Total length (no Ns): %.3e\n", $all_bases_noNs; print "\n";

printf "%-25s %0.2f %s\n", "As", $ALL_nucleotides{"A"}/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "Ts", $ALL_nucleotides{"T"}/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "Cs", $ALL_nucleotides{"C"}/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "Gs", $ALL_nucleotides{"G"}/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "(A + T)s", ($ALL_nucleotides{"A"} + $ALL_nucleotides{"T"})/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "(G + C)s", ($ALL_nucleotides{"G"} + $ALL_nucleotides{"C"})/$all_bases*100, "%";
printf "%-25s %0.2f %s\n", "Ns", $ALL_nucleotides{"N"}/$all_bases*100, "%";
printf "%-25s %i \n", "Capture Gaps", $gaps;
printf "%-25s %i \n", "Capture Gaps Length", $gap_length;
printf "%-25s %0.2f %s\n", "Capture Gaps Length/Total Length (%)", $total_gaps_length, "%";

open (CSV, ">$csv_file");
print CSV "Sequences,$total_Contigs_all_sets\n";
print CSV "Total Length (bp),$all_bases\n";
print CSV "Total length (no Ns),".$all_bases_noNs."\n";
print CSV "A,".( $ALL_nucleotides{"A"}/$all_bases*100)."\n";
print CSV "T,".($ALL_nucleotides{"T"}/$all_bases*100)."\n";
print CSV "C,".($ALL_nucleotides{"C"}/$all_bases*100)."\n";
print CSV "G,".($ALL_nucleotides{"G"}/$all_bases*100)."\n";
print CSV "A+T,".( ($ALL_nucleotides{"A"} + $ALL_nucleotides{"T"})/$all_bases*100)."\n";
print CSV "C+G,".( ($ALL_nucleotides{"G"} + $ALL_nucleotides{"C"})/$all_bases*100)."\n";
print CSV "N,".($ALL_nucleotides{"N"}/$all_bases*100)."\n";
print CSV "Capture Gaps,".$gaps."\n";
print CSV "Capture Gaps Length,".$gap_length."\n";
print CSV "Capture Gaps Length/Total Length (%),".$total_gaps_length."\n";

print CSV "Set,Number Seqs,% Seqs,Total Length (bp),% Bases,MinLen,MaxLen,Average Len,Median Len,n50,N_set,N_total,Length (no Ns)\n";
DOMINO::printHeader("","#"); 

foreach my $keys (sort {$a<=>$b} keys %parts_array) {
	my @array;
	my $total_pb_this_set = 0;
	my $total_N_this_set = 0;
	foreach my $keys2 (sort keys %parts_array) {
		if ($keys2 < $keys) {next;}
		push (@array, @{ $parts_array{$keys2}{"length"}});
		$total_pb_this_set += $parts_array{$keys2}{"TOTAL_COUNT"};
		$total_N_this_set += $parts_array{$keys2}{"N_COUNT"};
	}
	my $set = ">$keys bp"; 
	print "\nSet: $set\n";
	my $string2print_set = &get_stats(\@array, $all_bases, $total_Contigs_all_sets);

	my ($N_total_perc, $N_set_perc) = 0; 
	if ($total_N_this_set) {
		$N_total_perc = sprintf ("%0.2f",($total_N_this_set/$all_bases*100));
		$N_set_perc = sprintf ("%0.2f",($total_N_this_set/$total_pb_this_set*100));
		printf "%-25s %0.2f\n", "Ns (%set) ", $N_set_perc;
		printf "%-25s %0.2f\n", "Ns (%total) ", $N_total_perc;
	}
	my $all_bases_noNs_set = $total_pb_this_set - $total_N_this_set;
	if (!$N_set_perc) {$N_set_perc=0;}
	print CSV $set.",".$string2print_set.",".$N_set_perc.",".$N_total_perc.",".$all_bases_noNs_set."\n";	
	printf "%-25s %.3e\n", "Total length (no Ns) ", $all_bases_noNs_set; print "\n";	
	DOMINO::printHeader("","#"); 
}
close (CSV);


sub calcN50 {
	my @x = sort {$b<=>$a} @{$_[0]};
	my $n = $_[1];
	my $total = DOMINO::sum(\@x);
	my $total_seqs = scalar @x;
	my ($count, $n50) = (0,0);
	my $l50 = 0;
	
	for (my $j=0; $j<@x; $j++){
        $count += $x[$j];
        $l50++;
        if($count >= ($total*$n/100)){
            $n50=$x[$j];
            last;
	}}
	return ($n50, $l50);
}
	
sub calcMedian {
	my @arr = @_;
	my @sArr = sort{$a<=>$b} @arr;
	my $arrLen = @arr;
	my $median;
	if($arrLen % 2 == 0) {
		$median = ($sArr[$arrLen/2-1] + $sArr[$arrLen/2])/2;
	} else { $median = $sArr[$arrLen/2]; }
	return $median;
}

sub baseCount {
	my $seq = $_[0];
	my $tAs += $seq =~ s/A/A/gi;
	my $tTs += $seq =~ s/T/T/gi;
	my $tGs += $seq =~ s/G/G/gi;
	my $tCs += $seq =~ s/C/C/gi;
	my $Ns += (length $seq) - $tAs - $tTs - $tGs - $tCs;
	return ($tAs, $tTs, $tCs, $tGs, $Ns);	
}

sub get_stats {
	
	my $array_ref = $_[0];
	my $set_bases = $_[1];
	my $total = $_[2];
	my @array = @$array_ref;
	my $bases = DOMINO::sum($array_ref);
	
	my $percentage_pb_bases_this_set = ($bases/$set_bases)*100;
	my $percentage_pb_bases_this_set_print = sprintf "%0.5f",$percentage_pb_bases_this_set;
	my $totalContigs = scalar @array;
	my $percentage_contigs_this_set = ($totalContigs/$total)*100;
	my $percentage_contigs_returned = sprintf "%0.5f",$percentage_contigs_this_set;
	
	my ($minReadLen, $maxReadLen) = DOMINO::min_max($array_ref);
	my $avgReadLen = sprintf "%0.2f", $bases/$totalContigs;
	my $medianLen = calcMedian(@array);
	my ($n50, $L50) = calcN50($array_ref, 50);

	printf "%-25s %d\n" , "Total sequences", $totalContigs;
	printf "%-25s %0.2f\n" , "Total sequences (%)", $percentage_contigs_returned;
	printf "%-25s %.3e\n" , "Total bases", $bases;
	printf "%-25s %0.2f\n" , "Total bases(%)", $percentage_pb_bases_this_set_print;
	printf "%-25s %d\n" , "Min sequence length", $minReadLen;
	printf "%-25s %d\n" , "Max sequence length", $maxReadLen;
	printf "%-25s %0.2f\n", "Average sequence length", $avgReadLen;
	printf "%-25s %0.2f\n", "Median sequence length", $medianLen;
	printf "%-25s %0.2f\n", "N50: ", $n50;
	printf "%-25s %0.2f\n", "L50: ", $L50;
	my $string = $totalContigs.",".$percentage_contigs_returned.",".$bases.",".$percentage_pb_bases_this_set_print.",".$minReadLen.",".$maxReadLen.",".$avgReadLen.",".$medianLen.",".$n50.",".$L50;
	return $string;
}
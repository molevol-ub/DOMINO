#!/usr/bin/perl
#######################################################################################
###	DOMINO: Development of molecular markers in non-model organisms using NGS data  ###
###											###
###	Authors:									###
###	Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro 	###
###	Sánchez-Gracia, and Julio Rozas.					     	###
###											###
#######################################################################################

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin."/../lib";
my $scripts = $FindBin::Bin;
my $binaries = $FindBin::Bin."/../";

BEGIN {

	use Getopt::Long;
	use Pod::Usage;
	use Data::Dumper;
	use POSIX qw(strftime);

	use List::MoreUtils qw(firstidx);

	use DOMINO;
	use File::Copy;
	use File::Find qw(find);
	use List::Uniq qw(uniq);

	use File::Path qw(remove_tree);
	use Cwd qw(abs_path);  
	use Parallel::ForkManager;
	use Spreadsheet::WriteExcel;

	use Time::HiRes qw(usleep nanosleep);
}

my @modules = (
	"Getopt::Long", 
	"Pod::Usage",
	"Data::Dumper",
	"POSIX",
	"FindBin",
	"DOMINO",
	"File::Copy",
	"File::Find;",
	"List::Uniq",
	"File::Path",
	"Cwd",
	"Parallel::ForkManager",
	"Spreadsheet::WriteExcel",
	"Time::HiRes",
	"List::Util"
	);

$| = 1; ## Flush output after every print and do not wait until full

print "\n\n";
DOMINO::printHeader("","#"); DOMINO::printHeader(" MODULES ","#"); DOMINO::printHeader("","#");
print "\nChecking perl module dependencies...\n\n";
for (my $i=0; $i< scalar @modules; $i++) {
	print "\tChecking module: ".$modules[$i];
	for (my $i=0; $i < 20; $i++) { print "."; usleep 5000; }	
	my $new_module = $modules[$i] =~ s/\:\:/\//gr;
	if (exists $INC{$new_module.".pm"}) {
    	print $new_module.".pm [OK]\n";
	} else { print " [X]\n\n\tATTENTION: $new_module is missing but DOMINO might still work appropiate...]\n\n";}
}
print "\n\n";

## to check
DOMINO::printHeader("","#"); DOMINO::printHeader(" BINARIES ","#"); DOMINO::printHeader("","#");
print "\nChecking binary dependencies from other sources...\n\n";

my %binaries = (
	"samtools v1.3.1", $binaries."samtools-1.3.1/samtools",
	"bowtie2 v2.2.9", $binaries."bowtie2-2.2.9/",
	"BLAST", $binaries."NCBI_BLAST/",
	"mothur v1.32", $binaries."MOTHUR_v1.32.0/mothur",
	"MIRA", $binaries."mira_v4.0/bin/mira",
	"CAP3", $binaries."cap3/bin/cap3",
);

## check each binary
foreach my $path (keys %binaries) {
	print "\tChecking $path:... $binaries{$path}\n"
}

print "\n\n";

## to check
DOMINO::printHeader("","#"); DOMINO::printHeader(" UTILS ","#"); DOMINO::printHeader("","#");
print "\n\nChecking perl scripts...\n\n";

## check each perl script
my %script;
my @all_perl;
my @arrayfiles = @{ DOMINO::readDir($scripts) };
for (my $i=0; $i < scalar @arrayfiles; $i++) {
	if ($arrayfiles[$i] =~ /.*\.pl/) {
		push ( @all_perl, $scripts."/".$arrayfiles[$i]);
}}
my @arrayfiles2 = @{ DOMINO::readDir($binaries) };
for (my $i=0; $i < scalar @arrayfiles2; $i++) {
	if ($arrayfiles2[$i] =~ /.*\.pl/) {
		push ( @all_perl, $binaries.$arrayfiles2[$i]);
}}

#print Dumper \@all_perl;
for (my $j=0; $j < scalar @all_perl; $j++) {
	
	next if ($all_perl[$j] eq $scripts."/DM_MarkerSelection.pl");
	next if ($all_perl[$j] eq $scripts."/DM_ParseMSA_files.pl");
	
	print "\tChecking";
	for (my $i=0; $i < 20; $i++) { print "."; usleep 5000; }	
	
	my $command = "perl -c $all_perl[$j]";
	my $call = system($command);
	#print $command."\t"; print $call."\n";
}

print "\n\nEverything seems OK.\nFinish\n"






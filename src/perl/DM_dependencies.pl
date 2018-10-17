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

BEGIN {

	use Getopt::Long;
	use Pod::Usage;
	use Data::Dumper;
	use POSIX qw(strftime);
	use FindBin;
	use lib $FindBin::Bin."/lib";
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
	"Time::HiRes"
	);

$| = 1; ## Flush output after every print and do not wait until full

print "\nChecking perl module dependencies...\n\n";
for (my $i=0; $i< scalar @modules; $i++) {
	print "\tChecking module: ".$modules[$i];
	for (my $i=0; $i < 20; $i++) { print "."; usleep 5000; }	
	my $new_module = $modules[$i] =~ s/\:\:/\//gr;
	if (exists $INC{$new_module.".pm"}) {
    	print $new_module.".pm [OK]\n";
	} else { print " [X]\n\n\tATTENTION: $new_module is missing but DOMINO might still work appropiate...]\n\n";}
}

## to check
my $pipeline_path = abs_path($0);
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; } 

print "\nChecking binary dependencies from other sources...\n\n";

my %binaries = (
	"samtools v1.3.1", $scripts_path."samtools-1.3.1/samtools",
	"bowtie2 v2.2.9", $scripts_path."bowtie2-2.2.9/",
	"BLAST", $scripts_path."NCBI_BLAST/",
	"mothur v1.32", $scripts_path."MOTHUR_v1.32.0/mothur",
	
);

foreach my $path (keys %binaries) {
	print "\tChecking $path:...\n"
}


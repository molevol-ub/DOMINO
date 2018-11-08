#!/usr/bin/perl
###########################################################################################
###	DOMINO: Development of molecular markers in non-model organisms using NGS data	###
###											
###	Authors:								
###	Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro 
###	Sánchez-Gracia, and Julio Rozas.		
###########################################################################################
##	Usage:
##      perl DM_Clean_v1.1.pl
##
##    ###########################
##    ### General Information ###
##    ###########################
##      [-h|--help] [-man] [-v|--version] [-MoreInfo]
##
##    #########################
##    ### Mandatory options ###
##    #########################
##      [-type_file int_value] [-input_file file] [-b|--barcode_file file]
##      [-o|--outputFolder path]
##
##    ######################
##    ### Optional flags ###
##    ######################
##      ## Cleaning
##
##      [-l|--cut_Off_ReadLen int_value] [-s|--cut_Off_Quality int_value]
##      [-m|--minLen int_value] [-thr|--threshold_DUST int_value]
##
##      ## Database search
##
##      [-db|--blast_database file] [-no_db_search] [-only_user_db]
##
##      ## Others
##
##      [-bdiffs int_value] [-p|--number_cpu int_value] [-only_tag_files]
##      [-TempFiles]
##
use strict;
use warnings;
use Getopt::Long;		
use Pod::Usage;
use POSIX qw(strftime);
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/lib";
BEGIN {
	require DOMINO;
	require List::Uniq;
	require File::Copy;
	require File::Find; use File::Find qw(find);			
	require File::Path; use File::Path qw(remove_tree);
	require Cwd; use Cwd qw(abs_path);  
	require Parallel::ForkManager;
} 
my $pipeline_path = abs_path($0); 
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/";  }

##################################
## Initializing some variables	##
##################################
my (
## Options provided
$bdiffs, $read_length_GUI, $minimum_qual_GUI, $minimum_length_GUI, $outFolder, 
$manual, $version, $threshold, $noOfProcesses, $user_barcodes_file, $helpAsked, $helpAsked1,
$avoidDelTMPfiles, @user_blast_db, $skipping_BLAST, $user_option_file_type, $only_user_db, @blast_db,
$further_information, $onlyTagging_files, $debugger, @user_blast_db_tmp, @user_files,

## Others
$step_time, %domino_files, @file_abs_path, @file_names, %domino_params
);	

## Get input options
my %input_options = %{ DOMINO::input_option_hash() }; 

######################
## Get user options ##
######################
GetOptions(
	"h" => \$helpAsked1,
	"help" => \$helpAsked,
	"man" => \$manual,
	"v|version" => \$version,
	"MoreInfo" => \$further_information,	
	######################

	"type_file=i" => \$user_option_file_type,
	"input_file=s" => \@user_files,  
	"l|cut_Off_ReadLen=i" => \$read_length_GUI,
	"s|cut_Off_Quality=i" => \$minimum_qual_GUI,
	"m|minLen=i" => \$minimum_length_GUI,
	"thr|threshold_DUST=i" => \$threshold,
	"o|outputFolder=s" => \$outFolder, 
	"p|number_cpu=i" => \$noOfProcesses,
	"db|blast_database=s" => \@user_blast_db_tmp,			
	######################

	"bdiffs=i" => \$bdiffs,
	"b|barcode_file=s" => \$user_barcodes_file,			
	"TempFiles" => \$avoidDelTMPfiles,
	"no_db_search" => \$skipping_BLAST,
	######################

	"only_user_db" => \$only_user_db,
	"only_tag_files|OTG" => \$onlyTagging_files,			
	######################

	"Debug" => \$debugger,
);

## Help/manual/version
pod2usage( -exitstatus => 0, -verbose => 1 ) if ($helpAsked);
pod2usage( -exitstatus => 0, -verbose => 0 ) if ($helpAsked1);
pod2usage( -exitstatus => 0, -verbose => 2 ) if ($manual);
pod2usage( -exitstatus => 0, -verbose => 99, -sections => "VERSION") if ($version);
pod2usage( -exitstatus => 0, -verbose => 99, -sections => "NAME|VERSION|DESCRIPTION|AUTHOR|COPYRIGHT|LICENSE|DATE|CITATION") if ($further_information);

=pod

=over 2

=item B<>

=item B<############################################################>

=item B<######## DOMINO Input Data and Pre-processing phase ########>

=item B<############################################################>

=back

=head1 NAME

=over 2

=item B<>

DM_Clean_v1.1.pl 

=back
	
=head1 VERSION

=over 2

=item B<>

DOMINO v1.1 ## Revised 07-11-2018

=back

=head1 SYNOPSIS

=over 2

=item B<>

perl DM_Clean_v1.1.pl 

=item B<###########################>

=item B<### General Information ###>

=item B<###########################>

[-h|--help] [-man] [-v|--version] [-MoreInfo]

=item B<#########################>

=item B<### Mandatory options ###>

=item B<#########################>

[-type_file int_value] [-input_file file] [-b|--barcode_file file] [-o|--outputFolder path]

=item B<######################>

=item B<###	Optional flags ###>

=item B<######################>

## Cleaning

[-l|--cut_Off_ReadLen int_value] [-s|--cut_Off_Quality int_value] [-m|--minLen int_value] [-thr|--threshold_DUST int_value]

## Database search

[-db|--blast_database file] [-no_db_search] [-only_user_db]

## Others 

[-bdiffs int_value]  [-p|--number_cpu int_value] [-only_tag_files] [-TempFiles] 
	
=back
	
=head1 DESCRIPTION

=over 2

=item B<>

=item B<############################################>

=item B<######## GENERAL DOMINO DESCRIPTION ########>

=item B<############################################>

=item B<>

DOMINO is a new software application specifically designed for improving the development of DNA markers in non-model organisms using either NGS data or pre-computed multiple sequence alignments (MSA) in various formats (including MSA from RAD data). 

We used Perl and C++ scripting languages to combine some popular, open-source, bioinformatics tools for NGS data with a set of new developed functions in an integrated bioinformatics pipeline for custom-made marker discovery or selection. 

Using NGS data (raw or pre-processed reads) from typical genome partitioning or low depth sequencing experiments (including or not reference sequences, such as, for example, genome scaffolds), or pre-computed multiple sequence alignments (MSA) of multiple taxa (taxa panel, which may represent different populations or species), DOMINO searches for sequence regions with a minimum level of variation, which can be optionally flanked by stretches of conserved sequences (appropriate for further PCR primers design). 

The DOMINO workflow consists in four main phases: Input Data and Pre-processing, Assembly, Alignment/Mapping and Marker Discovery/Selection. Users can enter to the DOMINO workflow at different points in order to perform specific shortcut runs for marker discovery or selection according to the type of input data.

In a typical full run from raw NGS reads (using marker identification and marker selection modules), the program applies an assembly-based approach; the pipeline is therefore optimized to work with genome partitioning methods in which the size-selected or enriched fragments are long enough to permit a successful assembly and have been fully (or nearly fully) sequenced. For MSA inputs, as for example, restriction-site associated DNA (RAD) variation and related methods (e.g. RADseq, ddRAD, GBS) DOMINO can searh the MSA of each loci (in this case RAD loci) previously obtained with commonly used software packages for highly informative markers (marker selection module); current version accepts the PHYLIP and FASTA format as well as the pyRAD (*.loci) and STACKS (*.fa) output file formats. 

After the maker development phase, DOMINO provides i) the list the genomic regions (with their corresponding coordinates in each contig) exhibiting the minimum levels of nucleotide variation specified by the user among the selected taxa (or covering at least a pre-defined minimum number of them), and the marker MSAs for these taxa (DOMINO marker identification module) or ii) the list of markers exhibiting the same pre-definable features among a set of aligned sequences supplied by the user and previously obtained in other software (DOMINO marker selection module). 

Furthermore, if the taxa panel has been designed in a convenient phylogenetic context and the user asks for highly conserved regions flanking to be included in the developed markers, these markers should be suitable for further surveys of variation in an extended set of phylogenetically related taxa, i.e. focal taxa.

=item B<>

=item B<############################################################>

=item B<######## DOMINO Input Data and Pre-processing phase ########>

=item B<############################################################>

=item B<>

This script is a Perl pipeline for the Quality Control and the accommodation of NGS files and it is the first step in the package DOMINO.

DOMINO accepts input sequence data files in two different formats, the 454 Pyrosequencing Standard Flowgram Format (SFF), and FASTQ format. These input files can contain 454 or Illumina (single or paired-end) raw reads from m taxa (taxa panel). The sequences from each taxon should be properly identified with a specific barcode (aka, tag, MID or index), or loaded in separate files, also appropriately named (see the DOMINO manual in the DOMINO Web site for details). DOMINO is designed to filter low quality, low complexity, contaminant and very short reads using either default or user-specified filtering parameters. Mothur, PRINSEQ, NGS QC toolkit, BLAST, as well as new Perl functions specifically written for DOMINO (DM scripts) are used to perform these tasks.

DOMINO uses Mothur v1.32.0 to extract reads from SFF files and store them in FASTQ format, which are subsequently converted to FASTA and QUAL files. Low quality or very short reads are trimmed, or definitely removed, using NGS QC Toolkit v2.3.1. PRINSEQ v0.20.3 package is used to eliminate low complexity reads using the implemented DUST algorithm. 

Putative contaminant sequences, such as bacterial DNA frequently found in genomic samples, cloning vectors, adaptors, linkers, and high-throughput library preparation primers, can also be removed using a DOMINO function that performs a BLAST search (BLAST v2.2.28) against UniVec database (http://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) and/or against a user-supplied contaminant database (see the DOMINO manual).

=item B<>

=back

=head1 ARGUMENTS

=over 2

=item B<>

=item B<#####################################>

=item B<######## General Information ########>

=item B<#####################################>

=item B<>

=item B<-h>|B<--help>

Brief help message.

=item B<>

=item B<-man>
	
Longer help message.

=item B<>

=item B<-v|--version>

Program version; ignore other arguments.

=item B<>

=item B<-MoreInfo>

Further program information.

=item B<>

=item B<#####################################>

=item B<############### INPUT ###############>

=item B<#####################################>

=item B<>

=item B<-input_file [file]>

Path to the input file(s). Use multiple entries for multiple files.

If the file provided is a unique file containing all the taxa reads accordingly tag using MID tags, (type_file 1,2,4,6) the name for the file is not important as the taxa ids are within the barcodes files. (-b option)

Example 454 unique file: C<-input_file 454_reads.sff>

On the other hand, if providing a unique file per taxa (or a pair) provide files per taxa named as "[xxid-][yyy][_Rn].fastq" or "[yyy][_Rn].fastq". Where:

'[xxid-]' might be present or not. [xx] could be any character (or none). Please avoid using dots (.)
	
'[yyy]' taxon identifier. [Mandatory]

'[_Rn]' If paired-end data, R1 or R2, for the left and the right reads, respectively.

Example Single End: C<-input_file file1 -input_file file2 -input_file file3>

Example Pair End: C<-input_file user_file1_R1.fastq -input_file user_file1_R2.fastq -input_file user_file2_R1.fastq -input_file user_file2_R2.fastq>

The taxon identifier used for these files would be: user_file1, user_file2

Example 2 Pair End: C<-input_file user_file1_id-sp1_R1.fastq -input_file id-sp1_R2.fastq -input_file user_file2_id-sp2_R1.fastq -input_file user_file2_id-sp2_R2.fastq>

The taxon identifier used for these files would be: sp1,sp2

=item B<>

=item B<-type_file [int_value]>

Several types of input files are accepted. Select only one option:

1: Single Standard Flowgram Format (SFF) file, with taxa accordingly tagged.

2: Roche 454 reads FASTQ file with taxa accordingly tagged.

3: Roche 454 reads multiple FASTQ files, each file containing the reads from one taxon. 

4: Single Illumina (single-end) reads FASTQ file with taxa accordingly tagged.

5: Multiple Illumina (single-end) FASTQ files, each file containing the reads from one taxon. 

6: Single Illumina paired-end reads FASTQ file, each pair containing all reads of the different taxa accordingly tagged

7: Multiple Illumina paired-end reads FASTQ files, each pair from the same taxon containing left and right reads, respectively. 

=item B<>

=item B<-b|--barcode_file [file]>

CSV file containing barcodes.

The format accepted for the barcodes file is:
		
		barcode,ACACGACGACT,MID1,Dmelanogaster
		barcode,ACTACGTCTCT,MID2,Dsimulans
		barcode,ACTACGTATCT,MID3,Dyakuba


Example: 

C<-type_file 1 -input_file file.sff -b barcodes.txt>

C<-type_file 2 -input_file 454_file.fastq -b barcodes.txt>

C<-type_file 4 -input_file illumina_file.fastq -b barcodes.txt>

C<-type_file 6 -input_file illumina_pair-end_file_R1.fastq -input_file illumina_pair-end_file_R2.fastq -b barcodes.txt>

=item B<>

=item B<#####################################>

=item B<############### OUTPUT ##############>

=item B<#####################################>

=item B<>

=item B<-o|--outputFolder [path]>

Folder containing the complete DOMINO project.

DOMINO generates a number of subfolders within the main project folder during the process. DOMINO is intended to keep original folder names for further analysis.

=item B<>

=item B<#####################################>

=item B<############ PARAMETERS #############>

=item B<#####################################>

=item B<>

=item B<-bdiffs [int_value]>

Maximun number of errors accepted in the Multiplex Identifier (MID) barcode. [Mothur Default: 1, Max: 2]

=item B<-l|--cut_Off_ReadLen [int_value]>

Cut-off value for the percentage of read length of given quality. [NGSQCToolkit Default: 70]

=item B<-s|--cut_Off_Quality [int_value]>

The cut-off value for PHRED quality score for quality filtering. [NGSQCToolkit Default: 20]

=item B<-m|--minLen [int_value]>

Filter sequences shorter than the given minimum length. [NGSQCToolkit Default: 100]

=item B<-thr|--threshold_DUST [int_value]>

Threshold value for sequence complexity filtering (between 0 and 100). [PRINSEQ Default: 7]

=item B<-p|--number_cpu [int_value]>

Number of threads/cores to be used. [Default: 2]

=item B<>

=item B<#####################################>

=item B<############## OPTIONS ##############>

=item B<#####################################>

=item B<>

=item B<-only_user_db>

Only use user supplied contaminant databases. Use with the -db flag. 

Default DOMINO vector and contaminant databases:

	+ E.coli BL21

	+ Pseudomonas aeruginosa

	+ Saccharomyces cerevisae

	+ Staphyloccocus aureus
	
	+ UniVec Database

=item B<-db|--blast_database [file]>

File(s) with user-provided databases in fasta format (.fasta, .fa, .fna). 

Example: C<-db file1.fasta -db file2.fa -db file3.fna>
	
=item B<-no_db_search>

Use this flag to skip the search for contaminants. The rest of the Quality Control (QC) analysis would be performed.

=item B<-only_tag_files>

Use this flag to specify that Quality Control (QC) is not desired. 

=item B<-TempFiles>

Keep all intermediate files.

=item B<>

=item B<>

=item B<#####################################>

=item B<##### Command Line Examples #########>

=item B<#####################################>

=item B<>

=item B<454 SFF: Unique file>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 1 -o test_SFF_file 
 -input file.sff -b MID_barcodes.txt

=item B<>

=item B<454 Fastq: Unique file>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 2 -o test_454_fastq 
 -input file.fastq -b MID_barcodes.txt

=item B<>

=item B<454: Multiple fastq files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq 

=item B<>

=item B<Illumina: Unique file>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 4 -o test_Illumina 
 -input Illumina_Multiple-taxa.fastq -b MID_barcodes.txt 

=item B<>

=item B<Illumina: Multiple fastq files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 5 -o test_Illumina_multi_files 
 -input sp1_Illumina.fastq -input sp2_Illumina.fastq -input sp3_Illumina.fastq 
 
=item B<>

=item B<Illumina Paired-end Files: Unique file containing all taxa>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 6 -o test_Illumina_pair_end 
 -b MID_barcodes.txt -input reads_R1.fastq -input reads_R2.fastq 

=item B<>

=item B<Illumina Paired-end Files: Multiple Paired-end FASTQ files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 7 -o test_mult_pair_end 
 -input sp1_R1.fastq -input sp1_R2.fastq -input sp2_R1.fastq -input sp2_R2.fastq

=item B<>

=item B<Using User provided databases>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq -db file1.fasta -db file2.fasta 
 -only_user_db

=item B<>

=item B<No Contaminant Search>

 perl [DOMINO_scripts_path]/DM_Clean_v1.1.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq -no_db_search

=item B<>

=back	

=head1 AUTHOR

=over 2

Cristina Frias-Lopez, Jose F. Sanchez-Herrero, Miquel A. Arnedo, Alejandro Sanchez-Gracia and Julio Rozas.
	
Evolutionary Genomics and Bioinformatics Group, Departament de Genètica, Microbiologia i Estadística and Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Av. Diagonal 643, Barcelona 08028, Spain

Availability: 
DOMINO is freely available from www.ub.edu/softevol/domino or via github (https://github.com/molevol-ub/DOMINO)
	
Bugs & Comments

If you encounter bugs that are not listed here, please send a
report to the main authors.
  	
   	J.F.Sanchez: jfsanchezherrero@ub.edu
	A.Sanchez-Gracia: elsanchez@ub.edu
   	J.Rozas: jrozas@ub.edu 


=back

=head1 COPYRIGHT

=over 2

Copyright (C) 2016 Evolutionary Genomics and Bioinformatics Group, Universitat Barcelona.

=back

=head1 LICENSE

=over 2

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=back

=head1 DATE

=over 2

07 - 11 - 2016

=back	

=head1 CITATION

=over 2

Bioinformatics first published online August 16, 2016 
doi:10.1093/bioinformatics/btw534 

=back

=cut
##################################
##	Checking user options		##
##################################
if (!$outFolder || !$user_option_file_type) { DOMINO::dieNicely(); }
if (!$bdiffs) { $bdiffs = 1; }
if (!$read_length_GUI) { $read_length_GUI = 70; }
if (!$minimum_length_GUI) { $minimum_length_GUI = 100; }
if (!$minimum_qual_GUI) { $minimum_qual_GUI = 20; }
if (!$threshold) { $threshold = 7;}
unless ($noOfProcesses) { $noOfProcesses = 2; }

# Generating an output folder where the user decide it
my $random_number = int(rand(100));
if ($outFolder =~ /.*\/$/) { chop $outFolder;}
my $dirname = abs_path($outFolder);
if (!$dirname) { DOMINO::dieNicely(); }
mkdir $dirname, 0755;
my $datestring = strftime "%Y%m%d%H%M", localtime;
my $clean_folder = $dirname."/".$datestring."_DM_clean_data";
my $intermediate_folder = $clean_folder."/".$datestring."_DM_cleaning_intermediate_data";
my $error_log = $dirname."/".$datestring."_DM_Cleaning_ERROR.txt";
my $param_Detail_file = $dirname."/".$datestring."_DM_Cleaning_Parameters.txt";
if (-e $param_Detail_file) { File::Copy::move($param_Detail_file, $param_Detail_file."_old_".$random_number);
} elsif (-e $error_log) { File::Copy::move($error_log, $error_log."_old_".$random_number); }

### Lets start the pipeline
my $start_time = $step_time = time;
print "\n"; 
DOMINO::printHeader("","#"); DOMINO::printHeader(" DOMINO Cleaning Stage ","#"); DOMINO::printHeader("","#"); 
print "\n"; DOMINO::printHeader("","+"); DOMINO::printHeader(" Analysis Started ","+");  DOMINO::printHeader("","+"); 
DOMINO::printDetails("Starting the process: [ ".(localtime)." ]\n\n", $param_Detail_file);
my $localtime = localtime; push (@{ $domino_params{'clean_data'}{'start'} }, $localtime);

## Getting scripts path variable
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
my $BLAST = $scripts_path."NCBI_BLAST/";
my $db_dirname_default = $scripts_path."db_default"; # Default databases provided
my $db_dirname = $dirname."/db";

DOMINO::printHeader(" Input File and Parameter Preprocessing ","#");
DOMINO::printDetails("\n+ Output Directory: ".$dirname." ...OK\n\n", $param_Detail_file); push (@{ $domino_params{'clean_data'}{'folder'} }, $dirname);
if (scalar (@user_files) == 0 || !$user_option_file_type) {
	DOMINO::printError("No input files provided"); DOMINO::dieNicely();
} else {
	## Check files absolute path and name
	for (my $i = 0; $i < scalar @user_files; $i++) {
		my $tmp = abs_path($user_files[$i]);
		push (@file_abs_path, $tmp);
		## Get the name
		my @tmp_file_path = split ("/", $tmp);
		my @tmp_name;
		if ($tmp_file_path[-1] =~ /.*\.fastq/) {
			@tmp_name = split ("\.fastq", $tmp_file_path[-1]);
		} else {
			@tmp_name = split ("\.sff", $tmp_file_path[-1]); 		
		}
		push (@file_names, $tmp_name[0]);		
}}
chdir $dirname;

## Option 1: User barcodes provided
my @MID_species_array;
if ($user_option_file_type == 1 || $user_option_file_type == 2 || $user_option_file_type == 4 || $user_option_file_type == 6 ) {
	if ($user_barcodes_file) {
		my $user_barcodes_path_file = abs_path($user_barcodes_file);
		push (@{$domino_files{"original"}{"barcodes_file"}}, $user_barcodes_path_file);
		unless (-e $user_barcodes_path_file && -r $user_barcodes_path_file) {
			DOMINO::printError("File $user_barcodes_path_file does not exist or it is not readable\n Please provide a readable file");
			&FormatError_barcodes(); DOMINO::dieNicely();
		}
		my $ROCHE_oligos_file = $dirname."/ROCHE.oligos";
		open (ROCHE, ">$ROCHE_oligos_file");	
		open (BARCODE, "<$user_barcodes_path_file");
		while (<BARCODE>) {
			next if /^#/ || /^\s*$/;
			my $line = $_; chomp $line;
			$line =~ s/\t*/\ /g;
			$line =~ s/\s*//g;
			my @barcode_array = split ("\,", $line);
			if (scalar @barcode_array < 2) { &format_error(); DOMINO::dieNicely(); }
			print ROCHE "barcode\t\t$barcode_array[1]\t\t$barcode_array[3]\n";
		} close(BARCODE);
	} else {
		DOMINO::printError("\n\nERROR: No barcodes file provided...\n"); &FormatError_barcodes(); DOMINO::dieNicely();
}}

###################
##	Checking files
###################
my $file_type = $input_options{$user_option_file_type};
if ($input_options{$user_option_file_type}) {
	DOMINO::printDetails("+ Type of file(s): Option: $user_option_file_type -- $file_type ...OK\n", $param_Detail_file); 
	push (@{ $domino_params{'clean_data'}{'option'} }, "$user_option_file_type"."--".$file_type);
} else { 
	DOMINO::printError("\n\nERROR: Wrong type of file provided\nPlease provide a valid type of file:\n");
	DOMINO::printInput_type(); DOMINO::dieNicely();
}
DOMINO::printDetails("+ Checking file(s):\n", $param_Detail_file);

#### Type file: 1 
if ($file_type eq "454_sff") {
	if (scalar (@file_abs_path) != 1) {
		DOMINO::printError("Please provide a unique SFF file"); DOMINO::dieNicely();
	} elsif ($file_abs_path[0] !~ /(.*)\.sff/) {
		DOMINO::printError("Wrong 454 file provided.\nPlease provide a binary SFF extension file or provide a different type of file"); DOMINO::dieNicely();
	} elsif (-e -r -s $file_abs_path[0]) {
		my $id = $file_names[0];
		DOMINO::printDetails("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = SFF ...OK\n\t\tIt also contains an identifier ($id) in the name for later analysis...OK\n\n", $param_Detail_file);	
		push (@{$domino_files{"original"}{'SFF'}}, $file_abs_path[0]); 

#### Type file: 2
}} elsif ($file_type eq "454_fastq") {
	if (scalar (@file_abs_path) != 1) {
		DOMINO::printError("Please provide a unique 454 Roche FASTQ file containing all the taxa read accordingly tagged with MID sequences"); DOMINO::dieNicely();
	} elsif (scalar (@file_abs_path) == 1) {
		my $format = DOMINO::check_file_format($file_abs_path[0]);
		if ($file_abs_path[0] =~ /.*\.f.*q$/) {
			if ($format =~ /fastq/) {
				DOMINO::printDetails("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n", $param_Detail_file);		
				push (@{$domino_files{"original"}{'454_fastq'}}, $file_abs_path[0]); 

			} else { 
				DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[0]... Format: $format"); DOMINO::dieNicely();
		
#### Type file: 3
}}}} elsif ($file_type eq "454_multiple_fastq") {
	if (scalar (@file_abs_path) == 1) {
		DOMINO::printError("Please provide multiple 454 Roche FASTQ file containing each of the taxa reads"); DOMINO::printFormat_message; DOMINO::dieNicely();
	} else {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);

#### Type file: 4
}}} elsif ($file_type eq "Illumina") {
	if (scalar (@file_abs_path) != 1) {
		DOMINO::printError("Please provide a unique Illumina FASTQ file containing all the taxa read accordingly tagged with MID sequences"); DOMINO::dieNicely();
	} elsif (scalar (@file_abs_path) == 1) {
		my $format = DOMINO::check_file_format($file_abs_path[0]);
		if ($file_abs_path[0] =~ /.*\.f.*q$/) {
			if ($format =~ /fastq/) {
				DOMINO::printDetails("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n", $param_Detail_file);		
				push (@{$domino_files{"original"}{'illu_FASTQ'}}, $file_abs_path[0]); 
			} else { 
				DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[0]... Format: $format"); DOMINO::dieNicely();

#### Type file: 5
}}}} elsif ($file_type eq "Illumina_multiple_fastq") {  
	if (scalar (@file_abs_path) == 1) {
		DOMINO::printError("Please provide multiple Illumina FASTQ file containing containing each of the taxa reads"); DOMINO::printFormat_message; DOMINO::dieNicely();
	} elsif (scalar (@file_abs_path) != 1) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);

#### Type file: 6
}}} elsif ($file_type eq "Illumina_pair_end") {
	if (scalar (@file_abs_path) != 2) {
		DOMINO::printError("Please provide Illumina Paired-End FASTQ files. Provide the left (_R1) and then the right (_R2) file containing all the taxa read accordingly tagged with MID sequences\nPlease tagged left read file with xxx_R1.fastq and right read file as xxx_R2.fastq\n"); DOMINO::dieNicely();
	} elsif (scalar (@file_abs_path) == 2) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			my $format = DOMINO::check_file_format($file_abs_path[$i]);
			if ($file_abs_path[$i] =~ /.*\_R(\d+)\.f.*q$/) { ## .*_R[*].fastq
				my $pair_int = $1; my $pair;
				if ($pair_int == 1) { $pair = "Left reads file";	
				} elsif ($pair_int == 2) { $pair = "Right reads file";	
				} else { DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[$i]...\nFormat: $format\nPair: $pair_int\n"); DOMINO::dieNicely(); }
				if ($format =~ /fastq/) {
					DOMINO::printDetails("\t$file_abs_path[$i]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n\t\tPaired end file = ".$pair." ...OK\n", $param_Detail_file);
					push (@{$domino_files{"original"}{'illu_PE'}}, $file_abs_path[$i]);
				} else { DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[$i]...\nFormat: $format\nPair: $pair_int\n"); DOMINO::dieNicely();

#### Type file: 7
}}}}} elsif ($file_type eq "Illumina_pair_end_multiple_fastq") {  
	if (scalar (@file_abs_path) == 2 ) {
		DOMINO::printError("Please provide Illumina Paired-End FASTQ files for each taxa. Provide the left (_R1) and then the right (_R2) file containing each taxa reads"); DOMINO::printFormat_message; DOMINO::dieNicely();
	} elsif (scalar (@file_abs_path) != 2) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);
}}}

## Check if paired end files provided, all pairs have been correctly provided
if ($file_type eq "Illumina_pair_end_multiple_fastq") {
	my %species_names = %{ $domino_files{'original'} };
	foreach my $keys (keys %species_names) {
		my $sum_pairs; my $files;
		##@{ $domino_files{"original"}{$name} }
		my @array = @{ $species_names{$keys} };
		if (scalar @array != 2) {
			DOMINO::printError("Wrong pair of FASTQ files provided. Please check pair of files corresponding to $keys:\n$files\n"); DOMINO::dieNicely();
		} else { &debugger_print("OK: $keys\n\n"); }
}}

&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

##################################################
##	Printing user options 			##
##################################################
if ($onlyTagging_files) {
	DOMINO::printDetails("+ Cleaning process would be skipped...OK\n", $param_Detail_file);
	DOMINO::printDetails("+ Only extracting and tagging of the reads would be done...OK\n", $param_Detail_file);	
	push (@{ $domino_params{'clean_data'}{'onlyTagging_files'} }, "1");
} else {
	## Print information of the input parameters
	DOMINO::printDetails("+ Number of Processor: ".$noOfProcesses." ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'clean_data'}{'CPU'} }, $noOfProcesses);

	DOMINO::printDetails("+ Number of mismatches allowed in the barcode sequence: $bdiffs ...OK\n", $param_Detail_file);	
	push (@{ $domino_params{'clean_data'}{'mismatches'} }, $bdiffs);

	DOMINO::printDetails("+ Minimum read length: $minimum_length_GUI pb ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'clean_data'}{'read_length'} }, $minimum_length_GUI);

	DOMINO::printDetails("+ Minimum QUAL: $minimum_qual_GUI ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'clean_data'}{'QUAL'} }, $minimum_qual_GUI);

	DOMINO::printDetails("+ Minimum length of a read satisfying QUAL cutoff: $read_length_GUI % ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'clean_data'}{'QUAL_cutoff'} }, $read_length_GUI);

	DOMINO::printDetails("+ Threshold for complexity/entropy of a read: $threshold ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'clean_data'}{'entropy'} }, $threshold);

	unless ($avoidDelTMPfiles) {
		DOMINO::printDetails("+ Deleting of temporary files would be done ...OK\n", $param_Detail_file);
	} else { 
		DOMINO::printDetails("+ Deleting temporary files would be avoid ...OK\n", $param_Detail_file);
	}

	## Check databases to use in BLAST contaminations search
	unless ($skipping_BLAST) {
		DOMINO::printDetails("+ Checking Databases to use during BLAST search of contaminants:\n\tA folder named as 'db' would be generated to place the databases we are going to use\n", $param_Detail_file);
		mkdir $db_dirname, 0755; chdir $db_dirname;
		if ($only_user_db) {
			DOMINO::printDetails("\t\tNo default databases would be used and only user provided databases would be used instead\n", $param_Detail_file);
			push (@{ $domino_params{'clean_data'}{'only_user_DB'} }, "1");
		} else {
			DOMINO::printDetails("\tDefault databases provided in the package would be used\n\tCopying...\n", $param_Detail_file);
			push (@{ $domino_params{'clean_data'}{'default_DB'} }, "1");
			&download_db(); 
		}

		## Check if User provide any database
		if (@user_blast_db_tmp) {
			DOMINO::printDetails("\t+ Additional user provided databases would be used:\n", $param_Detail_file);
			push (@{ $domino_params{'clean_data'}{'additional'} }, "1");
			for (my $i = 0; $i < scalar @user_blast_db_tmp; $i++) {
				print "\t$user_blast_db_tmp[$i]\n";
				push (@{ $domino_params{'clean_data'}{'additional_db'} }, $user_blast_db_tmp[$i]);
				my $tmp_abs_name = abs_path($user_blast_db_tmp[$i]);
				my $format_returned = DOMINO::check_file_format($tmp_abs_name);
				if ($format_returned !~ /fasta/) { 
					DOMINO::printError("Please, enter a valid FASTA file, $user_blast_db[$i] is not valid"); DOMINO::dieNicely();
				} else {
					push (@user_blast_db, $tmp_abs_name);
		}}}
		DOMINO::printDetails("\n\t+ Database(s):\n", $param_Detail_file); &looking4db();
	} else {
		DOMINO::printDetails("+ No BLAST search of contaminants ...OK\n", $param_Detail_file);
		push (@{ $domino_params{'clean_data'}{'NO_search'} }, "1");

}}
## Print info about where to print info and error
DOMINO::printDetails("\n+ Parameters details would be print into file: ".$clean_folder."/".$datestring."_DM_Cleaning_Parameters.txt...\n", $param_Detail_file);
DOMINO::printDetails("+ Errors occurred during the process would be print into file: ".$clean_folder."/".$datestring."_DM_Cleaning_ERROR.txt...\n\n", $param_Detail_file);
chdir $dirname; print "\n"; &time_log(); print "\n";

###########################
##		Start the analysis 
###########################
if ($user_option_file_type == 1 || $user_option_file_type == 2 || $user_option_file_type == 4 || $user_option_file_type == 6) {
	
	##################################################################################################
	## If the type of file is SFF, 454_FASTQ, Illumina and Illumina paired-end and it is a unique   ##
	## file containing information of the different taxa, tagged with MID tags we would extract  	##
	## and classify 										##
	##################################################################################################
	
	print "\n"; DOMINO::printHeader(" Extracting taxa sequences from the file provided ","%"); 
	my $file_tmp_barcodes = $domino_files{"original"}{"barcodes_file"}[0];
	push (@{ $domino_params{'clean_data'}{'barcodes'} },$file_tmp_barcodes);

	DOMINO::printDetails("+ Barcodes file provided: $file_tmp_barcodes ...OK\n", $param_Detail_file);
	
	for (my $i=0; $i < scalar @file_abs_path; $i++) {
		if ($file_type eq "454_sff") { 
			# Extracting the sff file given
			print "+ Extracting information from the SFF file provided as input\n\t-> This might take a while for large files\n- Entering to MOTHUR executable";
			&mothur_sffinfo($file_abs_path[$i], $file_names[$i]); 	print "\n"; &time_log(); print "\n";
			print "\n- Exiting MOTHUR executable\n- Job done: sff file has been extracted and sequences has been trimmed according to their tag\n";
			&splitting_fastq($file_names[$i], $i, $dirname); 	print "\n"; &time_log(); print "\n";
		} else {
			&extracting_fastq($file_abs_path[$i], $file_names[$i]); 
			&mothur_trim_fastq($file_names[$i]); &splitting_fastq($file_names[$i], $i, $dirname);
	}}
	&time_log(); &debugger_print("domino_files"); &debugger_print("Ref", \%domino_files);
	my $array_ref = DOMINO::readDir($dirname);
	
	foreach my $taxa (keys %domino_files) {
		next if ($taxa eq 'original');
		my $dump_file = $domino_files{$taxa}{"DIR"}[0]."/dumper_extracted.txt";
		open (DUMP_IN, "$dump_file");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_IN);
	}
} 
&debugger_print("domino_files"); &debugger_print("Ref", \%domino_files); 
print "\n"; &time_log(); print "\n";

######################################################################################################
## If user only wants to tag files accordingly to the names input and do no cleaning step at all.	##
######################################################################################################
if ($onlyTagging_files) { 
	if ($file_type eq "454_multiple_fastq" || $file_type eq "Illumina_multiple_fastq" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		my $pm =  new Parallel::ForkManager($noOfProcesses);
		$pm->run_on_finish( 
			sub { my ($pid, $exit_code, $ident) = @_; 
			print "\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
		} );
		$pm->run_on_start( sub { my ($pid,$ident)=@_; print "** Rename ids for $ident with PID $pid\n\n"; } );
		my $int=0;
		foreach my $keys (keys %{ $domino_files{'original'} }) {
    		$int++;
			my $pid = $pm->start($keys) and next;
    		print "\nSending child process for renaming fastq for $keys\n\n";
			my %domino_files_threads_OTG;
			my $dir = $dirname."/".$keys;
			unless (-d $dir) { mkdir $dir, 0755; push (@{ $domino_files_threads_OTG{$keys}{"DIR"} }, $dir); }
			my @array = @{ $domino_files{'original'}{$keys}};
			for (my $j=0; $j < scalar @array; $j++) {	
				my $copy_name = &rename_seqs($domino_files{'original'}{$keys}[$j], $keys, "fastq");
				File::Copy::copy($domino_files{'original'}{$keys}[$j], $copy_name);
				push (@{ $domino_files_threads_OTG{$keys}{"FINAL_fq"} }, $copy_name);			
			}
			my $domino_files_threads_OTG_text = $dir."/dumper-hash_OTG.txt";
			DOMINO::printDump(\%domino_files_threads_OTG, $domino_files_threads_OTG_text);			
    		$pm->finish($int); # pass an exit code to finish
		}
		$pm->wait_all_children;
		print "\n** All child processes have finished...\n\n";		

		## Check each taxa dump file conainting file info
		&debugger_print("Retrieve info from files");
		## read dirs in hash domino_files and get dumper file
		my %tmp_hash = %domino_files;
		foreach my $taxa (keys %{ $tmp_hash{'original'} }) {
			my $dump_file = $dirname."/".$taxa."/dumper-hash_OTG.txt";
			open (DUMP_IN, "$dump_file");
			while (<DUMP_IN>) {
				my $line = $_; chomp $line;
				my @array = split("\t", $line);
				&debugger_print($line);
				push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
			} close (DUMP_IN);
	}} # else unique file?
	&debugger_print("domino_files"); &debugger_print("Ref", \%domino_files);	
	&sort_files_and_folders();
	print "\n+ Exiting the script...\n"; DOMINO::finish_time_stamp($start_time); exit();	 
		
} else {

	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";
	## Copy files and get ready for QC
	if ($file_type eq "454_multiple_fastq" || $file_type eq "Illumina_multiple_fastq" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		my $pm_copy =  new Parallel::ForkManager($noOfProcesses);
		$pm_copy->run_on_finish( 
			sub { my ($pid, $exit_code, $ident) = @_; 
			print "\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
		} );
		$pm_copy->run_on_start( sub { my ($pid,$ident)=@_; print "** Copy file for $ident with PID $pid\n\n"; } );
		my $int=0;
		foreach my $keys (keys %{ $domino_files{'original'} }) {
			$int++;
			my $pid = $pm_copy->start($keys) and next;
    		print "\nSending child process for copying fastq for $keys\n\n";
			my %domino_files_cp;
			#$species_names{$keys} = $keys;
			my $dir = $dirname."/".$keys;
			unless (-d $dir) { mkdir $dir, 0755; push (@{ $domino_files_cp{$keys}{"DIR"} }, $dir); }			
			for (my $i=0; $i < scalar @{ $domino_files{'original'}{$keys} }; $i++) {
				my @path_name = split("/", $domino_files{'original'}{$keys}[$i]);
				my $copy_name;
				if ($path_name[-1] =~ /.*\_id\-.*/) { $copy_name = $dir."/".$path_name[-1];
				} else { $copy_name = $dir."/reads_id-".$path_name[-1];
				}
				File::Copy::copy($domino_files{'original'}{$keys}[$i], $copy_name);
				push (@{ $domino_files_cp{$keys}{"EXTRACTED"} }, $copy_name);
			}
			my $domino_copy_files_threads_text = $dirname."/dumper-hash_cp.txt";
			DOMINO::printDump(\%domino_files_cp, $domino_copy_files_threads_text);
	 		$pm_copy->finish($int); # pass an exit code to finish
		}
		$pm_copy->wait_all_children; print "\n** All child processes have finished...\n\n";

		## Check each taxa dump file conainting file info
		&debugger_print("Retrieve info from dumper-hash_cp.txt");
		my $dump_file = $dirname."/dumper-hash_cp.txt";
		open (DUMP_IN, "$dump_file");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_IN);
	}
	print "\n"; &time_log(); print "\n";
}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";

##############################
# Concatenating databases
my $merge_fasta_database = $db_dirname."/merge_seqs.fasta";
open (DB, ">$merge_fasta_database");
unless ($skipping_BLAST) {
	for (my $i = 0; $i < @blast_db; $i++) {	
		my $file = $blast_db[$i];
		push (@{ $domino_files{'original'}{'db'} }, $file);
		push (@{ $domino_params{'clean_data'}{'db'} }, $file);
		open(FILE, $file) || die "Could not open the $file ...\n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
    		my ($titleline, $sequence) = split(/\n/,$_,2);
    		next unless ($sequence && $titleline);
    		print DB ">".$titleline."\n".$sequence;
    	} close(FILE); $/ = "\n";
}} close (DB);

##############################
## QC analysis for each taxa
DOMINO::printHeader(" Quality Control Analysis ","#");
## Send threads for each taxa
my $int_taxa = 0;
my $pm =  new Parallel::ForkManager($noOfProcesses); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
$pm->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
	print "\n\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
} );
$pm->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** QC Analysis started for $ident with PID $pid\n\n"; } );
foreach my $keys (sort keys %domino_files) {
	next if ($keys eq 'original');	
	$int_taxa++;
	my $pid = $pm->start($keys) and next; print "\nSending child process for renaming fastq for $keys\n\n";
	
	my %domino_files_threads_QC;
	my $taxa_dir = $domino_files{$keys}{"DIR"}[0];
	chdir $taxa_dir; print "Changing dir to $taxa_dir\n";

	##################################################################
	## 	Quality control step 1.1 -- PRINSEQ Dust algorithm	##
	##################################################################
	print "\n+ Quality control 1.1: low complexity reads: DUST algorithm: $keys\n";
	my $dust_clean_reads; my $dust_reads; 
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {		
		($dust_clean_reads, $dust_reads) = &prinseq_dust($domino_files{$keys}{'EXTRACTED'}[0], $domino_files{$keys}{'EXTRACTED'}[1], $threshold, $taxa_dir);
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads."_1.fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads."_2.fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads."_1.fastq");
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads."_2.fastq");
	} else {
		($dust_clean_reads, $dust_reads) = &prinseq_dust($domino_files{$keys}{'EXTRACTED'}[0], $domino_files{$keys}{'EXTRACTED'}[0], $threshold, $taxa_dir);
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads.".fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads.".fastq");
	}
	print "\n"; &time_log(); print "\n";
	&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

	##################################################
	## Quality control step 1.2 -- NGS_QC_toolkit	##
	##################################################
	print "- Quality control 1.2: short reads and poor quality: \n";
	print "+ Cleaning poor quality sequences and short sequences. Please check the stats file for results or the html report\n\n";
	my $NGS_QC_call = "perl ".$scripts_path."NGSQCToolkit_v2.3.1/";
	if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
		$NGS_QC_call .= "454QC.pl -i ".$dust_clean_reads.".fasta ".$dust_clean_reads.".qual N";
		$NGS_QC_call .= " -m ".$minimum_length_GUI; # Filter sequences shorter than
		push (@{ $domino_files_threads_QC{$keys}{"NGS_QC"} }, $dust_clean_reads.".fasta_filtered");
		push (@{ $domino_files_threads_QC{$keys}{"NGS_QC"} }, $dust_clean_reads.".qual_filtered");
	} else {
		$NGS_QC_call .= "IlluQC.pl ";		
		if ($file_type eq "Illumina_multiple_fastq" || $file_type eq "Illumina") { 
			$NGS_QC_call .= "-se ".$dust_clean_reads.".fastq N A"; 
			push (@{ $domino_files_threads_QC{$keys}{"NGS_QC"} }, $dust_clean_reads.".fastq_filtered");
		} else { ## Pair end
			$NGS_QC_call .= "-pe ".$dust_clean_reads."_1.fastq ".$dust_clean_reads."_2.fastq N A";
			push (@{ $domino_files_threads_QC{$keys}{"NGS_QC"} }, $dust_clean_reads."_1.fastq_filtered");
			push (@{ $domino_files_threads_QC{$keys}{"NGS_QC"} }, $dust_clean_reads."_2.fastq_filtered");
	}}
	$NGS_QC_call .= " -t 2 -o $taxa_dir -p 1 -l ".$read_length_GUI." -s ".$minimum_qual_GUI." 2> $error_log";
		# Number of processor to use (ONLY 1 because this would be already parallelized)
		# Setting the cutoff percentage of read length that should be of given quality [Default 70]
		# Setting the minimum cutoff for PHRED quality, default 20 
	
	## Sending command
	&debugger_print("NGS QC Toolkit command: $NGS_QC_call\n");
	my $NGSQC_result = system($NGS_QC_call);
	if ($NGSQC_result != 0) { DOMINO::printError("Cleaning step failed when calling NGSQC for file $dust_clean_reads..."); exit(); }	
	print "\n"; &time_log();
	&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

	##########################################################################
	## Quality control step 2 -- BLAST contaminant and vector search	##
	##########################################################################
	if ($skipping_BLAST) {
		#If user does not want to search against any database, we would skip the whole step
		if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
			my @array_files = @{ $domino_files_threads_QC{$keys}{"NGS_QC"} };
			for (my $i=0; $i < scalar @array_files; $i++) {
				my $final_name;
				my @path_name = split("/", $domino_files_threads_QC{$keys}{"NGS_QC"}[$i]);
				if ($path_name[-1] =~ /.*\_id\-(.*)\_R\d+\.clean\_dust\_(\d+).*/) {
					for (my $j=0; $j < scalar $#path_name; $j++) {
						$final_name .= $path_name[$j]."/";
					} $final_name .= "QC-filtered_id-".$1."_R".$2.".fastq";
				}
				push (@{ $domino_files_threads_QC{$keys}{"FINAL_fq"} }, $final_name);
				File::Copy::move( $domino_files_threads_QC{$keys}{"NGS_QC"}[$i], $domino_files_threads_QC{$keys}{"FINAL_fq"}[$i]);
			}
		} else {			
			my $final_name = $taxa_dir."/QC-filtered_id-".$keys.".fastq";
			my $final_file = DOMINO::qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"NGS_QC"}[0], $domino_files_threads_QC{$keys}{"NGS_QC"}[1], $keys, $final_name, $file_type);
			push (@{ $domino_files_threads_QC{$keys}{"FINAL_fq"} }, $final_file);
		}			
		DOMINO::printHeader(" Cleaning Step 2: BLAST and Vector screen ","#"); 
		print "+ Skipping the BLAST step for searching against databases for putative contaminants\n";
	} else {
		DOMINO::printHeader(" Cleaning Step 2: BLAST and Vector screen ","#"); 
		print "+ Eliminating that reads that could be potential contaminations\n+ Default databases would be used plus the one(s) provided by the user (if any)\n\n\n";
		DOMINO::printHeader(" Database Generation ","%"); print "+ Generation of Databases: \n";
		print "Please be patient, it could take a while for large files...\n\n";
		my @reads_db;
		if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
			for (my $i=0; $i < 2; $i++) {
				my @name = split("\.fastq_filtered", $domino_files_threads_QC{$keys}{"NGS_QC"}[$i]);
				&extracting_fastq($domino_files_threads_QC{$keys}{"NGS_QC"}[$i], $name[0], "YES");
				push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $name[0].".filtered.fasta");
				push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $name[0].".filtered.qual");
				my ($blast_DB, $blast_DB_message) = DOMINO::makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$i], $BLAST, $error_log);
				&debugger_print($blast_DB_message); push (@reads_db, $blast_DB);
			}		
		} elsif ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
			push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $domino_files_threads_QC{$keys}{"NGS_QC"}[0]);
			push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $domino_files_threads_QC{$keys}{"NGS_QC"}[1]);
			my ($blast_DB, $blast_DB_message) = DOMINO::makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0], $BLAST, $error_log);
			&debugger_print($blast_DB_message); push (@reads_db, $blast_DB);

		} else {
			my @name = split("\.fastq_filtered", $domino_files_threads_QC{$keys}{"NGS_QC"}[0]);
			&extracting_fastq($domino_files_threads_QC{$keys}{"NGS_QC"}[0], $name[0], "YES");
			push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $name[0].".filtered.fasta");
			push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $name[0].".filtered.qual");
			my ($blast_DB, $blast_DB_message) = DOMINO::makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0], $BLAST, $error_log);
			&debugger_print($blast_DB_message); push (@reads_db, $blast_DB);
		}
		print "\n"; &time_log(); print "\n";
		&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

		# Blasting our reads against the databases generated
		print "\n"; DOMINO::printHeader(" Database Search: BLAST ","%"); 
		my $out_file = $taxa_dir."/blast_search.txt";
		my $db; for (my $j=0; $j < scalar @reads_db; $j++) { $db .= $reads_db[$j]." "; } chop $db;
		my ($blastn, $blastn_message) = DOMINO::blastn($merge_fasta_database, $db, $out_file, $BLAST);
		&debugger_print($blastn_message);
		if ($blastn != 0) { DOMINO::printError("BLASTN failed...\n"); exit(); } 
		print "\n"; &time_log(); print "\n"; print "+ Parsing the BLAST results for $keys...\n"; 

		## Check BLAST result
		my $no_contamination = 0;
		if (-s $out_file) {
		
			##############################################################################################
			## CAUTION: bear in mind, it could be two possible contamination sources: 					##	
			##			- Bacteria/human/other: would be nearly the whole read with few mismatches		##
			##			  against the reference database sequence										##
			##			- Univec contamination: these are primers, adaptors and some cloning vectors	##
			##			  used in wet labs. Not used in our experiment but could be present some small	##	
			##		      traces of other experiments.	 												##
			##	 Instead of using an evalue threshold we would use a similarity of 100% as long as 		##
			##	 evalue relies on query length, database size,....										##	
			##	 If we have a contamination read, the possibility of mismatch would be nearly 0.0		##
			##############################################################################################
		
			my $align_len_read = '90'; my $similarity_thr = '97'; my $var;
			my $message;
			if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
				($var, $message) = DOMINO::readFASTA_hashLength($reads_db[0].".fasta_filtered");
			} else { 
				($var, $message) = DOMINO::readFASTA_hashLength($reads_db[0].".fasta"); 
			}
			&debugger_print($message);
			my %hash_length = %$var; my $fixed_length;
			if ($hash_length{"UNDEF"}) { $fixed_length = $hash_length{"UNDEF"}; }			
			print "- Reading the BLAST result file $out_file ...\n";
			$/ = "\n";
			my $contamination_ReadIds = $taxa_dir."/contaminant_ReadIds_".$keys.".txt";
			my $contaminants = $taxa_dir."/contaminants_Identified_".$keys.".txt";
			
			open (READS, ">$contamination_ReadIds"); open (CONTAMINANTS, ">$contaminants");	
			open(BLAST_result, "<$out_file") || die "Could not open the $out_file file.\n";
			while(<BLAST_result>){
				my $line = $_;
				chomp $line;
				$line =~ s/\s+/\t/g; $line =~ s/\t+/\t/g;	
				my @array = split ("\t", $line);
				## gnl|uv|NGB00361.1:1-92	HWI-ST993:479:C3DBYACXX:4:1101:10015:3765	100.00		33		0	0	60	92	101	69	2e-11	62.1
				## query_id					## subjectID								#perc		#length							#evalue
				## array[0]					array[1]									array[2] 	array[3]						#array[10]
				if ($array[0] =~ /gnl|uv.*/) { ## Hits similar to vectors
					if ($array[2] > $similarity_thr) { 
						print READS $array[1]."\n"; print CONTAMINANTS $array[0]."\n";
					}
				} else { ## Rest of the reads
					if ($array[2] > $similarity_thr) { 
						if ($fixed_length) {
							my $match_align = $array[3];
							my $percent = (($match_align*100)/$fixed_length);
							if ($percent > $align_len_read) { 
								print READS $array[1]."\n"; print CONTAMINANTS $array[0]."\n";
						}} elsif ($hash_length{$array[1]}) {
							my $match_align = $array[3];
							my $percent = (($match_align*100)/$hash_length{$array[1]});
							if ($percent > $align_len_read) { 
								print READS $array[1]."\n"; print CONTAMINANTS $array[0]."\n";
						}} else { next;				
			}}}} 
			close(BLAST_result); close (READS); close(CONTAMINANTS);
			
			## TODO: Print contaminants ($array[0]) to a file (cat | uniq ) to show user not the reads but the contaminants found
			
			## Discard reads
			if (-s $contamination_ReadIds) {
				## Check this paired-end
				my @array = @{$domino_files_threads_QC{$keys}{"FASTA4BLAST"}};
				for (my $j=0; $j < scalar @array; $j++) {
					DOMINO::mothur_remove_seqs($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$j], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[$j], $taxa_dir, $contamination_ReadIds, $mothur_path);
					my @name;
					if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
						@name = split("\.fasta_filtered", $domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$j]);
						push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST_filtered"} }, $name[0].".pick.fasta_filtered");
						push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST_filtered"} }, $name[0].".pick.qual_filtered");
					} else {
						@name = split("\.filtered\.fasta", $domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$j]);
						push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST_filtered"} }, $name[0].".filtered.pick.fasta");
						push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST_filtered"} }, $name[0].".filtered.pick.qual");
					}
					my $final_name;
					my @path_name = split("/", $domino_files_threads_QC{$keys}{"FASTA4BLAST_filtered"}[$j]);
					if (scalar @array > 1) {
						if ($path_name[-1] =~ /.*\_id\-(.*)\_R\d+\.clean\_dust\_(\d+).*/) {
							for (my $j=0; $j < scalar $#path_name; $j++) {
								$final_name .= $path_name[$j]."/";
							} $final_name .= "QC-filtered_id-".$1."_R".$2.".fastq";
					}} else {
						if ($path_name[-1] =~ /.*\_id\-(.*)\.clean\_dust.*/) {
							for (my $j=0; $j < scalar $#path_name; $j++) {
								$final_name .= $path_name[$j]."/";
							} $final_name .= "QC-filtered_id-".$1.".fastq";
					}}
					my $fastq_filtered = DOMINO::qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST_filtered"}[$j], $domino_files_threads_QC{$keys}{"QUAL4BLAST_filtered"}[$j], $keys, $final_name, $file_type);
					push (@{ $domino_files_threads_QC{$keys}{"FINAL_fq"} }, $fastq_filtered);
				} print "\n+ BLAST search and parsing of results for $keys done...\n";		
				
				## Get contaminants for each taxa
				push(@{$domino_files_threads_QC{$keys}{"CONTAMINANTS"} }, $contaminants);
				push(@{$domino_files_threads_QC{$keys}{"reads_CONTAMINANTS"} }, $contamination_ReadIds);
				
			} else { $no_contamination = 1; }
		} else { $no_contamination = 1; }
		if ($no_contamination == 1) {
			print "\n+ No contaminants were found for $keys...\n";
			my @array_files = @{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} };
			if (scalar @array_files > 1) {
				for (my $i=0; $i < scalar @array_files; $i++) {
					my $final_name;
					my @path_name = split("/", $domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$i]);
					if ($path_name[-1] =~ /.*\_id\-(.*)\_R\d+\.clean\_dust\_(\d+).*/) {
						for (my $j=0; $j < scalar $#path_name; $j++) {
							$final_name .= $path_name[$j]."/";
						} $final_name .= "QC-filtered_id-".$1."_R".$2.".fastq";
					}
					my $fastq_filtered = DOMINO::qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$i], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[$i], $keys, $final_name, $file_type);
					push (@{ $domino_files_threads_QC{$keys}{"FINAL_fq"} }, $fastq_filtered);
			}} else {
				my $final_name;
				my @path_name = split("/", $domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0]);
				if ($path_name[-1] =~ /.*\_id\-(.*)\.clean\_dust.*/) {
					for (my $j=0; $j < scalar $#path_name; $j++) {
						$final_name .= $path_name[$j]."/";
					} $final_name .= "QC-filtered_id-".$1.".fastq";
				}
				my $fastq_filtered = DOMINO::qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[0], $keys, $final_name, $file_type);
				push (@{ $domino_files_threads_QC{$keys}{"FINAL_fq"} }, $fastq_filtered);
	}}} 
	print "\n"; &time_log(); print "\n";
	
	&debugger_print("DOMINO Files threads QC"); &debugger_print("Ref", \%domino_files_threads_QC); print "\n";
	
	my $domino_files_threads_QC_text = $taxa_dir."/dumper-hash_QC.txt";
	DOMINO::printDump(\%domino_files_threads_QC, $domino_files_threads_QC_text);

	########################################
	## 	Getting files ready for finishing	
	########################################
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") { ##adjust pairs?
	} else { ## What is missing?? 
	}
	$pm->finish($int_taxa); # pass an exit code to finish
}
$pm->wait_all_children; print "\n** All QC child processes have finished...\n\n";		

## Check each taxa dump file conainting file info
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";

&debugger_print("Retrieve info from files");
foreach my $taxa (keys %domino_files) {
	next if ($taxa eq 'original');
	if ($domino_files{$taxa}{'DIR'}) {
		my $dump_file = $domino_files{$taxa}{'DIR'}[0]."/dumper-hash_QC.txt";
		open (DUMP_IN, "$dump_file");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
		}
		close (DUMP_IN);
}}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";

##########################################################################
## 		Deleting temporary files and sorting the output files			##
##########################################################################
chdir $dirname; 
print "\n\n"; DOMINO::printHeader(" Deleting Temporary Files ","#");
print "+ Deleting or renaming some files...\n"; 
&sort_files_and_folders(); print "\n";

##################################################
## 	Generating output folders					##
##################################################
#print "DOMINO Files\n"; print Dumper \%domino_files;
my $dump_info_DOMINO = $clean_folder."/DOMINO_dump_information.txt"; DOMINO::printDump(\%domino_files, $dump_info_DOMINO);
my $dump_param_DOMINO = $clean_folder."/DOMINO_dump_param.txt"; DOMINO::printDump(\%domino_params, $dump_param_DOMINO);
print "\n Job done succesfully, exiting the script\n"; 
DOMINO::finish_time_stamp($start_time); print "\n";

exit(0);


#######################################
#######	     SUBROUTINES	#######
#######################################
sub check_file {
	## Populate %species_names with the ids of each file
	my $file_to_check = $_[0];

	my $format_returned;
	if (-e -r -s $file_to_check) { 
		my $name; my $pair_int;
		if ($file_to_check =~ /.*id\-.*/) {
			if ($file_to_check =~ /.*id\-(.*)\_R(\d+)\.f.*q$/) { ## [xxx]id-[yyy](_R[*]).fastq
				$name = $1; $pair_int = $2;
			} elsif ($file_to_check =~ /.*id\-(.*)\.f.*q$/) {			
				$name = $1;
		}} else {
			my $tmp_name;
			if ($file_to_check =~ /(.*)\_R(\d+)\.f.*q$/) { ## [xxx]id-[yyy](_R[*]).fastq
				$tmp_name = $1; $pair_int = $2;
			} elsif ($file_to_check =~ /(.*)\.f.*q$/) {			
				$tmp_name = $1;
			}
			my @tmp_file_path = split ("/", $tmp_name);
			$name = $tmp_file_path[-1];		
		}
		$format_returned = DOMINO::check_file_format($file_to_check);
		my $pair;
		if ($pair_int) {
			if ($pair_int == 1) { $pair = "Left reads file";	
			} elsif ($pair_int == 2) { $pair = "Right reads file";	
			} else { DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\nPair: $pair_int"); DOMINO::printFormat_message; DOMINO::dieNicely();
		}}
		
		if ($format_returned =~ /fastq/) {
			DOMINO::printDetails("\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n\t\tIt also contains an identifier ($name) in the name for later analysis...OK\n", $param_Detail_file);
			push (@{ $domino_files{"original"}{$name} }, $file_to_check); 
			if ($pair) {  DOMINO::printDetails("\t\tPaired end file = ".$pair." ...OK\n", $param_Detail_file);	 }
		} else { DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\n"); DOMINO::printFormat_message; DOMINO::dieNicely();
	}} else { DOMINO::printError("Please provide a valid FASTQ file: $file_to_check...It is not readable or writable or it does not exist. "); DOMINO::printFormat_message; DOMINO::dieNicely();} 
	DOMINO::printDetails("\n", $param_Detail_file);
	&debugger_print("File: $file_to_check -- $format_returned -- OK;\n");
}

sub checking_names {
	
	##########################################################################################
	##	 																					##
	##  This function checks the names of the provided databases							##
	## 	It uses the module File::Find, looks for the files we want to generate the database ##
	##	to do the search. 																	##
	## 																						##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $dir = $_[0];
	my $id;
	opendir(DIR, $dir);
	my @files = readdir(DIR);
	for (my $i = 0; $i <scalar @files; $i++) {
		if ($files[$i] =~ /.*gz/) {
			DOMINO::gunzipping($db_dirname."/".$files[$i]);
		}
	}
	close(DIR);
	opendir(DIR, $dir);
	@files = readdir(DIR);
		
	for (my $i = 0; $i <scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		if ($files[$i] =~ /UniVec\.fasta/) { next;
		} elsif ($files[$i] =~ /(chr.*)\.fa/) {
			File::Copy::move($files[$i], "./".$1.".fasta"); next;
		} elsif ($files[$i] =~ /NC_.*/) {
			$id = &get_id($files[$i]);
		}
		if ($files[$i] =~ /(.*)\.fna/) {
			File::Copy::move($files[$i], "./".$id.".fasta");
		} elsif ($files[$i] =~ /\.fasta/) {
			next;
	}}
	close(DIR);
}

sub sort_files_and_folders {

	##########################################################################################
	##	 																					##
	##  This function checks the files in the directory and place them in a clean or 		##
	## 	original folder, just to make things look nicer.									##
	## 																						##
	##	Jose Fco. Sanchez Herrero, 07/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	opendir(DIR, $dirname);
	my $random_number = int(rand(100));
	
	## Checking if the Directory already exists because a previous analysis
	if (-d $clean_folder) { 
		File::Copy::move($clean_folder, $clean_folder."_old_".$random_number);
	} 
	if (-d $intermediate_folder) {
		File::Copy::move($intermediate_folder, $intermediate_folder."_old_".$random_number);
	} 
	mkdir $clean_folder, 0755; mkdir $intermediate_folder, 0755;
	push (@{ $domino_files{"main"}{"tmp_clean_data"} }, $intermediate_folder);
	push (@{ $domino_files{"main"}{"clean_data"} }, $clean_folder);
		
	#my $domino_files_general = $intermediate_folder."/dumper-hash_DOMINO_files.txt"; DOMINO::printDump(\%domino_files, $domino_files_general);

	print "+ Deleting temporary files...\n";	
	my $files_ref = DOMINO::readDir($dirname); my @files = @$files_ref;
	my $temp_dir = $dirname."/Temporary_files_Cleaning";
	mkdir $temp_dir, 0755;
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next;
		## Folder for each taxa
		} elsif (-d $files[$i]) {
			if ($files[$i] =~ /db/) { remove_tree ($files[$i]); }		
		## Intermediate
		} elsif ($files[$i] =~ /reads_renamed.txt/){
			File::Copy::move($dirname."/".$files[$i], $intermediate_folder);
		## Logs
		} elsif ($files[$i] =~ /.*ERROR.*/) {  next;
		} elsif ($files[$i] =~ /DM_Cleaning_log.txt/) { next;
		} elsif ($files[$i] =~ /DM_command-line_log.txt/) { next; 
		} elsif ($files[$i] =~ /.*Param.*/) {
			File::Copy::move($dirname."/".$files[$i], $clean_folder); 		
		## To discard
		} elsif (-z $files[$i]) { ## Check if file is empty
			File::Copy::move($dirname."/".$files[$i], $temp_dir); 		
		} else { File::Copy::move($dirname."/".$files[$i], $temp_dir); }
	}

	foreach my $tags (keys %domino_files) {
		if ($tags eq 'original' | $tags eq 'main') {
			next;	
		} else {
			my @array_final;		
			if ($domino_files{$tags}{'FINAL_fq'}) {
				@array_final = @{$domino_files{$tags}{'FINAL_fq'}};
			} else { 
				@array_final = @{ $domino_files{$tags}{'EXTRACTED'} };
			}
			for (my $i=0; $i < scalar @array_final; $i++) { File::Copy::move($array_final[$i], $clean_folder); }
			if ($domino_files{$tags}{"CONTAMINANTS"}[0]) { File::Copy::move($domino_files{$tags}{"CONTAMINANTS"}[0], $intermediate_folder); }
			if ($domino_files{$tags}{"reads_CONTAMINANTS"}[0]) { File::Copy::move($domino_files{$tags}{"reads_CONTAMINANTS"}[0], $intermediate_folder); }
			
			## Check each taxa folder
			my $dir = $domino_files{$tags}{'DIR'}[0];
			my $files_taxa_ref = DOMINO::readDir($dir);
			my @files_taxa_ref = @$files_taxa_ref;
			for (my $f=0; $f < scalar @files_taxa_ref; $f++) {
				if ($files_taxa_ref[$f] eq ".DS_Store" || $files_taxa_ref[$f] eq "." || $files_taxa_ref[$f] eq ".." ) { next;
				} elsif ($files_taxa_ref[$f] =~ /.*_stat/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder."/Quality_Filtering_statistics_".$tags.".txt");
					push (@{ $domino_files{$tags}{"QC_stats"} }, $intermediate_folder."/Quality_Filtering_statistics_".$tags.".txt");
				} elsif ($files_taxa_ref[$f] =~ /.*unPair.*/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} elsif ($files_taxa_ref[$f] =~ /.*singleton.*/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} elsif ($files_taxa_ref[$f] =~ /.*html/) { 
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder."/Cleaning_step_report_".$tags."_QC-Analysis.html"); 		
					push (@{ $domino_files{$tags}{"QC_stats_html"} }, $intermediate_folder."/Cleaning_step_report_".$tags."_QC-Analysis.html");
										
				} elsif ($files_taxa_ref[$f] =~ /.*png/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} else {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $temp_dir);
			}}
			remove_tree($dir);
	}}
	unless ($avoidDelTMPfiles) { remove_tree($temp_dir); }
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";	
}

sub debugger_print {
	my $string = $_[0];
	my $ref = $_[1];
	## Print to a file
	if ($debugger) {
	if ($string eq "Ref") {
		print "DEBUG:\t"; print Dumper $ref; print "\n"; ## if filehandle OUT: print OUT Dumper $ref
	} else {
		print "DEBUG:\t".$string."\n";
}}}

sub download_db {
	
	##########################################################################################
	##	 																					##
	##  This function gets the default databases distributed with the package				##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	opendir(DIR, $db_dirname_default);
	my @files = readdir(DIR);	
	for (my $i = 0; $i <scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		File::Copy::copy($db_dirname_default."/".$files[$i], $db_dirname);
	}	
}

sub extracting_fastq {
	
	##########################################################################################
	##	 																					##
	##  This function generates a FASTA and QUAL file given a FASTQ file 					##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 28/04/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $file = $_[0]; my $file_name = $_[1]; my $option = $_[2];
	print "+ Extracting FASTQ file: $file\n";
	my ($fasta_file,$qual_file);
	if (!$option) {
		$fasta_file = $file_name.".fasta";
		$qual_file = $file_name.".qual";
	} else {
		$fasta_file = $file_name.".filtered.fasta";
		$qual_file = $file_name.".filtered.qual";
	}
	open (FASTA, ">$fasta_file");
	open (QUAL, ">$qual_file");
	open (F1, "$file") or DOMINO::printError("Could not open file $file") and exit();
	while(<F1>) {
		my @Read = ();
		chomp(my $id = $_);
		last if($id=~ /^\n$/);

		## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
		for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>;}
		chomp(my $QualLine = $Read[2]); chomp(my $SeqLine = $Read[0]);
		
		#print $id."\n".$SeqLine."\n".$QualLine."\n";
		
		## Print into fasta and qual files		
		my @seq_id = split ("\@", $id);
		my $qual_string = DOMINO::convert_ASCII_to_number($QualLine);
		print FASTA ">".$seq_id[1]."\n".$SeqLine."\n";
		print QUAL ">".$seq_id[1]."\n".$qual_string."\n";
	} close(F1); close (FASTA); close (QUAL);
}

sub FormatError_barcodes { DOMINO::printError("The only format accepted for the Barcodes file is:\n\n##########\nbarcode, ACACGACGACT, MID1, Dmelanogaster\nbarcode, ACTACGTCTCT, MID2, Dsimulans\nbarcode, ACTACGTATCT, MID3, Dyakuba\n##########\nPlease note this is a comma-separated value file (.csv)"); }
	
sub looking4db {

	##########################################################################################
	##	 																					##
	##  This function looks for the databases in the db folder and checks the format		##
	## 	It uses the module File::Find, looks for the files we want to generate the database ##
	##	to do the search. 																	##
	## 																						##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 	jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $dir = $db_dirname;
	for (my $i = 0; $i < scalar @user_blast_db; $i++) { File::Copy::copy($user_blast_db[$i], $db_dirname); }
	&checking_names($dir);
	find(\&print_file, $dir);

	for (my $i = 0; $i < scalar @blast_db; $i++) { DOMINO::printDetails("\t\t- $blast_db[$i]\n", $param_Detail_file); }
	sub print_file{
		my $element = $_;
		if(-f $element && $element =~ /.fa/){
			push (@blast_db, "$File::Find::name");
}}}

sub mothur_sffinfo {
	
	###########################################################################################
	##	 																					 ##
	##  This function uses the MOTHUR standalone version to extract and trim tags of the 	 ##
	##	reads of a "sff" file format provided.												 ##
	##	Jose Fco. Sanchez Herrero, 28/04/2014 jfsanchezherrero@ub.edu						 ##
	## 																						 ##
	###########################################################################################
	
	# This subroutine takes as input a sff file and extracts it, classifies according to the tags provided, by Roche and trims the seqs
	my $file = $_[0];
	my $file_name = $_[1];
	my $directory = $dirname;

	## The '#' is necessary to be in the first place when calling mothur with a given set of commands
	my $mothur_call = $mothur_path." '#set.dir(output=$directory); sffinfo(sff='".$file."'); trim.seqs(fasta="."$file_name".".fasta, qfile=".$file_name.".qual, oligos=ROCHE.oligos, bdiffs=$bdiffs, processors=$noOfProcesses)'";
	&debugger_print("MOTHUR command (sff info): $mothur_call\n");
	my $mothur_result = system($mothur_call);
	if ($mothur_result != 0) { DOMINO::printError("MOTHUR failed when trying to proccess the file..."); exit(); } else { print "Done...OK\n\n"; }	
}

sub mothur_trim_fastq {
	
	# This subroutine takes as input a FASTQ file AND classifies according to the tags provided, by Roche and trims the seqs
	my $file_name = $_[0];
	my $directory = $dirname;

	## The '#' is necessary to be in the first place when calling mothur with a given set of commands
	my $mothur_call = $mothur_path." '#set.dir(output=".$directory."); trim.seqs(fasta=".$file_name.".fasta, qfile=".$file_name.".qual, oligos=ROCHE.oligos, bdiffs=".$bdiffs.", processors=".$noOfProcesses.")'";
	&debugger_print("MOTHUR command (trim fastq): $mothur_call\n");
	print "+ Triming the reads according to MID tag\n";
	my $mothur_result = system($mothur_call);
	if ($mothur_result != 0) { 	DOMINO::printError("MOTHUR failed when trying to proccess the file..."); exit(); } else { print "Done...OK\n\n"; }
}

sub prinseq_dust {
	
	###########################################################################################
	##	 																					 ##
	##  This function uses the script prinseq-lite.pl distributed in PRINSEQ 				 ##
	##  It would process reads containing low complexity under a threshold given             ##
	##	and the reads containing ambiguous bases and N's								     ##
	## 																						 ##
	###########################################################################################
	
	my $file1 = $_[0]; my $file2 = $_[1]; my $thr = $_[2]; my $dir = $_[3];
	
	my $prinseq_dust_command = "perl ".$scripts_path."PRINSEQ-lite_0.20.4/prinseq-lite.pl";
	my @name = split("/", $file1);
	my @id = split(".fastq", $name[-1]); 
	my $out_good = $dir."/".$id[0].".clean_dust";
	my $out_bad = $dir."/".$id[0].".dust";
	
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		$prinseq_dust_command .= " -fastq ".$file1." -fastq2 ".$file2;
	} else { $prinseq_dust_command .= " -fastq ".$file1; }

	if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
		$prinseq_dust_command .= " -out_format 2"; ## Generate fasta +  qual for NGSQC step
	}
	$prinseq_dust_command .= " -out_good ".$out_good." -out_bad ".$out_bad." -lc_method dust -lc_threshold ".$thr; 
	$prinseq_dust_command .= " -ns_max_n 1 -noniupac 2> $error_log";
	&debugger_print("PRINSEQ command:\n$prinseq_dust_command\n");
	my $dust_result = system($prinseq_dust_command);	
	if ($dust_result != 0) { DOMINO::printError("DUST cleaning step failed when trying to proccess the file... DOMINO would not stop in this step...");  }
	return ($out_good, $out_bad);	
}

sub rename_seqs {

	my $file = $_[0];
	my $id_Taxa = $_[1];
	my $typeFile = $_[2];
	
	my @file_name = split("/", $file);
	my $dir = $dirname."/".$id_Taxa;
	unless (-d $dir) { mkdir $dir, 0755; }
	my $file_path = $dir."/renamed_".$file_name[-1];
	my ($codeReturn, $messageDebugger) = DOMINO::check_ID_length($file, "fastq");
	&debugger_print($messageDebugger);
	my ($pairReturn, $type);
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		($pairReturn, $type) = DOMINO::check_paired_file($_[0], $typeFile);
		&debugger_print("Pair: ".$pairReturn."\t"."Type: ".$type);
		if ($pairReturn eq "L" ) { 
			$codeReturn = 0; ## Force to rewrite seq					
			&debugger_print("Lets rename $file : pair was $pairReturn");
		} elsif ($pairReturn eq "R") { 
			$codeReturn = 0; ## Force to rewrite seq
			&debugger_print("Lets rename $file : pair was $pairReturn");
		} # else {
			## There is no need to change file according to pair
	}	
	
	if ($codeReturn == 1) {
		## Length < 40 char
		my $name = $dir."/".$file_name[-1];
		return ($name);
		&debugger_print("No need to change file if of $file neither because of length nor type of pair");		
	
	} else {
		&debugger_print("Rename file $file");		
		open (F1, "$file") or DOMINO::printError("Could not open file $file") and exit();
		open (FASTQ_out, ">$file_path");
		my $int = 0;
		while(<F1>) {
			my @Read = ();
			chomp(my $id = $_);
			last if($id=~ /^\n$/);
			## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
			for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>; }
			chomp(my $qual = $Read[2]);
			chomp(my $seq = $Read[0]);
			$int++;			
			my $id_new = "Seq_".$int."_".$id_Taxa;
			my $new_pair;
			if ($pairReturn) { 
				if ($pairReturn eq "L" ) { $new_pair = "1";					
				} elsif ($pairReturn eq "R") { $new_pair = "2";					
				} else { $new_pair = $pairReturn; }
				# /1|/2 or 1:N|2:N
				$id_new .= "/".$new_pair;
			}
			if (length($id_new) > 40) {
				$id_new = "Seq_".$int;
				if ($pairReturn) { $id_new .= "/".$new_pair; }
			}
			print FASTQ_out "\@".$id_new."\n".$seq."\n+\n".$qual."\n";
		}
		close(F1); close (FASTQ_out);
		return $file_path;
	}	
}

sub splitting_fastq {

	my $file_name = $_[0];
	my $int = $_[1];
	my $general_dir = $_[2];
	my $fasta_file = $general_dir."/".$file_name.".trim.fasta";
	my $qual_file = $general_dir."/".$file_name.".trim.qual";
	my $group_file = $general_dir."/".$file_name.".groups";
	my %files_generated;

	print "- Checking files generated and splitting into multiple threads...\n";
	
	$/ = "\n"; ## Telling Perl where a new line starts
	open (GROUPS, $group_file) or DOMINO::printError("Cannot open file groups: $group_file") and exit();
	# Generate a hash with the identifier for each sequence
	while (<GROUPS>) {
		chomp; my $line = $_;
		$line =~ s/\s+/\t/g; $line =~ s/\t+/\t/g;
		my @array = split ("\t", $line);
		my $filehandle;
		if ($domino_files{$array[1]}{'IDS'}) {
			$filehandle = $domino_files{$array[1]}{'IDS'}[0];
		} else {
			$filehandle = $general_dir."/".$file_name."-".$array[1]."_ids.txt";
			push (@{ $domino_files{$array[1]}{'IDS'} }, $filehandle);
			push (@{ $domino_files{$array[1]}{'DIR'} }, $general_dir."/".$array[1]);
		}
		open (OUT, ">>$filehandle"); print OUT $array[0]."\n"; close(OUT);
	} close (GROUPS);
	&debugger_print("Taxa ids reads:\n"); &debugger_print("Ref", \%domino_files);
	
	my $pm_mothur =  new Parallel::ForkManager($noOfProcesses);
	$pm_mothur->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
		print "\n\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
	} );
	$pm_mothur->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** Mothur retrieving sequences started for $ident with PID $pid\n\n"; } );
	
	foreach my $files_ids (sort keys %domino_files) {
		next unless ( exists ($domino_files{$files_ids}{'IDS'}));
		my $pid = $pm_mothur->start($files_ids) and next;
		print "\nSending child process for retriving reads for $files_ids\n\n";
		my $dir = $general_dir."/".$files_ids; unless (-d $dir) { mkdir $dir, 0755; }
		DOMINO::mothur_retrieve_seqs($fasta_file, $qual_file, $dir, $domino_files{$files_ids}{'IDS'}[0], $mothur_path);
		push (@{ $files_generated{$files_ids}{'FASTA_extracted'} }, $dir."/".$file_name.".trim.pick.fasta");
		push (@{ $files_generated{$files_ids}{'QUAL_extracted'} }, $dir."/".$file_name.".trim.pick.qual");			
		my $name;
		if ($file_type eq "Illumina_pair_end") {
			my $pair; my $type;
			($pair, $type) = DOMINO::check_paired_file($dir."/".$file_name.".trim.pick.fasta", "fasta");
			$name = $dir."/reads_id-".$files_ids."_R".$pair.".fastq";
		} else { $name = $dir."/reads_id-".$files_ids.".fastq"; }

		my $fastq_file = DOMINO::qualfa2fq_modified_bwa($dir."/".$file_name.".trim.pick.fasta", $dir."/".$file_name.".trim.pick.qual", $files_ids, $name, $file_type);
		push (@{ $files_generated{$files_ids}{'EXTRACTED'} }, $fastq_file);			
		my $extracted = $dir."/dumper_extracted.txt";
		&debugger_print("####### Extracted dumper ".$extracted); 
		DOMINO::printDump(\%files_generated, $extracted);
		$pm_mothur->finish($int_taxa); # pass an exit code to finish
	}
	$pm_mothur->wait_all_children;
	print "\n** All Mothur retrieve processes have finished...\n\n";	
}

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}


__END__
sub printing_statistics {

	## TODO
	my $hash_Ref = $_[0];
	my %reads = %{$hash_Ref};
	my $original_reads;
	## Print into console output, text file and excel file
	if ($file_type eq "Illumina_multiple_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "Illumina_pair_end_multiple_fastq") { 
		print "+ Several original files provided:\n";
		foreach my $files (keys %reads) {
			if ($files =~ /original.*/) { print "\t".$reads{$files}[0]." containing ".$reads{$files}[1]." reads\n"; }
		}	
		print "\n";
		if ($file_type eq "Illumina_pair_end_multiple_fastq") {
			print "+ Files were merged into unique paired-end read merged files:\n\t$reads{Merged_multiple_fastq}[0] containing $reads{Merged_multiple_fastq}[1] reads\n\t$reads{Merged_multiple_fastq_1}[0] containing $reads{Merged_multiple_fastq_1}[1] reads\n";
		} else {
			print "+ Files were merged into a unique merged file:\n\t$reads{Merged_multiple_fastq}[0] containing $reads{Merged_multiple_fastq}[1] reads\n";
		}
		$original_reads = $reads{Merged_multiple_fastq}[1];
	} elsif ($file_type eq "Illumina_pair_end") {
		print "+ Original files:\n\t$reads{original}[0] contains $reads{original}[1]\n\t$reads{original_1}[0] contains $reads{original_1}[1]\n";
		$original_reads = $reads{original}[1];
	} else {
		print "+ Original file:\n\t$reads{original}[0] contains $reads{original}[1]\n";
		$original_reads = $reads{original}[1];
	}
	my $cleaned_Reads_step1_1 = $original_reads - $reads{"NGS_QC_filtered"}[1];
	print "\n+ Cleaning Step 1:\n";
	print "\tQuality control 1.1: short reads and poor quality reads:\t$cleaned_Reads_step1_1 reads were discarded\n";
	my $dust_reads;
	if ($file_type eq "Illumina_pair_end_multiple_fastq" || $file_type eq "Illumina_pair_end") {
		$dust_reads = $reads{"NGS_QC_filtered"}[1] - $reads{"DUST_filtered"}[1];
	} else {
		$dust_reads = $reads{"DUST"}[1];
	}
	
	print "\tQuality control 1.2: low complexity reads: DUST algorithm\t".$dust_reads." reads were discarded\n\n";

	if ($skipping_BLAST) {
		print "+ No BLAST step for searching against databases for putative contaminants\n";
	} else {
		print "+ Cleaning Step 2: BLAST and Vector screen\n";
		my $deleted_reads = $reads{"DUST_filtered"}[1] - $reads{"BLAST_filtered"}[1];
		print "\t $deleted_reads reads were discarded\n";
	}

	print "\n";
	unless ($file_type eq "454_fastq" || $file_type eq "Illumina" || $file_type eq "454_sff" || $file_type eq "Illumina_multiple_fastq" || $file_type eq "454_multiple_fastq") { 
	## Paired-end reads files
		print "+ Due to the cleaning process, some pairs remain unpaired:\n\tUnpaired reads are in file: ".$reads{unPaired_Reads}[0]." containing ".$reads{unPaired_Reads}[1]." reads\n";		
	}
	
	print "+ Summary:\n";
	
	my $total_discarded_reads = $original_reads - $reads{"BLAST_filtered"}[1];
	my $percentage = ($total_discarded_reads/$original_reads)*100;
	my $h = sprintf ("%.3f", $percentage);
	print "\tA total amount of $total_discarded_reads reads were discarded during the cleaning process acounting for: $h %...\n";
	print "\tCleaned reads are in file: ".$reads{"BLAST_filtered"}[0]." containing ".$reads{"BLAST_filtered"}[1]." reads: ".(100 - $h)."\n";
}

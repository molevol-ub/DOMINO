#!/usr/bin/perl
#######################################################################################
###	DOMINO: Development of molecular markers in non-model organisms using NGS data 	###
###																					###
###	Authors:																		###
###	Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro 		###
###	Sánchez-Gracia, and Julio Rozas.					     						###
###																					###
#######################################################################################
##	Usage:
##      perl DM_Clean_v1.0.0.pl
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
	require List::Uniq;
	require File::Copy;
	require File::Find;
	use File::Find qw(find);			
	require File::Path;
	use File::Path qw(remove_tree);
	require Cwd;
	use Cwd qw(abs_path);  
	require Parallel::ForkManager;
}
## TODO: Write a debug message if importing any module fails

my $pipeline_path = abs_path($0); 
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/";  }

## TODO
## print each QC thread analysis into a file not to show so many confusing messages

##################################
##	Initializing some variables	##
##################################
my (
## Options provided
$bdiffs, $read_length_GUI, $minimum_qual_GUI, $minimum_length_GUI, $outFolder, 
$manual, $version, $threshold, $noOfProcesses, $user_barcodes_file, $helpAsked, 
$avoidDelTMPfiles, @user_blast_db, $skipping_BLAST, $user_option_file_type, $only_user_db, @blast_db,
$further_information, $onlyTagging_files, $debugger, @user_blast_db_tmp, @user_files,

## Others
%taxa_dir, %species_names, $step_time, %domino_files, @file_abs_path, @file_names
);	

my %input_options = (
	1 =>'454_sff', 2 =>'454_fastq', 3=>'454_multiple_fastq', 4 => 'Illumina', 
	5 => 'Illumina_multiple_fastq', 6 => 'Illumina_pair_end', 
	7 => 'Illumina_pair_end_multiple_fastq');

######################
## Get user options	##
######################
GetOptions(
	"type_file=i" => \$user_option_file_type,
	"input_file=s" => \@user_files,  
	"bdiffs=i" => \$bdiffs,
	"l|cut_Off_ReadLen=i" => \$read_length_GUI,
	"s|cut_Off_Quality=i" => \$minimum_qual_GUI,
	"m|minLen=i" => \$minimum_length_GUI,
	"thr|threshold_DUST=i" => \$threshold,
	"o|outputFolder=s" => \$outFolder, 
	"p|number_cpu=i" => \$noOfProcesses,
	"b|barcode_file=s" => \$user_barcodes_file,
			
	"db|blast_database=s" => \@user_blast_db_tmp,
			
	"TempFiles" => \$avoidDelTMPfiles,
	"no_db_search" => \$skipping_BLAST,
	"only_user_db" => \$only_user_db,
	"only_tag_files|OTG" => \$onlyTagging_files,
			
	"h|help" => \$helpAsked,
	"man" => \$manual,
	"v|version" => \$version,
	"MoreInfo" => \$further_information,
	
	"Debug" => \$debugger,
);

## Help/manual/version
pod2usage( -exitstatus => 0, -verbose => 1 ) if ($helpAsked);
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

DM_Clean_v1.0.0.pl 

=back
	
=head1 VERSION

=over 2

=item B<>

DOMINO v1.0.0

=back

=head1 SYNOPSIS

=over 2

=item B<>

perl DM_Clean_v1.0.0.pl 

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

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 1 -o test_SFF_file 
 -input file.sff -b MID_barcodes.txt

=item B<>

=item B<454 Fastq: Unique file>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 2 -o test_454_fastq 
 -input file.fastq -b MID_barcodes.txt

=item B<>

=item B<454: Multiple fastq files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq 

=item B<>

=item B<Illumina: Unique file>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 4 -o test_Illumina 
 -input Illumina_Multiple-taxa.fastq -b MID_barcodes.txt 

=item B<>

=item B<Illumina: Multiple fastq files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 5 -o test_Illumina_multi_files 
 -input sp1_Illumina.fastq -input sp2_Illumina.fastq -input sp3_Illumina.fastq 
 
=item B<>

=item B<Illumina Paired-end Files: Unique file containing all taxa>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 6 -o test_Illumina_pair_end 
 -b MID_barcodes.txt -input reads_R1.fastq -input reads_R2.fastq 

=item B<>

=item B<Illumina Paired-end Files: Multiple Paired-end FASTQ files>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 7 -o test_mult_pair_end 
 -input sp1_R1.fastq -input sp1_R2.fastq -input sp2_R1.fastq -input sp2_R2.fastq

=item B<>

=item B<Using User provided databases>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq -db file1.fasta -db file2.fasta 
 -only_user_db

=item B<>

=item B<No Contaminant Search>

 perl [DOMINO_scripts_path]/DM_Clean_v1.0.0.pl -type_file 3 -o test_454_multiple 
 -input sp1.fastq -input sp2.fastq -input sp3.fastq -no_db_search

=item B<>

=back	

=head1 AUTHOR

=over 2

Cristina Frias-Lopez, Jose F. Sanchez-Herrero, Miquel A. Arnedo, Alejandro Sanchez-Gracia and Julio Rozas.
	
Evolutionary Genomics and Bioinformatics Group, Departament de Genètica, Microbiologia i Estadística and Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Av. Diagonal 643, Barcelona 08028, Spain
	
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

1 - 04 - 2016

=back	

=head1 CITATION

=over 2

To add Reference when available

=back

=cut
##########################################################################################
##					Checking user options												##
##########################################################################################
if (!$outFolder || !$user_option_file_type) { 
	Pod::Usage::pod2usage("Try 'perl $0 -man' for more information\n", -exitstatus => 0, -verbose => 0 ); 
}
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
if (!$dirname) { exit(); }
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
print "\n"; &print_Header("","#"); &print_Header(" DOMINO Cleaning Stage ","#"); 
&print_Header("","#"); print "\n"; &print_Header("","+");  
&print_Header(" Analysis Started ","+");  &print_Header("","+"); 
&print_DOMINO_details("Starting the process: [ ".(localtime)." ]\n\n");

## Getting scripts path variable
my $BLAST = $scripts_path."NCBI_BLAST_v2.2.28/";
my $db_dirname_default = $scripts_path."db_default"; # Default databases provided
my $db_dirname = $dirname."/db";
if ($user_barcodes_file) { 
	my $tmp_user_barcodes_file_abs_path = abs_path($user_barcodes_file);
	push (@{$domino_files{"original"}{"barcodes_file"}}, $tmp_user_barcodes_file_abs_path); 
}

&print_Header(" Input File and Parameter Preprocessing ","#");
&print_DOMINO_details("\n+ Output Directory: ".$dirname." ...OK\n\n");
if (scalar (@user_files) == 0 || !$user_option_file_type) {
	&printError("No input files provided"); &dieNicely();
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

## Check taxa names provided
## User provides the taxa names according to the type of file:
### 1) File containing multiple taxa accordingly tag by manufacturer using MID barcodes
##	- User provides a unique fastq/sff file and barcode file
## 	"
##	barcode,ATTGCTACGAT,MID2,Dmelanogaster
##	barcode,ACAGCTACGAT,MID6,Dsimulans
##	barcode,ATTGCAGCGAT,MID1,Dyakuba
##  "	

## Option 1: User barcodes provided
my @MID_species_array; my %species_names_tmp;
if ($user_option_file_type == 1 || $user_option_file_type == 2 || $user_option_file_type == 4 || $user_option_file_type == 6 ) {
	if ($user_barcodes_file) {
		my %hash_Return = &check_barcode_user_file($domino_files{"original"}{"barcodes_file"}[0]);
		%species_names_tmp = %hash_Return;
	} else {
		&printError("\n\nERROR: No barcodes file provided...\n"); &dieNicely();
}}

##########################################################################################
##				Checking files															##
##########################################################################################
my $file_type = $input_options{$user_option_file_type};
if (!$input_options{$user_option_file_type}) {
	&printError("\n\nERROR: Wrong type of file provided\nPlease provide a valid type of file:\n");
	&printError("According to the type of NGS files provided several options are available but only one option would be provided: -type_input [int]
 1: A single file in Standard Flowgram Format (SFF), 454 file, containing all the reads of the different taxa accordingly tagged
 2: FASTQ files coming from 454 Roche. A single file containing all the reads of the different taxa accordingly tagged
 3: Multiple FASTQ files coming from 454 Roche. Each file contains each taxa reads. 
 4: FASTQ file from Illumina single end, containing all the reads of the different taxa accordingly tagged
 5: Multiple FASTQ file from Illumina single end. Each file contains each taxa reads. 
 6: A single pair of FASTQ files from Illumina paired-end sequencing: Each pair would contain all the reads of the different taxa accordingly tagged. Please tagged left read file with xxx_R1.fastq and right read file as xxx_R2.fastq
 7: Multiple FASTQ files from Illumina paired-end: Each pair of files containing each taxa left and right reads respectively\n\n");
	Pod::Usage::pod2usage("Try 'perl $0 -man' for more information\n");	
	exit();
} else {
	&print_DOMINO_details("+ Type of file(s): Option: $user_option_file_type -- $file_type ...OK\n");
}
&print_DOMINO_details("+ Checking file(s):\n");

#### Type file: 1 
if ($file_type eq "454_sff") {
	if (scalar (@file_abs_path) != 1) {
		&printError("Please provide a unique SFF file"); &dieNicely();
	} elsif ($file_abs_path[0] !~ /(.*)\.sff/) {
		&printError("Wrong 454 file provided.\nPlease provide a binary SFF extension file or provide a different type of file"); &dieNicely();
	} elsif (-e -r -s $file_abs_path[0]) {
		my $id = $file_names[0];
		&print_DOMINO_details("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = SFF ...OK\n\t\tIt also contains an identifier ($id) in the name for later analysis...OK\n\n");	
		push (@{$domino_files{"original"}{'SFF'}}, $file_abs_path[0]); 

#### Type file: 2
}} elsif ($file_type eq "454_fastq") {
	if (scalar (@file_abs_path) != 1) {
		&printError("Please provide a unique 454 Roche FASTQ file containing all the taxa read accordingly tagged with MID sequences"); &dieNicely();
	} elsif (scalar (@file_abs_path) == 1) {
		my $format = &check_file_format($file_abs_path[0]);
		if ($file_abs_path[0] =~ /.*\.f.*q$/) {
			if ($format =~ /fastq/) {
				&print_DOMINO_details("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n");		
				push (@{$domino_files{"original"}{'454_fastq'}}, $file_abs_path[0]); 

			} else { 
				&printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[0]... Format: $format"); &dieNicely();
		
#### Type file: 3
}}}} elsif ($file_type eq "454_multiple_fastq") {
	if (scalar (@file_abs_path) == 1) {
		&printError("Please provide multiple 454 Roche FASTQ file containing each of the taxa reads"); &printFormat_message(); &dieNicely();
	} else {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);

#### Type file: 4
}}} elsif ($file_type eq "Illumina") {
	if (scalar (@file_abs_path) != 1) {
		&printError("Please provide a unique Illumina FASTQ file containing all the taxa read accordingly tagged with MID sequences"); &dieNicely();
	} elsif (scalar (@file_abs_path) == 1) {
		my $format = &check_file_format($file_abs_path[0]);
		if ($file_abs_path[0] =~ /.*\.f.*q$/) {
			if ($format =~ /fastq/) {
				&print_DOMINO_details("\t$file_abs_path[0]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n");		
				push (@{$domino_files{"original"}{'illu_FASTQ'}}, $file_abs_path[0]); 
			} else { 
				&printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[0]... Format: $format"); &dieNicely();

#### Type file: 5
}}}} elsif ($file_type eq "Illumina_multiple_fastq") {  
	if (scalar (@file_abs_path) == 1) {
		&printError("Please provide multiple Illumina FASTQ file containing containing each of the taxa reads"); &printFormat_message(); &dieNicely();
	} elsif (scalar (@file_abs_path) != 1) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);

#### Type file: 6
}}} elsif ($file_type eq "Illumina_pair_end") {
	if (scalar (@file_abs_path) != 2) {
		&printError("Please provide Illumina Paired-End FASTQ files. Provide the left (_R1) and then the right (_R2) file containing all the taxa read accordingly tagged with MID sequences\nPlease tagged left read file with xxx_R1.fastq and right read file as xxx_R2.fastq\n"); &dieNicely();
	} elsif (scalar (@file_abs_path) == 2) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			my $format = &check_file_format($file_abs_path[$i]);
			if ($file_abs_path[$i] =~ /.*\_R(\d+)\.f.*q$/) { ## .*_R[*].fastq
				my $pair_int = $1; my $pair;
				if ($pair_int == 1) { $pair = "Left reads file";	
				} elsif ($pair_int == 2) { $pair = "Right reads file";	
				} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[$i]...\nFormat: $format\nPair: $pair_int\n"); &dieNicely(); }
				if ($format =~ /fastq/) {
					&print_DOMINO_details("\t$file_abs_path[$i]\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format." ...OK\n\t\tPaired end file = ".$pair." ...OK\n");
					push (@{$domino_files{"original"}{'illu_PE'}}, $file_abs_path[$i]);
				} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_abs_path[$i]...\nFormat: $format\nPair: $pair_int\n"); &dieNicely();

#### Type file: 7
}}}}} elsif ($file_type eq "Illumina_pair_end_multiple_fastq") {  
	if (scalar (@file_abs_path) == 2 ) {
		&printError("Please provide Illumina Paired-End FASTQ files for each taxa. Provide the left (_R1) and then the right (_R2) file containing each taxa reads"); &printFormat_message(); &dieNicely();
	} elsif (scalar (@file_abs_path) != 2) {
		for (my $i = 0; $i < scalar @file_abs_path; $i++) {
			&check_file($file_abs_path[$i]);
}}}

## Check if paired end files provided, all pairs have been correctly provided
if ($file_type eq "Illumina_pair_end_multiple_fastq") {
	foreach my $keys (keys %species_names) {
		my $sum_pairs; my $files;
		##@{ $domino_files{"original"}{$name} }
		my @array = @{ $domino_files{"original"}{$keys} };
		if (scalar @array != 2) {
			&printError("Wrong pair of FASTQ files provided. Please check pair of files corresponding to $keys:\n$files\n"); &dieNicely();
		} else { &debugger_print("OK: $keys\n\n"); }
}}

&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

##########################################################################################
##					Printing user options												##
##########################################################################################
if ($onlyTagging_files) {
	&print_DOMINO_details("+ Cleaning process would be skipped...OK\n");
	&print_DOMINO_details("+ Only extracting and tagging of the reads would be done...OK\n");	
} else {
	## Print information of the input parameters
	&print_DOMINO_details("+ Number of Processor: ".$noOfProcesses." ...OK\n");
	&print_DOMINO_details("+ Number of mismatches allowed in the barcode sequence: $bdiffs ...OK\n");	
	&print_DOMINO_details("+ Minimum read length: $minimum_length_GUI pb ...OK\n");
	&print_DOMINO_details("+ Minimum QUAL: $minimum_qual_GUI ...OK\n");
	&print_DOMINO_details("+ Minimum length of a read satisfying QUAL cutoff: $read_length_GUI % ...OK\n");
	&print_DOMINO_details("+ Threshold for complexity/entropy of a read: $threshold ...OK\n");
	unless ($avoidDelTMPfiles) {
		&print_DOMINO_details("+ Deleting of temporary files would be done ...OK\n");
	} else { 
		&print_DOMINO_details("+ Deleting temporary files would be avoid ...OK\n");
	}

	## Check databases to use in BLAST contaminations search
	unless ($skipping_BLAST) {
		&print_DOMINO_details("+ Checking Databases to use during BLAST search of contaminants:\n\tA folder named as 'db' would be generated to place the databases we are going to use\n");
		mkdir $db_dirname, 0755; chdir $db_dirname;
		if ($only_user_db) {
			&print_DOMINO_details("\t\tNo default databases would be used and only user provided databases would be used instead\n");
		} else {
			&print_DOMINO_details("\tDefault databases provided in the package would be used\n\tCopying...\n");
			&download_db(); 
		}

		## Check if User provide any database
		if (@user_blast_db_tmp) {
			&print_DOMINO_details("\t+ Additional user provided databases would be used:\n");
			for (my $i = 0; $i < scalar @user_blast_db_tmp; $i++) {
				print "\t$user_blast_db_tmp[$i]\n";
				my $tmp_abs_name = abs_path($user_blast_db_tmp[$i]);
				my $format_returned = &check_file_format($tmp_abs_name);
				if ($format_returned !~ /fasta/) { 
					&printError("Please, enter a valid FASTA file, $user_blast_db[$i] is not valid"); &dieNicely();
				} else {
					push (@user_blast_db, $tmp_abs_name);
		}}}
		&print_DOMINO_details("\n\t+ Database(s):\n"); &looking4db();
	} else {
		&print_DOMINO_details("+ No BLAST search of contaminants ...OK\n");
}}
## Print info about where to print info and error
&print_DOMINO_details("\n+ Parameters details would be print into file: $param_Detail_file...\n");
&print_DOMINO_details("+ Errors occurred during the process would be print into file: $error_log...\n\n");
chdir $dirname;

##########################################################################################
##					Start the analysis													##
##########################################################################################
if ($user_option_file_type == 1 || $user_option_file_type == 2 || $user_option_file_type == 4 || $user_option_file_type == 6) {
	
	##################################################################################################
	## If the type of file is SFF, 454_FASTQ, Illumina and Illumina paired-end and it is a unique   ##
	## file containing information of the different taxa, tagged with MID tags we would extract  	##
	## and classify 																  				##
	##################################################################################################
	
	print "\n"; &print_Header(" Extracting taxa sequences from the file provided ","%"); 
	my $file_tmp_barcodes = $domino_files{"original"}{"barcodes_file"}[0];
	&print_DOMINO_details("+ Barcodes file provided: $file_tmp_barcodes ...OK\n");
	
	for (my $i=0; $i < scalar @file_abs_path; $i++) {
		if ($file_type eq "454_sff") { 
			# Extracting the sff file given
			print "+ Extracting information from the SFF file provided as input\n\t-> This might take a while for large files\n- Entering to MOTHUR executable";
			&mothur_sffinfo($file_abs_path[$i], $file_names[$i]);
			print "\n- Exiting MOTHUR executable\n- Job done: sff file has been extracted and sequences has been trimmed according to their tag\n";
			&splitting_fastq($file_names[$i], $i);
		} else {
			&extracting_fastq($file_abs_path[$i], $file_names[$i]); 
			&mothur_trim_fastq($file_names[$i]); &splitting_fastq($file_names[$i], $i);
	}}
			
	## Exit the script and wait for user to input the taxa names according to tags identified
	&time_stamp();	
	&debugger_print("domino_files"); &debugger_print("Ref", \%domino_files);
	my $array_ref = &read_dir($dirname);
	for (my $i=0; $i < scalar @$array_ref; $i++) {
		if ($$array_ref[$i] eq "." || $$array_ref[$i] eq ".." || $$array_ref[$i] eq ".DS_Store") { next;}
		if ($$array_ref[$i] eq "db") {next;}
		if (-d $$array_ref[$i]) {
			$taxa_dir{$$array_ref[$i]} = $dirname."/".$$array_ref[$i];				
			$species_names{$$array_ref[$i]} = $$array_ref[$i]; # get species tag
			my $array_ref2 = &read_dir($$array_ref[$i]);
			push (@{ $domino_files{$$array_ref[$i]}{"DIR"}}, $taxa_dir{$$array_ref[$i]});
			my $dump_file = $taxa_dir{$$array_ref[$i]}."/dumper_extracted.txt";
			open (DUMP_IN, "$dump_file");
			while (<DUMP_IN>) {
				my $line = $_; chomp $line;
				my @array = split("\t", $line);
				&debugger_print($line);
				push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
			} close (DUMP_IN);
	}}
	&debugger_print("Taxa names:\n"); &debugger_print("Ref", \%species_names);
} 
&debugger_print("domino_files"); &debugger_print("Ref", \%domino_files);

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
				push (@{ $domino_files_threads_OTG{$keys}{"FINAL"} }, $copy_name);			
			}
			my $domino_files_threads_OTG_text = $dir."/dumper-hash_OTG.txt";
			&print_dump(\%domino_files_threads_OTG, $domino_files_threads_OTG_text);			
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
	&create_folder_cleaning_step(); &delete_temp_file(); 
	print "\n+ Exiting the script...\n"; &finish_time_stamp(); exit();	 
		
} else {
	
	## Maybe copying big files could take a while ## TODO: implement threads here
	if ($file_type eq "454_multiple_fastq" || $file_type eq "Illumina_multiple_fastq" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		## Copy files and get ready for QC
		foreach my $keys (keys %{ $domino_files{'original'} }) {
			$species_names{$keys} = $keys;
			my $dir = $dirname."/".$keys;
			unless (-d $dir) { 
				mkdir $dir, 0755; push (@{ $domino_files{$keys}{"DIR"} }, $dir);
				$taxa_dir{$keys} = $dir;
			}			
			for (my $i=0; $i < scalar @{ $domino_files{'original'}{$keys} }; $i++) {
				my @path_name = split("/", $domino_files{'original'}{$keys}[$i]);
				my $copy_name;
				if ($path_name[-1] =~ /.*\_id\-.*/) {
					$copy_name = $taxa_dir{$keys}."/".$path_name[-1];
				} else {
					$copy_name = $taxa_dir{$keys}."/reads_id-".$path_name[-1];
				}
				File::Copy::copy($domino_files{'original'}{$keys}[$i], $copy_name);
				push (@{ $domino_files{$keys}{"EXTRACTED"} }, $copy_name);
	}}}
	print "\n"; &time_stamp(); print "\n";
}
&debugger_print("domino_files"); &debugger_print("Ref", \%domino_files);
&print_Header(" Quality Control Analysis ","#");
# Concatenating databases
my $merge_fasta_database = $db_dirname."/merge_seqs.fasta";
open (DB, ">$merge_fasta_database");
unless ($skipping_BLAST) {
	for (my $i = 0; $i < @blast_db; $i++) {	
		my $file = $blast_db[$i];
		push (@{ $domino_files{'original'}{'db'} }, $file);
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

## Send threads for each taxa
my $int_taxa = 0;
## Sent child process
my $pm =  new Parallel::ForkManager($noOfProcesses);
$pm->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
	print "\n\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
} );
$pm->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** QC Analysis started for $ident with PID $pid\n\n"; } );

## QC analysis for each taxa
foreach my $keys (sort keys %species_names) {	
	$int_taxa++;
	my $pid = $pm->start($keys) and next; print "\nSending child process for renaming fastq for $keys\n\n";
	
	my %domino_files_threads_QC;
	chdir $taxa_dir{$keys}; print "Changing dir to $taxa_dir{$keys}\n";

	##############################################################################
	## 		Quality control step 1.1 -- PRINSEQ Dust algorithm					##
	##############################################################################
	print "\n+ Quality control 1.1: low complexity reads: DUST algorithm: $keys\n";
	my $dust_clean_reads; my $dust_reads; 
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {		
		($dust_clean_reads, $dust_reads) = &prinseq_dust($domino_files{$keys}{'EXTRACTED'}[0], $domino_files{$keys}{'EXTRACTED'}[1], $threshold, $taxa_dir{$keys});
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads."_1.fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads."_2.fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads."_1.fastq");
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads."_2.fastq");
	} else {
		($dust_clean_reads, $dust_reads) = &prinseq_dust($domino_files{$keys}{'EXTRACTED'}[0], $domino_files{$keys}{'EXTRACTED'}[0], $threshold, $taxa_dir{$keys});
		push (@{$domino_files_threads_QC{$keys}{"DUST_filtered"}}, $dust_clean_reads.".fastq");		
		push (@{$domino_files_threads_QC{$keys}{"DUST"}}, $dust_reads.".fastq");
	}
	print "\n"; &time_stamp(); print "\n";

	&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

	##############################################################################
	## 		Quality control step 1.2 -- NGS_QC_toolkit							##
	##############################################################################
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
	$NGS_QC_call .= " -t 2 -o ".$taxa_dir{$keys}." -p 1"." -l ".$read_length_GUI." -s ".$minimum_qual_GUI." 2> $error_log";
		# Number of processor to use (ONLY 1 because this would be already parallelized)
		# Setting the cutoff percentage of read length that should be of given quality [Default 70]
		# Setting the minimum cutoff for PHRED quality, default 20 
	
	## Sending command
	print "NGS QC Toolkit command: ".$NGS_QC_call."\n";
	my $NGSQC_result = system($NGS_QC_call);
	if ($NGSQC_result != 0) { &printError("Cleaning step failed when calling NGSQC for file $dust_clean_reads..."); &dieNicely(); }	
	print "\n"; &time_stamp();

	&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

	##############################################################################
	## 		Quality control step 2 -- BLAST contaminant and vector search		##
	##############################################################################
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
				push (@{ $domino_files_threads_QC{$keys}{"FINAL"} }, $final_name);
				File::Copy::move( $domino_files_threads_QC{$keys}{"NGS_QC"}[$i], $domino_files_threads_QC{$keys}{"FINAL"}[$i]);
			}
		} else {			
			my $final_name;
			my @path_name = split("/", $domino_files_threads_QC{$keys}{"NGS_QC"}[0]);
			if ($path_name[-1] =~ /.*\_id\-(.*)\.clean\_dust.*/) {
				for (my $j=0; $j < scalar $#path_name; $j++) {
					$final_name .= $path_name[$j]."/";
				} $final_name .= "QC-filtered_id-".$1.".fastq";
			}
			push (@{ $domino_files_threads_QC{$keys}{"FINAL"} }, $final_name);
			File::Copy::move( $domino_files_threads_QC{$keys}{"NGS_QC"}[0], $final_name);
		}			
		&print_Header(" Cleaning Step 2: BLAST and Vector screen ","#"); 
		print "+ Skipping the BLAST step for searching against databases for putative contaminants\n";
	} else {
		&print_Header(" Cleaning Step 2: BLAST and Vector screen ","#"); 
		print "+ Eliminating that reads that could be potential contaminations\n+ Default databases would be used plus the one(s) provided by the user (if any)\n\n\n";
		&print_Header(" Database Generation ","%"); print "+ Generation of Databases: \n";
		print "Please be patient, it could take a while for large files...\n\n";
		my @reads_db;
		if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
			for (my $i=0; $i < 2; $i++) {
				my @name = split("\.fastq_filtered", $domino_files_threads_QC{$keys}{"NGS_QC"}[$i]);
				&extracting_fastq($domino_files_threads_QC{$keys}{"NGS_QC"}[$i], $name[0], "YES");
				push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $name[0].".filtered.fasta");
				push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $name[0].".filtered.qual");
				push (@reads_db, &makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$i]));
		}} elsif ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
			push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $domino_files_threads_QC{$keys}{"NGS_QC"}[0]);
			push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $domino_files_threads_QC{$keys}{"NGS_QC"}[1]);
			push (@reads_db, &makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0]));
		} else {
			my @name = split("\.fastq_filtered", $domino_files_threads_QC{$keys}{"NGS_QC"}[0]);
			&extracting_fastq($domino_files_threads_QC{$keys}{"NGS_QC"}[0], $name[0], "YES");
			push (@{ $domino_files_threads_QC{$keys}{"FASTA4BLAST"} }, $name[0].".filtered.fasta");
			push (@{ $domino_files_threads_QC{$keys}{"QUAL4BLAST"} }, $name[0].".filtered.qual");
			push (@reads_db, &makeblastdb($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0]));
		}
		print "\n"; &time_stamp(); print "\n";
		&debugger_print("domino_files_threads_QC"); &debugger_print("Ref", \%domino_files_threads_QC);

		# Blasting our reads against the databases generated
		print "\n"; &print_Header(" Database Search: BLAST ","%"); 
		my $filter = $BLAST."blastn -query ".$merge_fasta_database." -evalue 1e-10 -db '";
		#Appending databases
		for (my $j=0; $j < scalar @reads_db; $j++) { $filter .= $reads_db[$j]." "; }		
		my $out_file = $taxa_dir{$keys}."/blast_search.txt";
		$filter .= "' -out $out_file -outfmt 6 2> $error_log";
		print "BLASTN command: $filter\n"; my $blastn = system($filter);
		if ($blastn != 0) { &printError("BLASTN failed...\n"); &dieNicely(); } 
		print "\n"; &time_stamp(); print "\n";
		print "+ Parsing the BLAST results for $keys...\n"; 

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
			if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
				$var = &read_FASTA_hash_length($reads_db[0].".fasta_filtered");
			} else { $var = &read_FASTA_hash_length($reads_db[0].".fasta"); }
			my %hash_length = %$var; my $fixed_length;
			if ($hash_length{"UNDEF"}) { $fixed_length = $hash_length{"UNDEF"}; }			
			print "- Reading the BLAST result file $out_file ...\n";
			$/ = "\n";
			my $contamination_ids = $taxa_dir{$keys}."/contaminant_ids.txt";
			open (CONTAMINATION, ">$contamination_ids");	
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
					if ($array[2] > $similarity_thr) { print CONTAMINATION $array[1]."\n"; }
				} else { ## Rest of the reads
					if ($array[2] > $similarity_thr) { 
						if ($fixed_length) {
							my $match_align = $array[3];
							my $percent = (($match_align*100)/$fixed_length);
							if ($percent > $align_len_read) { print CONTAMINATION $array[1]."\n"; }
						} elsif ($hash_length{$array[1]}) {
							my $match_align = $array[3];
							my $percent = (($match_align*100)/$hash_length{$array[1]});
							if ($percent > $align_len_read) { print CONTAMINATION $array[1]."\n"; }
						} else { next;				
			}}}} close(BLAST_result); close (CONTAMINATION);
			
			## TODO: Print contaminants ($array[0]) to a file (cat | uniq ) to show user not the reads but the contaminants found

			## Discard reads
			if (-s $contamination_ids) {
				## Check this paired-end
				my @array = @{$domino_files_threads_QC{$keys}{"FASTA4BLAST"}};
				for (my $j=0; $j < scalar @array; $j++) {
					&mothur_remove_seqs($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$j], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[$j], $taxa_dir{$keys}, $contamination_ids);
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
					my $fastq_filtered = &qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST_filtered"}[$j], $domino_files_threads_QC{$keys}{"QUAL4BLAST_filtered"}[$j], $keys, $final_name);
					push (@{ $domino_files_threads_QC{$keys}{"FINAL"} }, $fastq_filtered);
				} print "\n+ BLAST search and parsing of results for $keys done...\n";		
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
					my $fastq_filtered = &qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[$i], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[$i], $keys, $final_name);
					push (@{ $domino_files_threads_QC{$keys}{"FINAL"} }, $fastq_filtered);
				}
			} else {
				my $final_name;
				my @path_name = split("/", $domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0]);
				if ($path_name[-1] =~ /.*\_id\-(.*)\.clean\_dust.*/) {
					for (my $j=0; $j < scalar $#path_name; $j++) {
						$final_name .= $path_name[$j]."/";
					} $final_name .= "QC-filtered_id-".$1.".fastq";
				}
				my $fastq_filtered = &qualfa2fq_modified_bwa($domino_files_threads_QC{$keys}{"FASTA4BLAST"}[0], $domino_files_threads_QC{$keys}{"QUAL4BLAST"}[0], $keys, $final_name);
				push (@{ $domino_files_threads_QC{$keys}{"FINAL"} }, $fastq_filtered);
	}}} 
	
	print "\n"; &time_stamp(); print "\n";

	my $domino_files_threads_QC_text = $taxa_dir{$keys}."/dumper-hash_QC.txt";
	&print_dump(\%domino_files_threads_QC, $domino_files_threads_QC_text);

	##########################################################################
	## 		Getting files ready for finishing								##
	##########################################################################
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {	
		##adjust pairs?
	} else {
		## What is missing??
	}
	$pm->finish($int_taxa); # pass an exit code to finish
}
$pm->wait_all_children; print "\n** All QC child processes have finished...\n\n";		

## Check each taxa dump file conainting file info
&debugger_print("Retrieve info from files");
foreach my $taxa (keys %taxa_dir) {
	my $dump_file = $taxa_dir{$taxa}."/dumper-hash_QC.txt";
	open (DUMP_IN, "$dump_file");
	while (<DUMP_IN>) {
		my $line = $_;
		chomp $line;
		my @array = split("\t", $line);
		&debugger_print($line);
		push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
	}
	close (DUMP_IN);
}

&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
print "\n";

##############################################################
## 	Check the amount of cleaning done						##
##############################################################
#my $hash_reads_Ref = &checking_discarded();
chdir $dirname;
&create_folder_cleaning_step();
print "\n\n"; &print_Header(" Deleting Temporary Files ","#");

##########################################################################
## 		Deleting temporary files and sorting the output files			##
##########################################################################
print "- Deleting or renaming some files...\n"; 
&delete_temp_file(); print "\n";

##########################################################################
## 	Printing some STATISTICS, generating output folders					##
##########################################################################
&print_Header("","#"); &print_Header(" STATISTICS ","#"); &print_Header("","#");
#&printing_statistics($hash_reads_Ref);  
print "\n Job done succesfully, exiting the script\n";
&finish_time_stamp(); print "\n";
exit(0);


######################################################################################
##																					##	
##									SUBROUTINES										##
##																					##
######################################################################################
sub check_barcode_user_file {
	
	my $user_barcodes_path_file = $_[0];
	unless (-e $user_barcodes_path_file && -r $user_barcodes_path_file) {
		&printError("File $user_barcodes_path_file does not exist or it is not readable\n Please provide a readable file");
		&format_error(); &dieNicely();
	}
	my %species_names_hash;
	my $ROCHE_oligos_file = $dirname."/ROCHE.oligos";
	open (ROCHE, ">$ROCHE_oligos_file");	
	open (BARCODE, "<$user_barcodes_path_file");
	while (<BARCODE>) {
			next if /^#/ || /^\s*$/;
			my $line = $_;
			chomp $line;
			$line =~ s/\t*/\ /g;
			$line =~ s/\s*//g;
			my @barcode_array = split ("\,", $line);
			if (scalar @barcode_array < 2) { &format_error(); &dieNicely(); }
			print ROCHE "barcode\t\t$barcode_array[1]\t\t$barcode_array[3]\n";
			$species_names_hash{$barcode_array[3]} = $barcode_array[2];
	}
	close(BARCODE);

	sub format_error {
		&printError("The only format accepted for the Barcodes file is:
			barcode, ACACGACGACT, MID1, Dmelanogaster                                                                   
			barcode, ACTACGTCTCT, MID2, Dsimulans                                                                  
			barcode, ACTACGTATCT, MID3, Dyakuba
	Please note this is a comma-separated value file (.csv)"); 
	}
	return %species_names_hash;
}

sub check_id_length {
	
	my $file = $_[0];
	my $option = $_[1];
	open (F1, "$file") or &printError("Could not open file $file") and &dieNicely();
	
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
			&debugger_print("Length of ID: $file");
			if ($length_id < 40) {
				&debugger_print("Length < 40 char"); return (1);
			} else {
				&debugger_print("Length > 40 char. Lets rename the files"); 
				return (0);
		}} close (F1);
	} else {	
		$/ = ">"; ## Telling perl where a new line starts
		while (<F1>) {		
			next if /^#/ || /^\s*$/;
			chomp;
	    	my ($titleline, $sequence) = split(/\n/,$_,2);
	    	next unless ($sequence && $titleline);
			my $length_id = length($titleline);
	    	if ($length_id < 40) {
				&debugger_print("Length < 40 char"); return (1);
			} else {
				&debugger_print("Length > 40 char. Lets rename the files"); 
				return (0);
		}} close(F1); $/ = "\n";
}}

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
		$format_returned = &check_file_format($file_to_check);
		my $pair;
		if ($pair_int) {
			if ($pair_int == 1) { $pair = "Left reads file";	
			} elsif ($pair_int == 2) { $pair = "Right reads file";	
			} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\nPair: $pair_int"); &printFormat_message(); &dieNicely();
		}}
		
		if ($format_returned =~ /fastq/) {
			&print_DOMINO_details("\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n\t\tIt also contains an identifier ($name) in the name for later analysis...OK\n");
			push (@{ $domino_files{"original"}{$name} }, $file_to_check); 
			if ($pair) {  &print_DOMINO_details("\t\tPaired end file = ".$pair." ...OK\n");	 }
		} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\n"); &printFormat_message(); &dieNicely();
	}} else { &printError("Please provide a valid FASTQ file: $file_to_check...It is not readable or writable or it does not exist. "); &printFormat_message(); &dieNicely();} 
	&print_DOMINO_details("\n");
	&debugger_print("File: $file_to_check -- $format_returned -- OK;\n");
}

sub check_file_format {    
    my $file = $_[0];
    my $id; my $count = 3;
    my $fasta = my $fastq = my $qual = 0;
    my $format = 'unknown';
    open(FILE, $file) or &printError("Could not open file $file") and &dieNicely();
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

sub check_paired_file_fastq {

	my $file = $_[0];
	my ($tmp_id, $pair, $pair_return, @pair, @types, $type);
	my $count = 0;
	open (F1, "$file") or &printError("Could not open file $file") and &dieNicely();
	while(<F1>) {
		my @Read = ();
		chomp(my $id = $_);
		last if($id=~ /^\n$/);
		## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
		for(my $i=0; $i<3; $i++) { $Read[$i] = <F1>; }
		chomp(my $QualLine = $Read[2]);
		chomp(my $SeqLine = $Read[0]);

		if ($id =~ /(.*)\/(\d+)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/2
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "/1,/2");

		} elsif ($id =~ /(\S*)\s+(\d+)\:N\:\.*/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918 2:N:0:CCGTCC
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "1:,2:");
		} elsif ($id =~ /(.*)\/(\R|\L)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/L
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "L,R");
		}
		if ($count == 5) {last;}
	} close(F1);
	for (my $i = 1; $i < scalar @pair; $i++) {
		if ($pair[0] == $pair[$i]) {
			$pair_return = $pair[$i];
		} else {
			&printError("The reads provided within the file $file do not belong to the same pair..."); &dieNicely();
		}
		if ($types[0] eq $types[$i]) {
			$type = $types[$i];
		} else {
			&printError("The reads provided within the file $file do not belong to the same pair..."); &dieNicely();
		}
	}
	return ($pair_return, $type);
}

sub check_paired_file_fasta {

	my $file = $_[0];
	my ($count, $tmp_id, $pair, $pair_return, @pair, @types, $type);
	open (FILE, "$file") or &printError("Could not open file $file") and &dieNicely();
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($id, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $id);
   		if ($id =~ /(.*)\/(\d+)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/2
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "/1,/2");
		} elsif ($id =~ /(\S*)\s+(\d+)\:N\:\.*/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918 2:N:0:CCGTCC
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "1:,2:");
		} elsif ($id =~ /(.*)\/(\R|\L)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/L
			$tmp_id = $1; $pair = $2; $count++;
			push(@pair, $pair);
			push(@types, "L,R");
		} else {
			print $id."\n";
		}
		if ($count == 5) { last; }
	}
	close(FILE);
	$/ = "\n";
	for (my $i=1; $i < scalar @pair; $i++) {
		if ($pair[0] == $pair[$i]) {
			$pair_return = $pair[$i];
		} else { &printError("Some error ocurred when parsing the paired ends files..."); &dieNicely(); }
		if ($types[0] eq $types[$i]) {
			$type = $types[$i];
		} else { &printError("The reads provided within the file $file do not belong to the same pair..."); &dieNicely();
	}}
	return ($pair_return, $type);
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
			&gunzipping($db_dirname."/".$files[$i]);
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

sub create_folder_cleaning_step {

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
		
	my $domino_files_general = $intermediate_folder."/dumper-hash_DOMINO_files.txt";
	&print_dump(\%domino_files, $domino_files_general);

}

sub debugger_print {
	my $string = $_[0];
	my $ref = $_[1];
	## Print to a file
	if ($debugger) {
		if ($string eq "Ref") {
			print "DEBUG:\t";
			print Dumper $ref; ## if filehandle OUT: print OUT Dumper $ref;
			print "\n";
		} else {
			print "DEBUG:\t".$string."\n";
}}}

sub delete_temp_file {

	##########################################################################################
	##	 																					##
	##  This function deletes the temporary files generated 								##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################

	print "+ Deleting temporary files...\n";	
	my $files_ref = &read_dir($dirname); my @files = @$files_ref;
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
			if ($domino_files{$tags}{'FINAL'}) {
				@array_final = @{$domino_files{$tags}{'FINAL'}};
			} else { 
				@array_final = @{ $domino_files{$tags}{'EXTRACTED'} };
			}
			for (my $i=0; $i < scalar @array_final; $i++) { File::Copy::move($array_final[$i], $clean_folder); }
			
			## Check each taxa folder
			my $dir = $domino_files{$tags}{'DIR'}[0];
			my $files_taxa_ref = &read_dir($dir);
			my @files_taxa_ref = @$files_taxa_ref;
			for (my $f=0; $f < scalar @files_taxa_ref; $f++) {
				if ($files_taxa_ref[$f] eq ".DS_Store" || $files_taxa_ref[$f] eq "." || $files_taxa_ref[$f] eq ".." ) { next;
				} elsif ($files_taxa_ref[$f] =~ /.*_stat/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder."/Quality_Filtering_statistics_".$tags.".txt");
				} elsif ($files_taxa_ref[$f] =~ /.*unPair.*/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} elsif ($files_taxa_ref[$f] =~ /.*singleton.*/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} elsif ($files_taxa_ref[$f] =~ /.*html/) { 
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder."/Cleaning_step_report_".$tags."_QC-Analysis.html"); 		
				} elsif ($files_taxa_ref[$f] =~ /.*png/) {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $intermediate_folder);
				} else {
					File::Copy::move($dir."/".$files_taxa_ref[$f], $temp_dir);
			}}
			remove_tree($dir);
	}}
	unless ($avoidDelTMPfiles) { remove_tree($temp_dir); }
}
	
sub dieNicely { exit(); }

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
	open (F1, "$file") or &printError("Could not open file $file") and &dieNicely();
	while(<F1>) {
		my @Read = ();
		chomp(my $id = $_);
		last if($id=~ /^\n$/);

		## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
		for(my $i=0; $i<3; $i++) { 
			$Read[$i] = <F1>;
		}
		chomp(my $QualLine = $Read[2]);
		chomp(my $SeqLine = $Read[0]);
		
		#print $id."\n".$SeqLine."\n".$QualLine."\n";
		
		## Print into fasta and qual files		
		my @seq_id = split ("\@", $id);
		my $qual_string = &convert_ASCII_to_number($QualLine);
		print FASTA ">".$seq_id[1]."\n".$SeqLine."\n";
		print QUAL ">".$seq_id[1]."\n".$qual_string."\n";
	} close(F1); close (FASTA); close (QUAL);
}

sub finish_time_stamp {
	my $finish_time = time;
	print "\n\n"; &print_Header("","+"); 
	&print_Header(" ANALYSIS FINISHED ","+"); 
	&print_Header("","+");  print "[".(localtime)." ]\t";
	my $secs = $finish_time - $start_time; 
	my $hours = int($secs/3600); $secs %= 3600; 	
	my $mins = int($secs/60); $secs %= 60; 
	printf ("Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs);
}

sub gunzipping {

	##########################################################################################
	##	 																					##
	##  This function unzips the files zipped. Only valid for UNIX and if gunzip installed	##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $gunzip_file = $_[0];
	my $process;
	if ($gunzip_file =~ /.*tar.*/) {
		my $tar_gz = "tar -zxvf $gunzip_file";
		$process = system($tar_gz);
	} else {
		my $gunzip = "gunzip -d ".$gunzip_file;
		$process = system($gunzip);
	}
	if ($process != 0) {
		&printError("Gunzipping the file failed when trying to process the file...\n");
		&printError("DOMINO would not die here...\n");
}}
 
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

	for (my $i = 0; $i < scalar @blast_db; $i++) { &print_DOMINO_details("\t\t- $blast_db[$i]\n"); }
	sub print_file{
		my $element = $_;
		if(-f $element && $element =~ /.fa/){
			push (@blast_db, "$File::Find::name");
}}}

sub makeblastdb {

	##########################################################################################
	##	 																					##
	##  This function makes a valid BLAST db for each file given using makeblastdb			##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
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
	return $db;
}

sub mothur_sffinfo {
	
	###########################################################################################
	##	 																					 ##
	##  This function uses the MOTHUR standalone version to extract and trim tags of the 	 ##
	##	reads of a "sff" file format provided.												 ##
	##	Jose Fco. Sanchez Herrero, 28/04/2014 jfsanchezherrero@ub.edu						 ##
	## 																						 ##
	###########################################################################################
	
	# This subroutine takes as input a sff file and extracts it, classifies according to the tags provided, by Roche and trims the seqs
	my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
	my $file = $_[0];
	my $file_name = $_[1];
	my $directory = $dirname;

	## The '#' is necessary to be in the first place when calling mothur with a given set of commands
	$mothur_path .= " '#set.dir(output=$directory); sffinfo(sff='".$file."'); trim.seqs(fasta="."$file_name".".fasta, qfile=".$file_name.".qual, oligos=ROCHE.oligos, bdiffs=$bdiffs, processors=$noOfProcesses)'";
	my $mothur_result = system($mothur_path);
	if ($mothur_result != 0) { &printError("MOTHUR failed when trying to proccess the file..."); &dieNicely(); } else {
		print "Done...OK\n\n";
	}	
}

sub mothur_trim_fastq {
	
	# This subroutine takes as input a FASTQ file AND classifies according to the tags provided, by Roche and trims the seqs
	my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";
	my $file_name = $_[0];
	my $directory = $dirname;

	## The '#' is necessary to be in the first place when calling mothur with a given set of commands
	$mothur_path .= " '#set.dir(output=".$directory."); trim.seqs(fasta=".$file_name.".fasta, qfile=".$file_name.".qual, oligos=ROCHE.oligos, bdiffs=".$bdiffs.", processors=".$noOfProcesses.")'";
	print "+ Triming the reads according to MID tag\n$mothur_path\n";

	my $mothur_result = system($mothur_path);
	if ($mothur_result != 0) { 	&printError("MOTHUR failed when trying to proccess the file..."); &dieNicely(); } else { print "Done...OK\n\n"; }
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

sub prinseq_dust {
	
	###########################################################################################
	##	 																					 ##
	##  This function uses the script prinseq-lite.pl distributed in PRINSEQ 				 ##
	##  It would process reads containing low complexity under a threshold given             ##
	##	and the reads containing ambiguous bases and N's								     ##
	## 																						 ##
	###########################################################################################
	
	my $file1 = $_[0];
	my $file2 = $_[1];
	my $thr = $_[2];
	my $dir = $_[3];
	
	my $prinseq_dust_command = "perl ".$scripts_path."PRINSEQ-lite_0.20.4/prinseq-lite.pl";
	my @name = split("/", $file1);
	my @id = split(".fastq", $name[-1]); 
	my $out_good = $dir."/".$id[0].".clean_dust";
	my $out_bad = $dir."/".$id[0].".dust";
	
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		$prinseq_dust_command .= " -fastq ".$file1." -fastq2 ".$file2;
	} else {
		$prinseq_dust_command .= " -fastq ".$file1;
	}

	if ($file_type eq "454_fastq" || $file_type eq "454_multiple_fastq" || $file_type eq "454_sff" ) {
		$prinseq_dust_command .= " -out_format 2"; ## Generate fasta +  qual for NGSQC step
	}

	$prinseq_dust_command .= " -out_good ".$out_good;
	$prinseq_dust_command .= " -out_bad ".$out_bad;
	$prinseq_dust_command .= " -lc_method dust -lc_threshold ".$thr; 
	$prinseq_dust_command .= " -ns_max_n 1 -noniupac"; 
	$prinseq_dust_command .= " 2> $error_log";

	print "PRINSEQ command:\n".$prinseq_dust_command."\n";

	my $dust_result = system ($prinseq_dust_command);	
	if ($dust_result != 0) { &printError("DUST cleaning step failed when trying to proccess the file... DOMINO would not stop in this step...");  }
	return ($out_good, $out_bad);	
}

sub print_DOMINO_details {
	
	my $string = $_[0];
	open (PARAM, ">>$param_Detail_file");
	print PARAM $string;
	print $string;
	close(PARAM);
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
	print STDERR $message."\n";
}

sub printFormat_message {
	print "\n\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\nWhere:\n\txxx: any character or none. Please avoid using dots (.)\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n\n";
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

sub print_dump {
	my $ref = $_[0]; 
	my $out_file = $_[1];
	my %tmp_hash = %$ref;
	
	open (DUMP, ">>$out_file"); 
	&debugger_print("Printing to a file $out_file");
	foreach my $names (sort keys %tmp_hash) {
		my %hash = %{$tmp_hash{$names}};
		foreach my $tags (sort keys %hash) {
			my @array = @{ $hash{$tags} };
			for (my $i=0; $i < scalar @array; $i++) {
				print DUMP $names."\t".$tags."\t".$hash{$tags}[$i]."\n";
				&debugger_print($names."\t".$tags."\t".$hash{$tags}[$i]);
			}
		}
	} close (DUMP);
}
 
sub qualfa2fq_modified_bwa {

	##########################################################################################
	##  This is a modified version of the perl script provided in the bwa-0.7.8 package		##	
	## 	to convert fasta+qual files into fastq files										##
	##	Jose Fco. Sanchez Herrero, 06/05/2014 jfsanchezherrero@ub.edu						##
	##########################################################################################
	
	my ($fhs, $fhq, $q, $a);
	my $taxa_name = $_[2];
	my $name = $_[3];
	my @temp_name = split ("fasta", $_[0]);
	my $fastq_file; my $bool_change_name = 0; my $pair;

	my $codeReturn = &check_id_length($_[0], "fasta");
	my ($pairReturn, $type);
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		($pairReturn, $type) = &check_paired_file_fasta($_[0]);
		&debugger_print("Pair: ".$pairReturn."\t"."Type: ".$type);
		if ($pairReturn eq "L" ) { 
			$codeReturn = 0; ## Force to rewrite seq					
			&debugger_print("Lets change id $_[0] : pair was $pairReturn");
		} elsif ($pairReturn eq "R") { 
			$codeReturn = 0; ## Force to rewrite seq
			&debugger_print("Lets change id $_[0] : pair was $pairReturn");
		} # else {
			## There is no need to change file according to pair
	}
	if ($codeReturn == 0) { $bool_change_name = 1; } 
	
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		($pair, $type) = &check_paired_file_fasta($_[0]);
		&debugger_print("Pair: ".$pair."\t"."Type: ".$type);
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
 
sub read_FASTQ_ids_file {
	
	my $file = $_[0];
	my %hash;		
	open (F1, "$file") or &printError("Could not open file $file") and &dieNicely();
	while(<F1>) {
		my @Read = ();
		chomp(my $id = $_);
		last if($id=~ /^\n$/);

		## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
		for(my $i=0; $i<3; $i++) { 
			$Read[$i] = <F1>;
		}
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

sub read_FASTA_hash {

	my $file = $_[0];
	my %hash;
	open(FILE, $file) || die "Could not open the $file ...\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	$hash{$titleline} = $sequence;
	}
	close(FILE);
	$/ = "\n";
	my $hashRef = \%hash;
	return $hashRef;
}

sub read_FASTA_hash_length {

	my $file = $_[0];
	my %hash;
	my $counter;
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
			&debugger_print("Fixed size for reads in $file -- $length");
			my %return = ("UNDEF" => $length); 
			return \%return; 
		} 
	}
	close(FILE); $/ = "\n"; 
	&debugger_print("Different sequence lengths");
	&debugger_print("Ref", \%hash);
	return \%hash;
}

sub read_dir {
	my $dir = $_[0];
	opendir(DIR, $dir) or &printError("Can not open folder $dir...") and &dieNicely;
	my @dir_files = readdir(DIR);
	my $array_ref = \@dir_files;
	return $array_ref;
}

sub rename_seqs {

	my $file = $_[0];
	my $id_Taxa = $_[1];
	my $typeFile = $_[2];
	
	my @file_name = split("/", $file);
	my $dir = $dirname."/".$id_Taxa;
	unless (-d $dir) { mkdir $dir, 0755; }
	my $file_path = $dir."/renamed_".$file_name[-1];
	
	my $codeReturn = &check_id_length($file, "fastq");

	my ($pairReturn, $type);
	if ($file_type eq "Illumina_pair_end" || $file_type eq "Illumina_pair_end_multiple_fastq") {
		if ($typeFile eq 'fastq' ) {
			($pairReturn, $type) = &check_paired_file_fastq($_[0]);
		} else {
			($pairReturn, $type) = &check_paired_file_fasta($_[0]);
		}
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
		
		open (F1, "$file") or &printError("Could not open file $file") and &dieNicely();
		open (FASTQ_out, ">$file_path");
		my $int = 0;
		while(<F1>) {
			my @Read = ();
			chomp(my $id = $_);
			last if($id=~ /^\n$/);

			## Read each entry of the FASTQ file; 1) @Seq_id 2) Sequence 3) +Seq_id  4) Qual
			for(my $i=0; $i<3; $i++) { 
				$Read[$i] = <F1>;
			}
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

sub seq_counter {
	
	my $file = $_[0];	
	my $option_format = &check_file_format($file);
	my ($nSequences, $nLines);
	open (F1, "$file") or &printError("Could not open file $file") and &dieNicely(); 
	while (<F1>) { $nLines++; } close(F1);
	if ($option_format eq "fastq") { 
		$nSequences = int ($nLines/4);
	} elsif ($option_format eq "fasta") {
		$nSequences = int ($nLines/2);
	}
	return $nSequences;
}

sub splitting_fastq {

	my $file_name = $_[0];
	my $int = $_[1];
	my $fasta_file = $file_name.".trim.fasta";
	my $qual_file = $file_name.".trim.qual";
	my $group_file = $dirname."/".$file_name.".groups";
	my %files_generated;

	$/ = "\n"; ## Telling Perl where a new line starts
	open (GROUPS, $group_file) or &printError("Cannot open file groups: $group_file") and &dieNicely;
	my %files;
	# Generate a hash with the identifier for each sequence
	while (<GROUPS>) {
		chomp; my $line = $_;
		$line =~ s/\s+/\t/g; $line =~ s/\t+/\t/g;
		my @array = split ("\t", $line);
		my $filehandle;
		if ($files{$array[1]}) {
			$filehandle = $files{$array[1]}
		} else {
			$filehandle = $dirname."/".$file_name."-".$array[1]."_ids.txt";
			$files{$array[1]} = $filehandle;
		}
		open (OUT, ">>$filehandle"); print OUT $array[0]."\n"; close(OUT);
	} close (GROUPS);
	&debugger_print("Taxa ids reads:\n"); &debugger_print("Ref", \%files);
	
	my $pm_mothur =  new Parallel::ForkManager($noOfProcesses);
	$pm_mothur->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
		print "\n\n** Child process finished for $ident with PID $pid and exit code: $exit_code\n\n"; 
	} );
	$pm_mothur->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** Mothur retrieving sequences started for $ident with PID $pid\n\n"; } );
	
	foreach my $files_ids (sort keys %files) {
		my $pid = $pm_mothur->start($files_ids) and next;
		print "\nSending child process for retriving reads for $files_ids\n\n";
		my $dir = $dirname."/".$files_ids; unless (-d $dir) { mkdir $dir, 0755; }
		&mothur_retrieve_seqs($fasta_file, $qual_file, $dir, $files{$files_ids});

		push (@{ $files_generated{$files_ids}{'reads_ids'} }, $files{$files_ids});
		push (@{ $files_generated{$files_ids}{'FASTA_extracted'} }, $dir."/".$file_name.".trim.pick.fasta");
		push (@{ $files_generated{$files_ids}{'QUAL_extracted'} }, $dir."/".$file_name.".trim.pick.qual");			
		my $name;
		if ($file_type eq "Illumina_pair_end") {
			my $pair; my $type;
			($pair, $type) = &check_paired_file_fasta($dir."/".$file_name.".trim.pick.fasta");
			$name = $dir."/reads_id-".$files_ids."_R".$pair.".fastq";
		} else { $name = $dir."/reads_id-".$files_ids.".fastq"; }

		my $fastq_file = &qualfa2fq_modified_bwa($dir."/".$file_name.".trim.pick.fasta", $dir."/".$file_name.".trim.pick.qual", $files_ids, $name);
		push (@{ $files_generated{$files_ids}{'EXTRACTED'} }, $fastq_file);			
		my $extracted = $dir."/dumper_extracted.txt";
		&debugger_print("####### Extracted dumper ".$extracted); &print_dump(\%files_generated, $extracted);
		$pm_mothur->finish($int_taxa); # pass an exit code to finish
	}
	$pm_mothur->wait_all_children;
	print "\n** All Mothur retrieve processes have finished...\n\n";	
}

sub time_stamp {	
	my $current_time = time;
	print "[ ".(localtime)." ]\t";
	my $secs = $current_time - $step_time; 
	my $hours = int($secs/3600); $secs %= 3600; 
	my $mins = int($secs/60); $secs %= 60; 
	$step_time = $current_time;
	printf ("Step took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

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
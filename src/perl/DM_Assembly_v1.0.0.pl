#!/usr/bin/perl
###########################################################################################
### DOMINO: Development of molecular markers in non-model organisms using NGS data 	###
###											###
### Authors:										###
### Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro 		###
### Sánchez-Gracia, and Julio Rozas.					     		###
###											###
###########################################################################################
##	Usage:
##      perl DM_Assembly_v1.0.0.pl
##
##    ###########################
##    ### General Information ###
##    ###########################
##      [-h|--help] [-man] [-v|--version] [-MoreInfo]
##
##    #########################
##    ### Mandatory options ###
##    #########################
##      [-type_file int_value] [-DOMINO_files || -user_files]
##      [-o|--outputFolder path]
##
##    ######################
##    ### Optional flags ###
##    ######################
##      ## CAP3
##
##      [-use_CAP3] [-overCAP3|--overalpping_CAP3 int_value]
##      [-simCAP3|--similarity_CAP3 int_value]
##
##      ## Others
##
##      [-p|--processes int_value] [-mrs|--min_relative_score int_value]
##      [-input_files file] [-TempFiles]
##
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use POSIX qw(strftime);
use FindBin;
use lib $FindBin::Bin."/lib";

BEGIN {
	require List::Uniq;
	use List::Uniq qw(uniq);
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
## TODO:: use debugger_print for printing Debugging messages

##################################
##	Initializing some variables	##
##################################
my ($helpAsked, $job, $type_job, $reference_fasta_file, $hirep, @fastq_files, 
$avoidDelTMPfiles, $file_type, $cap3flag, %cap3files, %MID_species_hash, @user_files, 
$noOfProcesses, $version, %files_used, $manual, $noClip, $abs_folder, $mrs, $debugger, $flagSpades,
@CAP3_directories, $DOMINO_files, $overlap_CAP3, $similar_CAP3, $user_files, $further_information,
@file_abs_path, @file_names, $step_time, %nucleotides, $total_Contigs_all_sets);

my %bases = ("A" => 0, "C" => 0, "T" => 0, "G" => 0, "N" => 0);
my %input_options = (1 =>'454_sff', 2 =>'454_fastq', 3=>'454_multiple_fastq', 
	4 => 'Illumina', 5 => 'Illumina_multiple_fastq', 6 => 'Illumina_pair_end', 
	7 => 'Illumina_pair_end_multiple_fastq');

######################
## Get user options ##
######################
GetOptions(
	"h|help" => \$helpAsked,
	"man" => \$manual,
	"v|version" => \$version, 
	"MoreInfo" => \$further_information,
	
	"type_file=i" => \$file_type,
	"DOMINO_files" => \$DOMINO_files,

	"user_files" => \$user_files,
	"input_files=s" => \@user_files,

	"o|outputFolder=s" => \$abs_folder,
	"mrs|min_relative_score=i" => \$mrs,
	
	"use_CAP3" => \$cap3flag,
	"overCAP3|overalpping_CAP3=i" => \$overlap_CAP3,
	"simCAP3|similarity_CAP3=i" => \$similar_CAP3,
	
	"TempFiles" => \$avoidDelTMPfiles,
	"p|number_cpu=i" => \$noOfProcesses,
	
	"SPAdes" => \$flagSpades,
	
	"Debug" => \$debugger,	
);

## Help/manual/version asked
pod2usage( -exitstatus => 0, -verbose => 1 ) if ($helpAsked);
pod2usage( -exitstatus => 0, -verbose => 2 ) if ($manual);
pod2usage( -exitstatus => 0, -verbose => 99, -sections => "VERSION") if ($version);
pod2usage( -exitstatus => 0, -verbose => 99, -sections => "NAME|VERSION|DESCRIPTION|AUTHOR|COPYRIGHT|LICENSE|DATE|CITATION") if ($further_information);

=pod

=over 2

=item B<>

=item B<###############################################################>

=item B<################ DOMINO Assembly pipeline #####################>

=item B<###############################################################>

=back

=head1 NAME

=over 2

DM_Assembly_v1.0.0.pl  

=back	

=head1 VERSION

=over 2

DOMINO v1.0.0

=back	
	
=head1 SYNOPSIS

=over 2

=item B<>

perl DM_Assembly_v1.0.0.pl  

=item B<###########################>
	
=item B<### General Information ###>

=item B<###########################>

[-h|--help] [-man] [-v|--version] [-MoreInfo]

=item B<#########################>

=item B<### Mandatory options ###>

=item B<#########################>

[-type_file int_value] [-DOMINO_files || -user_files] [-o|--outputFolder path]

=item B<######################>

=item B<###	Optional flags ###>

=item B<######################>

## CAP3

[-use_CAP3] [-overCAP3|--overalpping_CAP3 int_value] [-simCAP3|--similarity_CAP3 int_value]

## Others

[-p|--processes int_value] [-mrs|--min_relative_score int_value] [-input_files file] [-TempFiles] 

=item B<>

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

=item B<#######################################>

=item B<######## DOMINO Assembly phase ########>

=item B<#######################################>

=item B<>
	
DOMINO performs separate assemblies, one for each panel taxon, using MIRA v4.0.2 with the pre-processed reads from the previous step or with those supplied by the user. 

Although the default parameter values vary in function of the particular sequencing technology, the majority of them are shared (see the DOMINO manual). 

In order to avoid including repetitive and chimeric regions, all contigs (and the corresponding reads) identified as HAF6, HAF7 and MNRr by the MIRA algorithm are discarded from the mapping/alignment phase. 

Since MIRA can generate redundant contigs because of polymorphic and paralogous regions, we have implemented a specific DOMINO function that performs a similarity clustering of all contigs (i.e., an all versus all contigs BLAST search) to identify and remove such redundancies. All positive hits (E-value < 10-50), with an aligned region length higher or equal than the 85% of the shorter contig, and a minimum percentage of similarity of the 85%, are collapsed to the longest contig of the pair for the mapping phase. 

The DOMINO command line version also includes an option to perform a second iterative assembly step using the software CAP3. If selected, this option uses MIRA output sequences (contigs and singletons) as input for CAP3 under a relaxed parameter scheme: overlapping consensus similarity of 80 bp and percent identity of 98%. Users can modify these default parameters in the command prompt. 
	
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

=item B<################ INPUT ##############>

=item B<#####################################>

=item B<>

=item B<-type_file [string]>

Several types of input files are accepted. Select only one option:

3: Roche 454 reads multiple FASTQ files, each file containing the reads from one taxon. 

5: Multiple Illumina (single-end) FASTQ files, each file containing the reads from one taxon. 

7: Multiple Illumina paired-end reads FASTQ files, each pair from the same taxon containing left and right reads, respectively. 

=item B<>

=item B<-DOMINO_files || -user_files>

=item B<>

B<-DOMINO_files [Default]>

Default option. DOMINO uses the pre-processed reads in /DM_clean_data folder for the assembly.

If QC is not desired, please prepare the files using '-only_tag_files' option using the DM_Clean_v1.0.0.pl script or provide appropriately tagged files.

B<-user_files>

Use this option jointly with the flag '-input_files' to provide pre-processed reads (454 Roche or Illumina (single/pair-end) other than the generated in the previous phase. 

=item B<-input_files [files]> 

Provide input files to run DOMINO starting from this point.

Provide one input file per taxa named as "[xxid-][yyy][_Rn].fastq" or "[yyy][_Rn].fastq". Where:

'[xxid-]' might be present or not. [xx] could be any character (or none). Please avoid using dots (.)
	
'[yyy]' taxon identifier. [Mandatory]

'[_Rn]' If paired-end data, R1 or R2, for the left and the right reads, respectively.

Single End:

Example: C<-input_files path_to_file1/reads.id-Dmelanogaster.fastq -input_files path_to_file2/reads.id-Dsimulans.fastq -input_files path_to_file3/reads.id-Dyakuba.fastq>
	
Paired-End:

Example C<-input_files path_file1_left/reads.id-Dmelanogaster_R1.fastq -input_files path_file1_right/reads.id-Dmelanogaster_R2.fastq -input_files path_file2_left/reads.id-Dsimulans_R1.fastq -input_files path_file2_right/reads.id-Dsimulans_R2.fastq -input_files path_file3_left/reads.id-Dyakuba_R1.fastq -input_files path_file3_right/reads.id-Dyakuba_R2.fastq>

=item B<>

=item B<#####################################>

=item B<############# OUTPUT ################>

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

=item B<-mrs|--min_relative_score [int_value]>

Minimum % of matching between two reads to be considered for the assembly (MIRA assembly package). Increasing this value will save memory, but one might lose sensitivity. [MIRA Default: 80]

=item B<-overCAP3|--overalpping_CAP3 [int_value]> 

Minimum Overlap length cutoff (>15) for CAP3 assembly. [Default: 80]

=item B<-simCAP3|--similarity_CAP3 [int_value]>

Minimum Overlap percent identity cutoff (>65) for CAP3 assembly. [Default: 97]

=item B<>

=item B<#####################################>

=item B<############### OPTIONS #############>

=item B<#####################################>

=item B<>

=item B<-use_CAP3> 

Use CAP3 for a second assembly round. [Default: Off].

=item B<>

=item B<-TempFiles> 

Keep all intermediate files.

=item B<-p|--number_cpu [int_value]>

Number of threads/cores to be used. [Default: 2]

=item B<>

=item B<--SPAdes>

Use SPAdes Genome Assembler for a better assembly. It is mandatory a Linux server and at least 30GiB of free Memory RAM.

=item B<>

=item B<#####################################>

=item B<##### Command Line Examples #########>

=item B<#####################################>

=item B<>

=item B<454: DOMINO Multiple fastq files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -DOMINO_files -type_file 3 -p 3 -mrs 80

=item B<454: DOMINO Multiple fastq files. Use CAP3 for scaffolding>

 perl DM_Assembly_v1.0.0.pl -o test_folder -DOMINO_files -type_file 3 -p 3 -mrs 80 
 -useCAP3 -overCAP3 80 -simCAP3 97 -TempFiles

=item B<Illumina: DOMINO Multiple fastq files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -DOMINO_files -type_file 5

=item B<Illumina Paired-end Files: DOMINO Multiple Paired-end FASTQ files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -DOMINO_files -type_file 7

=item B<454: User Multiple fastq files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -type_file 3 -user_files 
 -input_file 454_user_file1.fastq -input_file 454_user_file2.fastq 
 -input_file 454_user_file3.fastq

=item B<Illumina: User Multiple fastq files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -type_file 5 -user_files 
 -input_file illumina_user_file1.fastq -input_file illumina_user_file2.fastq 
 -input_file illumina_user_file3.fastq

=item B<Illumina Paired-end Files: User Multiple Paired-end FASTQ files>

 perl DM_Assembly_v1.0.0.pl -o test_folder -type_file 7 -user_files 
 -input_file user_file_left_1.fastq -input_file user_file_right_1.fastq 
 -input_file user_file_left_2.fastq -input_file user_file_right_2.fastq

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

###########################
## Checking user options ##
###########################
pod2usage( -exitstatus => 0, -verbose => 1 ) if ($helpAsked);
pod2usage( -exitstatus => 0, -verbose => 2 ) if ($manual);
if (!$abs_folder || !$file_type) {
	Pod::Usage::pod2usage(-exitstatus => 0, -verbose => 0 );
} elsif ($user_files) {
	if (!$file_type) {
		print "Some parameters are missing\tPlease provide type of file\n"; &dieNicely();
}}

## Using default options if not provided
if (!$mrs) { $mrs = 80; }
if (!$overlap_CAP3) { $overlap_CAP3 = 80; }
if (!$similar_CAP3) { $similar_CAP3 = 97; }
if (!$noOfProcesses) { $noOfProcesses = 2; }

## Getting some PATH variables
my $pipeline_path = abs_path($0); ## abs_path($0) es donde esta el script!!!
my $path_abs_folder = abs_path($abs_folder);
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; }
my $MIRA_exec = $scripts_path."mira_v4.0/bin/mira";
my $CAP3_exec = $scripts_path."cap3/bin/cap3";
my $BLAST = $scripts_path."NCBI_BLAST_v2.2.28/";
my $mothur_path = $scripts_path."MOTHUR_v1.32.0/mothur";

## Checking if the Directory already exists because a previous analysis
unless (-d $path_abs_folder) { mkdir $path_abs_folder; } 
my $random_number = int(rand(1000));
my $datestring = strftime "%Y%m%d%H%M", localtime;

my $dirname = $path_abs_folder."/".$datestring."_DM_assembly";
my $dirname_tmp = $dirname."/".$datestring."_DM_intermediate_assembly";
my $error_log = $dirname."/".$datestring."_DM_Assembly_ERROR.txt";
my $param_Detail_file = $dirname."/".$datestring."_DM_Assembly_Parameters.txt";
if (-e $param_Detail_file) {
	File::Copy::move($param_Detail_file, $param_Detail_file."_old_".$random_number);
} elsif (-e $error_log) {
	File::Copy::move($error_log, $error_log."_old_".$random_number);
}	
																							  
if ($user_files) {
	## Check files absolute path and name
	for (my $i = 0; $i < scalar @user_files; $i++) {
		my $tmp = abs_path($user_files[$i]);
		push (@file_abs_path, $tmp);
		my @tmp_file_path = split ("/", $tmp);
		push (@file_names, $tmp_file_path[-1]);		
}}

## Generate directories
if (-d $dirname) { File::Copy::move($dirname, $dirname."_old_".$random_number); }
mkdir $dirname, 0755; chdir $dirname;
if (-d $dirname_tmp) { File::Copy::move($dirname_tmp, $dirname_tmp."_old_".$random_number); }
mkdir $dirname_tmp, 0755;

## Start the Analysis
print "\n"; &print_Header("","+");  &print_Header(" Analysis Started ","+"); &print_Header("","+"); print "\n";
my $start_time = $step_time = time;
&print_DOMINO_details("Starting the process: [ ".(localtime)." ]\n\n");
# Checking user options
&print_Header(" Input File and Parameter Preprocessing ","#");
&print_DOMINO_details("\n+ Output Directory: ".$dirname." ....OK\n");

## If using SPADES for assembly instead of MIRA
my ($assembly_directory, $assembly_directory_abs_path, $CAP3_directory, $CAP3_directory_abs_path);
if ($flagSpades) {
	$assembly_directory = "spades_assemblies";
} else {
	my $assembly_directory = "MIRA_assemblies"; 
	if ($cap3flag) {
		$CAP3_directory = "CAP3_assemblies"; mkdir $CAP3_directory, 0755;
		$CAP3_directory_abs_path = abs_path($CAP3_directory); 
}}
$assembly_directory_abs_path = abs_path($assembly_directory);
mkdir $assembly_directory_abs_path, 0755; chdir $assembly_directory_abs_path;

##################################
##	Checking files 		##
##################################
if ($input_options{$file_type}) {
	if ($file_type == 1 || $file_type == 2 || $file_type == 4|| $file_type == 6) {
		&printError("Input file type provided not supported, please read the documentation..."); &dieNicely();	
	}
} else {
	&printError("\n\nERROR: Wrong type of file provided\nPlease provide a valid type of file:\n");
	&printError("According to the type of NGS files provided several options are available but only one option would be provided: -type_input [int]
 1: A single file in Standard Flowgram Format (SFF), 454 file, containing all the reads of the different taxa accordingly tagged
 2: FASTQ files coming from 454 Roche. A single file containing all the reads of the different taxa accordingly tagged
 3: Multiple FASTQ files coming from 454 Roche. Each file contains each taxa reads. 
 4: FASTQ file from Illumina single end, containing all the reads of the different taxa accordingly tagged
 5: Multiple FASTQ file from Illumina single end. Each file contains each taxa reads. 
 6: A single pair of FASTQ files from Illumina paired-end sequencing: Each pair would contain all the reads of the different taxa accordingly tagged
 7: Multiple FASTQ files from Illumina paired-end: Each pair of files containing each taxa left and right reads respectively
 
 During the assembly stage, only separate files for each taxa is valid, so only options 3, 5 and 7 could be provided\n\n");
	Pod::Usage::pod2usage("Try 'perl $0 -man' for more information\n");	exit();
}

if ($user_files) { ## If users provides files for the assembly either than the DOMINO_clean_data
	if (scalar (@user_files) == 0 || !$file_type) {
		&printError("No input files provided"); &dieNicely();		
	} else {
		## Checking files
		&print_DOMINO_details("+ Type of file(s): Option $file_type : $input_options{$file_type} ...OK\n");
		&print_DOMINO_details("+ User files option provided: ....OK\n");
		&print_DOMINO_details("+ Checking file(s) user provided:\n");
		
		if ($file_type == 1) {
			&printError("Please provide a 454 Roche SFF file extracted into multiple FASTQ files tagged with MID ids like 'Seq_56-328-1_MID1_sp1'"); &dieNicely();
		} elsif ($file_type == 2) {
			&printError("Please provide a multiple FASTQ files containing each taxa read accordingly tagged with MID ids like 'Seq_56-328-1_MID1_sp1'"); &dieNicely();
		} elsif ($file_type == 3) {
			if (scalar (@file_abs_path) == 1) {
				&printError("Please provide multiple 454 Roche FASTQ file containing each of the taxa reads"); &dieNicely();
			} else {
				for (my $i = 0; $i < scalar @file_abs_path; $i++) {
					&check_file($file_abs_path[$i]);
					system("ln -s $file_abs_path[$i]");
					push (@fastq_files, $file_names[$i]);
			}}
		} elsif ($file_type == 4) {
			&printError("Please provide multiple Illumina FASTQ files containing each taxa reads accordingly tagged with MID ids like 'Seq_56-328-1_MID1_sp1'"); &dieNicely();
		} elsif ($file_type == 5) { 
			if (scalar (@file_abs_path) == 1) {
				&printError("Please provide multiple Illumina FASTQ files containing containing each of the taxa reads"); &dieNicely();
			} elsif (scalar (@file_abs_path) != 1) {
				for (my $i = 0; $i < scalar @file_abs_path; $i++) {
					&check_file($file_abs_path[$i]);
					system("ln -s $file_abs_path[$i]");
					push (@fastq_files, $file_names[$i]);
		}}} elsif ($file_type == 6) {
			&printError("type_file == 6 is not valid. Please provide multiple Illumina Paired-End FASTQ files for each taxa (type_file 7) (2 files for specie)");
			&printFormat_message(); &dieNicely();
		} elsif ($file_type == 7) {  ## Specify user to input left-right for each taxa: -input sp1_left.fastq -input sp1_right.fastq -input sp2_left.fastq -input sp2_right.fastq
			if (scalar (@file_abs_path) == 2 ) {
				&printError("type_file == 6 is not valid. Please provide multiple Illumina Paired-End FASTQ files for each taxa (type_file 7) (2 files for specie)");
				&printFormat_message(); &dieNicely();
			} elsif (scalar (@file_abs_path) != 2) {
				for (my $i = 0; $i < scalar @file_abs_path; $i++) {
					&check_file($file_abs_path[$i]);
					system("ln -s $file_abs_path[$i]");
					push (@fastq_files, $file_names[$i]);
		}}} else { ## Option provided is not recognized or mispelled
			&printError("Please provide a valid file type option"); &dieNicely();
}}} else {
	## Go to DOMINO_clean_data and get the files
	&print_DOMINO_details("+ Type of file(s): Option $file_type : $input_options{$file_type} ...OK\n+ No user files provided: Default DOMINO cleaning files ...OK\n+ Get the clean FASTQ files generated in the cleaning step ...OK\n");
	&get_fastq_files();	
	for (my $i = 0; $i < scalar @file_abs_path; $i++) { &check_file($file_abs_path[$i]); }
}
print "\n"; &time_stamp();
if (scalar @file_abs_path == 0) { &printError("No clean files were provided for the assembly...\n"); &printFormat_message(); &dieNicely();}

## Check if paired end files provided, all pairs have been correctly provided
if ($file_type == 7) {
	foreach my $keys (keys %MID_species_hash) {
		my $sum_pairs; my $files;
		for (my $i=0; $i < scalar @{ $MID_species_hash{$keys} }; $i++) {
			my $string2split = $MID_species_hash{$keys}[$i];
			my @array = split("----", $string2split);
			my $pair = $array[0]; my $file = $array[1];
			$files .= $file."\n"; $sum_pairs += $pair;
		}
		if ($sum_pairs ne 3) { &printError("Wrong pair of FASTQ files provided. Please check pair of files corresponding to $keys:\n$files\n"); &dieNicely(); }
}}

###########################
## Printing user options ##
###########################
## Print information of the input parameters
&print_DOMINO_details("+ Threads: ".$noOfProcesses." ...OK\n");
&print_DOMINO_details("+ Minimum relative Score (mrs) to use in MIRA assembly: $mrs ...OK\n");
if ($cap3flag) {
	&print_DOMINO_details("+ CAP3 would be used in a second round of assembly ...OK\n");
	&print_DOMINO_details("+ Similarity between reads for CAP3 assembly: $similar_CAP3 ...OK\n");
	&print_DOMINO_details("+ Overlapping between reads for CAP3 assembly: $overlap_CAP3 ...OK\n");
} else { 
if ($flagSpades) { &print_DOMINO_details("+ MIRA has been disabled, SPAdes would perform the assembly ...OK\n");	
} else { &print_DOMINO_details("+ CAP3 has been disabled, only MIRA would perform the assembly ...OK\n");
}}
unless ($avoidDelTMPfiles) {
	&print_DOMINO_details("+ Deleting of temporary files would be done ...OK\n");
} else { &print_DOMINO_details("+ Deleting temporary files would be avoid ...OK\n"); }

## Print info about where to print info and error
&print_DOMINO_details("\n+ Parameters details would be print into file: $param_Detail_file...\n");
&print_DOMINO_details("+ Errors occurred during the process would be print into file: $error_log...\n\n");

## Assembly
if ($flagSpades) {

	## If user provides Spades flag, we will assume, in this version, he is using a Linux server with
	## enough CPUs and RAM for assembly	
	## We would use "cat /proc/meminfo " in order to check the available memory, if greater than 30GiB
	## we would split in a reasonable amount of processes to proceed in parallel with the assembly 
	
	# Get memory
	system("cat /proc/meminfo > memory_server.txt");	
	my $memory_server_file = $assembly_directory_abs_path."/memory_server.txt";
	my $total_available;
	if (-e -r -s $memory_server_file) {
		&print_Header(" Memory RAM Usage retrieval ","#");
		my %memory_hash;
		open (MEM, $memory_server_file);
		while (<MEM>) {
			my $line = $_;
			chomp $line;
			$line =~ s/\s*//g;			
			&debugger_print($line);			
			my @array = split("\:", $line);
			if ($array[1] =~ /(\d+)(.*)/) {
				push (@{ $memory_hash{$array[0]} }, $1);
				push (@{ $memory_hash{$array[0]} }, $2);
			}
			if ($array[0] eq "Cached") { last; }
		} close (MEM);
		&debugger_print("Ref", \%memory_hash);
		foreach my $keys (keys %memory_hash) {
			if ($memory_hash{$keys}[1] eq 'kB') {
				$memory_hash{$keys}[0] = $memory_hash{$keys}[0]/1000000;
		}}
		print "+ Total Memory: ".$memory_hash{"MemTotal"}[0]." GiB\n";
		print "+ Free Memory: ".$memory_hash{"MemFree"}[0]." GiB\n";
		print "+ Cached Memory: ".$memory_hash{"Cached"}[0]." GiB\n";
		print "+ Buffers Memory: ".$memory_hash{"Buffers"}[0]." GiB\n";
		$total_available = $memory_hash{"MemFree"}[0]+$memory_hash{"Cached"}[0]; ## Cache Memory would be release once we send a command
		my $p = ($total_available/$memory_hash{"MemTotal"}[0])*100;
		print "Total Available Memory (Cached +  Free): ".$total_available."GiB\n";
		print "+ Percentage of Available Memory: $p %\n";
		
		## Optimize CPUs/taxa
		my $noOfProcesses_SPAdes;
		my $amount_taxa = scalar keys %MID_species_hash;
		if ($total_available > 30) { ## We expect at least 30 GiB of RAM free
			unless ($noOfProcesses > 8) { &printError("To make the most of your server and SPAdes please select more CPUs using -p option") and &dieNicely();}
			# Get number of taxa to assembly
			if ($total_available > 500) {
				$noOfProcesses_SPAdes = int($noOfProcesses/$amount_taxa);
			} elsif ($total_available > 200) {
				$noOfProcesses_SPAdes = int($noOfProcesses/4);
			} elsif ($total_available > 100) {
				$noOfProcesses_SPAdes = int($noOfProcesses/3);
			} elsif ($total_available > 50) {
				$noOfProcesses_SPAdes = $noOfProcesses;
			}

			print "\n\n+ Given the characteristics of the server and memory RAM available, DOMINO has decided to split the $amount_taxa taxa 
			into different subprocesses and assign to each one, a total amount of $noOfProcesses_SPAdes CPUs out of $noOfProcesses CPUs\n\n";
				
			&print_Header(" SPAdes Assembly of each taxa ","#");
			## Call spades for the assembly and send threads for each taxa
	
			my $int_taxa = 0;
			my $pm =  new Parallel::ForkManager($noOfProcesses); 			## Sent child process
			$pm->run_on_finish( 
			sub { my ($pid, $exit_code, $ident) = @_; 
				print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; 
			} );
			$pm->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** SPAdes assembly started with PID $pid\n\n"; } );
			&debugger_print("Ref", \%MID_species_hash);
			foreach my $keys (keys %MID_species_hash) {
				my $pid = $pm->start($int_taxa) and next; print "\nSending child process for SPAdes assembling $keys\n\n";
				$int_taxa++;
				my @files = ("", "");
				if ($file_type == 7) { ## illumina_PE
					for (my $i=0; $i < scalar @{ $MID_species_hash{$keys} }; $i++) {
						my $string2split = $MID_species_hash{$keys}[$i];
						my @array = split("----", $string2split);
						if ($array[0] == 1) { $files[0] = $array[1];
						} else { $files[1] = $array[1]; }  
				}}
	
				## Send SPAdes command				
				my $assembly_dir = $keys."_assembly";
				my $spades_path = "python ".$scripts_path."SPAdes-3.8.1-Linux/bin/spades.py -o ".$assembly_dir." ";
				if ($file_type == 7) { ## illumina_PE
					$spades_path .= "-1 $files[0] -2 $files[1] ";				
				} else { $spades_path .= "-s $MID_species_hash{$keys} ";}
				$spades_path .= "-t $noOfProcesses_SPAdes";
				print "Sending command:\n".$spades_path."\n";
				unless ($debugger) {
					$spades_path .= " > /dev/null"; ## discarding SPAdes output		
				}
				my $system_call = system($spades_path);
				if ($system_call != 0) {
					&printError("Something happened when calling SPAdes for assembly reads..."); &dieNicely();
				} print "\n"; &time_stamp();
	
				## Get contig file
				my $contigs_file = $assembly_dir."/contigs.fasta";
				my $new_contigs_file = $dirname."/assembly_id-".$keys.".contigs.fasta";
				File::Copy::move($contigs_file, $new_contigs_file);
	
				## Finish assembly and generate statistics
				&Contig_Stats($new_contigs_file);
				$pm->finish($int_taxa); # pass an exit code to finish
			}
			$pm->wait_all_children; print "\n** All Assembly child processes have finished...\n\n";		
	
			#######################################
			### Cleaning or renaming some files	###
			#######################################
			print "\n\n";
			unless ($avoidDelTMPfiles) {
				&print_Header(" Cleaning all the intermediary files generated ","#"); 
				&clean_assembling_folders(); print "\n\n";
			}
	
			#######################
			### Finish the job	###
			#######################
			&finish_time_stamp(); print "\n Job done succesfully, exiting the script\n"; exit();
		} else { &printError("Only $total_available GiB available of memory RAM, please bear in mind SPAdes would not complete the task...") and &dieNicely();}	
	} else { &printError("There was an error when retrieving Memory RAM information...\n\nAre you sure this is a linux sever?...") and &dieNicely();}
} else {
	
	#########################################################################
	# 	MIRA ASSEMBLY STEP OF THE READS OF EACH TAXA			#
	#########################################################################
	print "\n"; &print_Header("","#"); &print_Header(" MIRA Assembly Step for each taxa ","#"); &print_Header("","#"); print "\n";
	print "\n"; &print_Header(" Generating an assembly for each taxa ", "%"); print "\n";
	&debugger_print("Ref", \%MID_species_hash);
	foreach my $keys (keys %MID_species_hash) {
		my ($MIRA_manifest_file, @fastq_files_this_taxa, $name_of_project);

		if ($file_type == 7) {
			my @files = ("", "");
			for (my $i=0; $i < scalar @{ $MID_species_hash{$keys} }; $i++) {
				my $string2split = $MID_species_hash{$keys}[$i];
				my @array = split("----", $string2split);
				if ($array[0] == 1) { $files[0] = $array[1];
				} else { $files[1] = $array[1]; }  
			}
			($MIRA_manifest_file, $name_of_project)  = &Generate_manifest_file_pair_end( $files[0], $files[1]);
		} elsif ($file_type == 5) {		
			$MIRA_manifest_file = &Generate_manifest_file($MID_species_hash{$keys}, "Illumina");
			push (@fastq_files_this_taxa, $MID_species_hash{$keys});
		} elsif ($file_type == 3) {
			$MIRA_manifest_file = &Generate_manifest_file($MID_species_hash{$keys}, "454");		
			push (@fastq_files_this_taxa, $MID_species_hash{$keys});
		} 
			
		my $abs_path_manifest_file = $assembly_directory_abs_path."/".$MIRA_manifest_file;
		print $abs_path_manifest_file."\n- Calling MIRA now for $keys assembly...\n- It might take a while...\n";
		my $mira_exe = $MIRA_exec." -t $noOfProcesses ".$abs_path_manifest_file;
		if ($debugger) {
			$mira_exe .= " 2> $error_log"; ## show on screen or maybe print to a file
		} else {
			$mira_exe .= " > /dev/null 2> $error_log"; ## discarding MIRA output		
		}
	
		print $mira_exe."\n\n";
		my $system_call = system($mira_exe);
		if ($system_call != 0) {
			&printError("Something happened when calling MIRA for assembly reads..."); &dieNicely();
		}print "\n"; &time_stamp();
		
		unless ($file_type == 7) {
			my @tmp = split ("\.fastq", $fastq_files_this_taxa[0]);
			$name_of_project = $tmp[0];
		}
		my $folder = $name_of_project."_assembly";
		my $MID_identifier;
		if ($folder =~ /.*id\-(.*)\_assembly/) { $MID_identifier = $1; }
		my $fasta_file = $assembly_directory_abs_path."/".$folder."/".$name_of_project."_d_results/".$name_of_project."_out.unpadded.fasta";
		my $qual_file = $assembly_directory_abs_path."/".$folder."/".$name_of_project."_d_results/".$name_of_project."_out.unpadded.fasta.qual";
		my $contigs_file = $assembly_directory_abs_path."/assembly_id-".$MID_identifier.".contigs-MIRA.fasta";
		my $contigs_file_name = "assembly_id-".$MID_identifier.".contigs-MIRA.fasta";
		my $read_tag_list_file = $assembly_directory_abs_path."/".$folder."/".$name_of_project."_d_info/".$name_of_project."_info_readtaglist.txt";
	
		# Get path for qual and contigs files for CAP3 scaffolding if specified
		$cap3files{$contigs_file_name} = $qual_file;	
	
		## Obtain reads/contigs identified as repeats
		my $array_ref = \@fastq_files_this_taxa;
		&extractReadRepeats($read_tag_list_file, $array_ref, $fasta_file, $contigs_file, $assembly_directory_abs_path."/".$folder);
		chdir $assembly_directory_abs_path;
	}
	print "\n\n";
	&print_Header("","#"); &print_Header(" MIRA Assembly Step finished ","#");  &print_Header("","#"); print "\n\n";
	
	if ($cap3flag) { print "\n\nGetting ready for scaffolding step using CAP3...\n"; }
	## Generates folders and generates link for files of each taxa into them
	
	opendir(DIR, $assembly_directory_abs_path);
	my @files = readdir(DIR);
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		unless (-d $files[$i]) { ## if a directory
			if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
			my $MID_identifier;                     
			if ($files[$i] =~ /(.*)id\-(.*)\.contigs\-MIRA\.fasta/) {
				$MID_identifier = $2; 
				my $fasta_file = $assembly_directory_abs_path."/".$files[$i];
				my $new_fasta_file = $dirname."/".$1."id\-".$2.".contigs.fasta";
					
				if ($cap3flag) { # If cap3 is used for scaffolding move files and generate folders
					my $qual_tmp_file = $cap3files{$files[$i]};
					my $qual_file = $fasta_file.".qual";
					print "Fetching qual values for $fasta_file\n";
					my $hash_identifiers_ref = &read_FASTA_hash_length($fasta_file);
					my %hash_identifiers = %{$hash_identifiers_ref};
					$/ = ">"; ## Telling perl where a new line starts
					open (QUAL, $qual_tmp_file);
					open (OUT, ">$qual_file"); while (<QUAL>) {
					next if /^#/ || /^\s*$/;
					chomp;
					my ($seq_id, $sequence) = split(/\n/,$_,2);
					next unless ($sequence && $seq_id);
					if ($hash_identifiers{$seq_id}) {
						$seq_id = ">".$seq_id;
					print OUT $seq_id."\n".$sequence."\n";                  
					}} $/ = "\n";	close (QUAL); close (OUT);
																	
					## Generate directories
					my $MID_directory = $CAP3_directory_abs_path."/".$MID_identifier;
					unless (-d $MID_directory) {
						mkdir $MID_directory, 0755;
						push (@CAP3_directories, $MID_directory);
					}                                               
					system("ln -s $fasta_file $MID_directory/");
					system("ln -s $qual_file $MID_directory/");                             
				} else { &change_seq_names($fasta_file, $new_fasta_file, $MID_identifier);
}}}}}




###########################################################################
###	cap3 Scaffolding of the contigs MIRA generated for each taxa	###
###########################################################################
if ($cap3flag) { 

	print "\n\n"; &print_Header("","#"); 
	&print_Header(" CAP3 Assembly Step Started ","#"); &print_Header("","#"); print "\n\n";
	for (my $j = 0; $j < scalar @CAP3_directories; $j++) {
		my ($fasta_file, $qual_file); 
		opendir(MID_DIR, $CAP3_directories[$j]);
		my @files = readdir(MID_DIR);
		for (my $i = 0; $i < scalar @files; $i++) {
			if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
			if ($files[$i] =~ /(.*)fasta\.qual/) {
				$qual_file = $files[$i]; next;
			} else {
				$fasta_file = $CAP3_directories[$j]."/".$files[$i];
		}}
	
		###################################
		### System call for cap3 	###
		###################################
		my $command_CAP3 = $CAP3_exec." ".$fasta_file." -o ".$overlap_CAP3." -p ".$similar_CAP3;
		$command_CAP3 .= " 2> $error_log";
		print $command_CAP3."\n";
		my $system_CAP3 = system($command_CAP3);
		if ($system_CAP3 != 0) { &printError("Some error happened when calling CAP3 for assembly reads..."); &dieNicely(); }
	
		###########################################
		### Get contigs and singlets assembled 	###
		###########################################
		&fetch_qual_singlets($fasta_file);
		my $fasta = &get_contigs_Assembly($fasta_file); 
		File::Copy::move($fasta, $dirname);
		print "\n"; &time_stamp();
}}

###########################################
### Cleaning or renaming some files	###
###########################################
print "\n\n";
unless ($avoidDelTMPfiles) {
	&print_Header(" Cleaning all the intermediary files generated ","#"); 
	&clean_assembling_folders(); print "\n\n";
}

###############################
## Generating some statistics #
###############################
chdir $dirname;
my $array_Ref = &read_dir($dirname);
my @dirname_files = @$array_Ref;
for (my $i = 0; $i < scalar @dirname_files; $i++) {
	if ($dirname_files[$i] eq ".DS_Store" || $dirname_files[$i] eq "." || $dirname_files[$i] eq ".." ) { next; }
	if ($dirname_files[$i] =~ /.*\.fasta/) {
		print "+ Generating some statistics for Assembly file: $dirname_files[$i]\n";
		&Contig_Stats($dirname_files[$i]);
		$total_Contigs_all_sets = "";
		undef %nucleotides;
		print "\n\n";
}}

###########################
### 	Finish the job	###
###########################
&finish_time_stamp(); print "\n Job done succesfully, exiting the script\n"; exit();

##########################
##	SUBROUTINES	##
##########################

sub baseCount {
	my $seq = $_[0];
	my $tAs += $seq =~ s/A/A/gi; my $tTs += $seq =~ s/T/T/gi;
	my $tGs += $seq =~ s/G/G/gi; my $tCs += $seq =~ s/C/C/gi;
	my $Ns += (length $seq) - $tAs - $tTs - $tGs - $tCs;
	$nucleotides{"A"} += $tAs; $nucleotides{"T"} += $tTs;
	$nucleotides{"C"} += $tCs; $nucleotides{"G"} += $tGs;
	$nucleotides{"N"} += $Ns;
}

sub blastn {

	##########################################################################################
	##											##
	##  This function uses BLASTN to check for the putative contaminants			##
	##	Jose Fco. Sanchez Herrero, 20/02/2014 	jfsanchezherrero@ub.edu			##
	##########################################################################################

	my $file = $_[0];
	my $db;
	if ($file =~ /(.*)\.fasta/ ) { $db = $1 ; } 
	my $filter = $BLAST."blastn -query ".$file." -evalue 1e-5 -db ".$db." -out blast_search.txt -outfmt 6";
	my $blastn = system($filter);
	if ($blastn != 0) { die "BLASTN failed...\n"; } 	
}

sub calcN50 {
	my @x = @{$_[0]};
	my $n = $_[1];
	my $total;
	for (my $j=0; $j<@x; $j++){ $total += $x[$j]; }
	my ($count, $n50) = (0,0);
	for (my $j=0; $j<@x; $j++){
        $count += $x[$j];
        if($count >= ($total*$n/100)){
            $n50=$x[$j]; last;
    }}
	return $n50;
}
	
sub calcMedian {
	my @arr = @_;
	my @sArr = sort{$a<=>$b} @arr;
	my $arrLen = @arr;
	my $median;
	if($arrLen % 2 == 0) {
		$median = ($sArr[$arrLen/2-1] + $sArr[$arrLen/2])/2;
	} else {
		$median = $sArr[$arrLen/2];
	}
	return $median;
}

sub change_seq_names {
	
	my $file = $_[0];
	my $name_file = $_[1];
	my $identifier = $_[2];
	
	$/ = ">"; ## Telling perl where a new line starts
	open (FILE, $file) or die "Can not open file $file\n";
	open (OUT, ">", "$name_file");
	my $count;
	while (<FILE>) {
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($seq_id, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $seq_id);
    	$count++;
    	$seq_id = ">Contig_".$count."_".$identifier;
    	
    	my @array_seq = split ("\n", $sequence);
    	my $string_seq;
    	for (my $i = 0; $i < scalar @array_seq; $i++) { $string_seq .= $array_seq[$i]; }
		print OUT $seq_id."\n".uc($string_seq)."\n";    	
	}
	close (FILE); 
	close (OUT);	
	$/ = "\n";
}

sub check_file_format {
    
    my $file = $_[0];
    my $id;
    my $count = 3;
    my $fasta = my $fastq = my $qual = 0;
    my $format = 'unknown';

    open(FILE, $file) or &printError("Could not open file $file");
    while (<FILE>) {
        if($count-- == 0) { last;
        } elsif(!$fasta && /^\>\S+\s*/o) { $fasta = 1; $qual = 1;
        } elsif($fasta == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) { $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) { $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) { $id = $1; $fastq = 1;
        } elsif($fastq == 1 && (/^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o)) { $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/o) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/o);
    }}

    if($fasta == 2) { $format = 'fasta';
    } elsif($qual == 2) { $format = 'qual';
    } elsif($fastq == 3) { $format = 'fastq'; }
    return $format;
}

sub check_file {
	## Populate %MID_species_hash with the ids of each file
	my $file_to_check = $_[0];
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
		my $format_returned = &check_file_format($file_to_check);
		my $pair;
		if ($pair_int) {
			if ($pair_int == 1) { $pair = "Left reads file";	
			} elsif ($pair_int == 2) { $pair = "Right reads file";	
			} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\nPair: $pair_int"); &printFormat_message(); &dieNicely();
		}}
		
		if ($format_returned =~ /fastq/) {
			&print_DOMINO_details("\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n\t\tIt also contains an identifier ($name) in the name for later analysis...OK\n");
			if ($pair) { 
				&print_DOMINO_details("\t\tPaired end file = ".$pair." ...OK\n");	
				my $string2push = $pair_int."----".$file_to_check;
				push (@{ $MID_species_hash{$name} }, $string2push);
			} else { 
				$MID_species_hash{$name} = $file_to_check;
		}} else { &printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\n"); &printFormat_message(); &dieNicely();
	}} else { &printError("Please provide a valid FASTQ file: $file_to_check...It is not readable or writable or it does not exist. "); &printFormat_message(); &dieNicely();
	}
	&print_DOMINO_details("\n");
}

sub check_paired_file {

	my $file = $_[0];
	my ($count, $tmp_id, $pair, $pair_return, @pair);
	open (F1, "$file") or &printError("Could not open file $file");;
	my $nLines;
	while (<F1>) { $nLines++; }
	close(F1);
	my $isEOF = 1;
	if($nLines/4 > 0) { $isEOF = 0; }
	my $lineCount = 0;
	open (F1, "$file");
	while(!$isEOF) {
		my @Read = ();
		for(my $i=0; $i<4; $i++) { $Read[$i] = <F1>; }
		chomp(my $QualLine = $Read[3]);
		chomp(my $SeqLine = $Read[1]);
		chomp(my $id = $Read[0]);
		## TODO: add the different types of PE reads if MIRA allows them
		if ($id =~ /(.*)\/(\d+|\w+)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/2
			$tmp_id = $1; $pair = $2;
			$count++;
			push(@pair, $pair);
		} elsif ($id =~ /(\S*)(\d+)\:N\:\.*/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918 2:N:0:CCGTCC
			$tmp_id = $1; $pair = $2;
			$count++;
			push(@pair, $pair);
		} elsif ($id =~ /(.*)\/(\R|\L)/) {
			## @MG00HS16:404:C3W8UACXX:7:1101:1572:1918/L
			$tmp_id = $1; $pair = $2;
			$count++;
			push(@pair, $pair);
		}
		if ($count == 5) { last; }
	}
	close(F1);
	for (my $i = 1; $i < scalar @pair; $i++) {
		if ($pair[0] == $pair[$i]) {
			$pair_return = $pair[$i];
		} else {
			&printError("Some error ocurred when parsing the paired ends files..."); &dieNicely();
	}}
	return $pair_return;
}

sub clean_assembling_folders {
		
	##########################################################################
	##									##
	##  This function generates a cleaning of all the temporary files 	##
	##	generated during the assembly					##
	## 		       							##
	##########################################################################
	
	chdir $dirname;
	my $files_ref = &read_dir($dirname);
	my @files = @$files_ref;
	print "+ Removing temporary files and folders...\n";
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." || $files[$i] eq "tmp_dir") { next; }
		if (-z $files[$i]) { remove_tree($dirname."/".$files[$i]); next;}
		if ($files[$i] =~ /.*ERROR.*/) { next;
		} elsif ($files[$i] =~ /.*Paramete.*/) { next;
		} elsif ($files[$i] eq "CAP3_assemblies") {
			remove_tree($dirname."/".$files[$i]);
		} elsif ($files[$i] eq "MIRA_assemblies") {
			remove_tree($dirname."/".$files[$i]);
}}}

sub Contig_Stats { 
	
	my $fasta_file = $_[0];

	### This is a modification of the script provided in NGS QC toolkit for
	### the analysis of N50, average lentgh etc.
	my (@len, %contig_length, $all_bases);
	if(!defined($fasta_file)) { print "ERROR: No input files are provided\nDOMINO would not die here but not perform any statistics on for the assembly contigs\n"; return; }
	my @name = split("\.fasta", $fasta_file);
	my $outFile = $name[0]."-statistics.txt";
	open(OUT, ">$outFile");
	my @parts = (500, 1000, 2000, 5000, 10000, 50000); my %parts_array;
	open(FILE, $fasta_file) or &printError("Could not open the $fasta_file ...\n") and &dieNicely();
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
		my ($titleline, $sequence) = split(/\n/,$_,2);
		next unless ($sequence && $titleline);
		my @split_titleline = split(" ",$titleline);
		chomp $sequence;
		my $len = length $sequence; ## Get total bases A+T+C+G
		&baseCount($sequence);
		$all_bases += $len; $total_Contigs_all_sets++;
		for (my $j = 0; $j < scalar @parts; $j++) {
			if ($len <= $parts[$j] ) {
				push (@{ $parts_array{$parts[$j]}}, $len); last;
			} elsif ($len > $parts[5]) { ## Get contigs >50kb (if any)
				push (@{ $parts_array{">50Kb"}}, $len); last;			
		}}
	}
	close(FILE);
	$/ = "\n";
	my $nucl_A = $nucleotides{"A"}; my $nucl_T = $nucleotides{"T"};
	my $nucl_C = $nucleotides{"C"}; my $nucl_G = $nucleotides{"G"};
	my $nucl_N = $nucleotides{"N"};
	
	print "## Assembly Statistics ##\n"; 											print OUT "## Assembly Statistics ##\n";
	print "Assembly Statisitcs for file: $fasta_file\n\n";							print OUT "Assembly Statisitcs for file: $fasta_file\n\n";
	print "## General Statistics ##\n";												print OUT "## General Statistics ##\n";
	printf "%-25s %0.2f %s\n", "As:", $nucl_A/$all_bases*100, "%";					printf OUT "%-25s %0.2f %s\n", "As:", $nucl_A/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "Ts:", $nucl_T/$all_bases*100, "%"; 					printf OUT "%-25s %0.2f %s\n", "Ts:", $nucl_T/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "Cs:", $nucl_C/$all_bases*100, "%";					printf OUT "%-25s %0.2f %s\n", "Cs:", $nucl_C/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "Gs:", $nucl_G/$all_bases*100, "%"; 					printf OUT "%-25s %0.2f %s\n", "Gs:", $nucl_G/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "(A + T)s:", ($nucl_A+$nucl_T)/$all_bases*100, "%"; 	printf OUT "%-25s %0.2f %s\n", "(A + T)s:", ($nucl_A+$nucl_T)/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "(G + C)s:", ($nucl_G+$nucl_C)/$all_bases*100, "%"; 	printf OUT "%-25s %0.2f %s\n", "(G + C)s:", ($nucl_G+$nucl_C)/$all_bases*100, "%";
	printf "%-25s %0.2f %s\n", "Ns:", $nucl_N/$all_bases*100, "%";					printf OUT "%-25s %0.2f %s\n", "Ns:", $nucl_N/$all_bases*100, "%";
	print "----------------------------------------\n";								print OUT "----------------------------------------\n";

	foreach my $types (keys %parts_array) {
		my @array = @{ $parts_array{$types} };
		my $total_pb_this_set;
		my $set = $types;
		my @tmp;
		for (my $i = 0; $i < scalar @array; $i++) {
			$total_pb_this_set += $array[$i];
			push (@tmp, $array[$i]);
		}
		print "\nAssembly Statistics for Sequences: <".$set."pb\n"; 		print OUT "\nAssembly Statistics for Sequences: <".$set."pb\n";
		my @sort_array = sort @tmp; my $array_ref_2 = \@sort_array;
		print "Set: <".$set."\n"; &get_stats($array_ref_2,$all_bases, $total_pb_this_set);
		print "----------------------------------------\n"; 			print OUT "----------------------------------------\n";
	}
	
	sub get_stats {
		my $array_ref = $_[0];
		my $all_bases = $_[1];
		my $total_pb_this_set = $_[2];

		my @array = @$array_ref;
		my $bases = $total_pb_this_set;
		my $percentage_pb_bases_this_set = ($bases/$all_bases)*100;
		my $percentage_pb_bases_this_set_print = sprintf ("%0.3f",$percentage_pb_bases_this_set);
		my $totalContigs = scalar @array;
		my $percentage_contigs_this_set = ($totalContigs/$total_Contigs_all_sets)*100;
		my $percentage_contigs_returned = sprintf ("%0.3f", $percentage_contigs_this_set);
		my $avgReadLen = sprintf ("%0.2f", $bases/$totalContigs);
		my $medianLen = calcMedian(@array);
		my $n25 = calcN50($array_ref, 25); my $n50 = calcN50($array_ref, 50); 
		my $n75 = calcN50($array_ref, 75); my $n90 = calcN50($array_ref, 90);
		my $n95 = calcN50($array_ref, 95);
	
		printf "%-25s %d\n" , "Total sequences:", $totalContigs; 						printf OUT "%-25s %d\n" , "Total sequences:", $totalContigs;
		printf "%-25s %0.2f\n" , "Total sequences (%):", $percentage_contigs_this_set; 	printf OUT "%-25s %0.2f\n" , "Total sequences (%):", $percentage_contigs_this_set;
		printf "%-25s %d\n" , "Total bases:", $bases; 									printf OUT "%-25s %d\n" , "Total bases:", $bases;
		printf "%-25s %0.2f\n" , "Total bases(%):", $percentage_pb_bases_this_set_print; printf OUT "%-25s %0.2f\n" , "Total bases(%):", $percentage_pb_bases_this_set_print;
		printf "%-25s %0.2f\n", "Average sequence length:", $avgReadLen; 				printf OUT "%-25s %0.2f\n", "Average sequence length:", $avgReadLen;
		printf "%-25s %0.2f\n", "Median sequence length:", $medianLen; 					printf OUT "%-25s %0.2f\n", "Median sequence length:", $medianLen;
		printf "%-25s %0.2f\n", "N25: ", $n25;											printf OUT "%-25s %0.2f\n", "N25: ", $n25;	
		printf "%-25s %0.2f\n", "N50: ", $n50;											printf OUT "%-25s %0.2f\n", "N50: ", $n50;	
		printf "%-25s %0.2f\n", "N75: ", $n75;											printf OUT "%-25s %0.2f\n", "N75: ", $n75;	
		printf "%-25s %0.2f\n", "N90: ", $n90;											printf OUT "%-25s %0.2f\n", "N90: ", $n90;	
		printf "%-25s %0.2f\n", "N95: ", $n95;											printf OUT "%-25s %0.2f\n", "N95: ", $n95;	
	}
	close(OUT);
}

sub dieNicely { pod2usage(-exitstatus => 1, -verbose => 0); }

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

sub extractReadRepeats {

	my $read_taglist_file = $_[0];
	my $array_files_ref = $_[1];
	my @array_files = @$array_files_ref;
	my $assembled_contigs = $_[2];
    my $clean_contigs = $_[3];
    my $folder = $_[4];

    chdir $folder;
	my $clean_folder_path = &get_earliest("clean_data", $path_abs_folder);
	unless (-d $clean_folder_path) {mkdir $clean_folder_path, 0755;}
	my $preclean_folder = $clean_folder_path."/preRepeatScanner_DOMINO_clean_data";
	unless (-d $preclean_folder) { mkdir $preclean_folder, 0755; }
	
#file ...readtaglist.txt
##
# conName       cFromPadded     cToPadded       cFromUnpadded   cToUnpadded     type    rName   rFromPadded     rToPadded       rFromUnpadded   rToUnpadded     comment
#
##clean_reads.id-sp3_c1   546     453     546     453     HAF3    Seq652_sp3      52      145     52      145     
##clean_reads.id-sp3_c1   452     452     452     452     HAF2    Seq652_sp3      146     146     146     146     
##clean_reads.id-sp3_c1   451     396     451     396     HAF6    Seq653_sp3      147     202     147     202     
	
	my @array_repetitive_reads;
	my @array_repetitive_contigs;
	open(IN, $read_taglist_file) or die( "Couldn't open $read_taglist_file for reading" );
   	while(<IN>) {
     	# do whatever needs to be done
		chomp;
		my $line=$_;
		next if ($line=~/^#/); next if ($line=~/^\s$/);
		my @fields = split(/\s+/, $line);
		my ($contig, $flag, $read, @array_reads);
		$contig = $fields[0]; $flag = $fields[5]; $read = $fields[6];
		if (($flag eq 'HAF6') || ($flag eq 'HAF7') || ($flag eq 'MNRr')){ ## Discard reads identified as highly repetitive
			push (@array_repetitive_reads, $read);
			push (@array_repetitive_contigs, $contig);
	}}
    close(IN);

	my $accns_file_reads = "ids2remove_reads.tab";
	my $accns_file_contigs =  "ids2remove_contigs.tab";
	open(ACCNS_R, ">$accns_file_reads");
	open(ACCNS_C, ">$accns_file_contigs");
	
	my @array_repetitive_contigs_sort = sort(@array_repetitive_contigs);
	my @array_repetitive_contigs_sort_uniq = uniq(@array_repetitive_contigs_sort);
	for (my $i=0; $i < scalar @array_repetitive_contigs_sort_uniq; $i++) { 
		print ACCNS_C $array_repetitive_contigs_sort_uniq[$i]."\n"; 
	} 
	close(ACCNS_C);
	
	my @array_repetitive_reads_sort = sort(@array_repetitive_reads);
	my @array_repetitive_reads_sort_uniq = uniq(@array_repetitive_reads_sort);
	for (my $j=0; $j < scalar @array_repetitive_reads_sort_uniq; $j++) { 
		print ACCNS_R $array_repetitive_reads_sort_uniq[$j]."\n"; 
	}
	close(ACCNS_R);

	## Discard reads and contigs using mothur
	my $abs_path = abs_path();
	for (my $j=0; $j < scalar @array_files; $j++) {
		my $file = $array_files[$j];
		my $file_abs_path;
		for (my $h=0; $h < scalar @file_abs_path; $h++) {
			if ($file_abs_path[$h] =~ /.*$file/) {
				$file_abs_path = $file_abs_path[$h];
				last;	
		}}
		if (-z $accns_file_reads) {
			print "\n+ No reads discarded as there were no repeats identified during MIRA assembly...\n";		
		} else {
			my $line = $mothur_path." '#set.dir(output=$abs_path); remove.seqs(accnos=$accns_file_reads, fastq=$file_abs_path)'";
			print "\n+ Calling mothur executable for discarding reads in file $file...\n\n";
			my $system_call = system($line);
			my @array = split(".fastq", $file); my $filtered_file = $array[0].".pick.fastq";
			File::Copy::move($file_abs_path, $preclean_folder); File::Copy::move($filtered_file, $file_abs_path);
	}}

	my $contigs2cluster;
	if (-z $accns_file_contigs) {
		print "+ No contigs discarded as there were no repeats identified during MIRA assembly...\n";
		$contigs2cluster = $assembled_contigs;		
	} else {
		my $line = $mothur_path." '#set.dir(output=$abs_path); remove.seqs(accnos=$accns_file_contigs, fasta=$assembled_contigs)'";
		print "\n+ Calling mothur executable for discarding contigs in file $assembled_contigs...\n";
		my $system_call = system($line);
		my @array_name = split("/", $assembled_contigs);
		my @array = split(".fasta", $array_name[-1]);
		my $filtered_file = $array[0].".pick.fasta";
		$contigs2cluster = $filtered_file;
	}	
	
	## Use BLAST for clustering sequences
	print "\n\n+ Generate a BLAST database for $contigs2cluster...\n"; &makeblastdb($contigs2cluster);
	print "+ BLAST search now...\n"; &blastn($contigs2cluster);
	
	my $contig_length_Ref = &read_FASTA_hash($contigs2cluster);
	my $perc_aln_desired = 0.85; my $iden_desired = 85;
	
	## Filter BLAST results
	print "+ Filtering BLAST search now...\n";
	my $blast_search = "blast_search.txt";
	my (%contigs_keep, @contigs_seen);
	my $first_hit = 0;
	open (BLAST, $blast_search); while (<BLAST>) {
		my $line = $_;
		chomp $line;
		my ($query, $subject, $perc_iden, $aln, $mismatch, $gap_open, $query_start, $query_end, $subject_start, $subject_end, $eval, $score ) = split("\t", $line);
		my $length_sub = length($$contig_length_Ref{$subject});
		my $perc_aln = $length_sub*$perc_aln_desired;	
		if ($query eq $subject) { next;  }
		if ($eval < 1e-50 && $aln >= $perc_aln && $perc_iden >= $iden_desired) { 
			if ($first_hit == 0) {
				$first_hit++;
				push(@{$contigs_keep{$query}}, $subject); push (@contigs_seen, $subject);
			} else {
				my $flag = 0;
				foreach my $keys (keys %contigs_keep) {
					if (grep /.*$subject.*/, @{$contigs_keep{$keys}}) { $flag = 1; last; }
					if (grep /.*$query.*/, @{$contigs_keep{$keys}}) { $flag = 1; last; }
				}
				if ($flag == 0) { 
					push(@{$contigs_keep{$query}}, $subject); push (@contigs_seen, $subject);	
	}}}}
	close (BLAST);
	
	foreach my $keys (keys %$contig_length_Ref) { unless (grep /$keys/, @contigs_seen) { push (@{ $contigs_keep{$keys}}, $keys); } }
	open (CLN, ">$clean_contigs");
	foreach my $keys (sort keys %contigs_keep) {
		my $largest = length($$contig_length_Ref{$keys});
		my $largest_id = $keys;	
		for (my $i=0; $i < scalar @{$contigs_keep{$keys}}; $i++) {
			my $seq = $contigs_keep{$keys}[$i];
			my $len = length( $$contig_length_Ref{$seq} );
			if ($len > $largest) { $largest = $len; $largest_id = $seq; }
		} 
		print CLN ">$largest_id\n$$contig_length_Ref{$largest_id}\n"; 
	}
	close(CLN);	
	print "+ Done...\n\n";
}

sub fetch_qual_singlets {
	
	my $fasta = $_[0];
	my $singlets_file = $fasta.".cap.singlets";
	my $qual_file = $fasta.".qual";
	my $singlets_qual_file = $singlets_file.".qual";
	
	print "Fetching qual singlets for $fasta\n";
	print $singlets_file."\n";	
	
	my %hash_identifiers;
	$/ = ">"; ## Telling perl where a new line starts
	open (QUAL, $qual_file);
	open (OUT, ">$singlets_qual_file");
	while (<QUAL>) {
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($seq_id, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $seq_id);
    	$seq_id = ">".$seq_id;		
		$hash_identifiers{$seq_id} = $sequence;
	} #while
	close (QUAL);
	
	open (FILE, $singlets_file) or &printError("Could not open file.") and &dieNicely;
	# Get the idenfitiers of the singlet file, all of them
	while (<FILE>) {
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($seq_id, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $seq_id);
		$seq_id = ">".$seq_id;		
		if ($hash_identifiers{$seq_id}) {
			print OUT $seq_id."\n".$hash_identifiers{$seq_id}."\n";
		}
	}
	close (FILE); close (OUT);
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

sub get_fastq_files {

	##########################################################################################
	##	 																					##
	##  This function gets the FASTQ files generated in the previous cleaning step			##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################

	my $clean_folder = &get_earliest("clean_data", $path_abs_folder);
	if ($clean_folder eq 'NO') {
		&printError("No Clean Data folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); &dieNicely();
	}
	opendir(DIR, $clean_folder) or &printError("$clean_folder is not readdable or not accesible...\n") and &dieNicely();
	my @files = readdir(DIR);
	my $flag=0;
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		if ($files[$i] eq "preRepeatScanner_DOMINO_clean_data") { $flag = 1;}
	}
	if ($flag == 1) {
		##
		print "+ A pre assembly folder correcting the repeats have been, using the original clean files...\n";
		my $pre_assembly_dir = $clean_folder."/preRepeatScanner_DOMINO_clean_data";
		opendir(PREDIR, $pre_assembly_dir) or &printError("$pre_assembly_dir is not readdable or not accesible...\n") and &dieNicely();
		my @pre_files = readdir(PREDIR);
		for (my $h=0; $h < scalar @pre_files; $h++) {
			if ($pre_files[$h] eq ".DS_Store" || $pre_files[$h] eq "." || $pre_files[$h] eq ".." ) { next; }
			File::Copy::move($pre_assembly_dir."/".$pre_files[$h], $clean_folder);
		}
		remove_tree($pre_assembly_dir);
	}
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		if ($files[$i] =~ /.*id-.*fastq/) {
			if ($files[$i] =~ /.*MIDtag\.fastq/) { next; }
			push (@fastq_files, $files[$i]);  ## Push only the names			
			push (@file_abs_path, $clean_folder."/".$files[$i]); ## push the whole file path			
			my $file_path = $clean_folder."/".$files[$i];	
			system("ln -s $file_path");		
	}}
}

sub Generate_manifest_file {
	
	##########################################################################################
	##  This function generates a manifest file for Illumina or 454 for MIRA assembler		##
	##	Jose Fco. Sanchez Herrero, 10/06/2014 jfsanchezherrero@ub.edu						##
	##########################################################################################

	my $individual_fastq_file = $_[0];
	my $technology = $_[1];
	my $name_of_project;
	if ($individual_fastq_file =~ /.*(id\-.*)\.f*/) { $name_of_project = $1; }
	&print_Header(" Generate MIRA manifest file ","%"); 
	print "- Generating the manifest file now for $individual_fastq_file...\n";
	
	my $manifest_file = "manifest_$name_of_project.txt";
	open (OUT, ">$manifest_file");
	print OUT "# Manifest file for $technology Project $name_of_project to use MIRA\n";
	print OUT "# First part: defining some basic things\n";
	print OUT "project = $name_of_project\n";
	print OUT "job = genome, denovo, accurate\n";
	my $parameter = "parameters = --hirep_something -NW:cnfs=warn\\\n";  
	$parameter .= "\tCOMMON_SETTINGS -GE:not=$noOfProcesses:amm=on:kpmf=20 -OUT:ors=no:orc=no:otc=no:orw=no:rtd=yes -CL:ascdc=no\\\n";
	if ($technology eq "454") { $parameter .= "\t454_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	} else { $parameter .= "\tSOLEXA_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";}
	print OUT $parameter."\nreadgroup = ".$individual_fastq_file."_reads\ndata = $individual_fastq_file\n";
	if ($technology eq "454") { print OUT "technology = 454\n";
	} else { print OUT "technology = solexa\n"; }
	close (OUT);
	return $manifest_file;	
}

sub Generate_manifest_file_pair_end {
	
	##########################################################################################
	##  This function generates a manifest file for Illumina or 454 for MIRA assembler	##
	##	Jose Fco. Sanchez Herrero, 10/06/2014 jfsanchezherrero@ub.edu			##
	##########################################################################################

	my $left_pair_fastq_file = $_[0];
	my $right_pair_fastq_file = $_[1];
	my $tmp;
	if ($left_pair_fastq_file =~ /.*(id\-.*)\_R\d+\.f*/) { $tmp = $1; }
	my $name_of_project = $tmp;
	&print_Header(" Generate MIRA manifest file ","%"); 
	print "- Generating the manifest file now for $left_pair_fastq_file and $right_pair_fastq_file...\n";
	my $manifest_file = "manifest_$name_of_project.txt";
	open (OUT, ">$manifest_file");
	print OUT "# Manifest file for Solexa Paired End Project $name_of_project to use MIRA\n";
	print OUT "# First part: defining some basic things\n";
	print OUT "project = $name_of_project\n";
	print OUT "job = genome, denovo, accurate\n";
	my $parameter = "parameters = --hirep_something -NW:cnfs=warn\\\n";  
	$parameter .= "\tCOMMON_SETTINGS -GE:not=$noOfProcesses:amm=on:kpmf=20 -OUT:ors=no:orc=no:otc=no:orw=no:rtd=yes -CL:ascdc=no\\\n";
	$parameter .= "\tSOLEXA_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	print OUT $parameter."\n";
	print OUT "readgroup = ".$name_of_project."_reads\nautopairing\n";
	print OUT "data = $left_pair_fastq_file $right_pair_fastq_file\n";
	print OUT "technology = solexa\n";
	close (OUT);
	return ($manifest_file, $name_of_project);	
}

sub get_earliest {
	
	my $option = $_[0];
	my $folder = $_[1];
	#$option == "clean_data, assembly, mapping" 
	
	my $array_files_ref=&read_dir($folder);
	my @array_files = @$array_files_ref;
	my (%mapping_dirs, $earliest);
	for (my $i=0; $i<scalar @array_files;$i++) {
		if ($array_files[$i] eq "." || $array_files[$i] eq ".." || $array_files[$i] eq ".DS_Store") {next;}
		if ($array_files[$i] =~ /(\d+)\_DM\_$option.*/) {
			my $time_stamp=$1;
			if (!$earliest) { $earliest=$time_stamp; } else { if ($time_stamp > $earliest) {$earliest=$time_stamp}}
			$mapping_dirs{$time_stamp} = $folder."/".$array_files[$i];
	}}
	if (!exists $mapping_dirs{$earliest}) { return 'NO';
	} else {  return $mapping_dirs{$earliest}; }	
}

sub get_contigs_Assembly {

	my $fasta = $_[0];
	my $singlets_file = $fasta.".cap.singlets";
	my $contigs_file = $fasta.".cap.contigs";
	my $qual_contigs_file = $fasta.".cap.contigs.qual";
	my $singlets_qual_file = $singlets_file.".qual";
	
	my @tmp = split ("\.contigs\-MIRA\.fasta.*", $fasta);
	my $contigs_singlets_file_fasta = $tmp[0].".contigs.fasta";
	my $contigs_singlets_file_qual = $tmp[0].".contigs.qual";
	
	my $tmp_fasta = "tmp.fasta";
	my $tmp_qual = "tmp.qual";
	open (OUT_fasta, ">$tmp_fasta"); open (OUT_qual, ">$tmp_qual");
	open (CONTIGS, $contigs_file); while (<CONTIGS>) { print OUT_fasta $_; } close (CONTIGS);
	open (SINGLETS, $singlets_file); while (<SINGLETS>) { print OUT_fasta $_; } close (SINGLETS);close (OUT_fasta);
	open (CONTIGS_qual, $qual_contigs_file); while (<CONTIGS_qual>) { print OUT_qual $_; } close (CONTIGS_qual);
	open (SINGLETS_qual, $singlets_qual_file); while (<SINGLETS_qual>) { print OUT_qual $_; } close (SINGLETS_qual); close (OUT_qual);
	
	my $taxa_id;
	if ($contigs_singlets_file_fasta =~ /.*id\-(.*)\.contigs.fasta/) { $taxa_id = $1; }
	&change_seq_names($tmp_fasta, $contigs_singlets_file_fasta, $taxa_id);
	return ($contigs_singlets_file_fasta);
}

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
	my $make_blast_db = $BLAST."makeblastdb";
	$make_blast_db .= " -in ".$file;
	$make_blast_db .= " -dbtype nucl";
	if ($file =~ /(.*)\.fasta/ ) { $db = $1 ; } 
	$make_blast_db .= " -out ".$db;
	$make_blast_db .= " -logfile ".$db.".log";
	
	my $makeblastresult = system($make_blast_db);
	if ($makeblastresult != 0) {
		die "Generating the database failed when trying to proccess the file...\n";
	}
}

sub qualfa2fq_modified_bwa {

	##########################################################################################
	##	 																					##
	##  This is a modified version of the perl script provided in the bwa-0.7.8 package		##	
	## 	to convert fasta+qual files into fastq files										##
	## 																						##
	##	Jose Fco. Sanchez Herrero, 06/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my ($fhs, $fhq, $q, $a);
	my $name = $_[2];
	open($fhs, ($_[0])); #fasta_file
	open($fhq, ($_[1])); #qual_file
	my @temp_name = split ("fasta", $_[0]);
	my $fastq_file;
	if ($name) {
		$fastq_file = $temp_name[0]."fastq_".$name;
	} else {
		$fastq_file = $temp_name[0]."fastq";
	}
	open (FASTQ_out, ">$fastq_file");	
	print "+ Generating a fastq file $fastq_file now...\n";
	
	$/ = ">"; <$fhs>; <$fhq>; $/ = "\n";
	while (<$fhs>) {
		my $id = <$fhq>;
  		print FASTQ_out "\@".$id;
  		$/ = ">";
  		$a = <$fhs>; $q = <$fhq>;
  		chomp($q); chomp($a);
		my @array_seq = split ("\n", $a);
		my $seq;
		for (my $i = 0; $i < scalar @array_seq; $i++) { $seq .= $array_seq[$i]; }
		print FASTQ_out $seq."\n+".$id;
		my @array_qual = split (" ", $q);
		my $qual;
		for (my $j = 0; $j < scalar @array_qual; $j++) {
			if ($array_qual[$j] > 80) {
				$array_qual[$j] = $array_qual[$j] - 30;
			}
			$qual .= $array_qual[$j]." ";
		}
  		$qual =~ s/\s*(\d+)\s*/chr($1+33)/eg; # Convert ASCII characters
  		print FASTQ_out $qual."\n";
  		$/ = "\n";
	}
	close($fhs); close($fhq); close(FASTQ_out);
	$/ = "\n"; ## Telling Perl where a new line starts
	return $fastq_file;
}

sub print_DOMINO_details {
	my $string = $_[0];
	open (PARAM, ">>$param_Detail_file");
	print PARAM $string;
	print STDOUT $string;
	close(PARAM);
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

sub printFormat_message {
	print "\n\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\nWhere:\n\txxx: any character or none.Please avoid using dots (.)\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n\n";
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

sub read_dir {
	my $dir = $_[0];
	opendir(DIR, $dir);
	my @dir_files = readdir(DIR);
	my $array_ref = \@dir_files;
	return $array_ref;
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
   		$sequence =~ s/\n//g;
    	$hash{$titleline} = uc($sequence);
	}
	close(FILE);
	$/ = "\n";
	my $hashRef = \%hash;
	return $hashRef;
}

sub read_FASTA_hash_length {

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
    	for (my $i=0; $i < scalar @array_titleline; $i++) {
    		chomp $sequence;
    		$sequence =~ s/\n//g;
    		my $size = length($sequence);
		   	$hash{$array_titleline[$i]} = $size;
		   	last;
		}
	}
	close(FILE);
	$/ = "\n";
	my $hashRef = \%hash;
	return $hashRef;
}

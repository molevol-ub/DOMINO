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
##      [-input_files file] [-TempFiles] [-SPAdes]
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
	require DOMINO;
	require List::Uniq; use List::Uniq qw(uniq);
	require File::Copy; 
	require File::Find; use File::Find qw(find);			
	require File::Path; use File::Path qw(remove_tree);
	require Cwd; use Cwd qw(abs_path);  
	require Parallel::ForkManager;
}
## TODO: Write a debug message if importing any module fails
## TODO:: use debugger_print for printing Debugging messages

##################################
##	Initializing some variables	##
##################################
my ($helpAsked, %domino_files, $avoidDelTMPfiles, $file_type, $cap3flag,
$manual, $abs_folder, $mrs, $debugger, $flagSpades, $noOfProcesses, $version,
@user_files, $DOMINO_files, $overlap_CAP3, $similar_CAP3, $user_files, $further_information,
@file_abs_path, $step_time, %nucleotides, $check_threads,

$total_Contigs_all_sets);

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
	"check_Threads" => \$check_threads,
	
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

DOMINO v1.0.1

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

[-SPAdes] [-p|--processes int_value] [-mrs|--min_relative_score int_value] [-input_files file] [-TempFiles] 

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

There are two strategies for this assembly phase. User can use as default MIRA and a second scaffolding step using CAP3 (optional) or user can choose to use SPAdes, but we will assume, in this version, the user is using a Linux server with enough CPUs and RAM for assembly.

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

Use CAP3 for a second assembly round. Not recommended. [Default: Off].

=item B<>

=item B<-TempFiles> 

Keep all intermediate files.

=item B<-p|--number_cpu [int_value]>

Number of threads/cores to be used. [Default: 2]

=item B<>

=item B<--SPAdes>

Use SPAdes Genome Assembler for a better assembly. It is mandatory a Linux server and at least 30GiB of free Memory RAM.

Available option for Illumina reads (single or paired-end reads). 454 is not supported anymore. 

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

=item B<Illumina: DOMINO Multiple fastq files -- SPAdes >

 perl DM_Assembly_v1.0.0.pl -o test_folder -DOMINO_files -type_file 5 -SPAdes

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

21 - 09 - 2016

=back

=head1 CITATION

=over 2

Bioinformatics first published online August 16, 2016 
doi:10.1093/bioinformatics/btw534 

=back
=cut

###########################
## Checking user options ##
###########################
if (!$abs_folder || !$file_type) { DOMINO::dieNicely();
} elsif ($user_files) {
	if (!$file_type) {
		print "Some parameters are missing\tPlease provide type of file\n"; 
		DOMINO::printInput_type(); DOMINO::dieNicely();
}}

## Using default options if not provided
if (!$mrs) { $mrs = 80; } #default
if (!$overlap_CAP3) { $overlap_CAP3 = 80; } #default
if (!$similar_CAP3) { $similar_CAP3 = 97; } #default
if (!$noOfProcesses) { $noOfProcesses = 2; } #default

## Getting some PATH variables
my $pipeline_path = abs_path($0); ## abs_path($0) is where the perl script is
my $path_abs_folder = abs_path($abs_folder);
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; }
my $MIRA_exec = $scripts_path."mira_v4.0/bin/mira";
my $CAP3_exec = $scripts_path."cap3/bin/cap3";
my $BLAST = $scripts_path."NCBI_BLAST/";
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
		push (@{ $domino_files{"USER_FILES"} }, $tmp);
}}

## Generate directories
if (-d $dirname) { File::Copy::move($dirname, $dirname."_old_".$random_number); }
mkdir $dirname, 0755; chdir $dirname;
if (-d $dirname_tmp) { File::Copy::move($dirname_tmp, $dirname_tmp."_old_".$random_number); }
mkdir $dirname_tmp, 0755;

## Start the Analysis
print "\n"; DOMINO::printHeader("","+");  DOMINO::printHeader(" Analysis Started ","+"); DOMINO::printHeader("","+"); print "\n";
my $start_time = $step_time = time;
DOMINO::printDetails("Starting the process: [ ".(localtime)." ]\n\n", $param_Detail_file);
# Checking user options
DOMINO::printHeader(" Input File and Parameter Preprocessing ","#");
DOMINO::printDetails("\n+ Output Directory: ".$dirname." ....OK\n", $param_Detail_file);

## If using SPADES for assembly instead of MIRA
my ($assembly_directory, $CAP3_directory);
if ($flagSpades) {
	$assembly_directory = $dirname."/spades_assemblies";
} else {
	$assembly_directory = $dirname."/MIRA_assemblies"; 
	if ($cap3flag) {
		$CAP3_directory = $dirname."/CAP3_assemblies"; 
	}
}
mkdir $assembly_directory, 0755; chdir $assembly_directory;

####################
##	Checking files 	
####################
if ($input_options{$file_type}) {
	if ($file_type == 1 || $file_type == 2 || $file_type == 4|| $file_type == 6) {
		&printError("Input file type provided not supported, please read the documentation..."); 
		DOMINO::printInput_type(); DOMINO::dieNicely();	
	}
} else {
	&printError("\n\nERROR: Wrong type of file provided\nPlease provide a valid type of file:\n");
	DOMINO::printInput_type(); DOMINO::dieNicely();
}

if ($user_files) { ## If users provides files for the assembly either than the DOMINO_clean_data
	if (scalar (@user_files) == 0 || !$file_type) {
		&printError("No input files provided"); DOMINO::printInput_type(); DOMINO::dieNicely();		
	} else {
		## Checking files
		DOMINO::printDetails("+ Type of file(s): Option $file_type : $input_options{$file_type} ...OK\n", $param_Detail_file);
		DOMINO::printDetails("+ User files option provided: ....OK\n", $param_Detail_file);
		DOMINO::printDetails("+ Checking file(s) user provided:\n", $param_Detail_file);
		
		my @file_abs_path_user = @{$domino_files{"USER_FILES"}};
		if (scalar (@file_abs_path_user) == 1) {
			if ($file_type == 3) { 
				&printError("Please provide multiple 454 Roche FASTQ file containing each of the taxa reads");
			} elsif ($file_type == 5) { 
				&printError("Please provide multiple Illumina FASTQ files containing containing each of the taxa reads");
			} elsif ($file_type == 7) {  ## Specify user to input left-right for each taxa: -input sp1_left.fastq -input sp1_right.fastq -input sp2_left.fastq -input sp2_right.fastq
				&printError("Please provide multiple Illumina Paired-End FASTQ files for each taxa (type_file 7) (2 files for specie)");
			}
			DOMINO::printInput_type(); DOMINO::dieNicely();	
		} else {
			for (my $i = 0; $i < scalar @file_abs_path_user; $i++) {
				&check_file($file_abs_path[$i]);
				&debugger_print("ln -s $file_abs_path[$i]"); system("ln -s $file_abs_path[$i]");
	}}}
} else {
	## Go to DOMINO_clean_data and get the files
	DOMINO::printDetails("+ Type of file(s): Option $file_type : $input_options{$file_type} ...OK\n+ No user files provided: Default DOMINO cleaning files ...OK\n+ Get the clean FASTQ files generated in the cleaning step ...OK\n", $param_Detail_file);
	my $clean_folder = DOMINO::get_earliest("clean_data", $path_abs_folder);
	if ($clean_folder eq 'NO') { &printError("No Clean Data folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely(); }
	my $files_ref = DOMINO::readDir($clean_folder);
	my $flag=0;
	if (grep /preRepeatScanner_DOMINO_clean_data/, @$files_ref) {
		print "+ A pre assembly folder correcting the repeats have been, using the original clean files...\n";
		my $pre_assembly_dir = $clean_folder."/preRepeatScanner_DOMINO_clean_data";
		my $pre_files_ref =  DOMINO::readDir($pre_assembly_dir);
		for (my $h=0; $h < scalar @$pre_files_ref; $h++) {
			if ($$pre_files_ref[$h] eq ".DS_Store" || $$pre_files_ref[$h] eq "." || $$pre_files_ref[$h] eq ".." ) { next; }
			File::Copy::move($pre_assembly_dir."/".$$pre_files_ref[$h], $clean_folder);
		}
		remove_tree($pre_assembly_dir);
	}
	for (my $i = 0; $i < scalar @$files_ref; $i++) {
		if ($$files_ref[$i] eq ".DS_Store" || $$files_ref[$i] eq "." || $$files_ref[$i] eq ".." ) { next; }
		if ($$files_ref[$i] =~ /.*MIDtag\.fastq/) { next; }
		if ($$files_ref[$i] =~ /.*id\-.*fastq/) {
			my $file_path = $clean_folder."/".$$files_ref[$i];	
			&check_file($file_path); &debugger_print("ln -s $file_path"); system("ln -s $file_path");		
	}}
}
print "\n"; &time_log();
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

## Check if paired end files provided, all pairs have been correctly provided
if ($file_type eq 7) {
	foreach my $keys (keys %domino_files) {
		my @species_files = @{$domino_files{$keys}{'original'}};
		if (scalar @species_files != 2) {
			&printError("Wrong pair of FASTQ files provided. Please check pair of files corresponding to $keys:\n"); DOMINO::dieNicely();
		} else { &debugger_print("OK: $keys\n\n"); }
}}

foreach my $keys (keys %domino_files) {
	my $dir = $assembly_directory."/".$keys."_assembly";
	push (@{ $domino_files{$keys}{'DIR'} }, $dir);
}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

###########################
## Printing user options ##
###########################
## Print information of the input parameters
DOMINO::printDetails("+ Threads: ".$noOfProcesses." ...OK\n", $param_Detail_file);
DOMINO::printDetails("+ Minimum relative Score (mrs) to use in MIRA assembly: $mrs ...OK\n", $param_Detail_file);
if ($cap3flag) {
	DOMINO::printDetails("+ CAP3 would be used in a second round of assembly ...OK\n", $param_Detail_file);
	DOMINO::printDetails("+ Similarity between reads for CAP3 assembly: $similar_CAP3 ...OK\n", $param_Detail_file);
	DOMINO::printDetails("+ Overlapping between reads for CAP3 assembly: $overlap_CAP3 ...OK\n", $param_Detail_file);
} else { 
if ($flagSpades) { DOMINO::printDetails("+ MIRA has been disabled, SPAdes would perform the assembly ...OK\n", $param_Detail_file);	
} else { DOMINO::printDetails("+ CAP3 has been disabled, only MIRA would perform the assembly ...OK\n", $param_Detail_file);
}}
unless ($avoidDelTMPfiles) {
	DOMINO::printDetails("+ Deleting of temporary files would be done ...OK\n", $param_Detail_file);
} else { DOMINO::printDetails("+ Deleting temporary files would be avoid ...OK\n", $param_Detail_file); }

## Print info about where to print info and error
DOMINO::printDetails("\n+ Parameters details would be print into file: $param_Detail_file...\n", $param_Detail_file);
DOMINO::printDetails("+ Errors occurred during the process would be print into file: $error_log...\n\n", $param_Detail_file);

## Assembly
if ($flagSpades) {

	## If Spades flag, we will assume, in this version, it is using a Linux server with enough CPUs and RAM for assembly	
	## We would use "cat /proc/meminfo " in order to check the available memory, if greater than 30GiB
	## we would split in a reasonable amount of processes to proceed in parallel with the assembly 
	
	# Get memory
	&debugger_print("cat /proc/meminfo > memory_server.txt");
	system("cat /proc/meminfo > memory_server.txt");	
	my $memory_server_file = $assembly_directory."/memory_server.txt";
	my $total_available;
	if (-e -r -s $memory_server_file) {
		DOMINO::printHeader(" Memory RAM Usage retrieval ","#");
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
		my $subProcesses;
		my $amount_taxa = scalar keys %domino_files;
		if ($total_available > 30) { ## We expect at least 30 GiB of RAM free
			unless ($noOfProcesses >= 8) { &printError("To make the most of your server and SPAdes please select more CPUs using -p option") and DOMINO::dieNicely();}
			# Get number of taxa to assembly
			if ($total_available > 500) {
				$noOfProcesses_SPAdes = int($noOfProcesses/$amount_taxa);
				$subProcesses = $amount_taxa;
			} elsif ($total_available > 200) {
				$noOfProcesses_SPAdes = int($noOfProcesses/4);
				$subProcesses = 4;
			} elsif ($total_available > 100) {
				$noOfProcesses_SPAdes = int($noOfProcesses/3);
				$subProcesses = 3;
			} else {
				$noOfProcesses_SPAdes = $noOfProcesses;
				$subProcesses = 1;
			} 
			print "\n\n+ Given the characteristics of the server and memory RAM available, DOMINO has decided to split the $amount_taxa taxa\ninto different subprocesses and assign to each one, a total amount of $noOfProcesses_SPAdes CPUs out of $noOfProcesses CPUs\n\n";
			if ($check_threads) { &debugger_print("Exiting as -check_Threads flag provided, only checking availability of server...\n"); 	exit(); }
				
			DOMINO::printHeader(" SPAdes Assembly of each taxa ","#");
			## Call spades for the assembly and send threads for each taxa
	
			my $int_taxa = 0;
			my $pm =  new Parallel::ForkManager($subProcesses); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
			$pm->run_on_finish( 
			sub { my ($pid, $exit_code, $ident) = @_; 
				print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; 
			} );
			$pm->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** SPAdes assembly started with PID $pid\n\n"; } );

			foreach my $keys (sort keys %domino_files) {
				next if ($keys eq 'original'); 
				my %files_threads;
				my $pid = $pm->start($int_taxa) and next; print "\nSending child process for SPAdes assembling $keys\n\n";
				$int_taxa++;
		
				## Send SPAdes command				
				my $assembly_dir = $domino_files{$keys}{'DIR'}[0];
				my $spades_path = "python ".$scripts_path."SPAdes-3.8.1-Linux/bin/spades.py -o ".$assembly_dir." ";
				if ($file_type == 7) { ## illumina_PE
					$spades_path .= "-1 ".$domino_files{$keys}{'original'}[0]." -2 ".$domino_files{$keys}{'original'}[1];				
				} else { 
					$spades_path .= "-s ".$domino_files{$keys}{'original'}[0];
				}
				$spades_path .= " -t $noOfProcesses_SPAdes";
				unless ($debugger) {
					$spades_path .= " > /dev/null"; ## discarding SPAdes output		
				}
				&debugger_print("Sending command: $spades_path\n");
				print "Sending SPAdes command...\n";
				my $system_call = system($spades_path);
				if ($system_call != 0) {
					&printError("Something happened when calling SPAdes for assembly reads..."); DOMINO::dieNicely();
				} print "\n"; &time_log();
	
				## Get contig file
				my $contigs_file = $assembly_dir."/contigs.fasta";
				my $new_contigs_file = $dirname."/assembly_id-".$keys.".contigs.fasta";
				$files_threads{$keys}{'FINAL'} = $new_contigs_file;
				File::Copy::move($contigs_file, $new_contigs_file);
	
				## Finish assembly and generate statistics
				my $stats = &Contig_Stats($new_contigs_file); 
				$files_threads{$keys}{'stats'} = $stats;

				&debugger_print("DOMINO threads Assembly files"); &debugger_print("Ref", \%files_threads);
				my $spades_dump_hash = $assembly_dir."/dumper_files_threads.txt";
				DOMINO::printDump(\%files_threads, $spades_dump_hash);
				
				$pm->finish($int_taxa); # pass an exit code to finish
			}
			$pm->wait_all_children; print "\n** All Assembly child processes have finished...\n\n";		
	
			## Check each taxa dump file conainting file info
			&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";
			&debugger_print("Retrieve info from files");
			foreach my $taxa (keys %domino_files) {
				if ($domino_files{$taxa}{'DIR'}) {
					my $dump_file = $domino_files{$taxa}{'DIR'}[0]."/dumper_files_threads.txt";
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
	
			#######################################
			### Cleaning or renaming some files	###
			#######################################
			print "\n\n";
			unless ($avoidDelTMPfiles) {
				DOMINO::printHeader(" Cleaning all the intermediary files generated ","#"); 
				&clean_assembling_folders(); print "\n\n";
			}
	
			#######################
			### Finish the job	###
			#######################
			&finish_time_log(); print "\n Job done succesfully, exiting the script\n"; exit();
		} else { &printError("Only $total_available GiB available of memory RAM, please bear in mind SPAdes would not complete the task...") and DOMINO::dieNicely();}	
	} else { &printError("There was an error when retrieving Memory RAM information...\n\nAre you sure this is a linux sever?...") and DOMINO::dieNicely();}
} else {
	
	#########################################################################
	# 	MIRA ASSEMBLY STEP OF THE READS OF EACH TAXA			#
	#########################################################################
	print "\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" MIRA Assembly Step for each taxa ","#"); DOMINO::printHeader("","#"); print "\n";
	print "\n"; DOMINO::printHeader(" Generating an assembly for each taxa ", "%"); print "\n";
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";
	foreach my $keys (sort keys %domino_files) {
		my $MIRA_manifest_file;
		my @species_files = @{$domino_files{$keys}{'original'}};
		
		if (scalar @species_files == 2) {
			$MIRA_manifest_file = &Generate_manifest_file($keys, 'Illumina', $species_files[0], $species_files[1]);	
		} elsif ($file_type == 5) {		
			$MIRA_manifest_file = &Generate_manifest_file($keys, 'Illumina', $species_files[0]);	
		} elsif ($file_type == 3) {
			$MIRA_manifest_file = &Generate_manifest_file($keys, '454', $species_files[0]);	
		}
		
		my $abs_path_manifest_file = $assembly_directory."/".$MIRA_manifest_file;
		push (@{ $domino_files{$keys}{'manifest'}}, $abs_path_manifest_file);
		&debugger_print("MIRA manifest: $abs_path_manifest_file\n");
		print "- Calling MIRA now for $keys assembly...\n- It might take a while...\n";
		my $mira_exe = $MIRA_exec." -t $noOfProcesses ".$abs_path_manifest_file;
		if ($debugger) { $mira_exe .= " 2> $error_log"; ## show on screen or maybe print to a file
		} else { $mira_exe .= " > /dev/null 2> $error_log"; ## discarding MIRA output
		}
	
		print "- Sending MIRA command for $keys...\n"; &debugger_print("MIRA command: $mira_exe\n\n"); 
		my $system_call = system($mira_exe);
		if ($system_call != 0) { &printError("Something happened when calling MIRA for assembly reads..."); DOMINO::dieNicely(); }
		print "\n"; &time_log();
		push(@{ $domino_files{$keys}{'contigs_fasta'}}, $domino_files{$keys}{'DIR'}[0]."/".$keys."_d_results/".$keys."_out.unpadded.fasta");
		push(@{ $domino_files{$keys}{'contigs_qual'}}, $domino_files{$keys}{'DIR'}[0]."/".$keys."_d_results/".$keys."_out.unpadded.fasta.qual");
		push(@{ $domino_files{$keys}{'readTagList'}}, $domino_files{$keys}{'DIR'}[0]."/".$keys."_d_info/".$keys."_info_readtaglist.txt");
		push(@{ $domino_files{$keys}{'contigsMIRA'}}, $domino_files{$keys}{'DIR'}[0]."/assembly_id-".$keys.".contigs-MIRA.fasta");

		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";

		## Obtain reads/contigs identified as repeats
		## Extract read repeats
		
		chdir $domino_files{$keys}{'DIR'}[0];
		my $clean_folder_path = DOMINO::get_earliest("clean_data", $path_abs_folder);
		unless (-d $clean_folder_path) {mkdir $clean_folder_path, 0755;}
		my $preclean_folder = $clean_folder_path."/preRepeatScanner_DOMINO_clean_data";
		unless (-d $preclean_folder) { mkdir $preclean_folder, 0755; }
		
		my (@array_repetitive_reads, @array_repetitive_contigs);
		open(IN, $domino_files{$keys}{'readTagList'}[0]) or &printError( "Couldn't open read tag list for $keys taxa for reading and discarding repeats");
		while(<IN>) {
			# do whatever needs to be done
			chomp;
			my $line=$_;
			next if ($line=~/^#/); next if ($line=~/^\s$/);
	#file ...readtaglist.txt
	##
	# conName       cFromPadded     cToPadded       cFromUnpadded   cToUnpadded     type    rName   rFromPadded     rToPadded       rFromUnpadded   rToUnpadded     comment
	#
	##clean_reads.id-sp3_c1   546     453     546     453     HAF3    Seq652_sp3      52      145     52      145     
	##clean_reads.id-sp3_c1   452     452     452     452     HAF2    Seq652_sp3      146     146     146     146     
	##clean_reads.id-sp3_c1   451     396     451     396     HAF6    Seq653_sp3      147     202     147     202     

			my @fields = split(/\s+/, $line);
			my ($contig, $flag, $read, @array_reads);
			$contig = $fields[0]; $flag = $fields[5]; $read = $fields[6];
			if (($flag eq 'HAF6') || ($flag eq 'HAF7') || ($flag eq 'MNRr')){ ## Discard reads identified as highly repetitive
				push (@array_repetitive_reads, $read);
				push (@array_repetitive_contigs, $contig);
		}}
		close(IN);
	
		## Print contigs to discard
		my $accns_file_contigs =  "ids2remove_contigs_".$keys."tab"; open(ACCNS_C, ">$accns_file_contigs");
		my @array_repetitive_contigs_sort = sort(@array_repetitive_contigs); my @array_repetitive_contigs_sort_uniq = uniq(@array_repetitive_contigs_sort);
		for (my $i=0; $i < scalar @array_repetitive_contigs_sort_uniq; $i++) { print ACCNS_C $array_repetitive_contigs_sort_uniq[$i]."\n"; } close(ACCNS_C);
		
		## Print reads to discard
		my $accns_file_reads = "ids2remove_reads_".$keys.".tab"; open(ACCNS_R, ">$accns_file_reads"); 
		my @array_repetitive_reads_sort = sort(@array_repetitive_reads); my @array_repetitive_reads_sort_uniq = uniq(@array_repetitive_reads_sort);
		for (my $j=0; $j < scalar @array_repetitive_reads_sort_uniq; $j++) { print ACCNS_R $array_repetitive_reads_sort_uniq[$j]."\n"; } close(ACCNS_R);
	
		## Discard reads and contigs using mothur
		if (-z $accns_file_reads) { print "\n+ No reads discarded as there were no repeats identified during MIRA assembly...\n";		
		} else {
			for (my $j=0; $j < scalar @species_files; $j++) {
				print "\n+ Calling mothur executable for discarding reads in file $species_files[$j]...\n\n";
				DOMINO::mothur_remove_seqs($domino_files{$keys}{'original'}[$j], 'YES', $domino_files{$keys}{'DIR'}[0], $accns_file_reads, $mothur_path);
				my $folder_files = DOMINO::readDir($domino_files{$keys}{'DIR'}[0]);
				my @file = grep /.*pick\.fastq/, @{$folder_files};
				File::Copy::move($domino_files{$keys}{'original'}[$j], $preclean_folder);
				File::Copy::move($file[0], $domino_files{$keys}{'original'}[$j]);
		}}
		my $contigs2cluster;
		if (-z $accns_file_contigs) {
			print "+ No contigs discarded as there were no repeats identified during MIRA assembly...\n";
			$contigs2cluster = $domino_files{$keys}{'contigs_fasta'}[0];		
		} else {
			print "\n+ Calling mothur executable for discarding contigs in file ".$domino_files{$keys}{'contigs_fasta'}[0]."...\n";
			DOMINO::mothur_remove_seqs($domino_files{$keys}{'contigs_fasta'}[0], 'NO', $domino_files{$keys}{'DIR'}[0], $accns_file_contigs, $mothur_path);
			my $folder_files = DOMINO::readDir($domino_files{$keys}{'DIR'}[0]);
			my @file = grep /.*pick\.fasta/, @{$folder_files};
			$contigs2cluster = $file[0];
		}	
		
		## Use BLAST for clustering sequences
		print "\n\n+ Generate a BLAST database for $contigs2cluster...\n"; 
		my ($db_generated, $db_generated_message) = DOMINO::makeblastdb($contigs2cluster, $BLAST, $error_log);
		&debugger_print($db_generated_message);
		my $blast_search = "blast_search.txt";
		print "+ BLAST search now...\n"; 
		my ($blastn, $blastn_message) = DOMINO::blastn($contigs2cluster, $db_generated, $blast_search, $BLAST);
		&debugger_print($blastn);
		if ($blastn != 0) { &printError("BLASTN failed...\n"); exit(); } 
		my $contig_length_Ref = DOMINO::readFASTA_hash($contigs2cluster);
		my $perc_aln_desired = 0.85; my $iden_desired = 85;
		
		## Filter BLAST results
		print "+ Filtering BLAST search now...\n";
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
		
		foreach my $keys (sort keys %$contig_length_Ref) { unless (grep /$keys/, @contigs_seen) { push (@{ $contigs_keep{$keys}}, $keys); } }
		my $clean_contigs = $domino_files{$keys}{'contigsMIRA'}[0];
		open (CLN, ">$clean_contigs");
		foreach my $keys (sort keys %contigs_keep) {
			my $largest = length($$contig_length_Ref{$keys});
			my $largest_id = $keys;	
			for (my $i=0; $i < scalar @{$contigs_keep{$keys}}; $i++) {
				my $seq = $contigs_keep{$keys}[$i];
				my $len = length( $$contig_length_Ref{$seq} );
				if ($len > $largest) { $largest = $len; $largest_id = $seq; }
			} 
			print CLN ">$largest_id\n$$contig_length_Ref{$largest_id}"; 
		}
		close(CLN);	
		print "+ Done...\n\n";
		chdir $assembly_directory;
		DOMINO::printHeader(" Assembly finished for $keys ","#"); &time_log(); print "\n\n";
	}
	print "\n\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" MIRA Assembly Step finished ","#");  DOMINO::printHeader("","#"); print "\n\n";
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";
	
	if ($cap3flag) { print "\n\nGetting ready for scaffolding step using CAP3...\n"; mkdir $CAP3_directory, 0755; }
	## Generates folders and generates link for files of each taxa into them
	foreach my $taxa (keys %domino_files) {
		
		my $contigsMIRA_Fasta_file = $domino_files{$taxa}{'contigsMIRA'}[0];
		my $FINAL_fasta_file = $dirname."/assembly_id-".$taxa.".contigs.fasta";			
		push( @{ $domino_files{$taxa}{'FINAL'}}, $FINAL_fasta_file);
					
		if ($cap3flag) { # If cap3 is used for scaffolding move files and generate folders
			&debugger_print("DOMINO files"); &debugger_print("Ref", \%domino_files);
			chdir $CAP3_directory;
			my $qual_tmp_file = $domino_files{$taxa}{'contigs_qual'}[0];
			my $contigsMIRA_qual_file = $domino_files{$taxa}{'contigsMIRA'}[0].".qual";
			push (@{$domino_files{$taxa}{'contigsMIRA_qual'}}, $contigsMIRA_qual_file);
			print "Fetching qual values for $contigsMIRA_Fasta_file\n";
			my $hash_identifiers_ref = DOMINO::readFASTA_IDSfile($contigsMIRA_Fasta_file);
			my %hash_identifiers = %{$hash_identifiers_ref};
			$/ = ">"; ## Telling perl where a new line starts
			open (QUAL, $qual_tmp_file);
			open (OUT, ">$contigsMIRA_qual_file");
			while (<QUAL>) {
				next if /^#/ || /^\s*$/;
				chomp;
				my ($seq_id, $sequence) = split(/\n/,$_,2);
				next unless ($sequence && $seq_id);
				if ($hash_identifiers{$seq_id}) { 
					print OUT ">".$seq_id."\n".$sequence;
					#print ">".$seq_id."\n";
					}
			} $/ = "\n";
			close (QUAL); close (OUT);
																	
			## Generate directories
			my $MID_directory = $CAP3_directory."/".$taxa;
			unless (-d $MID_directory) { mkdir $MID_directory, 0755; }                                               
			&debugger_print("ln -s $contigsMIRA_Fasta_file $MID_directory/"); system("ln -s $contigsMIRA_Fasta_file $MID_directory/"); 
			&debugger_print("ln -s $contigsMIRA_qual_file $MID_directory/"); system("ln -s $contigsMIRA_qual_file $MID_directory/"); 
			push(@{$domino_files{$taxa}{'CAP3_dir'}}, $MID_directory);
			chdir $MID_directory;
				
			###########################################################################
			###	cap3 Scaffolding of the contigs MIRA generated for each taxa	###
			###########################################################################
			&debugger_print("DOMINO files"); &debugger_print("Ref", \%domino_files);
			print "\n\n"; DOMINO::printHeader("","#"); 
			DOMINO::printHeader(" CAP3 Assembly Step Started ","#"); DOMINO::printHeader("","#"); print "\n\n";
			my @tmp_array = split("/", $contigsMIRA_Fasta_file);
			my $fasta_name_contigsMIRA = $tmp_array[-1];
			my $command_CAP3 = $CAP3_exec." ".$fasta_name_contigsMIRA." -o ".$overlap_CAP3." -p ".$similar_CAP3." 2> $error_log";
			&debugger_print("CAP3 command: $command_CAP3\n"); my $system_CAP3 = system($command_CAP3);
			if ($system_CAP3 != 0) { &printError("Some error happened when calling CAP3 for assembly reads..."); DOMINO::dieNicely(); }
			
			###########################################
			### Get contigs and singlets assembled 	###
			###########################################
			my $singlets_file_CAP3 = $fasta_name_contigsMIRA.".cap.singlets"; push (@{$domino_files{$taxa}{'singletsCAP3'}}, $singlets_file_CAP3);
			my $contigs_file_CAP3 = $fasta_name_contigsMIRA.".cap.contigs"; push (@{$domino_files{$taxa}{'contigsCAP3'}}, $contigs_file_CAP3);

			my $tmp_fasta = "tmp.fasta"; open (OUT_fasta, ">$tmp_fasta"); 
			open (CONTIGS, $contigs_file_CAP3); while (<CONTIGS>) { print OUT_fasta $_; } close (CONTIGS);
			open (SINGLETS, $singlets_file_CAP3); while (<SINGLETS>) { print OUT_fasta $_; } close (SINGLETS);close (OUT_fasta);

			&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";
			&change_seq_names($tmp_fasta, $FINAL_fasta_file, $taxa);
			File::Copy::move($FINAL_fasta_file, $dirname);
			print "\n"; &time_log();			
		} else { 
			&change_seq_names($contigsMIRA_Fasta_file, $FINAL_fasta_file, $taxa);
}}} print "\n"; &time_log();

&debugger_print("DOMINO files"); &debugger_print("Ref", \%domino_files);
my $dump_hash = $dirname_tmp."/dumper_assembly_files.txt";
DOMINO::printDump(\%domino_files, $dump_hash);

###############################
## Generating some statistics #
###############################
chdir $dirname;
my $array_Ref = DOMINO::readDir($dirname);
my @dirname_files = @$array_Ref;
foreach my $taxa (keys %domino_files) {
	print "+ Generating some statistics for Assembly file: $domino_files{$taxa}{'FINAL'}[0]\n";
	if ($cap3flag) {
		print "+ Statistics for MIRA assembly:\n";
		my $stats_file_MIRA = &Contig_Stats($domino_files{$taxa}{'contigsMIRA'}[0]);
		$domino_files{$taxa}{'MIRA_stats'} = $stats_file_MIRA;
		print "+ Statistics for CAP3 scaffolding:\n";
	}
	my $stats_file = &Contig_Stats($domino_files{$taxa}{'FINAL'}[0]);
	$domino_files{$taxa}{'FINAL_stats'} = $stats_file;	
	print "\n\n";
}

###########################################
### Cleaning or renaming some files	###
###########################################
print "\n\n";
unless ($avoidDelTMPfiles) {
	DOMINO::printHeader(" Cleaning all the intermediary files generated ","#"); 
	&clean_assembling_folders(); print "\n\n";
}

###########################
### 	Finish the job	###
###########################
&finish_time_log(); print "\n Job done succesfully, exiting the script\n"; exit(0);

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

sub check_file {
	## Populate %domino_files with the ids of each file
	my $file_to_check = $_[0];
	if (-e -r -s $file_to_check) { 
		my $name; my $pair_int; my $typePair;
		if ($file_type == 7) {
			if ($file_to_check =~ /.*id\-(.*)\_(R\d+)\.(fastq|fq)$/) { ## [xxx]id-[yyy](_R[*]).fastq
				($pair_int, $typePair) = DOMINO::check_paired_file($file_to_check, "fastq");
				$name = $1;
				push (@{ $domino_files{$name}{'original'}}, $file_to_check); 
			}		
		} else {
			if ($file_to_check =~ /.*id\-(.*)\.(fastq|fq)$/) {			
				$name = $1;
				push (@{ $domino_files{$name}{'original'}}, $file_to_check); 
			}
		}	
		my $format_returned = DOMINO::check_file_format($file_to_check);
		if ($format_returned =~ /fastq/) {
			DOMINO::printDetails("\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n\t\tIt also contains an identifier ($name) in the name for later analysis...OK\n", $param_Detail_file);
			if ($file_type == 7) { 
				DOMINO::printDetails("\t\tPaired end file = ".$pair_int." ...OK\n", $param_Detail_file);	
			}
		} else { 
			&printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\n"); DOMINO::printFormat_message(); DOMINO::dieNicely();
		}
	} else { 
		&printError("Please provide a valid FASTQ file: $file_to_check...It is not readable or writable or it does not exist. "); DOMINO::printFormat_message(); DOMINO::dieNicely();
	}
	DOMINO::printDetails("\n", $param_Detail_file);
}

sub clean_assembling_folders {
		
	##########################################################################
	##									##
	##  This function generates a cleaning of all the temporary files 	##
	##	generated during the assembly					##
	## 		       							##
	##########################################################################
	
	chdir $dirname;
	my $files_ref = DOMINO::readDir($dirname);
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
		}
	}
}

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
	open(FILE, $fasta_file) or &printError("Could not open the $fasta_file ...\n") and DOMINO::dieNicely();
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
	
	return $outFile;
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

sub finish_time_log {

	my $finish_time = time;
	print "\n\n"; DOMINO::printHeader("","+"); 
	DOMINO::printHeader(" ANALYSIS FINISHED ","+"); 
	DOMINO::printHeader("","+"); 
	print DOMINO::time_stamp."\t";
	my $secs = $finish_time - $start_time; 
	my $hours = int($secs/3600); $secs %= 3600; 	
	my $mins = int($secs/60); $secs %= 60; 
	printf ("Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub Generate_manifest_file {
	
	##########################################################################################
	##  This function generates a manifest file for Illumina or 454 for MIRA assembler		##
	##	Jose Fco. Sanchez Herrero, 10/06/2014 jfsanchezherrero@ub.edu						##
	##########################################################################################

	my $name_of_project = $_[0]; my $technology = $_[1];
	my $fastq_file = $_[2];	my $fastq_file2 = $_[3];
	
	DOMINO::printHeader(" Generate MIRA manifest file for $name_of_project ","%"); 
	my $manifest_file = "manifest_$name_of_project.txt";
	print "- Generating the manifest file $manifest_file...\n";
	
	open (OUT, ">$manifest_file");
	print OUT "# Manifest file for $technology Project $name_of_project to use MIRA\n";
	print OUT "# First part: defining some basic things\n";
	print OUT "project = $name_of_project\n";
	print OUT "job = genome, denovo, accurate\n";
	
	my $parameter = "parameters = --hirep_something -NW:cnfs=warn\\\n";  
	$parameter .= "\tCOMMON_SETTINGS -GE:not=$noOfProcesses:amm=on:kpmf=20 -OUT:ors=no:orc=no:otc=no:orw=no:rtd=yes -CL:ascdc=no\\\n";

	if ($technology eq "454") { 
		$parameter .= "\t454_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	} else { 
		$parameter .= "\tSOLEXA_SETTINGS -AS:mrpc=1 -AL:mrs=$mrs -OUT:sssip=yes";
	}

	if ($fastq_file2) { ## PE reads
		print OUT "readgroup = ".$name_of_project."_reads\nautopairing\n";
		print OUT "data = $fastq_file $fastq_file2\n";
	} else {
		print OUT $parameter."\nreadgroup = ".$name_of_project."_reads\ndata = $fastq_file\n";
	}

	if ($technology eq "454") { 
		print OUT "technology = 454\n";
	} else { 
		print OUT "technology = solexa\n"; }
	close (OUT);
	
	return $manifest_file;	
}

sub printError {
    my $msg = $_[0];
	print "\n\n";DOMINO::printHeader(" ERROR ","!!"); print "\n";
    print $msg."\n\nTry \'perl $0 -h|--help or -man\' for more information.\nExit program.\n";
	print "\n\n"; DOMINO::printHeader("","!!"); DOMINO::printHeader("","!!"); 
    DOMINO::printError_log($msg, $error_log);
}

sub time_log {	
	my $current_time = time;
	print DOMINO::time_stamp."\t";
	my $secs = $current_time - $step_time; 
	my $hours = int($secs/3600); $secs %= 3600; 
	my $mins = int($secs/60); $secs %= 60; 
	$step_time = $current_time;
	printf ("Step took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}
#!/usr/bin/perl
###########################################################################################
### DOMINO: Development of molecular markers in non-model organisms using NGS data 	###
###											
### Authors:								
### Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro
### Sánchez-Gracia, and Julio Rozas.					     		
###########################################################################################
##	Usage:
##      perl DM_Assembly_v1.1.pl
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

##################################
##	Initializing some variables	##
##################################
my ($helpAsked, %domino_files, $avoidDelTMPfiles, $file_type, $cap3flag,
$manual, $abs_folder, $mrs, $debugger, $flagSpades, $noOfProcesses, $version,
@user_files, $DOMINO_files, $overlap_CAP3, $similar_CAP3, $user_files, $further_information,
@file_abs_path, $step_time, $check_threads, $helpAsked1, %domino_params);

######################
## Get user options ##
######################
GetOptions(
	"h" => \$helpAsked1,
	"help" => \$helpAsked,
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
pod2usage( -exitstatus => 0, -verbose => 0 ) if ($helpAsked1);
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

DM_Assembly_v1.1.pl  

=back	

=head1 VERSION

=over 2

DOMINO v1.1 ## Revised 07-11-2018

=back	
	
=head1 SYNOPSIS

=over 2

=item B<>

perl DM_Assembly_v1.1.pl  

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

If QC is not desired, please prepare the files using '-only_tag_files' option using the DM_Clean_v1.1.pl script or provide appropriately tagged files.

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

 perl DM_Assembly_v1.1.pl -o test_folder -DOMINO_files -type_file 3 -p 3 -mrs 80

=item B<454: DOMINO Multiple fastq files. Use CAP3 for scaffolding>

 perl DM_Assembly_v1.1.pl -o test_folder -DOMINO_files -type_file 3 -p 3 -mrs 80 
 -useCAP3 -overCAP3 80 -simCAP3 97 -TempFiles

=item B<Illumina: DOMINO Multiple fastq files>

 perl DM_Assembly_v1.1.pl -o test_folder -DOMINO_files -type_file 5

=item B<Illumina: DOMINO Multiple fastq files -- SPAdes >

 perl DM_Assembly_v1.1.pl -o test_folder -DOMINO_files -type_file 5 -SPAdes

=item B<Illumina Paired-end Files: DOMINO Multiple Paired-end FASTQ files>

 perl DM_Assembly_v1.1.pl -o test_folder -DOMINO_files -type_file 7

=item B<454: User Multiple fastq files>

 perl DM_Assembly_v1.1.pl -o test_folder -type_file 3 -user_files 
 -input_file 454_user_file1.fastq -input_file 454_user_file2.fastq 
 -input_file 454_user_file3.fastq

=item B<Illumina: User Multiple fastq files>

 perl DM_Assembly_v1.1.pl -o test_folder -type_file 5 -user_files 
 -input_file illumina_user_file1.fastq -input_file illumina_user_file2.fastq 
 -input_file illumina_user_file3.fastq

=item B<Illumina Paired-end Files: User Multiple Paired-end FASTQ files>

 perl DM_Assembly_v1.1.pl -o test_folder -type_file 7 -user_files 
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

07 - 11 - 2016

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

## Get input options
my %input_options = %{ DOMINO::input_option_hash() }; 

## Using default options if not provided
if (!$mrs) { $mrs = 80; } #default
if (!$overlap_CAP3) { $overlap_CAP3 = 80; } #default
if (!$similar_CAP3) { $similar_CAP3 = 97; } #default
if (!$noOfProcesses) { $noOfProcesses = 2; } #default
push (@{ $domino_params{'assembly'}{'CPU'} }, $noOfProcesses);

## Getting some PATH variables
my $pipeline_path = abs_path($0); ## abs_path($0) is where the perl script is
my $path_abs_folder = abs_path($abs_folder);
my @script_path_array = split ("/", $pipeline_path);
my $scripts_path;
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; }
my $domino_Scripts = $scripts_path."scripts";

## Checking if the Directory already exists because a previous analysis
unless (-d $path_abs_folder) { mkdir $path_abs_folder; } 
my $random_number = int(rand(1000));
my $datestring = strftime "%Y%m%d%H%M", localtime;

my $dirname = $path_abs_folder."/".$datestring."_DM_assembly";
push (@{ $domino_params{'assembly'}{'folder'} }, $dirname);

my $dirname_tmp = $dirname."/".$datestring."_DM_intermediate_assembly";
my $error_log = $dirname."/".$datestring."_DM_Assembly_ERROR.txt"; push (@{ $domino_params{'assembly'}{'errorLog'} }, $error_log);
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
my $localtime = localtime; push (@{ $domino_params{'assembly'}{'start'} }, $localtime);

# Checking user options
DOMINO::printHeader(" Input File and Parameter Preprocessing ","#");
DOMINO::printDetails("\n+ Output Directory: ".$dirname." ....OK\n", $param_Detail_file);

## If using SPADES for assembly instead of MIRA
my ($assembly_directory, $CAP3_directory);
if ($flagSpades) {
	$assembly_directory = $dirname."/spades_assemblies";
	push (@{ $domino_params{'assembly'}{'spades'} }, "YES");
} else {
	$assembly_directory = $dirname."/MIRA_assemblies"; 
	push (@{ $domino_params{'assembly'}{'MIRA'} }, "YES");
	if ($cap3flag) {
		$CAP3_directory = $dirname."/CAP3_assemblies"; 
		push (@{ $domino_params{'assembly'}{'CAP3'} }, "YES");
		push (@{ $domino_params{'assembly'}{'CAP3_directory'} }, $CAP3_directory);
	}
}
mkdir $assembly_directory, 0755; chdir $assembly_directory;
push (@{ $domino_params{'assembly'}{'dir'} }, $assembly_directory);

####################
##	Checking files 	
####################
if ($input_options{$file_type}) {
	if ($file_type == 1 || $file_type == 2 || $file_type == 4|| $file_type == 6) {
		DOMINO::printError("Input file type provided not supported, please read the documentation..."); 
		DOMINO::printInput_type(); DOMINO::dieNicely();	
	}
} else {
	DOMINO::printError("\n\nERROR: Wrong type of file provided\nPlease provide a valid type of file:\n");
	DOMINO::printInput_type(); DOMINO::dieNicely();
}

if ($user_files) { ## If users provides files for the assembly either than the DOMINO_clean_data
	if (scalar (@user_files) == 0 || !$file_type) {
		DOMINO::printError("No input files provided"); DOMINO::printInput_type(); DOMINO::dieNicely();		
	} else {
		## Checking files
		DOMINO::printDetails("+ Type of file(s): Option $file_type : $input_options{$file_type} ...OK\n", $param_Detail_file);
		push (@{ $domino_params{'assembly'}{'input_type'} }, "$file_type"."--".$input_options{$file_type});

		DOMINO::printDetails("+ User files option provided: ....OK\n", $param_Detail_file);
		DOMINO::printDetails("+ Checking file(s) user provided:\n", $param_Detail_file);
		
		my @file_abs_path_user = @{$domino_files{"USER_FILES"}};
		if (scalar (@file_abs_path_user) == 1) {
			if ($file_type == 3) { 
				DOMINO::printError("Please provide multiple 454 Roche FASTQ file containing each of the taxa reads");
			} elsif ($file_type == 5) { 
				DOMINO::printError("Please provide multiple Illumina FASTQ files containing containing each of the taxa reads");
			} elsif ($file_type == 7) {  ## Specify user to input left-right for each taxa: -input sp1_left.fastq -input sp1_right.fastq -input sp2_left.fastq -input sp2_right.fastq
				DOMINO::printError("Please provide multiple Illumina Paired-End FASTQ files for each taxa (type_file 7) (2 files for specie)");
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
	if ($clean_folder eq 'NO') { DOMINO::printError("No Clean Data folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely(); }
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
			DOMINO::printError("Wrong pair of FASTQ files provided. Please check pair of files corresponding to $keys:\n"); DOMINO::dieNicely();
		} else { &debugger_print("OK: $keys\n\n"); }
}}

foreach my $keys (keys %domino_files) {
	my $dir = $assembly_directory."/".$keys."_assembly";
	push (@{ $domino_files{$keys}{'DIR'} }, $dir);
	push (@{ $domino_params{"assembly"}{'taxa'} }, $keys);
}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
push (@{ $domino_params{"assembly"}{'file_type'} }, $file_type);
push (@{ $domino_params{"assembly"}{'file_type_option'} }, $input_options{$file_type});

###########################
## Printing user options ##
###########################
## Print information of the input parameters
DOMINO::printDetails("+ Threads: ".$noOfProcesses." ...OK\n", $param_Detail_file);
DOMINO::printDetails("+ Minimum relative Score (mrs) to use in MIRA assembly: $mrs ...OK\n", $param_Detail_file);
push (@{ $domino_params{'assembly'}{'mrs'}}, $mrs);

if ($cap3flag) {
	DOMINO::printDetails("+ CAP3 would be used in a second round of assembly ...OK\n", $param_Detail_file);
	DOMINO::printDetails("+ Similarity between reads for CAP3 assembly: $similar_CAP3 ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'assembly'}{'simCAP3'} }, $similar_CAP3);

	DOMINO::printDetails("+ Overlapping between reads for CAP3 assembly: $overlap_CAP3 ...OK\n", $param_Detail_file);
	push (@{ $domino_params{'assembly'}{'overCAP3'} }, $overlap_CAP3);

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

#Print parameters
my $dump_file = $dirname."/DOMINO_dump_information.txt"; DOMINO::printDump(\%domino_files, $dump_file);
my $dump_param = $dirname."/DOMINO_dump_param.txt"; DOMINO::printDump(\%domino_params, $dump_param);

# Lets go!
##############
## Assembly ##
##############
if ($flagSpades) {

	## If Spades flag, we will assume, in this version, it is using a Linux server with enough CPUs and RAM for assembly	
	## We would use "cat /proc/meminfo " in order to check the available memory, if greater than 30GiB
	## we would split in a reasonable amount of processes to proceed in parallel with the assembly 

	my $domino_Scripts_Spades = $domino_Scripts."/DM_runSPAdes.pl";
	my $command = "perl $domino_Scripts_Spades $path_abs_folder"."/"." $step_time"; print $command."\n"; 
	system($command);
	## when finished: check spades.success or spaces.failed in the folder assembly_directory

	my $succesul_run = $dirname."/spades.success"; 
	my $failed_run = $dirname."/spades.fail";
	
	if ($succesul_run) {
		print "SPADes finished successfully...\n";
	} elsif ($failed_run) {
		DOMINO::printError("SPADes failed...\n"); DOMINO::dieNicely();
	} else { DOMINO::printError("SPADes failed...\n"); DOMINO::dieNicely(); }
} else {
	#############################################################
	# 	MIRA ASSEMBLY STEP OF THE READS OF EACH TAXA			#
	#############################################################
	print "\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" MIRA Assembly Step for each taxa ","#"); DOMINO::printHeader("","#"); print "\n";
	print "\n"; DOMINO::printHeader(" Generating an assembly for each taxa ", "%"); print "\n";
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n";

	my $domino_Scripts_MIRA = $domino_Scripts."/DM_runMIRA.pl";
	my $command = "perl $domino_Scripts_MIRA $path_abs_folder"."/"." $step_time"; 
	#print $command."\n"; 
	system($command);
	
	my $succesul_run = $dirname."/MIRA.success"; 
	my $failed_run = $dirname."/MIRA.fail";
	
	if (-e -r -s $succesul_run) {
		print "MIRA finished successfully...\n";
	} elsif (-e -r -s $failed_run) {
		DOMINO::printError("MIRA failed...\n"); DOMINO::dieNicely();
	} else {
		DOMINO::printError("MIRA was interrupted...\n"); DOMINO::dieNicely();
	}
} print "\n"; &time_log();

## Check each taxa dump file conainting file info
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); print "\n"; 
&debugger_print("Retrieve info from files");
my @array_files;
foreach my $taxa (keys %domino_files) {
	if ($domino_files{$taxa}{'DIR'}) {
		my $dump_file = $domino_files{$taxa}{'DIR'}[0]."/dumper_files_threads.txt";
		push (@array_files, $dump_file);
	}	
}
my $dump_file_hash = DOMINO::retrieve_info(\@array_files, %domino_files);
&debugger_print("DOMINO Files"); &debugger_print("Ref", $dump_file_hash); print "\n";

###########################################
### Cleaning or renaming some files	###
###########################################
print "\n\n";
unless ($avoidDelTMPfiles) {
	DOMINO::printHeader(" Cleaning all the intermediary files generated ","#"); 
	&clean_assembling_folders(); print "\n\n";
}
if (-r -e -s $dump_param) { remove_tree($dump_param); $dump_param = $dirname."/DOMINO_dump_param.txt"; DOMINO::printDump(\%domino_params, $dump_param); }
if (-r -e -s $dump_file) { remove_tree($dump_file);   $dump_file = $dirname."/DOMINO_dump_information.txt"; DOMINO::printDump($dump_file_hash, $dump_file);}

###########################
### 	Finish the job	###
###########################
DOMINO::finish_time_stamp($start_time); print "\n Job done succesfully, exiting the script\n"; exit(0);

##########################
##	SUBROUTINES	##
##########################

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
			DOMINO::printError("Wrong FASTQ file provided. Please provide a valid FASTQ file: $file_to_check...\nFormat: $format_returned\n"); DOMINO::printFormat_message(); DOMINO::dieNicely();
		}
	} else { 
		DOMINO::printError("Please provide a valid FASTQ file: $file_to_check...It is not readable or writable or it does not exist. "); DOMINO::printFormat_message(); DOMINO::dieNicely();
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

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}
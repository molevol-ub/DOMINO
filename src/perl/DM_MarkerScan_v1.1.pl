#!/usr/bin/perl
#######################################################################################
###	DOMINO: Development of molecular markers in non-model organisms using NGS data  ###
###											###
###	Authors:									###
###	Cristina Frías-López, José F. Sánchez-Herrero, Miquel A. Arnedo, Alejandro 	###
###	Sánchez-Gracia, and Julio Rozas.					     	###
###											###
#######################################################################################
##	Usage:
##      perl DM_MarkerScan_1.1.pl
##
##    ###########################
##    ### General Information ###
##    ###########################
##      [-h] [--help] [-man] [-v|--version] [-MoreInfo]
##
##    #########################
##    ### Mandatory options ###
##    #########################
##      [-option string] [-type_input string] [-o|--outputFolder string]
##      [-taxa_names string] [-VD|--variable_divergence int_value]
##      [-CL|--conserved_length int_value] [-VL|--variable_length range]
##      [-DM|--development_module selection/discovery]
##
##    ######################
##    ### Optional flags ###
##    ######################
##      ## Input options
##
##      [-genome_fasta file] [-user_contig_files file] [-user_cleanRead_files
##      file] [-msa_file file] [-msa_folder directory] [-RADseq_file file]
##
##      ### Mapping
##
##      [-rdgopen int_value] [-rdgexten int_value] [-rfgopen int_value]
##      [-rfgexten int_value] [-mp|--mismatch_penalty int_value]
##      [-bowtie_local] [-map_contig_files] [-max_SoftClipping int_value]
##      [-SLCD|--significance_level_coverage_distribution float_value]
##      [-low_coverage_data]
##
##      ## Markers
##		[-MPA|--missing_perct_allowed float_value] 
##		[-MCT|--minimum_number_taxa_covered int_value] 
##		[-CD|--conserved_differences int_value] [-PV|--polymorphism] 
##		[-VP|--variable_positions range] 
##		[-SI|--sliding_interval int] [-dnaSP] 
##		[-V-SI_inc int_value] [-C-SI_inc int_value]
##		[-subset_offset_user int_value] [-totalContigs2use4markers int_value]
##		
##      ## Others
##
##      [-No_Profile_Generation|NPG] [-TempFiles] [-keep_bam_file]
##
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use POSIX qw(strftime);
use FindBin;
use lib $FindBin::Bin."/lib";
use List::MoreUtils qw(firstidx);
BEGIN {
	require DOMINO;
	require File::Copy;
	require File::Path; use File::Path qw(remove_tree);
	require Cwd; use Cwd qw(abs_path);  
	require Parallel::ForkManager;
	require Spreadsheet::WriteExcel;
}

##################################
##	Initializing some variables	##
##################################
my (
## User options
$folder, $helpAsked, $avoidDelete_tmp_files, $num_proc_user, $window_var_CONS, 
$window_size_CONS_range, $variable_divergence, $bowtie_local, $variable_positions_user_range,
$window_size_VARS_range, $input_type, $manual, $cigar_pct, $rdgopen, $rdgexten, $rfgexten, 
$rfgopen, $MID_taxa_names, $option, $mis_penalty, $msa_fasta_folder, $polymorphism_user,
$level_significance_coverage_distribution, $map_contig_files, $missing_allowed, $keepbam, 
$version, $DOMINO_simulations, $minimum_number_taxa_covered, $avoid_mapping, $further_information,
@user_cleanRead_files, @user_contig_files, $msa_file, $behaviour, $select_markers, $identify_markers,
$debugger, $helpAsked1, $VAR_inc, $CONS_inc, $option_all, $totalContigs2use4markers,$subset_offset_user,

## absolute path
@contigs_fasta_file_abs_path, @clean_fastq_file_abs_path, # toDiscard

## others
%domino_files, %domino_params, %domino_success_steps, $step_time, %discard_contigs, $pyRAD_file, $stacks_file, $radseq_like_data, 
$number_sp, $genome_fasta, $scripts_path, $dnaSP_flag, $SLIDING_user, %mapping_contigs, $file2dump_param
);


######################
## Get user options	##
######################
GetOptions(
	"h" => \$helpAsked1,
	"help" => \$helpAsked,
	"man" => \$manual,
	"v|version" => \$version, 
	"MoreInfo" => \$further_information,
	##########################

	"DM|development_module=s" => \$behaviour,
	"o|outputFolder=s" => \$folder,
	"p|number_cpu:i" => \$num_proc_user, ## default 2
	"option=s" => \$option,
	"type_input=s" => \$input_type,
	##########################

	"genome_fasta=s" => \$genome_fasta,
	"user_contig_files:s" => \@user_contig_files,
	"user_cleanRead_files:s" => \@user_cleanRead_files,
	"RADseq_file=s" => \$msa_file, ## Loci data enters as msa_file 
	"msa_file=s" => \$msa_file, ## Multiple alignments in one phylip file
	"msa_folder=s" => \$msa_fasta_folder, ## Multiple files in phylip or fasta
	##########################

	"rdgopen|read_gap_open_penalty=i" => \$rdgopen, #5
	"rdgexten|read_gap_extension_penalty=i" => \$rdgexten, #3
	"rfgopen|ref_gap_open_penalty=i" => \$rfgopen, #5
	"rfgexten|ref_gap_extension_penalty=i" => \$rfgexten, #3
 	"mp|mismatch_penalty=i" => \$mis_penalty, #4
	"PV|polymorphism" => \$polymorphism_user,
	"MPA|missing_perct_allowed:i" => \$missing_allowed,
 	"SLCD|significance_level_coverage_distribution=s" => \$level_significance_coverage_distribution,  ## -lscd 1e-05
 	"bowtie_local" => \$bowtie_local,
 	"keep_bam_file" => \$keepbam,
 	"dnaSP" => \$dnaSP_flag,
 	"SI|sliding_interval:i" => \$SLIDING_user,
 	"V-SI_inc:i" => \$VAR_inc,
 	"C-SI_inc:i" => \$CONS_inc,
 	"totalContigs2use4markers:i" => \$totalContigs2use4markers,
 	"subset_offset:i" => \$subset_offset_user,
	##########################

	"taxa_names=s" => \$MID_taxa_names,
	"CD|conserved_differences=i" => \$window_var_CONS, #1 ## variations_conserved_region
	"CL|conserved_length=s" => \$window_size_CONS_range, ##size_conserved_region
	"VD|variable_divergence=f" => \$variable_divergence, ## MINIMUN_variation_percentage: Valor 0-1
	"VL|variable_length=s" =>  \$window_size_VARS_range, #400::600 ## size_variable_region
	"VP|variable_positions=s" => \$variable_positions_user_range, ## 1::5
	"MCT|minimum_number_taxa_covered:i" => \$minimum_number_taxa_covered,
	##########################

	"NPG|No_Profile_Generation" => \$avoid_mapping,
	"TempFiles" => \$avoidDelete_tmp_files,
	"low_coverage_data" => \$DOMINO_simulations,
	"map_contig_files" => \$map_contig_files,
	"max_SoftClipping=i" => \$cigar_pct,
	##########################

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

=item B<#################################################################>

=item B<######## DOMINO Marker Discovery/Selection Pipeline ########>

=item B<#################################################################>

=back

=head1 NAME

=over 2

DM_MarkerScan_1.1.pl

=back
		
=head1 VERSION

=over 2

DOMINO v1.1 ## Revised 07-11-2018

=back
	
=head1 SYNOPSIS

=over 2

=item B<>
	
perl DM_MarkerScan_1.1.pl

=item B<###########################>
	
=item B<### General Information ###>

=item B<###########################>

[-h|--help] [-man] [-v|--version] [-MoreInfo]

=item B<#########################>

=item B<### Mandatory options ###>

=item B<#########################>

[-option string] [-type_input string] [-o|--outputFolder string] [-taxa_names string] [-VD|--variable_divergence int_value] [-CL|--conserved_length range] [-VL|--variable_length range] [-DM|--development_analysis selection/discovery]

=item B<######################>

=item B<###	Optional flags ###>

=item B<######################>

## Input options

[-genome_fasta file] [-user_contig_files file] [-user_cleanRead_files file] [-msa_file file] [-msa_folder directory] [-RADseq_file file]

### Mapping

[-rdgopen int_value] [-rdgexten int_value] [-rfgopen int_value] [-rfgexten int_value] [-mp|--mismatch_penalty int_value] [-bowtie_local] [-map_contig_files] [-max_SoftClipping int_value] [-SLCD|--significance_level_coverage_distribution float_value] [-low_coverage_data]

## Markers

[-MPA|--missing_perct_allowed float_value] [-MCT|--minimum_number_taxa_covered int_value] [-CD|--conserved_differences int_value] [-PV|--polymorphism] [-VP|--variable_positions range] [-SI|--sliding_interval int] [-dnaSP] [-subset_offset_user int_value] [-totalContigs2use4markers int_value]

## Others

[-No_Profile_Generation|NPG] [-TempFiles] [-keep_bam_file]

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

In a typical full run from raw NGS reads (using marker discovery and marker selection modules), the program applies an assembly-based approach; the pipeline is therefore optimized to work with genome partitioning methods in which the size-selected or enriched fragments are long enough to permit a successful assembly and have been fully (or nearly fully) sequenced. For MSA inputs, as for example, restriction-site associated DNA (RAD) variation and related methods (e.g. RADseq, ddRAD, GBS) DOMINO can searh the MSA of each loci (in this case RAD loci) previously obtained with commonly used software packages for highly informative markers (marker selection module); current version accepts the PHYLIP and FASTA format as well as the pyRAD (*.loci) and STACKS (*.fa) output file formats. 

After the maker development phase, DOMINO provides i) the list the genomic regions (with their corresponding coordinates in each contig) exhibiting the minimum levels of nucleotide variation specified by the user among the selected taxa (or covering at least a pre-defined minimum number of them), and the marker MSAs for these taxa (DOMINO marker discovery module) or ii) the list of markers exhibiting the same pre-definable features among a set of aligned sequences supplied by the user and previously obtained in other software (DOMINO marker selection module). 

Furthermore, if the taxa panel has been designed in a convenient phylogenetic context and the user asks for highly conserved regions flanking to be included in the developed markers, these markers should be suitable for further surveys of variation in an extended set of phylogenetically related taxa, i.e. focal taxa. 

=item B<>

=item B<#################################################################>

=item B<######## DOMINO Marker Discovery/Selection Pipeline ########>

=item B<#################################################################>

=item B<>

B<Mapping/Alignment phase>

DOMINO uses Bowtie2 to map the pre-processed reads from each taxon to the assembled contigs of the other m-1 taxa from the panel. Thus, in this step, DOMINO builds m(m-1) sequence alignment/map files (SAM/BAM files). Using the contigs from the assembly of taxon#1 as reference sequences, DOMINO will achieve 3 different alignment/mapping runs (reads of taxon#2 to taxon#1 references; reads of taxon#3 to taxon#1 references; reads of taxon#4 to taxon#1 references). The program will proceed in the same way with the assemblies of taxon#2, taxon#3 and taxon#4. In the case of a panel of m = 4 taxa, for example, DOMINO will build 4 x 3 = 12 SAM/BAM files during this step. The reason behind this particular mapping strategy lies in the dissimilar performance of alignment/mapping algorithms depending on the divergence between the reads and the reference sequences. 

C<Mapping errors and contigs with uneven sequencing depth>

Immediately after generating BAM files, DOMINO removes all unmapped contigs and multi-mapping reads. This step is critical to avoid alignment artifacts, which can create false positive markers (i.e., sequence regions with misleading high levels of nucleotide diversity). The contigs with an unusually large number of aligned reads, which can correspond to repetitive regions, are also removed (they are not suitable for designing single copy markers); in particular we discarded all contigs where the probability of observing the sequenced coverage in one or more positions was < 10-5 (from a Poisson with mean c, where c is the average per contig coverage estimated using all contigs in the sample of the same taxon); this cutoff values, nevertheless, could be changed in the command-line version. Later, DOMINO will build one pileup file per each BAM/SAM file using the SAMtools v0.1.19 suite (mpileup option).

C<Sequencing errors and ambiguity codes>

Since sequencing errors might have a great effect on the marker selection, DOMINO incorporates their own functions for detecting and masking putative sequencing errors, which apply a very conservative criterion for variant calling. First, to avoid the calling of spurious nucleotide variants in low sequencing coverage experiments (i.e. erroneously assigned variants fixed between the taxa from the panel), DOMINO mask the information from positions with only one read mapped to the reference. 

Furthermore, sequencing errors may also inflate the number of called polymorphisms under the Polymorphic Variants option in the marker Discovery/Selection phase. To avoid such undesirable effect, DOMINO incorporates a similar conservative criterion to use only highly credible polymorphisms. Under the Polymorphic Variant option, DOMINO will assume that each taxon represents a diploid individual; that is, in a true heterozygous position it expects approximately the same number of reads from each of the two segregating alleles. For positions with eight or more reads mapped, DOMINO discards those polymorphic variants in which the frequency of the minor allele is significantly lower than the expected (P < 0.05 in Binomial distribution with p = 0.5), likely corresponding to a sequencing error. For lower coverage values, DOMINO will use the information of a polymorphic variant only if the allele with the minor frequency is present in two or more reads. This testing procedure, applied independently for each position within each species, will likely discard many true polymorphic sites; this variant calling approach, however, makes DOMINO highly conservative in detecting true markers when including polymorphisms in the analysis (i.e., DOMINO will use only highly confident within-species segregating variants for the marker Discovery/Selection phase).

Ambiguity codes, either introduced by MIRA assembler in contig sequences or present in user-supplied reference sequences or MSA, are considered by DOMINO to decide whether a position is or not variable. For instance, if there is a Y (IUPAC ambiguity code for C/T) in the reference sequence and all aligned reads show an A at this position, the program computes a nucleotide change (because it is probably a fixed position between taxa). Instead, if all reads show a C or a T, DOMINO does not consider this position as variable in marker Discovery/Selection phase.

After applying all the above-mentioned post-mapping filters, DOMINO combines the variation profiles (arrays with the information about the state of each position, conserved or variable between taxa pairs) obtained from each of the m-1 pileup files including the same reference sequence (i.e., the same taxon), into a single multiple taxa variation profile (MTP). Since each of these references will be likely fragmented in i contigs, DOMINO will build i x m MTP per taxon. Each of these MTP will be independently scanned for regions containing candidate markers in the next phase.

If the user provides reference sequences from a single taxon (e.g. a genome draft), plus the reads from the m different taxa, the program builds only one MTP set (one per contig or scaffold in the supplied reference). This option allows mapping short reads from RAD sequencing (e.g., from RADSeq and related methods) to a reference, generating the variation profiles of each RAD region (RAD loci) for the next marker Discovery/Selection phase. 

On the other hand, if the input includes a single or multiple pre-computed MSA instead of NGS data, DOMINO skips the alignment/mapping phase and directly generates the single MTP set (one per aligned region). In this point, the program accepts MSA files in FASTA (multiple FASTA files, one per linked region), PHYLIP (multiple PHYLIP files, one per linked region, or one multi PHYLIP file with the alignment of all regions) and pyRAD LOCI (*.loci files generated by the program pyRAD) and STACKS fasta (batch_X.fa output files generated from the population analyses in the program STACKS) output files.

B<Marker Discovery/Selection phase>

Each MTP generated in the previous step is either scanned for the presence of candidate marker regions using a sliding window approach (DOMINO marker discovery module) or ii) used to select the markers with the desired features among a set of pre-computed MSA loaded in the previous TAB (DOMINO marker discovery module). In the first case, a specific DOMINO function searches for sequence regions of desired length (Variable region Length, VL), showing the minimum level of variation indicated by the user (Variable region Divergence, VD). DOMINO can also restrict that this variable region was flanked (or not) by highly conserved regions (Conserved region Divergence, CD) of a predefined length (Conserved region Length, CL); an information useful to further design PCR primers. Moreover, DOMINO can strictly restrict the search to a particular set of taxa (from the panel), or just specify the minimum number of taxa required to be covered by the marker (by changing the Minimum number of Covering Taxa parameter; MCT < m). As indicated, DOMINO can use or not the information from polymorphic sites. An appropriated combination of selected taxa and MCT and VD parameter values will allow the user select a large set of informative markers suitable for addressing a wide range of evolutionary questions. In the marker selection module, DOMINO allows directly selecting the most informative markers among the loaded by the user in the same way and with the same personalized features described above. For RAD loci, a particular range of variable positions (VP) between the closest taxa instead of the VD parameter must be specified. The latter option allows selecting informative RAD loci while excluding those exhibiting anomalous high levels of variation, which might reflect RAD tag clustering errors. 

Since DOMINO can work with more than one MTP set (m in a full DOMINO marker discovery run), some of the markers found in MTP based on different reference taxa may be redundant (they can cover the same genomic region, although with different coordinates; see Mapping/Alignment phase section), while other can be found only in one particular profile. To avoid reporting redundant information, we have implemented a function to collapse (using BLAST) these maker sequences, only reporting unique markers. In addition to the MSAs of each identified marker for the selected taxa, DOMINO reports the sequence of all contigs containing candidate markers, along with all relevant information, such as the exact positions delimiting the conserved and variable regions and divergence levels among the selected taxa in the variable region.

To maximize the probability of finding informative markers, the final list of candidates can include overlapped regions that fulfill the specified characteristics (not applicable to the DOMINO marker selection module). Operationally, all regions that meet the criteria for being considered a candidate marker (after moving the scanning window five or more base pairs) are listed as different markers in the final output. In this way, users can choose the best marker to be used directly for further analyses, or the most appropriated region of each contig to be PCR amplified and sequenced in additional focal species (i.e., the best marker from each linked block).

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

=item B<############# INPUT #################>

=item B<#####################################>

=item B<>

=item B<-taxa_names [csv_string]>

Names of the taxa to be included in the marker discovery/selection module.

Provide the names in a comma-separated string.

Example: -taxa_names name_sp1,name_sp2,name_sp3

=item B<>

=item B<-type_input [string]>

Type of NGS sequencing reads.

Two possibilities: single_end or pair_end 

=item B<>

=item B<-option [string]>

Available input options: 

B<>

B<-option DOMINO_files>

Use the files generated by DOMINO in previous phases.

Provide using the -o|--outputFolder the name of the main project folder used to perform previous steps.

B<>
	
B<-option user_assembly_contigs>

Provide one FASTA file per taxon using the option -user_contig_files (See below for further details about this option) and specify the reads to be used for mapping: pre-processed reads from DOMINO [Default], user supplied reads (use the option -user_cleanRead_files file) or the contigs themselves (use flag -map_contig_files, [Not recommended])

B<>

B<-option genome>

Provide a single FASTA file with a references sequence (use the option -genome_fasta file [See below for further details about this option]) and specify the reads to be used for mapping: pre-processed reads from DOMINO [Default] or user supplied reads (use the option -user_cleanRead_files file).

B<>

B<-option msa_alignment>

Provide MSA files in one of the two following formats:

1) A folder containing multiple MSA (each file containing the MSA of and independent region) in PHYLIP or FASTA formats. Use the option -msa_folder path.

2) A multi-MSA file with multiple aligned regions, each one in PHYLIP format. Use the option -msa_file filepath.

B<>

B<-option RADseq>

A multiple MSA file in pyRAD (.loci) or STACKS (.fa) output formats provided with the option -RADseq_file [file].

B<>

=item B<>

=item B<-user_contig_files [file]>

Provide one input FASTA file per species named as "[xxid-][yyy].contigs.fasta" or "[yyy].contigs.fasta". Where:

'[xx]' any character (or none). Please avoid using dots (.)
	
'[yyy]' taxon identifier. Make sure this identifier matches the name provided with -taxa_names and -user_cleanRead_files.
	
C<Example: -option user_assembly_contigs -user_contig_files path/id-Dmelanogaster.contigs.fasta -user_contig_files path/id-Dsimulans.contigs.fasta -user_contig_files path/id-Dyakuba.contigs.fasta -taxa_names Dmelanogaster,Dyakuba,Dsimulans>

=item B<-user_cleanRead_files [file]>

B<>

Provide one input file per taxa named as "[xxid-][yyy][_Rn].fastq" or "[yyy][_Rn].fastq". Where:

'[xxid-]' might be present or not. [xx] could be any character (or none). Please avoid using dots (.)
	
'[yyy]' taxon identifier. [Mandatory]

'[_Rn]' If paired-end data, R1 or R2, for the left and the right reads, respectively.

Single End:
	
C<Example: -user_cleanRead_files path/reads.id-Dmelanogaster.fastq -user_cleanRead_files path/reads.id-Dsimulans.fastq -user_cleanRead_files path/reads.id-Dyakuba.fastq>

Paired-End:

C<Example: -user_cleanRead_files path/reads.id-Dmelanogaster_R1.fastq -user_cleanRead_files path/reads.id-Dmelanogaster_R2.fastq -user_cleanRead_files path/reads.id-Dsimulans_R1.fastq -user_cleanRead_files path/reads.id-Dsimulans_R2.fastq -user_cleanRead_files path/reads.id-Dyakuba_R1.fastq -user_cleanRead_files path/reads.id-Dyakuba_R2.fastq>

=item B<-genome_fasta [file]>

Use this option to provide a a single reference (e.g. a genome, scaffolds, contigs) for the mapping phase, named as "[xxid-][yyy].fasta" or "[yyy].fasta". Where:

'[xxid-]' might be present or not. [xx] could be any character (or none). Please avoid using dots (.)
	
'[yyy]' taxon identifier of the reference genome.

C<Example: -genome_fasta Scaffolds-NCBI_id-Dmelanogaster_NCBI_v6.fasta>
 
=item B<-msa_folder [folder path]>

The folder containing the MSA files in PHYLIP or FASTA formats. Each file must include the MSA of a single aligned sequence region. Use in combination with -option msa_alignment

=item B<-msa_file [file]>

A multi-MSA file with multiple aligned regions, each one in PHYLIP format. Use in combination with -option msa_alignment

=item B<-RADseq_file [file]>

pyRAD (*.loci) or STACKS (xxx.fa) output MSA file.

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

=item B<-DM|--development_module [selection|discovery]>

Select the module for developing markers: 

-DM discovery 

Compatible with -option DOMINO_files, -option genome, -option user_assembly_contigs

-DM selection 

Compatible with -option msa_alignment and -option RADseq

=item B<>

=item B<-CL|--conserved_length [int_value]>

Sequence length for the two conserved regions flanking the marker.

Provide a single value or a range using the :: separator. 

If a unique value is provided, it would be used as the minimun for the conserved length. If two values provided, it would be interpreted as minimun and maximun for this conserved region. 

Example: -CL 30 or -CL 30::45

=item B<-CD|--conserved_differences [int_value]>

Maximun number of nucleotide differences across taxa allowed in the conserved region. [Default: 1]
	
=item B<-VL|--variable_length [range]>

Sequence length for the variable (marker) region.

Provide a single value or a range using the :: separator. 

Example: -VL 500 or -VL 300::500

=item B<-VD|--variable_divergence [float_value]>

Minimum obligatory level of variation between two taxa (among the selected taxa using option -taxa_names csv_string) in a region to be considered as a marker.

Example: -VD 0.01

To specify any level of variation (greater than 0) use a negative value. Example: -VD -1

=item B<-VP|--variable_positions [range]>

Range of variable positions between two taxa allowed in each RAD tag. For RADseq experiments only (-option RADseq).

Provide a range using the :: separator Example: -VP 3::6

Provide a unique value for loci contains exactly that amount of variable positions. Example: -VP 3

To specify only a lower bound for the range use the 999 code. Example: -VP 2::999.

=item B<-rdgopen|--read_gap_open_penalty [int_value]>

Bowtie2 read gap openning penalty. Please refer to Bowtie2 documentation for further information. [Bowtie2 Default: 5]
	
=item B<-rdgexten|--read_gap_extension_penalty [int_value]>

Bowtie2 read gap extension penalty parameter. Please refer to Bowtie2 documentation for further information. [Bowtie2 Default: 3]
	
=item B<-rfgopen|--ref_gap_open_penalty  [int_value]>

Bowtie2 reference gap openning penalty. Please refer to Bowtie2 documentation for further information. [Bowtie2 Default: 5]
	
=item B<-rfgexten|--ref_gap_extension_penalty [int_value]>

Bowtie2 reference gap extension penalty parameter. Please refer to Bowtie2 documentation for further information. [Bowtie2 Default: 3]
	
=item B<-mp|--mismatch_penalty [int_value]>
	
Bowtie2 mismatch penalty parameter. Please refer to Bowtie2 documentation for further information. [Bowtie2 Default: 4]
	
=item B<-SLCD|--significance_level_coverage_distribution [float_value]>

Significance level for the Poisson distribution in coverage-based filtering.

Provided it like scientific annotation: -SLCD 1e-07. [Default: 1e-05] 
		
=item B<-max_SoftClipping [int_value]>

Maximum percentage of soft clipping allowed when mapping. [Default: 5]

=item B<-p|--number_cpu [int_value]>

Number of threads/cores to be used. [Default: 2]
	
=item B<>

=item B<#####################################>

=item B<############ OPTIONS ################>

=item B<#####################################>

=item B<>

=item B<-MPA|--missing_perct_allowed [float_value] [Default 0.05]>

Maximun percentage of missing data per marker region.

=item B<-MCT|--minimum_number_taxa_covered [int_value] [Default Off]> 

Minimum number of taxa required in the alignment of a candidate region to be considered as a marker. 

Example: For a taxa panel of 4 species, specifying -taxa_names sp1,sp2,sp3,sp4 -MCT 3, DOMINO search sequence regions where any combination of 3 of the 4 taxa are present in the alignment.

=item B<-PV|--polymorphism [Default Off]>

Use polymorphic variants to estimate nucleotide VD and CD.

=item B<-low_coverage_data [Default Off]>

Use this option if you expect your data to have really low coverage. Sensibility and precision of the markers identified could be decrease.

=item B<-SI|--sliding_interval [int_value] [Default 10]>

When discovering new markers, DOMINO checks each region using a sliding window approach. This option SI is the increment to loop around this sequence. 

If too much markers are retrieved, please increment this value to obtain more spaced regions.

=item B<-C-SI_inc [int_value] [Default: Self-Adjustment]>

When discovering new markers, DOMINO checks each region using a sliding window approach. This options sets the interval to move if given a range for conserved region size.

If not given, the interval would be adjusted according to the range provided.

=item B<-V-SI_inc [int_value] [Default: Self-Adjustment]>

When discovering new markers, DOMINO checks each region using a sliding window approach. This options sets the interval to move if given a range for the variable region size.

If not given, the interval would be adjusted according to the range provided.

=item B<-dnaSP [Default Off]>

Use this option along with RADseq file or MSA [file|folder] to report any alignment with a minimun number of variations independently of the taxa given.

=item B<-keep_bam_file [Default Off]>

Keep the BAM files generated during the run.

=item B<-No_Profile_Generation|NPG [Default Off]>
	
Use this flag to skip the mapping phase and the generation of the variation profile. Uses data from a previous DOMINO run.

=item B<-totalContigs2use4markers [int]>

By default and in order to speed the computation, DOMINO would only use the largest 20000 contigs generated. Specify a diferent number or use -1 for all the contigs generated during the assembly. [Default: 20000]

=item B<-subset_offset [int]

By default and in order to speed the computation, DOMINO would split the set of contigs to use for markers [totalContigs2use4markers] into subsets of this given size. Specify the number according to the RAM available in your computer. [Default: 50]

=item B<-TempFiles [Default Off]>
	
Keep all intermediate files.

=item B<>

=item B<#####################################>

=item B<##### Command Line Examples #########>

=item B<#####################################>

=item B<DOMINO files: single end -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<DOMINO files: single end, No Mapping  -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 500 -CD 1 -NPG -MCT 2 -MPA 25 -DM discovery 

=item B<DOMINO files: paired-end -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option DOMINO_files
 -type_input pair_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 -DM discovery 

=item B<User provides contigs and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/Dmelanogaster.contigs.fasta -user_contig_files path_to_file2/Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/Dyakuba.contigs.fasta -user_cleanRead_files Dmelanogaster.clean.fastq 
 -user_cleanRead_files Dsimulans.clean.fastq -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<User provides contigs and reads (paired-end) -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option user_assembly_contigs -type_input pair_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -user_cleanRead_files reads_id-Dmelanogaster.clean.R1.fastq -user_cleanRead_files reads_id-Dmelanogaster.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dsimulans.clean.R1.fastq -user_cleanRead_files reads_id-Dsimulans.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dyakuba.clean.R1.fastq -user_cleanRead_files reads_id-Dyakuba.clean.R2.fastq 
 -DM discovery 

=item B<User provides contigs but no reads, and specifies to map contigs vs contigs  -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -map_contig_files -DM discovery 

=item B<User provides a reference genome and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option genome -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -genome_fasta path_to_genomes_folder/NCBI_id-Dpseudobscura.fasta
 -user_cleanRead_files Dmelanogaster.clean.fastq -user_cleanRead_files Dsimulans.clean.fastq 
 -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<MSA alignment: A single MSA in PHYLIP -- Single File -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/ -msa_file file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single MSA in FASTA -- Single File -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/ -msa_file file.fasta 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Discovery>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Selection>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in PHYLIP -- Folder -- Selection>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in FASTA MSA -- Folder -- Selection>

 perl DM_MarkerScan_v1.1.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in pyRAD format -- File -- Selection>

 perl DM_MarkerScan_v1.1.pl -option RADseq -o test/  -RADseq_file output.loci
 -taxa_names ind1,ind3,ind5,ind8,ind9 -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in STACKS format -- File -- Selection>

 perl DM_MarkerScan_v1.1.pl -option RADseq -o test/  -RADseq_file output.fa
 -taxa_names ind1,ind3,ind5,ind8,ind9 -VD 0.01 -DM selection

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

23 - 01 - 2018

=back

=head1 CITATION

=over 2

Bioinformatics first published online August 16, 2016 
doi:10.1093/bioinformatics/btw534 

=back

=cut

#############################################################################################################################
## Check mandatory options
if (!$folder) { DOMINO::printError("No folder provided...\n"); DOMINO::dieNicely(); }
my %type_input = ("single_end" => 1, "pair_end" => 1,);
if (!$type_input{$input_type}) {
	DOMINO::printError("Input type provided ($input_type) does not match single_end or pair_end\n");
	DOMINO::dieNicely();}
#############################################################################################################################

#######################################
###	Initialise some PATH variables	### 
#######################################
## Directory names
unless (-d $folder) { mkdir $folder, 0755; }
my $folder_abs_path = abs_path($folder);

## Dates & timestamps
my $random_number = int(rand(100));
my $datestring = strftime "%Y%m%d%H%M", localtime;
my $date = localtime; 
my $start_time = $step_time = time;

## mapping dirname
my $align_dirname = $folder_abs_path."/".$datestring."_DM_mapping"; 
my $mapping_parameters = $folder_abs_path."/".$datestring."_Mapping-Parameters.txt";
my $mapping_markers_errors_details = $folder_abs_path."/".$datestring."_Mapping_ERROR.txt";

## Markers dirname
my $marker_dirname = $folder_abs_path."/".$datestring."_DM_markers";
my $param_Detail_file_markers = $folder_abs_path."/".$datestring."_Markers-Parameters.txt";

## General binaries variables
my $pipeline_path = abs_path($0);
my @script_path_array = split ("/", $pipeline_path);
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; } 
my $domino_Scripts = $scripts_path."scripts";

############################################################################################################
## Check parameters provided, check older runs and print option
&check_options();

## Any previous run?
if ($avoid_mapping) { 
	my $answer = &check_previous(); 
	if ($answer eq "YES") { 
		$avoid_mapping = 1;
	} else { undef $avoid_mapping; }
}

## print options to screen
&print_options();
#######################################

## Debug print Dumper \%domino_files; print Dumper \%domino_params; exit();
my ($dump_file, $dump_param);
$dump_file = $align_dirname."/DOMINO_dump_information.txt"; 
$dump_param = $align_dirname."/DOMINO_dump_param.txt"; 

####################################################################################
## 								START
####################################################################################

if (!$avoid_mapping) {	
	####################################################################################
	##########	Mapping/Alignment of the contigs 		################################
	####################################################################################

	## dump information
	DOMINO::printDump(\%domino_files, $dump_file);
	DOMINO::printDump(\%domino_params, $dump_param);
	
	if ($option ne "msa_alignment") {
	
		#################################################################################################
		####
		####					MAPPING
		####
		## We would use Bowtie2 for mapping the reads in a separate script DM_MappingReads.pl	
		DOMINO::printHeader("", "#");	DOMINO::printHeader(" Mapping Process started ", "#"); DOMINO::printHeader("", "#"); print "\n";

		my $domino_Scripts_Mapping = $domino_Scripts."/DM_MappingReads.pl";
		my $command = "perl $domino_Scripts_Mapping ".$folder_abs_path." $step_time";
		print "\n[ System Call: ".$command." ]\n\n";
		system($command);
		#################################################################################################
		
		## retrieve information generated
		my @array = ($dump_file); %domino_files = %{ DOMINO::retrieve_info(\@array, \%domino_files) };
		
	} elsif ($option eq "msa_alignment") {

		#################################################################################################
		####
		####		Parse alignment
		####
		DOMINO::printHeader("", "#"); DOMINO::printHeader(" Parsing Alignment Process started ", "#"); DOMINO::printHeader("", "#"); print "\n";
		## to Debug
		## DM_ParseMSA_files.pl $folder_abs_path
		#################################################################################################
	}
	## Dump info up to now to a file if (-r -e -s $dump_file) { remove_tree($dump_file) }; DOMINO::printDump(\%domino_files, $dump_file);	
	## Dump parameters to a file if (-r -e -s $dump_param) { remove_tree($dump_param) }; DOMINO::printDump(\%domino_params, $dump_param);	

	## Move parameters and error file to folder
	File::Copy::move($mapping_parameters, $align_dirname."/");
	
} else {
	print "+ Files would be obtained...\n\n"; ## No_Profile_Generation|NPG: just get files
	## retrieve information generated
	my @array = ($dump_file); %domino_files = %{ DOMINO::retrieve_info(\@array, \%domino_files) };
	my @array_param = ($dump_param); %domino_params = %{ DOMINO::retrieve_info(\@array_param, \%domino_params) };
	## %domino_params
}

##########################################################################################
###################### MARKER DEVELOPMENT ################################################
##########################################################################################
## Print info about where to print info and error
DOMINO::printHeader("", "+"); DOMINO::printHeader(" Analysis of Molecular Markers started ", "+"); DOMINO::printHeader("", "+"); print "\n";
DOMINO::printDetails("\n+ Markers Directory: ".$marker_dirname." ...OK\n", $param_Detail_file_markers);
if (-d $marker_dirname) {  
	File::Copy::move($marker_dirname, $marker_dirname."_old_".$random_number);
	DOMINO::printDetails("+ Changing an old folder named as $marker_dirname to $marker_dirname"."_old_"."$random_number...OK\n", $param_Detail_file_markers);     						
} 
mkdir $marker_dirname, 0755; chdir $marker_dirname; &debugger_print("Changing dir to $marker_dirname");
if ($avoid_mapping) { 
	File::Copy::move($mapping_parameters, $marker_dirname."/"); 
}

## dump information
$dump_file = $marker_dirname."/DOMINO_dump_information.txt"; DOMINO::printDump(\%domino_files, $dump_file);
$dump_param = $marker_dirname."/DOMINO_dump_param.txt"; DOMINO::printDump(\%domino_params, $dump_param);

##########################################################
############# MARKER SELECTION ###########################
##########################################################
### MSA alignment
if ($option eq "msa_alignment") {
	#################################################################################################
	DOMINO::printHeader("", "#"); DOMINO::printHeader(" Parsing Alignment Process started ", "#"); DOMINO::printHeader("", "#"); print "\n";
	## to Debug
	## to Debug	## DM_MarkerSelection.pl $folder_abs_path
	## to Debug
	## to Debug	my $domino_Scripts_MarkerSelection = $domino_Scripts."/DM_MarkerSelection.pl";
	## to Debug	my $command = "perl $domino_Scripts_MarkerSelection ".$folder_abs_path." $step_time";
	## to Debug	print "\n[ System Call: ".$command." ]\n\n";
	## to Debug	system($command);
	#################################################################################################
}
&time_log(); print "\n";
## Debug print Dumper \%domino_files; print Dumper \%domino_params; exit();

## Other types of data
##########################################################
############# MARKER DISCOVERY ###########################
##########################################################
my $genome_marker_bool = 0;
my $all_markers_file = $marker_dirname."/markers.txt";
my $test=0; my $marker_success = 0;
foreach my $ref_taxa (sort keys %domino_files) {

	## For each taxa specified, obtain putative molecular markers
	unless ($domino_files{$ref_taxa}{'contigs'}) {next; }
	if ($genome_marker_bool == 1) {last;}
	print "\n";
	
	## Create a dir for each taxa
	DOMINO::printHeader(" Checking taxa files user specified ", "#"); 
	my $marker_dir;
	if ($genome_fasta) {
		print "Checking: \tGenome provided: $domino_files{$ref_taxa}{'contigs'}[0]\n\n";
		$genome_marker_bool = 1; $marker_dir = $marker_dirname."/DOMINO_markers_Genome";
	} else {
		print "Checking: \t$ref_taxa\n\n";
		$marker_dir = $marker_dirname."/markers_Ref_".$ref_taxa;
	}

	## DM_MarkerSelection.pl $folder_abs_path
	my $domino_Scripts_MarkerDiscovery = $domino_Scripts."/DM_MarkerDiscovery.pl";
	my $command = "perl $domino_Scripts_MarkerDiscovery ".$folder_abs_path." $step_time $ref_taxa $marker_dir";
	print "\n[ System Call: ".$command." ]\n\n";
	system($command);
	#################################################################################################

	## Check for marker_discovery.success/failed file
	if (-r -e -s $marker_dir."/marker_discovery.success") {
		$marker_success++;
	}
	$test++; if ($test == 2) { last; }	

} #each reference taxa

#################################################################################################
### LETS FINISH
#################################################################################################

if ($marker_success > 0) {
	## continue as there are some markers identified
} else {
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" No markers identified ", "#"); DOMINO::printHeader("", "#");
	## Finish and exit
	print "+ Termination of DOMINO marker identification\n+ No markers were identified using these parameters...\n";
	DOMINO::finish_time_stamp($start_time); print "\n\n Early termination, exiting the script\n\n\n"; 
	exit();
}

chdir $marker_dirname;
if ($genome_fasta) {
	## Print excel for clusterized results
	print "+ Generating an Excel file for DOMINO markers coordinates...\n"; 
	my $coordinates;
	foreach my $ref_taxa (sort keys %domino_files) { ## For each taxa specified, obtain putative molecular markers
		unless ($domino_files{$ref_taxa}{'contigs'}) {next; }
		$coordinates = $domino_files{$ref_taxa}{'coordinates'}[0];
	}
	my $excelbook = DOMINO::print_Excel(\$coordinates, \$marker_dirname);
	#################################################################################################

} else { ## multiple references

	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Clustering markers for unique results ", "#"); DOMINO::printHeader("", "#");

	### Clusterize markers for each taxa and provide unique results
	my $domino_Scripts_MarkerClusterize = $domino_Scripts."/DM_MarkerClusterize.pl";
	my $command = "perl $domino_Scripts_MarkerClusterize ".$folder_abs_path." $step_time";
	print "\n[ System Call: ".$command." ]\n\n";
	system($command);
	#################################################################################################
}

## Dumping info to a file
my $dump_info_DOMINO = $marker_dirname."/DOMINO_dump_information.txt"; if (-r -e -s $dump_info_DOMINO) { remove_tree($dump_info_DOMINO) }; DOMINO::printDump(\%domino_files, $dump_info_DOMINO);
my $dump_info_success = $marker_dirname."/DOMINO_dump_success.txt"; DOMINO::printDump(\%domino_success_steps, $dump_info_success);
unless ($avoidDelete_tmp_files) { 
	###########################
	## Delete temporary file ##
	###########################
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Deleting temporary files ", "#"); DOMINO::printHeader("", "#");
	## discard files and tmp folders
	foreach my $taxa (sort keys %domino_files) {
		if ($domino_files{$taxa}{'taxa'}) {
			remove_tree($domino_files{$taxa}{'PROFILE_FOLDER'}[0]);
			remove_tree($domino_files{$taxa}{'REF_DIR'}[0]);
			remove_tree($domino_files{$taxa}{'array_all_taxa'}[0]);
	}}
	#remove_tree($blast_dir);
}

## Move parameters and error file to folder
File::Copy::move($param_Detail_file_markers, $marker_dirname."/");
if (-z $mapping_markers_errors_details) { remove_tree($mapping_markers_errors_details); 
} else { File::Copy::move($mapping_markers_errors_details, $marker_dirname."/"); }

## Finish and exit
DOMINO::finish_time_stamp($start_time); print "\n\n Job done succesfully, exiting the script\n\n\n"; 
exit();


######################################################################################
##																					##	
##									SUBROUTINES										##
##																					##
######################################################################################

sub check_file {
	
	my $file_to_check = $_[0];
	my $id = $_[1];
	if (-e -r -s $file_to_check) { 
		if ($file_to_check =~ /(.*)\.fast.* || (.*)\.fa || (.*)\.fq/) { ## File is named fasta/fa
			my $format_returned = DOMINO::check_file_format($file_to_check);
			if ($format_returned eq "fasta") { ## fasta file is ok
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} elsif ($format_returned eq "fastq") { ##fastq file is ok
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} else { DOMINO::printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
		} else { ## File is not named fasta/fa or fastq/fq
			my $format_returned = DOMINO::check_file_format($file_to_check);
			if ($format_returned eq "fasta") { ## It is actually a fasta file
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} elsif ($format_returned eq "fastq") { ##fastq file is ok
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} else { DOMINO::printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
		}
		
		## Get file name
		my @file_to_check = split("/", $file_to_check);
		my $file_name = $file_to_check[-1];
		if ($id) {
			print "\t\tIt also contains an identifier ($id) in the name for later analysis...OK\n";
		} else { 
			if ($genome_fasta) {
				DOMINO::printError("$file_to_check...\n File does not contain any name provided using taxa_names option.\nMaybe it is under controlled but just bear it in mind...\nDOMINO would not die here...\n"); 
		}}
	} else { DOMINO::printError("Please provide several files for each taxa..."); DOMINO::printFormat_message(); DOMINO::dieNicely(); }
}

sub debugger_print {
	my $string = $_[0];
	my $ref = $_[1];
	## Print to a file
	if ($debugger) {
		if ($string eq "Ref") {
			print "\n******\nDEBUG:\n";
			print Dumper $ref; ## if filehandle OUT: print OUT Dumper $ref;
			print "******\n";
		} else {
			print "\n******\nDEBUG:\n".$string."\n******\n\n";
}}}

sub fastq_files {

	my $files_ref = $_[0];
	my @files = @$files_ref;
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		
		my @file_path_name = split("/",$files[$i]);
		my $file_name = $file_path_name[-1];
		
		if ($file_name =~ /.*id-(.*)(\_R\d+)\.fastq/g) {
			if ($domino_files{$1}{'taxa'}) { 
				push (@{ $domino_files{$1}{'reads'} }, $files[$i]); ## push the whole file path			
			} else { &printError("Please check the tag for the file $files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
		} elsif ($file_name =~ /.*id-(.*)\.fastq/g) {
			if ($domino_files{$1}{'taxa'}) { 
				push (@{ $domino_files{$1}{'reads'} }, $files[$i]); ## push the whole file path			
			} else { DOMINO::printError("Please check the tag for the file $files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
		} elsif ($files[$i] =~ /.*id-(.*)\.fastq/g) {
			if ($domino_files{$1}{'taxa'}) { 
				push (@{ $domino_files{$1}{'reads'} }, $files[$i]); ## push the whole file path			
			} else { DOMINO::printError("Please check the tag for the file $files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
		}
	}
}

sub get_clean_files {
	
	my $clean_folder = DOMINO::get_earliest("clean_data", $folder_abs_path);
	if ($clean_folder eq 'clean_data') {
		DOMINO::printError("No clean_data folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely();
	}
	DOMINO::printDetails("+ DOMINO clean reads would be retreived from: $clean_folder\n", $mapping_parameters, $param_Detail_file_markers);
	my $files_dir_ref = DOMINO::readDir($clean_folder);
	my @files = @$files_dir_ref; my @files2;
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		push (@files2, $clean_folder."/".$files[$i]);
	}
	&fastq_files(\@files2);
}

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}

sub user_cleanRead_files {	
	my @array;
	for (my $i=0; $i < scalar @user_cleanRead_files; $i++){
		push (@array, abs_path($user_cleanRead_files[$i]));
	}	
	my $user_cleanRead_files_ref = \@array;
	&fastq_files($user_cleanRead_files_ref);
}

sub check_options {
	
	##################################################
	## Get some info about the files names and tags ##
	##################################################
	
	## Check if a previous DOMINO parameters file exists
	if (-e $param_Detail_file_markers) { File::Copy::move($param_Detail_file_markers, $param_Detail_file_markers."_old_".$random_number); }
	if (-e $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $mapping_markers_errors_details."_old_".$random_number); }
	if (-e $mapping_parameters) { File::Copy::move($mapping_parameters, $mapping_parameters."_old_".$random_number); }
	
	## Behaviour DOMINO
	if (!$behaviour) { DOMINO::printError("\nPlease choose a development module for DOMINO between selection/discovery...\n"); DOMINO::dieNicely(); }
	if (!$option) { DOMINO::printError("\nPlease provide an option for DOMINO...\n"); DOMINO::dieNicely(); }
	if ($behaviour eq 'selection') {
		$select_markers=1;
		if ($option eq "DOMINO_files" || $option eq "user_assembly_contigs" || $option eq "genome") {
			DOMINO::printError("\nThe DOMINO development module SELECTION is not yet available for the option $option...\n"); DOMINO::dieNicely();
	}} elsif ($behaviour eq 'discovery') {
		if ($option eq "RADseq") { DOMINO::printError("\nThe DOMINO development module DISCOVERY is not suitable for the option $option...\n"); DOMINO::dieNicely(); }
		if (!$window_size_CONS_range || !$window_size_VARS_range) {
			unless ($option eq "msa_alignment") { DOMINO::printError("\nMandatory options are missing...\n"); DOMINO::dieNicely(); }
		}
		$identify_markers=1;
	} else { DOMINO::printError("\nPlease choose between selection/discovery...\n"); DOMINO::dieNicely(); }
	
	if ($option eq "DOMINO_files") {
		if (!$input_type) { DOMINO::printError("-type_input option is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
	} elsif ($option eq "user_assembly_contigs") {
		if (scalar @user_contig_files == 0) { DOMINO::printError("Contig assembled files were not provided...\nPlease provide several valid contig FASTA files for each taxa...."); DOMINO::dieNicely(); }
		if (scalar @user_cleanRead_files == 0) { DOMINO::printError("Clean Read files were not provided...\nPlease provide several FASTQ files for each taxa...."); DOMINO::dieNicely(); }
		if (!$input_type) { DOMINO::printError("input_type option is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
	} elsif ($option eq "genome") {
		if (!$genome_fasta) { DOMINO::printError("\nNo genome fasta file was provided...\nPlease provide a valid contig FASTA file to use as a reference using the option -genome_fasta [file]...."); DOMINO::dieNicely(); }
		if (!$input_type) { DOMINO::printError("Option -type_input is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
	} elsif ($option eq "msa_alignment") {
		if (!$msa_file and !$msa_fasta_folder) { DOMINO::printError("\nNo file or folder provided...\n"); DOMINO::dieNicely(); }
	} elsif ($option eq "RADseq") {
		$option = "msa_alignment"; $radseq_like_data = 1;
		if (!$msa_file) { DOMINO::printError("Exiting the script. No file provided. \nUse the option -RADseq_file [file]..\n"); DOMINO::dieNicely(); }
		push (@{ $domino_params{'marker'}{'RADseq'} }, "YES");
	} else { DOMINO::printError("\nOption provided is not known...\n"); DOMINO::dieNicely();}
	
	if (!$MID_taxa_names) {
	unless ($option eq "msa_alignment" || $option eq "RADseq") {
		DOMINO::printError("\nThe option -taxa_names option is missing...\n\nPlease provide it or DOMINO would not continue the process...\n"); DOMINO::dieNicely();
	}} else {
		$MID_taxa_names =~ s/\s+/\,/g;
		my @MID_name_array = split (",", $MID_taxa_names);
		push (@{$domino_params{"marker"}{'taxa_string'}}, $MID_taxa_names);
		for (my $j = 0; $j < scalar @MID_name_array; $j++) {
			push (@{$domino_files{'taxa'}{'user_Taxa'}}, $MID_name_array[$j]);
			push (@{$domino_files{$MID_name_array[$j]}{'taxa'}}, $MID_name_array[$j]);
			$number_sp++;
		}
		push (@{$domino_params{"mapping"}{'number_sp'}}, $number_sp);
	}
	if ($option eq "genome") {$number_sp++;} ## when reference genome provided, the ref taxa also counts
	
	## Start the Analysis
	print "\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" DOMINO Molecular Marker Development Stage ","#"); DOMINO::printHeader("","#"); print "\n"; DOMINO::printHeader("","+");  DOMINO::printHeader(" Analysis Started ","+");  DOMINO::printHeader("","+"); 
	DOMINO::printDetails("Starting the process: [ ".(localtime)." ]\n\n", $mapping_parameters, $param_Detail_file_markers);
	
	if ($keepbam) { print "\nKeepbam option not yet implemented....\nSorry for that...\n"; }
	
	##############################################################
	##		Checking and Printing user options					##
	##############################################################
	print "\n\n"; DOMINO::printHeader(" Input File and Parameter Preprocessing ","#"); print "\n";
	
	# Control if missing options
	if (!$variable_positions_user_range and !$variable_divergence) {
		DOMINO::printError("Exiting the script. A range for variable positions or a minimum divergence is missing.\n Use the option -VP|--variable_positions [min::max] or -VD|--variable_divergence [float number]..\n"); DOMINO::dieNicely();
	}
	unless (!$variable_divergence) { 
		if ($variable_divergence =~ /.*\,.*/) { $variable_divergence =~ s/\,/\./; }
		if ($variable_divergence < 0) {$variable_divergence = 0.000000000000000000000000000000001;} ## Set a very small value if -VD 0 
	}
	
	#############################################################################################################################
	## checking if any option is missing and using default
	if (!$num_proc_user) { $num_proc_user = 2; }
	if (!$rdgopen) { $rdgopen = 5; } 			## Bowtie Defaults
	if (!$rdgexten) { $rdgexten = 3; } 			## Bowtie Defaults
	if (!$rfgopen) { $rfgopen = 5; } 			## Bowtie Defaults
	if (!$rfgexten) { $rfgexten = 3; } 			## Bowtie Defaults
	if (!$mis_penalty) { $mis_penalty = 4;} 	## Bowtie Defaults
	if (!$SLIDING_user) { $SLIDING_user = 10; } ## Default DOMINO sliding interval for marker search
	if (!$cigar_pct) { $cigar_pct = 10; }
	if (!$window_var_CONS) { $window_var_CONS = 1;}
	if (!$level_significance_coverage_distribution) { $level_significance_coverage_distribution = 1e-05; }
	if (!$missing_allowed) { $missing_allowed = 0.1;} ## Def. 0.1: When looking for a marker if 10% of the length is missing for any specie, discard 
	if (!$totalContigs2use4markers) {$totalContigs2use4markers = 20000;} ## use by default the largest 20.000 contigs
	if (!$subset_offset_user) {$subset_offset_user = 50;} ## Split subsets into 50 contigs to avoid collapsing RAM
	
	## get optional
	my $answer_dnaSP = 0; if ($dnaSP_flag) { $answer_dnaSP++; $polymorphism_user=1;}
	my $answer_PV = 0; if ($polymorphism_user) { $answer_PV++; }
	my $BowtieLocal = 0; if ($bowtie_local) { $BowtieLocal++; } 
	my $mapContigFiles = 0; if ($map_contig_files) { $mapContigFiles++; }
	my $LowCoverageData = 0; if ($DOMINO_simulations) { $LowCoverageData++; }
	
	my ($variable_positions_user_min, $variable_positions_user_max);
	my ($window_size_CONS_min, $window_size_CONS_max);
	my ($window_size_VARS_min, $window_size_VARS_max);
	
	## Get ranges
	if ($variable_positions_user_range) {
		if ($variable_positions_user_range =~ m/.*\:\:.*/) {
			($variable_positions_user_min, $variable_positions_user_max) = split("::", $variable_positions_user_range);
		} elsif ($variable_positions_user_range =~ m/.*\:(\d+)/) {
			DOMINO::printError("\nPlease provide the range using 2 pair of dots like 2::7\n\n"); DOMINO::dieNicely();
		} else { $variable_positions_user_min = $variable_positions_user_max = $variable_positions_user_range; }
		push (@{ $domino_params{'marker'}{'variable_positions_user_range'} }, $variable_positions_user_range);
		push (@{ $domino_params{'marker'}{'variable_positions_user_max'} }, $variable_positions_user_max);
		push (@{ $domino_params{'marker'}{'variable_positions_user_min'} }, $variable_positions_user_min);
	}
	
	# Range Conserved size
	if ($window_size_CONS_range) {
		if ($window_size_CONS_range =~ m/.*\:\:.*/) {
			($window_size_CONS_min, $window_size_CONS_max) = split("::", $window_size_CONS_range);
		} elsif ($window_size_CONS_range =~ m/.*\:(\d+)/) {
			DOMINO::printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
		} else { $window_size_CONS_min = $window_size_CONS_max = $window_size_CONS_range;}
		push (@{ $domino_params{'marker'}{'window_size_CONS_range'} }, $window_size_CONS_range);
		push (@{ $domino_params{'marker'}{'window_size_CONS_max'} }, $window_size_CONS_max);
		push (@{ $domino_params{'marker'}{'window_size_CONS_min'} }, $window_size_CONS_min);
	}
	
	## Range Variable size
	if ($window_size_VARS_range) {
		if ($window_size_VARS_range =~ m/.*\:\:.*/) {
			($window_size_VARS_min, $window_size_VARS_max) = split("::", $window_size_VARS_range);
		} elsif ($window_size_VARS_range =~ m/.*\:(\d+)/) {
			DOMINO::printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
		} else { $window_size_VARS_min = $window_size_VARS_max = $window_size_VARS_range;}
		push (@{ $domino_params{'marker'}{'window_size_VARS_range'} }, $window_size_VARS_range);
		push (@{ $domino_params{'marker'}{'window_size_VARS_max'} }, $window_size_VARS_max);
		push (@{ $domino_params{'marker'}{'window_size_VARS_min'} }, $window_size_VARS_min);
	}
	## Variable Sliding window
	if ($behaviour eq "discovery") {
		if (!$VAR_inc) {
			my $difference_VAR = $window_size_VARS_max - $window_size_VARS_min;
			if ( $difference_VAR > 500) {  		$VAR_inc = 10; $CONS_inc = 5;
			} elsif ($difference_VAR > 300) { 	$VAR_inc = 5; $CONS_inc = 5;
			} elsif ($difference_VAR > 200) { 	$VAR_inc = 3; $CONS_inc = 3;
			} elsif ($difference_VAR > 100) { 	$VAR_inc = 2; $CONS_inc = 2;
			} else { $VAR_inc = 1; $CONS_inc = 1; }
		}
		## Conserved Sliding window
		if (!$CONS_inc) {
			# Set step for CONS range
			my $difference_CONS = $window_size_CONS_max - $window_size_CONS_min;
			if ( $difference_CONS >= 20) { $CONS_inc = 2; } else { $CONS_inc = 1; }
	}}
	
	# MCT
	if (!$minimum_number_taxa_covered) { $minimum_number_taxa_covered = $number_sp;  ## Force to be all the taxa
	} else { if ($minimum_number_taxa_covered > $number_sp) { DOMINO::printError("Minimum number of covered taxa (MCT) is bigger than the number of taxa provided...\n"); DOMINO::dieNicely(); }} 
	
	## Get all parameters
	## push parameters
	push (@{ $domino_params{'mapping'}{'option'} }, $option);
	push (@{ $domino_params{'mapping'}{'mapping_markers_errors_details'} }, $mapping_markers_errors_details);
	push (@{ $domino_params{'mapping'}{'mapping_parameters'} }, $mapping_parameters);		
	push (@{ $domino_params{'mapping'}{'cpu'} }, $num_proc_user);
	push (@{ $domino_params{'mapping'}{'type_input'} }, $input_type);
	push (@{ $domino_params{'mapping'}{'rdgopen'} }, $rdgopen);
	push (@{ $domino_params{'mapping'}{'rdgexten'} }, $rdgexten);
	push (@{ $domino_params{'mapping'}{'rfgopen'} }, $rfgopen);
	push (@{ $domino_params{'mapping'}{'rfgexten'} }, $rfgexten);
	push (@{ $domino_params{'mapping'}{'mis_penalty'} }, $mis_penalty);
	push (@{ $domino_params{'mapping'}{'level_significance_coverage_distribution'} }, $level_significance_coverage_distribution);
	push (@{ $domino_params{'mapping'}{'folder'} }, $align_dirname);	
	push (@{ $domino_params{'mapping'}{'dnaSP'} }, $answer_dnaSP);
	push (@{ $domino_params{'mapping'}{'polymorphism'} }, $answer_PV);
	push (@{ $domino_params{'mapping'}{'bowtie_local'} }, $BowtieLocal);
	push (@{ $domino_params{'mapping'}{'map_contig_files'} }, $mapContigFiles);
	push (@{ $domino_params{'mapping'}{'low_coverage_data'} }, $LowCoverageData);

	push (@{ $domino_params{'marker'}{'cpu'} }, $num_proc_user);
	push (@{ $domino_params{'marker'}{'variable_divergence'} }, $variable_divergence);
	push (@{ $domino_params{'marker'}{'behaviour'} }, $behaviour);
	push (@{ $domino_params{'marker'}{'option'} }, $option);
	push (@{ $domino_params{'marker'}{'folder'}}, $marker_dirname);
	push (@{ $domino_params{'marker'}{'SLIDING_user'} }, $SLIDING_user);
	push (@{ $domino_params{'marker'}{'cigar_pct'} }, $cigar_pct);
	push (@{ $domino_params{'marker'}{'V-SI_inc'} }, $VAR_inc);
	push (@{ $domino_params{'marker'}{'C-SI_inc'} }, $CONS_inc);
	push (@{ $domino_params{'marker'}{'subset_offset_user'} }, $subset_offset_user);
	push (@{ $domino_params{'marker'}{'totalContigs2use4markers'} }, $totalContigs2use4markers);
	push (@{ $domino_params{'marker'}{'missing_allowed'} }, $missing_allowed);
	push (@{ $domino_params{'marker'}{'window_var_CONS'} }, $window_var_CONS);
	push (@{ $domino_params{'marker'}{'totalContigs2use4markers'} }, $totalContigs2use4markers);
	push (@{ $domino_params{'marker'}{'subset_offset_user'} }, $subset_offset_user);
	push (@{ $domino_params{'marker'}{'MCT'} }, $minimum_number_taxa_covered);
	push (@{ $domino_params{'marker'}{'folder'} }, $marker_dirname);
	push (@{ $domino_params{'marker'}{'number_sp'} }, $number_sp);
	
	&debugger_print("DOMINO Parameters");&debugger_print("Ref", \%domino_params);
	#############################################################################################################################
}

sub check_previous {
	
	#############################################################################################################################
	## Check parameters previous run ## retrieve information generated
	my %domino_files_dump = %{ DOMINO::get_DOMINO_files($folder_abs_path."/", "mapping") };
	my %domino_params_dump = %{ DOMINO::get_parameters($folder_abs_path."/", "mapping") };
		#print Dumper \%domino_files_dump;
		#print Dumper \%domino_params_dump;
		
	my $undef_mapping=0;
	
	if (!$domino_params_dump{"mapping"}) { $undef_mapping++; } else {

	## Check parameters
	foreach my $keys (keys %{ $domino_params_dump{'mapping'} }) {
		next if ($keys eq "date"); next if ($keys eq "folder");
		next if ($keys eq "mapping_markers_errors_details");
		next if ($keys eq "mapping_parameters"); 
		next if ($keys eq "dump_file"); next if ($keys eq "cpu");
		next if ($keys eq "number_sp");
		
		my $prev = $domino_params_dump{'mapping'}{$keys}[0];
		my $curr = $domino_params{'mapping'}{$keys}[0];
		#print "Keys: $keys Curr: $curr Prev: $prev\n"; 
		unless ($prev eq $curr ) {
			$undef_mapping++; DOMINO::printError("There is difference: $keys $curr =/= $prev\n"); ## test
	}}
	&debugger_print("DOMINO params dump");&debugger_print("Ref", \%domino_files_dump);
	
	## Check files generated
	if ($genome_fasta) {		
		my $ref_genome_id;
		if ($genome_fasta =~/.*id-(.*)\.fasta/) {$ref_genome_id=$1;} else {$ref_genome_id="genome";}		
		my $profile = "PROFILE::Ref:$ref_genome_id";		
		foreach my $ref_taxa ( keys %domino_files ) {
			next if $ref_taxa eq 'taxa';
			next if $ref_taxa eq 'genome';	
			if ($domino_files_dump{$ref_taxa}{'taxa'}) {
				foreach my $taxa ( keys %domino_files ) {
					next if $ref_taxa eq $taxa; 
					next if $taxa eq 'taxa';
					unless ( $domino_files_dump{$ref_taxa}{$profile} ) {
						$undef_mapping++; &printError("There is not a profile folder for $ref_taxa vs $taxa ...\n");
		}}} else {$undef_mapping++; &printError("There is not a taxa name $ref_taxa in the previous run ...\n");
	}}} else {
		foreach my $ref_taxa ( keys %domino_files ) {
			next if $ref_taxa eq 'taxa';
			if ($domino_files_dump{$ref_taxa}{'taxa'}) {
				## Check for file
				unless ($domino_files_dump{$ref_taxa}{'contigs'}[0]){
					DOMINO::printError("There is not a contig file for $ref_taxa ...\n"); $undef_mapping++;
				}
				foreach my $taxa ( keys %domino_files ) {
					next if $ref_taxa eq $taxa; next if $taxa eq 'taxa';
					unless ( $domino_files_dump{$ref_taxa}{'PROFILE::Ref:'.$taxa} ) {
						$undef_mapping++; DOMINO::printError("There is not a profile folder for $ref_taxa vs $taxa ...\n");
		}}} else {$undef_mapping++; DOMINO::printError("There is not a taxa name $ref_taxa in the previous run ...\n"); 
	}}}}

	if ($undef_mapping > 0) {
		undef $avoid_mapping;
		DOMINO::printDetails("+ Although option -No_Profile_Generation was provided, it would be done again as parameters do not match with the available mapping folder...\n",$mapping_parameters, $param_Detail_file_markers);
		return ("NO");
	} else {
		DOMINO::printDetails("+ A previous profile has been generated with the same parameters and details...\n",$mapping_parameters, $param_Detail_file_markers);
		#%domino_files = %domino_files_dump;
		
		## Dump info up to now to a file 
		undef $domino_files_dump{"taxa"};
		foreach my $keys (keys %domino_files_dump) {
			foreach my $subkeys (keys %{ $domino_files_dump{$keys} }) {
			push (@{ $domino_files{$keys}{$subkeys} },  $domino_files_dump{$keys}{$subkeys}[0]);
		}}
		%domino_files = %{ DOMINO::get_uniq_hash(\%domino_files) };
		
		## Dump info up to now to a file 
		undef $domino_params{"mapping"}{"mapping_markers_errors_details"};
		$domino_params{"mapping"}{"mapping_markers_errors_details"}[0] = $mapping_markers_errors_details;
		
		undef $domino_params{"marker"}{"folder"};
		$domino_params{"marker"}{"folder"}[0] = $marker_dirname;

		undef $domino_params{"mapping"}{"folder"};
		$domino_params{"mapping"}{"folder"}[0] = $domino_params_dump{"mapping"}{"folder"}[0];

		undef %domino_params_dump;

		if (!$number_sp) {
			$number_sp = $domino_params{'marker'}{'number_taxa'}[0]; 
			$minimum_number_taxa_covered = $domino_params{'marker'}{'MCT'}[0];
			$MID_taxa_names = $domino_params{'marker'}{'taxa_string'}[0];
	}}
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
	#############################################################################################################################
	return ("YES");
}

sub print_options {
	#############################################################################################################################
	## Print Different options
	if (!$avoid_mapping) {
		if (-d $align_dirname) { ## Checking if the Directory already exists because a previous analysis
			File::Copy::move ($align_dirname, $align_dirname."_old_".$random_number);
			DOMINO::printDetails("+ Changing an old folder named as $align_dirname to $align_dirname"."_old_"."$random_number...OK\n", $mapping_parameters, $param_Detail_file_markers);     						
		}
		mkdir $align_dirname, 0755;
		if ($map_contig_files) { DOMINO::printDetails("+ Contig assembled files would be mapped...OK\n", $mapping_parameters, $param_Detail_file_markers); 
		} elsif ($genome_fasta) { DOMINO::printDetails("+ Clean reads files would be mapped to the genome provided...OK\n", $mapping_parameters, $param_Detail_file_markers); 
		} elsif ($msa_file) { DOMINO::printDetails("+ Parse the alingment file...OK\n", $mapping_parameters, $param_Detail_file_markers); 
		} elsif ($msa_fasta_folder) { DOMINO::printDetails("+ Parse the alingment folder provided...OK\n", $mapping_parameters);
		} else {  DOMINO::printDetails("+ Clean reads would be mapped...OK\n", $mapping_parameters, $param_Detail_file_markers); }
		DOMINO::printDetails("+ Alignment Directory: ".$align_dirname." ...OK\n", $mapping_parameters, $param_Detail_file_markers);
	}
	# Behaviour
	if ($behaviour eq 'selection') {		DOMINO::printDetails("+ DOMINO development module: Select informative markers ...OK\n", $mapping_parameters, $param_Detail_file_markers);
	} elsif ($behaviour eq 'discovery') { 	DOMINO::printDetails("+ DOMINO development module: Discover putative markers ...OK\n", $mapping_parameters, $param_Detail_file_markers); }
	
	# DOMINO option
	if ($radseq_like_data) { DOMINO::printDetails("+ Option: RADseq ...OK\n", $mapping_parameters, $param_Detail_file_markers);
	} else { DOMINO::printDetails("+ Option: $option ...OK\n", $mapping_parameters, $param_Detail_file_markers); }
	if ($map_contig_files) { 
		DOMINO::printDetails("+ Contig assembled files would be mapped...OK\n", $mapping_parameters, $param_Detail_file_markers);
		if ($keepbam) {DOMINO::printDetails("+ BAM files would be maintained for later visualization...\n", $mapping_parameters, $param_Detail_file_markers);}	
		DOMINO::printDetails("+ Type of file(s): $input_type ...OK\n+ Checking file(s):\n", $mapping_parameters, $param_Detail_file_markers);
	} elsif ($genome_fasta) { 
		DOMINO::printDetails("+ Clean reads files would be mapped to the genome provided...OK\n", $mapping_parameters, $param_Detail_file_markers); 
		if ($keepbam) {DOMINO::printDetails("+ SAM/BAM files would be maintained for later visualization...\n", $mapping_parameters, $param_Detail_file_markers);}	
		DOMINO::printDetails("+ Type of file(s): $input_type ...OK\n+ Checking file(s):\n", $mapping_parameters, $param_Detail_file_markers);
	} elsif ($option eq "msa_alignment") {
		if ($radseq_like_data) { 
			if ($msa_file =~ /.*\.loci/) { 
				DOMINO::printDetails("+ pyRAD loci data has been provided...OK\n", $mapping_parameters, $param_Detail_file_markers); 
				$pyRAD_file = 1;
			} elsif ($msa_file =~ /.*\.fa/) {
				DOMINO::printDetails("+ STACKS fasta file has been provided...OK\n", $mapping_parameters, $param_Detail_file_markers);  
				$stacks_file = 1;
		}} else { DOMINO::printDetails("+ Multiple sequence alignment has been provided...OK\n", $mapping_parameters, $param_Detail_file_markers); 
	}} else { 
		DOMINO::printDetails("+ Clean reads would be mapped to the contigs assembled...OK\n", $mapping_parameters, $param_Detail_file_markers); 
		if ($keepbam) {DOMINO::printDetails("+ SAM/BAM files would be maintained for later visualization...\n", $mapping_parameters, $param_Detail_file_markers);}	
		DOMINO::printDetails("+ Type of file(s): $input_type ...OK\n+ Checking file(s):\n", $mapping_parameters, $param_Detail_file_markers);
	}
	
	## Obtain the files for analysis: 
	if ($option eq 'user_assembly_contigs') {
		for (my $i = 0; $i < scalar @user_contig_files; $i++) {
			if ($user_contig_files[$i] eq ".DS_Store" || $user_contig_files[$i] eq "." || $user_contig_files[$i] eq ".." ) { next; }
			if ($user_contig_files[$i] =~ /.*id-(.*)\.contigs\.fasta/g) {
				if ($domino_files{$1}{'taxa'}) {
					&check_file($user_contig_files[$i], $1);
					push (@{ $domino_files{$1}{'contigs'} }, $user_contig_files[$i]); ## push the whole file path			
				} else { DOMINO::printError("Please check the tag for the file $user_contig_files[$i] \n...not matching any taxa name provided...\nDOMINO will not die here but take it into account..."); }
		}}
		if (!$map_contig_files) { &user_cleanRead_files();  push (@{ $domino_params{'mapping'}{'user_cleanRead_files'}}, 1);}
	} elsif ($option eq 'DOMINO_files') {
		unless ($avoid_mapping) {
			my $assembling_dirname = DOMINO::get_earliest("assembly", $folder_abs_path);
			if ($assembling_dirname eq 'assembly' || $assembling_dirname eq 'NO') {
				DOMINO::printError("No assembly folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely();
			}
			DOMINO::printDetails("+ DOMINO contigs would be retreived from: $assembling_dirname\n", $mapping_parameters, $param_Detail_file_markers);
			my $files_dir_ref = DOMINO::readDir($assembling_dirname);
			my @files = @$files_dir_ref;
			for (my $i = 0; $i < scalar @files; $i++) {
				if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
				if (-d $files[$i]) {next;}
				my $tmp_file_abs_path = $assembling_dirname."/".$files[$i];
				next if $files[$i] =~ /.*fasta\.fai$/;
				if ($files[$i] =~ /.*id\-(.*)\.contigs\.fasta/) {
					if ($domino_files{$1}{'taxa'}) {
						push (@{$domino_files{$1}{'contigs'}}, $tmp_file_abs_path);
			}}}
			## Obtain clean reads
			if (scalar @user_cleanRead_files == 0) { 
				&get_clean_files();	## use clean reads to map
			} else { &user_cleanRead_files();  push (@{ $domino_params{'mapping'}{'user_cleanRead_files'}}, 1);
			} ## user provides reads to map
	}} elsif ($option eq 'genome') {
		my $tmp = abs_path($genome_fasta);
	
		#push (@{ $domino_files{'genome'}{'contigs'}}, $tmp); 
		if ($genome_fasta =~/.*id-(.*)\.fasta/) {
			push (@{ $domino_files{$1}{'contigs'}}, $tmp);
			push (@{ $domino_files{$1}{'taxa'}}, "genome"); &check_file($tmp, $1);
		} else {
			push (@{ $domino_files{'genome'}{'contigs'}}, $tmp);
			push (@{ $domino_files{'genome'}{'taxa'}}, "1"); &check_file($tmp);
		}
		if (scalar @user_cleanRead_files == 0) {
			DOMINO::printError("Clean Read files were not provided...\nDOMINO would check in the output folder provided if there is a DOMINO_clean_data containing the FASTQ files for each taxa...."); 
			&get_clean_files();
		} else { 
			&user_cleanRead_files();
			push (@{ $domino_params{'mapping'}{'user_cleanRead_files'}}, 1);
		}	
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
		
	} elsif ($option eq "msa_alignment") {	
		if ($radseq_like_data) {
			my $rad_file_abs_path = abs_path($msa_file);
			push (@{$domino_files{'RADseq'}{'file'}}, $rad_file_abs_path);
			if ($pyRAD_file) { print "+ pyRAD data file: $rad_file_abs_path\n"; push (@{ $domino_params{'marker'}{'pyRAD'} }, 1);
			} elsif ($stacks_file) { print "+ STACKS file: $rad_file_abs_path\n"; 	push (@{ $domino_params{'marker'}{'STACKS'} }, 1);}		
			print "+ Checking file:\n";
			if (-f -e -r -s $rad_file_abs_path) {
				print "\tFile $rad_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
			} else { DOMINO::printError("File provided is not valid...\nPlease provide a valid file as specified in the DOMINO manual...."); DOMINO::dieNicely();
		}} elsif ($msa_file) {
			my $msa_file_abs_path = abs_path($msa_file);
			push (@{$domino_files{'MSA'}{'file'}}, $msa_file_abs_path);
			print "+ Multipe sequence alignment file provided: $msa_file_abs_path\n";		
			print "+ Checking file:\n";
			if (-f -e -r -s $msa_file_abs_path) {
				print "\t-File $msa_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
				chdir $align_dirname; system("ln -s $msa_file_abs_path");
				push (@{ $domino_params{'marker'}{'MSA'} }, 1);
			} else { DOMINO::printError("MSA file provided is not valid...\nPlease provide a valid contig MSA file as specified in the DOMINO manual...."); DOMINO::dieNicely();
		}} elsif ($msa_fasta_folder) {
			my $msa_folder_abs_path = abs_path($msa_fasta_folder);
			my @name_msa_folder = split("/",$msa_fasta_folder);
			push (@{$domino_files{'MSA_folder'}{'folder'}}, $msa_folder_abs_path);
			print "+ Multipe sequence alignment fasta folder provided: $msa_folder_abs_path\n+ Checking file(s):\n";
			if (-d $msa_folder_abs_path) {
				print "\t- Folder $msa_folder_abs_path\n\t\t- Folder exists, is readable and non-zero character...OK\n";
				chdir $align_dirname; 
				my @array = split("/", $msa_folder_abs_path);
				my $array_files_fasta_msa_ref = DOMINO::readDir($msa_folder_abs_path);
				print "\t- Checking files in folder provided...\n";
				for (my $i=0; $i < scalar @$array_files_fasta_msa_ref; $i++) {
					if ($$array_files_fasta_msa_ref[$i] eq "." || $$array_files_fasta_msa_ref[$i] eq ".." || $$array_files_fasta_msa_ref[$i] eq ".DS_Store"  ) {next;}
					my $file_path = $msa_folder_abs_path."/".$$array_files_fasta_msa_ref[$i];
					unless (-f -e -r -s $file_path) {
						DOMINO::printError("File $file_path is not readable or empty. Please discarded from the folder...\n"); DOMINO::dieNicely();
					} else { push (@{$domino_files{'MSA_folder'}{'files'}}, $$array_files_fasta_msa_ref[$i]);
				}} print "\t\t- Files checked and everything seems OK...\n\n";
				push (@{ $domino_params{'marker'}{'MSA_folder'} }, 1);
	
			} else { DOMINO::printError("MSA folder provided is not valid...\n"); DOMINO::dieNicely(); }
		} else { DOMINO::printError("MSA folder or file is missing...\n"); DOMINO::dieNicely(); }
		
		if (!$MID_taxa_names) {
			DOMINO::printDetails("+ No option -taxa_names provided.\n", $mapping_parameters, $param_Detail_file_markers);
			DOMINO::printDetails("+ DOMINO would verify all the taxa available...\n", $mapping_parameters, $param_Detail_file_markers);
			push (@{$domino_files{'taxa'}{'user_Taxa'}}, "all");		
	}}
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
	unless (!$MID_taxa_names) {
		DOMINO::printDetails("\n\n+ Taxa to use for the DOMINO development of molecular markers:\n", $mapping_parameters, $param_Detail_file_markers);
		my @split = split(",",$MID_taxa_names);
		for (my $i=0; $i < scalar @split; $i++) {
			DOMINO::printDetails("\tName: $split[$i]\n", $mapping_parameters, $param_Detail_file_markers);
		}
	}
	if ($behaviour eq 'selection') {
		DOMINO::printDetails("\n+ Parameters for the selection of molecular markers:\n", $param_Detail_file_markers);
		DOMINO::printDetails("\t- Variable Length (VL): All the available length would be used for each region\n", $param_Detail_file_markers); 
	} elsif ($behaviour eq 'discovery') {
		unless ($option eq "msa_alignment") {
			DOMINO::printDetails("\n\n+ Parameters for the mapping of molecular markers:\n", $mapping_parameters);
			DOMINO::printDetails("\t- Alignment of the reads using: Bowtie2\n", $mapping_parameters);
			if ($bowtie_local) {
				DOMINO::printDetails("\t- Local Bowtie: ON", $mapping_parameters); 
			}
			DOMINO::printDetails("\t- Read Gap Open penalty (rdgopen): ".$rdgopen."\n", $mapping_parameters);
			DOMINO::printDetails("\t- Read Gap Extension penalty (rdgexten): ".$rdgexten."\n", $mapping_parameters);
			DOMINO::printDetails("\t- Reference Gap Open penalty (rfgopen): ".$rfgopen."\n", $mapping_parameters);
			DOMINO::printDetails("\t- Reference Gap Open penalty (rfgexten): ".$rfgexten."\n", $mapping_parameters);
			DOMINO::printDetails("\t- Mismath penalty: ".$mis_penalty."\n", $mapping_parameters);
			DOMINO::printDetails("\t- Significance Level Coverage Distribution (SLCD): ".$level_significance_coverage_distribution."\n", $mapping_parameters);
		}
		DOMINO::printDetails("\n+ Parameters for the development of molecular markers:\n", $param_Detail_file_markers);
		my $tmp1 = $domino_params{'marker'}{'window_size_CONS_min'}[0];
		my $tmp2 = $domino_params{'marker'}{'window_size_CONS_max'}[0];
		DOMINO::printDetails("\t- Conserved Length (CL): $tmp1 -- $tmp2 (bp)\n", $param_Detail_file_markers);
		DOMINO::printDetails("\t- Conserved Differences (CD): $window_var_CONS\n", $param_Detail_file_markers);
		
		my $tmp3 = $domino_params{'marker'}{'window_size_VARS_min'}[0];
		my $tmp4 = $domino_params{'marker'}{'window_size_VARS_max'}[0];
		DOMINO::printDetails("\t- Variable Length (VL): $tmp3 -- $tmp4 (bp)\n", $param_Detail_file_markers);
		DOMINO::printDetails("\t- Sliding window increment (SI): $SLIDING_user (bp)\n", $param_Detail_file_markers);
		DOMINO::printDetails("\t- Sliding window increment Variable region (V-SI_inc): $VAR_inc (bp)\n", $param_Detail_file_markers);
		DOMINO::printDetails("\t- Sliding window increment Conserved region (C-SI_inc): $CONS_inc (bp)\n", $param_Detail_file_markers);	
	}
	
	if ($variable_divergence) {
		DOMINO::printDetails("\t- Variable Divergence (VD): $variable_divergence\n", $param_Detail_file_markers);	
	} else {
	if ($domino_params{'marker'}{'variable_positions_user_min'}[0] == $domino_params{'marker'}{'variable_positions_user_max'}[0]) {
		DOMINO::printDetails("\t- Variable Positions (VP): $domino_params{'marker'}{'variable_positions_user_min'}[0] (bp)\n", $param_Detail_file_markers);			
	} elsif ($domino_params{'marker'}{'variable_positions_user_max'}[0] == 999) {
		$domino_params{'marker'}{'variable_positions_user_max'}[0] = 99999999;
		DOMINO::printDetails("\t- Variable Positions (VP): > $domino_params{'marker'}{'variable_positions_user_min'}[0] (bp)\n", $param_Detail_file_markers);			
	} 
	
	if ($domino_params{'marker'}{'variable_positions_user_min'}[0] == 0) {
		$domino_params{'marker'}{'variable_positions_user_min'}[0] = 1;
	} else { DOMINO::printDetails("\t- Variable Positions (VP): $domino_params{'marker'}{'variable_positions_user_min'}[0] -- $domino_params{'marker'}{'variable_positions_user_max'}[0] (bp)\n", $param_Detail_file_markers);	
	}}
	
	## Common markers parameters
	unless (!$MID_taxa_names) { DOMINO::printDetails("\t- Minimum number of covered taxa (MCT): ".$minimum_number_taxa_covered."\n", $param_Detail_file_markers); }
	if ($polymorphism_user) { 
		DOMINO::printDetails("\t- Polymorphic variants would be detected (PV)...OK\n", $param_Detail_file_markers);
	}
	if ($dnaSP_flag) { DOMINO::printDetails("\t- dnaSP option [ON]...OK\n", $param_Detail_file_markers); }
	## Others parameters
	DOMINO::printDetails("\n+ Miscellaneous parameters:\n",$mapping_parameters, $param_Detail_file_markers);
	DOMINO::printDetails("\t- Threads: ".$num_proc_user."\n",$mapping_parameters, $param_Detail_file_markers);
	unless ($avoidDelete_tmp_files) {
		DOMINO::printDetails("\t- Deleting of temporary files would be done ...OK\n", $mapping_parameters, $param_Detail_file_markers);
	} else { DOMINO::printDetails("\t- Deleting temporary files would be avoid ...OK\n", $mapping_parameters, $param_Detail_file_markers); }
	if (!$avoid_mapping) {
		## Print info about where to print info and error
		DOMINO::printDetails("\t- DM mapping parameters details would be print into file:\n\t\t$mapping_parameters...\n", $mapping_parameters);
		DOMINO::printDetails("\t- DM mapping errors occurred during the process would be print into file:\n\t\t$mapping_markers_errors_details...\n", $mapping_parameters);
		chdir $align_dirname;
	}
	$mapping_markers_errors_details = $folder_abs_path."/".$datestring."_Markers_ERROR.txt";
	DOMINO::printDetails("\t- DM markers parameters details would be print into file:\n\t\t$param_Detail_file_markers...\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- DM markers errors occurred during the process would be print into file:\n\t\t$mapping_markers_errors_details...\n\n", $param_Detail_file_markers);
	print "\n"; &time_log(); print "\n";
	#############################################################################################################################
}
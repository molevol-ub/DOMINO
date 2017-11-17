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
##      perl DM_MarkerScan_1.0.3.pl
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
##		[-SI|--sliding_interval int] 
##		[-dnaSP]
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

BEGIN {
	require DOMINO;
	require List::Uniq; use List::Uniq qw(uniq);
	require File::Copy;
	require File::Path; use File::Path qw(remove_tree);
	require Cwd; use Cwd qw(abs_path);  
	require Parallel::ForkManager;
	require Spreadsheet::WriteExcel;
}

use List::MoreUtils qw(firstidx);

## TODO: Control if too many CPU (>6) provided or too many files (>20), do not print too much info

##################################
##	Initializing some variables	##
##################################
my $domino_version = "v1.0.3 ## Revised 15-11-2017";
my (
## User options
$folder, $helpAsked, $avoidDelete_tmp_files, $num_proc_user, $window_var_CONS, 
$window_size_CONS_range, $variable_divergence, $bowtie_local, $variable_positions_user_range,
$window_size_VARS_range, $input_type, $manual, $cigar_pct, $rdgopen, $rdgexten, $rfgexten, 
$rfgopen, $MID_taxa_names, $option, $mis_penalty, $msa_fasta_folder, $polymorphism_user,
$level_significance_coverage_distribution, $map_contig_files, $missing_allowed, $keepbam, 
$version, $DOMINO_simulations, $minimum_number_taxa_covered, $avoid_mapping, $further_information,
@user_cleanRead_files, @user_contig_files, $msa_file, $behaviour, $select_markers, $identify_markers,
$debugger, $helpAsked1, $VAR_inc, $CONS_inc, $option_all,

## absolute path
@contigs_fasta_file_abs_path, @clean_fastq_file_abs_path, # toDiscard

## others
%domino_files, %domino_params, $step_time, %discard_contigs, $pyRAD_file, $stacks_file, $radseq_like_data, 
$number_sp, $genome_fasta, $scripts_path, $dnaSP_flag, $SLIDING_user, %mapping_contigs
);

my %ambiguity_DNA_codes = (
"R" => ["A","G"], "Y" => ["C","T"], "K" => ["G","T"], "M" => ["A","C"], "S" => ["G","C"], 
"W" => ["A","T"], "B" => ["G","C","T"], "D" => ["A","G","T"], "H" => ["A","C","T"], 
"V" => ["A","C","G"], "N" => ["A","C","G","T"]);

my %ambiguity_DNA_codes_reverse;
foreach my $keys (sort keys %ambiguity_DNA_codes) {
	my @array = sort @{$ambiguity_DNA_codes{$keys}};
	my $string = join "",@array;
	$ambiguity_DNA_codes_reverse{$string} = $keys;	
}
	
######################
## Get user options	##
######################
GetOptions(
	"h" => \$helpAsked1,
	"help" => \$helpAsked,
	"man" => \$manual,
	"v|version" => \$version, 
	"MoreInfo" => \$further_information,

	"DM|development_module=s" => \$behaviour,

	"o|outputFolder=s" => \$folder,
	"p|number_cpu:i" => \$num_proc_user, ## default 2

	"option=s" => \$option,
	"type_input=s" => \$input_type,
	"genome_fasta=s" => \$genome_fasta,
	"user_contig_files:s" => \@user_contig_files,
	"user_cleanRead_files:s" => \@user_cleanRead_files,
	"RADseq_file=s" => \$msa_file, ## Loci data enters as msa_file 
	"msa_file=s" => \$msa_file, ## Multiple alignments in one phylip file
	"msa_folder=s" => \$msa_fasta_folder, ## Multiple files in phylip or fasta

	"rdgopen|read_gap_open_penalty=i" => \$rdgopen, #5
	"rdgexten|read_gap_extension_penalty=i" => \$rdgexten, #3
	"rfgopen|ref_gap_open_penalty=i" => \$rfgopen, #5
	"rfgexten|ref_gap_extension_penalty=i" => \$rfgexten, #3
 	"mp|mismatch_penalty=i" => \$mis_penalty, #4
	"PV|polymorphism" => \$polymorphism_user,

	"taxa_names=s" => \$MID_taxa_names,
	"CD|conserved_differences=i" => \$window_var_CONS, #1 ## variations_conserved_region
	"CL|conserved_length=s" => \$window_size_CONS_range, ##size_conserved_region
	"VD|variable_divergence=f" => \$variable_divergence, ## MINIMUN_variation_percentage: Valor 0-1
	"VL|variable_length=s" =>  \$window_size_VARS_range, #400::600 ## size_variable_region
	"VP|variable_positions=s" => \$variable_positions_user_range, ## 1::5
	"MCT|minimum_number_taxa_covered:i" => \$minimum_number_taxa_covered,

	"NPG|No_Profile_Generation" => \$avoid_mapping,
	"TempFiles" => \$avoidDelete_tmp_files,

	"low_coverage_data" => \$DOMINO_simulations,
	"map_contig_files" => \$map_contig_files,
	"max_SoftClipping=i" => \$cigar_pct,
	"MPA|missing_perct_allowed:i" => \$missing_allowed,
 	"SLCD|significance_level_coverage_distribution=s" => \$level_significance_coverage_distribution,  ## -lscd 1e-05
 	"bowtie_local" => \$bowtie_local,
 	"keep_bam_file" => \$keepbam,
 	"dnaSP" => \$dnaSP_flag,
 	"SI|sliding_interval:i" => \$SLIDING_user,
 	"V-SI_inc:i" => \$VAR_inc,
 	"C-SI_inc:i" => \$CONS_inc,

 	"all" => \$option_all, 	
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

DM_MarkerScan_1.0.2.pl

=back
		
=head1 VERSION

=over 2

DOMINO v1.0.2

=back
	
=head1 SYNOPSIS

=over 2

=item B<>
	
perl DM_MarkerScan_1.0.2.pl

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

[-MPA|--missing_perct_allowed float_value] [-MCT|--minimum_number_taxa_covered int_value] [-CD|--conserved_differences int_value] [-PV|--polymorphism] [-VP|--variable_positions range] [-SI|--sliding_interval int] [-dnaSP]

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

=item B<-all [Default Off]>

Use all the contigs generated during the assembly. By default and in order to speed the computation, DOMINO would only use the largest 50000 contigs generated.
	
=item B<-TempFiles [Default Off]>
	
Keep all intermediate files.

=item B<>

=item B<#####################################>

=item B<##### Command Line Examples #########>

=item B<#####################################>

=item B<DOMINO files: single end -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<DOMINO files: single end, No Mapping  -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 500 -CD 1 -NPG -MCT 2 -MPA 25 -DM discovery 

=item B<DOMINO files: paired-end -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option DOMINO_files
 -type_input pair_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 -DM discovery 

=item B<User provides contigs and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/Dmelanogaster.contigs.fasta -user_contig_files path_to_file2/Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/Dyakuba.contigs.fasta -user_cleanRead_files Dmelanogaster.clean.fastq 
 -user_cleanRead_files Dsimulans.clean.fastq -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<User provides contigs and reads (paired-end) -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option user_assembly_contigs -type_input pair_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -user_cleanRead_files reads_id-Dmelanogaster.clean.R1.fastq -user_cleanRead_files reads_id-Dmelanogaster.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dsimulans.clean.R1.fastq -user_cleanRead_files reads_id-Dsimulans.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dyakuba.clean.R1.fastq -user_cleanRead_files reads_id-Dyakuba.clean.R2.fastq 
 -DM discovery 

=item B<User provides contigs but no reads, and specifies to map contigs vs contigs  -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -map_contig_files -DM discovery 

=item B<User provides a reference genome and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option genome -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -genome_fasta path_to_genomes_folder/NCBI_id-Dpseudobscura.fasta
 -user_cleanRead_files Dmelanogaster.clean.fastq -user_cleanRead_files Dsimulans.clean.fastq 
 -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<MSA alignment: A single MSA in PHYLIP -- Single File -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/ -msa_file file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single MSA in FASTA -- Single File -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/ -msa_file file.fasta 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Discovery>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Selection>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in PHYLIP -- Folder -- Selection>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in FASTA MSA -- Folder -- Selection>

 perl DM_MarkerScan_v1.0.2.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in pyRAD format -- File -- Selection>

 perl DM_MarkerScan_v1.0.2.pl -option RADseq -o test/  -RADseq_file output.loci
 -taxa_names ind1,ind3,ind5,ind8,ind9 -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in STACKS format -- File -- Selection>

 perl DM_MarkerScan_v1.0.2.pl -option RADseq -o test/  -RADseq_file output.fa
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

07 - 11 - 2016

=back

=head1 CITATION

=over 2

Bioinformatics first published online August 16, 2016 
doi:10.1093/bioinformatics/btw534 

=back

=cut

## Check mandatory options
if (!$folder) { DOMINO::dieNicely(); } 

#######################################
###	Initialise some PATH variables	### 
#######################################
## Directory names
my $folder_abs_path = abs_path($folder);
unless (-d $folder_abs_path) { mkdir $folder_abs_path, 0755; }
## Setting a timestamp
my $random_number = int(rand(100));
my $datestring = strftime "%Y%m%d%H%M", localtime;
my $date = localtime;

## mapping dirname
my $align_dirname = $folder_abs_path."/".$datestring."_DM_mapping"; 
my $mapping_parameters = $folder_abs_path."/".$datestring."_Mapping-Parameters.txt";
my $mapping_markers_errors_details = $folder_abs_path."/".$datestring."_Mapping_ERROR.txt";
my $msa_dirname = $align_dirname."/MSA_files";

## Markers dirname
my $marker_dirname = $folder_abs_path."/".$datestring."_DM_markers";
my $param_Detail_file_markers = $folder_abs_path."/".$datestring."_Markers-Parameters.txt";

my $start_time = $step_time = time;
my $pipeline_path = abs_path($0);
my @script_path_array = split ("/", $pipeline_path);
for (my $j = 0; $j < $#script_path_array; $j++) { $scripts_path .= $script_path_array[$j]."/"; } 

## General variables
my $samtools_path = $scripts_path."samtools-1.3.1/samtools";
my $bowtie_path = $scripts_path."bowtie2-2.2.9/";
my $BLAST = $scripts_path."NCBI_BLAST/";

## Check if a previous DOMINO parameters file exists
if (-e $param_Detail_file_markers) { File::Copy::move($param_Detail_file_markers, $param_Detail_file_markers."_old_".$random_number); }
if (-e $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $mapping_markers_errors_details."_old_".$random_number); }
if (-e $mapping_parameters) { File::Copy::move($mapping_parameters, $mapping_parameters."_old_".$random_number); }

## Behaviour DOMINO
if (!$behaviour) { &printError("\nPlease choose a development module for DOMINO between selection/discovery...\n"); DOMINO::dieNicely(); }
if (!$option) { &printError("\nPlease provide an option for DOMINO...\n"); DOMINO::dieNicely(); }
if ($behaviour eq 'selection') {
	$select_markers=1;
	if ($option eq "DOMINO_files" || $option eq "user_assembly_contigs" || $option eq "genome") {
		&printError("\nThe DOMINO development module SELECTION is not yet available for the option $option...\n"); DOMINO::dieNicely();
}} elsif ($behaviour eq 'discovery') {
	if ($option eq "RADseq") { &printError("\nThe DOMINO development module DISCOVERY is not suitable for the option $option...\n"); DOMINO::dieNicely(); }
	if (!$window_size_CONS_range || !$window_size_VARS_range) {
		unless ($option eq "msa_alignment") { &printError("\nMandatory options are missing...\n"); DOMINO::dieNicely(); }
	}
	$identify_markers=1;
} else { &printError("\nPlease choose between selection/discovery...\n"); DOMINO::dieNicely(); }
push (@{ $domino_params{'general'}{'behaviour'} }, $behaviour);

# Control if missing options
if (!$variable_positions_user_range and !$variable_divergence) {
	&printError("Exiting the script. A range for variable positions or a minimum divergence is missing.\n Use the option -VP|--variable_positions [min::max] or -VD|--variable_divergence [float number]..\n"); DOMINO::dieNicely();
}
unless (!$variable_divergence) { 
	if ($variable_divergence =~ /.*\,.*/) { $variable_divergence =~ s/\,/\./; }
	if ($variable_divergence < 0) {$variable_divergence = 0.000000000000000000000000000000001;} ## Set a very small value if -VD 0 
}
push (@{ $domino_params{'marker'}{'variable_divergence'} }, $variable_divergence);


if ($option eq "DOMINO_files") {
	if (!$input_type) { &printError("-type_input option is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
} elsif ($option eq "user_assembly_contigs") {
	if (scalar @user_contig_files == 0) { &printError("Contig assembled files were not provided...\nPlease provide several valid contig FASTA files for each taxa...."); DOMINO::dieNicely(); }
	if (scalar @user_cleanRead_files == 0) { &printError("Clean Read files were not provided...\nPlease provide several FASTQ files for each taxa...."); DOMINO::dieNicely(); }
	if (!$input_type) { &printError("input_type option is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
} elsif ($option eq "genome") {
	if (!$genome_fasta) { &printError("\nNo genome fasta file was provided...\nPlease provide a valid contig FASTA file to use as a reference using the option -genome_fasta [file]...."); DOMINO::dieNicely(); }
	if (!$input_type) { &printError("Option -type_input is missing...\nPlease provide it in order to proceed with the computation...."); DOMINO::dieNicely(); }
} elsif ($option eq "msa_alignment") {
	if (!$msa_file and !$msa_fasta_folder) { &printError("\nNo file or folder provided...\n"); DOMINO::dieNicely(); }
} elsif ($option eq "RADseq") {
	$option = "msa_alignment";
	$radseq_like_data = 1;
	if (!$msa_file) { &printError("Exiting the script. No file provided. \nUse the option -RADseq_file [file]..\n"); DOMINO::dieNicely(); }
} else { &printError("\nOption provided is not known...\n"); DOMINO::dieNicely();}
push (@{ $domino_params{'general'}{'option'} }, $option);

##################################################
## Get some info about the files names and tags ##
##################################################
if (!$MID_taxa_names) {
unless ($option eq "msa_alignment" || $option eq "RADseq") {
	&printError("\nThe option -taxa_names option is missing...\n\nPlease provide it or DOMINO would not continue the process...\n"); DOMINO::dieNicely();
}} else {
	$MID_taxa_names =~ s/\s+/\,/g;
	my @MID_name_array = split (",", $MID_taxa_names);
	for (my $j = 0; $j < scalar @MID_name_array; $j++) {
		push (@{$domino_files{'taxa'}{'user_Taxa'}}, $MID_name_array[$j]);
		push (@{$domino_files{$MID_name_array[$j]}{'taxa'}}, 1);
		$number_sp++;
}}

## Start the Analysis
print "\n"; DOMINO::printHeader("","#"); DOMINO::printHeader(" DOMINO Molecular Marker Development Stage ","#"); DOMINO::printHeader("","#"); print "\n"; DOMINO::printHeader("","+");  DOMINO::printHeader(" Analysis Started ","+");  DOMINO::printHeader("","+"); 
DOMINO::printDetails("Starting the process: [ ".(localtime)." ]\n\n", $mapping_parameters, $param_Detail_file_markers);

##############################################################
##		Checking and Printing user options					##
##############################################################
print "\n\n"; DOMINO::printHeader(" Input File and Parameter Preprocessing ","#"); print "\n";

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

## Get ranges
my ($variable_positions_user_min, $variable_positions_user_max);
if ($variable_positions_user_range) {
	if ($variable_positions_user_range =~ m/.*\:\:.*/) {
		($variable_positions_user_min, $variable_positions_user_max) = split("::", $variable_positions_user_range);
	} elsif ($variable_positions_user_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 2::7\n\n"); DOMINO::dieNicely();
	} else { $variable_positions_user_min = $variable_positions_user_max = $variable_positions_user_range; }
	push (@{ $domino_params{'marker'}{'variable_positions_user_max'} }, $variable_positions_user_max);
	push (@{ $domino_params{'marker'}{'variable_positions_user_min'} }, $variable_positions_user_min);
}

# Range Conserved size
my ($window_size_CONS_min, $window_size_CONS_max);
if ($window_size_CONS_range) {
	if ($window_size_CONS_range =~ m/.*\:\:.*/) {
		($window_size_CONS_min, $window_size_CONS_max) = split("::", $window_size_CONS_range);
	} elsif ($window_size_CONS_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
	} else { $window_size_CONS_min = $window_size_CONS_max = $window_size_CONS_range;}
	push (@{ $domino_params{'marker'}{'window_size_CONS_max'} }, $window_size_CONS_max);
	push (@{ $domino_params{'marker'}{'window_size_CONS_min'} }, $window_size_CONS_min);
}

## Range Variable size
my ($window_size_VARS_min, $window_size_VARS_max);
if ($window_size_VARS_range) {
	if ($window_size_VARS_range =~ m/.*\:\:.*/) {
		($window_size_VARS_min, $window_size_VARS_max) = split("::", $window_size_VARS_range);
	} elsif ($window_size_VARS_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
	} else { $window_size_VARS_min = $window_size_VARS_max = $window_size_VARS_range;}
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
} else { if ($minimum_number_taxa_covered > $number_sp) { &printError("Minimum number of covered taxa (MCT) is bigger than the number of taxa provided...\n"); DOMINO::dieNicely(); }} 

## push default parameters
my $answer_PV; if ($polymorphism_user) { $answer_PV = 'YES'; } else { $answer_PV = 'NO';} 
my $BowtieLocal; if ($bowtie_local) { $BowtieLocal = 'YES'; } else { $BowtieLocal = 'NO';} 
my $mapContigFiles; if ($map_contig_files) { $mapContigFiles = 'YES'; } else { $mapContigFiles = 'NO';} 
my $LowCoverageData; if ($DOMINO_simulations) { $LowCoverageData = 'YES'; } else { $LowCoverageData = 'NO';} 
push (@{ $domino_params{'general'}{'cpu'} }, $num_proc_user);
push (@{ $domino_params{'general'}{'type_input'} }, $input_type);
push (@{ $domino_params{'mapping'}{'rdgopen'} }, $rdgopen);
push (@{ $domino_params{'mapping'}{'rdgexten'} }, $rdgexten);
push (@{ $domino_params{'mapping'}{'rfgopen'} }, $rfgopen);
push (@{ $domino_params{'mapping'}{'rfgexten'} }, $rfgexten);
push (@{ $domino_params{'mapping'}{'mis_penalty'} }, $mis_penalty);
push (@{ $domino_params{'mapping'}{'level_significance_coverage_distribution'} }, $level_significance_coverage_distribution);
push (@{ $domino_params{'mapping'}{'polymorphism'} }, $answer_PV);
push (@{ $domino_params{'mapping'}{'low_coverage_data'} }, $LowCoverageData);
push (@{ $domino_params{'mapping'}{'map_contig_files'} }, $mapContigFiles);
push (@{ $domino_params{'mapping'}{'bowtie_local'} }, $BowtieLocal);
push (@{ $domino_params{'marker'}{'SLIDING_user'} }, $SLIDING_user);
push (@{ $domino_params{'marker'}{'cigar_pct'} }, $cigar_pct);
push (@{ $domino_params{'marker'}{'V-SI_inc'} }, $VAR_inc);
push (@{ $domino_params{'marker'}{'C-SI_inc'} }, $CONS_inc);
&debugger_print("DOMINO Parameters");&debugger_print("Ref", \%domino_params);

## Check parameters previous run
if ($avoid_mapping) {
	my (%domino_files_dump, %domino_params_dump);
	my $undef_mapping=0;
	my $path_returned = DOMINO::get_earliest("mapping", $folder_abs_path);
	if ($path_returned eq 'NO') { $undef_mapping++; 
	} else { 		
		my ($file2dump, $file2dump_param);
		if ($path_returned =~ /(\d+)\_DM\_mapping$/) { 
			$file2dump = $path_returned."/".$1."_DUMP.txt";
			$file2dump_param = $path_returned."/".$1."_DUMP_param.txt";
		} 
		&debugger_print("Path: $path_returned");
		&debugger_print("DOMINO files: $file2dump");
		&debugger_print("DOMINO param: $file2dump_param");
		## DUMP.txt
		open (DUMP_IN, "$file2dump");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			push (@{ $domino_files_dump{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_IN);
		## DUMP_param.txt
		open (DUMP_PARAM, "$file2dump_param");
		while (<DUMP_PARAM>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			push (@{ $domino_params_dump{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_PARAM);

		## Check parameters
		foreach my $keys (keys %{ $domino_params_dump{'mapping'} }) {
			my $prev = $domino_params_dump{'mapping'}{$keys}[0];
			my $curr = $domino_params{'mapping'}{$keys}[0];
			unless ($prev eq $curr ) {
				$undef_mapping++; &printError("There is difference: $keys $curr =/= $prev\n"); ## test
		}}

		## Check files generated
		foreach my $ref_taxa ( keys %domino_files ) {
			next if $ref_taxa eq 'taxa';
			if ($domino_files_dump{$ref_taxa}{'taxa'}) {
				## Check for file
				unless ($domino_files_dump{$ref_taxa}{'contigs'}[0]){
					&printError("There is not a contig file for $ref_taxa ...\n"); $undef_mapping++;
				}
				foreach my $taxa ( keys %domino_files ) {
					next if $ref_taxa eq $taxa; next if $taxa eq 'taxa';
					unless ( $domino_files_dump{$ref_taxa}{'PROFILE::Ref:'.$taxa} ) {
						$undef_mapping++; &printError("There is not a profile folder for $ref_taxa vs $taxa ...\n");
		}}} else {$undef_mapping++; &printError("There is not a taxa name $ref_taxa in the previous run ...\n"); 
	}}}
	if ($undef_mapping > 0) {
		undef $avoid_mapping;
		DOMINO::printDetails("+ Although option -No_Profile_Generation was provided, it would be done again as parameters do not much with the available mapping folder...\n",$mapping_parameters, $param_Detail_file_markers);
	} else {
		DOMINO::printDetails("+ A previous profile has been generated with the same parameters and details...\n",$mapping_parameters, $param_Detail_file_markers);
		%domino_files = %domino_files_dump; 
		%domino_params = %domino_params_dump; 
		if (!$number_sp) {
			$number_sp = $domino_files{'number_taxa'}{'number'}[0]; 
			$minimum_number_taxa_covered = $domino_files{'number_taxa'}{'MCT'}[0];
			$MID_taxa_names = $domino_files{'taxa_string'}{'string'}[0];
}}}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

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
			} else { &printError("Please check the tag for the file $user_contig_files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
	}}
	if (!$map_contig_files) { &user_cleanRead_files(); }
} elsif ($option eq 'DOMINO_files') {
	unless ($avoid_mapping) {
		my $assembling_dirname = DOMINO::get_earliest("assembly", $folder_abs_path);
		if ($assembling_dirname eq 'assembly' || $assembling_dirname eq 'NO') {
			&printError("No assembly folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely();
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
		} else { &user_cleanRead_files(); 
		} ## user provides reads to map
}} elsif ($option eq 'genome') {
	my $tmp = abs_path($genome_fasta);
	push (@{ $domino_files{'genome'}{'contigs'}}, $tmp); 
	if ($genome_fasta =~/.*id-(.*)\.fasta/) {
		push (@{ $domino_files{'genome'}{'taxa'}}, $1); &check_file($tmp, $1);
	} else {
		push (@{ $domino_files{'genome'}{'taxa'}}, "1"); &check_file($tmp);
	}
	if (scalar @user_cleanRead_files == 0) {
		&printError("Clean Read files were not provided...\nDOMINO would check in the output folder provided if there is a DOMINO_clean_data containing the FASTQ files for each taxa...."); 
		&get_clean_files();
	} else { &user_cleanRead_files(); }	
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
	
} elsif ($option eq "msa_alignment") {	
	if ($radseq_like_data) {
		my $rad_file_abs_path = abs_path($msa_file);
		push (@{$domino_files{'RADseq'}{'file'}}, $rad_file_abs_path);
		if ($pyRAD_file) { print "+ pyRAD data file: $rad_file_abs_path\n"; 		
		} elsif ($stacks_file) { print "+ STACKS file: $rad_file_abs_path\n";}		
		print "+ Checking file:\n";
		if (-f -e -r -s $rad_file_abs_path) {
			print "\tFile $rad_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
		} else { &printError("File provided is not valid...\nPlease provide a valid file as specified in the DOMINO manual...."); DOMINO::dieNicely();
	}} elsif ($msa_file) {
		my $msa_file_abs_path = abs_path($msa_file);
		push (@{$domino_files{'MSA'}{'file'}}, $msa_file_abs_path);
		print "+ Multipe sequence alignment file provided: $msa_file_abs_path\n";		
		print "+ Checking file:\n";
		if (-f -e -r -s $msa_file_abs_path) {
			print "\t-File $msa_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
 			chdir $align_dirname; system("ln -s $msa_file_abs_path");
		} else { &printError("MSA file provided is not valid...\nPlease provide a valid contig MSA file as specified in the DOMINO manual...."); DOMINO::dieNicely();
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
					&printError("File $file_path is not readable or empty. Please discarded from the folder...\n"); DOMINO::dieNicely();
				} else { push (@{$domino_files{'MSA_folder'}{'files'}}, $$array_files_fasta_msa_ref[$i]);
			}} print "\t\t- Files checked and everything seems OK...\n\n";
		} else { &printError("MSA folder provided is not valid...\n"); DOMINO::dieNicely(); }
	} else { &printError("MSA folder or file is missing...\n"); DOMINO::dieNicely(); }
	
	if (!$MID_taxa_names) {
		DOMINO::printDetails("+ No option -taxa_names provided.\n", $mapping_parameters, $param_Detail_file_markers);
		DOMINO::printDetails("+ DOMINO would verify all the taxa available...\n", $mapping_parameters, $param_Detail_file_markers);
		push (@{$domino_files{'taxa'}{'user_Taxa'}}, "all");		
}}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
unless (!$MID_taxa_names) {
	DOMINO::printDetails("\n\n+ Taxa to use for the DOMINO development of molecular markers:\n", $mapping_parameters, $param_Detail_file_markers);
	foreach my $keys (sort keys %domino_files) { 
		if ($domino_files{$keys}{'taxa'}) {
			DOMINO::printDetails("\tName: $keys\n", $mapping_parameters, $param_Detail_file_markers);
}}}
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
	DOMINO::printDetails("\t- Conserved Length (CL): $window_size_CONS_min -- $window_size_CONS_max (bp)\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Conserved Differences (CD): $window_var_CONS\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Variable Length (VL): $window_size_VARS_min -- $window_size_VARS_max (bp)\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Sliding window increment (SI): $SLIDING_user (bp)\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Sliding window increment Variable region (V-SI_inc): $VAR_inc (bp)\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Sliding window increment Conserved region (C-SI_inc): $CONS_inc (bp)\n", $param_Detail_file_markers);	
}

if ($variable_divergence) {
	DOMINO::printDetails("\t- Variable Divergence (VD): $variable_divergence\n", $param_Detail_file_markers);	
} else {
if ($variable_positions_user_min == $variable_positions_user_max) {
	DOMINO::printDetails("\t- Variable Positions (VP): $variable_positions_user_min (bp)\n", $param_Detail_file_markers);			
} elsif ($variable_positions_user_max == 999) {
	$variable_positions_user_max = 99999999;
	DOMINO::printDetails("\t- Variable Positions (VP): > $variable_positions_user_min (bp)\n", $param_Detail_file_markers);			
} 

if ($variable_positions_user_min == 0) {
	$variable_positions_user_min = 1;
} else { DOMINO::printDetails("\t- Variable Positions (VP): $variable_positions_user_min -- $variable_positions_user_max (bp)\n", $param_Detail_file_markers);	
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
	
################################################################################################
################# 		Mapping/Alignment of the contigs 		################################
################################################################################################
if (!$avoid_mapping) {	
	if ($option ne "msa_alignment") {
		## We would use Bowtie2 for mapping the reads		
		DOMINO::printHeader("", "#");	DOMINO::printHeader(" Mapping Process started ", "#"); DOMINO::printHeader("", "#"); print "\n";
		
		#############################################################
		### Get Pre-assemble taxa read contigs of each taxa ### 
		#############################################################
		print "\n"; DOMINO::printHeader(" Get FASTQ files of the contigs generated ", "%"); print "\n";
		chdir $align_dirname; &debugger_print("Change dir to: ".$align_dirname);
		
		## Mapping of the reads, all taxa used as reference
		foreach my $reference_identifier (sort keys %domino_files) {
			unless ($domino_files{$reference_identifier}{'contigs'}) { next; }
			chdir $align_dirname; &debugger_print("Change dir to: ".$align_dirname);
	
			my (@sam_files, @clean_sam_files, @sorted_bam); 
				
			## Get the name and identifier of each reference fasta used	
			my @temp_contigs_name = split ("/", $domino_files{$reference_identifier}{'contigs'}[0]);
			my $contigs_fasta = $temp_contigs_name[$#temp_contigs_name];
			my $ref_Fasta = $reference_identifier.".fasta";
			
			###############################
			###	 Copy necessary files	### 
			###############################
			print "+ Generating symbolic links for necessary files\n"; 
			my $dir = $align_dirname."/".$reference_identifier; mkdir $dir, 0755; chdir $dir;
			system("ln -s $domino_files{$reference_identifier}{'contigs'}[0] $ref_Fasta");
			push (@{ $domino_files{$reference_identifier}{'dir'} }, $dir);
			## Generate a directory for each one
			print "+ Using as reference: $contigs_fasta\tID: $reference_identifier...OK\n";
			print "+ Generating a new directory $dir....OK\n\n";
			&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
		
			###################################
			###		 Index Contig file		### 
			###################################
			DOMINO::printHeader(" Indexing Contig File for mapping Reference ", "%");
			# Index contig reference file using Bowtie
			my $reference_tag = "reference_".$reference_identifier;
			print "- Reference: $ref_Fasta...\n";
			my $bowtie_index_call = $bowtie_path."bowtie2-build --threads $num_proc_user -f ".$ref_Fasta." ".$reference_tag;   
			&debugger_print("BOWTIE2 command: ".$bowtie_index_call);
			my $index_result = system ($bowtie_index_call);
			if ($index_result != 0) {
				&printError("Exiting the script. Some error happened when calling bowtie for indexing the file...\n"); DOMINO::dieNicely();
			} print "\n"; &time_log();	print "\n";
			
			###########################
			###	Align taxa Reads	### 
			###########################
			print "\n";	DOMINO::printHeader(" Aligning Reads Individually ", "%"); print "\n";
			my @clean_fastq_files;
			my @tmp_array;
			if (scalar @user_cleanRead_files > 0) { print "+ User clean reads files would be mapped\n";
			} elsif ($map_contig_files) { 			print "+ Contig files would be mapped\n";
			} else { 								print "+ Clean reads files would be mapped\n"; ## Map DOMINO clean reads
			}
			print "+ Obtain information of the reference sequences\n";
			my ($reference_hash_fasta_ref, $message) = DOMINO::readFASTA_hashLength($ref_Fasta); ## Obtain reference of a hash
			my $file2dump_seqs = $dir."/contigs_".$reference_identifier."_length.txt";
			push (@{ $domino_files{$reference_identifier}{"hash_reference_file"} }, $file2dump_seqs);
			DOMINO::printDump($reference_hash_fasta_ref,$file2dump_seqs,1);
			
			## Set DOMINO to use as many CPUs as provided
			## maximize process if > 10
			my ($split_CPU, $subprocesses); 
			if ($num_proc_user > 10) { 	$split_CPU = int($num_proc_user/$number_sp);  $subprocesses = $number_sp; 
			} else { 					$split_CPU = $num_proc_user; $subprocesses = 1; 			
			}			
			
			## Get files for later dump
			foreach my $reads_here (sort keys %domino_files) {
				unless ($domino_files{$reads_here}{'reads'}) { next; }
				my $file2dump = $dir."/dump_file_mapping_".$reads_here."_Ref_".$reference_identifier.".txt";
				#print "File to dump: ".$file2dump."\n";
				push (@{ $domino_files{$reads_here}{"DUMP_Mapping::Ref:".$reference_identifier} }, $file2dump);
			}
			
			my $pm_read_Reference =  new Parallel::ForkManager($subprocesses); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			my $total_subprocesses = $number_sp; my $count_subprocesses=0;
			print "+ Mapping process would be divided into $number_sp parts using up to $split_CPU/$num_proc_user CPUs provided\n";
			foreach my $reads (sort keys %domino_files) {
				unless ($domino_files{$reads}{'reads'}) { next; }
				$count_subprocesses++;
				print "\t- Mapping reads ($reads) vs reference ($reference_identifier): [$count_subprocesses/$total_subprocesses]\n";
				my $out_log_file = $dir."/reads_".$reads."-reference_".$reference_identifier."_logfile.txt";
				open (LOG, ">$out_log_file");
				chdir $domino_files{$reference_identifier}{'dir'}[0];
				print LOG "\nChange dir to: ".$domino_files{$reference_identifier}{'dir'}[0]."\n";
				my $pid = $pm_read_Reference->start() and next;
				my %domino_files_split_mapping;
				push(@{ $domino_files_split_mapping{$reads}{'LOG'}}, $out_log_file);

				## Mapping Parameters
				my $R_group_id = '--rg-id '.$reads;
				my $R_group_name = ' --rg '.$reads;
				my $threads = ' -p '.$split_CPU;
				my $mismatches = ' -N 1 --np 0'; ## Do not add penalty if read/ref got an ambiguous base
				my $read_gap_open = ' --rdg '.$rdgopen.','.$rdgexten;
				my $ref_gap_open = ' --rfg '.$rfgopen.','.$rfgexten;
				my $mismatch_penalty = ' --mp '.$mis_penalty;
				my $mapping_file = $domino_files{$reads}{'reads'}[0];
				my $botwie_system = $bowtie_path."bowtie2";
				if ($bowtie_local) { $botwie_system .= " --local"; }
				my $sam_name = $dir."/".$reference_tag."-taxa_".$reads.".sam";
				
				print LOG "+ Aligning reads for $reads against $reference_identifier as reference\n";
				if ($input_type eq 'pair_end') {
					my $second_Read_file = $domino_files{$reads}{'reads'}[1];
					print LOG "+ Reads 1: $mapping_file\n+ Reads 2: $second_Read_file\n";
					$botwie_system .= " -x ".$reference_tag." -q -1 $mapping_file -2 $second_Read_file -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
				} elsif ($input_type eq 'single_end') { ## Illumin single end, 454
					print LOG "+ Reads 1: $mapping_file...\n";
					if ($map_contig_files) { ## Mapping contigs
						$botwie_system .= " -x ".$reference_tag." -f -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
					} else { 
						$botwie_system .= " -x ".$reference_tag." -q -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
				}}
				print LOG "+ SAM file: $sam_name\n+ Mapping now...\n\n"; &debugger_print("BOWTIE2 command: ".$botwie_system); 
				
				### Map Reads
				my $error_bowtie = $dir."/reads_".$reads."-reference_".$reference_identifier."_mapping_logfile.txt";
				print LOG "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nMapping Statistics:\n\n"; 
				$botwie_system .= " 2> ".$error_bowtie;
				my $system_bowtie_call = system ($botwie_system);
				if ($system_bowtie_call != 0) {
					&printError("Exiting the script. Some error happened when calling bowtie for mapping the file $mapping_file...\n"); DOMINO::dieNicely();
				} 			
				push (@{$domino_files_split_mapping{$reads}{"SAM::Ref:".$reference_identifier}}, $sam_name);
				open (IN, $error_bowtie); while (<IN>) {print LOG $_; } close(IN);
				print LOG "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
				print LOG "\n\n+ Mapping finished for taxa $reads against $reference_identifier\n";		
				
				#################################
				## Get the the reference fasta ##
				#################################
				print LOG "\n+ Checking mapping reads in ".$sam_name."...\n";
				
				## Generate a sam for each contig
				my @temp = split ("\.sam", $sam_name);
				chdir $dir; my $dir_tmp = $temp[0]."_SPLIT"; mkdir $dir_tmp, 0755; chdir $dir_tmp;
				print LOG "\nChange dir to: ".$dir_tmp."\n";
				#print Dumper $reference_hash_fasta_ref;
				my @number_contigs = sort (keys %$reference_hash_fasta_ref);
				my $scalar = scalar @number_contigs;
				print LOG "+ This SAM file contains $scalar referense sequences...\n";
				system("ln -s $sam_name"); my @temp_name = split ("/", $sam_name);
				my $sam_base_name = $temp_name[-1];
				my $sorted_bam_file = &generate_bam($sam_base_name);
				&generate_index_bam($sorted_bam_file);
	
				if ($scalar == 1) {
					push (@{ $domino_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $temp_name[-1]);
				} else {
					print LOG "+ Using parallel threads ($split_CPU CPUs)...\n";
					print LOG "+ Splitting SAM file into several parts to speed the computation...\n"; 
					my $parts = int($scalar/$split_CPU); ## Split SAM into as many CPUs provided
					my @commands; my $iteration = 0; 
					while (1) {
						if (@number_contigs) {
							my @array1 = splice(@number_contigs, 0, $parts + 1);
							my $string = join(" ", @array1);
							my @temp_1 = split ("\.sorted.bam", $sorted_bam_file);
							my @temp_2 = split ("/", $temp_1[0]);
							my $sam_file_part = $dir_tmp."/".$temp_2[-1]."_part-".$iteration.".sam";
							#print "Contigs: ".$string."\n";
							my $command = $samtools_path." view -@ $split_CPU -Sh -o $sam_file_part $sorted_bam_file $string";
							push (@{ $domino_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $sam_file_part);
							push (@commands, $command); $iteration++; 
					} else { last; }}
					
					my $pm_SAM_split =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
					my $total_commands = scalar @commands;
					my $count_commands=0;
					for (my $a=0; $a < scalar @commands; $a++) {
						$count_commands++;
						print LOG "\t- Splitting SAM $sam_base_name: [$count_commands/$total_commands]\n";
						my $pid = $pm_SAM_split->start() and next; 
						&debugger_print("SAMTOOLS command: $commands[$a]");	
						my $system_call = system ($commands[$a]);
						$pm_SAM_split->finish(); # pass an exit code to finish
					}
					$pm_SAM_split->wait_all_children;
					&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files_split_mapping); 		
				}
				
				### Remove multimapping reads	### 
				##  DOMINO checks the SAM files generated and discards bad reads mapping, unmapping reads or multimapping reads. ##
				my @array_files_split = @{ $domino_files_split_mapping{$reads}{"SAM_Parts::Ref:".$reference_identifier}};				
				print LOG "\n+ Cleaning reads now in parallel threads ($split_CPU CPUs)...\n";
				## Get files for later dump
				for (my $i=0; $i < scalar @array_files_split; $i++) {
					my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part_".$i.".txt";
					#print "File to dump: ".$file2dump."\n";
					push (@{ $domino_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier} }, $file2dump);
				}		
				
				my $pm_SAM_parts =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				my $scalar_array_files_split = scalar @array_files_split; my $count_split=0;
				for (my $i=0; $i < scalar @array_files_split; $i++) {
					$count_split++;
					print LOG "\t- Checking each splitted file for $sam_base_name: [$count_split/$scalar_array_files_split]\n";	
					my $pid = $pm_SAM_parts->start() and next;
					my %domino_files_SAM_parts; my $discard_reads = 0; my $good_reads = 0; my $total_reads = 0;
					open (SAM, "<$array_files_split[$i]");
					my @temp = split ("\.sam", $array_files_split[$i]);
					my $output_sam = $temp[0]."_clean.sam"; 
					push (@{ $domino_files_SAM_parts{$reads}{"CLEAN_SAM_Parts::Ref:".$reference_identifier} }, $output_sam);
					open (SAM_OUT, ">$output_sam");	
					while (<SAM>) {
						chomp; my $line = $_;
						if ($line =~ /^\@.*/) { print SAM_OUT $line."\n"; next; }		
						my @sam = split ("\t", $line); 
						$total_reads++;				
						if ($line =~ /.*SA:Z:.*/) {   ## Discard Multimapping reads
							$discard_reads++; 
						} elsif ($sam[2] eq '*') { ## Discard unmapped reads
							## When calling Bowtie2 we use --no-unal not to show unaligned reads
							## but still we would check it, also for user input BAM files			
							$discard_reads++;
						} else { ## Check if read is efficiently mapping 
							## Even Bowtie2 is mapping end to end and allowing no clipping
							## we would check if any hard or soft clipping is introduced
							## (inherited from a previous version using BWA)
							
							my $cigar=$sam[5];
							my $xcent = 0; my $xcentmax = 0; my $pos = 0; my $NOpos = 0;
							my $cigar_pct_convert = $cigar_pct/100;
							
							#&debugger_print("\nLINE: ".$line."\nCIGAR: ".$cigar."\nMAX_CIGAR: ".$cigar_pct_convert."\n");
							while ($cigar !~ /^$/){         		
								if ($cigar =~ /^([0-9]+[MIDSH])/){
									my $cigar_part = $1;								
									#&debugger_print("CIGAR Part: $cigar_part\n");								
									if ($cigar_part =~ /(\d+)M/){ $pos += $1;
									} elsif ($cigar_part =~ /(\d+)I/){ $NOpos += $1;
									} elsif ($cigar_part =~ /(\d+)D/){ $NOpos += $1;
									} elsif ($cigar_part =~ /(\d+)S/){ $NOpos+= $1;
									} elsif ($cigar_part =~ /(\d+)H/){ $NOpos+= $1;
									} else { # die "Unexpected cigar: $cigar\n";
									}# close if
									$cigar =~ s/$cigar_part//;
							}}
							my $total = $pos + $NOpos; 				 ## total positions
							$xcentmax = $total * $cigar_pct_convert; ## Max bad CIGAR
							$xcent = $NOpos; ## Bad CIGAR
							if($xcent <= $xcentmax){
								print SAM_OUT $line."\n"; $good_reads++;
							} else {
								#&debugger_print("XCENT > XCENTMAX\nDISCARD READ!!\n"); 
								if ($bowtie_local) { print SAM_OUT $line."\n"; next; } 
								$discard_reads++;
					}}} close(SAM); close(SAM_OUT);	
					&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files_SAM_parts); 
				
					###################################################################
					## Generate sorted bam files in order to be able to get coverage ##
					###################################################################
					my $clean_sorted_bam = &generate_bam($array_files_split[$i]);
					push (@{ $domino_files_SAM_parts{$reads}{"clean_BAM_Parts::Ref:".$reference_identifier} }, $clean_sorted_bam);
	
					%discard_contigs = (); ## Initialize some hashes
					#print "\t- Generating coverage statistics for $clean_sorted_bam\n";
				
					## Generate Coverage statistics for the alignment file
					my @tmp_bam_name = split ("\.sorted.bam", $clean_sorted_bam);
					my $coverage_file = $tmp_bam_name[0]."_coverage_stats.txt";
					push (@{ $domino_files_SAM_parts{$reads}{"coverage_Parts::Ref:".$reference_identifier} }, $coverage_file);
	
					my $coverage_samtools_command = $samtools_path." depth ".$clean_sorted_bam." > ".$coverage_file;
					my $system_coverage_call = system ($coverage_samtools_command); 
					&debugger_print("SAMTOOLS command: $coverage_samtools_command");
					if ($system_coverage_call != 0) {
						&printError("Exiting the script. Some error happened when calling SAMtools for obtaining coverage of file $sorted_bam[$i]...\n"); DOMINO::dieNicely();
					}
					
					unless (-z $coverage_file) {
						#print "\t- Filtering Coverage Stats...\n\n";
						my ($sum_coverage_each, $total_positions_each);
						my %max_cov;
						open (COVERAGE, "<$coverage_file");
						while (<COVERAGE>) {
							my $line = $_;
							chomp $line;
							my @array = split (/\s+/,$line);
							if ($array[2]) {
								my $coverage_position = $array[2];
								$sum_coverage_each += $coverage_position;
								$total_positions_each++;
								my $contig = $array[0];
								if (defined($max_cov{$contig})) {
									if ($coverage_position > $max_cov{$contig}) { $max_cov{$contig} = $coverage_position; }
								} else { $max_cov{$contig} = $coverage_position; }
						}}
						close(COVERAGE);
						
						## Print max coverage to a file
						my $max_cov_file = $coverage_file."_max_coverage";
						push (@{$domino_files_SAM_parts{$reads}{"max_coverage_Parts::Ref:".$reference_identifier} }, $max_cov_file);
						open (OUT_COV, ">$max_cov_file");
						foreach my $contig (sort keys %max_cov) {
							print OUT_COV $contig.":".$max_cov{$contig}."\n";
						} close (OUT_COV);
						undef %max_cov;
						
						## Push useful info
						my $coverage_each = $sum_coverage_each.":".$total_positions_each;
						push (@{ $domino_files_SAM_parts{$reads}{"mean_coverage_Parts::Ref:".$reference_identifier} }, $coverage_each);
						push (@{$domino_files_SAM_parts{$reads}{"discard_reads::Ref:".$reference_identifier} }, $discard_reads);
						push (@{$domino_files_SAM_parts{$reads}{"good_reads::Ref:".$reference_identifier} }, $good_reads);
						push (@{$domino_files_SAM_parts{$reads}{"total_reads::Ref:".$reference_identifier} }, $total_reads);
		
						# Dump info into file
						DOMINO::printDump(\%domino_files_SAM_parts, $domino_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier}[$i]);
					}
					$pm_SAM_parts->finish(); # pass an exit code to finish
				}
				$pm_SAM_parts->wait_all_children;
				my @dump_files1 = @{ $domino_files_split_mapping{$reads}{"DUMP_Parts::Ref:".$reference_identifier} };
				my $hash_dump1 = &retrieve_info(\@dump_files1, \%domino_files_split_mapping);
				%domino_files_split_mapping = %{$hash_dump1};

				## Get total reads: discard, good and total
				my $discard_reads = 0; my $good_reads = 0; my $total_reads = 0; my @array;
				push (@array, $domino_files_split_mapping{$reads}{"discard_reads::Ref:".$reference_identifier});
				push (@array, $domino_files_split_mapping{$reads}{"good_reads::Ref:".$reference_identifier});
				push (@array, $domino_files_split_mapping{$reads}{"total_reads::Ref:".$reference_identifier});
				for (my $i=0; $i < scalar @array; $i++) {
					my @array_2 = @{$array[$i]};
					for (my $j=0; $j < scalar @array_2; $j++) {
						if ($i==0) { $discard_reads += $array_2[$j];
						} elsif ($i==1) { $good_reads += $array_2[$j];
						} else { $total_reads += $array_2[$j];
				}}}		
				
				##############################################################################################	
				## Calculate the probability of being under a poisson distribution with the mean of our data
				##############################################################################################
				my @stats = @{ $domino_files_split_mapping{$reads}{"mean_coverage_Parts::Ref:".$reference_identifier}};
				my ($sum_coverage, $total_positions);
				for (my $i=0; $i < scalar @stats; $i++) {
					my @stats_each = split(":", $stats[$i]);
					$sum_coverage += $stats_each[0];
					$total_positions += $stats_each[1];
				}
				my $mean_coverage = $sum_coverage/$total_positions;
				my $mean = sprintf ("%.3f", $mean_coverage);
	
				my @max_cov_files = @{ $domino_files_split_mapping{$reads}{"max_coverage_Parts::Ref:".$reference_identifier} };
				my $num_contigs; my %discard_contigs; my $contigs_discarded;
				for (my $h=0; $h < scalar @max_cov_files; $h++) {
					## Read Coverage file
					open (COV_READ, $max_cov_files[$h]);
					while (<COV_READ>) {
						chomp;
						my $line = $_;
						my @array = split(":", $line); 
						$num_contigs++; # $array[0] contig id in case we need					
						my $prob_poisson;
						if ($array[1] > 169) { ## Factorial would be out of range!
							$prob_poisson = 0;	
						} else { $prob_poisson = &Poisson_distribution($array[1], $mean_coverage); }		
						if ($prob_poisson < $level_significance_coverage_distribution) { # Discard
							$contigs_discarded++; $discard_contigs{$array[0]}++;
						} elsif ($prob_poisson eq 'nan') { # Discard
							$contigs_discarded++; $discard_contigs{$array[0]}++;
						} elsif ($array[1] == 0) {  # Discard
							$contigs_discarded++; $discard_contigs{$array[0]}++;
				}}}
	
				my $percentage_discard_Reads; my $h; my $perc_contig; my $h_cont;
				if ($discard_reads) { $percentage_discard_Reads = (1 - ($discard_reads/$total_reads))*100; $h = sprintf ("%.3f", $percentage_discard_Reads);} else { $h=100; $discard_reads=0;}
				if ($contigs_discarded) { $perc_contig = (1 - ($contigs_discarded/$num_contigs))*100; $h_cont = sprintf ("%.3f", $perc_contig); } else { $contigs_discarded = 0; $h_cont = 100; }
				
				## Adjust SAM files
				print LOG "\n+ Adjusting the SAM/BAM files using parallel threads ($split_CPU CPUs)\n+ Splitted files would be used...\n";
				my @parts_clean_sam = @{ $domino_files_split_mapping{$reads}{"CLEAN_SAM_Parts::Ref:".$reference_identifier}};
	
				## Get files for later dump
				for (my $i=0; $i < scalar @parts_clean_sam; $i++) {
					my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part2_".$i.".txt";
					#print "File to dump: ".$file2dump."\n";
					push (@{ $domino_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} }, $file2dump);
				}	
				
				my $pm_SAM_PILEUP =  new Parallel::ForkManager($split_CPU); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				my $total_parts_clean_sam = scalar @parts_clean_sam;
				my $count_totalParts=0;
				for (my $j=0; $j < scalar @parts_clean_sam; $j++) {
					$count_totalParts++;
					print LOG "\t- Adjusting SAM/BAM part file for $sam_base_name: [$count_totalParts/$total_parts_clean_sam]\n";
					my @basename = split("/", $parts_clean_sam[$j]);
					my @name = split(".sam", $basename[-1]);
					my $pid = $pm_SAM_PILEUP->start($name[0]) and next;
					my %domino_files_SAM_PILEUP;
					my @tmp_sam = split("\.clean.sam", $parts_clean_sam[$j]);
					my $sam_filter = $tmp_sam[0]."_filtered.sam";
					open (SAM_OUT, ">$sam_filter"); open (SAM, "<$parts_clean_sam[$j]");
					while (<SAM>) {
						chomp; my $line = $_;
						if ($line =~ /^@.*/ ) { print SAM_OUT $line."\n";	
						} else {
							my @array = split (/\s+/,$line);
							if (!$discard_contigs{$array[2]}) { print SAM_OUT $line."\n";	
					}}}
					close(SAM_OUT); close(SAM); undef %discard_contigs;
					#print LOG "- File checked: Contigs and Reads discarded...\n";
					push (@{ $domino_files_SAM_PILEUP{$reads}{"FILTERED_SAM_Parts::Ref:".$reference_identifier} }, $sam_filter);
					my $bam_filtered_returned = &generate_bam($sam_filter);
					push (@{ $domino_files_SAM_PILEUP{$reads}{"FILTERED_BAM_Parts::Ref:".$reference_identifier} }, $bam_filtered_returned);
					unless ($reads eq $reference_identifier) { ## DO NOT GENERATE FILTER PROFILE FOR REFERENCE
						print LOG "- Generate a PILEUP file for $sam_filter...\n";
						my $dir_returned = &generate_filter_PILEUP($bam_filtered_returned, $domino_files{$reference_identifier}{'contigs'}[0], $reference_hash_fasta_ref, $reference_identifier, $reads);
						print LOG "Finish PILEUP for $bam_filtered_returned\n";
						push (@{ $domino_files_SAM_PILEUP{$reads}{"PROFILE::Ref:".$reference_identifier} }, $dir_returned);
					}
					# Dump info into file
					DOMINO::printDump(\%domino_files_SAM_PILEUP, $domino_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier}[$j]);
					$pm_SAM_PILEUP->finish($name[0]); # pass an exit code to finish
				}
				$pm_SAM_PILEUP->wait_all_children;
				
				my @dump_files = @{ $domino_files_split_mapping{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} };
				my $hash_dump = &retrieve_info(\@dump_files, \%domino_files_split_mapping);
				%domino_files_split_mapping = %{$hash_dump};

				print LOG "\n\n\n";			
				my $stats_file = $dir."/mapping_statistics_Ref_".$reference_identifier."_Reads_".$reads.".txt";
				push (@{ $domino_files_split_mapping{$reference_identifier}{'stats'} }, $stats_file);
				print LOG "================ Filtering Statistics ======================\n"; 
				open (STATS, ">$stats_file");			
				print STATS "File: $sam_name\n"; 																print LOG "File: $sam_name\n";
				print STATS "Reference: ".$reference_identifier."\n"; 											print LOG "Reference: ".$reference_identifier."\n";
				print STATS "Reads: ".$reads."\n";																print LOG "Reads: ".$reads."\n";
				print STATS "==========================================================\n";						print LOG "==========================================================\n";
				print STATS "Total Contigs: $num_contigs\n";													print LOG "Total Contigs: $num_contigs\n";
				print STATS "Total Reads mapping: $total_reads\n";												print LOG "Total Reads mapping: $total_reads\n";
				print STATS "Coverage Mean: x".$mean."\n";														print LOG "Coverage Mean: x".$mean."\n";
				print STATS "\n******* CONTIGS *******\n";														print LOG "\n******* CONTIGS *******\n";
				print STATS "Contigs discarded: ".$contigs_discarded."\n";										print LOG "Contigs discarded: ".$contigs_discarded."\n";
				print STATS "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n"; 	print LOG "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n";
				print STATS "\n******* READS *******\n";														print LOG "\n******* READS *******\n";
				print STATS "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";	print LOG "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";
				print STATS "Reads remaining: $good_reads\t( ".$h." %)\n";										print LOG "Reads remaining: $good_reads\t( ".$h." %)\n";
				print LOG "============================================================\n"; 
				print LOG "\n\n\n";			

				DOMINO::printDump(\%domino_files_split_mapping, $domino_files{$reads}{"DUMP_Mapping::Ref:".$reference_identifier}[0]);
				$pm_read_Reference->finish(); # pass an exit code to finish
			} ## foreach reads
			$pm_read_Reference->wait_all_children;			
			
			## Get files for later dump
			foreach my $reads_here (sort keys %domino_files) {
				unless ($domino_files{$reads_here}{'reads'}) { next; }
				my @dump_files = @{ $domino_files{$reads_here}{"DUMP_Mapping::Ref:".$reference_identifier} };
				&retrieve_info(\@dump_files, \%domino_files);
			}
			&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
			print "\n\n"; DOMINO::printHeader("", "+");DOMINO::printHeader(" Mapping finished for Reference $reference_identifier ", "+");DOMINO::printHeader("", "+"); &time_log(); print "\n";		
		} # foreach reference
	} elsif ($option eq "msa_alignment") {	
		mkdir $msa_dirname, 0755; chdir $align_dirname;
		push (@{ $domino_files{'MSA_files'}{'folder'} }, $msa_dirname); 
	
		#####################################
		### Check MSA: file/folder/RADseq ###
		#####################################
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
		if ($radseq_like_data) {
			my @file_name = split("/", $domino_files{'RADseq'}{'file'}[0]);
			my $file_path = $align_dirname."/".$file_name[-1];
			system("ln -s $domino_files{'RADseq'}{'file'}[0] $file_path");
			&debugger_print("RADseq data file provided...");
			print "\n\n"; DOMINO::printHeader(" Checking RADseq file provided ", "%");
			
			## Split files
			print "+ Splitting file into multiple files to speed the computation...\n";
			my $size = DOMINO::get_size($file_path);
			my $chars = int($size/$num_proc_user);
				&debugger_print("Total Size: $size\nCharacters to split: $chars");
				&debugger_print("File: $file_path");		

			my $files_ref;
			if ($pyRAD_file) { $files_ref = DOMINO::loci_file_splitter($file_path, $chars,'loci');
			} else { $files_ref = DOMINO::fasta_file_splitter($file_path, $chars, 'fa'); }		
			push (@{ $domino_files{'RADseq'}{'parts'}}, @$files_ref);		
			&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
	
			for (my $i=0; $i < scalar @$files_ref; $i++) {
				my $file2dump = $align_dirname."/dump_file_split_Part_".$i.".txt";
				push (@{ $domino_files{'all'}{'dump_file_split'} }, $file2dump);
			}
	
			## Implement threads
			print "\n+ For each loci a MSA file would be generated...\n+ Parsing splitted files...\n+ Using parallel threads ($num_proc_user CPUs)...\n";
			if ($pyRAD_file) {		
				&debugger_print("pyRAD file provided...");
				### RADSEQ like data ####
				## Parse pyRAD loci file provided
				my $pm_pyRAD =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				$pm_pyRAD->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code\n"; } );
				$pm_pyRAD->run_on_start( sub { my ($pid,$ident)=@_; print "\t- pyRAD analysis for file $ident and PID=$pid started\n"; } );
				for (my $i=0; $i < scalar @$files_ref; $i++) {
					my @basename = split("/", $$files_ref[$i]);
					my @name = split(".loci", $basename[-1]);
			
					my $pid = $pm_pyRAD->start($name[-1]) and next; 
					my %domino_files_pyRAD_split;
					my $counter = 1; my %hash;
					open(FILE, $$files_ref[$i]) || die "Could not open the $$files_ref[$i] ...\n";
					while (<FILE>) {		
						next if /^#/ || /^\s*$/;
						my $line = $_;
						chomp $line;
						&debugger_print($line); &debugger_print("Counter: $counter");
						if ($line =~ /\/\//) { 
							my $file = $msa_dirname."/".$name[0]."_loci_".$counter.".fasta";
							foreach my $keys (sort keys %hash) { 
								if ($domino_files{$keys}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
									open (OUT, ">>$file"); print OUT ">".$keys."\n".$hash{$keys}."\n"; close(OUT);
							}} $counter++; undef %hash; next;
						}
						$line =~ s/\s+/\t/g; 
						$line =~ s/\%/>/g; ## Sometimes there is this symbol or at least in my test set
						my @array = split("\t", $line);
						$array[0] =~ s/\>//;
						$hash{$array[0]} = $array[1];
						if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
							$domino_files_pyRAD_split{$array[0]}{'taxa'}[0]++;
						}
					} close (FILE); undef %hash;
					#print Dumper \%domino_files_pyRAD_split;
					DOMINO::printDump(\%domino_files_pyRAD_split, $domino_files{'all'}{'dump_file_split'}[$i]);
					$pm_pyRAD->finish($name[-1]); # pass an exit code to finish
				}
				$pm_pyRAD->wait_all_children; 
				print "\n\n";
				print "***************************************************\n";
				print "**** All pyRAD parsing processes have finished ****\n";
				print "***************************************************\n\n";		
			} elsif ($stacks_file) {
				## Parse STACKS file provided
				&debugger_print("STACKS file provided...");
				my $pm_STACKS =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				$pm_STACKS->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code \n"; } );
				$pm_STACKS->run_on_start( sub { my ($pid,$ident)=@_; print "\t- STACKS analysis for file $ident and PID=$pid started\n"; } );
				for (my $i=0; $i < scalar @$files_ref; $i++) {
					my @basename = split("/", $$files_ref[$i]);
					my @name = split(".fa", $basename[-1]);
					my $pid = $pm_STACKS->start($name[-1]) and next; 
					
					my %domino_files_STACKS_split;
					my (%hash, $new_id); my $first = 0; my $previous;
					$/ = ">"; ## Telling perl where a new line starts
					open(FILE, $$files_ref[$i]) || die "Could not open the $$files_ref[$i] ...\n";
					while (<FILE>) {		
						next if /^#/ || /^\s*$/; chomp;
						my ($titleline, $sequence) = split(/\n/,$_,2);
						next unless ($sequence && $titleline);
						chomp $sequence;
						$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
						$titleline =~ s/\r//g;
						&debugger_print($titleline."\t".$sequence);
						if ($titleline =~ /(CLocus\_\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
							my $CLocus = $1; my $sample = $2;
							if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') { $domino_files_STACKS_split{$sample}{'taxa'}[0]++; }
							if ($first == 0) { 
								push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence); $first++;  
								$previous = $CLocus; next;
							}							
							if ($hash{$CLocus}) {
								push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence); $previous = $CLocus;
							} else {
								#&debugger_print("Ref", \%hash);
								my %Clocus_hash = %{ $hash{$previous} };
								#print "CLocus hash:\n"; print Dumper \%Clocus_hash;								
								## Parse
								my $out_file = $msa_dirname."/".$previous.".fasta";
								open (OUT, ">>$out_file");
								foreach my $keys (sort keys %Clocus_hash) {
									if ($domino_files{$keys}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
										my @array = @{$Clocus_hash{$keys}};
										for (my $i=0; $i < scalar @array; $i++) {
											my @split = split(":::", $array[$i]);
											print OUT ">".$split[0]."\n".$split[1]."\n";
								}}} close (OUT);
								push (@{ $hash{$CLocus}{$sample} }, $titleline.":::".$sequence);
								$previous = $CLocus;
					}}} close(FILE); $/ = "\n";
					
					## Parse last
					my %Clocus_hash = %{ $hash{$previous} };
					my $out_file = $msa_dirname."/".$previous.".fasta";
					open (OUT, ">>$out_file");
					foreach my $keys (sort keys %Clocus_hash) {
						if ($domino_files{$keys}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
							my @array = @{$Clocus_hash{$keys}};
							for (my $i=0; $i < scalar @array; $i++) {
								my @split = split(":::", $array[$i]);
								print OUT ">".$split[0]."\n".$split[1]."\n";
					}}} close (OUT);					
					DOMINO::printDump(\%domino_files_STACKS_split, $domino_files{'all'}{'dump_file_split'}[$i]);
					$pm_STACKS->finish($name[-1]); # pass an exit code to finish
				}
				$pm_STACKS->wait_all_children; 
				print "\n\n";
				print "****************************************************\n";
				print "**** All STACKS parsing processes have finished ****\n";
				print "****************************************************\n\n";		
			} 
			print "\n"; &time_log(); print "\n";
		} else {
			### MSA file or folder provided
			#mkdir $msa_dirname, 0755;
			if ($msa_file) {
				## Check the alignment format provided...
				print "\n\n"; DOMINO::printHeader(" Checking Alignment file provided ", "%");
				&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
				my $species_alignment_ref = &read_phylip_aln($domino_files{'MSA'}{'file'}[0]);
				if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
					undef $domino_files{'taxa'}{'user_Taxa'};		
					for (my $i=0; $i < scalar @$species_alignment_ref; $i++) {
						push (@{ $domino_files{$$species_alignment_ref[$i]}{'taxa'}}, 1);
						push (@{ $domino_files{'taxa'}{'user_Taxa'} }, $$species_alignment_ref[$i]);
				}}
				print "\n"; &time_log(); print "\n";		
			} elsif ($msa_fasta_folder) {				
				chdir $msa_dirname;
				print "\n\n"; DOMINO::printHeader(" Checking Alignment folder provided ", "%");
				my @array_files = @{ $domino_files{'MSA_folder'}{'files'} };
				my $tmp = $align_dirname."/tmp"; mkdir $tmp, 0755;
				for (my $i=0; $i < scalar @array_files; $i++) {
					my $file2dump = $tmp."/dump_file_split_Part_".$i.".txt";
					#print "File to dump: ".$file2dump."\n";
					push (@{ $domino_files{'all'}{'dump_file_split'} }, $file2dump);
				}
				&debugger_print("MSA folder provided...");
				print "+ Using parallel threads ($num_proc_user CPUs)...\n";			
				my $pm_MSA_folder =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				$pm_MSA_folder->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\t- Child process finished for file $ident; PID=$pid & ExitCode=$exit_code \n"; } );
				$pm_MSA_folder->run_on_start( sub { my ($pid,$ident)=@_; print "\t- Alignment analysis for file $ident and PID=$pid started \n"; } );
				for (my $i = 0; $i < scalar @array_files; $i++) {
					my @species_alignment;
					my $pid = $pm_MSA_folder->start($array_files[$i]) and next; 
					my $file_path = $domino_files{'MSA_folder'}{'folder'}[0]."/".$array_files[$i];
					my %alignment;			
					my %domino_files_MSA_folder;
					if ($array_files[$i] =~ /(.*)\.fasta/) {
						my $name = $1;
						open(FILE, $file_path) || die "Could not open the $file_path...\n";
						$/ = ">"; ## Telling perl where a new line starts
						while (<FILE>) {		
							next if /^#/ || /^\s*$/;
							chomp;
							my ($titleline, $sequence) = split(/\n/,$_,2);
							next unless ($sequence && $titleline);
							chomp $sequence; chomp $titleline;
							$titleline =~ s/\r//g;						
							if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
								unless (grep /$titleline/, @{ $domino_files_MSA_folder{'taxa'}{'user_Taxa'} }) { 
									$domino_files_MSA_folder{$titleline}{'taxa'}[0]++;
									push (@{ $domino_files_MSA_folder{'taxa'}{'user_Taxa'} }, $titleline);
								}
								$alignment{$titleline} = $sequence;
							} elsif ($domino_files{$titleline}{'taxa'}) { $alignment{$titleline} = $sequence; }
						} close(FILE); $/ = "\n";
						
						if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') { 
							system("ln -s $file_path");
						} else {
							my $out_file = $msa_dirname."/parsed_".$name.".fasta";
							open (OUT, ">$out_file");
							foreach my $seqs (sort keys %alignment) {
								print OUT ">".$seqs."\n".$alignment{$seqs}."\n";
							} close (OUT);
					}}

					DOMINO::printDump(\%domino_files_MSA_folder, $domino_files{'all'}{'dump_file_split'}[$i]);
					$pm_MSA_folder->finish($i); # pass an exit code to finish
				}
				$pm_MSA_folder->wait_all_children; 
				print "\n\n";
				print "*********************************************\n";
				print "**** All parsing processes have finished ****\n";
				print "*********************************************\n\n";
				my $files = scalar @array_files;
				print "\n\n+ Parsing of the $files files has been done...\n"; print "\n"; &time_log(); print "\n";
		}}
		if ($domino_files{'all'}{'dump_file_split'}) {
			my @dump_files = @{ $domino_files{'all'}{'dump_file_split'} };
			for (my $j=0; $j < scalar @dump_files; $j++) {
				if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
				open (DUMP_IN, "$dump_files[$j]");
				while (<DUMP_IN>) {
					my $line = $_; chomp $line;
					my @array = split("\t", $line);
					push (@{ $domino_files{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_IN); }}
			
		unless ($MID_taxa_names) {
			print "\n\n"; DOMINO::printHeader("","#"); print "NOTE:\n\n";
			print "\t+ No taxa names were provided so DOMINO have parsed and checked for the names...\n";	
			print "\t+ DOMINO would use all the taxa available...\n\n";
			my @taxa_sp;
			foreach my $keys (sort keys %domino_files) { 
				if ($domino_files{$keys}{'taxa'}) {
					DOMINO::printDetails("\tName: $keys\n", $mapping_parameters, $param_Detail_file_markers);
					$number_sp++;
					push (@taxa_sp, $keys);
			}}
			print "\n\n"; DOMINO::printHeader("","#"); print "\n\n";
			if (!$minimum_number_taxa_covered) { 
				$minimum_number_taxa_covered = 0;  ## Force to be any taxa
			}
			push (@{ $domino_files{'number_taxa'}{'number'}}, $number_sp);
			push (@{ $domino_files{'number_taxa'}{'MCT'}}, $minimum_number_taxa_covered);
			$MID_taxa_names = join(",", @taxa_sp);
			push (@{ $domino_files{'taxa_string'}{'string'}[0] }, $MID_taxa_names);
		}
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
	}
	
	## Dump info up to now to a file
	my $file2dump = $align_dirname."/".$datestring."_DUMP.txt";
	DOMINO::printDump(\%domino_files, $file2dump);	

	## Dump parameters to a file
	my $file2dump_param = $align_dirname."/".$datestring."_DUMP_param.txt";
	DOMINO::printDump(\%domino_params, $file2dump_param);	
	## Move parameters and error file to folder
	File::Copy::move($mapping_parameters, $align_dirname."/");
} else {
	print "+ Files would be obtained...\n\n"; ## No_Profile_Generation|NPG: just get files
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
### MSA alignment
if ($option eq "msa_alignment") {

	##########################################################
	############# MARKER SELECTION ###########################
	##########################################################
	my $profile_dir = $marker_dirname."/PROFILES";	mkdir $profile_dir, 0755;
	my $msa_dir_tmp = $marker_dirname."/MSA_markers_tmp"; mkdir $msa_dir_tmp, 0755;
	my $dir_Dump_file = $marker_dirname."/DUMP_files"; mkdir $dir_Dump_file, 0755;
	my $msa_dir = $marker_dirname."/MSA_markers"; mkdir $msa_dir, 0755;

	my $array_files_fasta_msa_ref = DOMINO::readDir($domino_files{'MSA_files'}{'folder'}[0]);
	my @array_files_fasta_msa = @$array_files_fasta_msa_ref;
	my $total_files = scalar @array_files_fasta_msa;
	print "+ Checking files in folder generated...\n";
	my $counter = 0;	
	my $pm_MARKER_MSA_files =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	$pm_MARKER_MSA_files->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
	$pm_MARKER_MSA_files->run_on_start( sub { my ($pid,$ident)=@_; } );
	for (my $i=0; $i < scalar @array_files_fasta_msa; $i++) {
		if ($array_files_fasta_msa[$i] eq "." || $array_files_fasta_msa[$i] eq ".." || $array_files_fasta_msa[$i] eq ".DS_Store"  ) {next;}
		$counter++; 
		if ($total_files > 100) {
			my $perc = sprintf( "%.3f", ( $counter/$total_files )*100 );
			print "\t- Checking each sequence: [ $perc % ]...\r";
		} else { print "\t- Checking each sequence: [$counter/$total_files]...\r";}	
	
		my $pid = $pm_MARKER_MSA_files->start($i) and next;
		my %domino_files_msa;
		my $file_path = $domino_files{'MSA_files'}{'folder'}[0]."/".$array_files_fasta_msa[$i];
		unless (-f -e -r -s $file_path) {
			&printError("File $file_path is not readable or empty. It would be skipped....\nPlease discarded from the folder...\n");
		} else {
			#print "\nChecking now: $array_files_fasta_msa[$i]\n";
			my ($region_id, $string2print_ref, $hash_ref_msa);
			if ($stacks_file) {
				#### STACKS
				$hash_ref_msa = DOMINO::readFASTA_hash($file_path);
				my $filename = $array_files_fasta_msa[$i];	
				my %hash = %$hash_ref_msa; my ( %hash2fill, %hash2return);
				foreach my $keys (sort keys %hash) {
				if ($keys =~ /(CLocus\_\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
					$region_id = $1; my $sample = $2;
					my $sequence = $hash{$keys};
					push (@{ $hash2fill{$sample} }, $keys.":::".$sequence);
				}}
				foreach my $keys (sort keys %hash2fill) {
					my $alleles = scalar @{ $hash2fill{$keys} };
					if ($alleles == 1) {
						my @array = split(":::", $hash2fill{$keys}[0]);
						$hash2return{$keys} = $array[1];					
					} elsif ($alleles == 2){					
						my @array = split(":::", $hash2fill{$keys}[0]);
						my @allele1 = split("", $array[1]);
						my @array2 = split(":::", $hash2fill{$keys}[1]);
						my @allele2 = split("", $array2[1]);
						my @string;
						for (my $i=0; $i < scalar @allele1; $i++) {
							my @tmp; push(@tmp, $allele1[$i]); push(@tmp, $allele2[$i]);
							my @tmp_uniq_sort = uniq(sort @tmp);
							if (scalar @tmp_uniq_sort == 1) {
								if ($tmp_uniq_sort[0] eq 'N') {  		push (@string, 'N');
								} elsif ($tmp_uniq_sort[0] eq '-') { 	push (@string, '-');
								} else { 								push (@string, $allele1[$i]);
							}} else {
								my $hash;
								$$hash{$allele1[$i]}++; $$hash{$allele2[$i]}++;
								push (@string, &get_amb_code($hash));
						}}
						my $tmp = join ("", @string);
						$hash2return{$keys} = $tmp;
						#print $array[1]."\n"; print $array2[1]."\n"; print $tmp."\n\n";
					} elsif ($alleles > 2) {
						&printError("ERROR: $region_id: [ $keys ] contains more than 2 alleles\nDOMINO would not die here but this locus would be skipped....\n");
						$pm_MARKER_MSA_files->finish;
				}}
				#print "\n\n"; print Dumper \%hash2return;

				undef %hash2fill;
				unless ($dnaSP_flag) {
					my $valueReturned = &check_marker_pairwise(\%hash2return);
					if ($valueReturned != 1) { $pm_MARKER_MSA_files->finish; }
				}
				$string2print_ref = &check_marker_ALL(\%hash2return, "Ref"); # Check there is a minimun variation
				if ($string2print_ref eq 'NO' ) { $pm_MARKER_MSA_files->finish;}
			} else {
				#### OTHER MSAs
				my @path_file = split("\.fasta", $array_files_fasta_msa[$i]);
				$region_id = $path_file[0];
				$hash_ref_msa = DOMINO::readFASTA_hash($file_path);				
				my $taxa4marker = scalar keys %{ $hash_ref_msa };
				if ($taxa4marker < $minimum_number_taxa_covered) { $pm_MARKER_MSA_files->finish; }
				unless ($dnaSP_flag) {
					my $valueReturned = &check_marker_pairwise($hash_ref_msa);
					if ($valueReturned != 1) { $pm_MARKER_MSA_files->finish; }
				}
				$string2print_ref = &check_marker_ALL($file_path); # Check there is a minimun variation
				if ($string2print_ref eq 'NO' ) { $pm_MARKER_MSA_files->finish;}
			}
			
			#print Dumper $string2print_ref;
			## taxa_included	var_sites	length	profile	effective_length
			##	0					1		2		3		4			
			
			my ($taxa, $var_sites, $length_string, $string_profile, $effective_length) = @{ $string2print_ref };
			my @array_taxa_split = split(",", $taxa);
			## Control steps
			unless ($variable_divergence) { 
				if ($var_sites > $variable_positions_user_max) {$pm_MARKER_MSA_files->finish;} 
				if ($var_sites < $variable_positions_user_min) {$pm_MARKER_MSA_files->finish;} 
			}
			
			if ($select_markers) {
				#### SELECTION MODE
				
				## Check marker
				my $variation_perc = ($var_sites/$effective_length)*100;
				my $h = sprintf ("%.3f", $variation_perc);
				my $string = $region_id."\t".$taxa."\t".$var_sites."\t".$effective_length."\t".$h;		
				push (@{ $domino_files_msa{$region_id}{'string'} }, $string);
											
				## Print profile
				my $profile_dir_file = $profile_dir."/".$region_id."_profile.txt";
				open (PRF, ">$profile_dir_file"); print PRF ">".$region_id."\n".$string_profile."\n"; close (PRF);
				#push (@{ $domino_files_msa{$region_id}{'profile'} }, $profile_dir_file);
			
				#print $file_path."\n"; print Dumper \%domino_files_msa;
				# Print format: same as input and *mmfas 				
				my $msa_fasta = $msa_dir."/".$region_id.".fasta";
				open (OUT_MSA, ">$msa_fasta");
				foreach my $keys (sort keys %{ $hash_ref_msa }) { print OUT_MSA ">".$keys."\n".$$hash_ref_msa{$keys}."\n"; }
				close (OUT_MSA);
				#push (@{ $domino_files_msa{$region_id}{'markers'} }, $msa_fasta);

				my $dump_folder_files = $dir_Dump_file."/dump_markers_".$region_id.".txt";
				# Dump into file # print Dumper \%domino_files_msa;
				DOMINO::printDump(\%domino_files_msa, $dump_folder_files);	
			} elsif ($identify_markers) { ## Identify markers in MSA alignments
				#### DISCOVERY MODE

				my $profile_dir_file = $profile_dir."/".$region_id.".txt";
				open (OUT, ">$profile_dir_file"); print OUT ">".$region_id."\n".$string_profile."\n"; close(OUT);
				push (@{ $domino_files_msa{$region_id}{'profile'} }, $profile_dir_file);
		

				my $mergeProfile = $profile_dir."/".$region_id."_merged_ARRAY.txt";
				my $string = $window_size_VARS_range;$string =~ s/\:\:/-/; my $string2 = $window_size_CONS_range; $string2 =~ s/\:\:/-/;	
				my $mergeCoord;

				if ($variable_divergence) { $mergeCoord = $profile_dir."/".$region_id."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
				} else { 					$mergeCoord = $profile_dir."/".$region_id."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
				}
				push (@{ $domino_files_msa{$region_id}{'mergeProfile'} }, $mergeProfile);
				push (@{ $domino_files_msa{$region_id}{'mergeCoord'} }, $mergeCoord);
				open (OUT_COORD, ">$mergeCoord");

				## Identify markers in MSA alignments
				my $infoReturned = &sliding_window_conserve_variable(\$region_id, \$string_profile); 
				#print $$fileReturned."\n";
				if (!$infoReturned) { $pm_MARKER_MSA_files->finish;
				} else { 
					my @array = @$infoReturned;
					for (my $j=0; $j < scalar @array; $j++) {
						print OUT_COORD $array[$j]."\n";
				}}
				close (OUT_COORD);
				
				my $file_markers_collapse = &check_overlapping_markers($mergeCoord, \$profile_dir_file);
				push (@{ $domino_files_msa{$region_id}{'markers_Merge'} }, $file_markers_collapse);
				
				# Retrieve fasta sequences...
				my $output_file = $msa_dir_tmp."/".$region_id."_markers_retrieved.txt";				
				my $array_Ref = &check_DOMINO_marker($output_file, $msa_dir_tmp, $file_markers_collapse, $file_path);
				unless (scalar @$array_Ref == 0) { 
					push (@{ $domino_files_msa{$region_id}{'markers_files'} }, @$array_Ref);
					my $dump_folder_files = $dir_Dump_file."/dump_markers_".$region_id.".txt";
					push (@{ $domino_files_msa{$region_id}{'markers'} }, $output_file);
					# Dump into file # print Dumper \%domino_files_msa;
					DOMINO::printDump(\%domino_files_msa, $dump_folder_files);	
	}}} $pm_MARKER_MSA_files->finish();
	}
	$pm_MARKER_MSA_files->wait_all_children;
	print "\n\n";
	print "**********************************************\n";
	print "**** All checking processes have finished ****\n";
	print "**********************************************\n\n";	
	print "\n"; &time_log(); print "\n"; chdir $marker_dirname;

	#################################################################################
	##	Once the coordinates are found, print different files with the information ##
	#################################################################################	
	### open Output and Error file
	my $output_file_coord = "DM_markers-summary.txt"; open (OUT_coord,">$output_file_coord")or die "Cannot write the file $output_file_coord";
	print OUT_coord "Region\t\tTaxa_included\tVariable_Positions\tEffective_length\tVariation(%)\n";
	my $out_file = "DM_markers";
	if ($pyRAD_file) { $out_file .= ".loci"; } elsif ($stacks_file) {$out_file .= ".fa";} else {$out_file .= ".fasta";}
	open (OUT, ">$out_file");
	print "+ Printing selected markers in $output_file_coord and $out_file...\n";
	
	if ($identify_markers) {
		my %hashRetrieve;
		my $array_files = DOMINO::readDir($dir_Dump_file);
		my @dump_files = @{ $array_files };
		for (my $j=0; $j < scalar @dump_files; $j++) {
			if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
			open (DUMP_IN, "$dir_Dump_file/$dump_files[$j]");
			while (<DUMP_IN>) {
				my $line = $_; chomp $line;
				my @array = split("\t", $line);
				&debugger_print($line);
				push (@{ $hashRetrieve{$array[0]}{$array[1]}}, $array[2]);
		} close (DUMP_IN); }

		foreach my $regions (sort keys %hashRetrieve) {
			if ($hashRetrieve{$regions}{'markers'}) {
				my @array_coord = @{$hashRetrieve{$regions}{'markers'}};
				for (my $i=0; $i < scalar @array_coord; $i++) {
					open (FILE, $array_coord[$i]); while (<FILE>) { print OUT_coord $_; } close(FILE);					
			}}
			if ($hashRetrieve{$regions}{'markers_files'}) {
				my @array_markers = @{$hashRetrieve{$regions}{'markers_files'}};
				for (my $j=0; $j < scalar @array_markers; $j++) {
					File::Copy::move($array_markers[$j], $msa_dir);
	}}}} else {		
		my %hashRetrieve;
		my $array_files_dump = DOMINO::readDir($dir_Dump_file);
		my @dump_files = @{ $array_files_dump };
		for (my $j=0; $j < scalar @dump_files; $j++) {
			if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
			open (DUMP_IN, "$dir_Dump_file/$dump_files[$j]"); while (<DUMP_IN>) { print OUT_coord $_; } close (DUMP_IN); 
		}
		my $array_files_msa = DOMINO::readDir($msa_dir);
		my @msa_files = @{ $array_files_msa };
		for (my $j=0; $j < scalar @msa_files; $j++) {
			if ($msa_files[$j] eq '.' || $msa_files[$j] eq '..' || $msa_files[$j] eq '.DS_Store') { next;}
			open (MSA, "$msa_dir/$msa_files[$j]"); while (<MSA>) { print OUT $_; } close (MSA); 
			if ($pyRAD_file) { print OUT "//\n"; }
	}} close(OUT_coord); close (OUT);

	## USE THE SUBROUTINE print_Excel and control if radseq_like_data
	print "+ Done...\n+ Retrieving informative locus has been done...\n+ Generating an Excel file for DOMINO markers identified...\n";
	my $excelBook = &print_Excel(\$output_file_coord, \$marker_dirname);

	unless ($avoidDelete_tmp_files) { 
		###########################
		## Delete temporary file ##
		###########################
		print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Deleting temporary files ", "#"); DOMINO::printHeader("", "#");
		remove_tree($msa_dir_tmp);
		remove_tree($dir_Dump_file);
	}	
	File::Copy::move($param_Detail_file_markers, $marker_dirname."/");
	if (-z $mapping_markers_errors_details) { remove_tree($mapping_markers_errors_details); 
	} else { File::Copy::move($mapping_markers_errors_details, $marker_dirname."/"); }	## Finish and exit
	&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; exit();
}
&time_log(); print "\n";

## Other types of data
##################################
###	Check taxa user specified  ### 
##################################
my $genome_marker_bool = 0;
my $all_markers_file = $marker_dirname."/markers.txt";
foreach my $ref_taxa (sort keys %domino_files) { ## For each taxa specified, obtain putative molecular markers
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
	mkdir $marker_dir, 0755; chdir $marker_dir; &debugger_print("Changing dir to $marker_dir");
	
	#######################################
	###		 Copy necessary files		### 
	#######################################
	print "+ Retrieve necessary files\n";
	if ($option eq "msa_alignment") { 
		push ( @{ $domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/merged.profile_ARRAY.txt");
	} else {
		my @taxa = sort @{ $domino_files{'taxa'}{'user_Taxa'} };
		my @uniq_sort_taxa = uniq(@taxa);
		my $name = join("_", @uniq_sort_taxa);
		push ( @{ $domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/$name.profile_ARRAY.txt");
	}

	##################################
	## Get the the reference fasta  ##
	##################################
	print "+ Checking the file specified as reference fasta... $domino_files{$ref_taxa}{'contigs'}[0]...OK\n";
	print "+ Generating symbolic links for necessary files\n"; 
	my $ref_Fasta = $ref_taxa.".fasta";
	system("ln -s $domino_files{$ref_taxa}{'contigs'}[0] $ref_Fasta");
	
	############################################################
	###	Generate an array information file for the reference ### 
	############################################################
	print "+ Reading the reference fasta file...\n";
	my $fasta_seqs = &retrieve_info(\@{ $domino_files{$ref_taxa}{"hash_reference_file"} }, 1);
	
	## Print each contig into a file
	print "+ Splitting file into multiple files to speed the computation...\n";
	print "+ Using parallel threads ($num_proc_user CPUs)...\n";			
	my $reference_dir = $marker_dir."/REF_DIR"; mkdir $reference_dir, 0755;
	push (@{ $domino_files{$ref_taxa}{'REF_DIR'}}, $reference_dir);
	my $size_file = DOMINO::get_size($ref_Fasta);
	my $parts2split = int($size_file/$num_proc_user);
	my $fasta_files_split = DOMINO::fasta_file_splitter($ref_Fasta, $parts2split, "fasta", $reference_dir);
		#&debugger_print("Total Size: $size_file\nCharacters to split: $parts2split"); #&debugger_print("Ref", $fasta_files_split);		
	
	## Generate a fasta for each contig
	my $pm_SPLIT_FASTA =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	my $total_fasta_files_split = scalar @{ $fasta_files_split };
	my $count_fasta_files_split=0;
	for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
		$count_fasta_files_split++;
		print "\t- Splitting file: [$count_fasta_files_split/$total_fasta_files_split]\n";
		my $pid = $pm_SPLIT_FASTA->start($i) and next;
		open(FILE, $$fasta_files_split[$i]) || die "Could not open the $$fasta_files_split[$i]...\n";
		my $size_file = $$fasta_files_split[$i]."_size";
		open (SIZE, ">$size_file");
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence; $sequence =~ s/\n//g; $titleline =~ s/\s/\t/g;
			my $file = $reference_dir."/".$titleline.".fasta";
			open(OUT, ">$file"); print OUT ">".$titleline."\n".uc($sequence)."\n"; close (OUT);
			print SIZE length($sequence)."\t".$titleline."\n";
		} close(FILE); $/ = "\n"; close (SIZE);
		$pm_SPLIT_FASTA->finish($i); # pass an exit code to finish
	}
	$pm_SPLIT_FASTA->wait_all_children;
	
	print "\n+ Concatenate contig size for each splitted file...\n";
	my $ref_Fasta_size = $ref_taxa.".size";
	my $tmp_Fasta_size = $ref_taxa.".size_tmp";
	for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
		my $size_file = $$fasta_files_split[$i]."_size";
		system ("cat $size_file >> $tmp_Fasta_size"); system ("rm $size_file");
	}
	print "\n+ Sorting contig size...\n";
	system ("sort -nr $tmp_Fasta_size >> $ref_Fasta_size"); system ("rm $tmp_Fasta_size");

	##########################################
	## 	Merge PILEUP information arrays     ##
	##########################################
	print "\n"; DOMINO::printHeader(" Fetching information from all the PROFILEs generated ", "#");
	my %pileup_files; 
	my (@clean_filtered_sam_files, @reference_bam_files, @pileup_Arrays);
	foreach my $reads (sort keys %domino_files) {
		unless ($domino_files{$reads}{'taxa'}) {next; }
		if ($domino_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa}) { push ( @reference_bam_files, @{ $domino_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa}});}		
		if ($domino_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa}) { push ( @clean_filtered_sam_files, @{ $domino_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa}});}
		next if ($reads eq $ref_taxa);
	}
	my $PILEUP_merged_folder_abs_path = $marker_dir."/PROFILE_merge_species";
	mkdir $PILEUP_merged_folder_abs_path, 0755; chdir $PILEUP_merged_folder_abs_path; 
	&debugger_print("Changing dir to $PILEUP_merged_folder_abs_path");
	push (@{ $domino_files{$ref_taxa}{'PROFILE_FOLDER'}}, $PILEUP_merged_folder_abs_path);
	
	######
	######
	###### FINDING THE GLITCH
	######	
	######
	######
	
	## get ordered size of contigs
	open (FILE, $marker_dir."/".$ref_Fasta_size);
	my @size_contigs;
	while (<FILE>) {
		chomp; my @array = split("\t", $_);
		push (@size_contigs, $array[1]);
	}
	close (FILE);	
	my $total_contigs = scalar @size_contigs; my $counter=0; my $stop=0;
	## check all or some
	if ($total_contigs > 50) { 
		$stop = 50; ## check largest ones
		print "\n\nATTENTION: There are too many contigs. DOMINO would only check for markers in the largest 50.000 ones\n";
		print "ATTENTION: Provide option -all for using the whole set\n\n";
	} else { $stop = $total_contigs; }
	if ($option_all) {	$stop = $total_contigs; }
	my $subset_offset = 5;
	
	print "+ Checking profiles of variation for each contig and merging information...\n";
	print "+ Using a sliding window approach...\n"; print "+ Using parallel threads ($num_proc_user CPUs)...\n";			

	## Check each markers using threads
	my $dir_Dump_file = $PILEUP_merged_folder_abs_path."/DUMP_files"; mkdir $dir_Dump_file, 0755;
	my $dir2print_markers = $PILEUP_merged_folder_abs_path."/MSA_fasta_tmp"; mkdir $dir2print_markers, 0755;
		
	## Check for markers: USING THREADS, ONE FOR EACH BLOCK OF CONTIGS
	my $pm_MARKER_PILEUP =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	my $SETS = int($total_contigs/$subset_offset) + 1;
	for (my $set=0; $set < $SETS; $set++) {
		my @subset_array; my $tmp=0;
		for ($counter=$counter; $counter < $total_contigs; $counter++) {
			push(@subset_array, $size_contigs[$counter]);
			$tmp++; if ($tmp eq $subset_offset) { $counter++; last; }
		}

		## SEND THREAD 
		my $pid = $pm_MARKER_PILEUP->start($set) and next;
		my (%pileup_files_threads, %contigs_pileup_fasta, @pileup_fasta);
		foreach my $reads (sort keys %domino_files) {
			next if ($reads eq $ref_taxa);
			next if ($reads eq 'taxa');
			for (my $j=0; $j < scalar @subset_array; $j++) {
				my $file = $domino_files{$reads}{"PROFILE::Ref:".$ref_taxa}[0]."/".$subset_array[$j]."_ARRAY.txt";
				unless (-e -r -s $file) { next; } # skip
				my $tmp_hash_reference = DOMINO::readFASTA_hash($file);
				my %tmp_fasta = %{$tmp_hash_reference};
				foreach my $seqs (sort keys %tmp_fasta) {
					push (@{ $contigs_pileup_fasta{$subset_array[$j]} }, $tmp_fasta{$seqs});
		}}}
		
		my $mergeProfile = $PILEUP_merged_folder_abs_path."/SET_$set"."_merged_ARRAY.txt";
		my $string = $window_size_VARS_range;$string =~ s/\:\:/-/; 
		my $string2 = $window_size_CONS_range; $string2 =~ s/\:\:/-/;	
		my $mergeCoord;
		if ($variable_divergence) { $mergeCoord = $PILEUP_merged_folder_abs_path."/SET_$set"."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
		} else { 					$mergeCoord = $PILEUP_merged_folder_abs_path."/SET_$set"."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
		}
		push (@{ $pileup_files_threads{"SET_$set"}{'mergeProfile'} }, $mergeProfile);
		push (@{ $pileup_files_threads{"SET_$set"}{'mergeCoord'} }, $mergeCoord);
		open (OUT_COORD, ">$mergeCoord");
		my $markers_shared_file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.tmp.txt";
		open (SHARED, ">$markers_shared_file");

		## NAME for output merged files
		my ($output_merged_file, $error_merged_file, $file);
		if ($variable_divergence) { 
			$file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
		} else {  
			$file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
		}		
		$output_merged_file = $file.".out"; $error_merged_file = $file.".err";
		push (@{ $pileup_files_threads{"SET_$set"}{'eachTaxaCoord'} }, $output_merged_file);

		## Merging variable and conserved information into a unique array, profile and generate coordinates
		foreach my $seqs (keys %contigs_pileup_fasta) {
			my $size = $$fasta_seqs{$seqs};
			my $tmp_string;			
			for (my $i = 0; $i < scalar $size; $i++) {
				my (@tmp, $pb);
				for (my $j = 0; $j < scalar @{ $contigs_pileup_fasta{$seqs} }; $j++){
					$pb = substr($contigs_pileup_fasta{$seqs}[$j], $i, 1);
					push (@tmp, $pb);					
				}
				my $flag = 0; my $missing_count = 0; my $missing_flag = 'N';
				my $missing_allowed_species = $minimum_number_taxa_covered;
				my $scalar_keys = scalar @tmp;
				for (my $keys = 0; $keys < scalar @tmp; $keys++) {
					if ($tmp[$keys] eq 1) { # One specie has a variable position
						$flag = 1; 
					} elsif ($tmp[$keys] eq 'N') { ## this position is missing for that specie
						$missing_count++;
					}}	
				my $tmp = $number_sp - $missing_count;
				if ($tmp < $missing_allowed_species) {
					$tmp_string .= $missing_flag; 
				} else { $tmp_string .= $flag; 
			}}		
			my $var_sites = $tmp_string =~ tr/1/1/; ## count variable sites
			my $cons_sites = $tmp_string =~ tr/0/0/; ## count conserved sites
			if ($var_sites != 0 && $cons_sites != 0) { 
				open (FH, ">>$mergeProfile"); print FH ">$seqs\n$tmp_string\n"; close(FH);
			} else { next; } 
			my $infoReturned = &sliding_window_conserve_variable(\$seqs, \$tmp_string);
			if ($infoReturned) { 
				my @array = @$infoReturned;
				my $tmp_file = "tmp_coord_set_".$set.".txt";
				open (TMP_COORD, ">$tmp_file");
				for (my $j=0; $j < scalar @array; $j++) {
					print OUT_COORD $array[$j]."\n";
					print TMP_COORD $array[$j]."\n";
					print SHARED $array[$j]."\t".$ref_taxa."\n";
				}
				close (TMP_COORD);
				} else { next; ## if empty next
			}
		
			######################################################################
			## Check the coordinates foreach taxa against the merge statistics  ##
			######################################################################
			foreach my $taxa (sort keys %domino_files) {
				unless ($domino_files{$taxa}{'taxa'}) { next; }
				if ($taxa eq $ref_taxa) {next;}	
				## For each taxa confirm profile
				my $pileup_each_taxa = $domino_files{$taxa}{"PROFILE::Ref:".$ref_taxa}[0]."/".$seqs."_ARRAY.txt";
				if (-f $pileup_each_taxa) {
					&get_coordinates_each_taxa(\$pileup_each_taxa, "tmp_coord_set_".$set.".txt", $taxa, \$output_merged_file, \$error_merged_file);
		}}}
		close (OUT_COORD);

		##########################################
		## Get Coordinates of Molecular Markers ##
		##########################################		
		my $shared_markers_others = $pileup_files_threads{"SET_$set"}{'eachTaxaCoord'}[0];
		my $markers_shared_file_sort = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.tmp_sort.txt";
		system("cat $shared_markers_others >> $markers_shared_file ; sort $markers_shared_file > $markers_shared_file_sort");
		my %coord_contig;
		open (MERGE_COORD,"<$markers_shared_file_sort") or die "Cannot open file $markers_shared_file_sort";
		while(<MERGE_COORD>){
			chomp;
			my $line = $_;
			next if ($line =~ m/^Contig\tCons1_ini:.*/o);
    		next if $line=~ m/^\s*$/o;
    		next if $line=~ /^\#/o;
    		my @array = split("\t",$line);
   			push (@{ $coord_contig{$array[0]}{$array[1].";".$array[2].";".$array[3]}}, $array[4]);
		}
		close (MERGE_COORD);	
		if (!%coord_contig) { $pm_MARKER_PILEUP->finish(); }
		
		## Print in tmp file for sorting and obtaining unique
		chdir $PILEUP_merged_folder_abs_path;
		my $tmp_file = $PILEUP_merged_folder_abs_path."/SET_$set"."_markers_shared.txt";

		open(TMP, ">$tmp_file");
		for my $scaffold (keys %coord_contig) {    		
			foreach my $marker (keys %{ $coord_contig{$scaffold} }) {
				my @taxa = @{ $coord_contig{$scaffold}{$marker} };
				## Check we find markers for all the taxa desired
				if (scalar @taxa < $minimum_number_taxa_covered) { next; }   		
				## Write DOMINO Markers Coordinates in tmp txt file
				my @string = split(";", $marker);
				my $string = join("\t", @string);
				my @sort_taxa = sort(@taxa);
				my $string_taxa = join(",", @sort_taxa); #taxa providing variation
				print TMP "$scaffold\t$string\t$string_taxa\n"; 
			}
		} close(TMP);	
		unless (-e -r -s $tmp_file) { $pm_MARKER_PILEUP->finish(); }
		
		## Collapse markers
		my $file_markers_collapse = &check_overlapping_markers($tmp_file, \$pileup_files_threads{"SET_$set"}{'mergeProfile'}[0]);

		# Retrieve fasta sequences...
		my $output_file = $PILEUP_merged_folder_abs_path."/SET_".$set."_markers_retrieved.txt";
		my $markers_print_ref = &check_DOMINO_marker($output_file, $dir2print_markers, $file_markers_collapse, $ref_taxa);
		
		unless (scalar @$markers_print_ref == 0) { 
			push (@{ $pileup_files_threads{"SET_$set"}{'markers'} }, $output_file);
			push (@{ $pileup_files_threads{"SET_$set"}{'markers_files'} }, @$markers_print_ref);
			my $dump_folder_files = $dir_Dump_file."/dump_markers_SET_".$set.".txt";
			DOMINO::printDump(\%pileup_files_threads, $dump_folder_files);	
			@pileup_fasta = (); %pileup_files_threads = ();
		}
		$pm_MARKER_PILEUP->finish(); # finish for each contig
	} 
	$pm_MARKER_PILEUP->wait_all_children; #each marker
	print "\n\n";
	print "******************************************************\n";
	print "**** All parallel parsing processes have finished ****\n";
	print "******************************************************\n\n";
	&time_log(); print "\n";

	######
	######
	###### FINDING THE GLITCH
	######	
	######	
	######	
	
	## Retrieve info of all markers identified...
	my $dump_files = DOMINO::readDir($dir_Dump_file);
	for (my $i=0; $i < scalar @$dump_files; $i++) {
		if ($$dump_files[$i] eq '.' || $$dump_files[$i] eq '..' || $$dump_files[$i] eq '.DS_Store') { next;}
		$$dump_files[$i] = $dir_Dump_file."/".$$dump_files[$i];
	}
	my %markers2retrieve;
	for (my $j=0; $j < scalar @{ $dump_files }; $j++) {
		if ($$dump_files[$j] eq '.' || $$dump_files[$j] eq '..' || $$dump_files[$j] eq '.DS_Store') { next;}
		open (DUMP_IN, "$$dump_files[$j]");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line; my @array = split("\t", $line);
			push (@{ $markers2retrieve{$array[0]}{$array[1]}}, $array[2]);
	} close (DUMP_IN); }
	# print Dumper \%markers2retrieve; #
	
	#################################################################
	## Get the information ready for the user to visualize contigs ##
	#################################################################
	print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Getting the information ready to present ", "#"); DOMINO::printHeader("", "#"); 
	my $markers_msa_folder = $marker_dir."/MSA_markers";  mkdir $markers_msa_folder, 0755; chdir $marker_dir;
	
	## Open output file for markers coordinates
	my $output_file;	
	if ($option eq "genome") { $output_file = $marker_dir."/DM_markers-coordinates.txt";		
	} else { $output_file = $marker_dir."/DM_markers-coordinates_ref_".$ref_taxa.".txt";	}
	open (OUT_markers, ">$output_file"); 
	
	## Open output file for contigs
	my $output_file_putative_contigs = $marker_dir."/DM_contigs.fasta";
	open (OUT, ">$output_file_putative_contigs"); 
	
	## Open output file for markers sequences
	my $file_coordinates = $marker_dir."/DM_sequence-markers.fasta"; 
	open (OUT_coord, ">$file_coordinates");

	push(@{ $domino_files{$ref_taxa}{'MSA_markers'}}, $markers_msa_folder);
	push(@{ $domino_files{$ref_taxa}{'markers'}}, $file_coordinates);
	push(@{ $domino_files{$ref_taxa}{'CONTIGS'}}, $output_file_putative_contigs);
	push(@{ $domino_files{$ref_taxa}{'coordinates'}}, $output_file);

	print "+ Copying reference fasta contigs...\n+ Printing reference sequences in $output_file_putative_contigs...\n";
	foreach my $subset (sort keys %markers2retrieve) {		
		## Move msa markers
		my @array_markers = @{ $markers2retrieve{$subset}{'markers_files'} };
		for (my $i=0; $i < scalar @array_markers; $i++) {
			File::Copy::move($array_markers[$i], $markers_msa_folder);	
		}
		
		my @ALL_contigs;
		# Get all contigs involved
		my $marker_contigs = $markers2retrieve{$subset}{'markers'}[0];
		open (IN_marker, "<$marker_contigs"); 
		while (<IN_marker>) { 
			print OUT_markers $_; ## print coordinates
			# Print sequence coordinates
			my $line = $_; chomp $line;
			my @array = split("\t", $line); 
			my @contig_name = split("_#", $array[0]);
			push (@ALL_contigs, $contig_name[0]);
		} close (IN_marker);
		
		# sort uniq
		my @sort_ALL_contigs = sort @ALL_contigs;
		my @uniq_contigs = uniq(@sort_ALL_contigs);
		
		## Printing contigs
		my %hash_contigs;
		for (my $c = 0; $c < scalar @uniq_contigs; $c++) {
			my $in_file = $reference_dir."/".$uniq_contigs[$c].".fasta";
			if (-e -r -s $in_file) { 
				my $hash_contigs = DOMINO::readFASTA_hash($in_file);
				$hash_contigs{$uniq_contigs[$c]} = $$hash_contigs{$uniq_contigs[$c]};
				print OUT ">".$uniq_contigs[$c]."\n".$$hash_contigs{$uniq_contigs[$c]}."\n";			
		}}
		## Get markers	
		open (IN_marker, "<$marker_contigs"); 
		while (<IN_marker>) { 
			# Print sequence coordinates
			my $line = $_; chomp $line;
			my @array = split("\t", $line); my @array1 = split(":", $array[1]); my @array3 = split(":", $array[3]);
			my @contig_name = split("_#", $array[0]);
			my ($seq_id, $seq, $array) = &fetching_range_seqs($contig_name[0], $array1[0], $array3[1], $hash_contigs{$contig_name[0]});
			print OUT_coord $seq_id."\n".uc($seq)."\n";
		} close (IN_marker);
	
	} close(OUT); close(OUT_markers); close (OUT_coord);
 	print "+ Marker development for $ref_taxa is finished here...\n\n"; &time_log(); print "\n";
} #each reference taxa

## Move parameters files
chdir $marker_dirname;

if ($genome_fasta) {
	## Print excel for clusterized results
	print "+ Generating an Excel file for DOMINO markers coordinates...\n"; 
	my $coordinates;
	foreach my $ref_taxa (sort keys %domino_files) { ## For each taxa specified, obtain putative molecular markers
		unless ($domino_files{$ref_taxa}{'contigs'}) {next; }
		$coordinates = $domino_files{$ref_taxa}{'coordinates'}[0];
	}
	my $excelbook = &print_Excel($coordinates, \$marker_dirname);
	## Move parameters and error file to folder
	File::Copy::move($param_Detail_file_markers, $marker_dirname."/");
	if (-z $mapping_markers_errors_details) { remove_tree($mapping_markers_errors_details); 
	} else { File::Copy::move($mapping_markers_errors_details, $marker_dirname."/"); }
	
	&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; exit();
}

#############################################################
## Clusterize markers using different taxa as reference ##
#############################################################
print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Clustering markers for unique results ", "#"); DOMINO::printHeader("", "#");
my $blast_dir = $marker_dirname."/clustering";
mkdir $blast_dir, 0755; chdir $blast_dir; 	&debugger_print("Changing dir to $blast_dir");

# read dir
my $files_dir_ref = DOMINO::readDir($marker_dirname);
my @markers_folders = @$files_dir_ref;

# output files
my $all_coordinates_file = $blast_dir."/all_coordinates.fasta"; open (ALL_coordinates, ">$all_coordinates_file");
my $all_contigs_file = $blast_dir."/all_contigs.fasta"; open (ALL_CONTIGS, ">$all_contigs_file"); 
print "+ Merging different DOMINO markers according to the taxa of reference...\n";

# get sequences
my @CoordMarker;
foreach my $keys (sort keys %domino_files) {
	if ($domino_files{$keys}{'taxa'}) {
		push (@CoordMarker, $domino_files{$keys}{'coordinates'}[0]);
		my $coordinate_file = $domino_files{$keys}{'markers'}[0];
		if (-e -r -s $coordinate_file) {
			my $hash = DOMINO::readFASTA_hash($coordinate_file);
			foreach my $seq (keys %{$hash}) {
				print ALL_coordinates ">".$seq."_taxa_".$keys."\n".$$hash{$seq}."\n";
		}}
		my $contig_file = $domino_files{$keys}{'CONTIGS'}[0];
		if (-e -r -s $contig_file) {
			my $size_file = DOMINO::get_size($contig_file);
			open (FH, $contig_file);
			my $chunk; read(FH, $chunk, $size_file); 
			print ALL_CONTIGS $chunk;
}}} close(ALL_CONTIGS); close(ALL_coordinates);

## Use BLAST for clustering sequences
print "+ Generate a BLAST database...\n"; 
my ($blast_DB, $blast_DB_message) = DOMINO::makeblastdb($all_coordinates_file, $BLAST, $mapping_markers_errors_details);
&debugger_print($blast_DB_message);
if ($blast_DB eq "1") {
	&printError("Early termination of the DOMINO Marker Scan...");
	my $msg= "\n\nPlease note that DOMINO could not find any markers for the parameters provided. Please Re-Run DOMINO using other parameters\n\n\n"; 
	print $msg; DOMINO::printError_log($msg); &finish_time_stamp();  exit();
} 
my $blast_search = "blast_search.txt"; print "+ BLAST search now...\n"; 
my ($blastn, $blastn_message) = DOMINO::blastn($all_coordinates_file, $blast_DB, $blast_search, $BLAST);
&debugger_print($blastn_message);

## Filter BLAST results
print "+ Filtering BLAST search now...\n";
my $contig_length_Ref = DOMINO::readFASTA_hash($all_coordinates_file);
my (%markers2keep, @markers_seen, %clusterized_contigs_keep);
my $first_hit = 0;
my $aln_overlapped = $window_size_VARS_max + $window_size_CONS_max;
print "+ Clustering markers with > $aln_overlapped bp overlapped...\n";
open (BLAST, $blast_search); while (<BLAST>) {
	my $line = $_;
	chomp $line;
	my @array = split("\t", $line);
	my $query = $array[0]; my $subject = $array[1];
	if ($query eq $subject) { next;} # same sequence
	my $query_string; my $subject_taxa;
	my $subject_string; my $taxa_query;
	if ($query =~ /(.*)\_taxa\_(.*)/) { $query_string = $1;   $taxa_query = $2; }
	if ($subject =~ /(.*)\_taxa\_(.*)/){$subject_string = $1; $subject_taxa = $2; }
	if ($subject_taxa eq $taxa_query) { next; } ## if same taxa
	if ($array[10] < 1e-20 && $array[3] > $aln_overlapped) {    ## how to clusterize...
		if ($query =~ /(.*)\_coord\_(.*)\_taxa\_(.*)/) { $clusterized_contigs_keep{$1}++;}
		if ($first_hit == 0) {
			$first_hit++;
			push( @{$markers2keep{$query_string} }, $subject_string); 
			push (@markers_seen, $subject_string);
		} else {
			my $flag_this = 0;
			foreach my $keys (sort keys %markers2keep) {
				if (grep /.*$subject_string.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
				if (grep /.*$query_string.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
			}
		if ($flag_this == 0) { 	push(@{$markers2keep{$query_string}}, $subject_string);  push (@markers_seen, $subject_string); }
}}} close (BLAST);
foreach my $keys (sort keys %$contig_length_Ref) { 
	if ($keys =~ /((.*)\_coord\_(.*))\_taxa\_(.*)/) { 
		unless (grep /$1/, @markers_seen) { 
			$markers2keep{$1}++; $clusterized_contigs_keep{$2}++;
}}}

## Printing definitely Results
my $definitely_results_dirname = $marker_dirname."/DOMINO_markers_Results";
mkdir $definitely_results_dirname, 0755; chdir $definitely_results_dirname; &debugger_print("Changing dir to $definitely_results_dirname");
print "+ Getting the information ready for the visualization of results...\n";
my $contig_def_results_sequences = "DM_contigs.fasta";
open (CONTIGS_DEF, ">$contig_def_results_sequences");
my $hash_contigs_ref = DOMINO::readFASTA_hash($all_contigs_file);
my %hash_contigs = %$hash_contigs_ref;
foreach my $keys (sort keys %hash_contigs) {
if ($clusterized_contigs_keep{$keys}) { print CONTIGS_DEF ">".$keys."\n".$hash_contigs{$keys}."\n";	}}
close(CONTIGS_DEF);

my %hash4markers2keep;
print "+ Filtering coordinates files...\n\n";
for (my $h = 0; $h < scalar @CoordMarker; $h++) {
	open (tmp_COOR, $CoordMarker[$h]);
	while (<tmp_COOR>) {
		my $line = $_;
		chomp $line;
		next if $line=~ m/^\s*$/o;
		if ($line =~ /.*Vari/) {next;}
		my @array_split = split("\t",$line);
		my @contig_split = split("\_#",$array_split[0]);
		if ($clusterized_contigs_keep{$contig_split[0]}) {
			my @array_start = split(":", $array_split[1]);
			my @array_end = split(":", $array_split[3]);
			my $coord_string = $array_start[0]."_".$array_end[1];
			my $coord2check = $contig_split[0]."_coord_".$coord_string;
			if ($markers2keep{$coord2check}) {
				$hash4markers2keep{$contig_split[0]}{$coord_string} = $line;
}}} close(tmp_COOR); }

my $coordinates_def_results = "DM_markers-coordinates.txt";
open (COOR, ">$coordinates_def_results");
print COOR "Contig\t\tConserved_Region\tVariable_Region\tConserved_Region\tMapping_Taxa\Length\tDivergence\n";
my %rename;

foreach my $contigs (sort keys %hash4markers2keep) {
	my $counter = 0;
	foreach my $markers (sort keys %{ $hash4markers2keep{$contigs} }) {
		$counter++;
		my $string2change = $hash4markers2keep{$contigs}{$markers};
		my $stringchanged;
		my @array = split("\t", $string2change);
		my $marker2search = $array[0];
		$marker2search =~ s/_#/_marker_/;
		my $tmp = $contigs."_#".$counter;
		$array[0] = $tmp;
		$rename{$marker2search} = $tmp;
		for (my $i=0; $i < scalar @array; $i++) { $stringchanged .= $array[$i]."\t"; }			
		print COOR $stringchanged."\n";
} print COOR "\n"; }
close(COOR); &time_log(); print "\n";

## Get MSA markers
my $markers_msa_folder = $definitely_results_dirname."/MSA_markers"; mkdir $markers_msa_folder, 0755;
foreach my $keys (sort keys %domino_files) {
	if ($domino_files{$keys}{'taxa'}) {
		my $MSA_markers_each_Taxa = $domino_files{$keys}{'MSA_markers'}[0];
		my $array_Ref = DOMINO::readDir($MSA_markers_each_Taxa);
		for (my $i=0; $i < scalar @{ $array_Ref }; $i++) {
			my @name = split(".fasta", $$array_Ref[$i]);
			if ( $rename{$name[0]} ) {
				File::Copy::copy($MSA_markers_each_Taxa."/".$$array_Ref[$i], $markers_msa_folder."/".$rename{$name[0]}.".fasta");	
}}}}

## Print excel for clusterized results
print "+ Generating an Excel file for DOMINO markers coordinates...\n"; 
&print_Excel(\$coordinates_def_results, \$definitely_results_dirname);

## Dumping info to a file
my $dump_info_DOMINO = $marker_dirname."/DOMINO_dump_information.txt"; DOMINO::printDump(\%domino_files, $dump_info_DOMINO);
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
	remove_tree($blast_dir);
}

# ToDo
if ($keepbam) { print "Keepbam option not yet implemented....\nSorry for that...\n"; }

## Move parameters and error file to folder
File::Copy::move($param_Detail_file_markers, $marker_dirname."/");
if (-z $mapping_markers_errors_details) { remove_tree($mapping_markers_errors_details); 
} else { File::Copy::move($mapping_markers_errors_details, $marker_dirname."/"); }

## Finish and exit
&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; 
exit();


######################################################################################
##																					##	
##									SUBROUTINES										##
##																					##
######################################################################################

sub binomial {
	## Returns the probability of an event ocurring $k times in $n attempts where the probability
	## of it occurring in a sample attempt is $p
    my $k = shift(@_); ## Occurrences
    my $n = shift(@_); ## Attempts
    my $p = shift(@_); ## Probability
    my $prob = ($p**$k) * ((1 - $p)**($n - $k)) * &functional_factorial($n) / (&functional_factorial($k) * &functional_factorial($n - $k));
    return $prob;
}

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
			} else { &printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
		} else { ## File is not named fasta/fa or fastq/fq
			my $format_returned = DOMINO::check_file_format($file_to_check);
			if ($format_returned eq "fasta") { ## It is actually a fasta file
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} elsif ($format_returned eq "fastq") { ##fastq file is ok
				print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
			} else { &printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
		}
		
		## Get file name
		my @file_to_check = split("/", $file_to_check);
		my $file_name = $file_to_check[-1];
		if ($id) {
			print "\t\tIt also contains an identifier ($id) in the name for later analysis...OK\n";
		} else { 
			if ($genome_fasta) {
				&printError("$file_to_check...\n File does not contain any name provided using taxa_names option.\nMaybe it is under controlled but just bear it in mind...\nDOMINO would not die here...\n"); 
		}}
	} else { &printError("Please provide several files for each taxa..."); DOMINO::printFormat_message(); DOMINO::dieNicely(); }
}

sub check_given_marker {
	
	my $array_Coord = $_[0];
	my $dna_seq = $_[1];
	my $total_length_sub = length($dna_seq);
	
	my @coordinates = @{ $array_Coord };
	my $coord_P1 = $coordinates[0]; my $coord_P2 = $coordinates[1];
	my $coord_P3 = $coordinates[2]; my $coord_P4 = $coordinates[3];
	my $coord_P5 = $coordinates[4]; my $coord_P6 = $coordinates[5];
	
	# Conserved
	my $length_string_P1_P2 = $coord_P2 - $coord_P1;
	my $string_P1_P2 = substr($dna_seq, $coord_P1, $length_string_P1_P2);
	my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//; ## Conserved 
	
	if ($length_string_P1_P2 < $window_size_CONS_min) {return 1;}
	if ($length_string_P1_P2 > $window_size_CONS_max) { 
		#print Dumper $array_Coord; print "ERROR length_string_P1_P2! $length_string_P1_P2 > $window_size_CONS_max\n"; return 1;
		if ($window_size_CONS_max != $window_size_CONS_min) { return 1; } ## If user specifies a range..
	}
	if ($count_string_P1_P2 > $window_var_CONS) { 
		#print Dumper $array_Coord; print "ERROR count_string_P1_P2! $count_string_P1_P2 > $window_var_CONS\n";
		return 1;
	} 
		
	# Variable
	my $length_string_P3_P4 = $coord_P4 - $coord_P3;
	if ($length_string_P3_P4 > $window_size_VARS_max) { 
		#print Dumper $array_Coord; print "ERROR! length_string_P3_P4 $length_string_P3_P4 > $window_size_VARS_max\n"; 
		return 1;
	}
	my $string_P3_P4 = substr($dna_seq, $coord_P3, $length_string_P3_P4);
	my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//; ## Variable
	
	if ($variable_divergence) {
		# If a minimun divergence, get the expected variable sites for the length
		my $expected_var_sites = int($variable_divergence * $total_length_sub);
		unless ($count_string_P3_P4 >= $expected_var_sites) { return 1; }			
	}

	# Conserved
	my $length_string_P5_P6 = $coord_P6 - $coord_P5;
	my $string_P5_P6 = substr($dna_seq, $coord_P5, $length_string_P5_P6);
	if (!$string_P5_P6) {
		#print Dumper $array_Coord; #print "Length: $length_string_P5_P6\n"; #print $dna_seq;
		return 1;
	}
	my $count_string_P5_P6 = $string_P5_P6  =~ tr/1//; ## Conserved
	if ($length_string_P5_P6 < $window_size_CONS_min) { return 1;}
	if ($length_string_P5_P6 > $window_size_CONS_max) { 
		#print Dumper $array_Coord; print "ERROR! length_string_P5_P6 $length_string_P5_P6 > $window_size_CONS_max\n"; return 1;
		if ($window_size_CONS_max != $window_size_CONS_min) { return 1; } ## If user specifies a range..
	}
	if ($count_string_P5_P6 > $window_var_CONS) {
		#print Dumper $array_Coord; print "ERROR! count_string_P5_P6 $count_string_P5_P6 > $window_var_CONS\n";
		return 1;
	}
	my $total_length = $coord_P6 - $coord_P1;
	return $total_length;
}

sub check_overlapping_markers {

	## Overlaps and maximizes domino markers obtained
	my $file = $_[0]; my $mergeArray_file = $_[1];
	&debugger_print("Checking file $file");
	my $contig_id; my %tmp_hash; my $marker_counter_tmp = 0;
	my @sequences;
	open (FILE, $file);
	while (<FILE>) {
		my $line = $_;
		chomp $line; 
		$line =~ s/ /\t/;
		my @array_lines = split ("\t", $line);		
		$contig_id = $array_lines[0];
		my @a = split(":", $array_lines[1]); ## conserved region1 coord
		my @b = split(":", $array_lines[2]); ## variable region coord
		my @c = split(":", $array_lines[3]); ## conserved region2 coord
		my @coordinates = ($a[0],$a[1],$b[0],$b[1],$c[0],$c[1]);
		my $string = join(",", @coordinates);
		my $taxa;
		if ($array_lines[4]) { $taxa = $array_lines[4];
		} else { $taxa = $MID_taxa_names;  }		
		push (@{ $tmp_hash{$contig_id}{$taxa}{$a[0]} }, $string); ## Keep record of taxa and coordinates
	} 
	close(FILE);

	# Debug 
	my %coord_seen;
	foreach my $contig (sort keys %tmp_hash) {
		foreach my $taxa (sort keys %{ $tmp_hash{$contig} }) {		
			foreach my $marker (keys %{ $tmp_hash{$contig}{$taxa} }) {			
				if ($coord_seen{$contig}{$taxa}{$marker}) {next;}
				my $bool = 1;
				my ($counter, $bad_counter) = 0;
				while ($bool) {
					$counter += $SLIDING_user;
					my $new_coord = $marker + $counter;
					if ($coord_seen{$contig}{$taxa}{$new_coord}) {next;}
					if ($tmp_hash{$contig}{$taxa}{$new_coord}) {
						push (@{ $tmp_hash{$contig}{$taxa}{$marker} }, @{ $tmp_hash{$contig}{$taxa}{$new_coord} });
						$coord_seen{$contig}{$taxa}{$new_coord}++;
						$tmp_hash{$contig}{$taxa}{$new_coord} = 1;
					} else {
						$bad_counter++;
						if ($bad_counter > 3) { ## We would consider the same marker if overlapping 3 SLIDING_user!!
							($bool,$counter,$bad_counter) = 0;
	}}}}}}
	my %tmp_coord;
	foreach my $contig (sort keys %tmp_hash) {
		foreach my $taxa (keys %{ $tmp_hash{$contig} }) {
			foreach my $marker (keys %{ $tmp_hash{$contig}{$taxa} }) {
				if ($tmp_hash{$contig}{$taxa}{$marker} == 1) {next;}
				my @array = sort(@{ $tmp_hash{$contig}{$taxa}{$marker} });
				my @sort_uniq = uniq(@array);
				push (@{ $tmp_coord{$contig}{$taxa}{$marker} }, @sort_uniq);
	}}}
	undef %coord_seen; undef %tmp_hash; ## release RAM

	## Set range values
	my $range = $window_size_VARS_max - $window_size_VARS_min; my @length;
	if ($range >= 500) { 	 @length = (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500); 
	} elsif ($range < 500) { @length = (50, 100, 200, 300, 400, 500); }

	# Debug print Dumper \%tmp_coord
	my @array = split(".txt", $file); my $hash_ref = DOMINO::readFASTA_hash($mergeArray_file); 
	my $file2return = $array[0]."_overlapped_Markers.txt"; open (OUT, ">$file2return"); # print $file2return."\n";	#

	foreach my $contig (sort keys %tmp_coord) {
		foreach my $taxa (keys %{ $tmp_coord{$contig} }) {
			foreach my $keys_markers (keys %{ $tmp_coord{$contig}{$taxa} }) {		
				my @array_coordinates = @{ $tmp_coord{$contig}{$taxa}{$keys_markers} };
				my (@good_ones, %hash2print, %marker_seen); 
				for (my $i=0; $i < scalar @array_coordinates; $i++) {
					my $coordinate = $array_coordinates[$i];
					my @array_coord = split (",", $coordinate);
					for (my $j = (scalar @array_coordinates) - 1; $j >= 0; $j--) {
						my $coordinate2 = $array_coordinates[$j];
						my @array_coord2 = split (",", $coordinate2);
						my @array2check = ($array_coord[0], $array_coord2[1], $array_coord2[2], $array_coord2[3], $array_coord2[4], $array_coord2[5]);
						my $string = join(":", @array2check).":".$taxa;
						#print $$hash_ref{$contig}."\n";
						if ($marker_seen{$string}) {next;}
						my $result = &check_given_marker(\@array2check, $$hash_ref{$contig});						
						if ($result ne 1) {
							my $id;
							for (my $j = 0; $j < scalar @length; $j++) {
								if ($result <= $length[$j] ) { 		$id = "$length[$j - 1] - $length[$j]"; last;
								} elsif ($result > $length[-1]) {	$id = "bigger"; last;
							}}
							push (@{ $hash2print{$array_coord[0]}{$id}}, $string);
							# print $keys_markers."\t".$array_coord[0]."\t".$array_coord2[5]."\t".$result."\t".$id."\t".$string."\n"; 
							$marker_seen{$string}++;
				}}}
				foreach my $keys (keys %hash2print) {
					foreach my $lent (sort keys %{ $hash2print{$keys} }) {
						my @array = @{ $hash2print{$keys}{$lent} };					
						for (my $i=0; $i < scalar @array; $i++) { 
							print OUT $contig."##".$array[$i]."\n";
						} print OUT "//\n";
	}}}}}
	close(OUT); 
	return $file2return;
}

sub check_DOMINO_marker {
	
	my $file = $_[0]; my $dir = $_[1];	
	my $DOMINO_markers_file = $_[2]; 
	my $ref_taxa_all = $_[3]; # if MSA alignment it is a file containing msa

	my @files; 
	
	## Check each group of overlapping markers
	open (MARKERS, "$DOMINO_markers_file") or die "Could not open file $DOMINO_markers_file";
	my $j = 1; my %hash_array;
	while (1) {
		my $chunk;
		if (!eof(MARKERS)) { 
			$/ = "\/\/"; ## Telling perl where a new line starts
			$chunk = <MARKERS>; 
			push(@{$hash_array{$j}}, $chunk);
		} 
		$j++; 
		last if eof(MARKERS);
	}
	$/ = "\n";	
	#print Dumper \%hash_array; 

	my $marker_number = 0;
	foreach my $group (sort {$a <=> $b} keys %hash_array) {
		## Split chunk push into array
		my @array_markers_tmp = @{ $hash_array{$group} };
		my @array_markers;
		for (my $i=0; $i < scalar @array_markers_tmp; $i++) {
			#print "CHUNK!\n"; print $array_markers_tmp[$i]."\n"; print "********\n";
			my @array = split("\n", $array_markers_tmp[$i]);
			for (my $j=0; $j < scalar @array; $j++) {
				chomp($array[$j]);
				if ($array[$j] =~ ".*:.*") {
					push(@array_markers, $array[$j]);
		}}}
		## For each group of markers, evaluate the features
		for (my $i=0; $i < scalar @array_markers; $i++) {
			my @contig_name = split("##", $array_markers[$i]);
			my @array = split(":", $contig_name[1]);
			my $coord_P1 = $array[0]; my $coord_P2 = $array[1];
			my $coord_P3 = $array[2]; my $coord_P4 = $array[3];
			my $coord_P5 = $array[4]; my $coord_P6 = $array[5];
			my $taxa = $array[6];
			# print $contig_name[0]."\t".join(",",@array)."\n"; exit();
			
			## Retrieve msa for this marker			
			my %hash; my %fasta_msa_sub;
			if ($option eq "msa_alignment") {
				my $fasta_msa_sub_ref = DOMINO::readFASTA_hash($ref_taxa_all);
				%fasta_msa_sub = %{ $fasta_msa_sub_ref };
			} else {
				my @desired_taxa = split(",", $taxa);
				for (my $j=0; $j < scalar @desired_taxa; $j++) {
					next if ($ref_taxa_all eq $desired_taxa[$j]);
					my $fasta_each_taxa = $domino_files{$desired_taxa[$j]}{"PROFILE::Ref:".$ref_taxa_all}[0]."/".$contig_name[0]."_sequence.fasta";
					if (-f $fasta_each_taxa) {
						$fasta_msa_sub{$desired_taxa[$j]} = $fasta_each_taxa;
					} else { 
						$fasta_msa_sub{$desired_taxa[$j]} = "missing"; 
				}}
				my $reference_file_contig = $domino_files{$ref_taxa_all}{'REF_DIR'}[0]."/".$contig_name[0].".fasta";
				$fasta_msa_sub{$ref_taxa_all} = $reference_file_contig;
			}

			foreach my $keys (sort keys %fasta_msa_sub ) {			
				if ($fasta_msa_sub{$keys} =~ /missing/) {next;}
				if ($fasta_msa_sub{$keys} =~ /.*fasta/) {
					open (FILE, $fasta_msa_sub{$keys});
					my ($titleline, $sequence);
					$/ = ">"; ## Telling perl where a new line starts
					while (<FILE>) {		
						next if /^#/ || /^\s*$/;
						chomp;
						($titleline, $sequence) = split(/\n/,$_,2);
						next unless ($sequence && $titleline);
						chomp $sequence;
						$sequence =~ s/\s+//g; $sequence =~ s/\r//g; $titleline =~ s/\r//g;
					}
					close(FILE); $/ = "\n";
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $coord_P1, $coord_P6, $sequence);
					$hash{$keys} = $seq;
				} else {
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $coord_P1, $coord_P6, $fasta_msa_sub{$keys});
					$hash{$keys} = $seq;
			}}
			my $seq_name;
			
			## Check pairwise MSA
			my $valueReturned = &check_marker_pairwise(\%hash);
			if ($valueReturned == 1) { ## if it is variable for each pairwise comparison
				## Get variable positions for the whole marker
				my $array_ref_returned = &check_marker_ALL(\%hash, "Ref");
				if ($array_ref_returned eq 'NO') { 
					remove_tree($msa_file);
				} else {
					$marker_number++;
					my $msa_file = $dir."/".$contig_name[0]."_marker_".$marker_number.".fasta";
					open (MSA_OUT, ">$msa_file");
					foreach my $seqs (sort keys %hash) { print MSA_OUT ">".$seqs."\n".$hash{$seqs}."\n"; }
					close (MSA_OUT);
					push (@files, $msa_file);
					my @arrayReturned = @$array_ref_returned; 	
					## 0: taxa join(::) ## 1: variable sites  ## 2: length of the marker ## 3: profile string ## 4: effective length
	
					my $msa_file_array = $dir."/".$contig_name[0]."_marker_".$marker_number."_ARRAY.txt";
					open (ARRAY, ">$msa_file_array"); print ARRAY ">".$contig_name[0]."_marker_".$marker_number."_ARRAY\n".$arrayReturned[3]."\n"; close (ARRAY);
	
					## Print into a file						
					my $percentage = ($arrayReturned[1]/$arrayReturned[4])*100;
					my $variation_sprintf = sprintf ("%.3f", $percentage);
					my $marker2write = $contig_name[0]."_#$marker_number";
					my @array2print = ($marker2write, $array[0].":".$array[1], $array[2].":".$array[3], $array[4].":".$array[5], $arrayReturned[2], $variation_sprintf, $array[6]); 
					my $string = join("\t", @array2print);
					open (OUT, ">>$file"); print OUT $string."\n"; close(OUT); last;		
				}
			}
		}
	} close (MARKERS);	
	return (\@files);
}

sub check_marker_pairwise {

	## Given a MSA file in a hash, check each taxa pairwise
	my $hash_ref = $_[0];
	my @taxa = keys %$hash_ref;
	my (%seen, %pairwise, %discard);
	for (my $j=0; $j < scalar @taxa; $j++) {
		my $reference = $taxa[$j];
		if (!$domino_files{$reference}{'taxa'}) {next;}
		foreach my $keys (sort keys %$hash_ref) {
			if (!$domino_files{$keys}{'taxa'}) {next;}
			if ($seen{$keys}) {next;}
			if ($discard{$keys}) {next;}
			if ($keys eq $reference) {next;}
			my $seq2check = $$hash_ref{$keys};			
			my @array_reference = split("",$$hash_ref{$reference});
			my @array2check = split("", $$hash_ref{$keys});
			my @array;
			for (my $i=0; $i < scalar @array_reference; $i++) {
				my $reference_nuc = $array_reference[$i];
				my $base2check = $array2check[$i];
				push(@array, &check_reference_bp($reference_nuc, $base2check));
			}
			my $string = join "", @array;
			my $var_sites_sub = $string =~ tr/1/1/;
			my $con_sites_sub = $string =~ tr/0/0/;
			my $total_sub = $con_sites_sub + $var_sites_sub;			
			if ($var_sites_sub == 0) { ## If does not fit the necessary divergence
				$seen{$reference}++; $discard{$keys}++; $pairwise{$reference}{$keys} = "NO!";
				if ($variable_divergence) { &debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t0\t".$variable_divergence);
				} else { 					&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t0\t".$variable_positions_user_range);
				}
				next;
			}
			my $percentage_sub = $var_sites_sub/$total_sub;
			if ($variable_divergence) { &debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t".$percentage_sub."\t".$variable_divergence);
			} else {					&debugger_print($reference."\t".$keys."\t".$var_sites_sub."\t".$con_sites_sub."\t".$percentage_sub."\t".$variable_positions_user_range);
			}
			if ($variable_positions_user_min) {
				if ($var_sites_sub > $variable_positions_user_min) { 
					if ($var_sites_sub < $variable_positions_user_max) { 
						$pairwise{$reference}{$keys} = "YES!";
						&debugger_print("YES");
					} else {&debugger_print("NO");}
				} else {
					## If does not fit the necessary divergence
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
					$pairwise{$reference}{$keys} = "NO!";
					&debugger_print("NO");
				}
			} elsif ($variable_divergence) {
				if ($percentage_sub > $variable_divergence) { 
					$pairwise{$reference}{$keys} = "YES!";
					&debugger_print("YES");
				} else {
					## If does not fit the necessary divergence
					&debugger_print("NO");
					$pairwise{$reference}{$keys} = "NO!";
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
		}}}
		$seen{$reference}++; ## avoid checking again the same pair
	}
	my $flag_fitting = 0;
	foreach my $keys (sort keys %pairwise) {
		if ($discard{$keys}) {next;}
		foreach my $k (sort keys %{$pairwise{$keys}}) {
			if ($discard{$k}) {next;}
			if ($pairwise{$keys}{$k} eq 'YES!') {
				$flag_fitting++;
	}}}

	#print Dumper \%pairwise;
	if ($number_sp == 2) {
		if ($flag_fitting == 1) {
			#print "YES!\n"; 
			return '1'; 
		} else {
			#print "NO!\n";
			return '0'; 
	}} else {
		if ($flag_fitting < $minimum_number_taxa_covered) {
			#print "NO!\n";
			return '0'; 
		} else {
			#print "YES!\n"; 
			return '1'; 
	}}	
}

sub check_marker_ALL {
	##########################################################################################
	##	This function checks each region using all taxa provided and generates the 			##	
	##	variation profile for the whole set													##
	## 		        																		##
	##	jfsanchezherrero@ub.edu 09/02/2016													##
	## 		        																		##
	##########################################################################################

	my $file = $_[0];	my $ref = $_[1];
	#print "check_marker_ALL $file\n";
	my (%hash, $length, @taxa, @length_seqs);
	if ($ref) { 
		foreach my $seqs (sort keys %{ $file }) {
			my @array = split("", $$file{$seqs});
			if (!$domino_files{$seqs}{'taxa'}) {next;}
			push (@{ $hash{$seqs}}, @array);
			$length = scalar @array;
			push (@length_seqs, $length);
			push (@taxa, $seqs);
		}	
	} else {
		open(FILE, $file) || die "Could not open the $file ...\n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence;
			$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
			$titleline =~ s/\r//g;
			#print $titleline."\n".$sequence."\n";
			my @array = split("", $sequence);
			if (!$domino_files{$titleline}{'taxa'}) {next;}
			push (@{ $hash{$titleline}}, @array);
			$length = scalar @array;
			push (@length_seqs, $length);
			push (@taxa, $titleline);
		}
		close(FILE); $/ = "\n";	
	}

	my @tmp_length = sort @length_seqs;
	my @tmp_length_uniq = uniq(@tmp_length);	
	if (scalar @tmp_length_uniq > 1) {
		&printError("There is problem: length of the markers do not match for $file..."); return "";
	} else { $length_seqs[0] = $length; }	
	
	my @profile;
	for (my $i=0; $i < $length; $i++) {
		my $flag_position = 0;
		my @tmp;
		foreach my $seqs (sort keys %hash) { 
			push (@tmp, $hash{$seqs}[$i]);
		}
		my @tmp_uniq = sort @tmp;
		my @tmp_uniq_sort = uniq(@tmp_uniq);
		if (scalar @tmp_uniq_sort == 1) {
			if ($tmp_uniq_sort[0] eq 'N') {
				push (@profile, 'N');
			} elsif ($tmp_uniq_sort[0] eq '-') {
				push (@profile, '-');
			} elsif ($ambiguity_DNA_codes{$tmp_uniq_sort[0]}) {
				push (@profile, '1');
			} else { push (@profile, '0'); }
		} else {
			## We are assuming the calling has been correctly done and
			## the ambiguity codes are due to polymorphism		
			my $escape_flag = 0;
			my (@tmp, @amb);
			for (my $j=0; $j < scalar @tmp_uniq_sort; $j++) {
				if ($tmp_uniq_sort[$j] eq '-') { ## Gaps would be codify as -
					push (@tmp, '-'); $escape_flag++;
				} elsif ($tmp_uniq_sort[$j] eq 'N') { ## Gaps would be codify as -
					push (@tmp, 'N'); $escape_flag++;
				} elsif ($ambiguity_DNA_codes{$tmp_uniq_sort[$j]}) {
					push(@amb, $tmp_uniq_sort[$j]);
				} else { push(@tmp, $tmp_uniq_sort[$j]); }
			}
			if ($escape_flag) { push (@profile, '-');			
			} else {
				if (scalar @amb == 0) { ## No ambiguous code
					push (@profile, '1');			
				} elsif (scalar @amb == 1) { ## 1 amb code
					for (my $i=0; $i < scalar @amb; $i++) {
						my $flag_yes = 0;
						for (my $k = 0; $k < scalar @{ $ambiguity_DNA_codes{$amb[$i]}}; $k++) {
							if (grep /$ambiguity_DNA_codes{$amb[$i]}[$k]/, @tmp) { $flag_yes++; }
						}
						if ($flag_yes > 0) {
							if ($polymorphism_user) { 	push (@profile, '1'); 	## if polymorphism
							} else { 					push (@profile, '0'); } ## The ambiguous is the present snps: 		YCT => [ Y > C/T ]
						} else { 						push (@profile, '1');	## The ambiguous is not the present snps  	MGG => [ M > A/C ]
				}}} elsif (scalar @amb > 1) {  			push (@profile, '1');   ## Several
	}}}}
	my $string = join ("", @profile); #print "\t\t\t  ".$string."\n";
	my $var_sites = $string =~ tr/1/1/; ## count variable sites
	if ($var_sites == 0) { return 'NO'; }
	my $con_sites = $string =~ tr/0/0/; ## count conserved sites
	my $count_length = $con_sites + $var_sites;
	my $missing = $length - $count_length;
	my $missing_allowed_length = $missing_allowed * $length;
	if ($missing > $missing_allowed_length) { return 'NO';}
	my $species = join (",", sort @taxa);
	my @array = ($species, $var_sites, $length, $string, $count_length);
	#print Dumper \@array; print "\n";
	return \@array;
}

sub check_reference_bp {

	my $reference_nuc = $_[0];
	my $base2check = $_[1];
	
	## For this taxa
	if ($base2check eq "-" || $reference_nuc eq "-") { 
		return '-';
	} elsif ($base2check eq "N" || $reference_nuc eq "N") { 
		return 'N';
	} elsif ($reference_nuc ne $base2check) { 
		## Check wether there is an ambiguity code or not
		## and also if user would like some polymorphism
		## and decide whether it is a variable or conserved position
		my $flag = 1;
		if ($ambiguity_DNA_codes{$reference_nuc} || $ambiguity_DNA_codes{$base2check}) { ## one or the other or both
			if ($ambiguity_DNA_codes{$reference_nuc} and $ambiguity_DNA_codes{$base2check}) {
				## Both bases are ambiguous
				for (my $h = 0; $h < scalar @{ $ambiguity_DNA_codes{$base2check}}; $h++) {
					for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$reference_nuc}}; $j++) {
						if ($ambiguity_DNA_codes{$reference_nuc}[$j] eq $ambiguity_DNA_codes{$base2check}[$h]) {
							$flag = 0;
						} else {
							if ($polymorphism_user) {
								$flag = 1;	last;
			}}}}
			} elsif ($ambiguity_DNA_codes{$reference_nuc}) {
				for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$reference_nuc}}; $j++) {
					if ($ambiguity_DNA_codes{$reference_nuc}[$j] eq $base2check) {
						$flag = 0;								
					} else {
						if ($polymorphism_user) {
							$flag = 1;	last;
				}}}
			} elsif ($ambiguity_DNA_codes{$base2check}) {
				for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$base2check}}; $j++) {
					if ($ambiguity_DNA_codes{$base2check}[$j] eq $reference_nuc) {
						$flag = 0;								
					} else {
						if ($polymorphism_user) {
							$flag = 1;	last;
				}}}}
		} else { 
			## If neither the reference or the base to check are ambiguous
			## and they are different, this would be a variable site
			$flag = 1;
		}
		if ($flag == 0) { return '0'; } else { return '1'; }

	} elsif ($reference_nuc eq $base2check) { 					## Both bases are the same, this is a conserved site
		return '0';
	} else {
		return '-';
	}
}

sub clean_tmp_files_marker_dir {
	
	##########################################################################################
	##	 																					##
	##  This function 																		##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $marker_folder = $_[0];
	my $files_dir_ref = DOMINO::readDir($marker_folder);
	my @tmp_files = @$files_dir_ref;
	my $tmp_dir = "intermediate_files"; mkdir $tmp_dir, 0755;
	my $pileup_merged;
	
	for (my $i = 0; $i < scalar @tmp_files; $i++) {
		## Skip
		if ($tmp_files[$i] eq ".DS_Store" || $tmp_files[$i] eq "." || $tmp_files[$i] eq ".." ) { next; }
		
		## Intermediate files 
		if ($tmp_files[$i] =~ /ARRAY\_files\_.*/)  { remove_tree ($tmp_files[$i]);  next;}
		if ($tmp_files[$i] =~ /PROFILE\_merge\_.*/) { $pileup_merged = $tmp_files[$i]; next;}

		## Output results
		if ($tmp_files[$i] =~ /DM_.*/) { next; } 
		if ($tmp_files[$i] =~ /Error_contigs.txt/) { next; } 

		## Temporary files 
		if ($tmp_files[$i] =~ /.*merge.*/ ) { File::Copy::move($tmp_files[$i], $tmp_dir); next;
 		} elsif ($tmp_files[$i] =~ /\.profile.*/ ) { File::Copy::move($tmp_files[$i], $tmp_dir); next;
 		} else { File::Copy::move($tmp_files[$i], $tmp_dir); next; }
	}

	## Delete temporary folder	
	unless ($avoidDelete_tmp_files) { 
		remove_tree($tmp_dir); remove_tree($pileup_merged);
	}
}

sub delete_files_mapping {  
	
	my $dir_to_delete = $_[0];
	my $ref_id = $_[1];
	my $files_dir_ref = DOMINO::readDir($dir_to_delete);
	my @mapping_files = @$files_dir_ref;
	my $temp_dir = "temp_dir";
	mkdir $temp_dir, 0755; 
	
	for (my $i=0; $i < scalar @mapping_files; $i++) {		
		if ($mapping_files[$i] eq ".DS_Store" || $mapping_files[$i] eq "." || $mapping_files[$i] eq ".." ) { next; }
		if (-z $mapping_files[$i]) {File::Copy::move($mapping_files[$i], $temp_dir)
		}elsif ($mapping_files[$i] =~ /$ref_id.*clean_filtered.sorted.bam$/) {next; 
		}elsif($mapping_files[$i] =~ /.*pileup\_ARRAY\.txt/) { File::Copy::move($mapping_files[$i], $temp_dir);
		}elsif($mapping_files[$i] =~ /ARRAY\_files\_.*/) { next; 
		}elsif($mapping_files[$i] =~ /.*clean\_filtered\.sam/) { next; 
		}elsif($mapping_files[$i] =~ /.*coverage\_stats\_After\_filtering\.txt/) { next; 
		}elsif($mapping_files[$i] =~ /tmp.*/) { remove_tree($mapping_files[$i]); 
		}else{ File::Copy::move($mapping_files[$i], $temp_dir);}
	}
	remove_tree($temp_dir);
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
		if ($files[$i] =~ /.*id-(.*)(\_R\d+)\.fastq/g) {
			if ($domino_files{$1}{'taxa'}) { 
				push (@{ $domino_files{$1}{'reads'} }, $files[$i]); ## push the whole file path			
			} else { &printError("Please check the tag for the file $files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
		} elsif ($files[$i] =~ /.*id-(.*)\.fastq/g) {
			if ($domino_files{$1}{'taxa'}) { 
				push (@{ $domino_files{$1}{'reads'} }, $files[$i]); ## push the whole file path			
			} else { &printError("Please check the tag for the file $files[$i] \n...not matching any taxa name provided..."); DOMINO::dieNicely(); }
		}
	}
}

sub fetching_range_seqs {
	
	my $contig_name = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $sequence = $_[3];
	
	my $length = $end - $start;	
	my $sub_seq = substr($sequence, $start, $length);
	my $seq_id = ">".$contig_name."_coord_".$start."_".$end;
	return($seq_id, $sub_seq);	
}

sub fetching_range_seqs_array_Pileup {
	my $contig_name = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $folder_ref = $_[3];
	my $folder_name = $_[4];
	
	my @array_folder = @$folder_ref;	
	for (my $i=0; $i < scalar @array_folder; $i++) {		
		if ($array_folder[$i] eq "." || $array_folder[$i] eq ".." || $array_folder[$i] eq ".DS_Store" ) {next;}
		if ($array_folder[$i] =~ /$contig_name\_ARRAY\.txt/) {
			my $ref_hash_array = DOMINO::readFASTA_hash($folder_name."/".$array_folder[$i]);
			my $array_seq_desired = $$ref_hash_array{$contig_name};
			my $real_start = $start - 1;
			my $length = $end - $real_start;	
			my $sub_array = substr ($array_seq_desired, $real_start, $length);
			return $sub_array;
}}}

sub finish_time_stamp {

	my $finish_time = time;
	print "\n\n"; DOMINO::printHeader("","+"); 
	DOMINO::printHeader(" ANALYSIS FINISHED ","+"); 
	DOMINO::printHeader("","+"); 
	print DOMINO::time_stamp();
	my $secs = $finish_time - $start_time; 
	my $hours = int($secs/3600); 
	$secs %= 3600; 	
	my $mins = int($secs/60); 
	$secs %= 60; 
	printf (" Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub functional_factorial {

    my $n = shift(@_);
    my $fact = 1;
    if (($n < 0) or (170 < $n)) { die "Factorial out of range"; }
    for (my $i = 1; $i <= $n; $i++) { $fact *= $i; }
    return $fact;
}

sub get_amb_code {
	my $ref_hash = $_[0];
	my @array = keys %{$ref_hash};
	my @array_sorted = sort @array;
	my $string = join "",@array_sorted;
	if ($ambiguity_DNA_codes_reverse{$string}) {
		return $ambiguity_DNA_codes_reverse{$string};
	} else { return "N"; }		
}

sub get_clean_files {
	
	my $clean_folder = DOMINO::get_earliest("clean_data", $folder_abs_path);
	if ($clean_folder eq 'clean_data') {
		&printError("No clean_data folder was found. Please Re-Run DOMINO mapping step or make sure you use the correct output directory..."); DOMINO::dieNicely();
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

sub get_coordinates_each_taxa {
	
	##########################################################################################
	##	 																					##
	##  This function checks the variation of each taxa against the merge PILEUP			##
	##  just to make sure not a unique specie is inflating the variation 					##
	## 		        																		##
	##	Cristina Frias 05/06/2014 	cristinafriaslopez@ub.edu								##
	## 		        																		##
	##########################################################################################
	
	my $file_MID_array = $_[0];
	my $merge_file_coord = $_[1];
	my $sp_id = $_[2];
	my $output_file = $_[3];
	my $error_file = $_[4];

	### declare vars script
	my ($contig_name, $seq_MID);
	my $coverage_var_pct = 0.75;

	### open FILE MID
	#print "\tChecking $sp_id coordinates...\n";	
	my $hash_ref_coordinates = DOMINO::readFASTA_hash($$file_MID_array);
	my %hash_seq = %{$hash_ref_coordinates};

	##### open output
	open (OUT, ">>$$output_file") or &printError("Could not open output file $$output_file"); 
	open (ERR, ">>$$error_file") or &printError("Could not open error file $$error_file");

	### open coord and check coord of file 1
	open (COORD,"<$merge_file_coord") or &printError("Could not open merge coordenates $$merge_file_coord");
	while(<COORD>){
		chomp;
		my $line = $_;
		next if ($line =~ m/^Contig\tCons1_ini.*/o);
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;

		###################################################################
		###		 Get the coordinates of the different regions 			### 
		###################################################################

		my @regions= split(/\t/,$line);
		my $contig_MID = $regions[0];
		my $cons1 = $regions[1];
		my $var = $regions[2];
		my $cons2 = $regions[3];

		### split data of regions
		my ($cons1_START, $cons1_END) = split(/:/,$cons1);
		my ($VAR_START, $VAR_END) = split(/:/,$var);
		my ($cons2_START, $cons2_END) = split(/:/,$cons2);
		
		my ($coord_P1, $coord_P2) = split(/:/, $cons1);
		my ($coord_P3, $coord_P4) = split(/:/, $var);
		my ($coord_P5, $coord_P6) = split(/:/, $cons2);

		## Check if exists
		if (exists ($hash_seq{$contig_MID}) ){ $seq_MID = $hash_seq{$contig_MID};
		} else { print ERR "$contig_MID is not mapping for file $file_MID_array\n"; next; }

		###################################################################
		###			 Get the sequence of the different regions 			### 
		###################################################################
	
		### CONS1 ###
		my $size_P1_P2 = $coord_P2 - $coord_P1;
		my $string_P1_P2 = substr ($seq_MID, $coord_P1, $size_P1_P2);
		my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//;

		### VAR ###
		my $size_P3_P4 = $coord_P4 - $coord_P3;
		my $string_P3_P4 = substr ($seq_MID, $coord_P3, $size_P3_P4);
		my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//;

		### CONS2 ###
		my $size_P5_P6 = $coord_P6 - $coord_P5;
		my $string_P5_P6 = substr ($seq_MID, $coord_P5, $size_P1_P2);
		my $count_string_P5_P6 = $string_P5_P6 =~ tr/1//;

		###################################################################
		### 		Check the composition of the regions				###
		###################################################################
		
		## More than variations allowed in conserved
		if ($count_string_P1_P2 > $window_var_CONS) { next; } 
		if ($count_string_P5_P6 > $window_var_CONS) { next; }

		#Get marker profile					
		my $total_length = $coord_P6 - $coord_P1;
		my $marker_profile = $string_P1_P2.$string_P3_P4.$string_P5_P6;
		
		### Check variation
		my $expected_var_sites;	my $flag_error=0;
		if ($variable_divergence) {
			# If a minimun divergence, get the expected variable sites for the length
			$expected_var_sites = int($variable_divergence * $total_length);
			unless ($count_string_P3_P4 >= $expected_var_sites) { $flag_error++; }			
		} else {
			# User provides a number of variable positions or a range
			if ($expected_var_sites > $variable_positions_user_min) {
				unless ($expected_var_sites < $variable_positions_user_max) {
					$flag_error++;
				}} else { $flag_error++; }
		}
		
		## Missing % of bases missing
		my $missing_count = $marker_profile =~ tr/N//;
		my $missing_count2 = $marker_profile =~ tr/-//;
		$missing_count += $missing_count2;
		my $missing_count_percent = ($missing_count/$total_length)*100;
		my $percent_total_length = $total_length * $missing_allowed;  ## Default 0.05

		#if ($debugger) {
		#	print "\n\n***********************************************\n";
		#	print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
		#	print "Contig: $contig_MID\n\nPositions:\n";
		#	print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
		#	print "\nSubsets\n";
		#	print "Conserved (P1-P2): ".length($string_P1_P2)."\n"; print "$string_P1_P2\n";
		#	print "Variable (P3-P4): ".length($string_P3_P4)."\n";  print $string_P3_P4."\n";
		#	print "Conserved (P5-P6): ".length($string_P5_P6)."\n"; print $string_P5_P6."\n";
		#	print "\nVariations:\n";
		#	print "Count_string_P1_P2: $count_string_P1_P2\n";
		#	print "Count_string_P3_P4: $count_string_P3_P4\n";
		#	print "Count_string_P5_P6: $count_string_P5_P6\n";
		#	print "\nMarker:\n";
		#	print "Total Length:$total_length\tVariable Positions:$count_string_P3_P4\tExpected:$expected_var_sites\n";
		#	print "Missing allowed: $percent_total_length %\tMissing:$missing_count_percent %\n";
		#	print "***********************************************\n";
		#}

		if ($flag_error > 0) {	
			print ERR "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\t$sp_id\n"; next; 
		} else {
			### Print coordinates if the meet requirements
			if ($missing_count_percent < $percent_total_length) {
				print OUT "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\t$sp_id\n";
	}}} close(COORD); close(OUT); close(ERR);
}

sub get_headers_sam {
	##########################################################################################
	##	 																					##
	##  This function checks the SAM files generated and discards bad reads mapping,		##
	##	unmapping reads or multimapping reads.												##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	##	Cristina Frias Lopez 				  cristinafriaslopez@ub.edu 					##
	## 		        																		##
	##########################################################################################

	my $file_sam = $_[0];
	my @array;
	open (SAM, "<$file_sam");
	while (<SAM>) {
		chomp;
		my $line = $_;
		if ($line =~ /^@/ ) {
			push (@array, $line);
	}}
	close(SAM);
	return \@array;
}

sub get_parameters {

	my %parameters;
	my $array_files_ref=DOMINO::readDir($folder_abs_path);
	my @array_files = @$array_files_ref;
	my (%dirs, $earliest);
	for (my $i=0; $i<scalar @array_files;$i++) {
		if ($array_files[$i] eq "." || $array_files[$i] eq ".." || $array_files[$i] eq ".DS_Store") {next;}
		if ($array_files[$i] =~ /(\d+)\_DM\_(.*)/) {
			my $time_stamp=$1; my $type = $2;
			$dirs{$type}{$time_stamp} = $folder_abs_path."/".$array_files[$i];
	}}
	
	foreach my $dir (sort keys %dirs) {
		my $last;
		foreach my $times (sort {$a<=>$b} keys %{$dirs{$dir}}) {
			$last = $times;	## Only used the last folder for each process			
		}
		if ($dir eq "assembly") {
			if ($dirs{$dir}{$last}) {
				my $file = $dirs{$dir}{$last}."/".$last."_DM_Assembly_Parameters.txt";
				open (FILE, $file);
				while (<FILE>) {
					my $line = $_; chomp $line;
					if ($line =~ /(.*):\s(.*)\s\.\.\.OK/) {
						my $string = $1; my $var = $2;
						if ($string =~ /.*provided/) { 				$parameters{$dir}{"files"} = $var;
						} elsif ($string =~ /Threads.*/) {			$parameters{$dir}{"proc"} = $var;
						} elsif ($string =~ /.*Score.*/) {			$parameters{$dir}{"mrs"} = $var;
						}  elsif ($string =~ /Similarity.*/) { 		$parameters{$dir}{"simCAP3"} = $var;
						}  elsif ($string =~ /Overlapping.*/) { 	$parameters{$dir}{"overCAP3"} = $var; }
					} elsif ($line =~ /CAP3 has been disabled.*/) { $parameters{$dir}{"CAP3"} = "NO";
					} elsif ($line =~ /CAP3 would be used.*/) {		$parameters{$dir}{"CAP3"} = "YES";
					} elsif ($line =~ /Starting the process:\s+\[(.*)\]/) {	$parameters{$dir}{"start"} = $1;
				}} close (FILE); }
				$parameters{$dir}{"folder"} = $dirs{$dir}{$last};
				
				my $files_assembly_ref = DOMINO::readDir($dirs{$dir}{$last});
				my @files_assembly = @$files_assembly_ref;
				for (my $j=0; $j < scalar @files_assembly; $j++) {
					if ($files_assembly[$j] eq "." || $files_assembly[$j] eq ".." || $files_assembly[$j] eq ".DS_Store") { next; }
					if ($files_assembly[$j] =~ /.*id\-(.*)\.contigs-statistics.txt/) {
						$parameters{$dir}{"stats"}{$1} = $dirs{$dir}{$last}."/".$files_assembly[$j];
				}}
		} elsif ($dir eq "clean_data") {
			if ($dirs{$dir}{$last}) {
				my $file = $dirs{$dir}{$last}."/".$last."_DM_Cleaning_Parameters.txt";
				open (FILE, $file);
				while (<FILE>) {
					my $line = $_; chomp $line;
					if ($line =~ /.*Checking file\(s\)\:/) {
						my %file_seqs;
						while (1) {
							my $new_line = <FILE>;
							if ($new_line =~ /\+.*/) { last; }
							if ($new_line =~ /\s+(\S+)/) {
								my $file_path = $1;
								my $counter_while = 0;
								while (1) {
									my $new_line2 = <FILE>;
									if (!$new_line2) {last;}
									last if $new_line2 =~ /^\s*$/;
									if ($new_line2 =~ /.*identifier.*\((\S+)\).*/) {
										$parameters{$dir}{"clean_files"}{$1} = $file_path;
										last; $counter_while++; if ($counter_while == 50) {last;}
					}}}}
					} elsif ($line =~ /Starting the process:\s+\[(.*)\]/) {	$parameters{$dir}{"start"} = $1;
					} elsif ($line =~ /(.*):\s(.*)\s\.\.\.OK/) {
						my $string = $1; my $var = $2;
						if ($string =~ /Type.*/) { 			 	$parameters{$dir}{"option"} = $var;
						} elsif ($string =~ /.*mismatches/) { 	$parameters{$dir}{"mismatches"} = $var;
						} elsif ($string =~ /.*read length.*/){	$parameters{$dir}{"read_length"} = $var;
						} elsif ($string =~ /Minimum QUAL.*/){	$parameters{$dir}{"QUAL"} = $var;
						} elsif ($string =~ /.*cutoff.*/) {		$parameters{$dir}{"QUAL_cutoff"} = $var;
						} elsif ($string =~ /.*entropy.*/) {	$parameters{$dir}{"entropy"} = $var;}
					} elsif ($line =~ /.*Database\(s\).*/) {
						my @db_seqs;
						my $counter_while = 0;
						while (1) {
							my $new_line = <FILE>;
							if ($new_line =~ /\s+\-\s+(\S+)/) {
								my @array = split("/", $1);
								push (@db_seqs, $array[-1]);
							} else { last; }
							$counter_while++; if ($counter_while == 50) {last;}
							}
						push (@{$parameters{$dir}{"dbs"}}, @db_seqs);
				}} close (FILE);				
			$parameters{$dir}{"folder"} = $dirs{$dir}{$last};
	
			# Retrieved Quality filtering results
			my $int_data = $dirs{$dir}{$last}."/".$last."_DM_cleaning_intermediate_data";
			my $ref = DOMINO::readDir($int_data);
			my @array_qc_stats = @$ref;
			for (my $i=0; $i < scalar @array_qc_stats; $i++) {
				if ($array_qc_stats[$i] eq "." || $array_qc_stats[$i] eq ".." || $array_qc_stats[$i] eq ".DS_Store") {next;}
				if ($array_qc_stats[$i] =~ /Quality\_Filtering\_statistics\_(.*)\.txt/) {
					my $QC_file = $int_data."/".$array_qc_stats[$i];
					$parameters{$dir}{"QC_analysis"}{$1} = $QC_file;
	}}}}}
	if (!%parameters) { return 0;	
	} else { return \%parameters;}
}

sub generate_bam {
	my $sam_file = $_[0];
	my $avoid = $_[1];
	my @temp = split ("\.sam", $sam_file);
	my $name = $temp[0]; my $bam_file = $name.".bam";
	&debugger_print("- Generating a BAM file for $sam_file\n"); 	
	my $system_samtools_sam2bam = $samtools_path." view -@ $num_proc_user -bS -o $bam_file $sam_file";
	&debugger_print("SAMTOOLS command: $system_samtools_sam2bam");	
	my $system_call = system ($system_samtools_sam2bam);
	if ($system_call != 0) {
		if (!$avoid) { &printError("Some error happened when calling SAMTOOLs for SAM -> BAM conversion $sam_file to $bam_file...."); DOMINO::dieNicely(); }
	}
	my $sorted;
	if ($avoid) { if ($avoid eq 'sam') { $sorted = &generate_sorted_bam($bam_file, "sam"); } }
	if (!$avoid) {$sorted = &generate_sorted_bam($bam_file);}
	return $sorted;
}

sub generate_sorted_bam {
	my $bam_file = $_[0];
	my $sam = $_[1];
	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	&debugger_print("- Sorting the BAM file: $bam_file\n"); 	
	my $sorted;	
	my $system_samtools_sort;
	if ($sam) {
		$sorted = $name.".sorted.sam";	
		$system_samtools_sort = $samtools_path." sort -@ $num_proc_user --output-fmt SAM -o $sorted $bam_file";
	} else {
		$sorted = $name.".sorted.bam";	
		$system_samtools_sort = $samtools_path." sort -@ $num_proc_user -o $sorted $bam_file";
	}	 
	&debugger_print("SAMTOOLS command: ".$system_samtools_sort);
	my $system_call_2 = system ($system_samtools_sort);
	if ($system_call_2 != 0) {
		&printError("Some error happened when calling SAMTOOLs for sorting BAM file...."); DOMINO::dieNicely();
	}
	return $sorted;
}

sub generate_sam {
	my $bam_file = $_[0];
	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	&debugger_print("- Generating a SAM file for $bam_file\n"); 	
	my $system_samtools_bam2sam = $samtools_path." view -@ $num_proc_user -h ".$bam_file." -o ".$name.".sam";
	&debugger_print("SAMTOOLS command: $system_samtools_bam2sam");
	my $system_call = system ($system_samtools_bam2sam);
	if ($system_call != 0) {
		&printError("Some error happened when calling SAMTOOLs for BAM -> SAM conversion....");  DOMINO::dieNicely();
	}
	return $name.".sam";
}

sub generate_index_bam {
	my $bam_file = $_[0];
	#print "\t- Generating an index bam file for $bam_file\n"; 	
	my $index_system = $samtools_path." index ".$bam_file;
	&debugger_print("SAMTOOLS command: ".$index_system);
	my $index_call = system($index_system);
	if ($index_call != 0) {
		&printError("Some error happened when calling SAMTOOLs for indexing the BAM file [$bam_file]...."); DOMINO::dieNicely();
	}
}

sub generate_filter_PILEUP {

	##########################################################################################
	##  This function generates a PILEUP and filters it										##
	##	Jose Fco. Sanchez Herrero, 08/06/2015 jfsanchezherrero@ub.edu						##
	##########################################################################################
	
	my $sorted_bam = $_[0]; my $contig_file = $_[1]; 
	my $reference_hash_fasta = $_[2]; my $reference_id = $_[3]; my $taxa = $_[4];

	my $dir_path = $domino_files{$reference_id}{'dir'}[0];
	my @temp_name = split ("\.sorted.bam", $sorted_bam);
	my ($ID, @sam);
	my $input_pileup = $temp_name[0].".profile";
	#push (@{ $domino_files{$reference_id}{'PILEUP_'.$taxa} }, $input_pileup);
	
	my $pileup_command = $samtools_path." mpileup -f ".$contig_file." -o ".$input_pileup." ".$sorted_bam." 2> ".$mapping_markers_errors_details;
	&debugger_print("SAMTOOLS command: ".$pileup_command);
	my $sytem_command_pileup = system ($pileup_command);
	if ($sytem_command_pileup != 0) {
		&printError("Exiting the script. Some error happened when calling SAMtools for generating the PILEUP for the file $contig_file...\n"); DOMINO::dieNicely();
	}
	my $tmp = $dir_path."/ARRAY_files_".$taxa."_PROFILE"; unless (-d $tmp) { mkdir $tmp, 0755; } 
	&debugger_print("Changing dir to $tmp");

  	#print "\t- Filtering the PILEUP generated\n";
	my ($previous_contig, $previous_fasta_contig, @array_positions, @fasta_positions);
	open (PILEUP,"<$input_pileup"); while (<PILEUP>){
		my $line = $_; chomp $line;
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;
		next if $line=~ m/\[REDUCE RESULT\]/;
		$line =~ s/\s/\t/g;
		my @pileup = split(/\t/,$line); ##	HKU61D301.clean.MID4_Nem_098_c8679	161	t	3	,,.	FA=
		my $contig = $pileup[0];
		my $pos_base; my $num_pos_base = $pileup[1]; my $num_pos_array = $num_pos_base -1;

		if (!$previous_contig) {
			my $array_positions_ref;
			($previous_contig, $array_positions_ref) = &initilize_contig($contig, $reference_hash_fasta);
			@array_positions = @$array_positions_ref;
			@fasta_positions = @array_positions;			
		} else {			
			if ($previous_contig ne $contig) {
				# Debug print Dumper \@array_positions;
				&print_coordinates(\@array_positions, \$previous_contig, $reference_id, $tmp); ## Print array into file $previous_contig
				# Debug	print Dumper \@fasta_positions;
				&print_fasta_coordinates(\@fasta_positions, \$previous_contig, $reference_id, $tmp); ## Print array into file $previous_contig
				my $array_positions_ref;
				($previous_contig, $array_positions_ref) = &initilize_contig($contig, $reference_hash_fasta);
				@array_positions = @$array_positions_ref;
				@fasta_positions = @array_positions;
		}}
		
		## Get array for each position
		if ($pileup[3] != 0) {
			my $read_base = $pileup[4]; my $ref_base = $pileup[2];
			my @base_record = split(//, $read_base);
			my (%posibilities, %polymorphism);
			my @base_parse;
			my $base_counter=0;
			for (my $i = 0; $i < scalar @base_record; $i++) {  
				if ($base_record[$i] =~ m/\^/) { $i++; next; } ## Starting to map a new read
				$base_record[$i]= uc($base_record[$i]);

				if ($base_record[$i] =~ m/A|G|C|T|a|g|c|t/) {
					$polymorphism{$base_record[$i]}++;
					$base_counter++;
					push (@base_parse, $base_record[$i]);
				} elsif ($base_record[$i] eq "." || $base_record[$i] eq ",") {
					$polymorphism{$ref_base}++;
					$base_counter++;
					push (@base_parse, $ref_base);
				} elsif ($base_record[$i] =~ m/(\+|\-)/) {# INDEL
					##	Example of INDEL: 
					##	HKU61D301.clean.MID4_Nem_098_c8679	280	g	2	,.+1T	I? 
					##  Contig_47_MID1	4	T	1	a+47cggatcgatcaaagtaagatatcatacttggaaggcaacatgcacgt	1	1	Position: 1
					my $length_base_record = length($read_base);
					my $indel;
					my $out_of_range = 0;
					my $indel_2 = $i+2;
					my $indel_3 = $i+3;
					if ($indel_2 >= $length_base_record) { $out_of_range = 1; }
					if ($indel_3 >= $length_base_record) { $out_of_range = 1; }
					if ($out_of_range == 0) {
						if ($base_record[$i+2] =~ /\d+/) { ## +/-12ATGAGACGATCC
							$indel = $base_record[$i+1].$base_record[$i+2];
							$i++; $i++;
						} elsif ($base_record[$i+3] =~ /\d+/) { ## +/-121ATGAGACGATCC...ATGAGACGATCC
							$indel = $base_record[$i+1].$base_record[$i+2].$base_record[$i+2];
							$i++; $i++;	$i++;
						} else { ## +/-1A
							$indel = $base_record[$i+1];
							$i++;
					}} else { ## +/-1A
						$indel = $base_record[$i+1];
						$i++;
					}
					$i = $i+$indel;
					next;
				} else { next; }
			}
			# Debug			
			#print Dumper @base_record; #print Dumper @base_parse; print Dumper %polymorphism;
			
			if ($ref_base eq "N") {
				$array_positions[$num_pos_array] = 'N'; ## Not informative enough  
				my @array_keys = keys %polymorphism;
				my $position;
				if (scalar @array_keys == 1) { ## a unique base is mapping
					$position = $array_keys[0];
				} else { ## get ambiguous code
					$position = &get_amb_code(\%polymorphism);
				}
				$fasta_positions[$num_pos_array] = $position; 
				# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
				next;	
			}
			
			if ($base_counter == 1) { ## There is only base mapping...
				unless ($map_contig_files || $DOMINO_simulations) {
					$array_positions[$num_pos_array] = 'N'; ## Not informative enough
					$fasta_positions[$num_pos_array] = 'N'; ## Not informative enough
					# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
					next;
			}} ## Let it be informative if DOMINO simulations or mapping contig files				
			
			my @array_keys = keys %polymorphism;
			if (scalar @array_keys == 1) { ## a unique base is mapping
				if ($ambiguity_DNA_codes{$ref_base}) {
					my $flag = 0; my @bases;
					for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$ref_base}}; $j++) {
						foreach my $keys (sort keys %polymorphism) {
							if ($keys eq $ambiguity_DNA_codes{$ref_base}[$j]) { $flag = 1; }
					}}
					if ($flag == 1) { 
						## Contig1	4	M	10	cccccccccc	## type 4
						$array_positions[$num_pos_array] = '0'; 
						$fasta_positions[$num_pos_array] = $ref_base;
					} else { 
						## Contig1	5	R	10	cccccccccc	## type 5
						$array_positions[$num_pos_array] = '1';
						$fasta_positions[$num_pos_array] = $array_keys[0];											
				}} else {
					if ($base_parse[0] eq $ref_base) {  
						##Contig1	2	T	10	..........	## type 2
						$array_positions[$num_pos_array] = '0';
						$fasta_positions[$num_pos_array] = $ref_base;						
					} else {  
						## Contig1	3	T	10	cccccccccc	## type 3
						$array_positions[$num_pos_array] = '1';
					    $fasta_positions[$num_pos_array] = $base_parse[0];
			}}} elsif (scalar @array_keys == 2) { ## Maybe true polymorphism or NGS error
				my ($value_return, $base_return) = &check_array($ref_base, \%polymorphism);
				$array_positions[$num_pos_array] = $value_return; 
			    $fasta_positions[$num_pos_array] = $base_return;
			} elsif (scalar @array_keys == 3) { ## Something odd...
				my ($last_key, $small_key, $highest_value, $smallest_value);
				my $flag=0;				
				my @array = sort values %polymorphism;
				for (my $i=1; $i < scalar @array; $i++) {
					if ($array[0] == $array[$i]) { $flag++; }
				}
				if ($flag != 0) { ## There are two with the same frequence! DISCARD!
					$array_positions[$num_pos_array] = 'N';
					$fasta_positions[$num_pos_array] = 'N';
				} else {
					## Get lowest value, discard it and check the rest.
					foreach my $keys (sort keys %polymorphism) {
						if ($polymorphism{$keys} eq $array[0]) {
							delete $polymorphism{$keys}; last;
					}}
					my ($value_return, $base_return) = &check_array($ref_base, \%polymorphism);
					$array_positions[$num_pos_array] = $value_return; 
			    	$fasta_positions[$num_pos_array] = $base_return;
			}} elsif (scalar @array_keys > 3) { 
				$array_positions[$num_pos_array] = 'N';
				$fasta_positions[$num_pos_array] = 'N';
		}} else { ## No base is mapping this reference
			##Contig1	1	A	0	## type 1
			$array_positions[$num_pos_array] = 'N';
			$fasta_positions[$num_pos_array] = 'N';
		}
		# Debug	print $array_positions[$num_pos_array]."\n"; print $fasta_positions[$num_pos_array]."\n";
	}
	close(PILEUP);
	&print_coordinates(\@array_positions, \$previous_contig, $reference_id, $tmp);
	&print_fasta_coordinates(\@fasta_positions, \$previous_contig, $reference_id, $tmp); ## Print array into file $previous_contig
	chdir $dir_path; &debugger_print("Changing dir to $dir_path");
	return ($tmp);
	
	sub check_array {
		my $ref_base = $_[0];
		my $ref_poly_hash = $_[1];
		my %polymorphism = %$ref_poly_hash;
		my $total=0;
		my ($last_key, $highest_value, $smallest_value);
		foreach my $keys (sort keys %polymorphism) {
			$total += $polymorphism{$keys};
			if (!$highest_value) {
				$highest_value = $polymorphism{$keys};
				$smallest_value = $polymorphism{$keys};
				$last_key = $keys;							
			} elsif ($polymorphism{$keys} > $highest_value) {
				$highest_value = $polymorphism{$keys};
				$last_key = $keys;							
			} elsif ($polymorphism{$keys} < $smallest_value) {
				$smallest_value = $polymorphism{$keys};
		}}		
		
		if ($highest_value >= 170) { return ("N","N");}
		## Check wether there are more than 8 positions mapping
		## and if not if there at least for the smallest two
		## bases.
		
		my $prob = &binomial($highest_value, $total, 0.5);
		my $significance_level = 1 - 0.95;
		if ($ambiguity_DNA_codes{$ref_base}) {
			## There is an ambiguous code in the reference
			my %ref_base_match;
			my $flag = 0;
			for (my $j = 0; $j < scalar @{ $ambiguity_DNA_codes{$ref_base}}; $j++) {
				foreach my $keys (sort keys %polymorphism) {
					if ($keys eq $ambiguity_DNA_codes{$ref_base}[$j]) {
						$flag++; $ref_base_match{$ambiguity_DNA_codes{$ref_base}[$j]}++;
			}}}
			if ($flag == 2) { 
				## Only if both possibilities are the same as the ambiguous code
				if ($prob < $significance_level) { 
					return ('0', $last_key); ##Contig1	10	R	10	AAAAAAAAAG
				} else { 
					## Contig1	9	R	10	AAAAAAGGGG
					## Contig1	9	R	10	AAGG
					## Contig1	9	R	10	AG
					if ($polymorphism_user) {
						return ('1', $ref_base);
					} else {
						return ('0', $ref_base);
			}}} else {
				if ($flag == 1) { ## Only one type of read is within the ambiguous code
					if ($prob < $significance_level) { 
						if ($ref_base_match{$last_key}) { 
							return ('0', $last_key); ## Contig1	11	R	10	AAAAAAAAAC
						} else { 
							return ('1', $last_key); ## Contig1	13	R	10	TTTTTTTTTTA
					}} else { 
						## Contig1	12	R	10	AAAAACCCCC
						## Contig1	12	R	10	AACC
						## Contig1	12	R	10	AC
						my $pos = &get_amb_code(\%polymorphism);
						if ($polymorphism_user) {							
							return ('1', $pos);
						} else {
							return ('0', $pos);
				}}} else { ## None of the reads mapping are similar to the ambiguous reference
					if ($prob < $significance_level) { 	
						return ('1', $last_key); ## Contig1	13	R	10	TTTTTTTTTTC
					} else { 
						## Contig1	13	R	10	TTTTTTCCCCC
						## Contig1	13	R	10	TTCC
						## Contig1	13	R	10	TC
						my $pos = &get_amb_code(\%polymorphism);
						return ('1', $pos);
		}}}} else { 
			## There is no ambiguity code in the reference base
			my $flag = 0;
			if ($polymorphism{$ref_base}) { ## if any mapping reads are the same as reference
				if ($total >= 8) {
					if ($prob < $significance_level) {
						if ($last_key eq $ref_base) { 
							return ('0', $last_key); ## Contig1	6	T	10	.........a
						} else { 
							return ('1', $last_key); ## Contig1	7	T	10	ccccccccc.
					}} else { 
						## Contig1	8	T	10	....cc..cc
						my $pos = &get_amb_code(\%polymorphism);
						if ($polymorphism_user) {
							return ('1', $pos);
						} else {
							return ('0', $pos);
				}}} else { ## less than 8 reads mapping here
					if ($smallest_value >= 2) {
						## Contig1	8	T	5	...cc
						my $pos = &get_amb_code(\%polymorphism);
						if ($polymorphism_user) {
							return ('1', $pos);
						} else {
							return ('0', $pos);
					}} else {
						if ($highest_value == $smallest_value) {
							return ('N', 'N'); ## Contig1	8	T	2	.c
						} elsif ($last_key eq $ref_base) { 
							return ('0', $ref_base); ## Contig1	8	T	4	...c
						} else {
							return ('1', $last_key); ## Contig1	8	T	4	.c
			}}}} else {
				## None of the reads are similar to reference
				if ($total >= 8) { 
					if ($prob < $significance_level) { 	
						return ('1', $last_key); ## Contig1	15	T	10	ccccccccccca
					} else { 
						my $pos = &get_amb_code(\%polymorphism);
						return ('1', $pos); ## Contig1	15	T	10	ccccccgggggg
				}} else {
					if ($smallest_value >= 2) {
						my $pos = &get_amb_code(\%polymorphism);
						return ('1', $pos); ## Contig1	8	T	5	GGGcc
					} else {
						if ($highest_value == $smallest_value) {
							return ('1', 'N'); ## Contig1	8	T	2	Gc
						} else {
							return ('1', $last_key); ## Contig1	8	T	3	GGc
	}}}}}}	

	sub initilize_contig {		
		my $current_contig = $_[0];
		my $reference_hash_fasta = $_[1];
		
		## Initialize new contig
		my @array_positions_sub;
		if (${$reference_hash_fasta}{$current_contig}) {
			my $tmp_size = ${$reference_hash_fasta}{$current_contig};
			@array_positions_sub = ("-") x $tmp_size;
		}
		my $ref = \@array_positions_sub;		
		return ($current_contig, $ref);
	}
	
	sub print_coordinates {
		## Print array into file $previous_contig
		my $coord_array_ref = $_[0];
		my $contig_name = $_[1];
		my $ref_id = $_[2];
		my $dir_tmp = $_[3];

		my @coord_array = @$coord_array_ref;
		my $seq_contig = join "", @coord_array;
		my $var_sites = $seq_contig =~ tr/1/1/;

		if ($var_sites != 0) {
			my $array_file = $dir_tmp."/".$$contig_name."_ARRAY.txt";
			open (FH, ">$array_file"); print FH ">$$contig_name\n$seq_contig\n"; close(FH);	
		}
	}
	
	sub print_fasta_coordinates {
		## Print array into file $previous_contig
		my $coord_array_ref = $_[0];
		my $contig_name = $_[1];
		my $ref_id = $_[2];
		my $dir_tmp = $_[3];
		
		my @coord_array_sub = @$coord_array_ref;
		my $seq_contig = join "", @coord_array_sub;
		my $array_file = $dir_tmp."/".$$contig_name."_sequence.fasta";
		open (FH, ">$array_file"); print FH ">$$contig_name\n$seq_contig\n"; close(FH);	
	}
}

sub Poisson_distribution {
	my ($x, $a) = @_;
	return unless $a >= 0 && $x >= 0 && $x == int($x); 
	return (($a ** $x) * exp(-$a))/&functional_factorial($x);	
}

sub printError {
    my $msg = $_[0];
	print "\n\n";DOMINO::printHeader(" ERROR ","!!"); print "\n";
    print $msg."\n\nTry \'perl $0 -h|--help or -man\' for more information.\nExit program.\n";
	print "\n\n"; DOMINO::printHeader("","!!"); DOMINO::printHeader("","!!"); 
    DOMINO::printError_log($msg, $mapping_markers_errors_details);
}

sub print_Excel {

	my $markers_file = $_[0]; my $path = $_[1];
	my @array_markers;
	open (FILE, "$$markers_file");
	while (<FILE>) {
		chomp; if ($_ =~ /.*Vari.*/) {next;}
		push (@array_markers, $_);
	} close (FILE);
	my $hash_parameters = &get_parameters();
	my $no_parameters; if ($hash_parameters == 0) { $no_parameters = 1; }
	
	#################################################################################
	##	Once the coordinates are found, print different files with the information ##
	#################################################################################	
	### open Output and Error file
	my $excel_woorkbook_name = $$path."/DM_markers-summary.xls";
	print "+ Printing putative Molecular markers in the XLS format file $excel_woorkbook_name...\n";
	
	## Generate an Excel file to write results
	my $workbook = Spreadsheet::WriteExcel->new($excel_woorkbook_name);
	my $worksheet_parameters = $workbook->add_worksheet("Parameters - Input");
    my $col = my $first_col = 0; 
    my $row = my $second_col = 1;

    my $format_main_heading = $workbook->add_format();
    $format_main_heading->set_size(14);
    $format_main_heading->set_font('Arial Black');

	my $format = $workbook->add_format(); # Add a format
	$format->set_size(12);
	
	my $format_left = $workbook->add_format(); # Add a format
	$format_left->set_size(12);
	$format_left->set_align('left');
	
	my $format_right = $workbook->add_format(); # Add a format
	$format_right->set_size(12);
	$format_right->set_align('right');
	
	my $format_bold = $workbook->add_format(); # Add a format
	$format_bold->set_size(12);
	$format_bold->set_font('Arial Black');
	$format_bold->set_align('left');

   	## Write Excel Heading and parameters
    $worksheet_parameters->write($row, $col, "DOMINO version:", $format_main_heading); $col++;
    $worksheet_parameters->write($row, $col, $domino_version, $format); $col++;
    #$worksheet_parameters->write($row, $col, $date, $format_left); 
    $row++; $row++;
    $worksheet_parameters->write($row, $first_col, "PROJECT:", $format_main_heading); 
    $col = $first_col; $row++; $row++;
    $worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
    $worksheet_parameters->write($row, $col, $folder_abs_path, $format); $row++;
    $col = $first_col;
    $worksheet_parameters->write($row, $col, "OPTION:", $format); $col++;
    if ($no_parameters) {
    	if ($radseq_like_data) {
			$worksheet_parameters->write($row, $col, "RADseq", $format); $row++;
		} else { $worksheet_parameters->write($row, $col, $option, $format); $row++; }
    } else {
    	if ($$hash_parameters{'clean_data'}{'option'}) {
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'option'}, $format); $row++;
    }} $row++; $row++;

	### Writing markers coordinates
	my $worksheet_markers = $workbook->add_worksheet("Output - Markers");
    my $col_markers = $first_col = 0; my $row_markers = $second_col = 1;
    $worksheet_markers->write($row_markers, $col_markers, "OUTPUT", $format_main_heading); $row_markers++; $col_markers++; $row_markers++;
    
    my $markers = 0;
    if ($select_markers) {
		$col_markers = $second_col;
		$worksheet_markers->write($row_markers, $col_markers, "Region ID", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Taxa included", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Variable sites", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Effective length", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Variation", $format_bold); $row_markers++;
    	
    	for (my $i = 0; $i < scalar @array_markers; $i++) {
			if ($array_markers[$i] eq "undef") { $row_markers++; next;}
			my @split = split("\t", $array_markers[$i]);
			$col_markers = $second_col;
			$worksheet_markers->write($row_markers, $col_markers, $split[2], $format_left); $col_markers++;	# Contig name
			$worksheet_markers->write($row_markers, $col_markers, $split[3], $format_left); $col_markers++; 	## taxa names
			$worksheet_markers->write($row_markers, $col_markers, $split[4], $format_right); $col_markers++;	## Variable sites
			$worksheet_markers->write($row_markers, $col_markers, $split[5], $format_right); $col_markers++;	## effective length
			$worksheet_markers->write($row_markers, $col_markers, $split[6], $format_right); $row_markers++;	## Variation percentage
			$markers++;
		}
    } else {
		$worksheet_markers->write($row_markers, $col_markers, "Conserved region: Left", $format_bold); $col_markers++; 
		$worksheet_markers->write($row_markers, $col_markers, "Variable region", $format_bold); $col_markers++; 
		$worksheet_markers->write($row_markers, $col_markers, "Conserved region: Right", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Taxa included", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Marker length", $format_bold);  $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Variation", $format_bold);  $row_markers++;
		$col_markers = $first_col;
		if ($option eq "msa_alignment") { 
			$worksheet_markers->write($row_markers, $col_markers, "Multiple alignment ID", $format_bold); $col_markers++;
		} else {
			$worksheet_markers->write($row_markers, $col_markers, "Reference contigs ID", $format_bold); $col_markers++;
		}    
		$worksheet_markers->write($row_markers, $col_markers, "Start - End", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Start - End", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Start - End", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "IDs", $format_bold); $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "Length (bp)", $format_bold);  $col_markers++;
		$worksheet_markers->write($row_markers, $col_markers, "% Variable positions", $format_bold);  $row_markers++;
	
		for (my $i = 0; $i < scalar @array_markers; $i++) {
			if ($array_markers[$i] eq "") { $row_markers++; next;}
			
			my @split = split("\t", $array_markers[$i]);
			$worksheet_markers->write($row_markers, $first_col, $split[0], $format); # Contig name
			my $position = $first_col + 1;
			$worksheet_markers->write($row_markers, $position, $split[1], $format_right);	$position++; ## Conserved left
			$worksheet_markers->write($row_markers, $position, $split[2], $format_right); $position++; ## Variable
			$worksheet_markers->write($row_markers, $position, $split[3], $format_right); $position++; ## Conserved Right
			$worksheet_markers->write($row_markers, $position, $split[6], $format_left); $position++; ## taxa names
			$worksheet_markers->write($row_markers, $position, $split[4], $format_right); $position++; ## variable region length
			$worksheet_markers->write($row_markers, $position, $split[5], $format_right); $row_markers++;	 	## Divergence
			$markers++;
	}} 

	## Print results and instructions
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" RESULTS ","#"); DOMINO::printHeader("", "#");
	print "\n+ DOMINO has retrieved $markers markers\n";
	my $instructions_txt = $marker_dirname."/Instructions.txt";
	my $MSA = $$path."/MSA_markers";
	open (OUT, ">$instructions_txt");
	my $string = "+ Several files and folders has been generated:
\t+ $marker_dirname: contains DOMINO markers detected for each taxa as a reference and the clusterized results.
\t+ $$path: contains the clusterized and definitive results for DOMINO markers.
\t+ $excel_woorkbook_name: contains information of the markers identified and parameters used by DOMINO.
\t+ $MSA folder contains a single file for each marker identified.\n\n";	
	print OUT $string; print $string."\n"; close(OUT);

	unless ($no_parameters) {
		if ($$hash_parameters{'clean_data'}) {
			$worksheet_parameters->write($row, $first_col, "PRE-PROCESSING PHASE:", $format_main_heading); 
			$col = $first_col; $row++; $row++;
			$worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'start'}, $format_left); $row++; $row++;
			$col = $first_col; 
			$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'folder'}, $format); $row++; $row++;
			$col = $first_col;
			$worksheet_parameters->write($row, $col, "INPUT DATA:", $format); $col++;
			$worksheet_parameters->write($row, $col, "Original Files included", $format_bold); $col++;
			$worksheet_parameters->write($row, $col, "DOMINO Label", $format_bold); $col++;
			foreach my $keys (sort keys %{$$hash_parameters{'clean_data'}{'clean_files'}}) {
				$col = $first_col; $col++; $row++;
				$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'clean_files'}{$keys}, $format_left);
				$col++; $worksheet_parameters->write($row, $col, $keys, $format_left);
			}
			$row++; $row++; $row++;
			$col = $first_col;
			$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col++; 
			$worksheet_parameters->write($row, $col, "Mismatches allowed in the barcode sequence:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'mismatches'}, $format_right); 
			$col = $first_col + 1; $row++;
			$worksheet_parameters->write($row, $col, "Minimum read length:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'read_length'}, $format_right); 
			$col = $first_col + 1; $row++;
			$worksheet_parameters->write($row, $col, "Minimum QUAL:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'QUAL'}, $format_right); 
			$col = $first_col + 1; $row++;
			$worksheet_parameters->write($row, $col, "Minimum length of a read satisfying QUAL cutoff:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'QUAL_cutoff'}, $format_right); 
			$col = $first_col + 1; $row++;
			$worksheet_parameters->write($row, $col, "Threshold for complexity/entropy of a read:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'clean_data'}{'entropy'}, $format_right); 
			$col = $first_col + 1; $row++; $row++;
			$worksheet_parameters->write($row, $col, "Database(s)", $format_bold); $row++;
			my @array_db = @{$$hash_parameters{'clean_data'}{'dbs'}};
			for (my $i=0; $i < scalar @array_db; $i++) {
				$worksheet_parameters->write($row, $col, $array_db[$i], $format_left); $row++;
			} $row++; $row++;
		
			if ($$hash_parameters{'clean_data'}{"QC_analysis"}) {
				my %QC_hash = %{$$hash_parameters{'clean_data'}{'QC_analysis'}};
				foreach my $taxa (sort keys %QC_hash) {
					my $worksheet_name = "Quality Filtering Stats $taxa";
					if (length($worksheet_name) >= 31) {
						$worksheet_name = "QC Stats $taxa";
						if (length($worksheet_name) >= 31) {
							$worksheet_name = "QC $taxa";
							if (length($worksheet_name) >= 31) {
								$worksheet_name = "$taxa";
					}}}
					
					my $worksheet_stats_QC = $workbook->add_worksheet($worksheet_name);
					my $col_stats_QC = 1; my $row_stats_QC = 2;
					my $file = $QC_hash{$taxa};
					open (FILE, "<$file");
					while (<FILE>) {
						my $line = $_;
						chomp $line;
						$col_stats_QC = 1;
						if ($line =~ /Input files/) {next;}
						if ($line =~ /File name/) {next;}
						my @array = split("\t", $line);
						for (my $i=0; $i < scalar @array; $i++) {
							if ($i > 1) { $array[$i] =~ s/\s//g;}
							if ($array[$i] =~ /Detailed QC statistics/) {
								$worksheet_stats_QC->write($row_stats_QC, $col_stats_QC, $array[$i], $format_left);
								$row_stats_QC++; $col_stats_QC++; $col_stats_QC++;
								$worksheet_stats_QC->write($row_stats_QC, $col_stats_QC, "Original", $format_left);
								$col_stats_QC++; $worksheet_stats_QC->write($row_stats_QC, $col_stats_QC, "Clean", $format_left);
							} else {
								$worksheet_stats_QC->write($row_stats_QC, $col_stats_QC, $array[$i], $format_left);
								$col_stats_QC++;
						}} $row_stats_QC++;
						} close (FILE); 
		}}}
	
		if ($$hash_parameters{'assembly'}) {
			$worksheet_parameters->write($row, $first_col, "ASSEMBLY PHASE:", $format_main_heading); 
			$row++; $row++;	$col = $first_col; 
			$worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'assembly'}{'start'}, $format_left); 
			$row++; $row++; $col = $first_col; 
			$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'assembly'}{'folder'}, $format_left); 
			$row++; $row++; $col = $first_col;
			if ($$hash_parameters{'assembly'}{'files'} eq "Default DOMINO cleaning files") {
				$worksheet_parameters->write($row, $col, "INPUT DATA", $format); $col++;
				$worksheet_parameters->write($row, $col, "Default DOMINO cleaned files", $format_left); $row++; $col++;
			} else { 
				## Control if user provides different files
			}			
			$row++; $row++; $col = $first_col;
			$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col = $first_col + 1; 
			$worksheet_parameters->write($row, $col, "Minimum read score (MIRA):", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters{'assembly'}{'mrs'}, $format_left); $col++;
			
			if ($$hash_parameters{'assembly'}{'CAP3'} eq 'YES') {
				$col = $second_col; $row++;
				$worksheet_parameters->write($row, $col, "CAP3:", $format); $col++;
				$worksheet_parameters->write($row, $col, "Enabled", $format_left); $col++;			
				$col = $second_col; $row++;
				$worksheet_parameters->write($row, $col, "Similarity:", $format); $col++;
				$worksheet_parameters->write($row, $col, $$hash_parameters{'assembly'}{'simCAP3'}, $format_left); $col++;
				$col = $second_col; $row++;
				$worksheet_parameters->write($row, $col, "Overlapping:", $format); $col++;
				$worksheet_parameters->write($row, $col, $$hash_parameters{'assembly'}{'overCAP3'}, $format_left); $col++;
			} else {
				$col = $second_col; $row++;
				$worksheet_parameters->write($row, $col, "CAP3: Disabled", $format);
			} $row++; $row++; $row++; 
			
			if ($$hash_parameters{'assembly'}{'stats'}) {
				my %hash = %{$$hash_parameters{'assembly'}{'stats'}};
				foreach my $taxa (sort keys %hash) {
					my $worksheet_name = "Assembly Stats $taxa";
					if (length($worksheet_name) >= 31) {
						$worksheet_name = "Ass. Stats $taxa";
						if (length($worksheet_name) >= 31) {
							$worksheet_name = "Stats $taxa";
							if (length($worksheet_name) >= 31) {
								$worksheet_name = "$taxa";
					}}}
					
					my $worksheet_stats = $workbook->add_worksheet($worksheet_name);
					my $col_stats = 1; my $row_stats = 2;
					open (STATS, $hash{$taxa});
					while (<STATS>) {
						my $line = $_;
						chomp $line;
						if ($line =~ /.*-.*/) { $row_stats++; next; }
						my @array = split("\:", $line);
						$col_stats = 1;
						for (my $h=0; $h < scalar @array; $h++) {
							if ($h > 0) { $array[$h] =~ s/\s//g; }
							if ($array[$h] =~ /\#\# Assembly Statistics \#\#/ ) { 
								$worksheet_stats->write($row_stats, $col_stats, $array[$h], $format_left); $col_stats++;
								$worksheet_stats->write($row_stats, $col_stats, $taxa, $format_left);
							} else {
								$worksheet_stats->write($row_stats, $col_stats, $array[$h], $format_left);
							} $col_stats++;
						} $row_stats++;
					} close (STATS);
	}}}}

	$worksheet_parameters->write($row, $first_col, "ALIGNMENT/MAPPING PHASE:", $format_main_heading);
    $col = $first_col; $row++; $row++;
    $worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
   	$worksheet_parameters->write($row, $col, $date, $format_left); 
   	$col = $first_col; $row++; $row++;
  	$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
   	$worksheet_parameters->write($row, $col, $align_dirname, $format);
   	$col = $first_col; $row++; $row++;
   	$worksheet_parameters->write($row, $col, "DATA:", $format); $col++;
	$worksheet_parameters->write($row, $col, "Option:", $format); 
   	my @option_position = ($row, ($col+1)); $row++;
   	
   	if ($option eq "msa_alignment") {
		my @row_typefile = ($row, $col+1);
   		$worksheet_parameters->write($row, $col, "Type of file:", $format); $row++;
		if ($radseq_like_data) {
   		   	$worksheet_parameters->write($option_position[0], $option_position[1], "RADseq", $format_left); 
			if ($pyRAD_file) {
			   	$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "pyRAD file provided", $format_left); 
			} elsif ($stacks_file) {
			   	$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "STACKS file provided", $format_left); 
		}} elsif ($msa_file) {
		   	$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple Sequence Alignment", $format_left); 
		   	$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "Alignment file provided", $format_left); 
		} elsif ($msa_fasta_folder) {
		   	$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple Sequence Alignment", $format_left); 
		   	$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "Alignment folder provided", $format_left); 
	}} else {
	   	if ($option eq 'user_assembly_contigs') {
		   	$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple references provided", $format_left); 
		} elsif ($option eq 'DOMINO_files') {
		   	$worksheet_parameters->write($option_position[0], $option_position[1], "DOMINO files ", $format_left); 
		} elsif ($option eq 'genome') {
		   	$worksheet_parameters->write($option_position[0], $option_position[1], "Single reference provided", $format_left); 
		} $row++;
		my $col_tag = $option_position[1] - 1;
	   	$worksheet_parameters->write($row, $col_tag, "Files used as reference(s):", $format_bold); $row++;
		for (my $i=0; $i < scalar @contigs_fasta_file_abs_path; $i++) {
		   	$worksheet_parameters->write($row, $col_tag, $contigs_fasta_file_abs_path[$i], $format_left); $row++;
		} $row++; 
	   	$worksheet_parameters->write($row, $col_tag, "Files used for mapping:", $format_bold); $row++;
		for (my $i=0; $i < scalar @clean_fastq_file_abs_path; $i++) {
		   	$worksheet_parameters->write($row, $col_tag, $clean_fastq_file_abs_path[$i], $format_left); $row++;
		} 
	   	$row++; $col = $first_col; 
		$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col++; 
   		my $tmp_col = $col;
		$worksheet_parameters->write($row, $tmp_col, "Read Gap Open penalty ", $format);
    	$worksheet_parameters->write($row, ($tmp_col + 1), $rdgopen, $format_right); $row++;
		$worksheet_parameters->write($row, $tmp_col, "Read Gap Extension penalty ", $format);
    	$worksheet_parameters->write($row, ($tmp_col + 1), $rdgexten, $format_right); $row++;
		$worksheet_parameters->write($row, $tmp_col, "Reference Gap Open penalty ", $format);
    	$worksheet_parameters->write($row, ($tmp_col + 1), $rfgopen, $format_right); $row++;
		$worksheet_parameters->write($row, $tmp_col, "Reference Gap Extension penalty ", $format);
    	$worksheet_parameters->write($row, ($tmp_col + 1), $rfgexten, $format_right); $row++;
		$worksheet_parameters->write($row, $tmp_col, "Mismath penalty ", $format);
    	$worksheet_parameters->write($row, ($tmp_col + 1), $mis_penalty, $format_right); $row++;
   	}

   	$row++; $row++;
	$worksheet_parameters->write($row, $first_col, "MARKER DISCOVERY PHASE:", $format_main_heading);  
    $col = $first_col; 
    $worksheet_parameters->write($row, $col, "STARTED:", $format); $col++;
   	$worksheet_parameters->write($row, $col, $date, $format_left); $row++; $row++;
   	$col = $first_col;
  	$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
   	$worksheet_parameters->write($row, $col, $marker_dirname, $format); $row++; $row++;
   	$col = $first_col;
  	$worksheet_parameters->write($row, $col, "DEVELOPMENT MODULE:", $format); $col++;
   	if ($behaviour eq 'selection') {
	   	$worksheet_parameters->write($row, $col, "Selection", $format); 
	} elsif ($behaviour eq 'discovery') {
   		$worksheet_parameters->write($row, $col, "Discovery", $format);
	} 
	$row++; $row++;	$col = $first_col;
   	$worksheet_parameters->write($row, $col, "PARAMETERS", $format); $col++;
	if ($polymorphism_user) {
		$worksheet_parameters->write($row, $col, "Polymorphic variants were used ", $format); $col++;
	} else {
		$worksheet_parameters->write($row, $col, "Polymorphic variants were not used ", $format); $col++;
	} $row++; $row++;

	unless ($select_markers) {
		$col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, "Parameters that can be modified only via the command-line version", $format_bold); $row++;
		$worksheet_parameters->write($row, $col, "Maximum percentage missing allowed (%) (MPA)", $format); $col++;
		my $x = $missing_allowed * 100; $worksheet_parameters->write($row, $col, $x, $format_right); $row++;
		$col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, "Significance Level Coverage Distribution (SLCD)", $format); $col++;
    	$worksheet_parameters->write($row, $col, $level_significance_coverage_distribution, $format_right); $row++; 
    	
    	$row++; $row++; $col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, "Markers Features",$format_bold); $row++;
		$col = $first_col + 1; 
    	$worksheet_parameters->write($row, $col, "Conserved Length (CL):", $format); $col++;
    	$worksheet_parameters->write($row, $col, $window_size_CONS_range, $format_right); $row++;
		$col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, "Differences in the conserved region (CD):", $format); $col++;
   		$worksheet_parameters->write($row, $col, $window_var_CONS, $format_right); $row++;    
		$col = $first_col + 1; 
    	$worksheet_parameters->write($row, $col, "Variable Length (VL):", $format); $col++;
		$worksheet_parameters->write($row, $col, $window_size_VARS_range, $format_right); $row++;        
		$col = $first_col + 1; 
	}

	if ($variable_divergence) {
		$col = $first_col + 1; 
    	$worksheet_parameters->write($row, $col, "Variable Divergence (%) (VD):", $format); $col++;
		$worksheet_parameters->write($row, $col, $variable_divergence, $format_right); $row++;	
	} else {
		$col = $first_col + 1; 
    	$worksheet_parameters->write($row, $col, "Variable Positions (VP):", $format); $col++;
    	$worksheet_parameters->write($row, $col, "$variable_positions_user_min -- $variable_positions_user_max", $format_right); $row++;    
	}

	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Minimum number of taxa covered (MCT) ", $format); $col++;
	$worksheet_parameters->write($row, $col, $minimum_number_taxa_covered, $format_right); $row++; $row++;
	$col = $first_col + 1;
	$worksheet_parameters->write($row, $col, "Taxa used for marker discovery:", $format_bold); $col++;
	$worksheet_parameters->write($row, $col, "ID", $format_bold); $row++;
	my $counter = 1;
	foreach my $taxa (sort keys %domino_files) {
		if ($domino_files{$taxa}{'taxa'}) {
			$col = $first_col + 1; 
			$worksheet_parameters->write($row, $col, $taxa, $format_left); $col++;	
			$worksheet_parameters->write($row, $col, $counter, $format_right); $counter++;	$row++;
	}} $workbook->close();
	return $excel_woorkbook_name;
}

sub read_phylip_aln {
	
	######################################################################
	##  This function process a MSA into a hash of sequences and it 	##
	## 	calls another function to generate the variation profiles.		##
	##	Jose F. Sanchez 02/2/2016										##
	######################################################################
	
	my $file = $_[0];
	my $region_provided = $_[1];
	
	open (INPUT, $file);
	my ($length, $species, %aln, $counter, $region);
	while (<INPUT>) {
		next if /^#/;
		next if /^\s+$/;
		my $line = $_;
		chomp $line;
		if ($line =~ /\s*(\d+)\s+(\d+)/) {
			$species = $1;
			$length = $2;
			$counter++;
			$region = "msa_".$counter;	
			next;
		}			
		if ($line =~ /(\S+)\s+(.*)/) {			
			$aln{$region}{0}{"name"} = $1;
			my $seq = $2;
			$seq =~ s/\s//g;
			$aln{$region}{0}{"seq"} = $seq;
			my $left_species = $species - 1;
			for (my $i=0; $i < $left_species; $i++) {
				my $next_line= <INPUT>;
				chomp $next_line;
				if ($next_line =~ /(\S+)\s+(.*)/) {
					my $id = $i+1;
					$aln{$region}{$id}{"name"} = $1;
					my $seq = $2;
					$seq =~ s/\s//g;
					$aln{$region}{$id}{"seq"} = $seq;
		}}} else {
			my @array = split("\t", $line);
			my $sequence = $aln{$region}{0}{"seq"};
			$sequence .= $array[0];
			$sequence =~ s/\n//;
			$aln{$region}{0}{"seq"} = $sequence;
			my $left_species = $species - 1;
			for (my $i=0; $i < $left_species; $i++) {
				my $next_line= <INPUT>;
				chomp $next_line;
				my $id = $i+1;
				my @array = split("\t", $next_line);
				my $sequence = $aln{$region}{$id}{"seq"};
				$sequence .= $array[0];
				$sequence =~ s/\n//;
				$aln{$region}{$id}{"seq"} = $sequence;
	}}}
	close(INPUT);
	my @array_taxa;
	foreach my $keys (sort keys %aln) {
		my %alignment;
		if (!$region_provided) { $region_provided = $keys; }
		my $region_Fasta = $msa_dirname."/".$region_provided.".fasta";
		open (OUT_MSA, ">$region_Fasta");
		foreach my $taxa (sort keys %{ $aln{$keys} }) {
			if ($domino_files{$aln{$keys}{$taxa}{"name"}}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
				print OUT_MSA ">".$aln{$keys}{$taxa}{"name"}."\n".$aln{$keys}{$taxa}{"seq"}."\n";
			}
			unless (grep /$aln{$keys}{$taxa}{"name"}/, @array_taxa ) {
				push (@array_taxa, $aln{$keys}{$taxa}{"name"});
		}}
		close(OUT_MSA);	
		$region_provided = "";
	}
	return \@array_taxa;
}

sub retrieve_info {
	
	my $dump_files_ref = $_[0];
	my $hash_Ref = $_[1];
	
	my %hash;
	unless ($hash_Ref eq 1) { %hash = %{$hash_Ref}; }
	&debugger_print("HASH to retrieve"); &debugger_print("Ref", \%hash); &debugger_print("\n");
	&debugger_print("Retrieve info from files");
	my @dump_files = @{ $dump_files_ref };
	for (my $j=0; $j < scalar @dump_files; $j++) {
		if ($dump_files[$j] eq '.' || $dump_files[$j] eq '..' || $dump_files[$j] eq '.DS_Store') { next;}
		unless (-e -r -s $dump_files[$j]) { next; }
		open (DUMP_IN, "$dump_files[$j]");
		while (<DUMP_IN>) {
			my $line = $_; chomp $line;
			my @array = split("\t", $line);
			&debugger_print($line);
			if ($hash_Ref eq 1) { $hash{$array[0]} = $array[1];
			} else { push (@{ $hash{$array[0]}{$array[1]}}, $array[2]); }
		} close (DUMP_IN);
	} 
	&debugger_print("HASH to retrieve"); &debugger_print("Ref", \%hash); &debugger_print("\n");
	return \%hash;
}

sub sliding_window_conserve_variable {

	my $id = $_[0]; my $seq = $_[1]; 
	my @output_file_info;
	my $dna_seq = $$seq; my $seqlen = length($$seq);
	
	## Loop
	for (my $i=0; $i < $seqlen; $i += $SLIDING_user) {
		# AAATATGACTATCGATCGATCATGCTACGATCGATCGATCGTACTACTACGACTGATCGATCGATCGACGACTGAC
		# 		P1		P2|P3							P4|P5						P6
		my $coord_P1 = $i; #print "P1: $coord_P1\n";
		for (my $h=$window_size_CONS_min; $h <= $window_size_CONS_max; $h += $CONS_inc) {
			my $coord_P2 = $i + $h; 			#print "P2: $coord_P2\n";
			my $coord_P3 = $coord_P2 + 1; 		#print "P3: $coord_P3\n";
			if ($coord_P3 > $seqlen) {next;} 

			for (my $j = $window_size_VARS_min; $j <= $window_size_VARS_max; $j += $VAR_inc) {
				my $coord_P4 = $coord_P3 + $j; #print "P4: $coord_P4\n";
				my $coord_P5 = $coord_P4 + 1; #print "P5: $coord_P5\n";
				my $coord_P6 = $coord_P5 + $h; #print "P6: $coord_P6\n";

				if ($coord_P4 > $seqlen){next;} 
				if ($coord_P5 > $seqlen){next;} 
				if ($coord_P6 > $seqlen){next;}				

				#print "\n******\n";
				#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\nPositions:\n";
				#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
				#print "******\n\n";

				## Get substrings
				my $length_string_P1_P2 = $coord_P2 - $coord_P1;
				my $string_P1_P2 = substr($dna_seq, $coord_P1, $length_string_P1_P2);

				my $length_string_P3_P4 = $coord_P4 - $coord_P3;
				my $string_P3_P4 = substr($dna_seq, $coord_P3, $length_string_P3_P4);

				my $length_string_P5_P6 = $coord_P6 - $coord_P5;
				my $string_P5_P6 = substr($dna_seq, $coord_P5, $length_string_P5_P6);
				
				## Count the variable positions in each 
				my $count_string_P1_P2 = $string_P1_P2 =~ tr/1//; ## Conserved 
				my $count_string_P3_P4 = $string_P3_P4 =~ tr/1//; ## Variable
				my $count_string_P5_P6 = $string_P5_P6  =~ tr/1//; ## Conserved
				my $count_N_P1_P2 = $string_P1_P2 =~ tr/N//; ## Conserved 
				my $count_N_P5_P6 = $string_P5_P6  =~ tr/N//; ## Conserved
				
				## More than variations allowed in conserved
				if ($count_string_P1_P2 > $window_var_CONS) { next; } 
				if ($count_string_P5_P6 > $window_var_CONS) { next; }
				if ($count_N_P1_P2 > $window_var_CONS) { next; } 
				if ($count_N_P5_P6 > $window_var_CONS) { next; } 

				#Get marker profile					
				my $total_length = $coord_P6 - $coord_P1;
				my $marker_profile = $string_P1_P2.$string_P3_P4.$string_P5_P6;
				
				### Check variation
				my $expected_var_sites;	my $flag_error=0;
				if ($variable_divergence) {
					# If a minimun divergence, get the expected variable sites for the length
					$expected_var_sites = int($variable_divergence * $total_length);
					unless ($count_string_P3_P4 >= $expected_var_sites) { $flag_error++; }			
				} else {
					# User provides a number of variable positions or a range
					if ($expected_var_sites > $variable_positions_user_min) {
						unless ($expected_var_sites < $variable_positions_user_max) {
							$flag_error++;
						}} else { $flag_error++; }
				}
				
				## Missing % of bases missing
				my $missing_count = $marker_profile =~ tr/N//;
				my $missing_count2 = $marker_profile =~ tr/-//;
				$missing_count += $missing_count2;
				my $missing_count_percent = ($missing_count/$total_length)*100;
				my $percent_total_length = $total_length * $missing_allowed;  ## Default 0.05

				if ($flag_error > 0) {	next; } #if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; }
				if ($missing_count_percent < $percent_total_length) {
					push (@output_file_info, "$$id\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6");
=head DEBUGGER
					#if ($debugger) {
						#print "\n\n***********************************************\n";
						#print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
						#print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\n\nPositions:\n";
						#print "P1:$coord_P1\nP2:$coord_P2\nP3:$coord_P3\nP4:$coord_P5\nP6: $coord_P6\n";
						#print "\nSubsets\n";
						#print "Conserved (P1-P2): ".length($string_P1_P2)."\n"; print "$string_P1_P2\n";
						#print "Variable (P3-P4): ".length($string_P3_P4)."\n";  print $string_P3_P4."\n";
						#print "Conserved (P5-P6): ".length($string_P5_P6)."\n"; print $string_P5_P6."\n";
						#print "\nVariations:\n";
						#print "Count_string_P1_P2: $count_string_P1_P2\n";
						#print "Count_string_P3_P4: $count_string_P3_P4\n";
						#print "Count_string_P5_P6: $count_string_P5_P6\n";
						#print "\nMarker:\n";
						#print "Total Length:$total_length\tVariable Positions:$count_string_P3_P4\tExpected:$expected_var_sites\n";
						#print "Missing allowed: $percent_total_length %\tMissing:$missing_count_percent %\n";
						#print "***********************************************\n";
					#}
=cut

	} else { next;	}}}}
	return \@output_file_info;	
}

sub time_log {
	
	my $current_time = time;
	print DOMINO::time_stamp();
	my $secs = $current_time - $step_time; 
	my $hours = int($secs/3600); 
	$secs %= 3600; 
	my $mins = int($secs/60); 
	$secs %= 60; 
	$step_time = $current_time;
	printf ("\tStep took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
}

sub user_cleanRead_files {
	my $user_cleanRead_files_ref = \@user_cleanRead_files;
	&fastq_files($user_cleanRead_files_ref);
}


__END__

################################################################################################################################################
PILEUP Explanation

	Each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities. 
	At the read base column, 
	+ a dot stands for a match to the reference base on the forward strand, 
	+ a comma for a match on the reverse strand, 
	+ `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. 
	+ A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. 
		The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:
		seq2 156 A 11  .$......+2AG.+2AG.+2AGGG    <975;:<<<<<

		Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. Here is an exmaple of a 4bp deletions from the reference, supported by two reads:
		seq3 200 A 20 ,,,,,..,.-4CACC.-4CACC....,.,,.^~. ==<<<<<<<<<<<::<;2<<
	
	+ A symbol `^' marks the start of a read segment which is a contiguous subsequence on the read separated by `N/S/H' CIGAR operations. 
	+ The ASCII of the character following `^' minus 33 gives the mapping quality. 
	+ A symbol `$' marks the end of a read segment. 
	
	Start and end markers of a read are largely inspired by Phil Green's CALF format. These markers make it possible to reconstruct the read sequences from pileup.
	SAMtools can optionally append mapping qualities to each line of the output. This makes the output much larger, but is necessary when a subset of sites are selected.
	
	Contig_25_MID3  489     C       19      ,.,.,.,...t,.,.,,..     AC9F9F;FIH86I;351II
	Contig_25_MID3  490     A       18      ,.,.,.,...,,.,,,..      AC7F9F:FIH76I851II
	Contig_25_MID3  491     G       18      cCcCcCc.CCccCccc.C      <C;F;F<FIF1:I5<6II
	Contig_25_MID3  492     G       18      ,.,.,.,...,,.,,,..      4E5E7F:FIF.9I3:9II
	Contig_25_MID3  493     A       19      ,.,.,.,...,.,,,..^$.^$. 4B5E7F9FIF9H3:6II55
	Contig_25_MID3  494     A       20      ,.,.,.,...,,.,,,....    4B6E7F9FIF.9H371II55
	Contig_25_MID3  495     C       20      tTtTtTtTGTttTtttTTTT    9B9F;F;FIF19H771II55
	Contig_25_MID3  496     A       20      ,.,.,.,...,,.,,,....    9A;F9B=FIF16F540II88
	
	
#####################################################################################################################


	Contig1	1	A	0	## type 1

	Contig1	2	T	10	..........	## type 2

	Contig1	3	T	10	cccccccccc	## type 3

	Contig1	4	M	10	cccccccccc	## type 4

	Contig1	5	R	10	cccccccccc	## type 5

	Contig1	6	T	10	.........a	## type 6

	Contig1	7	T	10	ccccccccc.	## type 7

	Contig1	8	T	10	....cc..cc	## type 8

	Contig1	9	T	10	ccccccccccca ## type 9

	Contig1	10	T	10	ccccccgggggg ## type 10

	Contig1	11	R	10	AAAAAAGGGG	## type 11

	Contig1	12	R	10	AAAAAAAAAG	## type 12

	Contig1	13	R	10	AAAAAAAAAC	## type 13

	Contig1	14	R	10	AAAAACCCCC	## type 14

	Contig1	15	R	10	TTTTTTTTTTA	## type 15

	Contig1	16	R	10	TTTTTTTTTTC	## type 16

	Contig1	17	R	10	TTTTTTCCCCC	## type 17

	Contig1	18	T	10	ACAAGGGGGG	## type 18

	Contig1	19	T	10	ACCCAAGGGGGG	## type 19

	Contig1	20	R	10	ACAAAAGGGGGG	## type 20

	Contig1	21	R	10	ACCCAAGGGGGG	## type 21

	Contig1	22	R	10	ACCCAAGGGGGGTTTTTT	## type 22

	Contig1	23	R	10	GCCCCCCCCCCCTTTTTT	## type 23
	
	
#################################################
## Merging the sam according to the user input ##
################################################
#print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Merging the SAM files according to user input ", "#"); DOMINO::printHeader("", "#"); print "\n";
#$merge_bam_all_sp = &merge_sam(\@sams_to_parse, \@sam_headers_line);
#print "\n"; &time_log(); print "\n";
#my @tmp = split ("_sorted\.bam", $merge_bam_all_sp);
	
sub merge_sam {
	##########################################################################################
	##	 																					##
	##  This function generates a merge sam of the files user specified						##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################
	
	my $array_sams2parse = $_[0];
	my @sams_to_parse_sub = @$array_sams2parse;
	my $array_samsheaders2parse = $_[1];
	my @sam_headers_line_sub = @$array_samsheaders2parse;

	my @array_sorted_uniq = uniq(@sam_headers_line_sub);
	my $file_header = "header.sam";
	open (HEADER, ">$file_header");
	for (my $i = 0; $i < scalar @array_sorted_uniq; $i++) {
		unless ($array_sorted_uniq[$i] =~ /^\@PG/) {
			unless ($i == $#array_sorted_uniq ) {
				print HEADER $array_sorted_uniq[$i]."\n";
			} else {
				print HEADER $array_sorted_uniq[$i];
	}}}
	close(HEADER);	

	my $sorted_bam; my $name;	

	## Generate sorted bam files
	for (my $j = 0; $j < scalar @sams_to_parse_sub; $j++) {
		my $temp_bam = &generate_bam($sams_to_parse_sub[$j]);
		$sorted_bam .= $temp_bam." ";
	}
	## Generate a name for the merge bam file
	foreach my $taxa (sort keys %domino_files) {
		if ($domino_files{$taxa}{'taxa'}) {
		$name .= $taxa."_";
	}}	
	## Merge the different bam files
	DOMINO::printHeader(" Merging the BAM files ", "%");
	my $system_samtools_merge = $samtools_path." merge -r -@ $num_proc_user -h header.sam ".$name."merged.bam ".$sorted_bam;
	&debugger_print("SAMTOOLS command: ".$system_samtools_merge);
	my $merge_system_command = system ($system_samtools_merge);
	if ($merge_system_command != 0) {
		&printError("Exiting the script. Some error happened when calling SAMtools for merging the different BAM Files...\n"); DOMINO::dieNicely();
	}
	print "Merge...\nDone\n";
	my $file_merged_sorted = &generate_sorted_bam($name."merged.bam");
	return $file_merged_sorted;
}	 
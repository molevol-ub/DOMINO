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
##      perl DM_MarkerScan_1.0.1.pl
##
##    ###########################
##    ### General Information ###
##    ###########################
##      [-h|--help] [-man] [-v|--version] [-MoreInfo]
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
##
##      [-MPA|--missing_perct_allowed float_value]
##      [-MCT|--minimum_number_taxa_covered int_value]
##      [-CD|--conserved_differences int_value] [-PV|--polymorphism]
##      [-VP|--variable_positions range]
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
## TODO: Write a debug message if importing any module fails
## TODO:: use debugger_print for printing Debugging messages

##################################
##	Initializing some variables	##
##################################
my $domino_version = "v1.0.1";
my (
## User options
$folder, $helpAsked, $avoidDelete_tmp_files, $num_proc_user, $window_var_CONS, 
$window_size_CONS_range, $variable_divergence, $bowtie_local, $variable_positions_user_range,
$window_size_VARS_range, $input_type, $manual, $cigar_pct, $rdgopen, $rdgexten, $rfgexten, 
$rfgopen, $MID_taxa_names, $option, $mis_penalty, $msa_fasta_folder, $polymorphism_user,
$level_significance_coverage_distribution, $map_contig_files, $missing_allowed, $keepbam, 
$version, $DOMINO_simulations, $minimum_number_taxa_covered, $avoid_mapping, $further_information,
@user_cleanRead_files, @user_contig_files, $msa_file, $behaviour, $select_markers, $identify_markers,
$debugger,
 
## absolute path
@contigs_fasta_file_abs_path, @clean_fastq_file_abs_path,

## others
%domino_files, $step_time, %discard_contigs, $pyRAD_file, $stacks_file, $radseq_like_data, 
@contigs_fasta_files, %putative_markers, 
%coord_markers, $merge_bam_all_sp, $number_sp, $genome_fasta,
$scripts_path, %msa_all_taxa_files, $dnaSP_flag,

%mapping_contigs
);

my %ambiguity_DNA_codes = (
"R" => ["A","G"], "Y" => ["C","T"], "K" => ["G","T"], "M" => ["A","C"], "S" => ["G","C"], 
"W" => ["A","T"], "B" => ["G","C","T"], "D" => ["A","G","T"], "H" => ["A","C","T"], 
"V" => ["A","C","G"], "N" => ["A","C","G","T"]);

my %ambiguity_DNA_codes_reverse;
foreach my $keys (keys %ambiguity_DNA_codes) {
	my @array = sort @{$ambiguity_DNA_codes{$keys}};
	my $string = join "",@array;
	$ambiguity_DNA_codes_reverse{$string} = $keys;	
}
	
######################
## Get user options	##
######################
GetOptions(
	"h|help" => \$helpAsked,
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

=item B<#################################################################>

=item B<######## DOMINO Marker Discovery/Selection Pipeline ########>

=item B<#################################################################>

=back

=head1 NAME

=over 2

DM_MarkerScan_1.0.0.pl

=back
		
=head1 VERSION

=over 2

DOMINO v1.0.1

=back
	
=head1 SYNOPSIS

=over 2

=item B<>
	
perl DM_MarkerScan_1.0.1.pl

=item B<###########################>
	
=item B<### General Information ###>

=item B<###########################>

[-h|--help] [-man] [-v|--version] [-MoreInfo]

=item B<#########################>

=item B<### Mandatory options ###>

=item B<#########################>

[-option string] [-type_input string] [-o|--outputFolder string] [-taxa_names string] [-VD|--variable_divergence int_value] [-CL|--conserved_length int_value] [-VL|--variable_length range] [-DM|--development_analysis selection/discovery]

=item B<######################>

=item B<###	Optional flags ###>

=item B<######################>

## Input options

[-genome_fasta file] [-user_contig_files file] [-user_cleanRead_files file] [-msa_file file] [-msa_folder directory] [-RADseq_file file]

### Mapping

[-rdgopen int_value] [-rdgexten int_value] [-rfgopen int_value] [-rfgexten int_value] [-mp|--mismatch_penalty int_value] [-bowtie_local] [-map_contig_files] [-max_SoftClipping int_value] [-SLCD|--significance_level_coverage_distribution float_value] [-low_coverage_data]

## Markers

[-MPA|--missing_perct_allowed float_value] [-MCT|--minimum_number_taxa_covered int_value] [-CD|--conserved_differences int_value] [-PV|--polymorphism] [-VP|--variable_positions range]

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

=item B<-PV|--polymorphism>

Use polymorphic variants to estimate nucleotide VD and CD.

=item B<-low_coverage_data>

Use this option if you expect your data to have really low coverage. Sensibility and precision of the markers identified could be decrease.

=item B<-dnaSP>

Use this option along with RADseq file or MSA [file|folder] to report any alignment with a minimun number of variations independently of the taxa given

=item B<-keep_bam_file>

Keep the BAM files generated during the run.

=item B<-No_Profile_Generation|NPG>
	
Use this flag to skip mapping phase and use data from a previous DOMINO run.
	
=item B<-TempFiles>
	
Keep all intermediate files.

=item B<>

=item B<#####################################>

=item B<##### Command Line Examples #########>

=item B<#####################################>

=item B<DOMINO files: single end -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<DOMINO files: single end, No Mapping  -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option DOMINO_files
 -type_input single_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 500 -CD 1 -NPG -MCT 2 -MPA 25 -DM discovery 

=item B<DOMINO files: paired-end -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option DOMINO_files
 -type_input pair_end -o test/ -taxa_names Dmelanogaster,Dsimulans,Dyakuba 
 -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 -DM discovery 

=item B<User provides contigs and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/Dmelanogaster.contigs.fasta -user_contig_files path_to_file2/Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/Dyakuba.contigs.fasta -user_cleanRead_files Dmelanogaster.clean.fastq 
 -user_cleanRead_files Dsimulans.clean.fastq -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<User provides contigs and reads (paired-end) -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option user_assembly_contigs -type_input pair_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -user_cleanRead_files reads_id-Dmelanogaster.clean.R1.fastq -user_cleanRead_files reads_id-Dmelanogaster.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dsimulans.clean.R1.fastq -user_cleanRead_files reads_id-Dsimulans.clean.R2.fastq 
 -user_cleanRead_files reads_id-Dyakuba.clean.R1.fastq -user_cleanRead_files reads_id-Dyakuba.clean.R2.fastq 
 -DM discovery 

=item B<User provides contigs but no reads, and specifies to map contigs vs contigs  -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option user_assembly_contigs -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -user_contig_files path_to_file1/clean_assembly_id-Dmelanogaster.contigs.fasta 
 -user_contig_files path_to_file2/clean_assembly_id-Dsimulans.contigs.fasta 
 -user_contig_files path_to_file3/clean_assembly_id-Dyakuba.contigs.fasta 
 -map_contig_files -DM discovery 

=item B<User provides a reference genome and reads (single end) -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option genome -type_input single_end -o test/ 
 -taxa_names Dmelanogaster,Dsimulans,Dyakuba -VD 0.01 -CL 40 -VL 400 -CD 1 -SLCD 1e-06 -mp 4 
 -genome_fasta path_to_genomes_folder/NCBI_id-Dpseudobscura.fasta
 -user_cleanRead_files Dmelanogaster.clean.fastq -user_cleanRead_files Dsimulans.clean.fastq 
 -user_cleanRead_files Dyakuba.clean.fastq -DM discovery 

=item B<MSA alignment: A single MSA in PHYLIP -- Single File -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/ -msa_file file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single MSA in FASTA -- Single File -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/ -msa_file file.fasta 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Discovery>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -CL 40 -VL 400 -CD 1 -DM discovery 

=item B<MSA alignment: A single PHYLIP file -- Multiple MSA -- Selection>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/ -msa_file multi_msa_file.phy 
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in PHYLIP -- Folder -- Selection>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection 

=item B<MSA alignment: Multiple files in FASTA MSA -- Folder -- Selection>

 perl DM_MarkerScan_v1.0.0.pl -option msa_alignment -o test/  -msa_folder /home/user/MSA
 -taxa_names cow,carp,horse,human -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in pyRAD format -- File -- Selection>

 perl DM_MarkerScan_v1.0.0.pl -option RADseq -o test/  -RADseq_file output.loci
 -taxa_names ind1,ind3,ind5,ind8,ind9 -VD 0.01 -DM selection

=item B<RAD-MSA alignment: single file in STACKS format -- File -- Selection>

 perl DM_MarkerScan_v1.0.0.pl -option RADseq -o test/  -RADseq_file output.fa
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

21 - 09 - 2016

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
my $mapping_parameters_short = $folder_abs_path."/".$datestring."_param.txt";
open (MP_SHORT, ">$mapping_parameters_short");

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
	} 
} elsif ($behaviour eq 'discovery') {
	if ($option eq "RADseq") { &printError("\nThe DOMINO development module DISCOVERY is not suitable for the option $option...\n"); DOMINO::dieNicely(); }
	if (!$window_size_CONS_range || !$window_size_VARS_range) {
		unless ($option eq "msa_alignment") {
			&printError("\nMandatory options are missing...\n"); DOMINO::dieNicely();
	}}
	$identify_markers=1;
} else { &printError("\nPlease choose between selection/discovery...\n"); DOMINO::dieNicely(); }

# Control if missing options
if (!$variable_positions_user_range and !$variable_divergence) {
	&printError("Exiting the script. A range for variable positions or a minimum divergence is missing.\n Use the option -VP|--variable_positions [min::max] or -VD|--variable_divergence [float number]..\n"); DOMINO::dieNicely();
}
unless (!$variable_divergence) { if ($variable_divergence < 0) {$variable_divergence = 0.000000000000000000000000000000001;} ## Set a very small value if -VD 0 
} 
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

##################################################
## Get some info about the files names and tags ##
##################################################
if (!$MID_taxa_names) {
	unless ($option eq "msa_alignment" || $option eq "RADseq") {
		&printError("\nThe option -taxa_names option is missing...\n\nPlease provide it or DOMINO would not continue the process...\n"); DOMINO::dieNicely();
	}
} else {
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
if (!$cigar_pct) { $cigar_pct = 10; }
if (!$window_var_CONS) { $window_var_CONS = 1;}
if (!$level_significance_coverage_distribution) { $level_significance_coverage_distribution = 1e-05; }
if (!$missing_allowed) { $missing_allowed = 0.1;} ## Def. 0.05: When looking for a marker if 5% of the length is missing for any specie, discard 
if ($avoid_mapping) { 
	## If user desires polymorphims, or modifies any other variable affecting the mapping,
	## force to do mapping again to make sure the arrays are correctly designed
	my %variables = (
		'rdgopen' => $rdgopen, 'rdgexten' => $rdgexten, 
		'rfgopen' => $rfgopen, 'rfgexten' => $rfgexten, 
		'mis_penalty' => $mis_penalty, 
		'bowtie_local' => 1, 'poly' => 1,
		'significance_level_coverage_distribution' => $level_significance_coverage_distribution,
	);
	## DOMINO would check for the latest mapping in order to find if parameters are the same.
	my $path_returned = DOMINO::get_earliest("mapping", $folder_abs_path);
	
	&debugger_print("DOMINO::get_earliest subroutine: $path_returned");
	
	if ($path_returned eq 'NO') {
		undef $avoid_mapping;
	} elsif ($path_returned eq 'mapping') {
		undef $avoid_mapping;
		DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
	} else {
		$align_dirname = $path_returned;
		my $time_folder;
		if ($align_dirname =~ /(\d+)\_DM\_mapping$/) { $time_folder = $1; }
		my $parameters_mapping = $align_dirname."/".$time_folder."_param.txt";
		my %tmp_hash;
		my $species_to_check; my $species_to_map;
		
		my $flag_mapping = my $flag_mapping_poly = my $flag_mapping_local = 0; my $taxa = 0;
		if (-e -r -s $parameters_mapping) {
			open(PA_MAP, $parameters_mapping);
			while (<PA_MAP>) {
				my $line = $_; chomp $line;
				&debugger_print("Previous Profile Generation file: $line");				
				my @array_split = split(":", $line);
				my $value = $variables{$array_split[0]};
				if ($array_split[0] eq "poly") { if (!$polymorphism_user) { $flag_mapping++; } $flag_mapping_poly = 1;		
				} elsif ($array_split[0] eq "bowtie_local") { if (!$bowtie_local) { $flag_mapping++; } $flag_mapping_local = 1;
				} elsif ($array_split[0] eq "Mapped") {
					unless ($array_split[1] eq "GenomeID") {
						if ($domino_files{$array_split[1]}{'taxa'}) { $species_to_map++; 
						} else {
							$flag_mapping++; &debugger_print("Taxa: $array_split[1] does not match the current parameters provided...");
							DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
							last;
			}}} elsif ($array_split[0] eq "taxa") {
					if ($domino_files{$array_split[1]}{'taxa'}) { $species_to_check++;
					} else {
						$flag_mapping++;
						&debugger_print("Taxa: $array_split[1] does not match the current parameters provided...");
						DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
						last;
			}} else {
					unless ($value == $array_split[1]) {  $flag_mapping++; }
			}} close(PA_MAP);

			if ($species_to_map ne $number_sp) { $taxa++; }
			if ($species_to_check ne $number_sp) { $taxa++; }
			if ($flag_mapping_poly == 0) { if ($polymorphism_user) {$flag_mapping++;} }
			if ($flag_mapping_local == 0) { if ($bowtie_local) {$flag_mapping++;} }
			if ($flag_mapping > 0) {
				undef $avoid_mapping;
				if ($taxa > 0) {
					DOMINO::printDetails("+ Although option -No_Profile_Generation | -NPG was provided, the profile would be generated again as the different taxa names provided do not much the ones available...\n",$mapping_parameters, $param_Detail_file_markers);
				} else {
					DOMINO::printDetails("+ Although option -No_Profile_Generation was provided, it would be done again as parameters do not much with the available mapping folder...\n",$mapping_parameters, $param_Detail_file_markers);
			}} else {
				DOMINO::printDetails("+ Generation of new profile of variation would be avoided as it has been previously done with the same parameters...OK\n",$mapping_parameters, $param_Detail_file_markers);
				DOMINO::printDetails("+ Alignment Directory: ".$align_dirname." ...OK\n", $mapping_parameters, $param_Detail_file_markers);
		}} else {
			undef $avoid_mapping;
			DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
}}}
&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

## Get ranges
my ($variable_positions_user_min, $variable_positions_user_max);
if ($variable_positions_user_range) {
	if ($variable_positions_user_range =~ m/.*\:\:.*/) {
		($variable_positions_user_min, $variable_positions_user_max) = split("::", $variable_positions_user_range);
	} elsif ($variable_positions_user_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 2::7\n\n"); DOMINO::dieNicely();
	} else { $variable_positions_user_min = $variable_positions_user_max = $variable_positions_user_range; 
}}

# Range Conserved size
my ($window_size_CONS_min, $window_size_CONS_max);
if ($window_size_CONS_range) {
	if ($window_size_CONS_range =~ m/.*\:\:.*/) {
		($window_size_CONS_min, $window_size_CONS_max) = split("::", $window_size_CONS_range);
	} elsif ($window_size_CONS_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
	} else { $window_size_CONS_min = $window_size_CONS_max = $window_size_CONS_range;
}}

## Range Variable size
my ($window_size_VARS_min, $window_size_VARS_max);
if ($window_size_VARS_range) {
	if ($window_size_VARS_range =~ m/.*\:\:.*/) {
		($window_size_VARS_min, $window_size_VARS_max) = split("::", $window_size_VARS_range);
	} elsif ($window_size_VARS_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
	} else { $window_size_VARS_min = $window_size_VARS_max = $window_size_VARS_range;
}}

# MCT
if (!$minimum_number_taxa_covered) { $minimum_number_taxa_covered = $number_sp;  ## Force to be all the taxa
} else {
	if ($minimum_number_taxa_covered > $number_sp) {
		&printError("Minimum number of covered taxa (MCT) is bigger than the number of taxa provided...\n"); DOMINO::dieNicely();
}} 

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
if ($behaviour eq 'selection') {
	DOMINO::printDetails("+ DOMINO development module: Select informative markers ...OK\n", $mapping_parameters, $param_Detail_file_markers);
} elsif ($behaviour eq 'discovery') {
	DOMINO::printDetails("+ DOMINO development module: Discover putative markers ...OK\n", $mapping_parameters, $param_Detail_file_markers);
}

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
	}} else {
		DOMINO::printDetails("+ Multiple sequence alignment has been provided...OK\n", $mapping_parameters, $param_Detail_file_markers); 
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
		if ($assembling_dirname eq 'assembly') {
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
		push (@{ $domino_files{'genome'}{'taxa'}}, $1); 
		&check_file($tmp, $1);
	} else {
		push (@{ $domino_files{'genome'}{'taxa'}}, "1"); 
		&check_file($tmp);
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
	foreach my $keys (keys %domino_files) { 
		if ($domino_files{$keys}{'taxa'}) {
			DOMINO::printDetails("\tName: $keys\n", $mapping_parameters, $param_Detail_file_markers);
			print MP_SHORT "taxa:$keys\n";
}}}
if ($behaviour eq 'selection') {
	DOMINO::printDetails("\n+ Parameters for the selection of molecular markers:\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Variable Length (VL): All the available length would be used for each region\n", $param_Detail_file_markers); 
} elsif ($behaviour eq 'discovery') {
	unless ($option eq "msa_alignment") {
		DOMINO::printDetails("\n\n+ Parameters for the mapping of molecular markers:\n", $mapping_parameters);
		DOMINO::printDetails("\t- Alignment of the reads using: Bowtie2\n", $mapping_parameters);
		if ($bowtie_local) {
			DOMINO::printDetails("\t- Local Bowtie: ON", $mapping_parameters); print MP_SHORT "bowtie_local:1\n";		
		}
		DOMINO::printDetails("\t- Read Gap Open penalty (rdgopen): ".$rdgopen."\n", $mapping_parameters); print MP_SHORT "rdgopen:$rdgopen\n";
		DOMINO::printDetails("\t- Read Gap Extension penalty (rdgexten): ".$rdgexten."\n", $mapping_parameters); print MP_SHORT "rdgexten:$rdgexten\n";
		DOMINO::printDetails("\t- Reference Gap Open penalty (rfgopen): ".$rfgopen."\n", $mapping_parameters); print MP_SHORT "rfgopen:$rfgopen\n";
		DOMINO::printDetails("\t- Reference Gap Open penalty (rfgexten): ".$rfgexten."\n", $mapping_parameters); print MP_SHORT "rfgexten:$rfgexten\n";
		DOMINO::printDetails("\t- Mismath penalty: ".$mis_penalty."\n", $mapping_parameters); print MP_SHORT "mis_penalty:$mis_penalty\n";
		DOMINO::printDetails("\t- Significance Level Coverage Distribution (SLCD): ".$level_significance_coverage_distribution."\n", $mapping_parameters);
		print MP_SHORT "significance_level_coverage_distribution:$level_significance_coverage_distribution\n";		
	}
	DOMINO::printDetails("\n+ Parameters for the development of molecular markers:\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Conserved Length (CL): $window_size_CONS_min -- $window_size_CONS_max (bp)\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Conserved Differences (CD): $window_var_CONS\n", $param_Detail_file_markers);
	DOMINO::printDetails("\t- Variable Length (VL): $window_size_VARS_min -- $window_size_VARS_max (bp)\n", $param_Detail_file_markers);
}

if ($variable_divergence) {
	DOMINO::printDetails("\t- Variable Divergence (VD): $variable_divergence\n", $param_Detail_file_markers);	
} else {
if ($variable_positions_user_min == $variable_positions_user_max) {
	DOMINO::printDetails("\t- Variable Positions (VP): $variable_positions_user_min (bp)\n", $param_Detail_file_markers);			
} elsif ($variable_positions_user_max == 999) {
	$variable_positions_user_max = 99999999;
	DOMINO::printDetails("\t- Variable Positions (VP): > $variable_positions_user_min (bp)\n", $param_Detail_file_markers);			
} else { DOMINO::printDetails("\t- Variable Positions (VP): $variable_positions_user_min -- $variable_positions_user_max (bp)\n", $param_Detail_file_markers);	
}}

## Common markers parameters
unless (!$MID_taxa_names) { DOMINO::printDetails("\t- Minimum number of covered taxa (MCT): ".$minimum_number_taxa_covered."\n", $param_Detail_file_markers); }
if ($polymorphism_user) { 
	DOMINO::printDetails("\t- Polymorphic variants would be detected (PV)...OK\n", $param_Detail_file_markers); print MP_SHORT "poly:1\n";
}
if ($dnaSP_flag) { DOMINO::printDetails("\t- dnaSP option [ON]...OK\n", $param_Detail_file_markers); print MP_SHORT "dnaSP:1\n"; }
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
close(MP_SHORT);
	
################################################################################################
################# 		Mapping/Alignment of the contigs 		################################
################################################################################################
if ($option ne "msa_alignment" and !$avoid_mapping) {

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
		
		## Generate a directory for each one
		my $dir = $align_dirname."/".$reference_identifier; mkdir $dir, 0755; chdir $dir;
		system("ln -s $domino_files{$reference_identifier}{'contigs'}[0] $ref_Fasta");
		push (@{ $domino_files{$reference_identifier}{'dir'} }, $dir);
		print "+ Using as reference: $contigs_fasta\tID: $reference_identifier...OK\n";
		print "+ Generating a new directory $dir....OK\n\n";
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
		
		###############################
		###	 Copy necessary files	### 
		###############################
		print "+ Generating symbolic links for necessary files\n"; system("ln -s $domino_files{$reference_identifier}{'contigs'}[0]");
		
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
		if (scalar @user_cleanRead_files > 0) { ## Already pushed into this array in line 403
			print "+ User clean reads files would be mapped\n";
		} elsif ($map_contig_files) {
			print "+ Contig files would be mapped\n";
		} else { ## Map DOMINO clean reads
			print "+ Clean reads files would be mapped\n";
		}
		print "+ Obtain information of the reference sequences\n";
		my ($reference_hash_fasta_ref, $message) = DOMINO::readFASTA_hashLength($contigs_fasta); ## Obtain reference of a hash
		my $file2dump_seqs = $dir."/contigs_".$reference_identifier."_length.txt";
		push (@{ $domino_files{$reference_identifier}{"hash_reference_file"} }, $file2dump_seqs);
		DOMINO::printDump($reference_hash_fasta_ref,$file2dump_seqs,1);

		foreach my $reads (sort keys %domino_files) {
			chdir $domino_files{$reference_identifier}{'dir'}[0];
			&debugger_print("Change dir to: ".$domino_files{$reference_identifier}{'dir'}[0]);
			unless ($domino_files{$reads}{'reads'}) { next; }
		
			## Mapping Parameters
			my $R_group_id = '--rg-id '.$reads;
			my $R_group_name = ' --rg '.$reads;
			my $threads = ' -p '.$num_proc_user;
			my $mismatches = ' -N 1 --np 0'; ## Do not add penalty if read/ref got an ambiguous base
			my $read_gap_open = ' --rdg '.$rdgopen.','.$rdgexten;
			my $ref_gap_open = ' --rfg '.$rfgopen.','.$rfgexten;
			my $mismatch_penalty = ' --mp '.$mis_penalty;
			my $mapping_file = $domino_files{$reads}{'reads'}[0];
			my $botwie_system = $bowtie_path."bowtie2";
			if ($bowtie_local) { $botwie_system .= " --local"; }
			my $sam_name = $dir."/".$reference_tag."-taxa_".$reads.".sam";
			
			print "+ Aligning reads for $reads against $reference_identifier as reference\n";
			if ($input_type eq 'pair_end') {
				my $second_Read_file = $domino_files{$reads}{'reads'}[1];
				print "+ Reads 1: $mapping_file\n+ Reads 2: $second_Read_file\n";
				$botwie_system .= " -x ".$reference_tag." -q -1 $mapping_file -2 $second_Read_file -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
			} elsif ($input_type eq 'single_end') { ## Illumin single end, 454
				print "+ Reads 1: $mapping_file...\n";
				if ($map_contig_files) { ## Mapping contigs
					$botwie_system .= " -x ".$reference_tag." -f -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
				} else {
					$botwie_system .= " -x ".$reference_tag." -q -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
				}
			}
			print "+ SAM file: $sam_name\n";
			print "+ Mapping now...\n\n";
			&debugger_print("BOWTIE2 command: ".$botwie_system); 
			
			### Map Reads
			print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nMapping Statistics:\n\n"; 
			my $system_bowtie_call = system ($botwie_system);
			if ($system_bowtie_call != 0) {
				&printError("Exiting the script. Some error happened when calling bowtie for mapping the file $mapping_file...\n"); DOMINO::dieNicely();
			} 			
			push (@{$domino_files{$reads}{"SAM::Ref:".$reference_identifier}}, $sam_name);
			print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
			print "\n\n+ Mapping finished for taxa $reads against $reference_identifier\n"; &time_log();		
			
			#################################
			## Get the the reference fasta ##
			#################################
			print "\n+ Checking mapping reads in ".$sam_name."...\n";
			
			## Generate a sam for each contig
			my @temp = split ("\.sam", $sam_name);
			chdir $dir;
			my $dir_tmp = $temp[0]."_SPLIT"; mkdir $dir_tmp, 0755; chdir $dir_tmp;

			#print Dumper $reference_hash_fasta_ref;
			my @number_contigs = sort (keys %$reference_hash_fasta_ref);
			my $scalar = scalar @number_contigs;
			print "+ This SAM file contains $scalar referense sequences...\n";
			system("ln -s $sam_name"); my @temp_name = split ("/", $sam_name);
			my $sorted_bam_file = &generate_bam($temp_name[-1]);
			&generate_index_bam($sorted_bam_file);

			if ($scalar == 1) {
				push (@{ $domino_files{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $temp_name[-1]);
			} else {
				print "+ Splitting SAM file into several parts to speed the computation...\n"; 
				print "+ Using parallel threads ($num_proc_user CPUs)...\n";
				my $parts = int($scalar/$num_proc_user); ## Split SAM into as many CPUs provided
				my @commands; my $iteration = 0; 
				while (1) {
					if (@number_contigs) {
						my @array1 = splice(@number_contigs, 0, $parts + 1);
						my $string = join(" ", @array1);
						my @temp_1 = split ("\.sorted.bam", $sorted_bam_file);
						my @temp_2 = split ("/", $temp_1[0]);
						my $sam_file_part = $dir_tmp."/".$temp_2[-1]."_part-".$iteration.".sam";
						my $command = $samtools_path." view -@ $num_proc_user -Sh -o $sam_file_part $sorted_bam_file $string";
						push (@{ $domino_files{$reads}{"SAM_Parts::Ref:".$reference_identifier} }, $sam_file_part);
						push (@commands, $command); $iteration++; 
				} else { last; }}
				
				my $pm_SAM_split =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
				$pm_SAM_split->run_on_finish( 
					sub { my ($pid, $exit_code, $ident) = @_; 
						print "\t** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n"; 
				} );
				$pm_SAM_split->run_on_start( sub { my ($pid,$ident)=@_; print "\t- SAMTOOLS command for file $ident and PID=$pid started\n"; } );
				for (my $a=0; $a < scalar @commands; $a++) {
					my $pid = $pm_SAM_split->start($a) and next; 
					&debugger_print("SAMTOOLS command: $commands[$a]");	
					my $system_call = system ($commands[$a]);
					$pm_SAM_split->finish($a); # pass an exit code to finish
				}
				$pm_SAM_split->wait_all_children;
				print "\n**********************************************\n";
				print "**** All SAMTOOLS child commands finished ****\n";
				print "**********************************************\n\n";
				&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files); 
	
			}
			
			### Remove multimapping reads	### 
			##  DOMINO checks the SAM files generated and discards bad reads mapping, unmapping reads or multimapping reads. ##
			my @array_files_split = @{ $domino_files{$reads}{"SAM_Parts::Ref:".$reference_identifier}};				
			print "+ Cleaning reads now in parallel threads ($num_proc_user CPUs)...\n";
			## Get files for later dump
			for (my $i=0; $i < scalar @array_files_split; $i++) {
				my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part_".$i.".txt";
				#print "File to dump: ".$file2dump."\n";
				push (@{ $domino_files{$reads}{"DUMP_Parts::Ref:".$reference_identifier} }, $file2dump);
			}		

			my $pm_SAM_parts =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			$pm_SAM_parts->run_on_finish( 
				sub { my ($pid, $exit_code, $ident) = @_; 
					print "\t** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n"; 
			} );
			$pm_SAM_parts->run_on_start( sub { my ($pid,$ident)=@_; print "+ Child process sent for checking file $ident with PID=$pid...\n"; } );
			for (my $i=0; $i < scalar @array_files_split; $i++) {
				my @basename = split("/", $array_files_split[$i]);
				my @name = split(".sam", $basename[-1]);
				my $pid = $pm_SAM_parts->start($name[0]) and next;
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
						#&debugger_print("TOTAL: ".$total."\nPOS: ".$pos."\nNOPOS = XCENT: ".$NOpos."\nXCENTMAX: ".$xcentmax."\n");
						if($xcent <= $xcentmax){
							#&debugger_print("XCENT < XCENTMAX\nGOOD READ!!\n");
							print SAM_OUT $line."\n"; 
							#open (OUT, $domino_files{$reference_identifier}{'mapping_'.$reads}[0]); print OUT $sam[2]."\n"; close (OUT);							
							$good_reads++;
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
				print "\t- Generating coverage statistics for $clean_sorted_bam\n";
			
				## Generate Coverage statistics for the alignment file
				my @tmp_bam_name = split ("\.sorted.bam", $clean_sorted_bam);
				my $coverage_file = $tmp_bam_name[0]."_coverage_stats.txt";
				push (@{ $domino_files_SAM_parts{$reads}{"coverage_Parts::Ref:".$reference_identifier} }, $coverage_file);

				my $coverage_samtools_command = $samtools_path." depth ".$clean_sorted_bam." > ".$coverage_file;
				my $system_coverage_call = system ($coverage_samtools_command); &debugger_print("SAMTOOLS command: $coverage_samtools_command");
				if ($system_coverage_call != 0) {
					&printError("Exiting the script. Some error happened when calling SAMtools for obtaining coverage of file $sorted_bam[$i]...\n"); DOMINO::dieNicely();
				}
				print "\t- Filtering Coverage Stats...\n\n";
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
				foreach my $contig (keys %max_cov) {
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
				DOMINO::printDump(\%domino_files_SAM_parts, $domino_files{$reads}{"DUMP_Parts::Ref:".$reference_identifier}[$i]);
				$pm_SAM_parts->finish($name[0]); # pass an exit code to finish
			}
			$pm_SAM_parts->wait_all_children;
			print "\n";
			print "************************************************\n";
			print "*** All SAM parsing child commands finished ****\n";
			print "************************************************\n\n";
			
			my @dump_files1 = @{ $domino_files{$reads}{"DUMP_Parts::Ref:".$reference_identifier} };
			&retrieve_info(\@dump_files1, \%domino_files);
						
			## Get total reads: discard, good and total
			my $discard_reads = 0; my $good_reads = 0; my $total_reads = 0; my @array;
			push (@array, $domino_files{$reads}{"discard_reads::Ref:".$reference_identifier});
			push (@array, $domino_files{$reads}{"good_reads::Ref:".$reference_identifier});
			push (@array, $domino_files{$reads}{"total_reads::Ref:".$reference_identifier});
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
			my @stats = @{ $domino_files{$reads}{"mean_coverage_Parts::Ref:".$reference_identifier}};
			my ($sum_coverage, $total_positions);
			for (my $i=0; $i < scalar @stats; $i++) {
				my @stats_each = split(":", $stats[$i]);
				$sum_coverage += $stats_each[0];
				$total_positions += $stats_each[1];
			}
			my $mean_coverage = $sum_coverage/$total_positions;
			my $mean = sprintf ("%.3f", $mean_coverage);

			my @max_cov_files = @{ $domino_files{$reads}{"max_coverage_Parts::Ref:".$reference_identifier} };
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
			print "+ Adjusting the SAM/BAM files using parallel threads ($num_proc_user CPUs)\n+ Splitted files would be used...\n";
			my @parts_clean_sam = @{ $domino_files{$reads}{"CLEAN_SAM_Parts::Ref:".$reference_identifier}};

			## Get files for later dump
			for (my $i=0; $i < scalar @parts_clean_sam; $i++) {
				my $file2dump = $dir_tmp."/dump_file_split_".$reads."_Part2_".$i.".txt";
				#print "File to dump: ".$file2dump."\n";
				push (@{ $domino_files{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} }, $file2dump);
			}	
			
			my $pm_SAM_PILEUP =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			$pm_SAM_PILEUP->run_on_finish( 
				sub { my ($pid, $exit_code, $ident) = @_; 
					print "\t** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n"; 
			} );
			$pm_SAM_PILEUP->run_on_start( sub { my ($pid,$ident)=@_; print "+ Child process sent for checking file $ident with PID=$pid...\n"; } );
			for (my $j=0; $j < scalar @parts_clean_sam; $j++) {
				my @basename = split("/", $parts_clean_sam[$j]);
				my @name = split(".sam", $basename[-1]);
				my $pid = $pm_SAM_PILEUP->start($name[0]) and next;
				my %domino_files_SAM_PILEUP;
				my @tmp_sam = split("\.clean.sam", $parts_clean_sam[$j]);
				my $sam_filter = $tmp_sam[0]."_filtered.sam";
				open (SAM_OUT, ">$sam_filter"); open (SAM, "<$parts_clean_sam[$j]");
				while (<SAM>) {
					chomp; my $line = $_;
					if ($line =~ /^@.*/ ) {
						print SAM_OUT $line."\n";	
					} else {
						my @array = split (/\s+/,$line);
						if (!$discard_contigs{$array[2]}) {
							print SAM_OUT $line."\n";	
				}}}
				close(SAM_OUT); close(SAM); undef %discard_contigs;
				print "\t- File checked: Contigs and Reads discarded...\n";
				push (@{ $domino_files_SAM_PILEUP{$reads}{"FILTERED_SAM_Parts::Ref:".$reference_identifier} }, $sam_filter);
				my $bam_filtered_returned = &generate_bam($sam_filter);
				push (@{ $domino_files_SAM_PILEUP{$reads}{"FILTERED_BAM_Parts::Ref:".$reference_identifier} }, $bam_filtered_returned);
				unless ($reads eq $reference_identifier) { ## DO NOT GENERATE FILTER PROFILE FOR REFERENCE
					print "\t- Generate a PILEUP file for $sam_filter...\n";
					my $dir_returned = &generate_filter_PILEUP($bam_filtered_returned, $domino_files{$reference_identifier}{'contigs'}[0], $reference_hash_fasta_ref, $reference_identifier, $reads);
					&debugger_print("Finish PILEUP for $bam_filtered_returned");
					push (@{ $domino_files_SAM_PILEUP{$reads}{"PROFILE::Ref:".$reference_identifier} }, $dir_returned);
				}
				# Dump info into file
				DOMINO::printDump(\%domino_files_SAM_PILEUP, $domino_files{$reads}{"DUMP2_Parts::Ref:".$reference_identifier}[$j]);
				$pm_SAM_PILEUP->finish($name[0]); # pass an exit code to finish
			}
			$pm_SAM_PILEUP->wait_all_children;
			print "\n";
			print "************************************************\n";
			print "*** All SAM parsing child commands finished ****\n";
			print "************************************************\n\n";
			
			my @dump_files = @{ $domino_files{$reads}{"DUMP2_Parts::Ref:".$reference_identifier} };
			&retrieve_info(\@dump_files, \%domino_files);
						
			my $stats_file = $dir."/mapping_statistics_Ref_".$reference_identifier."_Reads_".$reads.".txt";
			push (@{ $domino_files{$reference_identifier}{'stats'} }, $stats_file);
			DOMINO::printHeader(" Filtering Statistics ", "="); 
			open (STATS, ">$stats_file");			
			print STATS "File: $sam_name\n"; 																print "File: $sam_name\n";
			print STATS "Reference: ".$reference_identifier."\n"; 											print "Reference: ".$reference_identifier."\n";
			print STATS "Reads: ".$reads."\n";																print "Reads: ".$reads."\n";
			print STATS "==========================================================\n";						print "==========================================================\n";
			print STATS "Total Contigs: $num_contigs\n";													print "Total Contigs: $num_contigs\n";
			print STATS "Total Reads mapping: $total_reads\n";												print "Total Reads mapping: $total_reads\n";
			print STATS "Coverage Mean: x".$mean."\n";														print "Coverage Mean: x".$mean."\n";
			print STATS "\n******* CONTIGS *******\n";														print "\n******* CONTIGS *******\n";
			print STATS "Contigs discarded: ".$contigs_discarded."\n";										print "Contigs discarded: ".$contigs_discarded."\n";
			print STATS "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n"; 	print "Contigs remaining: ".($num_contigs - $contigs_discarded)."\t( ".$h_cont." %) \n";
			print STATS "\n******* READS *******\n";														print "\n******* READS *******\n";
			print STATS "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";	print "Reads discarded (multimapping, unmapped, low quality reads): $discard_reads\n";
			print STATS "Reads remaining: $good_reads\t( ".$h." %)\n";										print "Reads remaining: $good_reads\t( ".$h." %)\n";
			DOMINO::printHeader("", "="); print "\n";

			#unless ($avoidDelete_tmp_files) { &delete_files_mapping($dir, $reference_identifier); }	
		
		} ## foreach reads
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

		open (MP_SHORT, ">>$mapping_parameters_short");
		if ($genome_fasta) { print MP_SHORT "Mapped:GenomeID:".$reference_identifier."\n";
		} else { print MP_SHORT "Mapped:".$reference_identifier."\n"; }		
		close(MP_SHORT);
		print "\n\n"; DOMINO::printHeader("", "+");DOMINO::printHeader(" Mapping finished for Reference $reference_identifier ", "+");DOMINO::printHeader("", "+"); &time_log(); print "\n";		
		undef $reference_hash_fasta_ref;
	} # foreach reference
	my $file2dump = $align_dirname."/".$datestring."_DUMP.txt";
	DOMINO::printDump(\%domino_files, $file2dump);	

} elsif ($option eq "msa_alignment" && !$avoid_mapping) {

	mkdir $msa_dirname, 0755; chdir $align_dirname; 

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
		} else { $files_ref = DOMINO::file_splitter($file_path, $chars, 'fa'); }		
		push (@{ $domino_files{'RADseq'}{'parts'}}, @$files_ref);		
		&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);

		for (my $i=0; $i < scalar @$files_ref; $i++) {
			my $file2dump = $align_dirname."/dump_file_split_Part_".$i.".txt";
			#print "File to dump: ".$file2dump."\n";
			push (@{ $domino_files{'all'}{'dump_file_split'} }, $file2dump);
		}

		## Implement threads
		print "\n+ For each loci a MSA file would be generated...\n";
		print "+ Parsing splitted files...\n+ Using parallel threads ($num_proc_user CPUs)...\n";
		if ($pyRAD_file) {
	
			&debugger_print("pyRAD file provided...");
			### RADSEQ like data ####
			## Parse pyRAD loci file provided
			my $pm_pyRAD =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			$pm_pyRAD->run_on_finish( 
				sub { my ($pid, $exit_code, $ident) = @_; 
					print "\n** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n\n"; 
			} );
			$pm_pyRAD->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** pyRAD analysis for file $ident and PID=$pid started **\n"; } );
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
						foreach my $keys (keys %hash) { 
							if ($domino_files{$keys}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
								open (OUT, ">>$file"); print OUT ">".$keys."\n".$hash{$keys}."\n"; close(OUT);
							}
						} $counter++; undef %hash; next;
					}
					$line =~ s/\s+/\t/g; 
					$line =~ s/\%/>/g; ## Sometimes there is this symbol or at least in my test set
					my @array = split("\t", $line);
					$array[0] =~ s/\>//;
					$hash{$array[0]} = $array[1];
					if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
						push (@{ $domino_files_pyRAD_split{$array[0]}{'taxa'} }, 1);
					}
				} close (FILE); undef %hash;
				#print Dumper \%domino_files_pyRAD_split;
				DOMINO::printDump(\%domino_files_pyRAD_split, $domino_files{'all'}{'dump_file_split'}[$i]);
				$pm_pyRAD->finish($name[-1]); # pass an exit code to finish
			}
			$pm_pyRAD->wait_all_children; 
			print "***************************************************\n";
			print "**** All pyRAD parsing processes have finished ****\n";
			print "***************************************************\n\n";		
		} elsif ($stacks_file) {

			## Parse STACKS file provided
			&debugger_print("STACKS file provided...");
			my $pm_STACKS =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
			$pm_STACKS->run_on_finish( 
				sub { my ($pid, $exit_code, $ident) = @_; 
					print "\n** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n\n"; 
			} );
			$pm_STACKS->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** STACKS analysis for file $ident and PID=$pid started **\n"; } );
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
						if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') { 
							push (@{ $domino_files_STACKS_split{$sample}{'taxa'} }, 1);
						}
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
			$pm_MSA_folder->run_on_finish( 
				sub { my ($pid, $exit_code, $ident) = @_; 
					print "\n** Child process finished for file $ident; PID=$pid & ExitCode=$exit_code **\n\n"; 
			} );
			$pm_MSA_folder->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** Alignment analysis for file $ident and PID=$pid started **\n"; } );
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
								push (@{ $domino_files_MSA_folder{$titleline}{'taxa'}}, 1);
								push (@{ $domino_files_MSA_folder{'taxa'}{'user_Taxa'} }, $titleline);
							}
							$alignment{$titleline} = $sequence;
						} elsif ($domino_files{$titleline}{'taxa'}) {
							$alignment{$titleline} = $sequence;
					}} close(FILE); $/ = "\n";
					
					if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') { 
						system("ln -s $file_path");
					} else {
						my $out_file = $msa_dirname."/parsed_".$name.".fasta";
						open (OUT, ">$out_file");
						foreach my $seqs (keys %alignment) {
							print OUT ">".$seqs."\n".$alignment{$seqs}."\n";
						} close (OUT);
				}}
				if ($domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
					DOMINO::printDump(\%domino_files_MSA_folder, $domino_files{'all'}{'dump_file_split'}[$i]);
				}
				$pm_MSA_folder->finish($i); # pass an exit code to finish
			}
			$pm_MSA_folder->wait_all_children; 
			print "*********************************************\n";
			print "**** All parsing processes have finished ****\n";
			print "*********************************************\n\n";		
			
			my $files = scalar @array_files;
			print "\n\n+ Parsing of the $files files has been done...\n"; print "\n"; &time_log(); print "\n";
	}}
	
	exit();
	
	my @dump_files = @{ $domino_files{'all'}{'dump_file_split'} };
	&retrieve_info(\@dump_files, \%domino_files);

	unless ($MID_taxa_names) {
		print "\n\n"; DOMINO::printHeader("","#"); print "NOTE:\n\n";
		print "\t+ No taxa names were provided so DOMINO have parsed and checked for the names...\n";	
		print "\t+ DOMINO would use all the taxa available...\n\n";
		foreach my $keys (keys %domino_files) { 
			if ($domino_files{$keys}{'taxa'}) {
				DOMINO::printDetails("\tName: $keys\n", $mapping_parameters, $param_Detail_file_markers);
		}}
		print "\n\n"; DOMINO::printHeader("","#"); print "\n\n";
	}
	&debugger_print("DOMINO Files"); &debugger_print("Ref", \%domino_files);
}

## Move parameters and error file to folder
File::Copy::move($mapping_parameters, $align_dirname);
File::Copy::move($mapping_parameters_short, $align_dirname);
unless (-z $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $align_dirname); }


## TODO: Check previous	run and check if it is the same: MERGE with previous control steps
if ($avoid_mapping) { ## No_Profile_Generation|NPG: just get files
	########################################################################################
	###################### NoMapping: just get files #######################################
	########################################################################################

	print "+ Files would be obtained...\n\n";
	my $path_returned = DOMINO::get_earliest("mapping", $folder_abs_path);
	my $file2dump;
	if ($path_returned =~ /(\d+)\_DM\_mapping$/) { 
		$file2dump = $align_dirname."/".$1."_DUMP.txt";
	}
	&debugger_print($path_returned."\n".$file2dump);
	my @array_files;
	push (@array_files, $file2dump);

	## TO FIX: Check previous run and check if it is the same: MERGE with previous control steps
	my %hash; &retrieve_info(\@array_files, \%domino_files);
}

##########################################################################################
###################### MARKER DEVELOPMENT ################################################
##########################################################################################

### MSA alignment
if ($option eq "msa_alignment") {
=head
	my $array_files_fasta_msa_ref = DOMINO::readDir($msa_fasta_folder);
	@array_files_fasta_msa = @$array_files_fasta_msa_ref;
	my $fileRADmarkers = $align_dirname."/markers_parsed.txt";
	open (RAD, ">$fileRADmarkers");
	print "+ Checking files in folder generated...\n";
	for (my $i=0; $i < scalar @array_files_fasta_msa; $i++) {
		if ($array_files_fasta_msa[$i] eq "." || $array_files_fasta_msa[$i] eq ".." || $array_files_fasta_msa[$i] eq ".DS_Store"  ) {next;}
		my $file_path = $msa_fasta_folder."/".$array_files_fasta_msa[$i];
		unless (-f -e -r -s $file_path) {
			&printError("File $file_path is not readable or empty. Please discarded from the folder...\n"); DOMINO::dieNicely();
		}
		
		## Check marker
		my $hash_ref_msa = DOMINO::readFASTA_hash($file_path);
		my @taxa4marker = keys %$hash_ref_msa;
		if (scalar @taxa4marker < $minimum_number_taxa_covered) { next; }
		my $valueReturned = &check_marker_pairwise($hash_ref_msa);
		if ($valueReturned == 1) {
			my $string2print_ref = &check_marker_ALL($file_path);
			unless ($string2print_ref eq 'NO') {
				my $string2print = join("\t", @$string2print_ref);
				print RAD $file_path."\t".$string2print."\n";
	}}} 
	print "+ Files checked and everything seems OK...\n\n";
	close(RAD);	print "\n"; &time_log(); print "\n";

	##########################################################
	############# MARKER SELECTION ###########################
	##########################################################
	## Print info about where to print info and error
	DOMINO::printHeader("", "+"); DOMINO::printHeader(" Analysis of Molecular Markers started ", "+"); DOMINO::printHeader("", "+"); print "\n";
	DOMINO::printDetails("\n+ Markers Directory: ".$marker_dirname." ...OK\n", $param_Detail_file_markers);
	if (-d $marker_dirname) {  
		File::Copy::move($marker_dirname, $marker_dirname."_old_".$random_number);
		DOMINO::printDetails("+ Changing an old folder named as $marker_dirname to $marker_dirname"."_old_"."$random_number...OK\n", $param_Detail_file_markers);     						
	} 
	mkdir $marker_dirname, 0755; chdir $marker_dirname; &debugger_print("Changing dir to $marker_dirname");
	my $profile_dir = $marker_dirname."/PROFILES";	mkdir $profile_dir, 0755;
	my $msa_dir = $marker_dirname."/MSA_markers"; mkdir $msa_dir, 0755;
	my $tmp_dir;
	if ($select_markers) { $tmp_dir = $marker_dirname."/tmp_concatenate"; mkdir $tmp_dir, 0755; }
	
	#################################################################################
	##	Once the coordinates are found, print different files with the information ##
	#################################################################################	
	### open Output and Error file
	my $output_file = "DM_markers-summary.txt";
	open (OUT,">$output_file")or die "Cannot write the file $output_file";
	print OUT "Region\t\tTaxa_included\tVariable_Positions\tEffective_length\tVariation(%)\n";
	print "+ Printing selected markers in $output_file...\n";

	### USE THE SUBROUTINE print_Excel and control if radseq_like_data
	my @coord_markers_select;
	my %files;

	open (RAD_in, "<$fileRADmarkers");
	## msa_file_path	taxa_included	var_sites	length	profile	effective_length
	##	0					1				2			3		4			5
	while (<RAD_in>) {
		chomp;
		my ($file, $taxa, $var_sites, $length_string, $string_profile, $effective_length) = split("\t", $_);
		my @array_taxa_split = split(",", $taxa);
		my $region_id;
		my @path_file = split("/", $file);
		if ($path_file[-1] =~ /.*\_(loci\_\d+)\.fasta/) { 
			$region_id = $1; 
		} elsif ($path_file[-1] =~ /(.*)\.fasta/) {
			$region_id = $1; 
		}
		$files{$region_id} = $file;

		## Control steps
		unless ($variable_divergence) {
			if ($var_sites > $variable_positions_user_max) {next;}
		}
		if (scalar @array_taxa_split < $minimum_number_taxa_covered) { next; }   		

		## Marker seems ok...
		my %concatenate_tmp_files;
		if ($select_markers) {
			foreach my $keys (keys %MID_species_hash) {
				$concatenate_tmp_files{$keys} = $tmp_dir."/tmp_".$keys."_concatenate.fasta";
			}
			## print results into txt and xls
			my $variation_perc = ($var_sites/$effective_length)*100;
			my $h = sprintf ("%.3f", $variation_perc);
			my $string = $region_id."\t".$taxa."\t".$var_sites."\t".$effective_length."\t".$h;		
			push (@coord_markers_select, $string);
			print OUT $string."\n";
					
			## Print profile
			my $profile_dir_file = $profile_dir."/".$region_id."_profile.txt";
			open (PRF, ">$profile_dir_file");
			print PRF ">".$region_id."\n".$string_profile."\n";
			close (PRF);
	
			## Print sequence
			my $msa_fasta = $msa_dir."/".$region_id.".fasta";
			open (MSA, ">$msa_fasta");
			my $hash_file_Ref = DOMINO::readFASTA_hash($file);
			foreach my $keys (keys %MID_species_hash) {
				if ($$hash_file_Ref{$keys}) {
					my $file_concat = $concatenate_tmp_files{$keys};
					open (FILE, ">>$file_concat");
					print MSA ">".$keys."\n".$$hash_file_Ref{$keys}."\n";
					print FILE $$hash_file_Ref{$keys}." ";
					close (FILE);
				} else {
					my $file_concat = $concatenate_tmp_files{$keys};
					open (FILE, ">>$file_concat");
					my @array = ("-") x $length_string;
					my $string2print = join ("", @array);
					print FILE $string2print." ";
					close (FILE);
			}} close(MSA);
		} elsif ($identify_markers) {
			## Identify markers in MSA alignments
			#print Dumper \@array;			
			my $file = $profile_dir."/".$region_id.".txt";
			open (OUT, ">$file");
			print OUT ">".$region_id."-$taxa\n".$string_profile."\n";
			close(OUT);
	}} close (RAD_in); close (OUT);
	
	if ($identify_markers) {
		## Identify markers in MSA alignments
		print "- Checking each contig using a sliding window approach...\n";
		my $array_files_ref = &sliding_window_conserve_variable($profile_dir);
		my $merge_file = $profile_dir."/merge_coordinates.txt";
		open (MERGE, ">$merge_file");
		for (my $i=0; $i < scalar @$array_files_ref; $i++) {
			my $file = $$array_files_ref[$i];
			open (FILE, $file); while (<FILE>) { print MERGE $_; } close(FILE);
		} close(MERGE);
		my $hash_markers_overlap = &check_overlapping_markers($merge_file);
		#print Dumper \%DOMINO_markers;
		my %sequence_hash;
		foreach my $file (keys %files) {
			open(FILE, $files{$file}) || die "Could not open the $files{$file} ...\n";
			$/ = ">"; ## Telling perl where a new line starts
			while (<FILE>) {		
				next if /^#/ || /^\s*$/;
				chomp;
				my ($titleline, $sequence) = split(/\n/,$_,2);
				next unless ($sequence && $titleline);
				chomp $sequence;
				$sequence =~ s/\s+//g;
				$sequence =~ s/\r//g;
				$titleline =~ s/\r//g;
				if ($MID_species_hash{$titleline}) {
					$sequence_hash{$file}{$titleline} = $sequence;
				}} close(FILE); $/ = "\n";
		}
		my $array_Ref = &check_DOMINO_marker($output_file, \%sequence_hash, $msa_dir, $hash_markers_overlap);
		@coord_markers_select = @$array_Ref;
	} else {
		my $concatenate_markers = $marker_dirname."/concatenate_markers.fasta";
		print "+ Concatenating into file $concatenate_markers\n";
		
		my $array_ref_dir = DOMINO::readDir($tmp_dir);
		my @array_dir = @$array_ref_dir;
		open (OUT, ">$concatenate_markers");
		for (my $i=0; $i < scalar @array_dir; $i++) {
			if ($array_dir[$i] eq "." || $array_dir[$i] eq ".." || $array_dir[$i] eq ".DS_Store") { next; }
			if ($array_dir[$i] =~ /tmp\_(.*)\_concatenate\.fasta/) {
				my $string;
				my $file_in = $tmp_dir."/".$array_dir[$i];
				open (IN, "<$file_in");
				while (<IN>) {
					my $line = $_;
					chomp $line;
					$string .= $line;
				} close(IN); 
				print OUT ">$1\n$string\n";
		}} close (OUT);
	}
	print "+ Done...\n+ Retrieving informative locus has been done...\n+ Generating an Excel file for DOMINO markers identified...\n";
	&print_Excel(\@coord_markers_select);
	
	## Move parameters files
	File::Copy::move($param_Detail_file_markers, $marker_dirname);
	unless (-z $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $marker_dirname); }
	unless ($avoidDelete_tmp_files) { 
		remove_tree($profile_dir);
		if ($select_markers) { remove_tree($tmp_dir); }
	}

	## Finish and exit
	&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; exit(0);
=cut
}

## Other types of data
## Print info about where to print info and error
DOMINO::printHeader("", "+"); DOMINO::printHeader(" Analysis of Molecular Markers started ", "+"); DOMINO::printHeader("", "+"); print "\n";
DOMINO::printDetails("\n+ Markers Directory: ".$marker_dirname." ...OK\n", $param_Detail_file_markers);
if (-d $marker_dirname) {  
	File::Copy::move($marker_dirname, $marker_dirname."_old_".$random_number);
	DOMINO::printDetails("+ Changing an old folder named as $marker_dirname to $marker_dirname"."_old_"."$random_number...OK\n", $param_Detail_file_markers);     						
} 
mkdir $marker_dirname, 0755;

##################################
###	Check taxa user specified  ### 
##################################
my $genome_marker_bool = 0;
foreach my $ref_taxa (sort keys %domino_files) { ## For each taxa specified, obtain putative molecular markers
	unless ($domino_files{$ref_taxa}{'contigs'}) {next; }
	if ($genome_marker_bool == 1) {last;}
	print "\n";
	## Create a dir for each taxa
	DOMINO::printHeader(" Checking taxa files user specified ", "#"); 
	my $marker_dir;
	if ($genome_fasta) {
		print "Checking: \tGenome provided: $domino_files{$ref_taxa}{'contigs'}[0]\n\n";
		$genome_marker_bool = 1;
		$marker_dir = $marker_dirname."/DOMINO_markers_Genome";
	} else {
		print "Checking: \t$ref_taxa\n\n";
		$marker_dir = $marker_dirname."/markers_Ref_".$ref_taxa;
	}	
	mkdir $marker_dir, 0755; chdir $marker_dir; &debugger_print("Changing dir to $marker_dir");
	
	## Initialize some variables
	undef %coord_markers; 
	undef %putative_markers; 
	undef $merge_bam_all_sp;

	#######################################
	###		 Copy necessary files		### 
	#######################################
	print "+ Retrieve necessary files\n";
	if ($option eq "msa_alignment") { 
		push ( @{ $domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/merged.profile_ARRAY.txt");
	} else {
=head
		# DO NOT MERGE right now...
		#################################################
		## Merging the sam according to the user input ##
		#################################################
		#print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Merging the SAM files according to user input ", "#"); DOMINO::printHeader("", "#"); print "\n";
		#$merge_bam_all_sp = &merge_sam(\@sams_to_parse, \@sam_headers_line);
		#print "\n"; &time_log(); print "\n";
		#my @tmp = split ("_sorted\.bam", $merge_bam_all_sp);
=cut			
		my @taxa = sort @{ $domino_files{'taxa'}{'user_Taxa'} };
		my @uniq_sort_taxa = uniq(@taxa);
		my $name = join("_", @uniq_sort_taxa);
		push ( @{ $domino_files{$ref_taxa}{'array_all_taxa'} }, $marker_dir."/$name.profile_ARRAY.txt");
	}

	##################################
	## Get the the reference fasta  ##
	##################################
	print "+ Checking the file specified as reference fasta... $domino_files{$ref_taxa}{'contigs'}[0]...OK\n";
	system("ln -s $domino_files{$ref_taxa}{'contigs'}[0]");
	
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
	my $reference_fasta_file = $domino_files{$ref_taxa}{'contigs'}[0];
	my $size_file = DOMINO::get_size($reference_fasta_file);
	my $parts2split = int($size_file/$num_proc_user);
	my $fasta_files_split = DOMINO::fasta_file_splitter($reference_fasta_file, $parts2split, "fasta", $reference_dir);
		&debugger_print("Total Size: $size_file\nCharacters to split: $parts2split");
		&debugger_print("Ref", $fasta_files_split);		
	
	## IMPLEMENT THREADS!!!!!!!!
	my $pm_SPLIT_FASTA =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	$pm_SPLIT_FASTA->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
		print "\t** Child process finished with PID = $pid **\n"; 
	} );
	$pm_SPLIT_FASTA->run_on_start( sub { my ($pid,$ident)=@_; print "\t- Splitting process with PID = $pid started\n"; } );
	for (my $i=0; $i < scalar @{ $fasta_files_split }; $i++) {
		my $pid = $pm_SPLIT_FASTA->start($i) and next;
		open(FILE, $$fasta_files_split[$i]) || die "Could not open the $$fasta_files_split[$i]...\n";
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence;
			$sequence =~ s/\n//g; 
			$titleline =~ s/\s/\t/g;
			my $file = $reference_dir."/".$titleline.".fasta";
			open(OUT, ">$file");
			print OUT ">".$titleline."\n".uc($sequence)."\n";
			close (OUT);
		} close(FILE); $/ = "\n";
		$pm_SPLIT_FASTA->finish($i); # pass an exit code to finish
	}
	$pm_SPLIT_FASTA->wait_all_children; 
	print "\n\n";
	print "***********************************************\n";
	print "**** All splitting processes have finished ****\n";
	print "***********************************************\n\n";		
	
	##########################################
	## 	Merge PILEUP information arrays     ##
	##########################################
	print "\n"; DOMINO::printHeader(" Fetching information from all the PROFILEs generated ", "#");
	my %pileup_files; 
	my (@clean_filtered_sam_files, @reference_bam_files, @pileup_Arrays);
	foreach my $reads (sort keys %domino_files) {
		unless ($domino_files{$reads}{'taxa'}) {next; }
		if ($domino_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa}) {
			push ( @reference_bam_files, @{ $domino_files{$reads}{"FILTERED_BAM_Parts::Ref:".$ref_taxa } } );
		}		
		if ($domino_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa}) { 
			push ( @clean_filtered_sam_files, @{ $domino_files{$reads}{"FILTERED_SAM_Parts::Ref:".$ref_taxa } } );
		}
		
		next if ($reads eq $ref_taxa);
		my $PROFILE_dir = $domino_files{$reads}{"PROFILE::Ref:".$ref_taxa}[0];
		my $array_files_ref = DOMINO::readDir($PROFILE_dir);
		my @array_files = @$array_files_ref;
		for (my $y = 0; $y < scalar @array_files; $y++) {
			if ($array_files[$y] eq "." || $array_files[$y] eq ".." || $array_files[$y] eq ".DS_Store") {next;}
			my $tmp = "$PROFILE_dir/$array_files[$y]";		
			if ($array_files[$y] =~ /(.*)\_sequence\.fasta/) {
				$pileup_files{$1}{$reads."_FASTA"} = $tmp;
			} elsif ($array_files[$y] =~ /(.*)\_ARRAY\.txt/) {
				$pileup_files{$1}{$reads} = $tmp;
	}}}
	# Debug print	
	&debugger_print("Ref", \@reference_bam_files); &debugger_print("Ref", \@clean_filtered_sam_files); &debugger_print("Ref", \%pileup_files);

	my $PILEUP_merged_folder_abs_path = $marker_dir."/PROFILE_merge_species";
	mkdir $PILEUP_merged_folder_abs_path, 0755; chdir $PILEUP_merged_folder_abs_path; 
	&debugger_print("Changing dir to $PILEUP_merged_folder_abs_path");

	print "+ Checking profiles of variation for each contig and merging information...\n";
	print "+ Using a sliding window approach...\n"; print "+ Using parallel threads ($num_proc_user CPUs)...\n";			

	## Check each markers using threads
	my %pileup_files_threads;
	#my $pm_MARKER_PILEUP =  new Parallel::ForkManager($num_proc_user); ## Number of subprocesses equal to CPUs as CPU/subprocesses = 1;
	foreach my $contigs (sort keys %pileup_files) {
		#my $pid = $pm_MARKER_PILEUP->start($contigs) and next;
		my @pileup_fasta;
		foreach my $files (keys %{ $pileup_files{$contigs} }) {
			if ($files =~ /.*FASTA$/) {next;}
			my $tmp_hash_reference = DOMINO::readFASTA_hash($pileup_files{$contigs}{$files});
			my %tmp_fasta = %{$tmp_hash_reference};
			foreach my $seqs (keys %tmp_fasta) {
				push (@pileup_fasta, $tmp_fasta{$seqs});
		}}
		print ".";
		
		## Merging variable and conserved information into a unique array...
		my $size = $$fasta_seqs{$contigs};
		my $tmp_string;			
		for (my $i = 0; $i < scalar $size; $i++) {
			my (@tmp, $pb);
			for (my $j = 0; $j < scalar @pileup_fasta; $j++){
				$pb = substr($pileup_fasta[$j], $i, 1);
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
			} else {
				$tmp_string .= $flag; 
		}}
		my $array_all_taxa = $domino_files{$ref_taxa}{'array_all_taxa'}[0]; &debugger_print($array_all_taxa);
		open(OUT_PILEUP, ">>$array_all_taxa");
		my $var_sites = $tmp_string =~ tr/1/1/; ## count variable sites
		my $cons_sites = $tmp_string =~ tr/0/0/; ## count conserved sites
		if ($var_sites != 0 && $cons_sites != 0) { 
			#print ">$contigs\n$tmp_string\n"; 
			my $file = $PILEUP_merged_folder_abs_path."/".$contigs."_merged_ARRAY.txt";
			open (FH, ">$file"); print FH ">$contigs\n$tmp_string\n"; close(FH);
			push (@{ $pileup_files_threads{$contigs}{'mergeProfile'} }, $file);
			print OUT_PILEUP ">$contigs\n$tmp_string\n"; 
			my $fileReturned = &sliding_window_conserve_variable(\$contigs, \$tmp_string, $PILEUP_merged_folder_abs_path);
			if ($fileReturned eq 0) {  undef $pileup_files_threads{$contigs}; last;
			} else { push (@{ $pileup_files_threads{$contigs}{'mergeCoord'} }, $$fileReturned); }
		} 
		close (OUT_PILEUP); 
		print "."; 
		
		######################################################################
		## Check the coordinates foreach taxa against the merge statistics  ##
		######################################################################

		foreach my $taxa (sort keys %domino_files) {
			unless ($domino_files{$taxa}{'taxa'}) { next; }
			if ($taxa eq $ref_taxa) {next;}
				#$pileup_files_threads{$contigs}{$taxa}
				#$pileup_files_threads{$contigs}{$taxa."_FASTA"}
			my $string = $window_size_VARS_range;	$string =~ s/\:\:/-/;
			my $string2 = $window_size_CONS_range;  $string2 =~ s/\:\:/-/;
			my ($output_file, $error_file, $file);
			if ($variable_divergence) {
				$file = $PILEUP_merged_folder_abs_path."/".$contigs.".id_".$taxa."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
			} else { 
				$file = $PILEUP_merged_folder_abs_path."/".$contigs.".id_".$taxa."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
			}		
			$output_file = $file.".out"; $error_file = $file.".err";			
			if ($pileup_files{$contigs}{$taxa}) {
				my $pileup_each_taxa = $pileup_files{$contigs}{$taxa};
				&get_coordinates_each_taxa(\$pileup_each_taxa, \$pileup_files_threads{$contigs}{'mergeCoord'}[0], $taxa, \$output_file, \$error_file);
				push (@{ $pileup_files_threads{$contigs}{'eachTaxaCoord'} }, $output_file);
		}}
		if (!$pileup_files_threads{$contigs}{'mergeCoord'}[0]) { next; } ## If it is empty, jump to next...
		
		##########################################
		## Get Coordinates of Molecular Markers ##
		##########################################
		
		my $coord_retrieved_ref = &get_shared_coordinates(\@{$pileup_files_threads{$contigs}{'eachTaxaCoord'} }, \$pileup_files_threads{$contigs}{'mergeCoord'}[0], $ref_taxa);
		## Print in tmp file for sorting and obtaining unique
		chdir $PILEUP_merged_folder_abs_path;
		my $tmp_file = $PILEUP_merged_folder_abs_path."/".$contigs."_markers_shared.txt";
		open(TMP, ">$tmp_file");
		for my $scaffold (sort keys %{ $coord_retrieved_ref } ) {    		
			## Check we find markers for all the taxa desired
			my @string = split(";", $scaffold);
			my @array_taxa_split = @{ $$coord_retrieved_ref{$scaffold}{'taxa'} };
			my $species2find = $minimum_number_taxa_covered;
  			if (scalar @array_taxa_split < $species2find) { next; }   		
    		## Write DOMINO Markers Coordinates in tmp txt file
	   		my $string = join("\t", @string);
	   		my @sort_taxa = sort(@array_taxa_split);
	   		my $string_taxa = join(",", @sort_taxa);
    		print TMP "$string\t$string_taxa\n"; 
		} close(TMP);	

		## Collapse markers
		my $file_markers_collapse = &check_overlapping_markers($tmp_file); ## keep record of taxa: TODO

=head		
		# Retrieve fasta sequences...
		my %fasta_msa;
		foreach my $taxa (sort keys %domino_files) {
			unless ($domino_files{$taxa}{'taxa'}) { next; }
			if ($taxa eq $ref_taxa) {next;}
			$fasta_msa{$contigs}{$taxa} = $pileup_files{$contigs}{$taxa."_FASTA"};
		}
		my $refence_file_contig = $domino_files{$ref_taxa}{'REF_DIR'}[0]."/".$contigs.".fasta";
		if (-e -r -s $refence_file_contig) {
			$fasta_msa{$contigs}{$ref_taxa} = $refence_file_contig
		} else { print "ERROR: Something happenned when writting files to REF_DIR, skipping analysis for $contigs...\n"; next; }
		
		## Debug print Dumper \%fasta_msa;
		
		## Retrieve sequences for each marker:: TODO
		open (MARKER, $tmp_file);
		while (<MARKER>) {
			chomp;
			my $line = $_;
			#my $array_return = &check_DOMINO_marker($output_file, \%fasta_msa, $markers_msa_folder, $hash_markers_ref_collapse);
		}
		close(MARKER);
=cut		
		
		#$pm_MARKER_PILEUP->finish($contigs); # pass an exit code to finish
	} #each marker
	#$pm_MARKER_PILEUP->wait_all_children; 
	print "\n\n";
	print "******************************************************\n";
	print "**** All parallel parsing processes have finished ****\n";
	print "******************************************************\n\n";

	#################################################################
	## Get the information ready for the user to visualize contigs ##
	#################################################################
	print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Getting the information ready to present ", "#"); DOMINO::printHeader("", "#"); 

} #each reference taxa

###########################
## Delete temporary file ##
###########################
print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Deleting temporary files ", "#"); DOMINO::printHeader("", "#");
#&clean_tmp_files_marker_dir($dir); print "\n"; &time_log(); print "\n";
		

=head
## Move parameters files
File::Copy::move($param_Detail_file_markers, $marker_dirname);
unless (-z $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $marker_dirname); }

if ($genome_fasta) {
	## There is no need to clusterize markers
	&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n";
	&print_instructions(); exit();
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
my $all_coordinates_file = $blast_dir."/all_coordinates.fasta";
my $all_contigs_file = $blast_dir."/all_contigs.fasta";
open (ALL_CONTIGS, ">$all_contigs_file"); open (ALL_coordinates, ">$all_coordinates_file");
print "+ Merging different DOMINO markers according to the taxa of reference...\n";
# get sequences
for (my $h = 0; $h < scalar @markers_folders; $h++) {
	if ($markers_folders[$h] eq ".DS_Store" || $markers_folders[$h] eq "." || $markers_folders[$h] eq ".." ) { next; }
	if (-d $markers_folders[$h] || $markers_folders[$h] =~ /markers_Ref_(.*)/) {
		my $coordinate_file = $marker_dirname."/".$markers_folders[$h]."/DM_sequence-markers.fasta";
		my $contig_file = $marker_dirname."/".$markers_folders[$h]."/DM_contigs.fasta";
		open (CONTIG_FILE, $contig_file); while (<CONTIG_FILE>) { print ALL_CONTIGS $_; } close(CONTIG_FILE);
		open (COORD_FILE, $coordinate_file); while (<COORD_FILE>) { print ALL_coordinates $_; } close(COORD_FILE);
}} close(ALL_CONTIGS); close(ALL_coordinates);

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
open (BLAST, $blast_search); while (<BLAST>) {
	my $line = $_;
	chomp $line;
	my @array = split("\t", $line);
	my $query = $array[0]; my $subject = $array[1];
	if ($query eq $subject) { next;}
	if ($array[10] < 1e-20 && $array[3] > $window_size_VARS_min) {    ## how to clusterize...
		if ($first_hit == 0) {
			$first_hit++;
			push(@{$markers2keep{$query}}, $subject); push (@markers_seen, $subject);
		} else {
			my $flag_this = 0;
			foreach my $keys (keys %markers2keep) {
				if (grep /.*$subject.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
				if (grep /.*$query.*/, @{$markers2keep{$keys}}) { $flag_this = 1; last; }
			}
			if ($flag_this == 0) { push(@{$markers2keep{$query}}, $subject); push (@markers_seen, $subject);	}
}}}
close (BLAST);
foreach my $keys (keys %$contig_length_Ref) { unless (grep /$keys/, @markers_seen) { $markers2keep{$keys}++; } }

## Printing definitely Results
my $definitely_results_dirname = $marker_dirname."/DOMINO_markers_Results";
mkdir $definitely_results_dirname, 0755; chdir $definitely_results_dirname;
&debugger_print("Changing dir to $definitely_results_dirname");

my $sam_merged_DOMINO_markers;
unless ($option eq "msa_alignment") {
	my @merge_bam_all_sp = split ("/",$merge_bam_all_sp);
	my $tmp_basename = $merge_bam_all_sp[-1];
	my @tmp_name = split ("\_sorted.bam", $tmp_basename);
	$sam_merged_DOMINO_markers = $tmp_name[0]."_DOMINO_markers.sam";
}

## Check for Contigs
print "+ Clustering markers...\n+ Filtering contigs...\n";
my %hash = %$contig_length_Ref;
open (COORD_CONTIGS, ">DM_sequence-markers.fasta");
foreach my $keys (keys %hash) {
	my $contig_id;
	if ($keys =~ /(.*)_coord_(.*)/) { $contig_id = $1; }
	if ($markers2keep{$keys}) {
		print COORD_CONTIGS ">".$keys."\n".$hash{$keys}."\n";
		$clusterized_contigs_keep{$contig_id}{$keys} = $2;
}}
close(COORD_CONTIGS); undef %hash; undef %markers2keep;

print "+ Getting the information ready for the visualization of results...\n";
my $contig_def_results_sequences = "DM_contigs.fasta";
open (CONTIGS_DEF, ">$contig_def_results_sequences");
my $hash_contigs_ref = DOMINO::readFASTA_hash($all_contigs_file);
my %hash_contigs = %$hash_contigs_ref;
foreach my $keys (keys %hash_contigs) {
	if ($clusterized_contigs_keep{$keys}) {
		print CONTIGS_DEF ">".$keys."\n".$hash_contigs{$keys}."\n";	
}}
close(CONTIGS_DEF);

my %hash4markers2keep;
print "+ Filtering coordinates files...\n\n";
for (my $h = 0; $h < scalar @markers_folders; $h++) {
	if ($markers_folders[$h] eq ".DS_Store" || $markers_folders[$h] eq "." || $markers_folders[$h] eq ".." ) { next; }
	if (-d $markers_folders[$h] || $markers_folders[$h] =~ /markers_Ref_(.*)/) {
		my $id = $1;
		## Check for coordinates		
		my $coor = $marker_dirname."/".$markers_folders[$h]."/DM_markers-coordinates_ref_".$id.".txt";
		open (tmp_COOR, $coor);
		while (<tmp_COOR>) {
			my $line = $_;
			chomp $line;
			next if $line=~ m/^\s*$/o;
			if ($line =~ /.*Divergence/) {next;}
			my @array_split = split("\t",$line);
			my @contig_split = split("\_#",$array_split[0]);
			if ($clusterized_contigs_keep{$contig_split[0]}) {
				my @array_start = split(" - ", $array_split[1]);
				my @array_end = split(" - ", $array_split[3]);
				my $coord_string = $array_start[0]."_".$array_end[1];
				my $coord2check = $contig_split[0]."_coord_".$coord_string;
				if ($clusterized_contigs_keep{$contig_split[0]}{$coord2check}) {
					$hash4markers2keep{$contig_split[0]}{$array_start[0]} = $line;
		}}}
		close(tmp_COOR);		
}}

my $coordinates_def_results = "DM_markers-coordinates.txt";
open (COOR, ">$coordinates_def_results");
print COOR "Contig\t\tConserved_Region\tVariable_Region\tConserved_Region\tMapping_Taxa\tVariable-length\tDivergence\n";
my (@coord_cluster, %rename);
foreach my $contigs (keys %hash4markers2keep) {
	my $counter = 0;
	foreach my $markers (sort {$a <=> $b} keys %{ $hash4markers2keep{$contigs} }) {
		$counter++;
		my $string2change = $hash4markers2keep{$contigs}{$markers};
		my $stringchanged;
		my @array = split("\t",$string2change);
		my $tmp = $contigs."_#".$counter;
		$rename{$array[0]} = $tmp;
		$array[0] = $tmp;
		for (my $i=0; $i < scalar @array; $i++) { $stringchanged .= $array[$i]."\t"; }			
		print COOR $stringchanged."\n"; push (@coord_cluster, $stringchanged);
	} print COOR "\n"; push (@coord_cluster, "undef")
}
close(COOR); &time_log();	print "\n";

## Get MSA markers
my $markers_msa_folder = "MSA_markers"; mkdir $markers_msa_folder, 0755;
foreach my $markers ( keys %rename) {
	if ($msa_all_taxa_files{$markers}) {
		my $hash = DOMINO::readFASTA_hash($msa_all_taxa_files{$markers});
		my @file_name = split("/", $msa_all_taxa_files{$markers});
		my $file;
		if ($file_name[-1] =~ /(.*\_marker_)\d+\.fasta/) {
			$file = $1;			
			if ($rename{$markers} =~ /.*\_\#(\d+)/) {
				$file .= $1.".fasta";			
		}}
		my $file_path = $markers_msa_folder."/".$file;
		open (FILE, ">$file_path");
		foreach my $keys (keys %$hash) {
			print FILE ">".$keys."\n".$$hash{$keys}."\n";
		} close(FILE);
}}

## Print excel for clusterized results
print "+ Generating an Excel file for DOMINO markers coordinates...\n";
&print_Excel(\@coord_cluster);

if ($keepbam) {
	my (@header, @sams2merge);
	print "+ Printing the clusterized contigs alignments into a SAM file...\n";
	for (my $h = 0; $h < scalar @markers_folders; $h++) { 
		if ($markers_folders[$h] eq ".DS_Store" || $markers_folders[$h] eq "." || $markers_folders[$h] eq ".." ) { next; }
		if (-d $markers_folders[$h] || $markers_folders[$h] =~ /markers_Ref.*/) {
			my $files_markers_ref = DOMINO::readDir($marker_dirname."/".$markers_folders[$h]);
			my @markers_files = @$files_markers_ref;
			for (my $i=0; $i < scalar @markers_files; $i++) {
				if ($markers_files[$i] eq ".DS_Store" || $markers_files[$i] eq "." || $markers_files[$i] eq "..") { next; }
				if ($markers_files[$i] =~ /.*Ref.*\.sam/) {
					if (-z $marker_dirname."/".$markers_folders[$h]."/DM_sequence-markers.fasta") {next;}
					print "\t- Filtering $markers_files[$i]...\n";
					my $sam_file = $marker_dirname."/".$markers_folders[$h]."/".$markers_files[$i];
					my $output_sam = $markers_files[$i]."_tmp"; open (SAM_OUT, ">$output_sam"); 
					push (@sams2merge, $output_sam);
					open (SAM, "<$sam_file"); while (<SAM>) {
						chomp; my $line = $_;		
						if ($line =~ /^\@RG/) { push (@header, $line); 
						} elsif ($line =~ /^@.*SN:(\S+).*/ ) { 
							if (defined($clusterized_contigs_keep{$1})) { push (@header, $line); }
						} else {
							my @array = split (/\s+/,$line);
							if (defined($clusterized_contigs_keep{$array[2]})) { print SAM_OUT $line."\n";}
	}} close(SAM);	close (SAM_OUT); }}}}
	close(SAM_OUT);
	my @header_sort = sort @header;
	my @header_uniq = uniq(@header_sort);
	my $output_header = "DM_markers-mapped.sam"; 	open (HEADER, ">$output_header");
	for (my $i=0; $i < scalar @header_uniq; $i++) { print HEADER $header_uniq[$i]."\n"; }
	mkdir "tmp_files", 0755;
	for (my $j=0; $j< scalar @sams2merge; $j++) { 
		open (SAM, "<$sams2merge[$j]"); while (<SAM>) { print HEADER $_; } close(SAM);
		File::Copy::move($sams2merge[$j], "tmp_files")
	}
	close (HEADER); print " + Merge done\n";
	unless ($avoidDelete_tmp_files) { remove_tree("tmp_files"); }
	
	## Generate and Indexed and sorted BAM file
	my $bam_putative_markers = &generate_bam($output_header, 'yes');
	&generate_index_bam($bam_putative_markers);
}
unless ($avoidDelete_tmp_files) {
	print "\n+ Cleaning temporary files...\n";
	remove_tree($blast_dir);
}

## Finish and exit
&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; 
&print_instructions(); exit();

=cut



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

sub check_overlapping_markers {
	## Overlaps and maximizes domino markers obtained
	
	my $file = $_[0];
	my $contig_id;
	my %tmp_hash;
	#if ($debugger) { print "\n\nBefore overlapped...\n";}

	my $marker_counter_tmp = 0;
	open (FILE, $file);
	while (<FILE>) {
		my $line = $_;
		chomp $line;
	
		#print $line."\n";
	
		$line =~ s/ /\t/;
		my @array_lines = split ("\t", $line);		
		$contig_id = $array_lines[0];
		my @a = split(":", $array_lines[1]); ## conserved region1 coord
		my @b = split(":", $array_lines[2]); ## variable region coord
		my @c = split(":", $array_lines[3]); ## conserved region2 coord
		my $taxa = "taxa_".$array_lines[4]; 

		my @coordinates = ($a[0],$a[1],$b[0],$b[1],$c[0],$c[1]);
		my (@array, @array1, @array2, @array3, @array4, @array5);   
		 $array[0][0] = $a[0]; $array1[0][0] = $a[1];			
		$array2[0][0] = $b[0]; $array3[0][0] = $b[1];			
		$array4[0][0] = $c[0];$array5[0][0] = $c[1];
		my @arrayKey = (@array, @array1, @array2, @array3, @array4, @array5, $taxa);   

		if (!$tmp_hash{'marker_'.$marker_counter_tmp}) {
			push ( @{ $tmp_hash{ 'marker_'.$marker_counter_tmp } }, @arrayKey);
		} else {
			if ($tmp_hash{ 'marker_'.$marker_counter_tmp }[-1] =~ /taxa\_(.*)/) {
				unless ($1 eq $array_lines[4]) { 
					$marker_counter_tmp++;
					push ( @{ $tmp_hash{ 'marker_'.$marker_counter_tmp } }, @arrayKey); last;
			}}
			if ($tmp_hash{ 'marker_'.$marker_counter_tmp }[0][0] == $a[0]) {
				for (my $i=1; $i < scalar @coordinates; $i++) {
					push (@{$tmp_hash{ 'marker_'.$marker_counter_tmp }[$i]}, $coordinates[$i]); #keep all coordinates
			}} else {
				$marker_counter_tmp++; 	push ( @{ $tmp_hash{ 'marker_'.$marker_counter_tmp } }, @arrayKey);
	}}} close(FILE);

	my %tmp_coord;
	foreach my $marker (sort keys %tmp_hash ) {		
		my @array_arrays = @{ $tmp_hash{$marker} };
		#print $marker."\t";
		my $string2push_asValue;
		for (my $i = 0; $i < scalar @array_arrays; $i++) {
			if ($array_arrays[$i] =~ /taxa\_(.*)/) { next; }
			my @sub_array = sort {$a <=> $b} @{ $array_arrays[$i] };
			my @sub_array_uniq = uniq(@sub_array);
			$string2push_asValue .= $sub_array_uniq[-1].":";
		}
		my $string2push_asKey = $array_arrays[0][-1];			
		my @taxa_string = split("taxa\_", $tmp_hash{$marker}[-1]);
		$string2push_asValue .= $taxa_string[1];
		$tmp_coord{$string2push_asKey} = $string2push_asValue;
	}
	undef %tmp_hash;
	my %coord_seen; my $marker_counter = 0;

	my @array = split(".txt", $file);
	my $file2return = $array[0]."_overlapped_Markers.txt";
	open (OUT, ">$file2return");
	
	## TODO: control for different taxa in a given marker
	## TODO: control to not exceed min/max CONS or VAR provided 
	
	foreach my $keys_markers (sort keys %tmp_coord) {
		if ($coord_seen{$keys_markers}) { next; }
		$marker_counter++;
		my @array_coordinates;
		push (@array_coordinates, $keys_markers);
		my $bool = 1;
		my ($counter, $bad_counter) = 0;
		while ($bool) {
			$counter++;
			my $new_coord = $keys_markers + $counter;
			if ($tmp_coord{$new_coord}) {
				push (@array_coordinates, $new_coord); 	$coord_seen{$new_coord}++;
			} else {
				$bad_counter++;
				if ($bad_counter == 3) { ## We would consider the same marker if overlapping 3pb
					($bool,$counter,$bad_counter) = 0;
		}}}
		my $coordinate = $tmp_coord{$array_coordinates[-1]};
		my $id = $contig_id."_#".$marker_counter;
		my @array = split (":", $coordinate);
		my @cons1_coord_array = ($array[0], $array[1]);
		my @var_coord_array  = ($array[2], $array[3]);
		my @cons2_coord_array  = ($array[4], $array[5]);
		print OUT $id."\t".$array_coordinates[0].":".$cons1_coord_array[1]."\t".$var_coord_array[0].":".$var_coord_array[1]."\t".$cons2_coord_array[0].":".$cons2_coord_array[1]."\t".$array[-1]."\n";
		
		print $id."\t".$array_coordinates[0].":".$cons1_coord_array[1]."\t".$var_coord_array[0].":".$var_coord_array[1]."\t".$cons2_coord_array[0].":".$cons2_coord_array[1]."\t".$array[-1]."\n";
		
		#if ($debugger) {  print $id."\t".$array_coordinates[0].":".$cons1_coord_array[1]."\t".$var_coord_array[0].":".$var_coord_array[1]."\t".$cons2_coord_array[0].":".$cons2_coord_array[1]."\n"; }
	} close(OUT); undef %tmp_coord;
	#if ($debugger) {print "\n\n";}
	return $file2return;

}

sub check_DOMINO_marker {
	
	my $file = $_[0];
	my $fasta_seq_ref_hash = $_[1];
	my $dir = $_[2];
	my $DOMINO_markers_ref = $_[3];
	
	my %DOMINO_markers = %{$DOMINO_markers_ref};
	my %hash_seqs = %{$fasta_seq_ref_hash};
	my @array2return;	

	open (OUT, ">$file");	
	print OUT "Contig\t\tConserved_Region\tVariable_Region\tConserved_Region\tMapping_Taxa\tVariable-length\tDivergence\n";

	## Check each marker
   	my ($seq_name, $marker_number, $old_seq_name);
	for my $key (sort keys %DOMINO_markers ) {    
    	if ($DOMINO_markers{$key}[0]) {
			if ($key =~ /(.*)\_\#(.*)/) {
				$seq_name = $1;
				$marker_number = $2;
			}
			## Retrieve msa for this marker			
			my %hash;
			foreach my $keys (sort keys %{$hash_seqs{$seq_name}}) {			
				if ($hash_seqs{$seq_name}{$keys} =~ /.*fasta/) {
					open (FILE, $hash_seqs{$seq_name}{$keys});
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
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $DOMINO_markers{$key}[2], $DOMINO_markers{$key}[3], $sequence);
					$hash{$keys} = $seq;
				} else {
					my ($seq_id, $seq) = &fetching_range_seqs($keys, $DOMINO_markers{$key}[2], $DOMINO_markers{$key}[3], $hash_seqs{$seq_name}{$keys});
					$hash{$keys} = $seq;
			}}
			
			## Check pairwise MSA
			my $ref_hash = \%hash;
			my $valueReturned = &check_marker_pairwise($ref_hash);
			if ($valueReturned == 1) { ## if it is variable for each pairwise comparison
				## Print into MSA file
				my $msa_file = $dir."/".$seq_name."_marker_".$marker_number.".fasta";
				open (MSA_OUT, ">$msa_file");
				foreach my $seqs (keys %hash) {
					print MSA_OUT ">".$seqs."\n".$hash{$seqs}."\n";
				}
				close (MSA_OUT);
				my $tmp_abs = abs_path($msa_file);
				$msa_all_taxa_files{$seq_name."_#".$marker_number} = $tmp_abs;
				
				if (!$old_seq_name) {
					$old_seq_name = $seq_name;
				} elsif ($old_seq_name ne $seq_name) {
					print OUT "\n"; $old_seq_name = $seq_name;
					push (@array2return, "undef");
				}
				
				## Get variable positions for the whole marker
				my $array_ref_returned = &check_marker_ALL($msa_file);
				if ($array_ref_returned eq 'NO') { 
					remove_tree($msa_file);				
				} else {
					my @array = @$array_ref_returned; 	
					## 0: taxa join(::) ## 1: variable sites  ## 2: length of the marker ## 3: profile string ## 4: effective length
					# Conserved left
					my $string2write1 = $DOMINO_markers{$key}[0]." - ".$DOMINO_markers{$key}[1];
					# Variable
					my $string2write2 = $DOMINO_markers{$key}[2]." - ".$DOMINO_markers{$key}[3];
					# Conserved Right
					my $string2write3 = $DOMINO_markers{$key}[4]." - ".$DOMINO_markers{$key}[5];
					## Obtain divergence
					my $percentage = ($array[1]/$array[4])*100;
					my $h = sprintf ("%.3f", $percentage);
					my $string2return = "$key\t$string2write1\t$string2write2\t$string2write3\t$array[0]\t$array[4]\t$h";
					#print "Here2!\n".$key."\t".$array[1]."\t".$array[4]."\t".$array[3]."\t".$h."\n";
					print OUT $string2return."\n"; push (@array2return, $string2return);
	
					## Push coordinates into a hash for later analysis
					$coord_markers{$key} = [$DOMINO_markers{$key}[0], $DOMINO_markers{$key}[5], $DOMINO_markers{$key}[2], $DOMINO_markers{$key}[3]];

	}}}}
	close(OUT);	
	#print "Markers:\n"; #print Dumper \%DOMINO_markers;	
	#print "Coordinates:\n"; #print Dumper \%coord_markers;
	return \@array2return;
}

sub check_marker_pairwise {

	## Given a MSA file in a hash, check each taxa pairwise
	my $hash_ref = $_[0];
	my @taxa = keys %$hash_ref;
	my (%seen, %pairwise, %discard);
	for (my $j=0; $j < scalar @taxa; $j++) {
		my $reference = $taxa[$j];
		if (!$domino_files{$reference}{'taxa'}) {next;}
		foreach my $keys (keys %$hash_ref) {
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
			if ($var_sites_sub == 0) {
				## If does not fit the necessary divergence
				$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
				$pairwise{$reference}{$keys} = "NO!";
				next;
			}
			my $percentage_sub = $var_sites_sub/$total_sub;
			#print $string."\t".$var_sites_sub."\t".$con_sites_sub."\t".$percentage_sub."\t".$total_sub."\n";
			if ($variable_positions_user_min) {
				if ($var_sites_sub >= $variable_positions_user_min) { 
					$pairwise{$reference}{$keys} = "YES!";
				} else {
					## If does not fit the necessary divergence
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
					$pairwise{$reference}{$keys} = "NO!";
				}
			} elsif ($variable_divergence) {
				if ($percentage_sub >= $variable_divergence) { 
					$pairwise{$reference}{$keys} = "YES!";
				} else {
					## If does not fit the necessary divergence
					$pairwise{$reference}{$keys} = "NO!";
					$seen{$reference}++; $discard{$keys}++; ## avoid checking if it is not variable
		}}}
		$seen{$reference}++; ## avoid checking again the same pair
	}

	my $flag_fitting = 0;
	foreach my $keys (keys %pairwise) {
		if ($discard{$keys}) {next;}
		foreach my $k (keys %{$pairwise{$keys}}) {
			if ($discard{$k}) {next;}
			if ($pairwise{$keys}{$k} eq 'YES!') {
				$flag_fitting++;
	}}}
	if ($number_sp == 2) {
		if ($flag_fitting == 1) {
			return '1'; #print "YES!\n";
		} else {
			return '0'; #print "NO!\n";
	}} else {
		if ($flag_fitting < $minimum_number_taxa_covered) {
			return '0'; #print "NO!\n";
		} else {
			return '1'; #print "YES!\n"; #print Dumper \%pairwise;
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

	my $file = $_[0];
	my (%hash, $length, @taxa);
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
		my @array = split("", $sequence);
		if (!$domino_files{$titleline}{'taxa'}) {next;}
    	push (@{ $hash{$titleline}}, @array);
    	$length = scalar @array;
    	push (@taxa, $titleline);
	}
	close(FILE); $/ = "\n";
	
	my @profile;
	for (my $i=0; $i < $length; $i++) {
		my $flag_position = 0;
		my @tmp;
		foreach my $seqs (keys %hash) {
			push (@tmp, $hash{$seqs}[$i]);
		}
		my @tmp_uniq = sort @tmp;
		my @tmp_uniq_sort = uniq(@tmp_uniq);
		
		if (scalar @tmp_uniq_sort == 1) {
			if ($tmp_uniq_sort[0] eq 'N') {
				push (@profile, 'N');
			} elsif ($tmp_uniq_sort[0] eq '-') {
				push (@profile, '-');
			} else { push (@profile, '0'); }
		} else {
			## We are assuming the calling has been correctly done and
			## the ambiguity codes are due to polymorphism		
			my $escape_flag = 0;
			my (@tmp, @amb);
			for (my $j=0; $j < scalar @tmp_uniq_sort; $j++) {
				if ($tmp_uniq_sort[$j] eq '-') {
					push (@profile, '-'); ## Gaps would be codify as -
					$escape_flag++;
				} elsif ($tmp_uniq_sort[$j] eq 'N') {
					push (@profile, 'N'); ## Gaps would be codify as -
					$escape_flag++;
				} elsif ($ambiguity_DNA_codes{$tmp_uniq_sort[$j]}) {
					push(@amb, $tmp_uniq_sort[$j]);
				} else {
					push(@tmp, $tmp_uniq_sort[$j]);
			}}
			unless ($escape_flag) {
				if (scalar @amb == 0) { ## No ambiguous code
					push (@profile, '1');			
				} elsif (scalar @amb == 1) { ## 1 amb code
					for (my $i=0; $i < scalar @amb; $i++) {
						my $flag_yes = 0;
						for (my $k = 0; $k < scalar @{ $ambiguity_DNA_codes{$amb[$i]}}; $k++) {
							if (grep /$ambiguity_DNA_codes{$amb[$i]}[$k]/, @tmp) {
								$flag_yes++;
						}}
						if ($flag_yes > 0) {
							if ($polymorphism_user) {
								push (@profile, '1');
							} else { push (@profile, '0'); }
						} else {
							push (@profile, '1');
					}}
				} elsif (scalar @amb > 1) { ## Several
					push (@profile, '1');
	}}}}
	my $string = join ("", @profile);
	my $var_sites = $string =~ tr/1/1/; ## count variable sites
	if ($var_sites == 0) { return 'NO'; }
	
	my $con_sites = $string =~ tr/0/0/; ## count conserved sites
	my $count_length = $con_sites + $var_sites;
	my $missing = $length - $count_length;
	my $missing_allowed_length = $missing_allowed * $length;
	if ($missing > $missing_allowed_length) { return 'NO';}
	my $species = join (",", sort @taxa);
	my @array = ($species, $var_sites, $length, $string, $count_length);
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
	
	my $real_start = $start - 1;
	my $length = $end - $real_start;	
	my $contig_seq_desired = $sequence;
	my $sub_seq = substr ($contig_seq_desired, $real_start, $length);
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
	printf ("Whole process took %.2d hours, %.2d minutes, and %.2d seconds\n", $hours, $mins, $secs); 
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
	} else {
		return "N";
	}		
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
	open (OUT, ">$$output_file") or &printError("Could not open output file $$output_file"); 
	open (ERR, ">$$error_file") or &printError("Could not open error file $$error_file");

	### open coord and check coord of file 1
	open (COORD,"<$$merge_file_coord") or &printError("Could not open merge coordenates $$merge_file_coord");
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
			$expected_var_sites = ($variable_divergence * $total_length)*100;
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
			print ERR "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\n";
			next; 
		} else {
			### Print coordinates if the meet requirements
			if ($missing_count_percent < $percent_total_length) {
				print OUT "$contig_MID\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\n";
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
	
	foreach my $dir (keys %dirs) {
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

sub get_ready_to_view_putative_markers {
	
	my $reference = $_[0];
	my $reference_fasta = $_[1];
	my $folder_pileup_ref = $_[2];
	my $folder_pileup_name = $_[3];
	
	my %coord_contig;
	my $hash_markers_ref_collapse;
	
	my $pileup_Arrays_sub2_ref = $_[4];
	my @pileup_Arrays_sub2 = @$pileup_Arrays_sub2_ref;
	my (%fasta_msa, $id, $seq_name, $old_seq_name, $marker_number);

	my $output_file;	
	if ($option eq "genome") {
		$output_file = "DM_markers-coordinates.txt";		
	} else {
		$output_file = "DM_markers-coordinates_ref_".$reference.".txt";	
	}
	my $markers_msa_folder = "MSA_markers"; mkdir $markers_msa_folder, 0755;
	my $array_return = &check_DOMINO_marker($output_file, \%fasta_msa, $markers_msa_folder, $hash_markers_ref_collapse);
	if ($option eq "genome") { &print_Excel($array_return); }	
	my $output_file_putative_contigs = "DM_contigs.fasta";
	print "- Done...\n- Reading the reference fasta file...\n- Printing reference sequences in $output_file_putative_contigs...\n";
	
	open (OUT, ">$output_file_putative_contigs");
	## Printing coordinates
	my $file_coordinates = "DM_sequence-markers.fasta";
	open (OUT_coord, ">$file_coordinates");
	open(FILE, $reference_fasta) || die "Could not open the $reference_fasta ...\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	chomp $sequence;
    	$sequence =~ s/\n//g; $titleline =~ s/\s/\t/g;
		my @titleline_split = split("\t", $titleline);	
		my $match = $titleline_split[0];
		
		## Get the contig ids containing markers
		my @keys = keys %coord_markers;
		my @grep_keys = grep { /$match.*/ } @keys;
		#print "Match! $match\n";		
		#print scalar @grep_keys; print "\n";
		#print Dumper \@grep_keys;
		#print Dumper \%coord_markers;
		
		if (scalar @grep_keys > 0) {
			print OUT ">".$match."\n".$sequence."\n";
			for (my $h=0; $h < scalar @grep_keys; $h++) { 
				#print $grep_keys[$h]."\n";
				my ($seq_id, $seq, $array) = &fetching_range_seqs($titleline_split[0], $coord_markers{$grep_keys[$h]}[0], $coord_markers{$grep_keys[$h]}[1], $sequence);
				print OUT_coord $seq_id."\n".uc($seq)."\n";
	}}} close(FILE); $/ = "\n"; close(OUT_coord); close(OUT);
	
	if ($keepbam) {
		my @tmp_name = split ("\_sorted.bam", $merge_bam_all_sp);
		my $merge_sam_all_sp = $tmp_name[0].".sam";
		print "- Filtering $merge_sam_all_sp...\n";

		open (SAM, "<$merge_sam_all_sp");
		my $output_sam = "DM_markers-Ref_".$reference.".sam";
		open (SAM_OUT, ">$output_sam"); print "- Printing $output_sam...\n\n";
		while (<SAM>) {
			chomp; my $line = $_;		
			if ($line =~ /^\@RG/) { print SAM_OUT $line."\n"; }
			if ($line =~ /^@.*SN:(.*)\s+/ ) {
				if (defined($putative_markers{$1})) { print SAM_OUT $line."\n";	 }
			} else {
				my @array = split (/\s+/,$line);
				if (defined($putative_markers{$array[2]})) {
					print SAM_OUT $line."\n";	
		}}}
		close(SAM);	close(SAM_OUT);
		## Generate and Indexed and sorted BAM file
		my $bam_putative_markers = &generate_bam($output_sam, 'yes');
		DOMINO::printHeader(" Indexing the BAM file ", "%");
		&generate_index_bam($bam_putative_markers);
	}	
	print "Done...\n\n";	
	undef %fasta_msa;
}

sub get_shared_coordinates {

	##########################################################################################
	##	 																					##
	##  This function obtains the shared coordinates among the different PILEUPs generated,	##
	## 	and gets the information ready to present to the user.								##
	## 		        																		##
	##	Cristina Frias 05/06/2014 	cristinafriaslopez@ub.edu								##
	## 	Jose F. Sanchez 30/09/2014 	jfsanchezherrero@ub.edu									##
	## 		        																		##
	##########################################################################################

	my $file_coord_array_ref = $_[0];
	my $file_coord = $_[1];
	my $reference = $_[2];
	
	my %coord_contig;
	
	###################################
	### Read the merge PILEUP file  ###
	###################################
	open (MERGE_COORD,"<$$file_coord") or die "Cannot open file $$file_coord";
	while(<MERGE_COORD>){
		chomp;
		my $line = $_;
		next if ($line =~ m/^Contig\tCons1_ini:.*/o);
    	next if $line=~ m/^\s*$/o;
    	next if $line=~ /^\#/o;
    	$line =~ s/\s/;/g;
		push (@{ $coord_contig{$line}{'merged'} }, 1); # create empty array for each contig
		push (@{ $coord_contig{$line}{'taxa'} }, $reference); 
	}
	close (MERGE_COORD);	
	if (!%coord_contig) { return; } ## No coordinates were found
	
	###################################
	### Check each taxa coordinates ###
	###################################	
	my @array_file_coordinates = @{ $file_coord_array_ref };
	for (my $i = 0; $i < scalar @array_file_coordinates; $i++) {
		my $name;
		if ($array_file_coordinates[$i] =~ /.*id\_(.*)\-V(D|P).*/) { $name = $1; }
		open (FILE, "<$array_file_coordinates[$i]") || die "Cannot open file $array_file_coordinates[$i]";
		while(<FILE>){
			chomp;
			my $line= $_;
    	    next if ($line =~ m/^Contig\tCons1_ini:.*/o); next if $line=~ m/^\s*$/o; next if $line=~ /^\#/o;
	    	$line =~ s/\s/;/g;
	        if (exists $coord_contig{$line}){ 
       			push (@{ $coord_contig{$line}{'taxa'} }, $name); 
	    }}
		close (FILE);
	}
	#my $file = $$file_coord."_shared.txt";
	#DOMINO::printDump(\%coord_contig, $file);
	return \%coord_contig;
}

sub generate_bam {
	my $sam_file = $_[0];
	my $avoid = $_[1];
	my @temp = split ("\.sam", $sam_file);
	my $name = $temp[0]; my $bam_file = $name.".bam";
	print "\t- Generating a BAM file for $sam_file\n"; 	
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
	print "\t- Sorting the BAM file: $bam_file\n"; 	
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
	print "\t- Generating a SAM file for $bam_file\n"; 	
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
	print "\t- Generating an index bam file for $bam_file\n"; 	
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
						foreach my $keys (keys %polymorphism) {
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
					foreach my $keys (keys %polymorphism) {
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
		foreach my $keys (keys %polymorphism) {
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
				foreach my $keys (keys %polymorphism) {
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
	foreach my $taxa (keys %domino_files) {
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

sub parse_stacks_marker {
	
	my $ref_hash = $_[0];
	my $filename = $_[1];	
	
	my %hash = %$ref_hash;
	my %hash2return;
	foreach my $sample (keys %hash) {
		my @array = @{ $hash{$sample} };
		my $alleles = scalar @array;
		if ($alleles == 1) {
			$hash2return{$sample} = $hash{$sample}[0];
		} elsif ($alleles == 2){
			my @allele1 = split("", $hash{$sample}[0]);
			my @allele2 = split("", $hash{$sample}[1]);
			my @string;
			for (my $i=0; $i < scalar @allele1; $i++) {
				my @tmp;
				push(@tmp,$allele1[$i]);
				push(@tmp,$allele2[$i]);
				my @tmp_uniq = sort @tmp;
				my @tmp_uniq_sort = uniq(@tmp_uniq);
				if (scalar @tmp_uniq_sort == 1) {
					if ($tmp_uniq_sort[0] eq 'N') {
						push (@string, 'N');
					} elsif ($tmp_uniq_sort[0] eq '-') {
						push (@string, '-');
					} else {
						push (@string, $allele1[$i]);
				}} else {
					my $hash;
					$$hash{$allele1[$i]}++;
					$$hash{$allele2[$i]}++;
					my $value = &get_amb_code($hash);
					push (@string, $value);
			}}
			my $tmp = join ("", @string);
			$hash2return{$sample} = $tmp;
		} elsif ($alleles > 2) {
			DOMINO::printError_log("CLocus: ID: $sample contains more than 2 alleles");
	}}
	if (!%hash2return) {return;}
	my $file = $msa_dirname."/".$filename.".fasta";
	open (OUT, ">$file");
	foreach my $keys (keys %hash2return) {
		if ($domino_files{$keys}{'taxa'}) {
			print OUT ">".$keys."\n".$hash2return{$keys}."\n";
		}
	}
	close (OUT);
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

	my $array_markers_ref = $_[0];
	my @array_markers = @$array_markers_ref;
	my $hash_parameters = &get_parameters();
	my $no_parameters;
	if ($hash_parameters == 0) { $no_parameters = 1; }
	
	#################################################################################
	##	Once the coordinates are found, print different files with the information ##
	#################################################################################	
	### open Output and Error file
	my $excel_woorkbook_name = "DM_markers-summary.xls";
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
			$worksheet_markers->write($row_markers, $col_markers, $split[0], $format_left); $col_markers++;	# Contig name
			$worksheet_markers->write($row_markers, $col_markers, $split[1], $format_left); $col_markers++; 	## taxa names
			$worksheet_markers->write($row_markers, $col_markers, $split[2], $format_right); $col_markers++;	## Variable sites
			$worksheet_markers->write($row_markers, $col_markers, $split[3], $format_right); $col_markers++;	## effective length
			$worksheet_markers->write($row_markers, $col_markers, $split[4], $format_right); $row_markers++;	## Variation percentage
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
			if ($array_markers[$i] eq "undef") { $row_markers++; next;}
			my @split = split("\t", $array_markers[$i]);
			$worksheet_markers->write($row_markers, $first_col, $split[0], $format); # Contig name
			my $position = $first_col + 1;
			$worksheet_markers->write($row_markers, $position, $split[1], $format_right);	$position++; ## Conserved left
			$worksheet_markers->write($row_markers, $position, $split[2], $format_right); $position++; ## Variable
			$worksheet_markers->write($row_markers, $position, $split[3], $format_right); $position++; ## Conserved Right
			$worksheet_markers->write($row_markers, $position, $split[4], $format_left); $position++; ## taxa names
			$worksheet_markers->write($row_markers, $position, $split[5], $format_right); $position++; ## variable region length
			$worksheet_markers->write($row_markers, $position, $split[6], $format_right); $row_markers++;	 	## Divergence
	}} 

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
			foreach my $keys (keys %{$$hash_parameters{'clean_data'}{'clean_files'}}) {
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
								$worksheet_stats->write($row_stats, $col_stats, $array[$h], $format_left);
								$col_stats++;
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
}

sub print_instructions {
	
	my $instructions_txt = $marker_dirname."/Instructions.txt";
	open (OUT, ">$instructions_txt");
	
	my $string = "

INSTRUCTION FOR DOMINO FILES AND FOLDERS

Several files and folders has been generated:

+ $param_Detail_file_markers: contains information and parameters regarding the molecular markers developed with the characteristics desired.

+ $align_dirname Folder: contains the alignment/mapping information for each taxa.

+ $marker_dirname Folder and subfolders: contains DOMINO markers detected.

";	
	
	print OUT $string;
	print $string."\n";
	close(OUT);
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
	foreach my $keys (keys %aln) {
		my %alignment;
		if (!$region_provided) { $region_provided = $keys; }
		my $region_Fasta = $msa_dirname."/".$region_provided.".fasta";
		open (OUT_MSA, ">$region_Fasta");
		foreach my $taxa (keys %{ $aln{$keys} }) {
			if ($domino_files{$aln{$keys}{$taxa}{"name"}}{'taxa'} || $domino_files{'taxa'}{'user_Taxa'}[0] eq 'all') {
				print OUT_MSA ">".$aln{$keys}{$taxa}{"name"}."\n".$aln{$keys}{$taxa}{"seq"}."\n";
			}
			unless (grep /$aln{$keys}{$taxa}{"name"}/, @array_taxa ) {
				push (@array_taxa, $aln{$keys}{$taxa}{"name"});
			}
		}
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
	if ($hash_Ref eq 1) { return \%hash; }
}

sub sliding_window_conserve_variable {

	my $id = $_[0]; my $seq = $_[1]; my $dir = $_[2];
	
	#### open files	
	my $output_file_coordinates;
	my $string = $window_size_VARS_range;$string =~ s/\:\:/-/;
	my $string2 = $window_size_CONS_range; $string2 =~ s/\:\:/-/;

	if ($variable_divergence) {
		$output_file_coordinates = $dir."/".$$id."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
	} else {
		$output_file_coordinates = $dir."/".$$id."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
	}		
	open (OUTPUT,">$output_file_coordinates") or &printError("Could not open file $output_file_coordinates for writing results...");		
	my $dna_seq = $$seq; my $seqlen = length($$seq);
	
	## Loop
	for (my $i=0; $i < $seqlen; $i++) {
		# AAATATGACTATCGATCGATCATGCTACGATCGATCGATCGTACTACTACGACTGATCGATCGATCGACGACTGAC
		# 		P1		P2|P3							P4|P5						P6
		my $coord_P1 = $i; #print "P1: $coord_P1\n";
		
		# Set step for CONS range
		my $CONS_inc; ## Set as a variable? TODO
		my $difference_CONS = $window_size_CONS_max - $window_size_CONS_min;
		if ( $difference_CONS >= 20) { $CONS_inc = 2; } else { $CONS_inc = 1; }
		
		# Set step for CONS range
		my $VAR_inc; ## Set as a variable? TODO
		my $difference_VAR = $window_size_VARS_max - $window_size_VARS_min;
		if ( $difference_VAR > 500) {  		$VAR_inc = 10;
		} elsif ($difference_VAR > 300) { 	$VAR_inc = 5;
		} elsif ($difference_VAR > 200) { 	$VAR_inc = 3;
		} elsif ($difference_VAR > 100) { 	$VAR_inc = 2;
		} else { $VAR_inc = 1; }
		
		for (my $h=$window_size_CONS_min; $h <= $window_size_CONS_max; $h += $CONS_inc) {
			my $coord_P2 = $i + $h; #print "P2: $coord_P2\n";
			my $coord_P3 = $coord_P2 + 1; #print "P3: $coord_P3\n";
			if ($coord_P3 > $seqlen) { next; }				

			for (my $j = $window_size_VARS_min; $j <= $window_size_VARS_max; $j += $VAR_inc) {
				my $coord_P4 = $coord_P3 + $j; #print "P4: $coord_P4\n";
				my $coord_P5 = $coord_P4 + 1; #print "P5: $coord_P5\n";
				my $coord_P6 = $coord_P5 + $h; #print "P6: $coord_P6\n";
				if ($coord_P6 > $seqlen) { next; }				

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
					$expected_var_sites = ($variable_divergence * $total_length)*100;
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
					print OUTPUT "$$id\t$coord_P1:$coord_P2\t$coord_P3:$coord_P4\t$coord_P5:$coord_P6\n";
					
					#if ($debugger) {
					#	print "\n\n***********************************************\n";
					#	print "Marker: $coord_P1:$coord_P2 $coord_P3:$coord_P4 $coord_P5:$coord_P6\n";
					#	print "Contig: $$id\nMax: $seqlen\nVL: $j\nCL: $h\n\nPositions:\n";
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
					#	print "Missing allowed: $percent_total_length %\tMissing:$missing_count %\n";
					#	print "***********************************************\n";
					#}
				} else {
					# if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; }
					next;	
	}}}}
	close(OUTPUT);
	if (-e -r -s $output_file_coordinates) { return \$output_file_coordinates;	
	} else { return 0; }	
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

sub sliding_window_conserve_variable {

	my $id = $_[0]; my $seq = $_[1];

	### declare vars
	my ($contig_name, $degenerate_cons, $count, $count_VAR, $count_cons1, $count_cons2, $count_vars, $pos_cons1_ini, $pos_cons1_end, $pos_cons2_ini, $pos_cons2_end, $pos_VAR_ini, $pos_VAR_end);

	#### open files	
	my $output_file_coordinates;
	my $string = $window_size_VARS_range;$string =~ s/\:\:/-/;
	my $string2 = $window_size_CONS_range; $string2 =~ s/\:\:/-/;

	if ($variable_divergence) {
		$output_file_coordinates = $$id."-VD_".$variable_divergence."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
	} else {
		$output_file_coordinates = $$id."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$string2."-CD_".$window_var_CONS."-VL_".$string.".tab";
	}		
	open (OUTPUT,">$output_file_coordinates") or &Error($!);		
	my $dna_seq = $$seq; my $seqlen = length($$seq);
	
	### Loop through the sequence
	my $max = $seqlen - $window_size_CONS_max - $window_size_VARS_max; #	for (my $j = 0; $j < $seqlen - $window_size_CONS + 1; $j++){
	if ($max < 0) {next;}
	for (my $j = 0; $j < $max; $j++){			
		## Set a threshold for $window_size_VARS -VL 200::400				
		for (my $h=$window_size_VARS_min; $h <= $window_size_VARS_max; $h++) {
					
			my $pos_var_ini = $j + $window_size_CONS_min;
			my $pos_var_end = $j + $window_size_CONS_min + $h - 1;
			my $pos_cons2 = $pos_var_end + 1;			
						
			my $cons1_string = substr($dna_seq, $j, $window_size_CONS_min);
			my $VAR_string = substr($dna_seq, $pos_var_ini, $h);
			my $cons2_string = substr($dna_seq, $pos_cons2, $window_size_CONS_min);
			my $count_cons1 = $cons1_string  =~ tr/0//; ## Count the variable positions in each 
			my $count_vars =  $VAR_string 	 =~ tr/1//; ## Count the variable positions in each 
			my $count_cons2 = $cons2_string  =~ tr/0//; ## Count the variable positions in each 
						
			# Conserved Region1
			$pos_cons1_ini = $j + 1; 
			$pos_cons1_end = $j + $window_size_CONS_min;
			
			# Variable Region1
			$pos_VAR_ini = $pos_cons1_end + 1; 
			$pos_VAR_end = $pos_VAR_ini + $h - 1;
			
			# Conserved Region2
			$pos_cons2_ini = $pos_VAR_end + 1;  
			$pos_cons2_end = $pos_cons2_ini + $window_size_CONS_min - 1;
						
			if (($pos_cons2_ini > $seqlen) && ($pos_cons2_end > $seqlen)) { next; }					
			$degenerate_cons = $window_size_CONS_min - $window_var_CONS;
			
			#Get marker profile					
			my $total = $cons1_string.$VAR_string.$cons2_string;
			my $total_length = length($total);
	
			### print coord if ...
			## When parsing the results we are already getting those regions with at least a minimum variation user provided
			my $expected_var_sites;	my $flag_error=0;
			if ($variable_divergence) {
				$expected_var_sites = $variable_divergence * $total_length;
				unless ($count_vars > $expected_var_sites) { $flag_error++; }			
			} else {
				if ($count_vars > $variable_positions_user_min) {
					unless ($count_vars < $variable_positions_user_max) {
						$flag_error++;
					}} else { $flag_error++; }
			}
			my $missing_count = $total =~ tr/N//;
			my $missing_count2 = $total =~ tr/-//;
			$missing_count += $missing_count2;
			my $percent_total_length = $total_length * $missing_allowed; ## Def. 0.05: If 0.05% of the length is missing for any specie, discard
			
			if ($debugger) {
			#	print "Contig: $$id\nMax: $max\nWindow size variable: $h\n\nStart positions\nP0:$j\nP1:$pos_var_ini\nP2:$pos_var_end\nP3:$pos_cons2\n\nSubsets\n";
			#	print "Conserved (P1-P0): ".length($cons1_string)."\n".$cons1_string."\nVariable (P2-P1): ".length($VAR_string)."\n".$VAR_string."\n";
			#	print "Conserved (P3-P2): ".length($cons2_string)."\n".$cons2_string."\n";
			#	print "\nMarker:\nTotal Length:$total_length\tVariable Position:$count_vars\tExpected:$expected_var_sites\n";
			#	print "Missing allowed:$percent_total_length\tMissing:$missing_count\n";
			#	print "Conserved Region 1 CONSERVED:$count_cons1\tConserved Region 2 CONSERVED:$count_cons2\nMarker to Print:\n";
			#	print "$$id\t$pos_cons1_ini:$pos_cons1_end\t$pos_VAR_ini:$pos_VAR_end\t$pos_cons2_ini:$pos_cons2_end\n";
			}
		 	#if ($flag_error > 0) {if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; } next;}
			if ($missing_count < $percent_total_length) {
				if (($degenerate_cons <= $count_cons1 ) && ($degenerate_cons <= $count_cons2)) {
					print OUTPUT "$$id\t$pos_cons1_ini:$pos_cons1_end\t$pos_VAR_ini:$pos_VAR_end\t$pos_cons2_ini:$pos_cons2_end\n";
					#if ($debugger) { print "\n\nOK!\n\n######################################\n\n"; }
				} else {
					#if ($debugger) { print "\n\nERROR!\n\n######################################\n\n"; }
					##print STDERR "$contig\t$pos_cons1_ini:$pos_cons1_end\t$pos_VAR_ini:$pos_VAR_end\t$pos_cons2_ini:$pos_cons2_end\n";
					next;
				}
			}
		}
	}
	
	close(OUTPUT);	
	my $file = &check_overlapping_markers($output_file_coordinates);
	return $file;
	#if ($identify_markers) { return \@output_files; }
}

sub sliding_window_conserve_variable_old {
	
	##########################################################################################
	##	 																					##
	##  This function checks for conserved regions flanking variable regions in the 		##
	##	variation profiles (0/1)															##
	## 		        																		##
	##	Cristina Frias 05/06/2014 															##
	##	Jose F. Sanchez 02/2/2016															##
	## 		        																		##
	##########################################################################################
	
	my $pileup_folder = $_[0];
	my $array_ref_file = DOMINO::readDir($pileup_folder);
	my @array_files = @$array_ref_file;
	
	my @output_files;
	
	### declare vars
	my ($contig_name, $degenerate_cons, $count, $count_VAR, $string, $count_cons1, $count_cons2, $count_vars, $pos_cons1_ini, $pos_cons1_end, $pos_cons2_ini, $pos_cons2_end, $pos_VAR_ini, $pos_VAR_end);

	#### open files	
	for (my $u = 0; $u < scalar @array_files; $u++) {
		if ($array_files[$u] eq '.' || $array_files[$u] eq '..' || $array_files[$u] eq '.DS_Store') {next;}
		my @array_split = split("\.txt",$array_files[$u]);
		my $string = $window_size_VARS_range;
		$string =~ s/\:\:/-/;
		my $output_file_coordinates = $pileup_folder."/";
		if ($variable_divergence) {
			$output_file_coordinates .= $array_split[0]."-VD_".$variable_divergence."-CL_".$window_size_CONS_max."-CD_".$window_var_CONS."-VL_".$string.".tab";
		} else {
			$output_file_coordinates .= $array_split[0]."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$window_size_CONS_max."-CD_".$window_var_CONS."-VL_".$string.".tab";
		}
		
		if ($identify_markers) { push(@output_files, $output_file_coordinates); }
		
		open (OUTPUT,">$output_file_coordinates") or &Error($!);
		print OUTPUT "Contig\tCons1_ini:Cons1_end\tVAR_region_ini:VAR_region_end\tCons2_ini:Cons2_end\n";
		my $ref_hash = DOMINO::readFASTA_hash($pileup_folder."/".$array_files[$u]);
		foreach my $contig (keys %{$ref_hash}) {
			my $dna_seq = $$ref_hash{$contig};
			my $seqlen = length($dna_seq);
	
			### Loop through the sequence
			##$window_size_VARS_min, $window_size_VARS_max
			my $max = $seqlen - $window_size_CONS_max - $window_size_VARS_max; #	for (my $j = 0; $j < $seqlen - $window_size_CONS_max + 1; $j++){
			if ($max < 0) {next;}
			
			for (my $j = 0; $j < $max; $j++){			
				## Set a threshold for $window_size_VARS -VL 200::400				
				for (my $h=$window_size_VARS_min; $h <= $window_size_VARS_max; $h++) {
					
					my $pos_var_ini = $j + $window_size_CONS_max;
					my $pos_var_end = $j + $window_size_CONS_max + $h - 1;
					my $pos_cons2 = $pos_var_end + 1;			#$j + $window_size_CONS_max + $window_size_VARS;
					
					my $cons1_string = substr($dna_seq, $j, $window_size_CONS_max);
					my $VAR_string = substr($dna_seq, $pos_var_ini, $h);
					my $cons2_string = substr($dna_seq, $pos_cons2, $window_size_CONS_max);
					
					my $count_cons1 = $cons1_string  =~ tr/0//; ## Count the variable positions in each 
					my $count_vars =  $VAR_string 	 =~ tr/1//; ## Count the variable positions in each 
					my $count_cons2 = $cons2_string  =~ tr/0//; ## Count the variable positions in each 
					
					if ($count_vars == 0) {next;} ## It is not necessary to keep on calculating anything as there is no variation at all
					
					# Conserved Region1
					$pos_cons1_ini = $j + 1; 
					$pos_cons1_end = $j + $window_size_CONS_max;
		
					# Variable Region1
					$pos_VAR_ini = $pos_cons1_end + 1; 
					$pos_VAR_end = $pos_VAR_ini + $h - 1;
		
					# Conserved Region2
					$pos_cons2_ini = $pos_VAR_end + 1;  
					$pos_cons2_end = $pos_cons2_ini + $window_size_CONS_max - 1;
					
					if (($pos_cons2_ini > $seqlen) && ($pos_cons2_end > $seqlen)) { next; }					
					$degenerate_cons = $window_size_CONS_max - $window_var_CONS;
		
					### print coord if ...
					## When parsing the results we are already getting those regions with at least a minimum variation user provided
					## 			
					#if (($degenerate_cons <= $count_cons1 ) && ($degenerate_cons <= $count_cons2) && ($count_vars >= $expected_var)) {
					
					my $total = $cons1_string.$VAR_string.$cons2_string;
					my $missing_count = $total =~ tr/N//;
					my $missing_count2 = $total =~ tr/-//;
					$missing_count += $missing_count2;
					my $total_length = length($total);
					my $percent_total_length = $total_length * $missing_allowed; ## Def. 0.05: If 0.05% of the length is missing for any specie, discard 
					if ($missing_count < $percent_total_length) {
						if (($degenerate_cons <= $count_cons1 ) && ($degenerate_cons <= $count_cons2)) {
							print OUTPUT "$contig\t$pos_cons1_ini:$pos_cons1_end\t$pos_VAR_ini:$pos_VAR_end\t$pos_cons2_ini:$pos_cons2_end\n";
						} else {
							##print STDERR "$contig\t$pos_cons1_ini:$pos_cons1_end\t$pos_VAR_ini:$pos_VAR_end\t$pos_cons2_ini:$pos_cons2_end\n";
							next;
						}
					}
				}
			}
		}; # close for each contig inside the hash
		close(OUTPUT);
	}
	if ($identify_markers) { return \@output_files; }
}

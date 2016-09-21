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
my $domino_version = "v1.0.0";
my (
## User options
$folder, $helpAsked, $avoidDelete_tmp_files, $num_proc_user, $window_var_CONS, 
$window_size_CONS, $variable_divergence, $bowtie_local, $variable_positions_user_range,
$window_size_VARS_range, $input_type, $manual, $cigar_pct, $rdgopen, $rdgexten, $rfgexten, 
$rfgopen, $MID_taxa_names, $option, $mis_penalty, $msa_fasta_folder, $polymorphism_user,
$level_significance_coverage_distribution, $map_contig_files, $missing_allowed, $keepbam, 
$version, $DOMINO_simulations, $minimum_number_taxa_covered, $avoid_mapping, $further_information,
@user_cleanRead_files, @user_contig_files, $msa_file, $behaviour, $select_markers, $identify_markers,
$debugger,
 
## absolute path
$msa_folder_abs_path, $msa_file_abs_path, @contigs_fasta_file_abs_path, @clean_fastq_file_abs_path,

## others
$step_time, %discard_contigs, %MID_species_hash, $pyRAD_file, $stacks_file, $radseq_like_data, 
@contigs_fasta_files, @array_files_fasta_msa, %max_cov, %coord_contig, %putative_markers, 
%coord_markers, $merge_bam_all_sp, $number_sp, $genome_fasta, $genome_id, %mapping_contigs, 
$scripts_path, %msa_all_taxa_files);

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
	"CL|conserved_length=i" => \$window_size_CONS, #20 ##size_conserved_region
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
my $msa_dirname = $align_dirname."/MSA";

## Markers dirname
my $marker_dirname = $folder_abs_path."/".$datestring."_DM_markers";
my $param_Detail_file_markers = $folder_abs_path."/".$datestring."_Markers-Parameters.txt";

my $start_time = $step_time = time;
my $pipeline_path = abs_path($0);
my @script_path_array = split ("/", $pipeline_path);
for (my $j = 0; $j < $#script_path_array; $j++) {
	$scripts_path .= $script_path_array[$j]."/";
} 
my $mapping_parameters_short = $folder_abs_path."/".$datestring."_tmp.txt";
open (MP_SHORT, ">$mapping_parameters_short");
## General variables
my $samtools_path = $scripts_path."samtools-1.3.1/samtools";
my $bowtie_path = $scripts_path."bowtie2-2.2.9/";
my $BLAST = $scripts_path."NCBI_BLAST/";
my $CAP3_exec = $scripts_path."cap3/bin/cap3";

## Check if a previous DOMINO parameters file exists
if (-e $param_Detail_file_markers) { File::Copy::move($param_Detail_file_markers, $param_Detail_file_markers."_old_".$random_number); }
if (-e $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $mapping_markers_errors_details."_old_".$random_number); }
if (-e $mapping_parameters) { File::Copy::move($mapping_parameters, $mapping_parameters."_old_".$random_number); }

## Behaviour DOMINO
if (!$behaviour) {
	&printError("\nPlease choose a development module for DOMINO between selection/discovery...\n"); DOMINO::dieNicely();
}
if ($behaviour eq 'selection') {
	$select_markers=1;
	if ($option eq "DOMINO_files" || $option eq "user_assembly_contigs" || $option eq "genome") {
		&printError("\nThe DOMINO development module SELECTION is not yet available for the option $option...\n"); DOMINO::dieNicely();
	} 
} elsif ($behaviour eq 'discovery') {
	if ($option eq "RADseq") {
		&printError("\nThe DOMINO development module DISCOVERY is not suitable for the option $option...\n"); DOMINO::dieNicely();
	}
	if (!$window_size_CONS || !$window_size_VARS_range) {
		unless ($option eq "msa_alignment") {
			&printError("\nMandatory options are missing...\n"); DOMINO::dieNicely();
	}}
	$identify_markers=1;
} else { &printError("\nPlease choose between selection/discovery...\n"); DOMINO::dieNicely(); }

# Control if missing options
if (!$variable_positions_user_range and !$variable_divergence) {
	&printError("Exiting the script. A range for variable positions or a minimum divergence is missing.\n Use the option -VP|--variable_positions [min::max] or -VD|--variable_divergence [float number]..\n"); DOMINO::dieNicely();
}
if (!$option) {
	&printError("\nPlease provide an option for DOMINO...\n"); DOMINO::dieNicely();
}
unless (!$variable_divergence) { if ($variable_divergence < 0) {$variable_divergence = 0.000000000000000000000000000000001;} ## Set a very small value if -VD 0 
} 
if (!$option) { $option = "DOMINO_files";}
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
	if ($option eq "msa_alignment" || $option eq "RADseq") {
		if ($avoid_mapping) {
			&printError("\nThe option -taxa_names option is missing...\n\nPlease provide it or DOMINO would not continue the process...\n"); DOMINO::dieNicely();
		}
	} else {
		&printError("\nThe option -taxa_names option is missing...\n\nPlease provide it or DOMINO would not continue the process...\n"); DOMINO::dieNicely();
}} else {
	$MID_taxa_names =~ s/\s+/\,/g;
	my @MID_name_array = split (",", $MID_taxa_names);
	for (my $j = 0; $j < scalar @MID_name_array; $j++) {
		$MID_species_hash{$MID_name_array[$j]} = $MID_name_array[$j];
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
	
	&debugger_print("DOMINO::get_earliest subroutine: "); &debugger_print($path_returned);
	
	if ($path_returned eq 'NO') {
		undef $avoid_mapping;
	} elsif ($path_returned eq 'mapping') {
		undef $avoid_mapping;
		DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
	} else {
		$align_dirname = $path_returned;
		my $time_folder;
		if ($align_dirname =~ /(\d+)\_DM\_mapping/) { $time_folder = $1; }
		my $parameters_mapping = $align_dirname."/".$time_folder."_tmp.txt";
		my %tmp_hash;
		my $flag_mapping = my $flag_mapping_poly = my $flag_mapping_local = 0; my $taxa = 0;
		if (-e -r -s $parameters_mapping) {
			open(PA_MAP, $parameters_mapping);
			while (<PA_MAP>) {
				my $line = $_;
				chomp $line;
				my @array_split = split(":", $line);
				my $value = $variables{$array_split[0]};
				if ($array_split[0] eq "poly") {
					if (!$polymorphism_user) {
						$flag_mapping++;
					} $flag_mapping_poly = 1;		
				} elsif ($array_split[0] eq "bowtie_local") {
					if (!$bowtie_local) {
						$flag_mapping++;
					} $flag_mapping_local = 1;
				} elsif ($array_split[0] eq "Mapped") {
					unless ($array_split[1] eq "GenomeID") {
						if ($MID_species_hash{$array_split[1]}) { $tmp_hash{$array_split[1]}++;
				}}} else {
					unless ($value == $array_split[1]) { $flag_mapping++;}
			}}
			close(PA_MAP);
			foreach my $keys (keys %MID_species_hash) { unless ($tmp_hash{$keys}) { $flag_mapping++; $taxa++; } }
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
			}
		} else {
			undef $avoid_mapping;
			DOMINO::printDetails("+ Generation of new profile of variation would be done as it has been previously done with different parameters or it was incomplete...OK\n", $mapping_parameters, $param_Detail_file_markers);
}}}

## Get ranges
my ($variable_positions_user_min, $variable_positions_user_max);
if ($variable_positions_user_range) {
	if ($variable_positions_user_range =~ m/.*\:\:.*/) {
		($variable_positions_user_min, $variable_positions_user_max) = split("::", $variable_positions_user_range);
	} elsif ($variable_positions_user_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 2::7\n\n"); DOMINO::dieNicely();
	} else {
		$variable_positions_user_min = $variable_positions_user_max = $variable_positions_user_range; 
}}
my ($window_size_VARS_min, $window_size_VARS_max);
if ($window_size_VARS_range) {
	if ($window_size_VARS_range =~ m/.*\:\:.*/) {
		($window_size_VARS_min, $window_size_VARS_max) = split("::", $window_size_VARS_range);
	} elsif ($window_size_VARS_range =~ m/.*\:(\d+)/) {
		&printError("\nPlease provide the range using 2 pair of dots like 200::700\n\n"); DOMINO::dieNicely();
	} else {
		$window_size_VARS_min = $window_size_VARS_max = $window_size_VARS_range;
}}

if (!$minimum_number_taxa_covered) {
	$minimum_number_taxa_covered = $number_sp;  ## Force to be all the taxa
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
		my $tmp = abs_path($user_contig_files[$i]);
		push (@contigs_fasta_file_abs_path, $tmp);
	}
	if (!$map_contig_files) {
		for (my $i = 0; $i < scalar @user_cleanRead_files; $i++) {
			my $tmp = abs_path($user_cleanRead_files[$i]);
			push (@clean_fastq_file_abs_path, $tmp);
	}}
} elsif ($option eq 'DOMINO_files') {
	## Get reference contig files (FASTA & FASTQ) and push them into @contigs_fasta_files
	&get_MID_contigs();
	
	## Obtain clean reads
	if (scalar @user_cleanRead_files == 0) {
		## use clean reads to map
		&get_clean_files();
	} else {
		for (my $i = 0; $i < scalar @user_cleanRead_files; $i++) {
			my $tmp = abs_path($user_cleanRead_files[$i]);
			push (@clean_fastq_file_abs_path, $tmp);
	}}
} elsif ($option eq 'genome') {
	my $tmp = abs_path($genome_fasta);
	push (@contigs_fasta_file_abs_path, $tmp);
	if (scalar @user_cleanRead_files == 0) {
		&printError("Clean Read files were not provided...\nDOMINO would check in the output folder provided if there is a DOMINO_clean_data containing the FASTQ files for each taxa...."); 
		&get_clean_files();
		if (scalar @clean_fastq_file_abs_path == 0) {
			&printError("Clean Read files were not provided and DOMINO_clean_data folder seems to be missing...\nPlease provide a FASTQ files for each taxa using the -user_cleanRead_files option...."); DOMINO::dieNicely();
	}} else {
		for (my $i = 0; $i < scalar @user_cleanRead_files; $i++) {
			my $tmp = abs_path($user_cleanRead_files[$i]);
			push (@clean_fastq_file_abs_path, $tmp);
	}}
} elsif ($option eq "msa_alignment") {	
	if ($radseq_like_data) {
		$msa_file_abs_path = abs_path($msa_file);
		if ($pyRAD_file) { print "+ pyRAD data file: $msa_file_abs_path\n";
		} elsif ($stacks_file) { print "+ STACKS file: $msa_file_abs_path\n";}		
		print "+ Checking file:\n";
		if (-f -e -r -s $msa_file_abs_path) {
			print "\tFile $msa_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
 			chdir $align_dirname; system("ln -s $msa_file_abs_path");
		} else { 
			&printError("File provided is not valid...\nPlease provide a valid file as specified in the DOMINO manual...."); DOMINO::dieNicely();
	}} elsif ($msa_file) {
		$msa_file_abs_path = abs_path($msa_file);
		print "+ Multipe sequence alignment file provided: $msa_file_abs_path\n";		
		print "+ Checking file:\n";
		if (-f -e -r -s $msa_file_abs_path) {
			print "\tFile $msa_file_abs_path\n\t\tFile exists, is readable and non-zero character...OK\n";
 			chdir $align_dirname; system("ln -s $msa_file_abs_path");
		} else { 
			&printError("MSA file provided is not valid...\nPlease provide a valid contig MSA file as specified in the DOMINO manual...."); DOMINO::dieNicely();
	}} elsif ($msa_fasta_folder) {
		$msa_folder_abs_path = abs_path($msa_fasta_folder);
		print "+ Multipe sequence alignment fasta folder provided: $msa_folder_abs_path\n+ Checking file(s):\n";
		if (-d $msa_folder_abs_path) {
			print "\tFolder $msa_folder_abs_path\n\t\tFolder exists, is readable and non-zero character...OK\n";
 			chdir $align_dirname; system("ln -s $msa_folder_abs_path");
			my @array = split("/", $msa_folder_abs_path);
			my $array_files_fasta_msa_ref = DOMINO::readDir($msa_folder_abs_path);
			@array_files_fasta_msa = @$array_files_fasta_msa_ref;
			print "\t Checking files in folder provided...\n";
			for (my $i=0; $i < scalar @array_files_fasta_msa; $i++) {
				if ($array_files_fasta_msa[$i] eq "." || $array_files_fasta_msa[$i] eq ".." || $array_files_fasta_msa[$i] eq ".DS_Store"  ) {next;}
				my $file_path = $msa_folder_abs_path."/".$array_files_fasta_msa[$i];
				unless (-f -e -r -s $file_path) {
					&printError("File $file_path is not readable or empty. Please discarded from the folder...\n"); DOMINO::dieNicely();
			}} print "\t\t Files checked and everything seems OK..."
		} else { &printError("MSA folder provided is not valid...\n"); DOMINO::dieNicely(); }
	} else { &printError("MSA folder or file is missing...\n"); DOMINO::dieNicely(); }
}
unless ($option eq "msa_alignment") {
	## We would check the files for DOMINO detection of molecular markers
	if (scalar @contigs_fasta_file_abs_path == 0) {
		&printError("Contig assembled FASTA files were not able to be obtained...\nPlease keep DOMINO original names or re-run DOMINO Assembly step..."); DOMINO::dieNicely();
	} else {
		for (my $i = 0; $i < scalar @contigs_fasta_file_abs_path; $i++) {
			&check_file($contigs_fasta_file_abs_path[$i]);
}}}

unless (!$MID_taxa_names) {
	DOMINO::printDetails("\n\n+ Taxa to use for the DOMINO development of molecular markers:\n", $mapping_parameters, $param_Detail_file_markers);
	foreach my $keys (keys %MID_species_hash) { DOMINO::printDetails("\tName: $MID_species_hash{$keys}\n", $mapping_parameters, $param_Detail_file_markers); }
}
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
	DOMINO::printDetails("\t- Conserved Length (CL): $window_size_CONS (bp)\n", $param_Detail_file_markers);
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
	} else {
		DOMINO::printDetails("\t- Variable Positions (VP): $variable_positions_user_min -- $variable_positions_user_max (bp)\n", $param_Detail_file_markers);	
}}

## Common markers parameters
unless (!$MID_taxa_names) { DOMINO::printDetails("\t- Minimum number of covered taxa (MCT): ".$minimum_number_taxa_covered."\n", $param_Detail_file_markers); }
if ($polymorphism_user) { 
	DOMINO::printDetails("\t- Polymorphic variants would be detected (PV)...OK\n", $param_Detail_file_markers);
	print MP_SHORT "poly:1\n";
}
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
## Using simulated data we show the although Mugsy can align multiple taxa
## quite efficiently, it works really bad if low coverage data provided.
## Depending on the option provided, DOMINO would map reads against a reference (assembled or provided)
## or parse and alignment file provided.

if ($option ne "msa_alignment" and !$avoid_mapping) {	## We would use Bowtie2 for mapping the reads		
	DOMINO::printHeader("", "#");	DOMINO::printHeader(" Mapping Process started ", "#"); DOMINO::printHeader("", "#"); print "\n";
	
	#############################################################
	### Get Pre-assemble taxa read contigs of each taxa ### 
	#############################################################
	print "\n"; DOMINO::printHeader(" Get FASTQ files of the contigs generated ", "%"); print "\n";

	## Mapping of the reads, all taxa used as reference
	for (my $j = 0; $j < scalar @contigs_fasta_file_abs_path; $j++) {
		chdir $align_dirname;
		my (@sam_files, @clean_sam_files, @sorted_bam, $reference_identifier); 
			
		## Get the name and identifier of each reference fasta used	
		my @temp_contigs_name = split ("/", $contigs_fasta_file_abs_path[$j]);
		my $contigs_fasta = $temp_contigs_name[$#temp_contigs_name];

		## Check if reference file is another one containing not the names provided
		#foreach my $keys (keys %MID_species_hash) { if ($contigs_fasta =~ /.*($MID_species_hash{$keys}).*/g) { $reference_identifier = $1; } }
		if ($contigs_fasta =~ /.*id\-(.*)\.contigs\.fasta/) {
			$reference_identifier = $1;
			if ($reference_identifier =~ /(.*)\_R\d+/) {
				$reference_identifier = $1;
		}} elsif ($genome_fasta) {
			if ($contigs_fasta =~ /.*id\-(.*)\.fasta/) {
				$genome_id = $1; $reference_identifier = $1;			
			} else {
				&printError("Genome fasta file erroneously tagged...\n\nPlease provide a genome file tag as '[xx]id-[yyy].fasta'. Where:
'[xx]' any character (or none). Please avoid using dots (.)\n'[yyy]' taxon identifier of the reference genome....."); 
				DOMINO::dieNicely();
		}}
		
		my $ref_Fasta = $reference_identifier.".fasta";
		system("ln -s $contigs_fasta_file_abs_path[$j] $ref_Fasta");
		
		## Generate a directory for each one
		my $dir = $align_dirname."/".$reference_identifier; mkdir $dir, 0755; chdir $dir;
		print "+ Using as reference: $contigs_fasta\tID: $reference_identifier...OK\n";
		print "+ Generating a new directory $dir....OK\n\n";
			
		#######################################
		###		 Copy necessary files		### 
		#######################################
		print "+ Copying necessary files\n";
		system("ln -s $contigs_fasta_file_abs_path[$j]");
		#copy ($contigs_fasta_file_abs_path[$j], $dir."/".$contigs_fasta);
		
		my @clean_fastq_files;
		my @tmp_array;
		if (scalar @user_cleanRead_files > 0) { ## Already pushed into this array in line 403
			print "+ User clean reads files would be mapped\n+ Checking files and format\n";
			@tmp_array = @clean_fastq_file_abs_path;
		} elsif ($map_contig_files) {
			print "+ Contig files would be mapped\n+ Checking files and format\n";
			@tmp_array = @contigs_fasta_file_abs_path;
		} else { ## Map DOMINO clean reads
			print "+ Clean reads files would be mapped\n+ Checking files and format\n";
			@tmp_array = @clean_fastq_file_abs_path;
		}
		
		for (my $i = 0; $i < scalar @tmp_array; $i++) {
			my @temp_clean_Fastq_name = split ("/", $tmp_array[$i]);
			my $clean_fastq = $temp_clean_Fastq_name[-1];
			#copy ($tmp_array[$i], $clean_fastq_abs_path);
			system("ln -s $tmp_array[$i]");
			push (@clean_fastq_files, $clean_fastq);
			&check_file($clean_fastq);
		} print "Done...\n\n";
					
		###################################
		###		 Index Contig file		### 
		###################################
		DOMINO::printHeader(" Indexing Contig File for mapping Reference ", "%");
		# Index contig reference file using Bowtie
		my $reference = "reference_".$reference_identifier;
		print "- Reference: $contigs_fasta...\n";
		my $bowtie_index_call = $bowtie_path."bowtie2-build --threads $num_proc_user -f ".$contigs_fasta." ".$reference;   
		&debugger_print("BOWTIE2 command: ".$bowtie_index_call."\n");
		my $index_result = system ($bowtie_index_call);
		if ($index_result != 0) {
			&printError("Exiting the script. Some error happened when calling bowtie for indexing the file...\n"); DOMINO::dieNicely();
		} print "\n"; &time_log();	print "\n";

		###########################
		###	Align taxa Reads	### 
		###########################
		print "\n";	DOMINO::printHeader(" Aligning Reads Individually ", "%"); print "\n";
		my @mapping_files = @clean_fastq_files;
		my %files_used;
		for (my $i = 0; $i < scalar @mapping_files; $i++) {
			if (!$files_used{$mapping_files[$i]}) { 
				$files_used{$mapping_files[$i]}++;
			} else {
				next;
			}
			my $MID_key;
			foreach my $keys (keys %MID_species_hash) {
				if ($mapping_files[$i] =~ /.*($MID_species_hash{$keys}).*/g) {
					$MID_key = $1;	
					last;
			}}
			my $sam_name = "taxa_".$MID_species_hash{$MID_key}.".sam";
			push (@sam_files, $sam_name);
	
			## Map reads using bowtie
			my $R_group_id = '--rg-id '.$MID_species_hash{$MID_key};
			my $R_group_name = ' --rg '.$MID_key;
			my $threads = ' -p '.$num_proc_user;
			my $mismatches = ' -N 1 --np 0'; ## Do not add penalty if read/ref got an ambiguous base
			my $read_gap_open = ' --rdg '.$rdgopen.','.$rdgexten;
			my $ref_gap_open = ' --rfg '.$rfgopen.','.$rfgexten;
			my $mismatch_penalty = ' --mp '.$mis_penalty;
			my $mapping_file = $mapping_files[$i];
			my $botwie_system = $bowtie_path."bowtie2";
			if ($bowtie_local) { $botwie_system .= " --local"; }

			if ($input_type eq 'pair_end') {
				my ($name, $pair, $second_Read_file);
				for (my $k = 0; $k < scalar @mapping_files; $k++) { 				## Obtain mate for mapping file
					my $tmp = $mapping_file;
					$second_Read_file = $mapping_files[$k];
					if ($second_Read_file eq $mapping_file) { next;
					} elsif ($mapping_file =~ /.*id\-(.*)\_R(\d+).*/) { ##  reads.id-sp1_R1.fastq, reads.id-sp1_R2.fastq
						$name = $1; $pair = $2;
						if ($pair == 1) { $tmp =~ s/R1/R2/g;
						} else { $tmp =~ s/R2/R1/g; 
					}}					
					if ($tmp eq $second_Read_file) { 
						$files_used{$second_Read_file}++;
						last; 
				}}
				print "Aligning reads for $mapping_file and $second_Read_file ...\n";
				$botwie_system .= " -x ".$reference." -q -1 $mapping_file -2 $second_Read_file -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
			} elsif ($input_type eq 'single_end') { ## Illumin single end, 454
				print "Aligning reads for $mapping_file file...\n";
				if ($map_contig_files) { ## Mapping contigs
					$botwie_system .= " -x ".$reference." -f -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
				} else {
					$botwie_system .= " -x ".$reference." -q -U ".$mapping_file." -S ".$sam_name." ".$R_group_id.$R_group_name.$threads.$mismatches." --no-unal".$read_gap_open.$ref_gap_open.$mismatch_penalty;   
			}} 
			&debugger_print("BOWTIE2 command: ".$botwie_system."\n"); 
			my $system_bowtie_call = system ($botwie_system);
			if ($system_bowtie_call != 0) {
				&printError("Exiting the script. Some error happened when calling bowtie for mapping the file $mapping_file...\n"); DOMINO::dieNicely();
			} else { print "\n"; }
			print "\n"; &time_log(); print "\n";
		}			
		print "\n";	DOMINO::printHeader("", "#"); DOMINO::printHeader(" Mapping finished ", "#"); DOMINO::printHeader("", "#"); print "\n";
		&time_log(); print "\n";
		
		###################################
		###	Remove multimapping reads	### 
		###################################
		print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Cleaning the Mapping files generated ", "#"); DOMINO::printHeader("", "#");	
		print "\nCleaning reads now...\n";
		%mapping_contigs = ();
		for (my $i = 0; $i < scalar @sam_files; $i++) {
			print "\nChecking mapping reads in ".$sam_files[$i]."...\n";
			push (@clean_sam_files, &discard_reads_sam($sam_files[$i]));	
			## Check reads, print clean SAM file and push the name into array
		} print "\n"; &time_log(); print "\n";
		
		#################################
		## Get the the reference fasta ##
		#################################
		print "\n"; DOMINO::printHeader("", "%"); DOMINO::printHeader(" Obtain information of the Reference sequence ", "%"); DOMINO::printHeader("", "%"); print "\n";
		print "- Reading the reference fasta file...\n";
		my ($reference_hash_fasta_ref, $message) = DOMINO::readFASTA_hashLength($contigs_fasta); ## Obtain reference of a hash
		
		###################################################################
		## Generate sorted bam files in order to be able to get coverage ##
		###################################################################
		for (my $i = 0; $i < scalar @clean_sam_files; $i++) {
			push (@sorted_bam, &generate_bam($clean_sam_files[$i]));
		} print "\n"; &time_log(); print "\n";

		###################################
		###		 Index Contig file		### 
		###################################
		print "- Indexing the reference fasta file $contigs_fasta...\n";
		## Index contig file using samtools faidx
		my $samtools_index_system = $samtools_path." faidx ".$contigs_fasta;
		&debugger_print("SAMTOOLS command: ".$samtools_index_system."\n");
		my $samtools_index_system_call = system($samtools_index_system);
		if ($samtools_index_system_call != 0) {
			&printError("Exiting the script. Some error happened when calling SAMtools for indexing the file $contigs_fasta...\n"); DOMINO::dieNicely();
		}
		print "\n"; &time_log(); print "\n";

		DOMINO::printHeader("", "#"); DOMINO::printHeader(" Filter Contigs according to Coverage ", "#"); DOMINO::printHeader("", "#"); print "\n";
		for (my $i = 0; $i < scalar @sorted_bam; $i++) {
			%discard_contigs = (); %max_cov = (); ## Initialize some hashes
			print "- Generating coverage statistics for $sorted_bam[$i]\n- Obtaining the coverage of each pair base...\n";
			
			## Generate Coverage statistics for the alignment file
			my @tmp_bam_name = split ("\.sorted.bam", $sorted_bam[$i]);
			my $coverage_file = $tmp_bam_name[0]."_coverage_stats.txt";
			my $coverage_samtools_command = $samtools_path." depth ".$sorted_bam[$i]." > ".$coverage_file;
			&debugger_print("SAMTOOLS command: $coverage_samtools_command\n");
			my $system_coverage_call = system ($coverage_samtools_command);
			if ($system_coverage_call != 0) {
				&printError("Exiting the script. Some error happened when calling SAMtools for obtaining coverage of file $sorted_bam[$i]...\n"); DOMINO::dieNicely();
			} else {  print "\n"; }
			
			## Filter the file according to coverage					
			my $sam_filtered_returned = &filter_coverage_Adjust($coverage_file, $sorted_bam[$i]);
			my $bam_filtered_returned = &generate_bam($sam_filtered_returned);
			unless ($tmp_bam_name[0] =~ /.*$reference_identifier.*/) { ## DO NOT GENERATE FILTER PROFILE FOR REFERENCE
				&generate_filter_PILEUP($bam_filtered_returned, $contigs_fasta, $reference_hash_fasta_ref, $reference_identifier);
			}
			undef %discard_contigs; undef %max_cov; ## Initialize some hashes
			print "\n"; &time_log();	print "\n";
		}
		
		undef %mapping_contigs;
		unless ($avoidDelete_tmp_files) {
			############################
			## Delete Temporary Files ##
			############################
			DOMINO::printHeader(" Deleting Temporary Files of the Aligment Folder ", "%");
			&delete_files_mapping($dir, $reference_identifier);
		}	
		print "\n\n"; DOMINO::printHeader("", "+");DOMINO::printHeader(" Mapping finished ", "+");DOMINO::printHeader("", "+");print "Done...\n\n"; &time_log(); print "\n";		
		
		if ($genome_fasta) {
			print MP_SHORT "Mapped:GenomeID:".$reference_identifier."\n";
		} else {
			print MP_SHORT "Mapped:".$reference_identifier."\n";
		}
	}
	close(MP_SHORT);
} elsif ($option eq "msa_alignment" && !$avoid_mapping) {

	##########################################################################################
	###################### Check Alignment file/folder #######################################
	##########################################################################################

	my @species_alignment;
	## Get species to check and print into taxa_names.txt
	if ($radseq_like_data) {
		my $file = $msa_file_abs_path;
		if ($pyRAD_file) {
			## Parse pyRAD loci file provided
			my $counter = 1;
			open(FILE, $file) || die "Could not open the $file ...\n";
			my @name = split(".loci", $file);
			while (<FILE>) {		
				next if /^#/ || /^\s*$/;
				my $line = $_;
				chomp $line;
				if ($line =~ /\/\//) { next; }
				$line =~ s/\s+/\t/g; $line =~ s/\%/>/g;
				my @array = split("\t", $line);
				$array[0] =~ s/\>//;
				unless (grep /$array[0]/, @species_alignment) {
					push(@species_alignment, $array[0])							
			}} close (FILE);
		} elsif ($stacks_file) {
			
			## Parse STACKS file provided
			my (%hash, $new_id);
			open(FILE, $file) || die "Could not open the $file ...\n";
			my @name = split(".fa", $file);
			$/ = ">"; ## Telling perl where a new line starts
			while (<FILE>) {		
				next if /^#/ || /^\s*$/;
				chomp;
				my ($titleline, $sequence) = split(/\n/,$_,2);
				next unless ($sequence && $titleline);
				chomp $sequence;
				$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
				$titleline =~ s/\r//g;
				if ($titleline =~ /CLocus\_(\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
					unless (grep /$2/, @species_alignment) {
						push(@species_alignment, $2)							
			}}} close(FILE); $/ = "\n";
		} print "\n"; &time_log(); print "\n";
	} else {
		### MSA file or folder provided
		mkdir $msa_dirname, 0755;
		if ($msa_file) {
			## Check the alignment format provided...
			print "\n\n"; DOMINO::printHeader(" Checking Alignment file provided ", "%");
			my $species_alignment_ref = &read_phylip_aln($msa_file_abs_path);
			@species_alignment = @$species_alignment_ref;
			print "\n"; &time_log(); print "\n";
	
		} elsif ($msa_fasta_folder) {		
			print "\n\n"; DOMINO::printHeader(" Checking Alignment folder provided ", "%");
			print "+ Checking files...\n";
			my $files = 0;
			for (my $i = 0; $i < scalar @array_files_fasta_msa; $i++) {
				if ($array_files_fasta_msa[$i] eq "." || $array_files_fasta_msa[$i] eq ".." || $array_files_fasta_msa[$i] eq ".DS_Store"  ) {next;}
				my $file_path = $msa_folder_abs_path."/".$array_files_fasta_msa[$i];
				my %alignment;			
				if ($array_files_fasta_msa[$i] =~ /(.*)\.fasta/) {
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
						unless (grep /$titleline/, @species_alignment) {push(@species_alignment, $titleline);}
					} close(FILE); $/ = "\n";				
					system("ln -s $file_path $msa_dirname"); $files++;			
				} else {
					my $species_alignment_ref = &read_phylip_aln($msa_file, $array_files_fasta_msa[$i]);
					push (@species_alignment, @$species_alignment_ref);
			}}
			print "+ Parsing of the $files files has been done...\n"; print "\n"; &time_log(); print "\n";
	}}
	print "\n"; &time_log(); print "\n";
	unless ($MID_taxa_names) {
		#print into taxa_msa.txt the names of the taxa find
		my @species_alignment_sort = sort @species_alignment;
		my @species_alignment_uniq = uniq(@species_alignment_sort);
		my $out_file = $folder_abs_path."/taxa_names.txt";
		my $string;
		for (my $i=0; $i<scalar @species_alignment_uniq; $i++) {
			$string .= $species_alignment_uniq[$i].",";
		}
		
		if ($select_markers) {
			chdir $folder_abs_path; &debugger_print("Changing dir to $folder_abs_path");
			remove_tree($align_dirname);
			remove_tree($marker_dirname);
			remove_tree($mapping_parameters); 
			remove_tree($mapping_markers_errors_details); 
			remove_tree($param_Detail_file_markers);
		}
		chop $string; open (OUT, ">$out_file"); print OUT $string."\n"; close(OUT);
		print "\n\nExiting DOMINO...\n";
		&finish_time_stamp(); print "\n\n Job done succesfully, exiting the script\n\n\n"; 
		print "\n\n"; DOMINO::printHeader("","#"); print "NOTE:\n\n";
		print "\t+ No taxa names were provided so DOMINO have parsed and checked for the names...\n";	
		print "\t+ Please checked the file taxa_names.txt within the main project folder generated and\n\tprovide the taxa of interest using option -taxa_names with comma-separated values...\n";
		print "\t+ Re-run DOMINO providing the option -NPG for reducing the computational time...\n\n";
		DOMINO::printHeader("","#"); print "\n\n"; exit();
	}
	
} else { ## No_Profile_Generation|NPG: just get files
	########################################################################################
	###################### NoMapping: just get files #######################################
	########################################################################################
	if ($option eq 'DOMINO_files') {
		print "\n"; DOMINO::printHeader("", "+");DOMINO::printHeader(" Mapping has been avoided ", "+");
		DOMINO::printHeader("", "+"); print "+ Contig FASTA files would be obtained...\n\n";
		&get_MID_contigs(); &time_log(); print "\n";			
	} elsif ($option eq 'genome') {
		if ($contigs_fasta_file_abs_path[0] =~ /.*id\-(.*)\.fasta/) {
			$genome_id = $1;
}}}

## Move parameters and error file to folder
File::Copy::move($mapping_parameters, $align_dirname);
File::Copy::move($mapping_parameters_short, $align_dirname);
unless (-z $mapping_markers_errors_details) { File::Copy::move($mapping_markers_errors_details, $align_dirname); }

##########################################################################################
###################### MARKER DEVELOPMENT ################################################
##########################################################################################

### MSA alignment
if ($option eq "msa_alignment") {

	if ($pyRAD_file) {
		### RADSEQ like data ####
		my $file = $msa_file_abs_path;
		mkdir $msa_dirname, 0755;
		print "\n\n"; DOMINO::printHeader(" Checking pyRAD loci file provided ", "%");
		print "+ For each loci a MSA file would be generated...\n";
	
		## Parse pyRAD loci file provided
		my %hash;
		my $counter = 1;
		open(FILE, $file) || die "Could not open the $file ...\n";
		my @file_name = split("/", $file);
		my @name = split(".loci", $file_name[-1]);
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			my $line = $_;
			chomp $line;
			if ($line =~ /\/\//) {
				my $file = $msa_dirname."/".$name[0]."_loci_".$counter.".fasta";
				foreach my $keys (keys %hash) { 
					if ($MID_species_hash{$keys}) { 
						open (OUT, ">>$file");   		
						print OUT ">".$keys."\n".$hash{$keys}."\n";
						close(OUT);
				}}
				$counter++; undef %hash; next;
			}
			$line =~ s/\s+/\t/g; $line =~ s/\%/>/g;
			my @array = split("\t", $line);
			$array[0] =~ s/\>//;
			$hash{$array[0]} = $array[1];
		}
		close(FILE);
		undef %hash;
		$msa_fasta_folder = $msa_dirname;
		$msa_folder_abs_path = $msa_dirname;
		undef $msa_file;

	} elsif ($stacks_file) {
		### RADSEQ like data ####
		my $file = $msa_file_abs_path;
		mkdir $msa_dirname, 0755;
		print "\n\n"; DOMINO::printHeader(" Checking STACKS file provided ", "%");
		print "+ For each loci a MSA file would be generated...\n";
		
		## Parse STACKS file provided
		my (%hash, $new_id);
		open(FILE, $file) || die "Could not open the $file ...\n";
		my @file_name = split("/", $file);
		my @name = split(".fa", $file_name[-1]);
		$/ = ">"; ## Telling perl where a new line starts
		while (<FILE>) {		
			next if /^#/ || /^\s*$/;
			chomp;
			my ($titleline, $sequence) = split(/\n/,$_,2);
			next unless ($sequence && $titleline);
			chomp $sequence;
			$sequence =~ s/\s+//g; $sequence =~ s/\r//g;
			$titleline =~ s/\r//g;
			if ($titleline =~ /CLocus\_(\d+)\_(.*)\_Locus\_(\d+)\_Allele\_(\d+).*/){
				if (!$new_id) { $new_id = $1;}
				if ($new_id ne $1 ) {
					&parse_stacks_marker(\%hash, $new_id, $name[0]);
					undef %hash; $new_id = "";
				}
				push (@{$hash{$2}}, $sequence);
		}}
		close(FILE); $/ = "\n";
		&parse_stacks_marker(\%hash, $new_id, $name[0]);
		undef %hash;			
		$msa_fasta_folder = $msa_dirname;
		$msa_folder_abs_path = $msa_dirname;
		undef $msa_file;
	} else {
		$msa_fasta_folder = $msa_dirname;
		$msa_folder_abs_path = $msa_dirname;
	} print "\n"; &time_log(); print "\n";

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
foreach my $keys_hash (sort keys %MID_species_hash) { ## For each taxa specified, obtain putative molecular markers
	if ($genome_marker_bool == 1) {last;}
	
	print "\n";
	## Create a dir for each taxa
	DOMINO::printHeader(" Checking taxa files user specified ", "#"); 
	my ($dir, $ref_taxa);
	if ($genome_fasta) {
		$ref_taxa = $genome_id;
		print "Checking: \tGenome provided: $ref_taxa\n\n";
		$genome_marker_bool = 1;
		$dir = $marker_dirname."/DOMINO_markers_Genome";
	} elsif ($option eq "msa_alignment") {
		$ref_taxa = $keys_hash;
		print "Checking: \t$ref_taxa\n\n";
		$dir = $marker_dirname."/markers_Ref_".$ref_taxa;
	} else {
		$ref_taxa = $keys_hash;
		print "Checking: \t$ref_taxa\n\n";
		$dir = $marker_dirname."/markers_Ref_".$ref_taxa;
	}	
	mkdir $dir, 0755; chdir $dir; &debugger_print("Changing dir to $dir");
	
	## Initialize some variables
	undef %coord_contig; undef %coord_markers; undef %putative_markers; undef $merge_bam_all_sp;
	my (%pileup_files, $array_all_species, $contigs_fasta, %ref_fasta, @sams_to_parse, @sam_headers_line);

	#######################################
	###		 Copy necessary files		### 
	#######################################
	print "+ Copying necessary files\n";
	my ($reference_bam_file, $ref_arrays_pileups, $ref_clean_sam) = &get_filtered_sam_files($ref_taxa, $dir); ## filtered_sam_files and PILEUP arrays and $reference_bam_file
	my @pileup_Arrays = @$ref_arrays_pileups;
	my @clean_filtered_sam_files = @$ref_clean_sam;
	if ($genome_fasta) {$reference_bam_file = $ref_taxa.".sorted.bam";}

	unless ($option eq "msa_alignment") { print "\n";
		for (my $j = 0; $j < scalar @clean_filtered_sam_files; $j++) {
			my @temp = split ("/", $clean_filtered_sam_files[$j]);
			my $sam_file = $temp[$#temp];
			print "Checking sam file: ".$sam_file."...\n";
			push (@sams_to_parse, $clean_filtered_sam_files[$j]);
			my $array_ref = &get_headers_sam($sam_file);	## Get the header information of each file
			push (@sam_headers_line, @$array_ref);
		}

		#################################################
		## Merging the sam according to the user input ##
		#################################################
		print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Merging the SAM files according to user input ", "#"); DOMINO::printHeader("", "#"); print "\n";
		$merge_bam_all_sp = &merge_sam(\@sams_to_parse, \@sam_headers_line);
		print "\n"; &time_log(); print "\n";
		my @tmp = split ("_sorted\.bam", $merge_bam_all_sp);
		$array_all_species = $tmp[0].".profile_ARRAY.txt";	
	} else {
		$array_all_species = "./merged.profile_ARRAY.txt";
		my $array_contigs_aln = DOMINO::readDir($align_dirname);
		my @array_contigs = @$array_contigs_aln;
		for (my $i=0; $i < scalar @array_contigs; $i++) {
			if ($array_contigs[$i] eq "." || $array_contigs[$i] eq ".." || $array_contigs[$i] eq ".DS_Store" ) {next;}
			if ($array_contigs[$i] =~ /.*fasta/) { 
				my $file = $align_dirname."/".$array_contigs[$i];
				push (@contigs_fasta_file_abs_path, $file); 
	}}}
	
	##################################
	## Get the the reference fasta  ##
	##################################
	print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Obtaining information of the Reference sequence file provided ", "#"); DOMINO::printHeader("", "#");

	## Remove duplicates if any
	my @tmp = @contigs_fasta_file_abs_path;
	my @tmp_sort = sort(@tmp);
	my @tmp_sort_uniq = uniq(@tmp_sort);
	@contigs_fasta_file_abs_path = @tmp_sort_uniq;

	## Get the apropiate contig fasta file
	for (my $k = 0; $k < scalar @contigs_fasta_file_abs_path; $k++) {
		my @temp_contigs_name_1 = split ("/",$contigs_fasta_file_abs_path[$k]);
		my $contigs_fasta_name_1 = $temp_contigs_name_1[-1];
		if ($contigs_fasta_name_1 =~ /.*$ref_taxa.*/) {
			$contigs_fasta = $contigs_fasta_file_abs_path[$k];
			last;	
	}}	
	my @temp_contigs_name = split ("/", $contigs_fasta);
	my $contigs_fasta_name = $temp_contigs_name[-1];
	my $contigs_fasta_abs_path = $dir."/".$contigs_fasta_name;
	
	print "- Checking the file specified as reference fasta... $contigs_fasta_name...OK\n";
	system("ln -s $contigs_fasta");
	
	############################################################
	###	Generate an array information file for the reference ### 
	############################################################
	print "- Reading the reference fasta file...\n";
	open(FILE, $contigs_fasta_abs_path) || die "Could not open the $contigs_fasta_abs_path ...\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
		my ($titleline, $sequence) = split(/\n/,$_,2);
		next unless ($sequence && $titleline);
		$sequence =~ s/\n//g;
		my $tmp_size = length($sequence);
		my $size = $tmp_size - 1;
		$titleline =~ s/ /\t/g;
		my @array_title = split("\t", $titleline);
		push (@{ $ref_fasta{$array_title[0]} }, $size );
	}
	close(FILE);
	$/ = "\n";
	print "\n"; &time_log(); print "\n";
	
	##########################################
	## 	Merge PILEUP information arrays     ##
	##########################################
	print "\n"; DOMINO::printHeader(" Fetching information from all the PROFILEs generated ", "#");
	
	##########################################
	## Check conserved and variable regions ##
	##########################################
	for (my $i = 0; $i < scalar @pileup_Arrays; $i++) {
		my $array_files_ref = DOMINO::readDir($pileup_Arrays[$i]);
		my @array_files = @$array_files_ref;
		my $name;
		if ($pileup_Arrays[$i] =~ /.*\_taxa\_(.*)\_clean.*/) { $name = $1; }
		for (my $y = 0; $y < scalar @array_files; $y++) {
			if ($array_files[$y] eq "." || $array_files[$y] eq ".." || $array_files[$y] eq ".DS_Store") {next;}
			if ($array_files[$y] =~ /.*sequence\.fasta/) { next; } ## Maybe save path for sequence fasta file?
			my $tmp = "$pileup_Arrays[$i]/$array_files[$y]";
			$pileup_files{$array_files[$y]}{$name} = $tmp;
	}}

	print "- Merging variable and conserved information into a unique array...\n- Obtaining information from contigs\n";
	foreach my $contigs (keys %pileup_files) {
		my @pileup_fasta;
		foreach my $files (keys %{ $pileup_files{$contigs} }) {
			my $tmp_hash_reference = DOMINO::readFASTA_hash($pileup_files{$contigs}{$files});
			my %tmp_fasta = %{$tmp_hash_reference};
			foreach my $seqs (keys %tmp_fasta) {
				push (@pileup_fasta, $tmp_fasta{$seqs});
		}}
		## Merging variable and conserved information into a unique array...
		my $contig_name;
		if ($contigs =~ /.*\_ARRAY.*/) {
			my @contig_name = split("\_ARRAY", $contigs);
			$contig_name = $contig_name[0];
		} else { $contig_name = $contigs; }
		
		my $size = $ref_fasta{$contig_name}[0];
		my $tmp_string;			
		for (my $i = 0; $i <= scalar $size; $i++) {
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
		open(OUT_PILEUP, ">>$array_all_species");
		my $var_sites = $tmp_string =~ tr/1/1/; ## count variable sites
		my $cons_sites = $tmp_string =~ tr/0/0/; ## count conserved sites
		if ($var_sites != 0 && $cons_sites != 0) { 
			print OUT_PILEUP ">$contig_name\n$tmp_string\n"; 
		} 
		close (OUT_PILEUP);			
	}		
	print "\n"; &time_log(); print "\n"; DOMINO::printHeader(" Filtering the merged PROFILE array generated ", "%");
	my $folder = abs_path();
	my $tmp_PILEUP_merge_folder = "PROFILE_merge_species";
	my $PILEUP_merged_folder_abs_path = $folder."/".$tmp_PILEUP_merge_folder;
	mkdir $tmp_PILEUP_merge_folder, 0755; chdir $tmp_PILEUP_merge_folder; 
	&debugger_print("Changing dir to $tmp_PILEUP_merge_folder");
	my $array_all_species_abs_path = $folder."/".$array_all_species;
	open (ARRAY,"<$array_all_species_abs_path");
	$/ = ">"; ## Telling perl where a new line starts
	while (<ARRAY>) {		
		next if /^#/ || /^\s*$/; chomp;
		my ($titleline, $sequence) = split(/\n/,$_,2);
		next unless ($sequence && $titleline);
		chomp $sequence;
		$sequence =~ s/\n//g;
		my @array_name = split(/\s+/,$titleline); 
		my $contig = $array_name[0];
		my $file = $contig."_ARRAY.txt";
		open (FH, ">>$file"); print FH ">".$titleline."\n".$sequence."\n"; close(FH);
	}
	close(ARRAY); $/ = "\n";
	&sliding_window_conserve_variable($PILEUP_merged_folder_abs_path);
	print "\n"; &time_log(); print "\n";
	chdir $folder; &debugger_print("Changing dir to $folder");	
	my $ref_array_folder_files = DOMINO::readDir($tmp_PILEUP_merge_folder);
	my @array_folder_files = @$ref_array_folder_files;

	######################################################################
	## Check the coordinates of each taxa against the merge statistics  ##
	######################################################################
	print "\n"; DOMINO::printHeader(" Checking each taxa coordinates ", "#");
	for (my $contigs_merge =0; $contigs_merge < scalar @array_folder_files; $contigs_merge++) {
		if ($array_folder_files[$contigs_merge] eq '.' || $array_folder_files[$contigs_merge] eq '..' || $array_folder_files[$contigs_merge] eq '.DS_Store'  ) { next; }		
		my $array_file;
		if ($array_folder_files[$contigs_merge] =~ /(.*\_ARRAY\.txt)/) {
			$array_file = $1;
			$pileup_files{$array_file}{"merged"} = $PILEUP_merged_folder_abs_path."/".$array_folder_files[$contigs_merge];
		} elsif ($array_folder_files[$contigs_merge] =~ /(.*\_ARRAY)\-V(D|P).*/) {
			$array_file = $1.".txt";
			$pileup_files{$array_file}{"merged_coord"} = $PILEUP_merged_folder_abs_path."/".$array_folder_files[$contigs_merge];
	}}	

	foreach my $contigs (keys %pileup_files) {
		my ($contig_name, $coordinate_merge_file);
		if ($contigs =~ /(.*)\_ARRAY\.txt/) { 
			$contig_name = $1; 
		} else { next; }
		if ($pileup_files{$contigs}{"merged_coord"}) {
			print "Checking coordinates for sequence: ".$contig_name."\n";
			$coordinate_merge_file = $pileup_files{$contigs}{"merged_coord"};
		} else { next; }

		my @coordinates_each_contig;
		foreach my $files (sort keys %{ $pileup_files{$contigs}}) {
			if ($files =~ /.*merge.*/) {next;}
			my $string = $window_size_VARS_range;
			$string =~ s/\:\:/-/;
			my ($output_file, $error_file, $file);
			if ($variable_divergence) {
				$file = $PILEUP_merged_folder_abs_path."/".$contig_name.".id-$files"."-VD_".$variable_divergence."-CL_".$window_size_CONS."-CD_".$window_var_CONS."-VL_".$string;
			} else {
				$file = $PILEUP_merged_folder_abs_path."/".$contig_name.".id-$files"."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$window_size_CONS."-CD_".$window_var_CONS."-VL_".$string;
			}
			$output_file = $file.".out"; $error_file = $file.".err";
			
			my $array_pileup_each_taxa = $pileup_files{$contigs}{$files};
			&get_coordinates_each_taxa($array_pileup_each_taxa, $coordinate_merge_file, $files, $output_file, $error_file);
			push (@coordinates_each_contig, $output_file);
		}
		
		##########################################
		## Get Coordinates of Molecular Markers ##
		##########################################
		&get_shared_coordinates(\@coordinates_each_contig, $coordinate_merge_file, $ref_taxa);
	}
	chdir $folder; &debugger_print("Changing dir to $folder"); print "\n"; &time_log();	print "\n";

	#################################################################
	## Get the information ready for the user to visualize contigs ##
	#################################################################
	print "\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Getting the information ready to present ", "#"); DOMINO::printHeader("", "#"); 
	&get_ready_to_view_putative_markers($ref_taxa, $contigs_fasta, $ref_array_folder_files, $tmp_PILEUP_merge_folder, $ref_arrays_pileups);
	print "\n"; &time_log(); print "\n";
	chdir $folder; &debugger_print("Changing dir to $folder");

	###########################
	## Delete temporary file ##
	###########################
	print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" Deleting temporary files ", "#"); DOMINO::printHeader("", "#");
	&clean_tmp_files_marker_dir($dir); print "\n"; &time_log(); print "\n";
}

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

sub Calculate_poisson_distribution {  ## Check this level significance!
	
	my $mean = $_[0];
	my $cover_file = $_[1];
	my $h = sprintf ("%.3f", $mean);
	print "\t\t\tCoverage Mean: ".$h."\n";

	## Calculate the probability of being under a poisson distribution with the mean of our data 	
	my %total_contigs;
	open (COVERAGE, "<$cover_file");
	while (<COVERAGE>) {
		my $line = $_; chomp $line;
		my @array = split (/\s+/,$line);
		if ($array[2]) {
			my $coverage_position = $array[2];
			my $contig = $array[0];
			$total_contigs{$contig}++;
			if (defined($max_cov{$contig})) {
				if ($coverage_position > $max_cov{$contig}) {
					$max_cov{$contig} = $coverage_position;
			}} else {
				$max_cov{$contig} = $coverage_position;
	}}}
	close(COVERAGE);
	my $num_contigs = scalar keys %total_contigs;

	## If there is only one contig or a genome provided, it would be 
	## possible that it is discarded because of high coverage areas.
	unless ($num_contigs == 1) {
		foreach my $max_cov_contig (keys %max_cov) {
			my $prob_poisson;
			if ($max_cov{$max_cov_contig} > 169) { ## Factorial would be out of range!
				$prob_poisson = 0;	
			} else { $prob_poisson = &Poisson_distribution($max_cov{$max_cov_contig}, $mean); }

			if ($prob_poisson < $level_significance_coverage_distribution) { # Discard
				$discard_contigs{$max_cov_contig}++;
			} elsif ($prob_poisson eq 'nan') { # Discard
				$discard_contigs{$max_cov_contig}++;
			} elsif ($max_cov{$max_cov_contig} == 0) {  # Discard
				$discard_contigs{$max_cov_contig}++;
	}}}
	my $contigs_discarded = scalar keys %discard_contigs;
	print "\t\t\tTotal contigs: $num_contigs\n";
	print "\t\t\tContigs discarded: ".$contigs_discarded."\n";
	print "\t\t\tContigs remaining: ".($num_contigs-$contigs_discarded)."\n";
	undef %max_cov; undef %total_contigs;
	return $num_contigs;
}

sub check_file {
	
	my $file_to_check = $_[0];
	if (-e -r -s $file_to_check) { 
			if ($file_to_check =~ /(.*)\.fast.* || (.*)\.fa || (.*)\.fq/) { ## File is named fasta/fa
				my $format_returned = DOMINO::check_file_format($file_to_check);
				if ($format_returned eq "fasta") { ## fasta file is ok
					print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
				} elsif ($format_returned eq "fastq") { ##fastq file is ok
					print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
				} else { &printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
			} else { ## File is not named fasta/fa or fastq/fq
				my $format_returned = &check_file($file_to_check);
				if ($format_returned eq "fasta") { ## It is actually a fasta file
					print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
				} elsif ($format_returned eq "fastq") { ##fastq file is ok
					print "\t$file_to_check\n\t\tFile exists, is readable and non-zero character...OK\n\t\tFile format = ".$format_returned." ...OK\n";		
				} else { &printError("Wrong file provided $file_to_check..."); DOMINO::dieNicely(); }
			}
			
			## Get file name
			my @file_to_check = split("/", $file_to_check);
			my $file_name = $file_to_check[-1];
			my $flag_check = 0;
			foreach my $keys (keys %MID_species_hash) {
				if ( $file_name =~ /.*($keys).*/) {
					print "\t\tIt also contains an identifier ($1) in the name for later analysis...OK\n";
					$flag_check = 1; last;
			}}
			if ($flag_check == 0) { 
				if (!$genome_fasta) {
					&printError("$file_to_check...\n File does not contain any name provided using taxa_names option.\nMaybe it is under controlled but just bear it in mind...\nDOMINO would not die here...\n"); 
			}}
	} else { &printError("Please provide several files for each taxa..."); DOMINO::printFormat_message(); DOMINO::dieNicely(); }
}

sub check_overlapping_markers {

	my $file = $_[0];
	my (%tmp_hash, %hash_DM_markers);
	open (FILE, $file);
	print "- Checking overlapping DOMINO markers and merging into unique regions...\n";
	while (<FILE>) {
		my $line = $_;
		chomp $line;
		if ($line =~ /.*Region.*/) { next;}
		if ($line =~ /.*Cons1_ini.*/) { next;}
		
		$line =~ s/ /\t/;
		my @array = split ("\t", $line);		
		my $contig_id;
		my $tmp = $array[0];
		
		my $taxa;
		if ($option eq "msa_alignment") {
			my @tmp = split("-", $tmp);
			$taxa = $tmp[1]; $contig_id = $tmp[0];
		} else {
			$taxa = $array[4]; $contig_id = $array[0];
		}

		my @cons1_coord_array = split(":", $array[1]);
		my @var_coord_array = split(":", $array[2]);
		my @cons2_coord_array = split(":", $array[3]);
		my $string2push_asKey = $cons1_coord_array[0];			
		my $string2push_asValue = $cons1_coord_array[0].":".$cons1_coord_array[1]."_".$var_coord_array[0].":".$var_coord_array[1]."_".$cons2_coord_array[0].":".$cons2_coord_array[1];			
		push ( @{$tmp_hash{$contig_id}{$string2push_asKey}}, $string2push_asValue, $taxa);		
	}
	close(FILE);
	foreach my $keys_contigs (keys %tmp_hash) {	
		my %coord_seen;
		my $marker_counter = 0;
		foreach my $keys_markers (sort keys %{$tmp_hash{$keys_contigs}}) {
			if ($coord_seen{$keys_markers}) { next; }
			$marker_counter++;
			#print "Markers Coordinates:".$keys_markers."\n";
			my @array_coordinates;
			push (@array_coordinates, $keys_markers);
			my $bool = 1;
			my ($counter, $bad_counter) = 0;
			while ($bool) {
				$counter++;
				my $new_coord = $keys_markers + $counter;
				if ($tmp_hash{$keys_contigs}{$new_coord}) {
					push (@array_coordinates, $new_coord);
					$coord_seen{$new_coord}++;
				} else {
					$bad_counter++;
					if ($bad_counter == 3) { ## We would consider the same marker if overlapping 3pb
						($bool,$counter,$bad_counter) = 0;
			}}}
			
			my $coordinate = $tmp_hash{$keys_contigs}{$array_coordinates[-1]}[0];
			my $id = $keys_contigs."_#".$marker_counter;
			my @array = split ("_", $coordinate);		
			my @cons1_coord_array = split(":", $array[0]);
			my @var_coord_array = split(":", $array[1]);
			my @cons2_coord_array = split(":", $array[2]);
			
			$hash_DM_markers{$id}[0] = $array_coordinates[0]; # cons1_coord_array_start
			$hash_DM_markers{$id}[1] = $cons1_coord_array[1]; # cons1_coord_array_end
			$hash_DM_markers{$id}[2] = $var_coord_array[0]; # var_coord_array_start 
			$hash_DM_markers{$id}[3] = $var_coord_array[1]; # $var_coord_array_end 
			$hash_DM_markers{$id}[4] = $cons2_coord_array[0]; # cons2_coord_array_start
			$hash_DM_markers{$id}[5] = $cons2_coord_array[1]; # cons2_coord_array_end
			$hash_DM_markers{$id}[6] = $tmp_hash{$keys_contigs}{$array_coordinates[-1]}[1]; # taxa 		
		}
	}
	undef %tmp_hash; return \%hash_DM_markers;
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
		if (!$MID_species_hash{$reference}) {next;}
		foreach my $keys (keys %$hash_ref) {
			if (!$MID_species_hash{$keys}) {next;}
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
				}
			}
		}
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
		if (!$MID_species_hash{$titleline}) {next;}
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

sub discard_reads_sam {

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
	open (SAM, "<$file_sam");
	my @temp = split ("\.sam", $file_sam);
	my $output_sam = $temp[0]."_clean.sam";
	open (SAM_OUT, ">$output_sam");	
	my $discard_reads = 0;
	my $good_reads = 0;
	my $total_reads = 0;
	
	while (<SAM>) {
		chomp;
		my $line = $_;
		my $result = 1;
		if ($line =~ /^\@PG.*/ ) { next; }		
		if ($line =~ /^\@/ ) { print SAM_OUT $line."\n"; next; }		
		my @sam = split ("\\s", $line);		
		$total_reads++;
		
		if ($line =~ /SA:Z:/) {   ## Discard Multimapping reads
			$discard_reads++;   
			$result = 0;
			next;
		} elsif ($sam[2] eq '*') {  
			## Discard unmapped reads
			## When calling Bowtie2 we use -unal not to show unaligned reads
			## but still we would check it, also for user input BAM files			
			$discard_reads++;
			$result = 0;
			next;		
		} else { 			
			## Check if read is efficiently mapping
			## Even Bowtie2 is mapping end to end and allowing no clipping
			## we would check if any hard or soft clipping is introduced
			## (inherited from a previous version using BWA)
			
			my $cigar=$sam[5];
			my $cigar_record = $cigar;
			my ($xcent, $xcentmax, $pos, $NOpos);
			$xcent = $xcentmax = $pos = $NOpos = 0;
   			my $cigar_pct_convert = $cigar_pct/100;

			while ($cigar !~ /^$/){         		
         		if ($cigar =~ /^([0-9]+[MIDSH])/){
           			my $cigar_part = $1;
               		if ($cigar_part =~ /(\d+)M/){ $pos += $1;
            		}elsif($cigar_part =~ /(\d+)I/){ $NOpos += $1;
           			}elsif($cigar_part =~ /(\d+)D/){ $NOpos += $1;
           			}elsif ($cigar_part =~ /(\d+)S/){ $NOpos+= $1;
            		}elsif ($cigar_part =~ /(\d+)H/){ $NOpos+= $1;
            		}else{ # die "Unexpected cigar: $cigar\n";
         			}# close if
         			$cigar =~ s/$cigar_part//;
   			}}

			my $total = $pos + $NOpos; 				 ## total positions
			$xcentmax = $total * $cigar_pct_convert; ## Max bad CIGAR
    		$xcent = $NOpos; ## Bad CIGAR

    		if($xcent <= $xcentmax){
				print SAM_OUT $line."\n";
				$mapping_contigs{$sam[2]} = 1;
				$good_reads++;
			} else {
				if ($bowtie_local) {
					print SAM_OUT $line."\n"; next;
				} $discard_reads++;
	}}}
	close(SAM); close(SAM_OUT);
	
	my $percentage_discard = ($discard_reads/$total_reads)*100; 
	my $h = sprintf ("%.3f", $percentage_discard);
	print "This SAM file contains: $total_reads reads mapping the contigs provided\n";
	print "Because of multimapping reads, unmapped reads or low quality mapping reads some of them have been discarded: ".$h." %\n";
	print "There are ".$good_reads." reads remaining after the cleaning step\n";
	return $output_sam;
}

sub discard_coverage_sam {

	##########################################################################################
	##	 																					##
	##  This function checks the SAM files generated and discards reads with a high 		##
	##	heterogenity of coverage															##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	##	Cristina Frias Lopez 				  cristinafriaslopez@ub.edu 					##
	## 		        																		##
	##########################################################################################

	my $file_sam = $_[0];
	
	my @temp = split ("\.sam", $file_sam);
	my $output_sam = $temp[0]."_filtered.sam";
	open (SAM_OUT, ">$output_sam"); open (SAM, "<$file_sam");
	while (<SAM>) {
		chomp; my $line = $_;
		if ($line =~ /^@.*SN:(.*)\s+/ ) {
			print SAM_OUT $line."\n";	
		} else {
			my @array = split (/\s+/,$line);
			if (!$discard_contigs{$array[2]}) {
				print SAM_OUT $line."\n";	
	}}}
	close(SAM_OUT); close(SAM);
	return $output_sam;
	undef %discard_contigs;	
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

sub filter_coverage_Adjust {
	
	my $cov_file = $_[0];
	my $sorted_bam = $_[1];
		
	print "\n"; DOMINO::printHeader(" Filtering Coverage Stats for $cov_file ", "%"); print "- Reading the coverage file generated...\n\n";
	my ($mean_coverage) = &get_coverage_stats($cov_file);
	
	DOMINO::printHeader(" Filtering Statistics ", "="); 
	my $contigs = &Calculate_poisson_distribution($mean_coverage, $cov_file); DOMINO::printHeader("", "="); print "\n"; 
	my @tmp_sam = split("\.sorted.bam", $sorted_bam);
	my $sam_filter;
	if ($contigs == 1) {
		$sam_filter = $tmp_sam[0]."_filtered.sam";
		File::Copy::move($tmp_sam[0]."\.sam", $sam_filter);		
	} else {
		DOMINO::printHeader(" Adjusting the BAM file filtered ", "%"); print "- Reading the coverage file generated...\n\n";
		my $file_to_Adjust = $tmp_sam[0]."\.sam";
		$sam_filter = &discard_coverage_sam($file_to_Adjust);
	}
	return $sam_filter;
}

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
	my @files = @$files_dir_ref;
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		if ($files[$i] =~ /.*id-(.*)(\_R\d+)\.fastq/g) {
			my $MID_found = $1;
			if ($MID_species_hash{$MID_found}) { ## Make sure of only taking the taxa user specified
				push (@clean_fastq_file_abs_path, $clean_folder."/".$files[$i]); ## push the whole file path			
		}} elsif ($files[$i] =~ /.*id-(.*)\.fastq/g) {
			my $MID_found = $1;
			if ($MID_species_hash{$MID_found}) { ## Make sure of only taking the taxa user specified
				push (@clean_fastq_file_abs_path, $clean_folder."/".$files[$i]); ## push the whole file path			
}}}}

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
	my $hash_ref_coordinates = DOMINO::readFASTA_hash($file_MID_array);
	my %hash_seq = %{$hash_ref_coordinates};

	##### open output
	open (OUT, ">$output_file") or &Error($!);
	open (ERR, ">$error_file") or &Error($!);

	### open coord and check coord of file 1
	open (COORD,"<$merge_file_coord");
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

		### initialite var
		my ($count_mono_ref_cons1, $count_poli_ref_cons1) = 0;
		my ($count_mono_ref_VAR, $count_poli_ref_VAR) = 0;
		my ($count_mono_ref_cons2, $count_poli_ref_cons2) = 0;
		my $size_VAR_MID = 0;

		### split data of regions
		my ($cons1_START, $cons1_END) = split(/:/,$cons1);
		my ($VAR_START, $VAR_END) = split(/:/,$var);
		my ($cons2_START, $cons2_END) = split(/:/,$cons2);

		## Check if exists
		if (exists ($hash_seq{$contig_MID}) ){
			$seq_MID = $hash_seq{$contig_MID};
		} else {
			print ERR "$contig_MID is not mapping for file $file_MID_array\n";
			next;
		}

		###################################################################
		###			 Get the sequence of the different regions 			### 
		###################################################################
	
		### CONS1 ###
		my $size_CONS1 = $cons1_END - $cons1_START + 1;
		my $pileup_cons1 = substr ($seq_MID, $cons1_START-1, $size_CONS1);

		### VAR ###
		my $size_VAR = $VAR_END - $VAR_START +1; 
		my $pileup_var = substr ($seq_MID, $VAR_START-1, $size_VAR);	

		### CONS2 ###
		my $size_CONS2 = $cons2_END - $cons2_START + 1 ;
		my $pileup_cons2 = substr($seq_MID, $cons2_START-1, $size_CONS2);

		###################################################################
		### 		Check  the coverage of each region					###
		###################################################################
	
		### CONS1 REGION ###
		my $size_CONS_degenerate = $size_CONS1 - $window_var_CONS;
		$count_mono_ref_cons1 = $pileup_cons1 =~ tr/0//;
		$count_poli_ref_cons1 = $pileup_cons1 =~ tr/1//;

		### VARIABLE REGION ###
		$count_mono_ref_VAR = $pileup_var =~ tr/0//;
		$count_poli_ref_VAR = $pileup_var =~ tr/1//;
		
		my $total = $count_poli_ref_VAR + $count_mono_ref_VAR;

		### CONS2 REGION ###
		$count_mono_ref_cons2 = $pileup_cons2 =~ tr/0//; 
		$count_poli_ref_cons2 = $pileup_cons2 =~ tr/1//;
	
		###################################################################
		### 		Check the composition of the regions				###
		###################################################################
	
		my $delimiter = "\t";  # empty string
		my $string = join($delimiter, @regions);

		## Variability on variable region: User provides a minimum for variation within this region
		#my $expected_var_sites = ($variable_divergence/100) * $size_VAR; 0-100%
		#my $expected_var_sites = $variable_divergence * $size_VAR; # 0-1 
		
		###################################################################
		### 		Print coordinates if the meet requirements			###
		###################################################################
		my $flag_error;
		if(($count_mono_ref_cons1 >= $size_CONS_degenerate) && ($count_mono_ref_cons2 >= $size_CONS_degenerate)) {
			if ($variable_divergence) {
				my $expected_var_sites = $variable_divergence * $total;
				if ($count_poli_ref_VAR > $expected_var_sites) {
					print OUT "$string\n"; 
				} else { $flag_error++; }
			} else {
				if ($count_poli_ref_VAR > $variable_positions_user_min) {
					if ($count_poli_ref_VAR < $variable_positions_user_max) {
						print OUT "$string\n";
					} else { $flag_error++;}
				} else { $flag_error++; }
			}
		} else { $flag_error++; }
		if ($flag_error) {print ERR $string."\n";}
	}# close while
	close(COORD); close(OUT); close(ERR);
}

sub get_coverage_stats {
	
	my $cover_file = $_[0];
	my ($mean_coverage, $sum_coverage, $total_positions);
	
	open (COVERAGE, "<$cover_file");
	while (<COVERAGE>) {
		my $line = $_;
		chomp $line;
		my @array = split (/\s+/,$line);
		my $fake_pb = $array[1] - 1;
		if ($array[2]) {
			$sum_coverage += $array[2];
			$total_positions++;
	}}
	close(COVERAGE);
	$mean_coverage = $sum_coverage/$total_positions;
	
	# Return:
	return ($mean_coverage);
}

sub get_MID_contigs {
	
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
		if ($files[$i] =~ /.*id\-(.*)\.contigs\.fasta/) {
			my $MID_found = $1;
			if ($MID_found =~ /(.*)\_R\d+/) {
				my $MID_found_tmp = $1; $MID_found = $MID_found_tmp;				
			}			
			if ($MID_species_hash{$MID_found}) { ## Make sure of only taking the taxa user specified
				push (@contigs_fasta_file_abs_path, $tmp_file_abs_path); 
}}}}

sub get_filtered_sam_files {

	##########################################################################################
	##	 																					##
	##  This function gets the filtered SAM files generated in the previous cleaning step	##
	## 		        																		##
	##	Jose Fco. Sanchez Herrero, 08/05/2014 jfsanchezherrero@ub.edu						##
	## 		        																		##
	##########################################################################################

	my $MID_id = $_[0];
	my $folder = $_[1];
	my ($ref_bam_file, @pileup_Arrays_sub, @clean_filtered_sam_files_sub);
	my $files_ref = DOMINO::readDir($align_dirname);
	my @files = @$files_ref;
		
	for (my $i = 0; $i < scalar @files; $i++) {
		if ($files[$i] eq ".DS_Store" || $files[$i] eq "." || $files[$i] eq ".." ) { next; }
		if ($files[$i] =~ /.*fastq/) { next; }		
		my $abs_path_file = $align_dirname."/".$files[$i];
		if (-d $abs_path_file) {
			if ($files[$i] eq $MID_id) {
				my $dir_files_ref = DOMINO::readDir($abs_path_file);
				my @dir_files = @$dir_files_ref;
				my $scalar = scalar @dir_files;
				for (my $j = 0; $j < scalar @dir_files; $j++) {
					if ($dir_files[$j] =~ /taxa\_(.*)\_clean_filtered\.sam/) {
						if ($MID_species_hash{$1}) {
							push (@clean_filtered_sam_files_sub, $dir_files[$j]);
							system("ln -s $align_dirname/$files[$i]/$dir_files[$j] $folder");
						}
					} elsif ($dir_files[$j] =~ /ARRAY\_files\_taxa\_(.*)\_clean\_filtered\.profile/) {
						if ($MID_species_hash{$1}) {
							my $tmp = $abs_path_file."/".$dir_files[$j];
							push (@pileup_Arrays_sub, $tmp);
							system("ln -s $tmp $folder");
						}
					} elsif ($dir_files[$j] =~ /.*$MID_id.*sorted.bam/) {
						$ref_bam_file = $dir_files[$j];
						system("ln -s $align_dirname/$files[$i]/$dir_files[$j] $folder");						
					} else { 
						next; 
	}}}}}
	return ($ref_bam_file, \@pileup_Arrays_sub, \@clean_filtered_sam_files_sub);
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
	
	my $pileup_Arrays_sub2_ref = $_[4];
	my @pileup_Arrays_sub2 = @$pileup_Arrays_sub2_ref;
	my (%fasta_msa, $id, $seq_name, $old_seq_name, $marker_number);

	## Print in tmp file for sorting and obtaining unique
	my $tmp_file = "tmp_file.txt";
	open(TMP, ">$tmp_file");
	for my $scaffold (keys %coord_contig ) {    		
		## Check we find markers for all the taxa desired
		my $string = join(",", @{$coord_contig{$scaffold}});
		$string =~ s/merged\,/ /g;
		my @array_taxa_split = split(",", $string);
		my $species2find = $minimum_number_taxa_covered;
  		if (scalar @array_taxa_split < $species2find) { next; }   		
    	## Write DOMINO Markers Coordinates in tmp txt file
	   	my @contig_Array = split (/\s+/, $scaffold);
    	print TMP "$scaffold $string\n"; #Contig1161 322:341 342:741 742:761 taxa1,taxa2,taxa3
   		$putative_markers{$contig_Array[0]}++;
	} close(TMP);	

	#print Dumper \%putative_markers;
	my $hash_markers_ref_collapse = &check_overlapping_markers($tmp_file);

	### Get sequences for each taxa
	for (my $i = 0; $i < scalar @pileup_Arrays_sub2; $i++) {
		my $array_files_ref = DOMINO::readDir($pileup_Arrays_sub2[$i]);
		my @array_files = @$array_files_ref;
		my $name;
		if ($pileup_Arrays_sub2[$i] =~ /.*\_taxa\_(.*)\_clean.*/) { $name = $1; }
		for (my $y = 0; $y < scalar @array_files; $y++) {
			if ($array_files[$y] eq "." || $array_files[$y] eq ".." || $array_files[$y] eq ".DS_Store") {next;}
			if ($array_files[$y] =~ /.*ARRAY\.txt/) { }
			if ($array_files[$y] =~ /(.*)_sequence\.fasta/) { 
				my $tmp = "$pileup_Arrays_sub2[$i]/$array_files[$y]";
				$fasta_msa{$1}{$name} = $tmp;		
	}}}

	## Get the reference sequence
	my %reference_sequence;
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
		if ($putative_markers{$match}) { $fasta_msa{$match}{$reference} = uc($sequence); }
	} close(FILE); $/ = "\n";
	
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
	
	###################################
	### Read the merge PILEUP file  ###
	###################################
	open (MERGE_COORD,"<$file_coord") or die "Cannot open file $file_coord";
	while(<MERGE_COORD>){
		chomp;
		my $line = $_;
		next if ($line =~ m/^Contig\tCons1_ini:.*/o);
    	next if $line=~ m/^\s*$/o;
    	next if $line=~ /^\#/o;
		push (@{ $coord_contig{$line} }, "merged"); # create empty array for each contig
		push (@{ $coord_contig{$line} }, $reference); 
	}
	close (MERGE_COORD);	
	if (!%coord_contig) { return; } ## No coordinates were found
	
	###################################
	### Check each taxa coordinates ###
	###################################	
	my @array_file_coordinates = @$file_coord_array_ref;
	for (my $i = 0; $i < scalar @array_file_coordinates; $i++) {
		my $name;
		if ($array_file_coordinates[$i] =~ /.*id\-(.*)\-V(D|P).*/) { $name = $1; }
		open (FILE, "<$array_file_coordinates[$i]") || die "Cannot open file $array_file_coordinates[$i]";
		while(<FILE>){
			chomp;
			my $line= $_;
    	    next if ($line =~ m/^Contig\tCons1_ini:.*/o); next if $line=~ m/^\s*$/o; next if $line=~ /^\#/o;
	        if (exists $coord_contig{$line}){ 
	        	push (@{ $coord_contig{$line} }, $name); 
	    }}
		close (FILE);
}}

sub generate_bam {
	my $sam_file = $_[0];
	my $avoid = $_[1];
	my @temp = split ("\.sam", $sam_file);
	my $name = $temp[0]; my $bam_file = $name.".bam";
	print "- Generating a BAM file for $sam_file\n"; 	
	my $system_samtools_sam2bam = $samtools_path." view -@ $num_proc_user -bS -o $bam_file $sam_file";
	&debugger_print("SAMTOOLS command: $system_samtools_sam2bam\n");	
	my $system_call = system ($system_samtools_sam2bam);
	if ($system_call != 0) {
		if (!$avoid) { &printError("Some error happened when calling SAMTOOLs for SAM -> BAM conversion $sam_file to $bam_file...."); DOMINO::dieNicely(); }
	}
	my $sorted = &generate_sorted_bam($bam_file);	
	return $sorted;
}

sub generate_sorted_bam {
	my $bam_file = $_[0];
	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	print "- Sorting the BAM file: $bam_file\n"; 	
	my $sorted_bam = $name.".sorted.bam";	
	my $system_samtools_sort = $samtools_path." sort -@ $num_proc_user -o $sorted_bam $bam_file";
	&debugger_print("SAMTOOLS command: ".$system_samtools_sort."\n");
	my $system_call_2 = system ($system_samtools_sort);
	if ($system_call_2 != 0) {
		&printError("Some error happened when calling SAMTOOLs for sorting BAM file...."); DOMINO::dieNicely();
	}
	return $sorted_bam;
}

sub generate_sam {
	my $bam_file = $_[0];
	my @temp = split ("\.bam", $bam_file);
	my $name = $temp[0];
	print "- Generating a SAM file for $bam_file\n"; 	
	my $system_samtools_bam2sam = $samtools_path." view -@ $num_proc_user -h ".$bam_file." -o ".$name.".sam";
	&debugger_print("SAMTOOLS command: $system_samtools_bam2sam\n");
	my $system_call = system ($system_samtools_bam2sam);
	if ($system_call != 0) {
		&printError("Some error happened when calling SAMTOOLs for BAM -> SAM conversion....");  DOMINO::dieNicely();
	}
	return $name.".sam";
}

sub generate_index_bam {
	my $bam_file = $_[0];
	print "- Generating an index bam file for $bam_file\n"; 	
	my $index_system = $samtools_path." index ".$bam_file;
	&debugger_print("SAMTOOLS command: ".$index_system."\n");
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
	
	my $sorted_bam = $_[0];
	my $contig_file = $_[1];
	my $reference_hash_fasta = $_[2];
	my $reference_id = $_[3];
	
	my $dir_path = $align_dirname."/".$reference_id;
	
	my @temp_name = split ("\.sorted.bam", $sorted_bam);
	my ($ID, @sam);
	my $input_pileup = $temp_name[0].".profile";
	my $input_pileup_path = $dir_path."/".$input_pileup;

    DOMINO::printHeader(" Generate a PILEUP for $sorted_bam ", "%");
	my $pileup_command = $samtools_path." mpileup -f ".$contig_file." -o ".$input_pileup." ".$sorted_bam;
	&debugger_print("SAMTOOLS command: ".$pileup_command."\n");
	my $sytem_command_pileup = system ($pileup_command);
	if ($sytem_command_pileup != 0) {
		&printError("Exiting the script. Some error happened when calling SAMtools for generating the PILEUP for the file $contig_file...\n"); DOMINO::dieNicely();
	}
	print "Done...\n\n";
	
	my $tmp = $dir_path."/ARRAY_files_$input_pileup"; mkdir $tmp, 0755; chdir $tmp;
	&debugger_print("Changing dir to $tmp");
    DOMINO::printHeader(" Filtering the PILEUP generated ", "%");
	print "- Splitting the PILEUP file $input_pileup...\n";
	my ($previous_contig, $previous_fasta_contig, @array_positions, @fasta_positions);
	open (PILEUP,"<$input_pileup_path");
	while (<PILEUP>){
		my $line = $_;
		chomp $line;
		next if $line=~ m/^\s*$/o;
		next if $line=~ /^\#/o;
		next if $line=~ m/\[REDUCE RESULT\]/;
		$line =~ s/\s/\t/g;
		my @pileup = split(/\t/,$line); ##	HKU61D301.clean.MID4_Nem_098_c8679	161	t	3	,,.	FA=
		my $contig = $pileup[0];
		my ($num_pos_base, $pos_base);
		$num_pos_base = $pileup[1];
		my $num_pos_array = $num_pos_base -1;
		if (!$previous_contig) {
			my $array_positions_ref;
			($previous_contig, $array_positions_ref) = &initilize_contig($contig, $reference_hash_fasta);
			@array_positions = @$array_positions_ref;
			@fasta_positions = @array_positions;			
		} else {
			if ($previous_contig ne $contig) {
				&print_coordinates(\@array_positions, \$previous_contig, $reference_id); ## Print array into file $previous_contig
				&print_fasta_coordinates(\@fasta_positions, \$previous_contig, $reference_id); ## Print array into file $previous_contig
				my $array_positions_ref;
				($previous_contig, $array_positions_ref) = &initilize_contig($contig, $reference_hash_fasta);
				@array_positions = @$array_positions_ref;
				@fasta_positions = @array_positions;
		}}
		
		## Get array for each position
		if ($pileup[3] != 0) {
			my $read_base = $pileup[4];
			my $ref_base = $pileup[2];
			my @base_record = split(//, $read_base);
			my (%posibilities, %polymorphism);
			my @base_parse;
			my $base_counter=0;
			for (my $i = 0; $i < scalar @base_record; $i++) {  
				if ($base_record[$i] =~ m/\^/) { ## Starting to map a new read
					$i++; next; 
				}
				$base_record[$i]= uc($base_record[$i]);
				if ($base_record[$i] =~ m/A|G|C|T|a|g|c|t/) {
					$polymorphism{$base_record[$i]}++;
					$base_counter++;
					push (@base_parse, $base_record[$i]);
				} elsif ($base_record[$i] eq "." || $base_record[$i] eq ",") {
					$polymorphism{$ref_base}++;
					$base_counter++;
					push (@base_parse, $ref_base);
				} elsif ($base_record[$i] =~ m/(\+|\-)/) {	
					## There is an INDEL	
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
				} else {
					next;
			}}
			
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
				next;	
			}
			
			if ($base_counter == 1) { ## There is only base mapping...
				unless ($map_contig_files || $DOMINO_simulations) {
					$array_positions[$num_pos_array] = 'N'; ## Not informative enough
					$fasta_positions[$num_pos_array] = 'N'; ## Not informative enough
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
	}}
	close(PILEUP);
	&print_coordinates(\@array_positions, \$previous_contig, $reference_id);
	&print_fasta_coordinates(\@fasta_positions, \$previous_contig, $reference_id); ## Print array into file $previous_contig

	## Return
	#printDebug $dir_path."\n";
	chdir $dir_path;
	&debugger_print("Changing dir to $dir_path");

	
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
		my $current_previous_contig = $current_contig;
		my @array_positions;
		if ($mapping_contigs{$current_previous_contig}) {
			my $tmp_size = ${ $reference_hash_fasta}{$current_contig};
			my $size = $tmp_size;
			@array_positions = ("-") x $size;
		}
		my $ref = \@array_positions;		
		return ($current_previous_contig, $ref);
	}
	
	sub print_coordinates {
		## Print array into file $previous_contig
		my $coord_array_ref = $_[0];
		my $contig_name = $_[1];
		my $ref_id = $_[2];

		my @coord_array = @$coord_array_ref;
		my $seq_contig = join "", @coord_array;
		my $var_sites = $seq_contig =~ tr/1/1/;

		if ($var_sites != 0) {
			my $array_file = $$contig_name."_ARRAY.txt";
			open (FH, ">$array_file"); print FH ">$$contig_name\n$seq_contig\n"; close(FH);	
		}
	}
	
	sub print_fasta_coordinates {
		## Print array into file $previous_contig
		my $coord_array_ref = $_[0];
		my $contig_name = $_[1];
		my $ref_id = $_[2];
		my @coord_array = @$coord_array_ref;
		my $seq_contig = join "", @coord_array;
		my $array_file = $$contig_name."_sequence.fasta";
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
	foreach my $MID_barcodes (keys %MID_species_hash) {
		$name .= $MID_barcodes."_";
	}	
	## Merge the different bam files
	DOMINO::printHeader(" Merging the BAM files ", "%");
	my $system_samtools_merge = $samtools_path." merge -r -@ $num_proc_user -h header.sam ".$name."merged.bam ".$sorted_bam;
	&debugger_print("SAMTOOLS command: ".$system_samtools_merge."\n");
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
	my $id = $_[1];
	my $filename = $_[2];
	
	my %hash = %$ref_hash;
	my %hash2return;
	foreach my $sample (keys %hash) {
		my @array = @{ $hash{$sample} };
		my $alleles = scalar @array;
		if ($alleles == 1) {
			$hash2return{$sample} = $hash{$sample}[0];
		} elsif ($alleles == 2){
			my @allele1 = split("",$hash{$sample}[0]);
			my @allele2 = split("",$hash{$sample}[1]);
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
			DOMINO::printError_log("CLocus: $id -- ID: $sample contains more than 2 alleles");
	}}
	if (!%hash2return) {return;}
	my $file = $msa_dirname."/".$filename."_loci_".$id.".fasta";
	open (OUT, ">$file");
	foreach my $keys (keys %hash2return) {
		if ($MID_species_hash{$keys}) {
		print OUT ">".$keys."\n".$hash2return{$keys}."\n";
	}}
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
    	$worksheet_parameters->write($row, $col, $window_size_CONS, $format_right); $row++;
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
	foreach my $taxa (sort keys %MID_species_hash) {
		$col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, $taxa, $format_left); $col++;	
		$worksheet_parameters->write($row, $col, $counter, $format_right); $counter++;	$row++;
	}
	$workbook->close();
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
			print OUT_MSA ">".$aln{$keys}{$taxa}{"name"}."\n".$aln{$keys}{$taxa}{"seq"}."\n";
			unless (grep /$aln{$keys}{$taxa}{"name"}/, @array_taxa ) {
				push (@array_taxa, $aln{$keys}{$taxa}{"name"});
			}
		}
		close(OUT_MSA);	
		$region_provided = "";
	}
	return \@array_taxa;
}

sub sliding_window_conserve_variable {
	
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
	print "- Checking each contig using a sliding window approach...\n";
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
			$output_file_coordinates .= $array_split[0]."-VD_".$variable_divergence."-CL_".$window_size_CONS."-CD_".$window_var_CONS."-VL_".$string.".tab";
		} else {
			$output_file_coordinates .= $array_split[0]."-VPmin_".$variable_positions_user_min."-VPmax_".$variable_positions_user_max."-CL_".$window_size_CONS."-CD_".$window_var_CONS."-VL_".$string.".tab";
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
			my $max = $seqlen - $window_size_CONS - $window_size_VARS_max; #	for (my $j = 0; $j < $seqlen - $window_size_CONS + 1; $j++){
			if ($max < 0) {next;}
			
			for (my $j = 0; $j < $max; $j++){			
				## Set a threshold for $window_size_VARS -VL 200::400				
				for (my $h=$window_size_VARS_min; $h <= $window_size_VARS_max; $h++) {
					
					my $pos_var_ini = $j + $window_size_CONS;
					my $pos_var_end = $j + $window_size_CONS + $h - 1;
					my $pos_cons2 = $pos_var_end + 1;			#$j + $window_size_CONS + $window_size_VARS;
					
					my $cons1_string = substr($dna_seq, $j, $window_size_CONS);
					my $VAR_string = substr($dna_seq, $pos_var_ini, $h);
					my $cons2_string = substr($dna_seq, $pos_cons2, $window_size_CONS);
					
					my $count_cons1 = $cons1_string  =~ tr/0//; ## Count the variable positions in each 
					my $count_vars =  $VAR_string 	 =~ tr/1//; ## Count the variable positions in each 
					my $count_cons2 = $cons2_string  =~ tr/0//; ## Count the variable positions in each 
					
					if ($count_vars == 0) {next;} ## It is not necessary to keep on calculating anything as there is no variation at all
					
					# Conserved Region1
					$pos_cons1_ini = $j + 1; 
					$pos_cons1_end = $j + $window_size_CONS;
		
					# Variable Region1
					$pos_VAR_ini = $pos_cons1_end + 1; 
					$pos_VAR_end = $pos_VAR_ini + $h - 1;
		
					# Conserved Region2
					$pos_cons2_ini = $pos_VAR_end + 1;  
					$pos_cons2_end = $pos_cons2_ini + $window_size_CONS - 1;
					
					if (($pos_cons2_ini > $seqlen) && ($pos_cons2_end > $seqlen)) { next; }					
					$degenerate_cons = $window_size_CONS - $window_var_CONS;
		
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
		}}}}
		}; # close for each contig inside the hash
		close(OUTPUT);
	}
	if ($identify_markers) { return \@output_files; }
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


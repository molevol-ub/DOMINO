#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
	require Spreadsheet::WriteExcel;
}

#################################################################
my $absolute_path = $ARGV[0];
my $path = $ARGV[2];
my $markers_file = $ARGV[1];
if(!defined($markers_file)) { DOMINO::printError("No marker files are provided in DM_print_Excel module : [$0]\n"); exit; } ## give error message for DOMINO debug
#################################################################

#################################################################
my $domino_version ="DOMINO v1.1 ## Revised 07-11-2018";
my $hash_parameters = DOMINO::get_parameters($absolute_path."/", "markers");
my $hash_files = DOMINO::get_DOMINO_files($absolute_path."/", "markers");

my $hash_parameters2 = DOMINO::get_parameters($absolute_path."/", "assembly");
my $hash_files2 = DOMINO::get_DOMINO_files($absolute_path."/", "clean_data");

my $hash_parameters3 = DOMINO::get_parameters($absolute_path."/", "clean_data");

#print Dumper $hash_parameters;
#print Dumper $hash_parameters2;
#print Dumper $hash_files;
#print Dumper $hash_files2;

#################################################################

#################################################################
my $no_parameters; if ($hash_parameters == 0) { $no_parameters = 1; }
my @array_markers;
open (FILE, "$markers_file");
while (<FILE>) {
	chomp; if ($_ =~ /.*Vari.*/) {next;}
	push (@array_markers, $_);
} close (FILE);

#################################################################################
##	Once the coordinates are found, print different files with the information ##
#################################################################################	
### open Output and Error file
my $excel_woorkbook_name = $path."/DM_markers-summary.xls";
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
$worksheet_parameters->write($row, $col, $path, $format); $row++;
$col = $first_col;

$worksheet_parameters->write($row, $col, "OPTION:", $format); $col++;

if ($no_parameters) {
	if ($$hash_parameters{"markers"}{"RADseq"}[0]) {
		$worksheet_parameters->write($row, $col, "RADseq", $format); $row++;
	} else { $worksheet_parameters->write($row, $col, $$hash_parameters{"markers"}{"option"}, $format); $row++; }
} else {
	if ($$hash_parameters3{'clean_data'}{'option'}) {
		$worksheet_parameters->write($row, $col, $$hash_parameters2{'clean_data'}{'option'}, $format); $row++;
}} $row++; $row++;

### Writing markers coordinates
my $worksheet_markers = $workbook->add_worksheet("Output - Markers");
my $col_markers = $first_col = 0; my $row_markers = $second_col = 1;
$worksheet_markers->write($row_markers, $col_markers, "OUTPUT", $format_main_heading); $row_markers++; $col_markers++; $row_markers++;

if ($$hash_parameters{'marker'}{'behaviour'} eq 'selection') {
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
	}
} else {
	$worksheet_markers->write($row_markers, $col_markers, "Conserved region: Left", $format_bold); $col_markers++; 
	$worksheet_markers->write($row_markers, $col_markers, "Variable region", $format_bold); $col_markers++; 
	$worksheet_markers->write($row_markers, $col_markers, "Conserved region: Right", $format_bold); $col_markers++;
	$worksheet_markers->write($row_markers, $col_markers, "Taxa included", $format_bold); $col_markers++;
	$worksheet_markers->write($row_markers, $col_markers, "Marker length", $format_bold);  $col_markers++;
	$worksheet_markers->write($row_markers, $col_markers, "Variation", $format_bold);  $row_markers++;
	$col_markers = $first_col;
	if ($$hash_parameters{"marker"}{"option"} eq "msa_alignment") { 
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
}} 

unless ($no_parameters) {
	if (exists $$hash_parameters3{'clean_data'}{'start'}) {
		$worksheet_parameters->write($row, $first_col, "PRE-PROCESSING PHASE:", $format_main_heading); 
		$col = $first_col; $row++; $row++;
		$worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'start'}, $format_left); $row++; $row++;
		$col = $first_col; 
		$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'folder'}, $format); $row++; $row++;
		$col = $first_col;
		$worksheet_parameters->write($row, $col, "INPUT DATA:", $format); $col++;
		$worksheet_parameters->write($row, $col, "Original Files included", $format_bold); $col++;
		$worksheet_parameters->write($row, $col, "DOMINO Label", $format_bold); $col++;
		foreach my $keys (sort keys %{$$hash_parameters3{'clean_data'}{'clean_files'}}) {
			$col = $first_col; $col++; $row++;
			$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'clean_files'}{$keys}, $format_left);
			$col++; $worksheet_parameters->write($row, $col, $keys, $format_left);
		}
		$row++; $row++; $row++;
		$col = $first_col;
		$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col++; 
		$worksheet_parameters->write($row, $col, "Mismatches allowed in the barcode sequence:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'mismatches'}, $format_right); 
		$col = $first_col + 1; $row++;
		$worksheet_parameters->write($row, $col, "Minimum read length:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'read_length'}, $format_right); 
		$col = $first_col + 1; $row++;
		$worksheet_parameters->write($row, $col, "Minimum QUAL:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'QUAL'}, $format_right); 
		$col = $first_col + 1; $row++;
		$worksheet_parameters->write($row, $col, "Minimum length of a read satisfying QUAL cutoff:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'QUAL_cutoff'}, $format_right); 
		$col = $first_col + 1; $row++;
		$worksheet_parameters->write($row, $col, "Threshold for complexity/entropy of a read:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters3{'clean_data'}{'entropy'}, $format_right); 
		$col = $first_col + 1; $row++; $row++;
		$worksheet_parameters->write($row, $col, "Database(s)", $format_bold); $row++;
		my @array_db = @{ $$hash_parameters3{'clean_data'}{'db'}};
		for (my $i=0; $i < scalar @array_db; $i++) {
			$worksheet_parameters->write($row, $col, $array_db[$i], $format_left); $row++;
		} $row++; $row++;
	
		### Fix this
		if ($$hash_files2{'clean_data'}{"QC_analysis"}) {
			my %QC_hash = %{$$hash_parameters2{'clean_data'}{'QC_analysis'}};
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
		}}
		## up to here
	
	}

	if (exists $$hash_parameters2{'assembly'}{"date"}) {
		$worksheet_parameters->write($row, $first_col, "ASSEMBLY PHASE:", $format_main_heading); 
		$row++; $row++;	$col = $first_col; 
		$worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters2{'assembly'}{'start'}, $format_left); 
		$row++; $row++; $col = $first_col; 
		$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters2{'assembly'}{'folder'}, $format_left); 
		$row++; $row++; $col = $first_col;
		
		#if ($$hash_parameters{'assembly'}{'files'} eq "Default DOMINO cleaning files") {
		$worksheet_parameters->write($row, $col, "INPUT DATA", $format); $col++;
		$worksheet_parameters->write($row, $col, "Default DOMINO cleaned files", $format_left); $row++; $col++;
		## Control if user provides different files
		
		$row++; $row++; $col = $first_col;
		$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, "Minimum read score (MIRA):", $format); $col++;
		$worksheet_parameters->write($row, $col, $$hash_parameters2{'assembly'}{'mrs'}, $format_left); $col++;
		
		if ($$hash_parameters{'assembly'}{'CAP3'}) {
			$col = $second_col; $row++;
			$worksheet_parameters->write($row, $col, "CAP3:", $format); $col++;
			$worksheet_parameters->write($row, $col, "Enabled", $format_left); $col++;			
			$col = $second_col; $row++;
			$worksheet_parameters->write($row, $col, "Similarity:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters2{'assembly'}{'simCAP3'}, $format_left); $col++;
			$col = $second_col; $row++;
			$worksheet_parameters->write($row, $col, "Overlapping:", $format); $col++;
			$worksheet_parameters->write($row, $col, $$hash_parameters2{'assembly'}{'overCAP3'}, $format_left); $col++;
		} else {
			$col = $second_col; $row++;
			$worksheet_parameters->write($row, $col, "CAP3: Disabled", $format);
		} $row++; $row++; $row++; 
		
		## Fix
		if ($$hash_files2{'assembly'}{'stats'}) {
			my %hash = %{$$hash_parameters2{'assembly'}{'stats'}};
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
					my @array = split(",", $line);
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
		}}
		## Fix up to here
}}

$worksheet_parameters->write($row, $first_col, "ALIGNMENT/MAPPING PHASE:", $format_main_heading);
$col = $first_col; $row++; $row++;
$worksheet_parameters->write($row, $col, "DATE:", $format); $col++;
$worksheet_parameters->write($row, $col, $$hash_parameters{'mapping'}{'date'}[0], $format_left); 
$col = $first_col; $row++; $row++;
$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
$worksheet_parameters->write($row, $col, $$hash_parameters{'mapping'}{'folder'}[0], $format);
$col = $first_col; $row++; $row++;
$worksheet_parameters->write($row, $col, "DATA:", $format); $col++;
$worksheet_parameters->write($row, $col, "Option:", $format); 
my @option_position = ($row, ($col+1)); $row++;

if ($$hash_parameters{"marker"}{"option"} eq "msa_alignment") { 
	my @row_typefile = ($row, $col+1);
	$worksheet_parameters->write($row, $col, "Type of file:", $format); $row++;
	if ($$hash_parameters{"marker"}{"RADseq"}) {
		$worksheet_parameters->write($option_position[0], $option_position[1], "RADseq", $format_left); 
		if ($$hash_parameters{'marker'}{'pyRAD'}) {
			$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "pyRAD file provided", $format_left); 
		} elsif ($$hash_parameters{'marker'}{'STACKS'}) {
			$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "STACKS file provided", $format_left); 
	}} elsif ($$hash_parameters{'marker'}{'MSA'}) {
		$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple Sequence Alignment", $format_left); 
		$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "Alignment file provided", $format_left); 
	} elsif ($$hash_parameters{'marker'}{'MSA_folder'}) {
		$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple Sequence Alignment", $format_left); 
		$worksheet_parameters->write($row_typefile[0], $row_typefile[1], "Alignment folder provided", $format_left); 
}} else {
	if ($$hash_parameters{"marker"}{"option"} eq 'user_assembly_contigs') {
		$worksheet_parameters->write($option_position[0], $option_position[1], "Multiple references provided", $format_left); 
	} elsif ($$hash_parameters{"marker"}{"option"} eq 'DOMINO_files') {
		$worksheet_parameters->write($option_position[0], $option_position[1], "DOMINO files ", $format_left); 
	} elsif ($$hash_parameters{"marker"}{"option"} eq 'genome') {
		$worksheet_parameters->write($option_position[0], $option_position[1], "Single reference provided", $format_left); 
	} 
	$row++; $col = $first_col; 
	$worksheet_parameters->write($row, $col, "PARAMETERS:", $format); $col++; 
	my $tmp_col = $col;
	$worksheet_parameters->write($row, $tmp_col, "Read Gap Open penalty ", $format);
	$worksheet_parameters->write($row, ($tmp_col + 1), $$hash_parameters{"mapping"}{"rdgopen"}[0], $format_right); $row++;
	$worksheet_parameters->write($row, $tmp_col, "Read Gap Extension penalty ", $format);
	$worksheet_parameters->write($row, ($tmp_col + 1), $$hash_parameters{"mapping"}{"rdgexten"}[0], $format_right); $row++;
	$worksheet_parameters->write($row, $tmp_col, "Reference Gap Open penalty ", $format);
	$worksheet_parameters->write($row, ($tmp_col + 1), $$hash_parameters{"mapping"}{"rfgopen"}[0], $format_right); $row++;
	$worksheet_parameters->write($row, $tmp_col, "Reference Gap Extension penalty ", $format);
	$worksheet_parameters->write($row, ($tmp_col + 1), $$hash_parameters{"mapping"}{"rfgexten"}[0], $format_right); $row++;
	$worksheet_parameters->write($row, $tmp_col, "Mismath penalty ", $format);
	$worksheet_parameters->write($row, ($tmp_col + 1), $$hash_parameters{"mapping"}{"mis_penalty"}[0], $format_right); $row++;
}

$row++; $row++;
$worksheet_parameters->write($row, $first_col, "MARKER DISCOVERY PHASE:", $format_main_heading);  
$col = $first_col; 
$worksheet_parameters->write($row, $col, "STARTED:", $format); $col++;
$worksheet_parameters->write($row, $col, $$hash_parameters{'marker'}{'date'}[0], $format_left); $row++; $row++;
$col = $first_col;
$worksheet_parameters->write($row, $col, "FOLDER:", $format); $col++;
$worksheet_parameters->write($row, $col, $$hash_parameters{'marker'}{'folder'}[0], $format); $row++; $row++;
$col = $first_col;
$worksheet_parameters->write($row, $col, "DEVELOPMENT MODULE:", $format); $col++;
if ($$hash_parameters{'marker'}{'behaviour'}[0] eq 'selection') {
	$worksheet_parameters->write($row, $col, "Selection", $format); 
} elsif ($$hash_parameters{'marker'}{'behaviour'}[0] eq 'discovery') {
	$worksheet_parameters->write($row, $col, "Discovery", $format);
} 
$row++; $row++;	$col = $first_col;
$worksheet_parameters->write($row, $col, "PARAMETERS", $format); $col++;
if ($$hash_parameters{'marker'}{'polymorphism_user'}) {
	$worksheet_parameters->write($row, $col, "Polymorphic variants were used ", $format); $col++;
} else {
	$worksheet_parameters->write($row, $col, "Polymorphic variants were not used ", $format); $col++;
} $row++; $row++;

unless ($$hash_parameters{'marker'}{'behaviour'}[0] eq 'selection') {
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Parameters that can be modified only via the command-line version", $format_bold); $row++;
	$worksheet_parameters->write($row, $col, "Maximum percentage missing allowed (%) (MPA)", $format); $col++;
	my $x = $$hash_parameters{'marker'}{'missing_allowed'}[0] * 100; 
	$worksheet_parameters->write($row, $col, $x, $format_right); $row++;
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Significance Level Coverage Distribution (SLCD)", $format); $col++;
	$worksheet_parameters->write($row, $col, $$hash_parameters{"mapping"}{"level_significance_coverage_distribution"}[0], $format_right); $row++; 
	
	$row++; $row++; $col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Markers Features",$format_bold); $row++;
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Conserved Length (CL):", $format); $col++;
	my $merged_coord_CL = $$hash_parameters{'marker'}{'window_size_CONS_min'}[0]."--".$$hash_parameters{'marker'}{'window_size_CONS_max'}[0];
	$worksheet_parameters->write($row, $col, $merged_coord_CL, $format_right); $row++;
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Differences in the conserved region (CD):", $format); $col++;
	$worksheet_parameters->write($row, $col, $$hash_parameters{'marker'}{'window_var_CONS'}[0], $format_right); $row++;    
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Variable Length (VL):", $format); $col++;
	my $merged_coord_VL = $$hash_parameters{'marker'}{'window_size_VARS_min'}[0]."--".$$hash_parameters{'marker'}{'window_size_VARS_max'}[0];
	$worksheet_parameters->write($row, $col, $merged_coord_VL, $format_right); $row++;        
	$col = $first_col + 1; 
}

if ($$hash_parameters{'marker'}{'variable_divergence'}[0]) {
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Variable Divergence (%) (VD):", $format); $col++;
	$worksheet_parameters->write($row, $col, $$hash_parameters{'marker'}{'variable_divergence'}[0], $format_right); $row++;	
} else {
	$col = $first_col + 1; 
	$worksheet_parameters->write($row, $col, "Variable Positions (VP):", $format); $col++;
	my $merge = $$hash_parameters{'marker'}{'variable_positions_user_min'}[0]." -- ". $$hash_parameters{'marker'}{'variable_positions_user_max'}[0];
	$worksheet_parameters->write($row, $col, $merge, $format_right); $row++;    
}

$col = $first_col + 1; 
$worksheet_parameters->write($row, $col, "Minimum number of taxa covered (MCT) ", $format); $col++;
$worksheet_parameters->write($row, $col, $$hash_parameters{'marker'}{'MCT'}[0], $format_right); $row++; $row++;
$col = $first_col + 1;
$worksheet_parameters->write($row, $col, "Taxa used for marker discovery:", $format_bold); $col++;
$worksheet_parameters->write($row, $col, "ID", $format_bold); $row++;
my $counter = 1;
my $string = $$hash_parameters{'marker'}{'taxa_string'}[0];
if ($string) {
	my @array = split(",",$string);
	for (my $i=0; $i < scalar @array; $i++) {
		$col = $first_col + 1; 
		$worksheet_parameters->write($row, $col, $array[$i], $format_left); $col++;	
		my $j=$i+1;
		$worksheet_parameters->write($row, $col, $j, $format_right); $counter++;	$row++;
	} 
}
$workbook->close();

######
my $markers_count;
open (MARKER, "<$markers_file"); while (<MARKER>) { if ($_ =~ /.*\_\#.*/) {$markers_count++} } close (MARKER);
	
## Print results and instructions
print "\n\n"; DOMINO::printHeader("", "#"); DOMINO::printHeader(" RESULTS ","#"); DOMINO::printHeader("", "#");
if ($markers_count) { print "\n+ DOMINO has retrieved $markers_count markers\n"; }

my $marker_dirname = $$hash_parameters{"marker"}{"folder"}[0];
my $instructions_txt = $marker_dirname."/Instructions.txt";
my $MSA = $path."/MSA_markers";
open (OUT, ">$instructions_txt");
my $string2print = "+ Several files and folders has been generated:
\t+ $marker_dirname: contains DOMINO markers detected for each taxa as a reference and the clusterized results.
\t+ $path: contains the clusterized and definitive results for DOMINO markers.
\t+ $excel_woorkbook_name: contains information of the markers identified and parameters used by DOMINO.
\t+ $MSA folder contains a single file for each marker identified.\n\n";	
print OUT $string2print; print $string2print."\n"; close(OUT);

chdir $path; DOMINO::print_success_Step("excel");
exit();

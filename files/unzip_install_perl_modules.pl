#!/usr/bin/perl
use strict;
use warnings;
my $tar_folder = $ARGV[0];
my $abs_path_install = $ARGV[1];
if (scalar @ARGV != 2) {
	print "\nUsage:\n\tperl $0 modules-folder path\n\n";
	exit();
}

my $lib = $abs_path_install."/lib";
mkdir $lib, 0755; chdir $tar_folder;
my @exporter_lib_files;
my $ref_array_modules = &read_dir($tar_folder);
my $tmp_folder_copy = $tar_folder."/tmp_copy";
my @array_modules = @$ref_array_modules;
for (my $i=0; $i<scalar @array_modules; $i++) {
	if ($array_modules[$i] eq "." || $array_modules[$i] eq ".." || $array_modules[$i] eq ".DS_Store") {next;}
	chdir $tar_folder;
	if ($array_modules[$i] eq "Copy.pm") {
		mkdir $tmp_folder_copy, 0755;
		system("cp $array_modules[$i] $tmp_folder_copy"); 
	}
	if ($array_modules[$i] =~ /(.*)\.tar\.gz/) {
		print "\t- Installing module: ".$array_modules[$i]."\n";
		my $tmp_path = $tar_folder."/$1_tmp";
		mkdir $tmp_path, 0755; chdir $tmp_path;
		system ("ln -s $tar_folder/$array_modules[$i]");
		&untar_file($array_modules[$i]); 
		print "\n\n";

		# Compile and generate binaries of module
		my $ref_module = &read_dir($tmp_path);		
		my @module = @$ref_module;
		for (my $j=0; $j< scalar @module; $j++) {
			if ($module[$j] eq "." || $module[$j] eq ".." || $module[$j] eq ".DS_Store" ) {next;}			
			if ($module[$j] =~ /.*gz/) {
				next;
			} else {
				my $module_path = $tmp_path."/$module[$j]";
				chdir $module_path;
				my $ref_module_files = &read_dir($module_path);
				my @module_files = @$ref_module_files;
				my $make_bool=0;
				for (my $h=0; $h < scalar @module_files; $h++) {
					if ($module_files[$h] eq '.exists' || $module_files[$h] eq '.' || $module_files[$h] eq '..' || $module_files[$h] eq '.DS_Store') {next;}
					if ($module_files[$h] eq 'auto') {next;}
					if ($module_files[$h] eq "Makefile.PL") {
						system("perl Makefile.PL");	
						system("make > DOMINO_Error.log");
						$make_bool=1;	
				}}
				
				if ($make_bool == 1) {
					$ref_module_files = &read_dir($module_path."/blib/lib");
					@module_files = @$ref_module_files;
					for (my $h=0; $h < scalar @module_files; $h++) {
						if ($module_files[$h] eq '.exists' || $module_files[$h] eq '.' || $module_files[$h] eq '..' || $module_files[$h] eq '.DS_Store') {next;}
						if ($module_files[$h] eq 'auto') {next;}
						chdir $module_path;
						if ($module_files[$h] =~ /.*pm/) {
							system ("cp $module_files[$h] $lib");
						} else {
							my $lib_dir = $lib."/$module_files[$h]";
							my $sub_module_folder = $module_path."/blib/lib/".$module_files[$h];
							chdir $sub_module_folder;
							unless (-d $lib_dir) { mkdir $lib_dir, 0755; };
							my $ref_lib = &read_dir($sub_module_folder);
							my @ref_array = @$ref_lib;
							for (my $k=0; $k < scalar @ref_array; $k++) {
								if ($ref_array[$k] eq '.exists' || $ref_array[$k] eq '.' || $ref_array[$k] eq '..' || $ref_array[$k] eq '.DS_Store') { next; }	
								chdir $sub_module_folder;
								my $file2move = $sub_module_folder."/$ref_array[$k]";
								system ("cp -r $file2move $lib_dir");
					}}}
					chdir $tar_folder;
					system ("rm -rf $tmp_path");	
	}}}}
	print "\n";
}

my $file_module = $lib."/File";
my $copy_file = $tmp_folder_copy."/Copy.pm";
system ("cp $copy_file $file_module");

sub untar_file {
	my $file = $_[0];
	my $command = "tar -zvxf ".$file;
	system($command);
}

sub read_dir {
	my $dir = $_[0];
	opendir(DIR, $dir) or die "Could not open $dir";
	my @dir_files = readdir(DIR);
	my $array_ref = \@dir_files;
	return $array_ref;
}
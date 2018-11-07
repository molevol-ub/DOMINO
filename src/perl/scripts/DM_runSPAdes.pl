#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/../lib";
BEGIN {
	require DOMINO;
}

my $path = $ARGV[0];
my $step_time = $ARGV[1];

# Get memory
if(!defined($path)) { print "ERROR: No input files are provided: [$0]\n"; exit; } ## give error message for DOMINO debug
my $domino_version ="DOMINO v1.1 ## Revised 24-10-2018";

my $hash_parameters = DOMINO::get_parameters($path);
my $domino_files_ref = DOMINO::get_DOMINO_files($path);
my $assembly_directory = $$hash_parameters{"assembly"}{"dir"}[0];
my $noOfProcesses=$$hash_parameters{"assembly"}{"CPU"}[0];
my $memory_server_file = $assembly_directory."/memory_server.txt";
my $scripts_path = $FindBin::Bin."/../";

#&debugger_print("cat /proc/meminfo > memory_server.txt");
system("cat /proc/meminfo > memory_server.txt");	
my $total_available;
unless (-e -r -s $memory_server_file) {
	DOMINO::printErrir("There was an error when retrieving Memory RAM information...\n\nAre you sure this is a linux sever?...");
	DOMINO::print_fail_Step("spades"); exit();
}
DOMINO::printHeader(" Memory RAM Usage retrieval ","#");
my %memory_hash;
open (MEM, $memory_server_file);
while (<MEM>) {
	my $line = $_;
	chomp $line;
	$line =~ s/\s*//g;			
	#&debugger_print($line);			
	my @array = split("\:", $line);
	if ($array[1] =~ /(\d+)(.*)/) {
		push (@{ $memory_hash{$array[0]} }, $1);
		push (@{ $memory_hash{$array[0]} }, $2);
	}
	if ($array[0] eq "Cached") { last; }
} close (MEM);
#&debugger_print("Ref", \%memory_hash);
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
my $amount_taxa = scalar keys %{$domino_files_ref};
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
	DOMINO::printHeader(" SPAdes Assembly of each taxa ","#");
	## Call spades for the assembly and send threads for each taxa

	my $int_taxa = 0;
	my $pm =  new Parallel::ForkManager($subProcesses); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
	$pm->run_on_finish( 
	sub { my ($pid, $exit_code, $ident) = @_; 
		print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; 
	} );
	$pm->run_on_start( sub { my ($pid,$ident)=@_; print "\n\n** SPAdes assembly started with PID $pid\n\n"; } );

	foreach my $keys (sort keys %{$domino_files_ref}) {
		next if ($keys eq 'original'); 
		my %files_threads;
		my $pid = $pm->start($int_taxa) and next; print "\nSending child process for SPAdes assembling $keys\n\n";
		$int_taxa++;

		## Send SPAdes command				
		my $assembly_dir = $$domino_files_ref{$keys}{'DIR'}[0];
		my $spades_path = "python ".$scripts_path."SPAdes-3.8.1-Linux/bin/spades.py -o ".$assembly_dir." ";
		if ($$hash_parameters{"assembly"}{"file_type"}[0] == 7) { ## illumina_PE
			$spades_path .= "-1 ".$$domino_files_ref{$keys}{'original'}[0]." -2 ".$$domino_files_ref{$keys}{'original'}[1];				
		} else { 
			$spades_path .= "-s ".$$domino_files_ref{$keys}{'original'}[0];
		}
		$spades_path .= " -t $noOfProcesses_SPAdes";
		$spades_path .= " > spades.log"; ## discarding SPAdes output		
		#&debugger_print("Sending command: $spades_path\n");
		print "Sending SPAdes command...\n";
		my $system_call = system($spades_path);
		if ($system_call != 0) {
			&printError("Something happened when calling SPAdes for assembly reads..."); DOMINO::dieNicely();
		} print "\n"; &time_log();

		## Get contig file
		my $contigs_file = $assembly_dir."/contigs.fasta";
		my $new_contigs_file = $$hash_parameters{"assembly"}{"folder"}[0]."/assembly_id-".$keys.".contigs.fasta";
		$files_threads{$keys}{'FINAL'} = $new_contigs_file;
		File::Copy::move($contigs_file, $new_contigs_file);

		## Finish assembly and generate statistics
		my $stats = DOMINO::Contig_Stats($new_contigs_file); 
		$files_threads{$keys}{'stats'} = $stats;

		#&debugger_print("DOMINO threads Assembly files"); &debugger_print("Ref", \%files_threads);
		my $spades_dump_hash = $assembly_dir."/dumper_files_threads.txt";
		DOMINO::printDump(\%files_threads, $spades_dump_hash);
		
		$pm->finish($int_taxa); # pass an exit code to finish
	}
	$pm->wait_all_children; print "\n** All Assembly child processes have finished...\n\n";	
} else {
	DOMINO::printError("Only $total_available GiB available of memory RAM, please bear in mind SPAdes would not complete the task..."); 
	DOMINO::print_fail_Step("spades");
	exit();
}

DOMINO::print_success_Step("spades");

sub time_log {	
	my $step_time_tmp = DOMINO::time_log($step_time); print "\n"; 
	$step_time = $$step_time_tmp;
}

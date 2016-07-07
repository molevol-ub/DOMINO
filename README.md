=
DOMINO: Development Of Molecular markers In Non-model Organisms 
=

=
1. INTRODUCTION
=
DOMINO is a new software application specifically designed for improving the development of DNA markers in non-model organisms using either NGS data or pre-computed multiple sequence alignments (MSA) in various formats (including RAD loci). We used Perl and C++ scripting languages to combine some popular, open-source, bioinformatics tools for NGS data with a set of new developed functions in an integrated bioinformatics pipeline for custom-made marker discovery or selection based on user-defined criteria. 

Customizable features include the length of variable and conserved regions (when requested), the minimum levels (or a preferred range) of nucleotide variation, how to manage polymorphic variants, or which taxa (or what fraction of them) should be covered by the marker. All these criteria can be easily defined in a user-friendly GUI or under a command-line version that implements some extended options and that it is particularly useful for working with large NGS data sets in high performance computers (see also the DOMINO documentation). The regions identified or selected in DOMINO can be i) directly used as markers with a particular depth of taxonomic resolution, ii) utilized for their downstream PCR amplification in a broader taxonomic scope or iii) used as suitable templates to optimized bait design for target DNA enrichment techniques.

The DOMINO workflow consists in four main phases: “Input Data” and “Pre-processing”, “Assembly”, “Alignment/Mapping” and “Marker Discovery/Selection”. Users can enter to the DOMINO workflow at different points in order to perform specific shortcut runs for marker discovery or selection according to the type of input data. In a typical full run from raw NGS reads (using marker identification module), the program applies an assembly-based approach; the pipeline is therefore optimized to work with genome partitioning methods in which the size-selected or enriched fragments are long enough to permit a successful assembly and have been fully (or nearly fully) sequenced. For MSA inputs, as for example, restriction-site associated DNA (RAD) variation and related methods (e.g. RADseq, ddRAD, GBS) DOMINO can search the MSA of each loci (in this case RAD loci) previously obtained with commonly used software packages for highly informative markers (marker selection module); current version accepts the PHYLIP and FASTA format as well as the pyRAD (.loci) and STACKS (.fa) output file formats. 

After the marker development phase, DOMINO provides i) a list of the genomic regions (with their corresponding coordinates in each generated contig) that meet marker design criteria and the MSA of each of these markers MSAs for the specified taxa (DOMINO marker identification module outputs) or ii) the list of makers selected under the same criteria among a set of aligned sequences supplied by the user and previously obtained in other software (DOMINO marker selection module outputs). 

=
2. DOCUMENTATION
=
Documentation and additional information can be retrieved from the DOMINO oficial website: http://www.ub.edu/softevol/domino/

=
3. INSTALLATION
=
3.1. Operating Systems

DOMINO is developed and maintained under GNU/Unix using ANSI/ISO C and standard development tools (GCC, Autotools).

It can be compiled and run in Linux (Ubuntu 14.04 LTS) and Mac OSX systems (version > 10.10)

3.2. Requirements

- Perl Programming Language version v5.2+

- C++ compiler
	+ Mac OS users may need to install Xcode package via App Store
	+ Linux users may need to install build essential package o g++ compiler
		Ubuntu:
			sudo apt-get install build-essential 
			or 
			sudo apt-get install g++

- zlib compression library v1.2.3+ (SAMtools dependency) <http://www.zlib.net/>
	+ Mac OS: included and installed in the distribution.
	+ Linux users may need to install zlib1g-dev package
		Ubuntu/Linux Mint: sudo apt-get install zlib1g-dev
			
- Xterm (in some ubuntu-like distributions may not be included):
	+ Mac OS: included and installed in the distribution
	+ Ubuntu/Linux Mint: sudo apt-get install xterm


3.3. Installation

Two possible installation options:

3.3.1. DOMINO GUI:

In the DOMINO website, user can obtain installers that will guide the user through the complete installation process for Linux and Mac distributions.

Also, all the necessary installation files and folders (the user may need to compile and retrieve some Qt libraries) are available at the domino Git repository.

Once installed, the bin directory will contain the DOMINO perl scripts, the directories with the binaries of the external software included in the DOMINO pipeline and some mandatory Perl modules. In addition, a shortcut Desktop DOMINO icon will be generated.

3.3.2. Command-line:

In the Git directory, the user can find a pre-configured shell script (install.sh) written to handle the installation of command-line version (to install only the command-line version).

In a Unix-based bash console, change the directory to enter in the desired installation folder and run the script, with the command: sh install.sh

The bin directory will contain the DOMINO perl scripts, the directories with the binaries of the external software included in the DOMINO pipeline and some mandatory Perl modules.

=
4. BUG REPORTS
=

If you encounter bugs that are not listed here, please send comments and bug reports via GitHub


=
5. COMMUNITY
=

If you would like to see a new feature implemented, suggestions are welcome. Our aim is to have a "suggestion" box soon.

=
6. THANKS
=

The authors would like to thank all people who contributed (and are still contributing) in the creation and testing of this software. 

=
7. AUTHORS:
=
Cristina Frias-Lopez, Jose F. Sanchez-Herrero, Miquel A. Arnedo, Alejandro Sanchez-Gracia and Julio Rozas.
  	
   	C.Frias-Lopez: cristinafriaslopez@ub.edu
   	J.F.Sanchez-Herrero: jfsanchezherrero@ub.edu
	M.A. Arnedo: marnedo@ub.edu
	A.Sanchez-Gracia: elsanchez@ub.edu
   	J.Rozas: jrozas@ub.edu 
	
Evolutionary Genomics and Bioinformatics Group, Departament de Genètica, Microbiologia i Estadística and Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Av. Diagonal 643, Barcelona 08028, Spain

=
8. CITATION
=
*To add citation when available. Current status: under review*

=
9. COPYRIGHT AND LICENSE
=
Copyright (C) 2016 Evolutionary Genomics and Bioinformatics Group, University of Barcelona.

DOMINO is licensed under the GPLv3 license.  See `LICENSE' file for details. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>

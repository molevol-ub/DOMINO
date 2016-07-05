=
DOMINO: Development Of Molecular markers In Non-model Organisms 
=

=
1. INTRODUCTION
=
DOMINO is a new software application specifically designed for improving the development of DNA markers in non-model organisms using either NGS data or pre-computed multiple sequence alignments (MSA) in various formats (including MSA from RAD data). 

We used Perl and C++ scripting languages to combine some popular, open-source, bioinformatics tools for NGS data with a set of new developed functions in an integrated bioinformatics pipeline for custom-made marker discovery or selection. 

Using NGS data (raw or pre-processed reads) from typical genome partitioning or low depth sequencing experiments (including or not reference sequences, such as, for example, genome scaffolds), or pre-computed multiple sequence alignments (MSA) of multiple taxa 
(“taxa panel”, which may represent different populations or species), DOMINO searches for sequence regions with a minimum level of variation, which can be optionally flanked by stretches of conserved sequences (appropriate for further PCR primers design). 

The DOMINO workflow consists in four main phases: “Input Data” and “Pre-processing”, “Assembly”, “Alignment/Mapping” and “Marker Discovery/Selection”. Users can enter to the DOMINO workflow at different points in order to perform specific shortcut runs for marker 
discovery or selection according to the type of input data.

In a typical full run from raw NGS reads (using marker identification and marker selection modules), the program applies an assembly-based approach; the pipeline is therefore optimized to work with genome partitioning methods in which the size-selected or enriched fragments are long enough to permit a successful assembly and have been fully (or nearly fully) sequenced. For MSA inputs, as for example, restriction-site associated DNA (RAD) variation and related methods (e.g. RADseq, ddRAD, GBS) DOMINO can searh the MSA of each loci (in this case RAD loci) previously obtained with commonly used software packages for highly informative markers (marker selection module); current version accepts the PHYLIP and FASTA format as well as the pyRAD (*.loci) and STACKS (*.fa) output file formats. 

After the maker development phase, DOMINO provides i) the list the genomic regions (with their corresponding coordinates in each contig) exhibiting the minimum levels of nucleotide variation specified by the user among the selected taxa (or covering at least a pre-defined minimum number of them), and the marker MSAs for these taxa (DOMINO marker identification module) or ii) the list of markers exhibiting the same pre-definable features among a set of aligned sequences supplied by the user and previously obtained in other software (DOMINO marker selection module). 

Furthermore, if the taxa panel has been designed in a convenient phylogenetic context and the user asks for highly conserved regions flanking to be included in the developed markers, these markers should be suitable for further surveys of variation in an extended set of phylogenetically related taxa, i.e. “focal taxa”. 

=
2. DOCUMENTATION
=
Documentation and additional information can be retrieved from the DOMINO oficial website: http://www.ub.edu/softevol/domino/

=
3. INSTALLATION
=
3.1. Supported Operating Systems

DOMINO is developed and maintained under GNU/Unix using ANSI/ISO C and standard development tools (GCC, Autotools).

It can be compiled and run in Linux (it has been tested in Ubuntu 14.04 LTS) and Mac OSX systems (version > 10.10)

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
		Ubuntu/Linux Mint:
			sudo apt-get install zlib1g-dev
			
- Xterm (in some ubuntu-like distributions may not be included):
	+ Mac OS: included and installed in the distribution
	+ Ubuntu/Linux Mint:
		sudo apt-get install xterm


3.3. Installation

Two possible installation options:

3.3.1. DOMINO GUI:

Under the DOMINO oficial website, user can obtain installers that will guide the user through the complete installation process for Linux and Mac distributions. 

Also, all the necessary installation files and folders are available at the domino git repository but the user may need to compile and retrieve some Qt libraries. 

Once installed, the bin directory will contain the DOMINO perl scripts, the directories with the binaries of the external software included in the DOMINO pipeline and some mandatory Perl modules. In addition, a shortcut Desktop DOMINO icon will be generated.

3.3.2. Command-line:

In the git directory, the user can find a pre-configured shell script (install.sh) written to handle the installation of command-line version (please, notice that only the command-line version should be installed using this option).

In a Unix-based bash console, change the directory to fit the desired installation folder and run the script, with 'sh install.sh'
	
The bin directory will contain the DOMINO perl scripts, the directories with the binaries of the external software included in the DOMINO pipeline and some mandatory Perl modules.

------------------------------
4. BUG REPORTS
------------------------------

If you encounter bugs that are not listed here, please send a report to the main authors.
  	
   	J.F.Sanchez: jfsanchezherrero@ub.edu
	A.Sanchez-Gracia: elsanchez@ub.edu
   	J.Rozas: jrozas@ub.edu 


------------------------------
5. COMMUNITY
------------------------------

If you would like to see a new feature implemented, send suggestions to the main authors (See AUTHORS). Our aim is to have an "asked features" framework soon.

------------------------------
6. THANKS
------------------------------

The authors would like to thank all people who contributed (and are still contributing) in the creation and testing of this software. 

------------------------------
7. AUTHORS:
------------------------------
Cristina Frias-Lopez, Jose F. Sanchez-Herrero, Miquel A. Arnedo, Alejandro Sanchez-Gracia and Julio Rozas.
	
Evolutionary Genomics and Bioinformatics Group, Departament de Genètica, Microbiologia i Estadística and Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Av. Diagonal 643, Barcelona 08028, Spain


------------------------------
8. CITATION
------------------------------
*To add citation when available*

*Current status: under review*


------------------------------
9. COPYRIGHT AND LICENSE
------------------------------

Copyright (C) 2016 Evolutionary Genomics and Bioinformatics Group, University of Barcelona.

DOMINO is licensed under the GPLv3 license.  See `LICENSE' file for details. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>

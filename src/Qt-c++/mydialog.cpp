#include "mydialog.h"
#include "ui_mydialog.h"
#include <QtCore>
#include <iostream>
#include <QGraphicsPixmapItem>
#include <QGraphicsView>
#include <QLabel>

using namespace std;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ General Dialog Class: Standard Output and file information @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

MyDialog::MyDialog(QWidget *parent) : QDialog(parent), ui(new Ui::MyDialog) {
    // Set dialog window itself
    ui->setupUi(this);
    setModal(false);
}

MyDialog::~MyDialog() {
    delete ui;
}

void MyDialog::DOMINOlicense() {
    QString dir = QDir::currentPath();
    QString License;
    License = QString(dir + "/scripts/LICENSE.txt");
    QFile file(License);
    if(!file.open(QIODevice::ReadOnly)) {
        cout << "Error: LICENSE.txt file could not be read...\n";
    }
    QTextStream in(&file);
    while(!in.atEnd()) {
        QString line = in.readLine();
        MyDialog::textFileBrowser_dialog(line);
    }
    file.close();
}

void MyDialog::aboutDOMINO() {
    QString aboutDOMINO_guide = (

    "OVERVIEW\n\n"

    "DOMINO is a new software application specifically designed for improving the development of DNA markers in non-model organisms using either NGS data or pre-computed multiple sequence alignments (MSA) in various formats (including MSA from RAD data).\n\n"

    "We used Perl and C++ scripting languages to combine some popular, open-source, bioinformatics tools for NGS data with a set of new developed functions in an integrated bioinformatics pipeline for custom-made marker discovery or selection.\n\n"

    "Using NGS data (raw or pre-processed reads) from typical genome partitioning or low depth sequencing experiments (including or not reference sequences, such as, for example, genome scaffolds), or pre-computed multiple sequence alignments (MSA) of multiple taxa (“taxa panel”, which may represent different populations or species), DOMINO searches for sequence regions with a minimum level of variation, which can be optionally flanked by stretches of conserved sequences (appropriate for further PCR primers design).\n\n"

    "The DOMINO workflow consists in four main phases: “Input Data” and “Pre-processing”, “Assembly”, “Alignment/Mapping” and “Marker Discovery/Selection”. Users can enter to the DOMINO workflow at different points in order to perform specific shortcut runs for marker discovery or selection according to the type of input data.\n\n"

    "In a typical full run from raw NGS reads (using marker identification and marker selection modules), the program applies an assembly-based approach; the pipeline is therefore optimized to work with genome partitioning methods in which the size-selected or enriched fragments are long enough to permit a successful assembly and have been fully (or nearly fully) sequenced. For MSA inputs, as for example, restriction-site associated DNA (RAD) variation and related methods (e.g. RADseq, ddRAD, GBS) DOMINO can searh the MSA of each loci (in this case RAD loci) previously obtained with commonly used software packages for highly informative markers (marker selection module); current version accepts the PHYLIP and FASTA format as well as the pyRAD (*.loci) and STACKS (*.fa) output file formats.\n\n"

    "After the maker development phase, DOMINO provides i) the list the genomic regions (with their corresponding coordinates in each contig) exhibiting the minimum levels of nucleotide variation specified by the user among the selected taxa (or covering at least a pre-defined minimum number of them), and the marker MSAs for these taxa (DOMINO marker identification module) or ii) the list of markers exhibiting the same pre-definable features among a set of aligned sequences supplied by the user and previously obtained in other software (DOMINO marker selection module).\n\n"

    "Furthermore, if the taxa panel has been designed in a convenient phylogenetic context and the user asks for highly conserved regions flanking to be included in the developed markers, these markers should be suitable for further surveys of variation in an extended set of phylogenetically related taxa, i.e. “focal taxa”.\n\n"

    "AUTHORS\n\n"

     "Cristina Frias-Lopez, Jose F. Sanchez-Herrero, Miquel A. Arnedo, Alejandro Sanchez-Gracia and Julio Rozas.\n\n"

     "Evolutionary Genomics and Bioinformatics Group, Departament de Genètica, Microbiologia i Estadística and Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Av. Diagonal 643, Barcelona 08028, Spain\n\n"


     "BUGS & COMMENTS:\n\n"

     "If you encounter bugs that are not listed here, please send a report to the main authors.\n\n"

     "J.F.Sanchez: jfsanchezherrero@ub.edu\n"
     "A.Sanchez-Gracia: elsanchez@ub.edu\n"
     "J.Rozas: jrozas@ub.edu\n\n\n"

     "COPYRIGHT\n\n"
     "Copyright (C) 2016 Evolutionary Genomics and Bioinformatics Group, Universitat Barcelona.\n\n"

     "LICENSE\n\n"
     "DOMINO is licensed under the GPLv3 license.  See `LICENSE' file for details.\n\n"

      "This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\n\n"

      "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n"

      "See the GNU General Public License for more details.\n\n");

    MyDialog::show_helpguide(aboutDOMINO_guide);
}

void MyDialog::perl_dialog() {
    // Set window title
    MyDialog::QWidget::setWindowTitle("Log details of the process");
    //ui->label_DOMINO_workflow->hide();
}

void MyDialog::write_perl_output(QByteArray ba_outsring) {
    QString outstring;
    outstring.append(ba_outsring);
    ui->textBrowser_1->append(outstring);
}

void MyDialog::textFileBrowser_dialog(QString file_data) {
    ui->textBrowser_1->append(file_data);
}

void MyDialog::textFileBrowser_dialog_show(QString file_title) {
    // Set Window Title
    MyDialog::QWidget::setWindowTitle(file_title);
}

void MyDialog::show_helpguide(QString helpguide) {
    // Set sentence to show when asking for help in the cleaning process
    ui->textBrowser_1->setText(helpguide);
}

void MyDialog::about_cleaning_pipeline() {

   QString helpguide_cleaning = (

    "DOMINO accepts input sequence data files in two different formats, the 454 Pyrosequencing Standard Flowgram Format (SFF), and FASTQ format. These input files can contain 454 or Illumina (single or paired-end) raw reads from m taxa (taxa panel). The sequences from each taxon should be properly identified with a specific barcode (aka, tag, MID or index), or loaded in separate files, also appropriately named (see the DOMINO manual in the DOMINO Web site for details). DOMINO is designed to filter low quality, low complexity, contaminant and very short reads using either default or user-specified filtering parameters. Mothur, PRINSEQ, NGS QC toolkit, BLAST, as well as new Perl functions specifically written for DOMINO (DM scripts) are used to perform these tasks.\n\n"

    "DOMINO uses Mothur v1.32.0 to extract reads from SFF files and store them in FASTQ format, which are subsequently converted to FASTA and QUAL files. Low quality or very short reads are trimmed, or definitely removed, using NGS QC Toolkit v2.3.1. PRINSEQ v0.20.3 package is used to eliminate low complexity reads using the implemented DUST algorithm.\n\n"

    "Putative contaminant sequences, such as bacterial DNA frequently found in genomic samples, cloning vectors, adaptors, linkers, and high-throughput library preparation primers, can also be removed using a DOMINO function that performs a BLAST search (BLAST v2.2.28) against UniVec database (http://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) and/or against a user-supplied contaminant database (see the DOMINO manual).\n\n"
);
      MyDialog::show_helpguide(helpguide_cleaning);
}

void MyDialog::about_assembling_pipeline() {
    QString helpguide_assembling = (
    "DOMINO performs separate assemblies, one for each panel taxon, using MIRA v4.0.2 with the pre-processed reads from the previous step or with those supplied by the user.\n\n"

    "Although the default parameter values vary in function of the particular sequencing technology, the majority of them are shared (see the DOMINO manual).\n\n"

    "In order to avoid including repetitive and chimeric regions, all contigs (and the corresponding reads) identified as HAF6, HAF7 and MNRr by the MIRA algorithm are discarded from the mapping/alignment phase.\n\n"

    "Since MIRA can generate redundant contigs because of polymorphic and paralogous regions, we have implemented a specific DOMINO function that performs a similarity clustering of all contigs (i.e., an all versus all contigs BLAST search) to identify and remove such redundancies. All positive hits (E-value < 10-50), with an aligned region length higher or equal than the 85% of the shorter contig, and a minimum percentage of similarity of the 85%, are collapsed to the longest contig of the pair for the mapping phase.\n\n"

    "The DOMINO command line version also includes an option to perform a second iterative assembly step using the software CAP3. If selected, this option uses MIRA output sequences (contigs and singletons) as input for CAP3 under a relaxed parameter scheme: overlapping consensus similarity of 80 bp and percent identity of 98%. Users can modify these default parameters in the command prompt.\n\n"
    );
    MyDialog::show_helpguide(helpguide_assembling);
}

void MyDialog::about_DOMINO_marker_development() {
    QString helpguide_DOMINO_development = (

                "Mapping/Alignment phase\n\n"

                "DOMINO uses Bowtie2 to map the pre-processed reads from each taxon to the assembled contigs of the other m-1 taxa from the panel. Thus, in this step, DOMINO builds m(m-1) sequence alignment/map files (SAM/BAM files). Using the contigs from the assembly of taxon#1 as reference sequences, DOMINO will achieve 3 different alignment/mapping runs (reads of taxon#2 to taxon#1 references; reads of taxon#3 to taxon#1 references; reads of taxon#4 to taxon#1 references). The program will proceed in the same way with the assemblies of taxon#2, taxon#3 and taxon#4. In the case of a panel of m = 4 taxa, for example, DOMINO will build 4 x 3 = 12 SAM/BAM files during this step. The reason behind this particular mapping strategy lies in the dissimilar performance of alignment/mapping algorithms depending on the divergence between the reads and the reference sequences.\n\n"

                "Mapping errors and contigs with uneven sequencing depth\n\n"

                "Immediately after generating BAM files, DOMINO removes all unmapped contigs and multi-mapping reads. This step is critical to avoid alignment artifacts, which can create false positive markers (i.e., sequence regions with misleading high levels of nucleotide diversity). The contigs with an unusually large number of aligned reads, which can correspond to repetitive regions, are also removed (they are not suitable for designing single copy markers); in particular we discarded all contigs where the probability of observing the sequenced coverage in one or more positions was < 10-5 (from a Poisson with mean c, where c is the average per contig coverage estimated using all contigs in the sample of the same taxon); this cutoff values, nevertheless, could be changed in the command-line version. Later, DOMINO will build one pileup file per each BAM/SAM file using the SAMtools v0.1.19 suite (mpileup option).\n\n"

                "Sequencing errors and ambiguity codes\n\n"

                "Since sequencing errors might have a great effect on the marker selection, DOMINO incorporates their own functions for detecting and masking putative sequencing errors, which apply a very conservative criterion for variant calling. First, to avoid the calling of spurious nucleotide variants in low sequencing coverage experiments (i.e. erroneously assigned variants fixed between the taxa from the panel), DOMINO mask the information from positions with only one read mapped to the reference.\n\n"

                "Furthermore, sequencing errors may also inflate the number of called polymorphisms under the Polymorphic Variants option in the marker Discovery/Selection phase. To avoid such undesirable effect, DOMINO incorporates a similar conservative criterion to use only highly credible polymorphisms. Under the Polymorphic Variant option, DOMINO will assume that each taxon represents a diploid individual; that is, in a true heterozygous position it expects approximately the same number of reads from each of the two segregating alleles. For positions with eight or more reads mapped, DOMINO discards those polymorphic variants in which the frequency of the minor allele is significantly lower than the expected (P < 0.05 in Binomial distribution with p = 0.5), likely corresponding to a sequencing error. For lower coverage values, DOMINO will use the information of a polymorphic variant only if the allele with the minor frequency is present in two or more reads. This testing procedure, applied independently for each position within each species, will likely discard many true polymorphic sites; this variant calling approach, however, makes DOMINO highly conservative in detecting true markers when including polymorphisms in the analysis (i.e., DOMINO will use only highly confident within-species segregating variants for the marker Discovery/Selection phase).\n\n"

                "Ambiguity codes, either introduced by MIRA assembler in contig sequences or present in user-supplied reference sequences or MSA, are considered by DOMINO to decide whether a position is or not variable. For instance, if there is a Y (IUPAC ambiguity code for C/T) in the reference sequence and all aligned reads show an A at this position, the program computes a nucleotide change (because it is probably a fixed position between taxa). Instead, if all reads show a C or a T, DOMINO does not consider this position as variable in marker Discovery/Selection phase.\n\n"

                "After applying all the above-mentioned post-mapping filters, DOMINO combines the variation profiles (arrays with the information about the state of each position, conserved or variable between taxa pairs) obtained from each of the m-1 pileup files including the same reference sequence (i.e., the same taxon), into a single multiple taxa variation profile (MTP). Since each of these references will be likely fragmented in i contigs, DOMINO will build i x m MTP per taxon. Each of these MTP will be independently scanned for regions containing candidate markers in the next phase.\n\n"

                "If the user provides reference sequences from a single taxon (e.g. a genome draft), plus the reads from the m different taxa, the program builds only one MTP set (one per contig or scaffold in the supplied reference). This option allows mapping short reads from RAD sequencing (e.g., from RADSeq and related methods) to a reference, generating the variation profiles of each RAD region (RAD loci) for the next marker Discovery/Selection phase.\n\n"

                "On the other hand, if the input includes a single or multiple pre-computed MSA instead of NGS data, DOMINO skips the alignment/mapping phase and directly generates the single MTP set (one per aligned region). In this point, the program accepts MSA files in FASTA (multiple FASTA files, one per linked region), PHYLIP (multiple PHYLIP files, one per linked region, or one multi PHYLIP file with the alignment of all regions) and pyRAD LOCI (*.loci files generated by the program pyRAD) and STACKS fasta (batch_X.fa output files generated from the population analyses in the program STACKS) output files.\n\n"

                "Marker Discovery/Selection phase\n\n"

                "Each MTP generated in the previous step is either scanned for the presence of candidate marker regions using a sliding window approach (DOMINO marker identification module) or ii) used to select the markers with the desired features among a set of pre-computed MSA loaded in the previous TAB (DOMINO marker identification module). In the first case, a specific DOMINO function searches for sequence regions of desired length (Variable region Length, VL), showing the minimum level of variation indicated by the user (Variable region Divergence, VD). DOMINO can also restrict that this variable region was flanked (or not) by highly conserved regions (Conserved region Divergence, CD) of a predefined length (Conserved region Length, CL); an information useful to further design PCR primers. Moreover, DOMINO can strictly restrict the search to a particular set of taxa (from the panel), or just specify the minimum number of taxa required to be covered by the marker (by changing the Minimum number of Covering Taxa parameter; MCT < m). As indicated, DOMINO can use or not the information from polymorphic sites. An appropriated combination of selected taxa and MCT and VD parameter values will allow the user select a large set of informative markers suitable for addressing a wide range of evolutionary questions. In the marker selection module, DOMINO allows directly selecting the most informative markers among the loaded by the user in the same way and with the same personalized features described above. For RAD loci, a particular range of variable positions (VP) between the closest taxa instead of the VD parameter must be specified. The latter option allows selecting informative RAD loci while excluding those exhibiting anomalous high levels of variation, which might reflect RAD tag clustering errors.\n\n"

                "Since DOMINO can work with more than one MTP set (m in a full DOMINO marker identification run), some of the markers found in MTP based on different reference taxa may be redundant (they can cover the same genomic region, although with different coordinates; see Mapping/Alignment phase section), while other can be found only in one particular profile. To avoid reporting redundant information, we have implemented a function to collapse (using BLAST) these maker sequences, only reporting unique markers. In addition to the MSAs of each identified marker for the selected taxa, DOMINO reports the sequence of all contigs containing candidate markers, along with all relevant information, such as the exact positions delimiting the conserved and variable regions and divergence levels among the selected taxa in the variable region.\n\n"

                "To maximize the probability of finding informative markers, the final list of candidates can include overlapped regions that fulfill the specified characteristics (not applicable to the DOMINO marker selection module). Operationally, all regions that meet the criteria for being considered a candidate marker (after moving the scanning window five or more base pairs) are listed as different markers in the final output. In this way, users can choose the best marker to be used directly for further analyses, or the most appropriated region of each contig to be PCR amplified and sequenced in additional focal species (i.e., the best marker from each linked block).\n\n");
    MyDialog::show_helpguide(helpguide_DOMINO_development);
}

void MyDialog::on_pushButton_clicked() {
    close();
}

void MyDialog::printing_DOMINO_last_step_finished() {
    QString instructions_finish = (
    "DOMINO is finished here succesfully\n\n"
    "Several files and folders has been generated:\n"
    "+ DOMINO_markers_parameters_details.txt: contains information and parameters regarding the molecular markers developed with the characteristics desired.\n"
    "+ DOMINO_alignment_using_CleanReads Folder: contains the alignment/mapping information for each taxa.\n"
    "+ DOMINO_markers_using_CleanReads_... Folder and subfolders: contains DOMINO markers detected. For each taxa specified, DOMINO would search for putative molecular markers using each taxa as a reference each time. This folder would be named as 'markers_Reference_*'.\n"

    "\nOnce all the taxa have been used as reference, DOMINO would clusterized the information and markers found. In the directory 'clustering' several files regarding the BLAST search for clustering are presented although, the clusterized markers would appear in 'DOMINO_markers_results'.\n"
    "In each subfolder, the files 'DOMINO_shared_contigs' ('.txt' or '.xls') contain information of the coordinates of the different molecular markers identified. The sequences of these markers could be found in the file 'contig_sequences_DOMINO_markers.fasta', here we would find all the sequence and not just the trimmed markers identified. Those markers discarded because of not meeting the requirements specified, would be found in the file 'Error_contigs.txt'.\n"
    "Several other files named according to the species used in the development of molecular markers, would contain information necessary for obtaining a visual representation of the contigs using an appropiate software. So far, if we are using a contig viewer such as Geneious or Tablet, when importing an assembly, we would need to upload the file with 'sorted.bam' termination, bear in mind that the file 'sorted.bam.bai' should also be present in the same folder.\n");
    MyDialog::show_helpguide(instructions_finish);
}

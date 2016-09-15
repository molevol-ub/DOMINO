#include "mainwindow_amm.h"
#include "ui_mainwindow_amm.h"
#include "mydialog.h"
#include "change_progress_bar.h"
#include <QProcess>
#include <iostream>
#include <QString>
#include <QTextEdit>
#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QFileDialog>
#include <QByteArray>
#include <QTabWidget>
#include <QTreeView>
#include <QModelIndex>
#include <QDirModel>
#include <QDesktopServices>
#include <QInputDialog>
#include <QListView>
#include "results_dialog.h"
#include "checkOS.h"
#include <cmath>
#include <QMutex>
using namespace std;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%       DOMINO Main Window Class        %%
// %%                                       %%
// %%       Author:                         %%
// %%   Jose F. Sanchez Herrero             %%
// %%   jfsanchezherrero@ub.edu             %%
// %%                                       %%
// %%       Institution:                    %%
// %%   Evolutionary Genomics &             %%
// %%   Bioinformatics Group.               %%
// %%   University Barcelona                %%
// %%                                       %%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// @@@@@@@@@@@@@@@@@@@@@@@
// @@ General Variables @@
// @@@@@@@@@@@@@@@@@@@@@@@
extern QString APPDIR;
QString input_file_name_cleaning;
QString input_file_name_assembling;
QString input_file_name_DOMINO_marker;
QString input_file_name_DOMINO_marker_2;
QString outfolder;
QString barcode_user_file;
QString additional_perl_options_cleaning;
QString additional_perl_options_assembling;
QString additional_perl_options_markers;
QString additional_perl_options_mapping;
QStringList db_list;
QStringList type_files;
int num_proc; int thresh_dust; int min_length; int length_cutoff; int PHRED_sc; int bdiffs;
extern QString perl_path;
QString perl_abs_path;
bool debug_mode;
bool cleaning_process_killed;
bool assembling_process_killed;
bool marker_identification_killed;
bool marker_identification_finished;
bool type_file_cleaning;
bool type_file_assembling;
bool type_file_mappingMarker;
bool domino_short_cut;
bool cleaning_checked;
bool assembly_checked;
bool mapping_checked;
bool assembly_finished;
bool DOMINO_user_files;
bool names_checked;
int idealthreads (QThread::idealThreadCount()); // Get the ideal number of proccessors to use for this computer

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Initializer and destroyer @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
MainWindow_AMM::MainWindow_AMM(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow_AMM) {
    // ###############################################
    // ## Set some icons and the application itself ##
    // ###############################################
    ui->setupUi(this);

    // Set some buttons, tip and help information and icons
    ui->show_workflow_button->setToolTip("< img src = ':/images/GUI_icons/information.png'> Show DOMINO Workflow");
    ui->pushButton_showManual_intro_tab->setToolTip("< img src = ':/images/GUI_icons/information.png'> Show DOMINO Manual");
    ui->pushButton_showRef_intro_tab->setToolTip("< img src = ':/images/GUI_icons/information.png'> Show DOMINO Reference");
    ui->pushButton_exit_DOMINO->hide();

    // Perl path
    ui->Perl_comboBox->setToolTip("< img src = ':/images/GUI_icons/search_lense.png'> Get the default path for perl");
    ui->radioButton_default_perl_path->setToolTip("< img src = ':/images/GUI_icons/search_lense.png'> Use a default path for perl");
    ui->radioButton_custom_perl_path->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Specify a path for perl");
    ui->get_user_perl_path->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Use the Perl path specified in the text box");
    ui->textEdit_perl_path->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Specify a path for perl");

    // Folder
    ui->frame_3->setToolTip("< img src = ':/images/GUI_icons/folder.png'> DOMINO main folder for output results");
    ui->browse_project_dir_initial_tab->setToolTip("< img src = ':/images/GUI_icons/folder.png'> Browse DOMINO main folder for output results");

    // Start or Reset
    ui->continue_button_tab_about->setToolTip("< img src = ':/images/GUI_icons/start.png'> Start DOMINO");
    ui->continue_button_tab_about->setEnabled(false);
    ui->pushButton_Resume_project->setToolTip("< img src = ':/images/GUI_icons/forward.png'> Resume a previous project or start a new one entering to any tab desired");
    ui->pushButton_Resume_project->setEnabled(false);

    // Disable the different tabs of the application and set tips
    ui->tabWidget->setTabToolTip(0,"< img src = ':/images/GUI_icons/home.png'> Introductory tab");

    ui->tabWidget->setTabEnabled(1, false); // Cleaning tab
    ui->tabWidget->setTabToolTip(1,"Input your data and preprocess it");
    ui->tabWidget->setTabEnabled(2, false); // Assembling tab
    ui->tabWidget->setTabToolTip(2,"Assemble your data");
    ui->tabWidget->setTabEnabled(3, false); // Mapping tab
    ui->tabWidget->setTabToolTip(3,"Map or Retrieve alignment files");
    ui->tabWidget->setTabEnabled(4, false); // Molecular marker tab
    ui->tabWidget->setTabToolTip(4,"Discovery/Selection Molecular markers");
    ui->tabWidget->setTabEnabled(5, false); // File-Results Viewer tab
    ui->tabWidget->setTabToolTip(5,"Check files and folders generated");

    // Set processors
    ui->spinBox_num_proc->setValue(idealthreads);
    ui->spinBox_num_proc->setMaximum(idealthreads);
    ui->spinBox_num_proc->setMinimum(1);
    ui->frame_9->setToolTip("< img src = ':/images/GUI_icons/options.png'> Computational parameters to use in some parts of the pipeline");

    //Reference to several objects
    //perl_dialog_log = new MyDialog(this);
    //perl_dialog_log->setFixedSize(1206, 628);
    change_value = new change_progress_bar(this);
    STD_err_dialog = new Error_Dialog(this);
    STD_err_dialog->setFixedSize(637,470);
}

MainWindow_AMM::~MainWindow_AMM() {
    // ###############
    // ## Destroyer ##
    // ###############
    delete ui;
}

void MainWindow_AMM::on_actionDebug_Mode_triggered() { debug_mode = true; }

void MainWindow_AMM::show_exit_button() {
    ui->pushButton_exit_DOMINO->show();
    ui->pushButton_exit_DOMINO->setEnabled(true);
}

void MainWindow_AMM::on_pushButton_exit_DOMINO_clicked() { close(); }


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Setting different tabs @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::set_cleaning_tab() {
    // ########################################
    // ## Set some icons in the cleaning tab ##
    // ########################################

    // Enable buttons on Cleaning tab
    ui->tabWidget->setTabEnabled(1, true);
    ui->browse_Input_file->setFocus();

    // Disable buttons until user provides several options
    ui->pushButton_stop_cleaning->setEnabled(false);
    ui->run_Pipe_cleaning_button->setEnabled(false);
    ui->browse_barcode_file->setEnabled(false);
    ui->checkBox_barcode_user_file->setEnabled(false);
    ui->textEdit_user_barcode_file->setEnabled(false);
    ui->textBrowser_taxa_names->setEnabled(false);
    ui->label_8->setEnabled(false);

    // Set type of file combo box
    if (!type_file_cleaning) {
        //type_files << "Select a type of file" << "454 Roche SFF file (MID tags)" << "FASTQ file 454 Roche (MID tags)" << "Multiple FASTQ 454 Roche" << "Illumina Single End FASTQ file (MID tags)" << "Multiple Single End FASTQ files" << "Illumina Paired-end FASTQ file (MID tags)" << "Multiple Illumina Paired-end FASTQ files";
        type_files << "Select a type of file" << "Roche 454 SFF" << "Roche 454 FASTQ" << "Multiple Roche 454 FASTQ" << "Illumina single-end FASTQ" << "Multiple Illumina single-end FASTQ" << "Illumina paired-end" << "Multiple Illumina paired-end FASTQ";
        ui->comboBox_type_file_cleaning->addItem(type_files.at(0));
        for (int i = 1; i < type_files.size(); i++) {
            ui->comboBox_type_file_cleaning->addItem(QString::number(i) + ": " + type_files.at(i));
        }
        type_file_cleaning = true;
    }
    ui->comboBox_type_file_cleaning->setCurrentIndex(0);

    // Setting minimun, maximun and default values for the spin Boxes in the cleaning_step
    ui->spinBox_cutoff_PHRED->setMinimum(18);
    ui->spinBox_cutoff_PHRED->setValue(20);
    ui->spinBox_cutoff_PHRED->setMaximum(30);
    ui->spinBox_minimun_read_length->setMaximum(120);
    ui->spinBox_minimun_read_length->setValue(100);
    ui->spinBox_minimun_read_length->setMinimum(80);
    ui->spinBox_min_cutoff_length->setMinimum(65);
    ui->spinBox_min_cutoff_length->setValue(70);
    ui->spinBox_min_cutoff_length->setMaximum(80);
    ui->spinBox_threshold_DUST->setValue(7);
    ui->spinBox_threshold_DUST->setMinimum(7);
    ui->spinBox_threshold_DUST->setMaximum(15);

    //Set button icons and tips
    ui->run_Pipe_cleaning_button->setToolTip("< img src = ':/images/GUI_icons/start.png'> Run the computation");
    ui->browse_Input_file->setToolTip("< img src = ':/images/GUI_icons/document_file.png'> Browse input(s) file(s)...");
    ui->browse_Input_file->setEnabled(false);
    ui->browse_barcode_file->setToolTip("< img src = ':/images/GUI_icons/document_file.png'> Get the Barcodes file");
    ui->browse_databases->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Get the additional databases");
    ui->textEdit_user_databases_select->setToolTip("Databases selected");
    ui->Help_button_cleaning->setToolTip("< img src = ':/images/GUI_icons/help.png'> Show description of the cleaning process");
    ui->pushButton_stop_cleaning->setToolTip("< img src = ':/images/GUI_icons/close_delete'> Stop the analysis");
    ui->frame_8->setToolTip("< img src = ':/images/GUI_icons/options.png'> Cleaning-Trimming Parameters for Read Preprocessing");
    ui->frame_7->setToolTip("< img src = ':/images/GUI_icons/options.png'> Parameters for the BLAST search of contaminants");
    ui->frame_5->setToolTip("< img src = ':/images/GUI_icons/search_lense.png'> Input your data and type of databases");
    ui->comboBox_type_file_cleaning->setToolTip("< img src = ':/images/GUI_icons/options.png'> Select your type of file");
    ui->checkBox_no_BLAST_search->setToolTip("Choose wether to perform a BLAST search for putative contaminants or not");
    ui->checkBox_only_additional_database->setToolTip("Use only your own databases provided");
    ui->checkBox_Preprocess_NGS_files->setToolTip("Select wether to clean and trim reads or not");
    ui->textBrowser_taxa_names->setToolTip("Names provided to tag each taxa and reads");

    // Set default
    ui->checkBox_no_BLAST_search->setChecked(true);
    ui->checkBox_Preprocess_NGS_files->setChecked(true);
    ui->browse_databases->setEnabled(false);
}

void MainWindow_AMM::set_assembling_tab() {
    // ######################################################
    // ## Set some icons in the assembling and mapping tab ##
    // ######################################################

    ui->tabWidget->setTabEnabled(2, true); // Assembling tab
    ui->run_Assembling->setEnabled(false);
    ui->pushButton_stop_assembling->setEnabled(false);
    ui->textEdit_input_file_assembling->setEnabled(false);
    ui->browse_Input_file_assembling->setEnabled(false);
    ui->frame_11->setToolTip("Input data parameters for DOMINO assembly");

    if (!type_file_assembling) {
        ui->comboBox_type_file_assembling->addItem(type_files.at(0));
        for (int i = 1; i < type_files.size(); i++) {
            if (i == 3 || i == 5 || i == 7) {
                ui->comboBox_type_file_assembling->addItem(QString::number(i) + ": " + type_files.at(i));
        }}
        type_file_assembling = true;
        ui->comboBox_type_file_assembling->setCurrentIndex(0);
    }

    if (!DOMINO_user_files) {
        ui->comboBox_DOMINO_user_files->addItem("Select type of Data");  // index 0
        ui->comboBox_DOMINO_user_files->addItem("DOMINO pre-processed reads");  // index 1
        ui->comboBox_DOMINO_user_files->addItem("User files");          // index 2
        ui->comboBox_DOMINO_user_files->setCurrentIndex(0);
        DOMINO_user_files = true;
    }

    //Set button icons and tips
    ui->run_Assembling->setToolTip("< img src = ':/images/GUI_icons/start.png'> Run the computation");
    ui->Help_button_Assembling_tab->setToolTip("< img src = ':/images/GUI_icons/help.png'> Show description of the Assembling and Mapping process");
    ui->pushButton_stop_assembling->setToolTip("< img src = ':/images/GUI_icons/close_delete.png'> Stop the analysis");
    ui->textEdit_input_file_assembling->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Get new FASTQ files not used during the DOMINO cleaning step");
    ui->browse_Input_file_assembling->setToolTip("< img src = ':/images/GUI_icons/plus.png'> Get new FASTQ files not used during the DOMINO cleaning step");
    ui->comboBox_DOMINO_user_files->setToolTip("Choose wether to use DOMINO clean files or new FASTQ files");
    ui->frame_13->setToolTip("< img src = ':/images/GUI_icons/options.png'> Parameters for the assembly of reads using MIRA and CAP3");
    ui->comboBox_type_file_assembling->setToolTip("< img src = ':/images/GUI_icons/options.png'> Select your type of file");
}

void MainWindow_AMM::set_mapping_tab() {
    // #####################################################
    // ## Set some icons in the Marker Identification tab ##
    // #####################################################

    ui->tabWidget->setTabEnabled(3, true); // Mapping Tab
    ui->frame_15->setEnabled(true);
    ui->textEdit_input_file_DOMINO_marker->setEnabled(false);
    ui->textEdit_input_file_DOMINO_marker_2->setEnabled(false);

    // Read gap open/extension penalty
    ui->spinBox_pen_read_gap_open->setValue(5);
    ui->spinBox_pen_read_gap_open->setMaximum(10);
    ui->spinBox_pen_read_gap_open->setMinimum(1);
    ui->spinBox_pen_read_gap_ext->setValue(3);
    ui->spinBox_pen_read_gap_ext->setMaximum(10);
    ui->spinBox_pen_read_gap_ext->setMinimum(1);

    // Reference gap open/extension penalty
    ui->spinBox_pen_ref_gap_open->setValue(5);
    ui->spinBox_pen_ref_gap_open->setMaximum(10);
    ui->spinBox_pen_ref_gap_open->setMinimum(1);
    ui->spinBox_pen_ref_gap_ext->setValue(3);
    ui->spinBox_pen_ref_gap_ext->setMaximum(10);
    ui->spinBox_pen_ref_gap_ext->setMinimum(1);

    // Mismatch spinBox
    ui->spinBox_mismatch_penalty->setValue(4);
    ui->spinBox_mismatch_penalty->setMaximum(10);
    ui->spinBox_mismatch_penalty->setMinimum(1);

     // Behaviour
    ui->frame_17->setToolTip("Choose wether to discover new markers or select informative markers");
    ui->radioButton_identify->setEnabled(false);
    ui->radioButton_select->setEnabled(false);

    // Files to use
    if (!type_file_mappingMarker) {
        // option
        ui->comboBox_DOMINO_user_files_markers->addItem("Please Choose a type of file"); // 0
        ui->comboBox_DOMINO_user_files_markers->addItem("DOMINO Contigs"); // DOMINO files 1
        ui->comboBox_DOMINO_user_files_markers->addItem("Multiple taxa references"); // user_assembling_contigs 2
        ui->comboBox_DOMINO_user_files_markers->addItem("Single taxon reference(s)"); // genome 3
        ui->comboBox_DOMINO_user_files_markers->addItem("MSA files"); // msa_alignment 4
        ui->comboBox_DOMINO_user_files_markers->addItem("RAD-MSA files"); //  MSA_RADseq

        // type_input
        ui->comboBox_SINGLE_PAIR_end->addItem("Paired End");  //pair_end, index 0
        ui->comboBox_SINGLE_PAIR_end->addItem("Single End");//single_end, index 1
        type_file_mappingMarker = true;
    }

    //Set button icons and tips
    ui->comboBox_DOMINO_user_files_markers->setToolTip("Please choose wether to use DOMINO files, new contigs assemble and reads or a reference genome, a multiple alignment file or RADseq data");
    ui->comboBox_SINGLE_PAIR_end->setToolTip("Choose if Single or Paired-end reads provided");
    ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get new contig fasta files or genome fasta file");
    ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get new read FASTQ files");
    ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Multiple sequence alignment file or folder");
    ui->frame_14->setToolTip("Mapping and molecular markers file parameters");
    ui->frame_16->setToolTip("< img src = ':/images/GUI_icons/options.png'> Mapping Parameters");
    ui->Help_button_DOMINO_marker_Identif_tab_2->setToolTip("< img src = ':/images/GUI_icons/help.png'> Show description for the Mapping and molecular marker scan process");
    ui->pushButton_mapping->setToolTip("< img src = ':/images/GUI_icons/forward.png'> Proceed to the marker scan tab");
}

void MainWindow_AMM::set_MolecularMarker_identification_tab() {

    ui->tabWidget->setTabEnabled(4, true); // Molecular Marker Identification Tab
    ui->run_DOMINO_marker_identification_button->setEnabled(false);
    ui->pushButton_stop_DOMINO_marker_identification->setEnabled(false);

    // Setting minimun, maximun and default values for the spin Boxes in the ANMs tab
    // CRS: Conserved Region Size
    ui->spinBox_CRS->setValue(40);
    ui->spinBox_CRS->setMaximum(100);
    ui->spinBox_CRS->setMinimum(1);

    // mct
    ui->spinBox_msa_markers->setValue(0);
    ui->spinBox_msa_markers->setToolTip("Minimum number of taxa required for searching");
    ui->label_34->setToolTip("Minimum number of taxa required for searching");

    // Polymorphism
    ui->checkBox_polymorphism->setChecked(false);
    ui->checkBox_polymorphism->setToolTip("DOMINO will use polymorphisms as variable sites");

    // sam/bam checkbox
    ui->keepbam_checkbox->setToolTip("Keep SAM/BAM files generated for the visualization of results");

    // CD: Maximun variations conserved region
    ui->spinBox_MVCR->setValue(1);
    ui->spinBox_MVCR->setMinimum(0);
    ui->spinBox_MVCR->setMaximum(5);

    ui->run_DOMINO_marker_identification_button->setToolTip("< img src = ':/images/GUI_icons/start.png'> Run the computation");
    ui->pushButton_stop_DOMINO_marker_identification->setToolTip("< img src = ':/images/GUI_icons/close_delete'> Stop the analysis");
    ui->Help_button_DOMINO_marker_Identif_tab->setToolTip("< img src = ':/images/GUI_icons/help.png'> Show description for the Mapping and molecular marker scan process");
    ui->frame_18->setToolTip("< img src = ':/images/GUI_icons/options.png'> Select the taxa names for molecular marker scan");
    ui->frame_20->setToolTip("< img src = ':/images/GUI_icons/options.png'> Molecular marker Scan Parameters");

    if (!mapping_checked) { MainWindow_AMM::check_user_options("mapping"); }
}

void MainWindow_AMM::set_treeView_outputFolder() {
    // #########################################
    // ## Set some icons in the Tree View tab ##
    // #########################################

    ui->tabWidget->setTabEnabled(5, true);
    model_dir = new QDirModel(this);
    model_dir->setReadOnly(false);
    model_dir->setSorting(QDir::DirsFirst | QDir::IgnoreCase | QDir::Name);

    QModelIndex outfolder_index = model_dir->index(outfolder);
    ui->treeView->setModel(model_dir);
    ui->treeView->expand(outfolder_index);
    ui->treeView->scrollTo(outfolder_index);
    ui->treeView->setCurrentIndex(outfolder_index);
    ui->treeView->resizeColumnToContents(0);
    ui->pushButton_openFile->setToolTip("< img src = ':/images/GUI_icons/search_lense.png'> Open a plain text file in a Text Browser window");
    ui->pushButton_refresh_Dir->setToolTip("< img src = ':/images/GUI_icons/refresh.png'> Refresh the information");
}

void MainWindow_AMM::setLabelData() {

    int index_combo_box = ui->comboBox_DOMINO_user_files_markers->currentIndex();
    ui->label_data->clear();

    QString item;
    if (index_combo_box == 1) {
        item = QString("DOMINO Contigs");
    } else if (index_combo_box == 2) {
        item = QString("Multiple taxa references");
    } else if (index_combo_box == 3) {
        item = QString("Single taxon reference(s)");
    } else if (index_combo_box == 4) {
        item = QString("MSA files");
    } else if (index_combo_box == 5) {
        item = QString("RAD-MSA files");
    }

    QString develop_analysis("Development Module: ");
    if (ui->radioButton_identify->isChecked()) {
        develop_analysis.append(" Discovery");
    } else if (ui->radioButton_select->isChecked()) {
        develop_analysis.append(" Selection");
    }
    develop_analysis.append(" of informative markers\n");
    QString typeData("\nType of Data: "); typeData.append(item).append("\n");
    typeData.append(develop_analysis);

    //cout << typeData.toStdString() << endl;
    ui->label_data->setText(typeData);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Push Buttons in DOMINO Menu Bar @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_actionReference_triggered() {
    // ##################################################################
    // ## Open a dialog showing the reference of DOMINO paper and link ##
    // ##################################################################
    reference_dialog = new Error_Dialog(this);
    reference_dialog->setFixedSize(589,470);
    reference_dialog->setWindowTitle("Reference");
    reference_dialog->show_Reference();
    reference_dialog->show();
}

void MainWindow_AMM::on_actionManual_triggered() {
    // ########################################################
    // ## Open a web browser showing a PDF manual for DOMINO ##
    // ########################################################
    QString website_link = "http://www.ub.edu/softevol/domino/";
    QDesktopServices::openUrl(QUrl(website_link ));
}

void MainWindow_AMM::on_actionCite_DOMINO_triggered() {
    // ###############################################
    // ## Open a web browser showing a DOMINO paper ##
    // ###############################################
    QString website_link = "http://bioinformatics.oxfordjournals.org/content/early/2016/09/02/bioinformatics.btw534";
    QDesktopServices::openUrl(QUrl(website_link ));
}


void MainWindow_AMM::on_actionAbout_DOMINO_triggered(){
    about_dialog = new MyDialog(this);
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("About DOMINO");
    about_dialog->aboutDOMINO();
    about_dialog->show();
}

void MainWindow_AMM::on_actionLicense_triggered()
{
    about_dialog = new MyDialog(this);
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("DOMINO License");
    about_dialog->DOMINOlicense();
    about_dialog->show();
}


void MainWindow_AMM::on_actionShow_Workflow_triggered() { MainWindow_AMM::on_show_workflow_button_clicked(); }
void MainWindow_AMM::on_actionClose_triggered() { MainWindow_AMM::close(); }

void MainWindow_AMM::on_actionReset_GUI_triggered() {
    // Reset GUI
    qApp->quit();
    QProcess::startDetached(qApp->arguments()[0], qApp->arguments());
}

void MainWindow_AMM::on_actionAbout_Qt_triggered() {
    // Show info for Qt
    qApp->aboutQt();
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Push Buttons in Tree View Tab @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_pushButton_refresh_Dir_clicked() {
    // #################################################################
    // ## Refresh and uptades the tree view of the directory selected ##
    // #################################################################
    MainWindow_AMM::set_treeView_outputFolder();
}

void MainWindow_AMM::on_pushButton_openFile_clicked() {

    // ##########################################################
    // ## Open the file chosen in a text browser. Check if     ##
    // ## it is readable, a plain text file and not too large  ##
    // ##########################################################
    textFileBrowser = new MyDialog(this);
    textFileBrowser->setFixedSize(1206, 628);
    QModelIndex new_index = ui->treeView->currentIndex(); // Get the index of the file or dir
    if(!new_index.isValid()) return; // Check if it has valid index

    // Open if it is a file
    if(model_dir->fileInfo(new_index).isFile()) {
        // Get file path and file name
        QString File_name(model_dir->fileName(new_index));
        textFileBrowser->textFileBrowser_dialog_show(File_name);
        QString File_path(model_dir->filePath(new_index));
        QFile file(File_path);
        // Check if it is readable
        if (!file.open(QIODevice::ReadOnly)) {
            //cout << "File is not openned" << endl;
            STD_err_dialog->showMessage("Not able to open this file. It is not readable");
            STD_err_dialog->setWindowTitle("Error ocurred while opening the file");
            STD_err_dialog->show();
         } else {
            qint64 size_file(model_dir->fileInfo(new_index).size()); // Get file size
            if (size_file > 10000000) {
                // Show Error dialog if file is too big
                STD_err_dialog->showMessage("Not able to open this file because it is too large");
                STD_err_dialog->setWindowTitle("Error ocurred while opening the file");
                STD_err_dialog->show();
            } else {
                // Get file information into a string
                QTextStream file_info(&file);
                QString temp = file_info.readAll();

                // Check and only display plain text files
                if (temp.isSimpleText()) {
                    textFileBrowser->show();
                    textFileBrowser->textFileBrowser_dialog(temp);
                } else {
                    // Show Error dialog if file chosen is not plain text file
                    STD_err_dialog->showMessage("Not able to open this file as it is not a plain text file");
                    STD_err_dialog->setWindowTitle("Error ocurred while opening the file");
                    STD_err_dialog->show();
        }}}
    } else {
        // Show Error dialog if file chosen is a directory
        STD_err_dialog->showMessage("Not able to open this file as it is not a plain text file but a directory");
        STD_err_dialog->setWindowTitle("Error ocurred while opening the file");
        STD_err_dialog->show();
    }
}


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Introductory tab buttons @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_radioButton_default_perl_path_clicked() {
    ui->textEdit_perl_path->setEnabled(false);
    ui->get_user_perl_path->setEnabled(false);
    ui->Perl_comboBox->setEnabled(true);
    ui->continue_button_tab_about->setEnabled(true);
    ui->pushButton_Resume_project->setEnabled(true);

    // If several perl paths, split them and put them into combo_box
    QStringList perl_paths = perl_path.split("\n");
    for (int i = 0; i < perl_paths.size(); i++ ) {
        if (perl_paths.at(i).contains("gz")) continue;
        if (perl_paths.at(i).contains("share")) continue;
        if (perl_paths.at(i).contains("lib")) continue;
        if (perl_paths.at(i).contains(":")) continue;
        ui->Perl_comboBox->addItem(perl_paths.at(i));
    }
    ui->browse_project_dir_initial_tab->setFocus();
}

void MainWindow_AMM::on_radioButton_custom_perl_path_clicked() {
    ui->textEdit_perl_path->setEnabled(true);
    ui->get_user_perl_path->setEnabled(true);
    ui->Perl_comboBox->setEnabled(false);
}

void MainWindow_AMM::on_pushButton_Resume_project_clicked() {
    // ###################################################################################################
    // ## Pressing this check Box, it would be possible to enter any step and resume a previous project ##
    // ###################################################################################################
    if (outfolder.isEmpty()) {
        STD_err_dialog->showMessage("DOMINO output folder missing. Please provide it here before continuing with DOMINO");
        STD_err_dialog->setWindowTitle("Output folder missing");
        STD_err_dialog->show();
    } else {
        ui->tabWidget->setTabEnabled(1, true); // Enable and set Cleaning tab
        set_cleaning_tab();
        ui->tabWidget->setTabEnabled(2, true); // Enable and set Assembling & Mapping tab
        set_assembling_tab();
        ui->tabWidget->setTabEnabled(3, true); // Enable and set mapping tab
        set_mapping_tab();
        ui->tabWidget->setTabEnabled(4, false); // Disable Molecular Marker identification tab
        domino_short_cut = true;
    }
}

void MainWindow_AMM::on_get_user_perl_path_clicked() {
    // ######################################
    // ## Get the perl path user specified ##
    // ######################################
    perl_abs_path = ui->textEdit_perl_path->toPlainText();
    if (perl_abs_path == "") { perl_abs_path = perl_path; }
    ui->continue_button_tab_about->setEnabled(true);
    ui->pushButton_Resume_project->setEnabled(true);
    //cout << perl_path.toStdString() << endl;
}

void MainWindow_AMM::on_browse_project_dir_initial_tab_clicked() {
    // ##################################################################
    // ## Let the user browse project folder opening a window at $HOME ##
    // ##################################################################
    outfolder = QFileDialog::getExistingDirectory(this, tr("Project Folder for Output Results"), QDir::homePath());
    ui->textEdit_project_dir_initial_tab->setText(outfolder);
    MainWindow_AMM::set_treeView_outputFolder();
    ui->continue_button_tab_about->setFocus();
}

void MainWindow_AMM::on_continue_button_tab_about_clicked() {
    // ###########################################
    // ## Set the cleaning tab and start DOMINO ##
    // ###########################################
    if (outfolder.isEmpty()) {
        STD_err_dialog->showMessage("DOMINO output folder missing. Please provide it here before continuing with DOMINO");
        STD_err_dialog->setWindowTitle("Output folder missing");
        STD_err_dialog->show();
    } else {
        set_cleaning_tab();
        ui->tabWidget->setCurrentIndex(1);
        MainWindow_AMM::set_treeView_outputFolder();
    }
}

void MainWindow_AMM::on_show_workflow_button_clicked() {
    // ###########################################################
    // ## Shows a dialog with general information of DOMINO and ##
    // ## lets the user open workflow images of the three steps ##
    // ###########################################################
    workflow_dialog = new Dialog_workflow(this);
    workflow_dialog->setFixedSize(870,800);
    workflow_dialog->setToolTip("< img src = ':/images/GUI_icons/Help1.png'> Show DOMINO Workflow");
    workflow_dialog->setWindowTitle("General DOMINO process information");
    workflow_dialog->show();
}

void MainWindow_AMM::on_pushButton_showRef_intro_tab_clicked() { MainWindow_AMM::on_actionReference_triggered(); }
void MainWindow_AMM::on_pushButton_showManual_intro_tab_clicked() { MainWindow_AMM::on_actionManual_triggered(); }


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Input data and Cleaning Tab @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_Help_button_cleaning_clicked() {
    // ####################################################################
    // ## Sets Help button displaying information of the cleaning process ##
    // ####################################################################

    about_dialog = new MyDialog(this);
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("Information about the Input data and Cleaning step");
    about_dialog->about_cleaning_pipeline();
    about_dialog->show();
}

void MainWindow_AMM::on_browse_Input_file_clicked() {
    // ########################################################################
    // ## Let the user browse input file(s), opening a window at $HOME ##
    // ########################################################################
    QString example_path_tmp = APPDIR;
    QStringList example_path_list = example_path_tmp.split("/");
    QString example_path;
    int ind = example_path_list.size();
    int new_ind = ind-2;
    for (int i=0; i <= new_ind; i++) {
        example_path.append(example_path_list.at(i)).append("/");
    }
    example_path.append("example"); // example folder in DOMINO exec path

    QStringList input_cleaning_files = QFileDialog::getOpenFileNames(this,tr("Find file(s)"), example_path);
    QString input_cleaning_files_String = input_cleaning_files.join("\n");
    ui->textEdit_input_file->setText(input_cleaning_files_String);
    input_file_name_cleaning = input_cleaning_files_String;
    if (!cleaning_checked) { MainWindow_AMM::check_user_options("clean"); }

    int index = ui->comboBox_type_file_cleaning->currentIndex();
    if ((index == 3) | (index == 5) | (index == 7) ){
        MainWindow_AMM::check_species_names(input_cleaning_files, "cleaning");
    }
}

void MainWindow_AMM::on_browse_barcode_file_clicked() {
    // ###################################################################
    // ## Let the user browse a barcode file, opening a window at $HOME ##
    // ###################################################################
    barcode_user_file= QFileDialog::getOpenFileName(this,tr("Find Barcode file"), QDir::homePath());
    ui->checkBox_barcode_user_file->setChecked(true);
    ui->textEdit_user_barcode_file->setText(barcode_user_file);

    MainWindow_AMM::check_user_options("clean");
}

void MainWindow_AMM::on_browse_databases_clicked() {
    // ########################################################################
    // ## Let the user browse a several databases, opening a window at $HOME ##
    // ########################################################################
    db_list = QFileDialog::getOpenFileNames(this,tr("Find databases"), QDir::homePath(), tr("Fasta format(*.fasta)"));
    QString db_string = db_list.join("\n");
    ui->textEdit_user_databases_select->setText(db_string);
    //ui->checkBox_additional_database->setChecked(true);
}

void MainWindow_AMM::on_checkBox_no_BLAST_search_stateChanged(int flag) {

    if (flag == 2) {
        // checked
        ui->browse_databases->setEnabled(true);
        ui->checkBox_only_additional_database->setEnabled(true);
        ui->textEdit_user_databases_select->setEnabled(true);
    } else {
        // Unchecked
        ui->browse_databases->setEnabled(false);
        ui->checkBox_only_additional_database->setEnabled(false);
        ui->textEdit_user_databases_select->setEnabled(false);
    }
}

void MainWindow_AMM::on_checkBox_only_additional_database_clicked() { ui->browse_databases->setFocus(); }

void MainWindow_AMM::on_comboBox_type_file_cleaning_currentIndexChanged(int index) {

    if ((index == 1) | (index == 2) | (index == 4) | (index == 6)) {
        ui->browse_barcode_file->setEnabled(true);
        ui->browse_barcode_file->setFocus();
        ui->checkBox_barcode_user_file->setEnabled(true);
        ui->checkBox_barcode_user_file->setChecked(true);
        ui->textEdit_user_barcode_file->setEnabled(true);
        ui->textBrowser_taxa_names->setEnabled(false);
        ui->label_8->setEnabled(false);
        ui->browse_Input_file->setEnabled(true);
    } else if (index == 0) {
        ui->checkBox_barcode_user_file->setChecked(false);
        ui->checkBox_barcode_user_file->setEnabled(false);
        ui->checkBox_barcode_user_file->setChecked(false);
        ui->textBrowser_taxa_names->setEnabled(false);
        ui->label_8->setEnabled(false);
        ui->textEdit_user_barcode_file->setEnabled(false);
        ui->browse_Input_file->setEnabled(false);
        ui->browse_barcode_file->setEnabled(false);
    } else {
        ui->checkBox_barcode_user_file->setChecked(false);
        ui->checkBox_barcode_user_file->setEnabled(false);
        ui->checkBox_barcode_user_file->setChecked(false);
        ui->textEdit_user_barcode_file->setEnabled(false);
        ui->browse_barcode_file->setEnabled(false);
        ui->textBrowser_taxa_names->setEnabled(true);
        ui->label_8->setEnabled(true);
        ui->browse_Input_file->setEnabled(true);
    }
}

void MainWindow_AMM::on_checkBox_only_additional_database_stateChanged(int arg1) {
    if (arg1 == 2) {
        ui->browse_databases->setEnabled(true);
    } else {
        ui->browse_databases->setEnabled(false);
    }
}

void MainWindow_AMM::on_run_Pipe_cleaning_button_clicked(){
    // ##########################################################################################
    // ##  When user press Run Button, get some options and call the cleaning pipeline script  ##
    // ##########################################################################################

    bool run_process;
    // Control if user specified a Perl path if not use, the default given in the comboBox
    if (perl_abs_path == "") {
        perl_abs_path = ui->Perl_comboBox->currentText();
    }

    // Get the options the user input
    PHRED_sc = MainWindow_AMM::ui->spinBox_cutoff_PHRED->value(); //PHRED score
    length_cutoff = MainWindow_AMM::ui->spinBox_min_cutoff_length->value(); //cutoff length % satisfying PHRED score
    min_length = MainWindow_AMM::ui->spinBox_minimun_read_length->value(); //Minimum read length
    thresh_dust = MainWindow_AMM::ui->spinBox_threshold_DUST->value(); //threshold dust algorithm
    num_proc = MainWindow_AMM::ui->spinBox_num_proc->value(); //number CPU

    // Get type of file
    QString identification_string;
    int index = ui->comboBox_type_file_cleaning->currentIndex();
    if ((index == 1) | (index == 2) | (index == 4) | (index == 6)) {
        //check barcode file is provided
        run_process = true;
    } else {
        run_process = true;
    }

    // Check BLAST search radio button
    if (ui->checkBox_no_BLAST_search->isChecked()) {
        //
    } else {
        // Sets the option to avoid contaminant search ##
        QString no_BLAST("-no_db_search ");
        additional_perl_options_cleaning.append(no_BLAST);
    }

    // Check ONLY aditional databases radio button
    if (ui->checkBox_only_additional_database->isChecked()) {
        // Sets the option for BLAST search against only user databases provided ##
        QString only_db_user("-only_user_db ");
        additional_perl_options_cleaning.append(only_db_user);
    }

    // Check if only tag files checked
    if (ui->checkBox_Preprocess_NGS_files->isChecked()) {
        //Default behaviour
    } else {
        QString only_tagging("-only_tag_files ");
        additional_perl_options_cleaning.append(only_tagging);
    }

    // Check DO NOT delete tmp files radio button
    if (ui->checkBox_tmp_files->isChecked()) {
        // Sets the option to avoid deleting temporary files ##
        QString deleting_tmp(" -TempFiles ");
        additional_perl_options_cleaning.append(deleting_tmp);
    }

    // Check the barcode radio button
    if (ui->checkBox_barcode_user_file->isChecked()) {
        // Append the barcode file provided to the additional options when calling DOMINO cleaning pipeline
        if (barcode_user_file.isEmpty()) {
            STD_err_dialog->showMessage("Checkbox for Barcode file has been checked but no file has been selected.\nMake sure you choose a file or unchecked this option and use default Barcodes");
            STD_err_dialog->setWindowTitle("Barcode File ERROR");
            STD_err_dialog->show();
            run_process = false;
        } else {
            QString user_barcode_option("-b ");
            user_barcode_option.append(barcode_user_file);
            user_barcode_option.append(" ");
            additional_perl_options_cleaning.append(user_barcode_option);
        }
    }

    // Check aditional databases radio button
    if (ui->checkBox_only_additional_database->isChecked()) {
        if (db_list.isEmpty()) {
            STD_err_dialog->showMessage("Checkbox for Additional databases has been checked but no file has been selected.\nMake sure you choose a file or unchecked this option and use Default Contaminant Databases");
            STD_err_dialog->setWindowTitle("Additional Databases ERROR");
            STD_err_dialog->show();
            run_process = false;
        } else {
            QString db_option;
            for (int j=0; j < db_list.size(); ++j) {
                db_option.append("-db ");
                db_option.append(db_list.at(j));
                db_option.append(" ");
            }
            additional_perl_options_cleaning.append(db_option);
        }
    }

    //Run the pipeline
    if (run_process == true) {
        ui->run_Pipe_cleaning_button->setEnabled(false); // Avoid user can click twice and start another computation
        ui->pushButton_stop_cleaning->setEnabled(true); // Allow user to stop the computation
        MainWindow_AMM::cleaning_pipe(identification_string);
        //MainWindow_AMM::show_perl_dialog();
        MainWindow_AMM::start_status_bar();
    }
}

void MainWindow_AMM::on_pushButton_stop_cleaning_clicked() {
    // ###############################################
    // ## Stops and terminates the Cleaning Process ##
    // ###############################################

    cleaning_process_killed = true;
    proc_Perl_cleaning_pipe->kill();

    // Defatult tab parameters
    MainWindow_AMM::set_cleaning_tab();
    ui->textEdit_input_file->setText("");
    ui->textEdit_user_barcode_file->setText("");
    ui->textEdit_user_databases_select->setText("");
    input_file_name_cleaning = "";
}

void MainWindow_AMM::on_checkBox_Preprocess_NGS_files_stateChanged(int flag) {
    if (flag == 2) {
        // Preprocess of NGS files -- checked
        ui->label_4->setEnabled(true);
        ui->spinBox_cutoff_PHRED->setEnabled(true);
        ui->label_5->setEnabled(true);
        ui->label_6->setEnabled(true);
        ui->label_7->setEnabled(true);
        ui->spinBox_minimun_read_length->setEnabled(true);
        ui->spinBox_min_cutoff_length->setEnabled(true);
        ui->spinBox_threshold_DUST->setEnabled(true);
        ui->checkBox_no_BLAST_search->setEnabled(true);
        ui->checkBox_no_BLAST_search->setChecked(true);
    } else {
        // NO Preprocess of NGS files -- Unchecked
        ui->label_4->setEnabled(false);
        ui->spinBox_cutoff_PHRED->setEnabled(false);
        ui->label_5->setEnabled(false);
        ui->label_6->setEnabled(false);
        ui->label_7->setEnabled(false);
        ui->spinBox_minimun_read_length->setEnabled(false);
        ui->spinBox_min_cutoff_length->setEnabled(false);
        ui->spinBox_threshold_DUST->setEnabled(false);
        ui->checkBox_no_BLAST_search->setEnabled(false);
        ui->checkBox_no_BLAST_search->setChecked(false);
    }
    MainWindow_AMM::check_user_options("clean");
}


// @@@@@@@@@@@@@@@@@@@@
// @@ Assembling Tab @@
// @@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_browse_Input_file_assembling_clicked() {
    // ########################################################################
    // ## Let the user browse a SFF or FASTQ file, opening a window at $HOME ##
    // ########################################################################

    QStringList input_assembling_files = QFileDialog::getOpenFileNames(this,tr("Find file(s)"), QDir::homePath());
    QString input_assembling_files_String = input_assembling_files.join("\n");
    ui->textEdit_input_file_assembling->setText(input_assembling_files_String);
    input_file_name_assembling = input_assembling_files_String;

    // Set text in Text-browser
    ui->textEdit_input_file_assembling->setText(input_file_name_assembling);
    if (!assembly_checked) { MainWindow_AMM::check_user_options("assemble"); }
}

void MainWindow_AMM::on_run_Assembling_clicked() {
    // #############################################################
    // ##  When user press Run Button, get some options and call  ##
    // ##  the perl assembling-mapping pipeline script            ##
    // #############################################################

    ui->run_Assembling->setEnabled(false); // Avoid user can click twice and start another computation
    ui->pushButton_stop_assembling->setEnabled(true);
    MainWindow_AMM::run_assembling_pipeline();
}

void MainWindow_AMM::on_pushButton_stop_assembling_clicked() {
    // #########################################################
    // ## Stops and terminates the Assembling-Mapping Process ##
    // #########################################################

    int assembling_state = proc_perl_assembling_pipe->state();
    if (assembling_state == 2) {
        assembling_process_killed = true;
        proc_perl_assembling_pipe->kill();
    }
    // Default configuration for tab
    MainWindow_AMM::set_assembling_tab();
    ui->textEdit_input_file_assembling->setText("");
    input_file_name_assembling = "";
}

void MainWindow_AMM::on_Help_button_Assembling_tab_clicked()
{
    // ###############################################################################
    // ## Sets Help button displaying information of the assembling-mapping process ##
    // ###############################################################################

    about_dialog = new MyDialog(this);
    about_dialog->setToolTip("< img src = ':/images/GUI_icons/Help1.png'> Information about this Assembling and Mapping step");
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("Information about the Assembling step");
    about_dialog->about_assembling_pipeline();
    about_dialog->show();
}

void MainWindow_AMM::on_comboBox_DOMINO_user_files_currentIndexChanged(const QString &arg1) { // Assembly
    //cout << arg1.toStdString() << endl;
    if (arg1 == "User files") {
        ui->textEdit_input_file_assembling->setEnabled(true);
        ui->browse_Input_file_assembling->setEnabled(true);
        ui->run_Assembling->setEnabled(false);
    } else if (arg1 == "DOMINO pre-processed reads") {
        ui->textEdit_input_file_assembling->setEnabled(false);
        ui->browse_Input_file_assembling->setEnabled(false);
        ui->run_Assembling->setEnabled(true);
    } else {
        ui->textEdit_input_file_assembling->setEnabled(false);
        ui->browse_Input_file_assembling->setEnabled(false);
        ui->run_Assembling->setEnabled(false);
}}

void MainWindow_AMM::on_comboBox_type_file_assembling_currentIndexChanged(int index) {
    MainWindow_AMM::check_user_options("assemble"); //cout << QString::number(index).toStdString() << endl;
}


// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ DOMINO mapping Tab @@
// @@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_pushButton_mapping_clicked() {
    ui->tabWidget->setCurrentIndex(4);
    MainWindow_AMM::check_user_options("mapping");
    //MainWindow_AMM::on_pushButton_clicked();
}

void MainWindow_AMM::on_browse_file_DOMINO_marker_clicked() {
    // #########################################################################################
    // ## Let the user browse a file (contigs or reference genome), opening a window at $HOME ##
    // #########################################################################################

    int index_combo_box = ui->comboBox_DOMINO_user_files_markers->currentIndex();
    QString input_marker_files_String;
    if (index_combo_box == 2) { // multiple assembly contigs
        QStringList input_marker_files = QFileDialog::getOpenFileNames(this,tr("Find file(s)"), QDir::homePath());
        input_marker_files_String = input_marker_files.join("\n");
    } else if (index_combo_box == 3) { // genome
        input_marker_files_String = QFileDialog::getOpenFileName(this,tr("Find file"), QDir::homePath());
    }
    ui->textEdit_input_file_DOMINO_marker->setText(input_marker_files_String);
    input_file_name_DOMINO_marker = input_marker_files_String;
    //cout << input_file_name_DOMINO_marker.toStdString() << endl;
}

void MainWindow_AMM::on_browse_file_DOMINO_marker_2_clicked() {
    // #########################################################################################
    // ## Let the user browse a file (clean reads), opening a window at $HOME ##
    // #########################################################################################
    QStringList input_marker_files_2 = QFileDialog::getOpenFileNames(this,tr("Find file(s)"), QDir::homePath());
    QString input_marker_files_String_2 = input_marker_files_2.join("\n");
    ui->textEdit_input_file_DOMINO_marker_2->setText(input_marker_files_String_2);
    input_file_name_DOMINO_marker_2 = input_marker_files_String_2;
    MainWindow_AMM::check_species_names(input_marker_files_2, "marker");
}

void MainWindow_AMM::on_browse_file_DOMINO_marker_3_clicked() {
    // #########################################################################################
    // ## Let the user browse a file (msa file or folder), opening a window at $HOME ##
    // #########################################################################################

    ui->textEdit_input_file_DOMINO_marker_3->setText("");

    if (ui->radioButton_file->isChecked()) {
        // Browse file
        QString msa_file = QFileDialog::getOpenFileName(this,tr("Find file(s)"), QDir::homePath());
        ui->textEdit_input_file_DOMINO_marker_3->setText(msa_file);
    } else {
        // Browse folder
        QString msa_folder = QFileDialog::getExistingDirectory(this, tr("Find Alignment folder"), QDir::homePath());
        ui->textEdit_input_file_DOMINO_marker_3->setText(msa_folder);
    }
}

void MainWindow_AMM::on_comboBox_DOMINO_user_files_markers_currentIndexChanged(const QString &arg1) {

    if (arg1 == "DOMINO Contigs") {
        ui->frame_15->setEnabled(false);
        ui->frame_16->setEnabled(true);
        ui->pushButton_mapping->setEnabled(true);
        ui->comboBox_SINGLE_PAIR_end->setEnabled(true);
        ui->keepbam_checkbox->setEnabled(true);

        // Behaviour
        ui->radioButton_identify->setEnabled(true);
        ui->radioButton_identify->clicked();
        ui->radioButton_identify->setChecked(true);
        ui->radioButton_select->setEnabled(false);
        ui->radioButton_select->setToolTip("Selection Development analysis it is not available yet for this option of input data");

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->radioButton_file->setEnabled(false);
        ui->radioButton_folder->setEnabled(false);

        if (!assembly_finished) {
            // Retrieve clean fastq files from DOMINO_clean_data folder
            QDir dir = outfolder;
            dir.setFilter(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);
            QFileInfoList list_files_info_dir = dir.entryInfoList();
            QString DOMINO_clean_data_folder_string;
            for (int i=0; i < list_files_info_dir.size(); i++) {
                QFileInfo file(list_files_info_dir.at(i));
                QString name(file.fileName());
                if (name.contains("DM_clean_data")) {
                    DOMINO_clean_data_folder_string.append(name);
                }
            }
            DOMINO_clean_data_folder_string.append("/");
            DOMINO_clean_data_folder_string.prepend("/");
            DOMINO_clean_data_folder_string.prepend(outfolder);
            QDir DOMINO_clean_data_folder(DOMINO_clean_data_folder_string);
            DOMINO_clean_data_folder.setFilter(QDir::Files | QDir::NoSymLinks | QDir::NoDotAndDotDot);
            QFileInfoList list_files_info = DOMINO_clean_data_folder.entryInfoList();
            QStringList list_files;
            for (int i=0; i < list_files_info.size(); i++) {
                QFileInfo file(list_files_info.at(i));
                QString name(file.fileName());
                list_files << name;
            }
            MainWindow_AMM::check_species_names(list_files, "mapping");
        }
        //ui->tabWidget->setTabEnabled(4, true); // Enable and set Molecular Marker identification tab
        MainWindow_AMM::set_MolecularMarker_identification_tab();

    } else if (arg1 == "Multiple taxa references") {

        ui->frame_15->setEnabled(true);

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Multiple sequence alignment file or folder");
        ui->comboBox_SINGLE_PAIR_end->setEnabled(true);

        // Behaviour
        ui->radioButton_identify->setEnabled(true);
        ui->radioButton_identify->clicked();
        ui->radioButton_identify->setChecked(true);
        ui->radioButton_select->setEnabled(false);
        ui->radioButton_select->setToolTip("Selection Development analysis it is not available yet for this option of input data");

        // Reference
        ui->label_18->setEnabled(true);
        ui->textEdit_input_file_DOMINO_marker->setEnabled(true);
        ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Contig FASTA files");
        ui->browse_file_DOMINO_marker->setEnabled(true);

        // Clean reads
        ui->textEdit_input_file_DOMINO_marker_2->setEnabled(true);
        ui->label_16->setEnabled(true);
        ui->browse_file_DOMINO_marker_2->setEnabled(true);
        ui->textEdit_input_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");
        ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");

        // Mapping parameters
        ui->frame_16->setEnabled(true);
        ui->pushButton_mapping->setEnabled(true);
        ui->keepbam_checkbox->setEnabled(true);

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->radioButton_file->setEnabled(false);
        ui->radioButton_folder->setEnabled(false);

    } else if (arg1 == "Single taxon reference(s)") {

        ui->frame_15->setEnabled(true);

        // Behaviour
        ui->radioButton_identify->setEnabled(true);
        ui->radioButton_identify->clicked();
        ui->radioButton_identify->setChecked(true);
        ui->radioButton_select->setEnabled(false);
        ui->radioButton_select->setToolTip("Selection Development analysis it is not available yet for this option of input data");

        // MSA boxes
        ui->radioButton_file->setEnabled(false);
        ui->radioButton_folder->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Multiple sequence alignment file or folder");

        // Reference
        ui->label_18->setEnabled(true);
        ui->textEdit_input_file_DOMINO_marker->setEnabled(true);
        ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Genome FASTA files");
        ui->browse_file_DOMINO_marker->setEnabled(true);

        // Clean reads
        ui->textEdit_input_file_DOMINO_marker_2->setEnabled(true);
        ui->label_16->setEnabled(true);
        ui->browse_file_DOMINO_marker_2->setEnabled(true);
        ui->textEdit_input_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");
        ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");

        // Mapping parameters
        ui->frame_16->setEnabled(true);
        ui->pushButton_mapping->setEnabled(true);
        ui->keepbam_checkbox->setEnabled(true);
        ui->comboBox_SINGLE_PAIR_end->setEnabled(true);

    } else if (arg1 == "MSA files") {

        ui->frame_15->setEnabled(true);

        // Behaviour
        ui->radioButton_identify->setEnabled(true);
        ui->radioButton_select->setEnabled(true);
        ui->frame_17->setToolTip("Choose wether to discover new markers or select informative markers");

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(true);
        ui->browse_file_DOMINO_marker_3->setEnabled(true);
        ui->label_17->setEnabled(true);

        ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Multiple sequence alignment file or folder");
        ui->radioButton_file->setEnabled(true);
        ui->radioButton_folder->setEnabled(true);
        ui->radioButton_file->setChecked(true);
        ui->radioButton_folder->setChecked(false);

        // Reference
        ui->label_18->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker->setEnabled(false);
        ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get reference FASTA files");
        ui->browse_file_DOMINO_marker->setEnabled(false);

        // Clean reads
        ui->textEdit_input_file_DOMINO_marker_2->setEnabled(false);
        ui->label_16->setEnabled(false);
        ui->browse_file_DOMINO_marker_2->setEnabled(false);
        ui->comboBox_SINGLE_PAIR_end->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");
        ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");

        // Mapping parameters
        ui->frame_16->setEnabled(false);
        ui->pushButton_mapping->setEnabled(true);
        ui->keepbam_checkbox->setEnabled(false);

    } else if (arg1 == "RAD-MSA files") {

        ui->frame_15->setEnabled(true);

        // Behaviour
        ui->radioButton_identify->setEnabled(false);
        ui->radioButton_identify->setToolTip("Discovery Development analysis it is not available yet for this option of input data");
        ui->radioButton_select->setEnabled(true);
        ui->radioButton_select->clicked();
        ui->radioButton_select->setChecked(true);

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(true);
        ui->browse_file_DOMINO_marker_3->setEnabled(true);
        ui->label_17->setEnabled(true);
        ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get RADseq file: .loci (pyRAD) or .fa (STACKS)");
        ui->radioButton_file->setEnabled(true);
        ui->radioButton_file->setChecked(true);
        ui->radioButton_folder->setEnabled(false);

        // Reference
        ui->label_18->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker->setEnabled(false);
        ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get reference FASTA files");
        ui->browse_file_DOMINO_marker->setEnabled(false);

        // Clean reads
        ui->textEdit_input_file_DOMINO_marker_2->setEnabled(false);
        ui->label_16->setEnabled(false);
        ui->browse_file_DOMINO_marker_2->setEnabled(false);
        ui->comboBox_SINGLE_PAIR_end->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");
        ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");

        // Mapping parameters
        ui->frame_16->setEnabled(false);
        ui->pushButton_mapping->setEnabled(true);
        ui->keepbam_checkbox->setEnabled(false);

    } else {
        ui->frame_15->setEnabled(false);
        ui->keepbam_checkbox->setEnabled(true);

        // MSA boxes
        ui->textEdit_input_file_DOMINO_marker_3->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->browse_file_DOMINO_marker_3->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Multiple sequence alignment file or folder");

        // Reference
        ui->label_18->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker->setEnabled(false);
        ui->browse_file_DOMINO_marker->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get reference FASTA files");
        ui->browse_file_DOMINO_marker->setEnabled(false);

        // Clean reads
        ui->textEdit_input_file_DOMINO_marker_2->setEnabled(false);
        ui->label_16->setEnabled(false);
        ui->comboBox_SINGLE_PAIR_end->setEnabled(true);
        ui->browse_file_DOMINO_marker_2->setEnabled(false);
        ui->textEdit_input_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");
        ui->browse_file_DOMINO_marker_2->setToolTip("< img src = ':/images/GUI_icons/plus_small.png'> Get Clean FASTQ files");

        // Mapping parameters
        ui->frame_16->setEnabled(false);
    }

    MainWindow_AMM::setLabelData();
}

void MainWindow_AMM::on_radioButton_identify_clicked() {

    ui->label_25->setText("Length in bp (VL):");
    ui->label_25->setEnabled(true);
    ui->spinBox_VRS_max->setEnabled(true);

    // VL: Variable Region Size
    ui->spinBox_VRS_min->setValue(400);
    ui->spinBox_VRS_min->setMinimum(1);
    ui->spinBox_VRS_min->setMaximum(10000);
    ui->spinBox_VRS_min->hide();
    ui->label_11->hide();
    ui->label_10->hide();

    // VL: Variable Region Size
    ui->spinBox_VRS_max->setValue(400);
    ui->spinBox_VRS_max->setMinimum(1);
    ui->spinBox_VRS_max->setMaximum(10000);
    ui->spinBox_VRS_max->setEnabled(true);

    ui->frame_21->setEnabled(true);
    ui->label_28->setEnabled(true);
    ui->label_37->setEnabled(true);
    ui->doubleSpinBox_variation_percentage->setEnabled(true);
    ui->spinBox_CRS->setValue(20);
    ui->spinBox_MVCR->setValue(1);
}

void MainWindow_AMM::on_radioButton_select_clicked() {

    int index = ui->comboBox_DOMINO_user_files_markers->currentIndex();

    // Conserved
    ui->frame_21->setEnabled(false);
    ui->label_37->setEnabled(false);

    if (index == 4) {

        ui->label_25->setText("Length in bp (VL):");
        ui->label_25->setEnabled(false);
        ui->spinBox_VRS_min->setEnabled(false);
        ui->spinBox_VRS_max->setEnabled(false);

        // % variation
        ui->doubleSpinBox_variation_percentage->setEnabled(true);
        ui->label_28->setEnabled(true);

    } else {
        ui->spinBox_VRS_min->show();
        ui->label_11->show();
        ui->label_10->show();

        ui->spinBox_VRS_min->setEnabled(true);
        ui->spinBox_VRS_max->setEnabled(true);
        ui->label_25->setEnabled(true);
        ui->label_25->setText("Variable Positions (VP):");

        // VRS: Variable Region Size
        ui->spinBox_VRS_min->setValue(1);
        ui->spinBox_VRS_min->setMinimum(0);
        ui->spinBox_VRS_min->setMaximum(100);

        // VRS: Variable Region Size
        ui->spinBox_VRS_max->setValue(5);
        ui->spinBox_VRS_max->setMinimum(0);
        ui->spinBox_VRS_max->setMaximum(100);
        ui->frame_15->setEnabled(true);

        ui->doubleSpinBox_variation_percentage->setEnabled(false);
        ui->label_28->setEnabled(false);
    }
}


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ DOMINO marker identificiation Tab @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::on_run_DOMINO_marker_identification_button_clicked() {
    // ###############################################################
    // ##  When user press Run Button, get some options and call    ##
    // ##  the perl molecular marker identification pipeline script ##
    // ###############################################################
    ui->run_DOMINO_marker_identification_button->setEnabled(false); // Avoid user can click twice and start another computation
    // Check that user provided all the options
    MainWindow_AMM::run_DOMINO_marker_identification();
}

void MainWindow_AMM::on_pushButton_stop_DOMINO_marker_identification_clicked() {
    // ############################################################
    // ## Stops and terminates the MolMarkIdentification Process ##
    // ############################################################
    marker_identification_killed = true;
    proc_perl_DOMINO_marker_identification->kill();
    // Default tab parameters
    MainWindow_AMM::set_mapping_tab();
    ui->textEdit_input_file_DOMINO_marker->setText("");
    ui->textEdit_input_file_DOMINO_marker_2->setText("");
}

void MainWindow_AMM::on_Help_button_DOMINO_marker_Identif_tab_2_clicked() {
    // ############################################################################################
    // ## Sets Help button displaying information of the Molecular marker identification process ##
    // ############################################################################################
    about_dialog = new MyDialog(this);
    about_dialog->setToolTip("< img src = ':/images/GUI_icons/Help1.png'> Information about this Molecular Marker Identification step");
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("Information about the Molecular Marker Identification step");
    about_dialog->about_DOMINO_marker_development();
    about_dialog->show();
}

void MainWindow_AMM::on_Help_button_DOMINO_marker_Identif_tab_clicked() {
    // ############################################################################################
    // ## Sets Help button displaying information of the Molecular marker identification process ##
    // ############################################################################################
    about_dialog = new MyDialog(this);
    about_dialog->setToolTip("< img src = ':/images/GUI_icons/Help1.png'> Information about this Molecular Marker Identification step");
    about_dialog->setFixedSize(1206, 628);
    about_dialog->setWindowTitle("Information about the Molecular Marker Identification step");
    about_dialog->about_DOMINO_marker_development();
    about_dialog->show();
}

void MainWindow_AMM::checkbox_species_names_clicked() {
    int CBchecked = 0;
    QList<QCheckBox *> allbuttons = scroll_taxa_names->findChildren<QCheckBox *>();
    for (int i=0; i < allbuttons.size(); i++) {
        if (allbuttons.at(i)->isChecked()) {
            CBchecked++;
    }}
    ui->spinBox_msa_markers->setValue(CBchecked);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Actions: Connecting the different pipelines @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::onCleaningPipe_finished() {
    // #############################################################################
    // ## When CleaningPipe2 finished:                                            ##
    // ##   If user stopped it: set the tab for a new computation                 ##
    // ##   If finished properly, move to the next tab for assembling and mapping ##
    // #############################################################################

    if (cleaning_process_killed) {
        cleaning_process_killed = false;
        MainWindow_AMM::kill_CleaningProcess();
    } else {
        // Stop the progress bar and update tree View widget
        MainWindow_AMM::stop_status_bar();
        MainWindow_AMM::set_treeView_outputFolder();
        // Continue the workflow with the assembling step
        ui->tabWidget->setTabEnabled(2, true);
        ui->tabWidget->setCurrentIndex(2);
        MainWindow_AMM::set_assembling_tab();

        int type_cleaning = ui->comboBox_type_file_cleaning->currentIndex();
        if ((type_cleaning == 1) | (type_cleaning == 2)) {
            ui->comboBox_type_file_assembling->setCurrentIndex(1);
        } else if (type_cleaning == 3) {
            ui->comboBox_type_file_assembling->setCurrentIndex(1);
        } else if (type_cleaning == 4) {
            ui->comboBox_type_file_assembling->setCurrentIndex(2);
        } else if (type_cleaning == 5) {
            ui->comboBox_type_file_assembling->setCurrentIndex(2);
        } else if (type_cleaning == 6) {
            ui->comboBox_type_file_assembling->setCurrentIndex(3);
        } else if (type_cleaning == 7) {
            ui->comboBox_type_file_assembling->setCurrentIndex(3);
        }
        ui->comboBox_DOMINO_user_files->setCurrentIndex(1);
        ui->comboBox_type_file_assembling->setEnabled(false);
    }
}

void MainWindow_AMM::assembling_pipeline_finished() {
    // ##########################################################################################
    // ## When Assembling-Mapping finished properly, start the Molecular marker Identification ##
    // ##########################################################################################

    if (assembling_process_killed) {
        assembling_process_killed = false;
        MainWindow_AMM::kill_Assembling_Process();
    } else {
        MainWindow_AMM::stop_status_bar();
        MainWindow_AMM::set_treeView_outputFolder();
        ui->tabWidget->setEnabled(3);
        ui->tabWidget->setCurrentIndex(3);
        MainWindow_AMM::set_mapping_tab();
        assembly_finished = true;
        int type_assembling = ui->comboBox_type_file_assembling->currentIndex();
        if ((type_assembling == 1) | (type_assembling == 2)) {
            ui->comboBox_SINGLE_PAIR_end->setCurrentIndex(1);
        } else {
            ui->comboBox_SINGLE_PAIR_end->setCurrentIndex(0);
        }
        ui->comboBox_DOMINO_user_files_markers->setCurrentIndex(1); // DOMINO files
        ui->comboBox_SINGLE_PAIR_end->setEnabled(false);
    }
}

void MainWindow_AMM::onDOMINOmarker_identification_finished() {
    // ##############################################################################
    // ## When Molecular marker Identification finished, show the user the results ##
    // ##############################################################################

    MainWindow_AMM::stop_status_bar();
    if (marker_identification_killed) {
        marker_identification_killed = false;
        MainWindow_AMM::kill_MolMarkerIdentProcess();
    } else {
        // Show results to the user
        //cout << "I am finished here" << endl;
        marker_identification_finished = true;
        MainWindow_AMM::set_treeView_outputFolder();

        //textFileBrowser = new MyDialog(this);
        //textFileBrowser->setFixedSize(1206, 628);
        //textFileBrowser->textFileBrowser_dialog_show("DOMINO Computation is finished");
        //textFileBrowser->printing_DOMINO_last_step_finished();
        //textFileBrowser->show();

        ui->run_DOMINO_marker_identification_button->setEnabled(true);
        ui->run_DOMINO_marker_identification_button->setText("Re-Scan");
        ui->frame_12->setEnabled(true);
        ui->run_DOMINO_marker_identification_button->show();
        MainWindow_AMM::show_exit_button();
    }
}

void MainWindow_AMM::onDOMINOmarker_msa_finished() {
    MainWindow_AMM::stop_status_bar();
    if (marker_identification_killed) {
        marker_identification_killed = false;
        MainWindow_AMM::kill_MolMarkerIdentProcess();
    } else {
        // Show results to the user
        marker_identification_finished = true;
        MainWindow_AMM::set_treeView_outputFolder();
        ui->run_DOMINO_marker_identification_button->setText("Scan");
        ui->run_DOMINO_marker_identification_button->show();
        // Get names of species
        // Read file
        QString msa_taxa_file(outfolder);
        msa_taxa_file.append("/taxa_names.txt");
        QFile file(msa_taxa_file);
        if(!file.open(QIODevice::ReadOnly)) {
            cout << "Error: taxa_names.txt file could not be read...\n";
        }
        QTextStream in(&file);
        QStringList fields;
        while(!in.atEnd()) {
            QString line = in.readLine();
            fields = line.split(",");
        }
        file.close();       
        MainWindow_AMM::generate_layout(fields);
    }
}


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Actions: Writing standard output: files and dialogs @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::write_output_cleaning_perl_dialog() {
    // #####################################################################
    // ## Write output of the cleaning process into a dialog text browser ##
    // #####################################################################

    // Print Standard OUT to file and window
    QByteArray standard_out_byte = proc_Perl_cleaning_pipe->readAllStandardOutput(); // Get standard output
    //perl_dialog_log->write_perl_output(standard_out_byte);  // Redirect to Perl dialog

    // Generate a cleaning log file and append data
    QString pipe_log = QString("%1/DM_Cleaning_log.txt").arg(outfolder);
    QFile file_log_cleaning(pipe_log);
    file_log_cleaning.open(QIODevice::Append);
    if (!file_log_cleaning.isOpen()) {
        cout << "Error while opening the file for writing the output..." << endl;
    } else {
        QString standard_out_channel_finished_cleaning;
        standard_out_channel_finished_cleaning.append(standard_out_byte);
        QTextStream outstream_file_log(&file_log_cleaning);         //Point a Qtextstream object at the file
        outstream_file_log << standard_out_channel_finished_cleaning;
    }

    // Print error log too
    QByteArray standard_err_byte = proc_Perl_cleaning_pipe->readAllStandardError();
    QString standard_error_channel_finished_cleaning;
    standard_error_channel_finished_cleaning.append(standard_err_byte);
    if (!standard_error_channel_finished_cleaning.isEmpty()) {
        STD_err_dialog->showMessage(standard_error_channel_finished_cleaning);
        STD_err_dialog->setWindowTitle("DOMINO log details");
        STD_err_dialog->show();
    }
}

void MainWindow_AMM::write_output_assembling_perl_dialog() {
    // #######################################################################
    // ## Write output of the assembling process into a dialog text browser ##
    // #######################################################################

    QByteArray standard_out_byte = proc_perl_assembling_pipe->readAllStandardOutput(); // Get standard output
    //perl_dialog_log->write_perl_output(standard_out_byte);// Redirect to Perl dialog
    // Generate an assembling log file and append data
    QString pipe_log = QString("%1/DM_Assembly_log.txt").arg(outfolder);
    QFile file_log_assembling(pipe_log);
    file_log_assembling.open(QIODevice::Append);
    if (!file_log_assembling.isOpen()) {
        cout << "Error while opening the file for writing the output..." << endl;
    } else {
        QString standard_out_channel_finished_assembling;
        standard_out_channel_finished_assembling.append(standard_out_byte);
        QTextStream outstream_file_log(&file_log_assembling);         //Point a Qtextstream object at the file
        outstream_file_log << standard_out_channel_finished_assembling;
    }

    // Print error log too
    QByteArray standard_err_byte = proc_perl_assembling_pipe->readAllStandardError();
    QString standard_error_channel_finished_assembling;
    standard_error_channel_finished_assembling.append(standard_err_byte);
    if (!standard_error_channel_finished_assembling.isEmpty()) {
        STD_err_dialog->showMessage(standard_error_channel_finished_assembling);
        STD_err_dialog->setWindowTitle("DOMINO log details");
        STD_err_dialog->show();
    }
}

void MainWindow_AMM::write_output_DOMINO_marker_identification_perl_dialog() {
    // ##############################################################################################
    // ## Write output of the molecular marker identification process into a a dialog text browser ##
    // ##############################################################################################
    QByteArray standard_out_byte = proc_perl_DOMINO_marker_identification->readAllStandardOutput(); // Get standard output
    //perl_dialog_log->write_perl_output(standard_out_byte); // Redirect to Perl dialog
    // Generate an assembling log file and append data
    QString pipe_log = QString("%1/DM_MarkerScan_log.txt").arg(outfolder);
    QFile file_log_MolMarker(pipe_log);
    file_log_MolMarker.open(QIODevice::Append);
    if (!file_log_MolMarker.isOpen()) {
        cout << "Error while opening the file for writing the output..." << endl;
    } else {
        QString standard_out_channel_finished_MolMarker;
        standard_out_channel_finished_MolMarker.append(standard_out_byte);
        QTextStream outstream_file_log(&file_log_MolMarker);         //Point a Qtextstream object at the file
        outstream_file_log << standard_out_channel_finished_MolMarker;
    }

    // Print error log too
    QByteArray standard_err_byte = proc_perl_DOMINO_marker_identification->readAllStandardError();
    QString standard_error_channel_finished_marker;
    standard_error_channel_finished_marker.append(standard_err_byte);
    if (!standard_error_channel_finished_marker.isEmpty()) {
        STD_err_dialog->showMessage(standard_error_channel_finished_marker);
        STD_err_dialog->setWindowTitle("DOMINO log details");
        STD_err_dialog->show();
    }
}

//void MainWindow_AMM::show_perl_dialog() {
    // ###########################################
    // ## Shows the a dialog and a text browser ##
    // ###########################################

    //perl_dialog_log->perl_dialog();
    //perl_dialog_log->show();
//}

//void MainWindow_AMM::print_information_perl_dialog(QString command) {
    // ################################################################################
    // ## Print some information of the command sent and the project directory path  ##
    // ################################################################################
    //QString outfolder_user("Project Directory:");
    //QByteArray outfolder_user_byte = outfolder_user.toUtf8();
    //QByteArray outfolder_byte = outfolder.toUtf8();
    //QString command_perl_sent("Command line Instruction:");
    //QByteArray command_perl_sent_byte = command_perl_sent.toUtf8();
    //QString blank_line("\n");
    //QByteArray blank_line_byte = blank_line.toUtf8();
    //QByteArray command_byte = command.toUtf8();
    //perl_dialog_log->write_perl_output(outfolder_user_byte);
    //perl_dialog_log->write_perl_output(outfolder_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(command_perl_sent_byte);
    //perl_dialog_log->write_perl_output(command_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);
//}

void MainWindow_AMM::print_information_perl_command(QString command) {

    QString header("\n######################################################\n");
    QString outfolder_user("Project Directory:");
    QString command_perl_sent("Command line Instruction:");
    QString string2print;
    string2print.append(header).append(outfolder_user).append("\n").append(outfolder).append("\n").
            append(command_perl_sent).append("\n").append(command).append("\n").append(header);
    QString pipe_log = QString("%1/DM_command-line_log.txt").arg(outfolder);
    QFile file_log_MolMarker(pipe_log);
    file_log_MolMarker.open(QIODevice::Append);
    if (!file_log_MolMarker.isOpen()) {
        cout << "Error while opening the file for writing the output..." << endl;
    } else {
        QTextStream outstream_file_log(&file_log_MolMarker);         //Point a Qtextstream object at the file
        outstream_file_log << string2print;
    }
    cout << string2print.toStdString() << endl;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Actions: General    @@
// @@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::check_species_names(QStringList list_files, QString option) {

    QStringList taxa_names_provided_list;

    for (int j=0; j < list_files.size(); ++j) {
        QFileInfo File = list_files.at(j);
        QString name = File.baseName();
        QRegExp reg_id(".*id-");
        QRegExp reg_pair("_R");
        QRegExp param("arameters");
        QRegExp MID_tag("MIDtag");
        //cout << name.toStdString() << endl;
        if (name.contains(MID_tag)) {
            // skip this item
        } else if (name.contains(param) ) {
            // skip this item
        } else if (name.contains(reg_id)) {
            QStringList x;
            x = name.split(reg_id);
            QString taxa(x.at(1));
            if (taxa.contains(reg_pair)) {
                QStringList y;
                y = taxa.split(reg_pair);
                //cout << y.at(0).toStdString() << endl;
                taxa_names_provided_list << y.at(0);
            } else {
                taxa_names_provided_list << x.at(1);
        }} else {
            if (name.contains(reg_pair)) {
                QStringList y;
                y = name.split(reg_pair);
                taxa_names_provided_list << y.at(0);
            } else {
                taxa_names_provided_list << name;
    }}}

    if (option == "cleaning") {
        //taxa_names_provided_list
        taxa_names_provided_list.removeDuplicates();
        QString tmp2 = taxa_names_provided_list.join(", ");
        ui->textBrowser_taxa_names->setText(tmp2);
    }

    // Remove duplicates if any
    MainWindow_AMM::generate_layout(taxa_names_provided_list);

}

void MainWindow_AMM::generate_layout(QStringList list_ids) {

    dialog_grid = new QDialog(this);
    scroll_taxa_names = new QScrollArea(dialog_grid);
    scroll_taxa_names->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    scroll_taxa_names->setVerticalScrollBarPolicy((Qt::ScrollBarAlwaysOn));

    QWidget *view = new QWidget(dialog_grid);
    scroll_taxa_names->setWidget(view);
    scroll_taxa_names->setWidgetResizable(true);

    lay = new QGridLayout(this);
    view->setLayout(lay);

    // Remove duplicates if any
    list_ids.removeDuplicates();
    QStringList tmp(list_ids);
    int counter = 0; int row = 0; int col = 0;
    for (int i=0; i < tmp.size(); i++) {
        QCheckBox *dynamic = new QCheckBox(tmp.at(i));
        QObject::connect(dynamic, SIGNAL(clicked()), this, SLOT(checkbox_species_names_clicked()));
        if (counter == 2) {
            col = 0; counter = 0; row++;
           lay->addWidget(dynamic, row, col);
        } else {
            lay->addWidget(dynamic, row, col);
        }
        counter++; col++;
    }

    QGridLayout *grid = new QGridLayout(this);
    grid->addWidget(scroll_taxa_names);
    dialog_grid->setLayout(grid);
    dialog_grid->setToolTip("Please choose taxa to include in the analysis");
    dialog_grid->setWindowTitle("Taxa names");
    mapping_checked = true;
    ui->run_DOMINO_marker_identification_button->setEnabled(true);
}

void MainWindow_AMM::on_pushButton_clicked() {
    dialog_grid->show();
}

void MainWindow_AMM::on_spinBox_VRS_min_valueChanged(int arg1) {
    int max = ui->spinBox_VRS_max->value();
    if (arg1 >= max) {
        ui->spinBox_VRS_max->setValue(arg1);
}}

void MainWindow_AMM::on_spinBox_VRS_max_valueChanged(int arg1) {
    int min = ui->spinBox_VRS_min->value();
    if (arg1 <= min) {
        ui->spinBox_VRS_min->setValue(arg1);
}}

QString MainWindow_AMM::replace_whitespaces(QString string) {
    return string.replace(QRegExp(" "), "\\ ");
}

void MainWindow_AMM::check_user_options(QString step) {
    // ##############################################################
    // ## Checks if user has provided the minimum options required ##
    // ## for starting the computation. Then, enables Run button   ##
    // ##############################################################

    //cout << "Checking user options " + step.toStdString() << endl;

    if (step == "clean") {
        if (!input_file_name_cleaning.isEmpty()) {
            if (outfolder.isEmpty()) {
                // print error
                // Output folder is missing
                STD_err_dialog->showMessage("DOMINO output is folder missing. Please provide it in the main introductory tab");
                STD_err_dialog->setWindowTitle("Output folder missing");
                STD_err_dialog->show();
             } else {
                 if (ui->checkBox_barcode_user_file->isChecked()) {
                     QString barcode_user_file_string = ui->textEdit_user_barcode_file->toPlainText();
                     if (!barcode_user_file_string.isEmpty()) {
                         ui->run_Pipe_cleaning_button->setEnabled(true); //Enables Run button Cleaning Tab
                         cleaning_checked = true;
                     }
                 } else {
                     ui->run_Pipe_cleaning_button->setEnabled(true); //Enables Run button Cleaning Tab
                     cleaning_checked = true;
    }}}
    } else if (step == "assemble") {
         if (outfolder.isEmpty()) {
            // print error
            // Output folder is missing
            STD_err_dialog->showMessage("DOMINO output folder is missing. Please provide it in the main introductory tab");
            STD_err_dialog->setWindowTitle("Output folder missing");
            STD_err_dialog->show();
          } else {
             int index_type_file = ui->comboBox_type_file_assembling->currentIndex();   // 454 fastq, illumina or illumina pair end
             int index_combo_box = ui->comboBox_DOMINO_user_files->currentIndex();      // DOMINO clean files or user files
             if (index_type_file == 0) {
                ui->run_Assembling->setEnabled(false);
             } else {
                 if (index_combo_box == 2) {
                     if (input_file_name_assembling.isEmpty()) {
                         STD_err_dialog->showMessage("User files option has been checked but no files were provided...\n\n"
                                                     "\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\n"
                                                     "Where:\n\txxx: any character or none\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n"
                                                     "\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n");
                         STD_err_dialog->setWindowTitle("Assembling Input Files");
                         STD_err_dialog->show();
                     } else {
                         // Files provided and everything allrigth
                         QString additional_perl_options_assembling_tmp(" -user_files");
                         QString input_files_assembling(" -input_files ");
                         QStringList list_input_files_assembling = input_file_name_assembling.split("\n");
                         for (int i = 0; i < list_input_files_assembling.size(); i++ ) {
                             additional_perl_options_assembling_tmp.append(input_files_assembling).append(list_input_files_assembling.at(i));
                         }
                         additional_perl_options_assembling.append(additional_perl_options_assembling_tmp);
                         ui->run_Assembling->setEnabled(true); // Enables Run button Assembling Tab
                         assembly_checked = true;
                     }
                 } else if (index_combo_box == 1) {
                    // DOMINO user files
                     ui->run_Assembling->setEnabled(true); // Enables Run button Assembling Tab
                     assembly_checked = true;
                     ui->run_Assembling->setFocus();
          }}}
    } else if (step == "mapping") {
        additional_perl_options_mapping = "";
        if (outfolder.isEmpty()) {
           // print error
           // Output folder is missing
           STD_err_dialog->showMessage("DOMINO output folder is missing\n\n. Please provide it in the main introductory tab");
           STD_err_dialog->setWindowTitle("Output folder missing");
           STD_err_dialog->show();
        } else {
           MainWindow_AMM::setLabelData();
           int index_combo_box = ui->comboBox_DOMINO_user_files_markers->currentIndex();
           if (index_combo_box == 1) {                 // DOMINO files
               ui->run_DOMINO_marker_identification_button->setEnabled(true);
               mapping_checked = true;
           } else if (index_combo_box == 2) {                 // User contig files
                if (input_file_name_DOMINO_marker.isEmpty()) {
                    // Error message
                    STD_err_dialog->showMessage("User Contig files option has been checked but no files has been provided...\n\n"

                                                "Provide Contig FASTA file for each species.\n\n Please name your files like [xx]id-[yyy].contigs.fasta. Where:\n\n"

                                                "[xx] means it could be any character or none,\n\n"

                                                "[yyy] is any desired name for identifying the taxa contigs in these file.\n\n"

                                                );
                    STD_err_dialog->setWindowTitle("DOMINO marker User Input Files");
                    STD_err_dialog->show();
                    //ui->browse_file_DOMINO_marker->setFocus();
                } else {
                    // Check FASTQ user clean reads are provided
                    if (input_file_name_DOMINO_marker_2.isEmpty()) {
                        STD_err_dialog->showMessage(
                         "No Clean FASTQ files have been provided...\n"
                         "\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\n"
                         "Where:\n\txxx: any character or none\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n"
                         "\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n");
                        STD_err_dialog->setWindowTitle("DOMINO marker User Input Files");
                        STD_err_dialog->show();
                        //ui->browse_file_DOMINO_marker_2->setFocus();
                    } else {
                        // Get Contig files
                        QString additional_perl_options_marker_tmp("");
                        additional_perl_options_marker_tmp = QString(" -option user_assembly_contigs");
                        QString input_files_contigs_makers("");
                        input_files_contigs_makers = QString(" -user_contig_files ");
                        QStringList list_input_files_contigs_makers = input_file_name_DOMINO_marker.split("\n");
                        for (int i = 0; i < list_input_files_contigs_makers.size(); i++ ) {
                            additional_perl_options_marker_tmp.append(input_files_contigs_makers).append(list_input_files_contigs_makers.at(i));
                        }

                        // Using FASTQ user clean reads
                        QString user_cleanRead_files(" -user_cleanRead_files ");
                        QStringList list_input_files_contigs_makers_2 = input_file_name_DOMINO_marker_2.split("\n");
                        for (int i = 0; i < list_input_files_contigs_makers_2.size(); i++ ) {
                            additional_perl_options_marker_tmp.append(user_cleanRead_files).append(list_input_files_contigs_makers_2.at(i));
                        }
                        additional_perl_options_markers.append(additional_perl_options_marker_tmp);
                        ui->tabWidget->setTabEnabled(4, true); // Enable and set Molecular Marker identification tab
                        mapping_checked = true;
                        MainWindow_AMM::set_MolecularMarker_identification_tab();
                        ui->tabWidget->setCurrentIndex(4);
                        ui->run_DOMINO_marker_identification_button->setEnabled(true);
           }}
           } else if (index_combo_box == 3) { // genome
                // Check FASTQ user clean reads are provided
                if (input_file_name_DOMINO_marker.isEmpty()) {
                    STD_err_dialog->showMessage(

                     "No Genome has been provided...\n\n"
                     "Use a closed related genome for identifying molecular markers. \n\nThe name should be identified as with [xx]id-[yyy].fasta. Where:\n\n"
                     "[xx] means it could be any character or none\n\n"
                     "[yyy] is any desired name for identifying the reference genome provided. This name should not be provided within any other clean reads file.\n\n"
                                               );
                    STD_err_dialog->setWindowTitle("DOMINO marker User Input Files");
                    STD_err_dialog->show();
                    //ui->browse_file_DOMINO_marker_2->setFocus();

                } else if (input_file_name_DOMINO_marker_2.isEmpty()) {
                    STD_err_dialog->showMessage("User files option has been checked but no files were provided...\n\n"
                                                "\nPlease tag your files using: [xxx](id-)[yyy](_R[*]).fastq\n"
                                                "Where:\n\txxx: any character or none\n\tid-: Optional. If xxx is too long, please provide 'id-' to identify the name provided with [yyy]\n"
                                                "\tyyy: is any desired name for identifying the taxa reads in these file\n\t(_R[*]): if paired end files, please tag left file using R1 and right using R2\n\n");
                    STD_err_dialog->setWindowTitle("Assembling Input Files");
                    STD_err_dialog->show();
                } else {
                    // Get Contig files
                    QString additional_perl_options_marker_tmp("");
                    additional_perl_options_marker_tmp = QString(" -option genome");
                    QString input_files_contigs_makers("");
                    input_files_contigs_makers = QString(" -genome_fasta ");
                    QStringList list_input_files_contigs_makers = input_file_name_DOMINO_marker.split("\n");
                    for (int i = 0; i < list_input_files_contigs_makers.size(); i++ ) {
                        additional_perl_options_marker_tmp.append(input_files_contigs_makers).append(list_input_files_contigs_makers.at(i));
                    }
                    // Using FASTQ user clean reads
                    QString user_cleanRead_files(" -user_cleanRead_files ");
                    QStringList list_input_files_contigs_makers_2 = input_file_name_DOMINO_marker_2.split("\n");
                    for (int i = 0; i < list_input_files_contigs_makers_2.size(); i++ ) {
                        additional_perl_options_marker_tmp.append(user_cleanRead_files).append(list_input_files_contigs_makers_2.at(i));
                    }
                    additional_perl_options_markers.append(additional_perl_options_marker_tmp);

                    ui->tabWidget->setTabEnabled(4, true); // Enable and set Molecular Marker identification tab
                    mapping_checked = true;
                    MainWindow_AMM::set_MolecularMarker_identification_tab();
                    ui->tabWidget->setCurrentIndex(4);
                    ui->run_DOMINO_marker_identification_button->setEnabled(true);
                }
           } else if (index_combo_box == 4) { // Multiple alignment
               if (ui->textEdit_input_file_DOMINO_marker_3->toPlainText().isEmpty()) {
                   STD_err_dialog->showMessage("No Multiple alignment file has been provided...\n\n");
                   STD_err_dialog->setWindowTitle("DOMINO marker User Input Files");
                   STD_err_dialog->show();
               } else {
                   int counter= 0;
                   if (!ui->radioButton_select->isChecked() ) {counter++;}
                   if (!ui->radioButton_identify->isChecked()) {counter++;}
                   if (counter > 1) {
                       STD_err_dialog->showMessage("No Development Module for DOMINO is specified, please choose wether to select or to identify molecular markers");
                       STD_err_dialog->setWindowTitle("DOMINO Development Module");
                       STD_err_dialog->show();
                   } else {

                       ui->tab_markers->setEnabled(true);

                       QString additional_perl_options_marker_tmp = QString(" -option msa_alignment");
                       QString input_msa = ui->textEdit_input_file_DOMINO_marker_3->toPlainText();
                       QString input_files_option;
                       if (ui->radioButton_file->isChecked()) {
                           // msa_file
                           input_files_option = " -msa_file ";
                       } else {
                           // msa_folder
                           input_files_option = " -msa_folder ";
                       }
                       additional_perl_options_markers.append(additional_perl_options_marker_tmp).append(input_files_option).append(input_msa);
                       mapping_checked = true;

                       // Pre-comput msa_alignment
                       // Control if user specified a Perl path if not use, the default
                       if (perl_abs_path == "") { // Controlling if user specified a Perl path if not use, the default
                           perl_abs_path = ui->Perl_comboBox->currentText(); // Get current path for Perl script calling
                       }

                       QString behaviour(" -DM ");
                       if (ui->radioButton_identify->isChecked()) {
                           behaviour.append("discovery");
                           // Get some variables
                           int var_region_size_max = ui->spinBox_VRS_max->value();      // Variable Region Size max
                           //int var_region_size_min = ui->spinBox_VRS_min->value();      // Variable Region Size min
                           //QString var_region_size;
                           //var_region_size.append(QString::number(var_region_size_min)).append("::").append(QString::number(var_region_size_max));
                           //behaviour.append(scr_option).append(QString::number(con_region_size)).append(svr_option).append(var_region_size);

                           int con_region_size = ui->spinBox_CRS->value();                                                 // Conserved Region Size
                           QString scr_option(" -CL "); QString svr_option(" -VL ");
                           behaviour.append(scr_option).append(QString::number(con_region_size)).append(svr_option).append(QString::number(var_region_size_max));

                       } else if (ui->radioButton_select->isChecked()) {
                           behaviour.append("selection");
                       }

                       QString dir = QDir::currentPath();
                       proc_perl_DOMINO_marker_identification = new QProcess(this);
                       ui->pushButton_stop_DOMINO_marker_identification->setEnabled(true);
                       int num_processor = MainWindow_AMM::ui->spinBox_num_proc->value();               // Number of threads to use in the computation

                       // Set the perl command to send
                       QString MolMarker_command;
                       MolMarker_command = QString(perl_abs_path + " " + dir + "/scripts/DM_MarkerScan_v1.0.0.pl");

                       // Options
                       QString folder_option(" -o ");
                       QString folder_path = outfolder;
                       QString processor_option(" -p ");
                       double VarPercent_VR = ui->doubleSpinBox_variation_percentage->value();                         // Variability % in the Variable Region
                       QString vp_option(" -VD ");

                       if (ui->checkBox_tmp_files->isChecked()) {
                           // Do not delete temporary files
                           QString keep_files(" -TempFiles ");
                           additional_perl_options_markers.append(keep_files);
                       }

                       MolMarker_command.append(folder_option).append(folder_path).append(vp_option).append(QString::number(VarPercent_VR)).
                               append(processor_option).append(QString::number(num_processor)).
                               append(additional_perl_options_markers).append(behaviour);

                       if (debug_mode) { MolMarker_command.append(" -Debug"); }

                       // Connect the process signal finished with the next step of the pipeline
                       connect(proc_perl_DOMINO_marker_identification, SIGNAL(finished(int)), this, SLOT(onDOMINOmarker_msa_finished()));

                       // Connect the output of the process to print it in a dialog and into a file
                       connect(proc_perl_DOMINO_marker_identification, SIGNAL(readyReadStandardOutput()), this, SLOT(write_output_DOMINO_marker_identification_perl_dialog()));

                       // Start the process and print the output
                       //MainWindow_AMM::show_perl_dialog();
                       MainWindow_AMM::start_status_bar();
                       //MainWindow_AMM::print_information_perl_dialog(MolMarker_command);
                       MainWindow_AMM::print_information_perl_command(MolMarker_command);
                       proc_perl_DOMINO_marker_identification->start(MolMarker_command);
                       //cout << MolMarker_command.toStdString() << endl;
           }}
           } else if (index_combo_box == 5) { // RADseq
               if (ui->textEdit_input_file_DOMINO_marker_3->toPlainText().isEmpty()) {
                   STD_err_dialog->showMessage("No Multiple alignment file has been provided...\n\n");
                   STD_err_dialog->setWindowTitle("DOMINO marker User Input Files");
                   STD_err_dialog->show();
               } else {
                   QString additional_perl_options_marker_tmp = QString(" -option RADseq");
                   QString input_msa = ui->textEdit_input_file_DOMINO_marker_3->toPlainText();
                   QString input_files_option;
                   if (ui->radioButton_file->isChecked()) {
                       // msa_file
                       input_files_option = " -RADseq_file ";
                   }
                   additional_perl_options_markers.append(additional_perl_options_marker_tmp).append(input_files_option).append(input_msa);

                   // Control if user specified a Perl path if not use, the default
                   if (perl_abs_path == "") { // Controlling if user specified a Perl path if not use, the default
                       perl_abs_path = ui->Perl_comboBox->currentText(); // Get current path for Perl script calling
                   }

                   QString dir = QDir::currentPath();
                   proc_perl_DOMINO_marker_identification = new QProcess(this);
                   ui->pushButton_stop_DOMINO_marker_identification->setEnabled(true);

                   // Get some variables
                   int var_positions_min = ui->spinBox_VRS_min->value(); // Variable positions min
                   int var_positions_max = ui->spinBox_VRS_max->value(); // Variable positions max
                   QString var_positions_range;
                   var_positions_range.append(QString::number(var_positions_min)).append("::").append(QString::number(var_positions_max));
                   int num_processor = MainWindow_AMM::ui->spinBox_num_proc->value();               // Number of threads to use in the computation

                   // Set the perl command to send
                   QString MolMarker_command;
                   MolMarker_command = QString(perl_abs_path + " " + dir + "/scripts/DM_MarkerScan_v1.0.0.pl");

                   // Options
                   QString folder_option(" -o ");
                   QString folder_path = outfolder; QString vp_option(" -VP ");
                   QString processor_option(" -p ");


                   QString behaviour(" -DM ");
                   if (ui->radioButton_select->isChecked()) {
                       behaviour.append("selection");
                   }

                   if (ui->checkBox_tmp_files->isChecked()) {
                   // Do not delete temporary files
                        QString keep_files(" -TempFiles ");
                        additional_perl_options_markers.append(keep_files);
                   }

                    MolMarker_command.append(folder_option).append(folder_path).append(processor_option).append(vp_option).append(var_positions_range).
                           append(QString::number(num_processor)).append(additional_perl_options_markers).append(behaviour);

                    if (debug_mode) { MolMarker_command.append(" -Debug"); }


                    // Connect the process signal finished with the next step of the pipeline
                    connect(proc_perl_DOMINO_marker_identification, SIGNAL(finished(int)), this, SLOT(onDOMINOmarker_msa_finished()));

                    // Connect the output of the process to print it in a dialog and into a file
                    connect(proc_perl_DOMINO_marker_identification, SIGNAL(readyReadStandardOutput()), this, SLOT(write_output_DOMINO_marker_identification_perl_dialog()));

                    // Start the process and print the output
                    //MainWindow_AMM::show_perl_dialog();
                    MainWindow_AMM::start_status_bar();
                    //MainWindow_AMM::print_information_perl_dialog(MolMarker_command);
                    MainWindow_AMM::print_information_perl_command(MolMarker_command);
                    proc_perl_DOMINO_marker_identification->start(MolMarker_command);
                    //cout << MolMarker_command.toStdString() << endl;
}}}}}

void MainWindow_AMM::on_textEdit_perl_path_textChanged() { ui->get_user_perl_path->setEnabled(true); }

// @@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Actions: Status bar @@
// @@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::start_status_bar() {
    // ###############################################################################
    // ## Generate a status bar at the bottom to show the process is still running  ##
    // ###############################################################################

    ui->statusBar->showMessage("Running now...");
    statprogress = new QProgressBar(this);
    //ui->statusBar->addWidget(statprogress);
    ui->statusBar->addPermanentWidget(statprogress);
    statprogress->setMaximum(100);
    statprogress->setMinimum(0);
    statprogress->setTextVisible(false);
    change_value->start();
    connect(change_value, SIGNAL(value_change(int)), this, SLOT(onNumberchanged_progress_bar(int)));
}

void MainWindow_AMM::onNumberchanged_progress_bar(int percent) { statprogress->setValue(percent); }

void MainWindow_AMM::stop_status_bar() {
    // ##################################################################################
    // ## Remove the status bar at the bottom to show the process is already finished  ##
    // ##################################################################################
    ui->statusBar->showMessage("Process finished...");
    ui->statusBar->removeWidget(statprogress);
}


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Actions: Stopping a process by user input @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::kill_MolMarkerIdentProcess() {
    // ##########################################################
    // ## Sets the tab for a new MolMarkIdentification Process ##
    // ##########################################################

    MainWindow_AMM::stop_status_bar();
    MainWindow_AMM::set_treeView_outputFolder();
    ui->run_DOMINO_marker_identification_button->setEnabled(false);

    // Show a message saying Process has been aborted
    //QString command_stopped("Instruction to stop the process has been sent by the user\nAborting the process now...");
    //QByteArray command_stopped_byte = command_stopped.toUtf8();
    //QString blank_line("\n");
    //QByteArray blank_line_byte = blank_line.toUtf8();
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(command_stopped_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);

    STD_err_dialog->showMessage("Instruction to stop the process has been sent by the user\nAborting the process now...");
    STD_err_dialog->setWindowTitle("Stopping the computation");
    STD_err_dialog->show();
}

void MainWindow_AMM::kill_Assembling_Process() {
    // #######################################################
    // ## Sets the tab for a new Assembling-Mapping Process ##
    // #######################################################

    MainWindow_AMM::stop_status_bar();
    MainWindow_AMM::set_treeView_outputFolder();
    // Avoid when stop button pressed to go to next level, enable Run button.
    ui->tabWidget->setCurrentIndex(2);
    ui->run_Assembling->setEnabled(false);

    // Show a message saying Process has been aborted
    //QString command_stopped("Instruction to stop the process has been sent by the user\nAborting the process now...");
    //QByteArray command_stopped_byte = command_stopped.toUtf8();
    //QString blank_line("\n");
    //QByteArray blank_line_byte = blank_line.toUtf8();
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(command_stopped_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);

    STD_err_dialog->showMessage("Instruction to stop the process has been sent by the user\nAborting the process now...");
    STD_err_dialog->setWindowTitle("Stopping the computation");
    STD_err_dialog->show();
}

void MainWindow_AMM::kill_CleaningProcess() {
    // #############################################
    // ## Sets the tab for a new Cleaning Process ##
    // #############################################

    MainWindow_AMM::stop_status_bar();
    MainWindow_AMM::set_treeView_outputFolder();
    // Avoid when stop button pressed to go to next level, enable Run button.
    ui->tabWidget->setCurrentIndex(1);
    ui->run_Pipe_cleaning_button->setEnabled(false);
    // Show a message saying Process has been aborted
    //QString command_stopped("Instruction to stop the process has been sent by the user\nAborting the process now...");
    //QByteArray command_stopped_byte = command_stopped.toUtf8();
    //QString blank_line("\n");
    //QByteArray blank_line_byte = blank_line.toUtf8();
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);
    //perl_dialog_log->write_perl_output(command_stopped_byte);
    //perl_dialog_log->write_perl_output(blank_line_byte);

    STD_err_dialog->showMessage("Instruction to stop the process has been sent by the user\nAborting the process now...");
    STD_err_dialog->setWindowTitle("Stopping the computation");
    STD_err_dialog->show();
}


// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ Calling a QProcess @@
// @@@@@@@@@@@@@@@@@@@@@@@@
void MainWindow_AMM::cleaning_pipe(QString option) {
    // #############################################################################
    // ## Generate the command to call the perl cleaning command script pipeline  ##
    // #############################################################################

    proc_Perl_cleaning_pipe = new QProcess(this);
    proc_Perl_cleaning_pipe->setProcessChannelMode(QProcess::SeparateChannels); // Two channels: STD & STERR
    QString Species_names_mids = option;
    QString pipe_perl = QString(perl_abs_path + " " + APPDIR + "/scripts/DM_Clean_v1.0.0.pl");

    additional_perl_options_cleaning.append(Species_names_mids);
    connect(proc_Perl_cleaning_pipe, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(onCleaningPipe_finished()));

    QString len (" -l ");
    QString PH (" -s ");
    QString LENG (" -m ");
    QString proc (" -p ");
    QString out (" -outputFolder ");
    QString ad_op(" ");
    QString type_file(" -type_file ");
    int type_file_int = ui->comboBox_type_file_cleaning->currentIndex();
    QString input (" -input_file ");
    QStringList list_input_files = input_file_name_cleaning.split("\n");
    QString input_string;
    for (int i = 0; i < list_input_files.size(); i++ ) {
        // check if whitespaces anywhere on the input file name
        input_string.append(input).append(list_input_files.at(i));
    }

    QString cleaning_command = pipe_perl.append(input_string).append(type_file).append(QString::number(type_file_int)).
            append(len).append(QString::number(length_cutoff)).append(PH).append(QString::number(PHRED_sc)).append(LENG).
            append(QString::number(min_length)).append(proc).append(QString::number(num_proc)).append(out).append(outfolder).
            append(ad_op).append(additional_perl_options_cleaning);

    if (debug_mode) { cleaning_command.append(" -Debug"); }

    // Print information on the perl dialog and into a file
    connect(proc_Perl_cleaning_pipe, SIGNAL(readyReadStandardOutput()), this, SLOT(write_output_cleaning_perl_dialog())); // Every new data buffered would be appended due to readyReadStandardOutput
    //MainWindow_AMM::print_information_perl_dialog(cleaning_command);
    MainWindow_AMM::print_information_perl_command(cleaning_command);

    // Start the Process
    proc_Perl_cleaning_pipe->start(cleaning_command);
    //proc_Perl_cleaning_pipe->start("pwd");
    //QString command_sent("Command: ");
    //cout << command_sent.toStdString() << endl;
    //cout << cleaning_command.toStdString() << endl;
}

void MainWindow_AMM::run_assembling_pipeline() {
    // ######################################################################
    // ## Generate the command to call the perl assembling pipeline script ##
    // ######################################################################

    if (perl_abs_path == "") { // Controlling if user specified a Perl path if not use, the default
        perl_abs_path = ui->Perl_comboBox->currentText(); // Get current path for Perl script calling
    }
    proc_perl_assembling_pipe = new QProcess(this);                                     // Initialize a QProcess object

    // Get user options provided or default
    int min_ReadScore = MainWindow_AMM::ui->spinBox_min_read_score->value();            //Value Minimum Read Score MIRA assembly
    int processor_assembly = MainWindow_AMM::ui->spinBox_num_proc->value();  // Number of processors
    int type_file_assembling_int_tmp = ui->comboBox_type_file_assembling->currentIndex();
    int type_file_assembling_int;
    if (type_file_assembling_int_tmp == 1) {
        type_file_assembling_int = 3;
    } else if (type_file_assembling_int_tmp == 2) {
        type_file_assembling_int = 5;
    } else {
        type_file_assembling_int = 7;
    }

    // Generate perl command
    QString pipe_perl_Assembly = QString(perl_abs_path + " " + APPDIR + "/scripts/DM_Assembly_v1.0.0.pl");
    QString folder(" -o ");
    QString type_file_string(" -type_file ");
    QString min_read_score_string(" -mrs ");
    QString processors_assembly_string(" -p ");
    QString assembly_command = pipe_perl_Assembly.append(folder).append(outfolder).append(min_read_score_string).append(QString::number(min_ReadScore)).
            append(type_file_string).append(QString::number(type_file_assembling_int)).append(processors_assembly_string).append(QString::number(processor_assembly)).
            append(additional_perl_options_assembling);

    if (ui->checkBox_tmp_files->isChecked()) {
        QString tmp_files(" -TempFiles ");
        assembly_command.append(tmp_files);
    }

    if (debug_mode) { assembly_command.append(" -Debug"); }


    // Connect QProcess signal with different slots
    connect(proc_perl_assembling_pipe, SIGNAL(finished(int)), this, SLOT(assembling_pipeline_finished()));                      // Connect the process signal finished with the next step of the pipeline
    connect(proc_perl_assembling_pipe, SIGNAL(readyReadStandardOutput()), this, SLOT(write_output_assembling_perl_dialog()));   // Connect the output of the process to print it in a dialog and in a file

    // Start the process and print the output
    MainWindow_AMM::start_status_bar();
    //MainWindow_AMM::show_perl_dialog();
    //MainWindow_AMM::print_information_perl_dialog(assembly_command);
    MainWindow_AMM::print_information_perl_command(assembly_command);
    proc_perl_assembling_pipe->start(assembly_command);
    //cout << "Command sent: " << endl;
    //cout << assembly_command.toStdString() << endl;
}

void MainWindow_AMM::run_DOMINO_marker_identification() {
    // ###########################################################################################
    // ## Generate the command to call the perl Molecular marker identification pipeline script ##
    // ###########################################################################################

    QStringList list_names;
    int CBchecked = 0;
    QList<QCheckBox *> allbuttons = scroll_taxa_names->findChildren<QCheckBox *>();
    for (int i=0; i < allbuttons.size(); i++) {
        if (allbuttons.at(i)->isChecked()) {
            QString taxa_names(allbuttons.at(i)->text());
            list_names << taxa_names; CBchecked++;
    }}
    QString group_of_taxa_string = list_names.join(",");
    int option_int = ui->comboBox_DOMINO_user_files_markers->currentIndex();


    QString behaviour(" -DM ");
    if (ui->radioButton_identify->isChecked()) {
        behaviour.append("discovery");
    } else if (ui->radioButton_select->isChecked()) {
        behaviour.append("selection");
    }

    int mct_int = ui->spinBox_msa_markers->value();

    if (mct_int > CBchecked) {
        // Print error
        STD_err_dialog->showMessage("MCT Value is greater than the number of taxa specified");
        STD_err_dialog->setWindowTitle("DOMINO marker Taxa Input");
        STD_err_dialog->show();
    } else if (group_of_taxa_string.isEmpty()) {
        // Print error
        STD_err_dialog->showMessage("No taxa has been provided\nPlease provide several taxa for marker development");
        STD_err_dialog->setWindowTitle("DOMINO marker Taxa Input");
        STD_err_dialog->show();

    } else if(CBchecked == 1) {
        // Print error
        STD_err_dialog->showMessage("A single taxa has been provided\nPlease provide several taxa for marker development");
        STD_err_dialog->setWindowTitle("DOMINO marker Taxa Input");
        STD_err_dialog->show();

    } else {

        // Control if user specified a Perl path if not use, the default
        if (perl_abs_path == "") { // Controlling if user specified a Perl path if not use, the default
            perl_abs_path = ui->Perl_comboBox->currentText(); // Get current path for Perl script calling
        }
        proc_perl_DOMINO_marker_identification = new QProcess(this);
        ui->pushButton_stop_DOMINO_marker_identification->setEnabled(true);

        // Get some variables
        int num_processor = MainWindow_AMM::ui->spinBox_num_proc->value();               // Number of threads to use in the computation

        // Set the perl command to send
        QString MolMarker_command;
        MolMarker_command = QString(perl_abs_path + " " + APPDIR + "/scripts/DM_MarkerScan_v1.0.0.pl");

        // Options
        QString folder_option(" -o ");
        QString folder_path = outfolder;
        QString species_names(" -taxa_names ");
        QString processor_option(" -p ");
        QString additional_perl_options_markers_tmp;

        if (ui->keepbam_checkbox->isChecked()) {
            additional_perl_options_markers_tmp.append(" -keep_bam_file ");
        }

        if (ui->checkBox_polymorphism->isChecked()) { additional_perl_options_markers_tmp.append(" -PV "); }
        if (mct_int) {
            QString msa_option(" -MCT ");
            additional_perl_options_markers_tmp.append(msa_option).append(QString::number(mct_int));
        }

        if (ui->checkBox_tmp_files->isChecked()) {
            // Do not delete temporary files
            QString keep_files(" -TempFiles ");
            additional_perl_options_markers_tmp.append(keep_files);
        }
        if (option_int == 1) {
            QString input_option(" -option DOMINO_files ");
            additional_perl_options_markers_tmp.append(input_option);
        }
        QString tmp;
        if (option_int == 5) {
            // Get some variables
            int var_positions_min = ui->spinBox_VRS_min->value(); // Variable positions min
            int var_positions_max = ui->spinBox_VRS_max->value(); // Variable positions max
            QString vp_option(" -VP ");
            tmp.append(vp_option).append(QString::number(var_positions_min)).append("::").append(QString::number(var_positions_max));

        } else {
            QString vcr_option(" -CD ");
            QString scr_option(" -CL ");
            QString svr_option(" -VL ");
            QString divergence_option(" -VD ");
            QString read_gap_open(" -rdgopen ");
            QString read_gap_ext(" -rdgexten ");
            QString ref_gap_open(" -rfgopen ");
            QString ref_gap_ext(" -rfgexten ");
            QString mismatch(" -mp ");
            int var_region_size_max = ui->spinBox_VRS_max->value();                                                 // Variable Region Size
            //int var_region_size_min = ui->spinBox_VRS_min->value();      // Variable Region Size min
            //QString var_region_size;
            //var_region_size.append(QString::number(var_region_size_min)).append("::").append(QString::number(var_region_size_max));
            //behaviour.append(scr_option).append(QString::number(con_region_size)).append(svr_option).append(var_region_size);

            int con_region_size = ui->spinBox_CRS->value();                                                 // Conserved Region Size
            int max_variation_CR = ui->spinBox_MVCR->value();                                               // Maximun Variations in the Conserved Region
            double VarPercent_VR = ui->doubleSpinBox_variation_percentage->value();                         // Variability % in the Variable Region
            if (VarPercent_VR == 0.0) { VarPercent_VR = -1; }

            int read_gap_ext_value = ui->spinBox_pen_read_gap_ext->value();
            int read_gap_open_value = ui->spinBox_pen_read_gap_open->value();
            int ref_gap_ext_value = ui->spinBox_pen_ref_gap_ext->value();
            int ref_gap_open_value = ui->spinBox_pen_ref_gap_open->value();
            int mismatch_value = ui->spinBox_mismatch_penalty->value();                                    //
            int type = ui->comboBox_SINGLE_PAIR_end->currentIndex();                                       //Single or Pair end reads

            QString type_input(" -type_input ");
            if (type == 1) {
                type_input.append("single_end");
            } else {
                type_input.append("pair_end");
            }

            tmp.append(type_input).append(read_gap_ext).append(QString::number(read_gap_ext_value)).
                append(read_gap_open).append(QString::number(read_gap_open_value)).
                append(ref_gap_ext).append(QString::number(ref_gap_ext_value)).
                append(ref_gap_open).append(QString::number(ref_gap_open_value)).
                append(mismatch).append(QString::number(mismatch_value)).
                append(vcr_option).append(QString::number(max_variation_CR)).
                append(scr_option).append(QString::number(con_region_size)).
                //append(svr_option).append(var_region_size).
                append(svr_option).append(QString::number(var_region_size_max)).
                append(divergence_option).append(QString::number(VarPercent_VR));
        }

        MolMarker_command.append(folder_option).append(folder_path).
                append(processor_option).append(QString::number(num_processor)).append(QString(" -NPG")).
                append(species_names).append(group_of_taxa_string).append(behaviour).
                append(additional_perl_options_markers).append(tmp).append(additional_perl_options_markers_tmp);

        if (debug_mode) { MolMarker_command.append(" -Debug");}

        // Connect the process signal finished with the next step of the pipeline
        connect(proc_perl_DOMINO_marker_identification, SIGNAL(finished(int)), this, SLOT(onDOMINOmarker_identification_finished()));

        // Connect the output of the process to print it in a dialog and into a file
        connect(proc_perl_DOMINO_marker_identification, SIGNAL(readyReadStandardOutput()), this, SLOT(write_output_DOMINO_marker_identification_perl_dialog()));

        // Start the process and print the output
        //MainWindow_AMM::show_perl_dialog();
        MainWindow_AMM::start_status_bar();
        //MainWindow_AMM::print_information_perl_dialog(MolMarker_command);
        MainWindow_AMM::print_information_perl_command(MolMarker_command);
        proc_perl_DOMINO_marker_identification->start(MolMarker_command);
        //cout << MolMarker_command.toStdString() << endl;

        additional_perl_options_markers_tmp = "";
    }
    ui->run_DOMINO_marker_identification_button->setEnabled(false);
}

#ifndef MAINWINDOW_AMM_H
#define MAINWINDOW_AMM_H
#include <QMainWindow>
#include "mydialog.h"
#include "change_progress_bar.h"
#include <QProgressBar>
#include <QProcess>
#include <QDirModel>
#include "error_dialog.h"
#include "dialog_workflow.h"
#include "results_dialog.h"
#include <QScrollArea>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%       DOMINO Main Window Class        %%
// %%                                       %%
// %%       Author:                         %%
// %%   Jose F. Sanchez Herrero             %%
// %%   jfsanchezherrero@ub.edu             %%
// %%                                       %%
// %%       Institution:                    %%
// %%   Molecular Evolutionary Genetics     %%
// %%   Group, University Barcelona         %%
// %%                                       %%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace Ui {
class MainWindow_AMM;
}

class MainWindow_AMM : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow_AMM(QWidget *parent = 0);
    ~MainWindow_AMM();
    QProcess *proc_Perl_cleaning_pipe;
    QProcess *proc_perl_assembling_pipe;
    QProcess *proc_perl_DOMINO_marker_identification;
    //MyDialog *perl_dialog_log;
    change_progress_bar *change_value;
    QDirModel *model_dir;
    MyDialog *textFileBrowser;
    Error_Dialog *STD_err_dialog;
    Error_Dialog *reference_dialog;
    Dialog_workflow *workflow_dialog;
    QGridLayout *lay;
    //Results_Dialog *dialog_results;
    MyDialog *results_dialog;
    Error_Dialog *MyMessageDialog_results;
    QDialog *dialog_grid;
    QScrollArea *scroll_taxa_names;


private slots:

    // Set tabs settings
    void set_cleaning_tab();
    void set_assembling_tab();
    void set_mapping_tab();
    void set_MolecularMarker_identification_tab();
    void set_treeView_outputFolder();
    //void on_tabWidget_tabBarClicked(int index);

    // Actions
    void start_status_bar();
    void stop_status_bar();
    void onNumberchanged_progress_bar(int);
    //void show_perl_dialog();
    //void print_information_perl_dialog(QString command);
    void print_information_perl_command(QString command);
    void onCleaningPipe_finished();
    void assembling_pipeline_finished();
    void onDOMINOmarker_identification_finished();
    void onDOMINOmarker_msa_finished();
    void on_pushButton_refresh_Dir_clicked();
    void on_pushButton_openFile_clicked();
    void check_user_options(QString step);
    void kill_CleaningProcess();
    void kill_Assembling_Process();
    void kill_MolMarkerIdentProcess();
    QString replace_whitespaces(QString string);
    void show_exit_button();
    void on_pushButton_exit_DOMINO_clicked();
    void on_actionClose_triggered();
    void check_species_names(QStringList list_files, QString option);
    void checkbox_species_names_clicked();
    void generate_layout(QStringList list_ids);

    // Write log details of the process
    void write_output_cleaning_perl_dialog();
    void write_output_assembling_perl_dialog();
    void write_output_DOMINO_marker_identification_perl_dialog();

    // Introduction/About tab slots
    void on_radioButton_default_perl_path_clicked();
    void on_radioButton_custom_perl_path_clicked();
    void on_continue_button_tab_about_clicked();
    void on_browse_project_dir_initial_tab_clicked();
    void on_get_user_perl_path_clicked();
    void on_textEdit_perl_path_textChanged();
    void on_pushButton_Resume_project_clicked();
    void on_pushButton_showRef_intro_tab_clicked();
    void on_pushButton_showManual_intro_tab_clicked();

    // Cleaning tab
    void on_browse_barcode_file_clicked();
    void on_browse_databases_clicked();
    void on_show_workflow_button_clicked();
    void on_pushButton_stop_cleaning_clicked();
    void on_Help_button_cleaning_clicked();
    void on_browse_Input_file_clicked();
    void on_run_Pipe_cleaning_button_clicked();
    void on_comboBox_type_file_cleaning_currentIndexChanged(int index);
    void on_checkBox_Preprocess_NGS_files_stateChanged(int flag);
    void on_checkBox_no_BLAST_search_stateChanged(int flag);
    void on_checkBox_only_additional_database_clicked();
    void on_checkBox_only_additional_database_stateChanged(int arg1);

    // Assembling tab
    void on_run_Assembling_clicked();
    void on_pushButton_stop_assembling_clicked();
    void on_browse_Input_file_assembling_clicked();
    void on_comboBox_DOMINO_user_files_currentIndexChanged(const QString &arg1);
    void on_comboBox_type_file_assembling_currentIndexChanged(int index);
    void on_Help_button_Assembling_tab_clicked();

    // Molecular marker identification tab
    void on_run_DOMINO_marker_identification_button_clicked();
    void on_pushButton_stop_DOMINO_marker_identification_clicked();
    void on_Help_button_DOMINO_marker_Identif_tab_clicked();
    void on_Help_button_DOMINO_marker_Identif_tab_2_clicked();
    void on_browse_file_DOMINO_marker_3_clicked();
    void on_browse_file_DOMINO_marker_clicked();
    void on_comboBox_DOMINO_user_files_markers_currentIndexChanged(const QString &arg1);
    void on_browse_file_DOMINO_marker_2_clicked();
    void on_pushButton_mapping_clicked();
    void on_radioButton_identify_clicked();
    void on_radioButton_select_clicked();
    void on_spinBox_VRS_min_valueChanged(int arg1);
    void on_spinBox_VRS_max_valueChanged(int arg1);
    void setLabelData();

    // Calling a QProcess
    void cleaning_pipe(QString option);
    void run_assembling_pipeline();
    void run_DOMINO_marker_identification();

    // Menu Bar Actions
    void on_actionReference_triggered();
    void on_actionManual_triggered();
    void on_actionShow_Workflow_triggered();
    void on_actionReset_GUI_triggered();
    void on_actionAbout_Qt_triggered();
    void on_pushButton_clicked();
    void on_actionAbout_DOMINO_triggered();
    void on_actionLicense_triggered();
    void on_actionDebug_Mode_triggered();
    void on_actionCite_DOMINO_triggered();

private:
    Ui::MainWindow_AMM *ui;
    MyDialog *about_dialog;
    QProgressBar *statprogress;
};

#endif // MAINWINDOW_AMM_H

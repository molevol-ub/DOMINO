#include "results_dialog.h"
#include "ui_results_dialog.h"
#include "mydialog.h"
#include "mainwindow_amm.h"
#include <QDesktopServices>
#include <QTextStream>
#include <iostream>
using namespace std;

extern QString outfolder;
extern bool marker_identification_finished;

Results_Dialog::Results_Dialog(QWidget *parent):QDialog(parent),ui(new Ui::Results_Dialog) {
    ui->setupUi(this);
    Dialog_text_browser = new MyDialog(this);
    MyMessageDialog = new Error_Dialog(this);
    
    ui->pushButton_DOMINO_markers->setEnabled(false);
    ui->pushButton_DOMINO_markers->setToolTip("Open the molecular markers identified information in a new window");
    ui->pushButton_FilterStats->setToolTip("Open the Quality Filtering Statistics in a new window");

    if (marker_identification_finished) {
        ui->pushButton_DOMINO_markers->setEnabled(true);
    }
}

Results_Dialog::~Results_Dialog() {
    delete ui;
}

void Results_Dialog::on_pushButton_FilterStats_clicked() {

}

void Results_Dialog::on_pushButton_DOMINO_markers_clicked() {

}

void Results_Dialog::on_pushButton_Explanation_clicked() {

}

void Results_Dialog::on_pushButton_clicked() {


}

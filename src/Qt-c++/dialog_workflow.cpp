#include "dialog_workflow.h"
#include "ui_dialog_workflow.h"

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Dialog DOMINO Workflow Information Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Dialog_workflow::Dialog_workflow(QWidget *parent):QDialog(parent), ui(new Ui::Dialog_workflow) {
    ui->setupUi(this);
    Dialog_workflow::show_image("clean");
}

Dialog_workflow::~Dialog_workflow() {
    delete ui;
}

void Dialog_workflow::on_pushButton_cleaning_step_clicked() {
    Dialog_workflow::show_image("clean");
    ui->pushButton_cleaning_step->setFocus();
}

void Dialog_workflow::on_pushButton_assembling_step_clicked() {
    Dialog_workflow::show_image("assemble");
    ui->pushButton_assembling_step->setFocus();
}

void Dialog_workflow::on_pushButton_marker_step_clicked() {
    Dialog_workflow::show_image("marker");
    ui->pushButton_marker_step->setFocus();
}

void Dialog_workflow::on_pushButton_mapping_clicked() {
    Dialog_workflow::show_image("mapping");
    ui->pushButton_mapping->setFocus();
}

void Dialog_workflow::show_image(QString step) {
    QString image_path;

    if (step == "clean") {
        image_path = ":/images/GUI_icons/clean.jpg";
    }
    if (step == "assemble") {
        image_path = ":/images/GUI_icons/assembly.jpg";
    }
    if (step == "mapping") {
        image_path = ":/images/GUI_icons/mapping_alignment.jpg";
    }
    if (step == "marker") {
        image_path = ":/images/GUI_icons/marker_discovery.jpg";
    }

    //QImage image();
    QPixmap pixmap2show(image_path);
    //QSize desired_label_size = ui->label_workflow->size();
    ui->label_workflow->setPixmap(pixmap2show);
    //ui->label_workflow->setPixmap(pixmap2show.scaled(desired_label_size));
}

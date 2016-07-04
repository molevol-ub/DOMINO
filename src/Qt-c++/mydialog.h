#ifndef MYDIALOG_H
#define MYDIALOG_H
#include <QDialog>
#include "error_dialog.h"
#include <QLabel>
#include <QGridLayout>

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ General Dialog Class: Standard Output and file information @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

namespace Ui {
class MyDialog;
}

class MyDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MyDialog(QWidget *parent = 0);
    ~MyDialog();
    QLabel *imageLabel;

public slots:
    void about_cleaning_pipeline();
    void about_assembling_pipeline();
    void about_DOMINO_marker_development();
    void show_helpguide(QString helpguide);
    void write_perl_output(QByteArray ba_outsring);
    void perl_dialog();
    void textFileBrowser_dialog(QString file_data);
    void textFileBrowser_dialog_show(QString file_title);
    void printing_DOMINO_last_step_finished();
    void aboutDOMINO();
    void DOMINOlicense();

private slots:
    void on_pushButton_clicked();

private:
    Ui::MyDialog *ui;
};

#endif // MYDIALOG_H

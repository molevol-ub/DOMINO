#ifndef RESULTS_DIALOG_H
#define RESULTS_DIALOG_H

#include <QDialog>
#include "mydialog.h"
#include "error_dialog.h"
#include "mainwindow_amm.h"

namespace Ui {
class Results_Dialog;
}

class Results_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Results_Dialog(QWidget *parent = 0);
    ~Results_Dialog();



private slots:
    void on_pushButton_FilterStats_clicked();
    void on_pushButton_DOMINO_markers_clicked();
    void on_pushButton_Explanation_clicked();
    void on_pushButton_clicked();

private:
    Ui::Results_Dialog *ui;
};

#endif // RESULTS_DIALOG_H

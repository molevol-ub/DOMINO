#ifndef DIALOG_WORKFLOW_H
#define DIALOG_WORKFLOW_H
#include <QDialog>
#include <QLabel>

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Dialog DOMINO Workflow Information Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

namespace Ui {
class Dialog_workflow;
}

class Dialog_workflow : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog_workflow(QWidget *parent = 0);
    ~Dialog_workflow();

private slots:
    void on_pushButton_cleaning_step_clicked();
    void on_pushButton_assembling_step_clicked();
    void on_pushButton_marker_step_clicked();
    void show_image(QString step);

    void on_pushButton_mapping_clicked();

private:
    Ui::Dialog_workflow *ui;
};

#endif // DIALOG_WORKFLOW_H

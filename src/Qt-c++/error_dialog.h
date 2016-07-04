#ifndef ERROR_DIALOG_H
#define ERROR_DIALOG_H
#include <QDialog>

// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ Error Dialog Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@

namespace Ui {
class Error_Dialog;
}

class Error_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Error_Dialog(QWidget *parent = 0);
    ~Error_Dialog();

    void showMessage(QString Error_Message);
    void show_Reference();
    void show_BitRock_image();
    void show_GPL_image();
    void show_Qt_logo();
    void show_Perl_logo();
    void show_DOMINO_logo();

private slots:
    void on_pushButton_Error_Dialog_clicked();
    void on_pushButton_BitRock_clicked();
    void on_pushButton_GPL_clicked();
    void on_pushButton_Qt_logo_clicked();
    void on_pushButton_perl_clicked();
    void on_pushButton_DOMINO_clicked();

private:
    Ui::Error_Dialog *ui_error_dialog;
};

#endif // ERROR_DIALOG_H

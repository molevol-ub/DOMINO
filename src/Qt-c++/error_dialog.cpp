#include "error_dialog.h"
#include "ui_error_dialog.h"
#include <QDesktopServices>
#include <QMutex>
#include <QThread>

// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ Error Dialog Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@

Error_Dialog::Error_Dialog(QWidget *parent) : QDialog(parent), ui_error_dialog(new Ui::Error_Dialog) {
    ui_error_dialog->setupUi(this);
    setModal(true);
}

Error_Dialog::~Error_Dialog() {
    delete ui_error_dialog;
}

void Error_Dialog::showMessage(QString Error_Message) {
    setModal(false);
    ui_error_dialog->pushButton_DOMINO->hide();
    ui_error_dialog->pushButton_BitRock->hide();
    ui_error_dialog->pushButton_GPL->hide();
    ui_error_dialog->pushButton_Qt_logo->hide();
    ui_error_dialog->pushButton_perl->hide();
    ui_error_dialog->textBrowser_Error_Dialog->setText(Error_Message);
    ui_error_dialog->textBrowser_Error_Dialog->setFixedSize(621,401);
}

void Error_Dialog::on_pushButton_Error_Dialog_clicked() {
    close();
}

void Error_Dialog::show_Reference() {
    ui_error_dialog->textBrowser_Error_Dialog->hide();

    Error_Dialog::show_DOMINO_logo();
    ui_error_dialog->pushButton_DOMINO->setToolTip("DOMINO website");

    Error_Dialog::show_BitRock_image();
    ui_error_dialog->pushButton_BitRock->setToolTip("BitRock Installer Builder for Qt apps");

    Error_Dialog::show_GPL_image();
    ui_error_dialog->pushButton_GPL->setToolTip("General Public License");

    Error_Dialog::show_Qt_logo();
    ui_error_dialog->pushButton_Qt_logo->setToolTip("Qt Framework");

    Error_Dialog::show_Perl_logo();
    ui_error_dialog->pushButton_perl->setToolTip("Perl scripting language");

    ui_error_dialog->pushButton_Error_Dialog->move(500,430);
}

void Error_Dialog::show_BitRock_image() {
    QPixmap BitRock(":/images/GUI_icons/installersby_BitRock.png");
    QIcon icon(BitRock);
    ui_error_dialog->pushButton_BitRock->setIcon(icon);
    ui_error_dialog->pushButton_BitRock->setIconSize(BitRock.size());
    ui_error_dialog->pushButton_BitRock->show();
}

void Error_Dialog::show_GPL_image() {
    QPixmap GPL(":/images/GUI_icons/gplv3-127x51.png");
    QIcon icon(GPL);
    ui_error_dialog->pushButton_GPL->setIcon(icon);
    ui_error_dialog->pushButton_GPL->setIconSize(GPL.size());
    ui_error_dialog->pushButton_GPL->show();
}

void Error_Dialog::on_pushButton_BitRock_clicked() {
    QString website_link = "http://installbuilder.bitrock.com/index.html";
    QDesktopServices::openUrl(QUrl(website_link));
}

void Error_Dialog::on_pushButton_GPL_clicked() {
    QString website_link = "https://www.gnu.org/licenses/quick-guide-gplv3.html";
    QDesktopServices::openUrl(QUrl(website_link));
}

void Error_Dialog::show_Perl_logo() {
    QPixmap image(":/images/GUI_icons/Perl.gif");
    QIcon icon(image);
    ui_error_dialog->pushButton_perl->setIcon(icon);
    ui_error_dialog->pushButton_perl->setIconSize(image.size());
    ui_error_dialog->pushButton_perl->show();
}

void Error_Dialog::show_Qt_logo() {
    QPixmap image(":/images/GUI_icons/logo_qt.png");
    QIcon icon(image);
    ui_error_dialog->pushButton_Qt_logo->setIcon(icon);
    ui_error_dialog->pushButton_Qt_logo->setIconSize(image.size());
    ui_error_dialog->pushButton_Qt_logo->show();
}

void Error_Dialog::show_DOMINO_logo() {
    QPixmap image(":/images/GUI_icons/DOMINO_logo.png");
    QIcon icon(image);
    ui_error_dialog->pushButton_DOMINO->setIcon(icon);
    ui_error_dialog->pushButton_DOMINO->setIconSize(image.size());
    ui_error_dialog->pushButton_DOMINO->show();
}

void Error_Dialog::on_pushButton_Qt_logo_clicked() {
    QString website_link = "http://www.qt.io/";
    QDesktopServices::openUrl(QUrl(website_link));
}

void Error_Dialog::on_pushButton_perl_clicked()
{
    QString website_link = "https://www.perl.org/";
    QDesktopServices::openUrl(QUrl(website_link));
}

void Error_Dialog::on_pushButton_DOMINO_clicked() {
    QString website_link = "http://www.ub.edu/softevol/domino/";
    QDesktopServices::openUrl(QUrl(website_link ));
}

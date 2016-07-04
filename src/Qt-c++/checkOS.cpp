#include "checkOS.h"
#include <QProcess>
#include <QString>
#include <iostream>
#include <QCoreApplication>
#include <QDir>
#include "error_dialog.h"

using namespace std;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Check Operating Sytem Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

QString Opesystem;
QString APPDIR;
QString perl_path;

checkOS::checkOS() {

    /*
     * Just running the command change directory, "pwd", we would be able to know wether the user will be using Windows or Unix.
     * If the user is running a windows system, the directory will be starting with C:/ that would let us make the difference
     * between the Operating system we are interested and get where perl is.
     *
    */

    // First of all we would check the operating system we are running
    QString pwd;
    pwd = QCoreApplication::applicationDirPath();
    //pwd = QDir::currentPath();

    QString regexp_mac("/DOMINO.app");

    if (pwd.contains(regexp_mac)) {
        QStringList x;
        x = pwd.split(regexp_mac);
        pwd = x.at(0);
    }

    // We would get the directory where our application is running
    APPDIR = pwd.trimmed();
    //cout << "Working directory : " << APPDIR.toStdString() << endl;

    QString regexp("C:/");

    if (pwd.startsWith(regexp)) {
        //cout << "This is a Windows operating system" << endl;
        Opesystem = QString ("Windows");

    }else {

        // We have noticed that there is a difference between Fedora and Ubuntu or some other Linux distributions, if we run the command lets say, whereis perl,
        // it would print in Fedora perl: /usr/bin/perl but in Mac and the rest of Unix distributions it would just print /usr/bin/perl
        // If we want to take the path for perl or whatever it would be necessary to take this into account. Also, this happens in Debian

        QProcess Unix_dist;
        Unix_dist.setProcessChannelMode(QProcess::MergedChannels);  //captures both standard error and stdoutput in the same
        QString whereis("whereis perl");
        Unix_dist.start(whereis); //executing the command
        Unix_dist.waitForFinished(-1); //waiting until the process has finished
        QString result_dist = Unix_dist.readAllStandardOutput();
        QString distribution(result_dist); //converts bytearray on string
        QString regexp_perl("perl:");

        if (distribution.startsWith(regexp_perl)) {
            //cout << "This is a Fedora operating system" << endl;
            Opesystem = QString ("Fedora");
        } else {
            //cout << "This is a Unix operating system" << endl;
            Opesystem = QString ("Unix");
        }
    }

    //Get the perl absolute executable path from the Operating System
    checkOS::get_perl_path();
}

void checkOS::get_perl_path() {

    STD_err_dialog = new Error_Dialog(this);
    STD_err_dialog->resize(600,250);

    if (Opesystem == "Unix") {  //This would work for Linux, and most Unix distribution but not for fedora neither ubuntu
        //Once we know we are running a Unix operating system, we would check the paths for the languages we are interested
        QProcess envir_perl;
        envir_perl.setProcessChannelMode(QProcess::MergedChannels);  //captures both standard error and stdoutput in the same
        QString whereis_perl("whereis perl");
        envir_perl.start(whereis_perl); //executing the command
        envir_perl.waitForFinished(-1); //waiting until the process has finished
        QString result_perl = envir_perl.readAllStandardOutput();
        QString command_perl(result_perl); //converts bytearray on string

        if (command_perl == "") {
            STD_err_dialog->showMessage("Not able to find any default perl executable. Please specify it in the box");
            STD_err_dialog->setWindowTitle("Error ocurred while retrieving perl information");
            STD_err_dialog->show();
        } else {
            perl_path = QString(command_perl);
            perl_path = perl_path.trimmed();
            //cout << perl_path.toStdString() << endl;
        }
      } else if (Opesystem == "Fedora") { // or Ubuntu
         /* Fedora checking loop
          * perl: /usr/bin/perl /usr/share/man/man1/perl.1.gz
          * Check how to get the proper the path
          */
         //perl path
         // We would try to check if perl is already running
         QProcess perl_fedora;
         perl_fedora.setProcessChannelMode(QProcess::MergedChannels);  //captures both standard error and stdoutput in the same
         QString perl_version("perl -v");
         perl_fedora.start(perl_version); //executing the command
         perl_fedora.waitForFinished(-1); //waiting until the process has finished
         QString result_fedora = perl_fedora.readAllStandardOutput();
         QString command_fedora(result_fedora); //converts bytearray on string

         if (command_fedora == "") { // There is no perl installed in this computer, show error dialog
            STD_err_dialog->showMessage("Not able to find any default perl executable. Please specify it in the box");
            STD_err_dialog->setWindowTitle("Error ocurred while retrieving perl information");
            STD_err_dialog->show();
         } else {
             QProcess envir_perl_fedora;
             envir_perl_fedora.setProcessChannelMode(QProcess::MergedChannels);  //captures both standard error and stdoutput in the same
             QString where_perl_fedora("whereis perl");
             envir_perl_fedora.start(where_perl_fedora); //executing the command
             envir_perl_fedora.waitForFinished(-1); //waiting until the process has finished
             QString result_perl_fedora = envir_perl_fedora.readAllStandardOutput();
             QString command_perl_fedora(result_perl_fedora); //converts bytearray on string
             QString perl_path_fedora = QString(command_perl_fedora);
             QStringList perl_paths = perl_path_fedora.split(" ");
             perl_path = perl_paths.join("\n");

             perl_path = perl_path.trimmed();
             //cout << perl_path_fedora.toStdString() << endl;
             //cout << perl_path.toStdString() << endl;
             // perl: /usr/bin/perl /etc/perl /usr/lib/perl /usr/bin/X11/perl /usr/share/perl /usr/share/man/man1/perl.1.gz
             // If we have got something when calling perl -v means that there is already a version installed and running in the computer, so we would just call it like perl
             //perl_path = QString("perl");     // Maybe because of a List and just a name it doesnt work, but QProcess needs a list
         }
    } else if (Opesystem == "Windows") {
           //If we are not running a Unix computer, we would be running a Windows computer
           //perl path
           QProcess envir_perl_wind;
           envir_perl_wind.setProcessChannelMode(QProcess::MergedChannels);  //captures both standard error and stdoutput in the same
           QString where_perl_wind("where perl");
           envir_perl_wind.start(where_perl_wind); //executing the command
           envir_perl_wind.waitForFinished(-1); //waiting until the process has finished
           QString result_perl_wind = envir_perl_wind.readAllStandardOutput();
           QString command_perl_wind(result_perl_wind); //converts bytearray on string
           if (command_perl_wind == "") {
                STD_err_dialog->showMessage("Not able to find any default perl executable. Please specify it in the box");
                STD_err_dialog->setWindowTitle("Error ocurred while retrieving perl information");
                STD_err_dialog->show();
           } else {
               perl_path = QString(command_perl_wind);
               perl_path = perl_path.trimmed();
           }
    }
}

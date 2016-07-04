#include "mainwindow_amm.h"
#include <QApplication>
#include "checkOS.h"
#include <QCoreApplication>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%           DOMINO Main Class           %%
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

int main(int argc, char *argv[])
{
    QCoreApplication::addLibraryPath("./plugins"); //Load Plugins
    QApplication DOMINO(argc, argv); // Set application and main window name and properties

    QIcon icon(":/icon/GUI_icons/domino_pieza.ico");
    DOMINO.setWindowIcon(QIcon(icon));

    MainWindow_AMM main_w_DOMINO;
    checkOS();
    main_w_DOMINO.setFixedSize(823, 674);     //Set style for the main window of the application

    main_w_DOMINO.show();
    return DOMINO.exec();
}

#ifndef CHECKOS_H
#define CHECKOS_H
#include <mainwindow_amm.h>
#include <error_dialog.h>

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@ Check Operating Sytem Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class checkOS : public MainWindow_AMM {
public:
    checkOS();
    void get_perl_path();

    Error_Dialog *STD_err_dialog;
};

#endif // CHECKOS_H

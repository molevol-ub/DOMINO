#include "change_progress_bar.h"
#include "mainwindow_amm.h"
#include <QMutex>

// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ Progress Bar Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@
change_progress_bar::change_progress_bar(QObject *parent) : QThread(parent) {
}

void change_progress_bar::run() {
    int i = 0;
    for (i = 0; i < 100; i++) {
        QMutex mutex;
        mutex.lock();
        QThread::msleep(50);
        mutex.unlock();
        emit value_change(i);
        if (i == 99) { i = 0; }
    }
}


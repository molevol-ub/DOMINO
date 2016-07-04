#-------------------------------------------------
#
# Project created by QtCreator 2014-04-28T18:04:30
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QMAKE_MAC_SDK=macosx10.11

CONFIG += static

ICON = domino.icns

TARGET = DOMINO
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow_amm.cpp \
    mydialog.cpp \
    change_progress_bar.cpp \
    checkOS.cpp \
    error_dialog.cpp \
    dialog_workflow.cpp

HEADERS  += mainwindow_amm.h \
    mydialog.h \
    change_progress_bar.h \
    checkOS.h \
    error_dialog.h \
    dialog_workflow.h

FORMS    += mainwindow_amm.ui \
    mydialog.ui \
    error_dialog.ui \
    dialog_workflow.ui

RESOURCES += \
    Resource_file.qrc



#ifndef CHANGE_PROGRESS_BAR_H
#define CHANGE_PROGRESS_BAR_H
#include <QThread>

// @@@@@@@@@@@@@@@@@@@@@@@@
// @@ Progress Bar Class @@
// @@@@@@@@@@@@@@@@@@@@@@@@

class change_progress_bar : public QThread
{
    Q_OBJECT
public:
    explicit change_progress_bar(QObject *parent = 0);
    void run();

signals:
    void value_change(int percentage);

};

#endif // CHANGE_PROGRESS_BAR_H

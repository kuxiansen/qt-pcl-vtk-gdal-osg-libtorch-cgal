#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include <qfiledialog.h>
#include"ui_dpl_semseg.h"
class dpl_sem_seg :public QDialog
{
	Q_OBJECT
public:
	dpl_sem_seg(QWidget* parent = nullptr);
	~dpl_sem_seg();
signals:
	void send_block(QString blocksize, QString blockstride, QString modelpath);
public slots:
	void on_accepted();
	void on_canceled();
	void on_open_model();
private:
	Ui::dpl_semsegclass* ui;
};


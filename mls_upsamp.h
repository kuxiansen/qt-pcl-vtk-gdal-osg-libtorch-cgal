#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_mls_upsampling.h"

class mls_upsamp :public QDialog
{
	Q_OBJECT

public:
	mls_upsamp(QWidget* parent = nullptr);
	~mls_upsamp();
signals:
	void send_mlsup(QString radius, QString step);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::mls_upsamplingclass* ui;
};




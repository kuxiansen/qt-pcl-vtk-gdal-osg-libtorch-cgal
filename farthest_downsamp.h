#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_farthest_downsampling.h"


class farthest_downsamp :public QDialog
{
	Q_OBJECT

public:
	farthest_downsamp(QWidget* parent = nullptr);
	~farthest_downsamp();
signals:
	void send_farthestsamp(QString nums);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::farthest_downsampclass* ui;
};

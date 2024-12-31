#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include<QButtonGroup>
#include"ui_cur_downsampling.h"

class cur_downsamp :public QDialog
{
	Q_OBJECT

public:
	cur_downsamp(QWidget* parent = nullptr);
	~cur_downsamp();
signals:
	void send_cursamp(QString K_search, QString nums);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::cur_dowwnsamp* ui;
};

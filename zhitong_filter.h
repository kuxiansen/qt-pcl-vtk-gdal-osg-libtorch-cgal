#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_zhitong_filter.h"

class zhitong_filter  : public QDialog
{
	Q_OBJECT

public:
	zhitong_filter(QWidget*parent = nullptr);
	~zhitong_filter();
public slots:
	void on_accepted();
	void on_canceled();
signals:
	void sendcanshu(int data1, QString data2,QString data3);//参数，最小值，最大值
private:
	Ui::Dialog* ui;
};

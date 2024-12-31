#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_gauss_filter.h"



class gaussian_filter :public QDialog
{
	Q_OBJECT

public:
	gaussian_filter(QWidget* parent = nullptr);
	~gaussian_filter();
signals:
	void send_gaussian(QString sigma, QString sigmarela, QString threshold, QString radius);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::gauss_filter* ui;
};



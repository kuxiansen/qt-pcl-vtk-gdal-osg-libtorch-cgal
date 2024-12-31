#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_doubleline_filter.h"

class doubleline_filter :public QDialog
{
	Q_OBJECT

public:
	doubleline_filter(QWidget* parent = nullptr);
	~doubleline_filter();
signals:
	void send_doubleline(QString nums, QString sigmadis, QString sigmanormal);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::doubleline_filterclass* ui;
};

#pragma once
//ÉèÖÃÖĞÎÄ±àÂë
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include "ui_view_slider.h"


class View_slider :public QDialog
{
	Q_OBJECT
public:
	View_slider(QWidget* parent = nullptr);
	~View_slider();
signals:
	void senddata(QString data);
public slots:
	void on_buttonBox_rejected();
	void on_buttonBox_accepted();
private:
	Ui::View_sliderclass* ui;
};
#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_morphological_filter.h"
class morphological_filter :public QDialog
{
	Q_OBJECT

public:
	morphological_filter(QWidget* parent = nullptr);
	~morphological_filter();
signals:
	void send_morph(QString grid_max, QString slope, QString height_init, QString height_max);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::morph_filter* ui;
};


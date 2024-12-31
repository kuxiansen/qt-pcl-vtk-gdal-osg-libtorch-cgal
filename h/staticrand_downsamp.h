#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_staticrand_downsampling.h"

class staticrandownsamp :public QDialog
{
	Q_OBJECT

public:
	staticrandownsamp(QWidget* parent = nullptr);
	~staticrandownsamp();
signals:
	void send_nums(QString nums);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::staticrd_ds* ui;
};
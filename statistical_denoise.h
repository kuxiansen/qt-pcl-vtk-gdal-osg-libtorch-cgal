#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include<QButtonGroup>
#include"ui_statistical_denoising.h"

class statiscal_denoise :public QDialog
{
	Q_OBJECT

public:
	statiscal_denoise(QWidget* parent = nullptr);
	~statiscal_denoise();
signals:
	void send_std(QString d_max, QString nums, int id);
public slots:
	void on_accepted();
	void on_canceled();
	//void onRadioButtonClicked(int id);
private:
	Ui::stastical_denos* ui;
	QButtonGroup* buttonGroup;
};
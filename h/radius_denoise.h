#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include<QButtonGroup>
#include"ui_radius_denoising.h"

class radius_denoise :public QDialog
{
	Q_OBJECT

public:
	radius_denoise(QWidget* parent = nullptr);
	~radius_denoise();
signals:
	void send_radius_nois(QString radius, QString nums);
public slots:
	void on_accepted();
	void on_canceled();
	//void onRadioButtonClicked(int id);
private:
	Ui::radius_denos* ui;	
};


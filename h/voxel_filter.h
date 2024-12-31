#pragma once
#pragma execution_character_set("utf-8")
#include<QDialog>
#include <Qdebug>
#include"ui_voxel_filter.h"
class voxel_filter :public QDialog
{
	Q_OBJECT

public:
	voxel_filter(QWidget* parent = nullptr);
	~voxel_filter();
signals:
	void send_voxel(QString voxelnum, QString voxelsize);
public slots:
	void on_accepted();
	void on_canceled();
private:
	Ui::voxel_filterclass* ui;
};

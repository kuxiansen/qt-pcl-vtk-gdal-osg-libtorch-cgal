#include "voxel_filter.h"

voxel_filter::voxel_filter(QWidget* parent )
	: QDialog(parent), ui(new Ui::voxel_filterclass)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
}

voxel_filter::~voxel_filter()
{
	delete ui;
}
void voxel_filter::on_canceled()
{
	this->close();
}
void voxel_filter::on_accepted()
{
	emit send_voxel(ui->lineEdit_2->text(),ui->lineEdit->text());
	this->close();
}


#include "radius_denoise.h"

radius_denoise::radius_denoise(QWidget* parent)
	: QDialog(parent), ui(new Ui::radius_denos)
{

	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

radius_denoise::~radius_denoise()
{
	delete ui;
}

void radius_denoise::on_canceled()
{
	this->close();
}

void radius_denoise::on_accepted()
{
	emit send_radius_nois(ui->lineEdit_radius->text(), ui->lineEdit_minums->text());
	this->close();
}

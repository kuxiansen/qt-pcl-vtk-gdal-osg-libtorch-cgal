#include "farthest_downsamp.h"


farthest_downsamp::farthest_downsamp(QWidget* parent)
	: QDialog(parent), ui(new Ui::farthest_downsampclass)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

farthest_downsamp::~farthest_downsamp()
{
	delete ui;
}

void farthest_downsamp::on_canceled()
{
	this->close();
}

void farthest_downsamp::on_accepted()
{
	emit send_farthestsamp(ui->lineEdit->text());
	this->close();
}

#include "gaussian_filter.h"
gaussian_filter::gaussian_filter(QWidget* parent)
	: QDialog(parent), ui(new Ui::gauss_filter)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

gaussian_filter::~gaussian_filter()
{
	delete ui;
}

void gaussian_filter::on_canceled()
{
	this->close();
}

void gaussian_filter::on_accepted()
{
	emit send_gaussian(ui->lineEdit_sigma->text(), ui->lineEdit_sigmarela->text(), ui->lineEdit_threshold->text(),ui->lineEdit_radius->text());
	this->close();
}
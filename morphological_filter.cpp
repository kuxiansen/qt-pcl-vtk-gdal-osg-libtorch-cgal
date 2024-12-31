#include "morphological_filter.h"

morphological_filter::morphological_filter(QWidget* parent)
	: QDialog(parent), ui(new Ui::morph_filter)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

morphological_filter::~morphological_filter()
{
	delete ui;
}

void morphological_filter::on_canceled()
{
	this->close();
}

void morphological_filter::on_accepted()
{
	emit send_morph(ui->lineEdit_maxfrid->text(), ui->lineEdit_slope->text(), ui->lineEdit_indis->text(), ui->lineEdit_maxdis->text());
	this->close();
}
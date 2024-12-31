#include "doubleline_filter.h"
doubleline_filter::doubleline_filter(QWidget* parent)
	: QDialog(parent), ui(new Ui::doubleline_filterclass)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

doubleline_filter::~doubleline_filter()
{
	delete ui;
}

void doubleline_filter::on_canceled()
{
	this->close();
}

void doubleline_filter::on_accepted()
{
	emit send_doubleline(ui->lineEdit_num->text(), ui->lineEdit_dissigma->text(), ui->lineEdit_norsigma->text());
	this->close();
}
#include "cur_downsamp.h"

cur_downsamp::cur_downsamp(QWidget* parent)
	: QDialog(parent), ui(new Ui::cur_dowwnsamp)
{

	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

cur_downsamp::~cur_downsamp()
{
	delete ui;
}

void cur_downsamp::on_canceled()
{
	this->close();
}

void cur_downsamp::on_accepted()
{
	emit send_cursamp(ui->lineEdit_num->text(), ui->lineEdit_cur->text());
	this->close();
}

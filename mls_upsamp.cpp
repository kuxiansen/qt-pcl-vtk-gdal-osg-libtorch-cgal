#include "mls_upsamp.h"
mls_upsamp::mls_upsamp(QWidget* parent)
	: QDialog(parent), ui(new Ui::mls_upsamplingclass)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

mls_upsamp::~mls_upsamp()
{
	delete ui;
}

void mls_upsamp::on_canceled()
{
	this->close();
}

void mls_upsamp::on_accepted()
{
	emit send_mlsup(ui->lineEdit_radius->text(), ui->lineEdit_step->text());
	this->close();
}

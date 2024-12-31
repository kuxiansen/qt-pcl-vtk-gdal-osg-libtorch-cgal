#include "statistical_denoise.h"

statiscal_denoise::statiscal_denoise(QWidget* parent)
	: QDialog(parent), ui(new Ui::stastical_denos)
{
	
	ui->setupUi(this);
	buttonGroup = new QButtonGroup(this);
	buttonGroup->addButton(ui->single_thread, 1);
	buttonGroup->addButton(ui->many_thread, 2);
	ui->single_thread->setChecked(true);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	//connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(onRadioButtonClicked(int)));
}

statiscal_denoise::~statiscal_denoise()
{
	delete ui;
}

void statiscal_denoise::on_canceled()
{
	this->close();
}

void statiscal_denoise::on_accepted()
{
	emit send_std(ui->lineEdit->text(), ui->lineEdit_2->text(), buttonGroup->checkedId());
	this->close();
}

#include "zhitong_filter.h"

zhitong_filter::zhitong_filter(QWidget*parent)
	: QDialog(parent), ui(new Ui::Dialog)
{
	ui->setupUi(this);
	ui->comboBox->clear();
	QStringList list;
	list << "x" << "y" << "z" << "label";
	ui->comboBox->addItems(list);
	ui->comboBox->setCurrentIndex(0);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
}

zhitong_filter::~zhitong_filter()
{
	delete ui;
}

void zhitong_filter::on_canceled()
{
	this->close();
}

void zhitong_filter::on_accepted()
{
	emit sendcanshu(ui->comboBox->currentIndex(), ui->min_num->text(), ui->lineEdit->text());
	this->close();
}

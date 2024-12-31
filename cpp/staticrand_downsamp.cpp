#include"staticrand_downsamp.h"

staticrandownsamp::staticrandownsamp(QWidget* parent)
	: QDialog(parent), ui(new Ui::staticrd_ds)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
}

staticrandownsamp::~staticrandownsamp()
{
	delete ui;
}

void staticrandownsamp::on_canceled()
{
	this->close();
}
void staticrandownsamp::on_accepted()
{
	emit send_nums(ui->lineEdit->text());
	this->close();
}

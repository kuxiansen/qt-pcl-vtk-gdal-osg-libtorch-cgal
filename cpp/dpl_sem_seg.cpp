#include "dpl_sem_seg.h"


dpl_sem_seg::dpl_sem_seg(QWidget* parent) : QDialog(parent), ui(new Ui::dpl_semsegclass)
{
	ui->setupUi(this);
	connect(ui->okButton, SIGNAL(clicked()), this, SLOT(on_accepted()));
	connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(on_canceled()));
	connect(ui->open_model, SIGNAL(clicked()), this, SLOT(on_open_model()));
}

dpl_sem_seg::~dpl_sem_seg()
{
	delete ui;
}
void dpl_sem_seg::on_accepted()
{
	emit send_block(ui->block_size->text(), ui->block_stride->text(), ui->model_file->text());
	this->close();
}

void dpl_sem_seg::on_canceled()
{
	this->close();
}

void dpl_sem_seg::on_open_model()
{
	QString path = QFileDialog::getOpenFileName(this, "选择模型", ".//", "模型文件(*.pt);;所有文件(*.*)");
	ui->model_file->setText(path);
}

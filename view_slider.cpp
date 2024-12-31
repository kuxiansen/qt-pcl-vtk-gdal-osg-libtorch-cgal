#include"view_slider.h"
#include "ui_view_slider.h"
//#include <qmessagebox.h>

View_slider::View_slider(QWidget* parent) :
    QDialog(parent), ui (new Ui::View_sliderclass)
{
    
    ui->setupUi(this);
    ui->spinBox->setMinimum(0);
    ui->spinBox->setMaximum(10);
    ui->spinBox->setSingleStep(1);
    ui->spinBox->setValue(5);

    ui->ptsize_slider->setMinimum(0);
    ui->ptsize_slider->setMaximum(10);
    ui->ptsize_slider->setSingleStep(1);
    ui->ptsize_slider->setValue(5);

    connect(ui->spinBox, SIGNAL(valueChanged(int)), ui->ptsize_slider, SLOT(setValue(int)));
    connect(ui->ptsize_slider, SIGNAL(valueChanged(int)), ui->spinBox, SLOT(setValue(int)));
    //connect(ui->buttonBox, SIGNAL(triggered(bool)), this, SLOT(on_buttonBox_accepted()));//点云大小
}

View_slider::~View_slider()
{
    delete ui;
    //QMessageBox::warning(this, "提示", "修改点云尺寸成功");
}

void View_slider::on_buttonBox_accepted()
{
    QString value = QString("%1").arg(ui->ptsize_slider->value());
    emit senddata(value);
    this->close();
}

void View_slider::on_buttonBox_rejected()
{
    this->close();
}




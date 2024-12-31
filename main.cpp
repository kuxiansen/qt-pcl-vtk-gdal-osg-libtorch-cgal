#include "QTPCLVTK.h"
#include <QtWidgets/QApplication>
#include <vtkOutputWindow.h>
#include <QMovie>
#include <QPixmap>
#include <QSplashScreen>

int main(int argc, char *argv[])
{
    vtkOutputWindow::SetGlobalWarningDisplay(0);//������vtkOutputWindow����
    QApplication a(argc, argv);
    QPixmap pixmap(":/img/image/start.gif");
    QSplashScreen splash(pixmap); //
    splash.setWindowOpacity(0.8); // ���ô���͸����
    QLabel label(&splash);
    QMovie mv(":/img/image/start.gif");
    label.setMovie(&mv);
    mv.start();
    splash.show();
    //splash.showMessage("�������ڼ���......", Qt::AlignCenter, Qt::red); //��ʾ����
    
    for (int i = 0; i < 5000; i += mv.speed()) 
    {
        a.processEvents(); //ʹ��������ʾ���������ͬʱ������Ӧ���������¼�
        Sleep(mv.speed()); // ��ʱ
    }
    QTPCLVTK w;
    w.setWindowTitle("PCL_test");
    w.show();
    splash.finish(&w); //����������ʼ����ɺ������������
    return a.exec();
}

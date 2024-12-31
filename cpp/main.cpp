#include "QTPCLVTK.h"
#include <QtWidgets/QApplication>
#include <vtkOutputWindow.h>
#include <QMovie>
#include <QPixmap>
#include <QSplashScreen>

int main(int argc, char *argv[])
{
    vtkOutputWindow::SetGlobalWarningDisplay(0);//不弹出vtkOutputWindow窗口
    QApplication a(argc, argv);
    QPixmap pixmap(":/img/image/start.gif");
    QSplashScreen splash(pixmap); //
    splash.setWindowOpacity(0.8); // 设置窗口透明度
    QLabel label(&splash);
    QMovie mv(":/img/image/start.gif");
    label.setMovie(&mv);
    mv.start();
    splash.show();
    //splash.showMessage("程序正在加载......", Qt::AlignCenter, Qt::red); //显示文字
    
    for (int i = 0; i < 5000; i += mv.speed()) 
    {
        a.processEvents(); //使程序在显示启动画面的同时仍能响应鼠标等其他事件
        Sleep(mv.speed()); // 延时
    }
    QTPCLVTK w;
    w.setWindowTitle("PCL_test");
    w.show();
    splash.finish(&w); //在主体对象初始化完成后结束启动动画
    return a.exec();
}

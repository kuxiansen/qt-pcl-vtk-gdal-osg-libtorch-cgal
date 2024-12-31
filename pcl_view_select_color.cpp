#include "pcl_view_select_color.h"
#include <QHBoxLayout>
#include <qcolordialog.h>

pcl_view_select_color::pcl_view_select_color()
{
    QColor c = QColorDialog::getColor(Qt::white, this);

    if (c.isValid())
    {
        color = c;
    }
}

pcl_view_select_color::~pcl_view_select_color()
{
}

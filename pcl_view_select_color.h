//#ifndef PCL_VIEW_SELECT_COLOR_H
//#define PCL_VIEW_SELECT_COLOR_H
#pragma once
#include <QWidget>

class pcl_view_select_color : public QWidget
{
    Q_OBJECT

public:
    pcl_view_select_color();

    ~pcl_view_select_color();

    void setColor(const QColor& c)
    {
        if (c.isValid())
        {
            color = c;
        }

    }
    QColor getColor()
    {
        return color;
    }

private:

    QColor color;


};

// #endif // PCL_VIEW_SELECT_COLOR_H



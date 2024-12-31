//ֻ��һ��
#pragma once
//�������ı���
#pragma execution_character_set("utf-8")
#include <vtkAutoInit.h> 
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
//qtͷ�ļ�
#include <QtWidgets/QMainWindow>
#include <QFileDialog.h>
#include <QMessageBox.h>
#include <QMenuBar>
#include <QMenu>
#include <QToolBar>
#include <QStatusBar>
#include <QFileDialog.h>
#include <QMessageBox.h>
#include <QMenuBar>
#include <QMenu>
#include <QStatusBar>
#include <Qdebug>
#include <QString>
#include <QStandardItemModel>
#include <QProgressBar>
#include <QKeyEvent>

//pcl
#include <pcl/common/common.h>
#include <pcl/point_types.h>		
#include <pcl/point_cloud.h>				
#include <pcl/visualization/pcl_visualizer.h>	
#include <pcl/filters/filter.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h>
#include "pcl_view_select_color.h"
#include <vtkRenderWindow.h>
//ui
#include "ui_QTPCLVTK.h"
#include"view_slider.h"
#include"voxel_filter.h"
#include"zhitong_filter.h"
#include"staticrand_downsamp.h"
#include"filter_function.h"
#include"statistical_denoise.h"
#include"radius_denoise.h"
#include"cur_downsamp.h"
#include"gaussian_filter.h"
#include"doubleline_filter.h"
#include "mls_upsamp.h"
#include"farthest_downsamp.h"
#include"morphological_filter.h"
#include"dpl_sem_seg.h"
//las��osg��gdal
#include "lasreader.hpp"
#include <gdal_priv.h>
#include <osgQT/osgQOpenGLWidget>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Node>
#include <osgDB/ReadFile>
#include <osgTerrain/TerrainTile>
#include <osgTerrain/Layer>
#include <osg/StateSet>
#include <osg/Material>
#include<osg/ComputeBoundsVisitor>
#include <osgGA/TrackballManipulator>
typedef pcl::PointXYZL PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

class QTPCLVTK : public QMainWindow
{
    Q_OBJECT
public:
	QTPCLVTK(QWidget* parent = nullptr);
	~QTPCLVTK();
	//��ǰ�ĵ���
	PointCloudT::Ptr m_Cloud;
	//���ӻ�������
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	//�������꼫ֵ
	PointT p_min, p_max;
	//int p_size;
	double maxLen;
	static int color_set;
	//staticQ StandardItem* itemFolder;
	std::vector<int>cloud_index;
	std::vector<PointCloudT::Ptr>cloud_vec;
	void view_updata(std::vector<PointCloudT::Ptr>& cloud_vec, std::vector<int>& cloud_index);
	double getMinValue(PointT p1, PointT p2);
	double getMaxValue(PointT p1, PointT p2);
	static void CreateCloudFromTxt(const std::string& file_path, pcl::PointCloud<pcl::PointXYZL>::Ptr cloud);
	static void CreateCloudFromlas(const std::string& file_path, pcl::PointCloud<pcl::PointXYZL>::Ptr cloud);
	// ������ص�����
	static void pp_callback_PointsSelect(const pcl::visualization::PointPickingEvent& event, void* args );
private:
	//std::vector<QString>cloud_path;
	Ui::QTPCLVTKClass* ui;
	//osgQOpenGLWidget* osgWidget;
	View_slider* dialog_slider;
	pcl_view_select_color* dialog_colorselect;
signals:

public slots:
	//osg��ʼ��
	void initOSG();
	//osg�������
	void updataosgview();
	//��ģ��
	void on_actionOpendem_triggered();
	//����ѡȡ����
	void on_point_select_clicked();
	//�򿪵���
	void on_actionOpen_triggered();
	//�رյ���
	void on_actionClose_triggered();
	//�������
	void on_actionSave_triggered();
	//���οؼ���ѡ��
	void on_treeView_clicked(const QModelIndex& index);
	//����ͼ
	void on_actionUp_triggered();
	//����ͼ
	void on_actionBottom_triggered();
	//ǰ��ͼ
	void on_actionFront_triggered();
	//����ͼ
	void on_actionBack_triggered();
	//����ͼ
	void on_actionLeft_triggered();
	//����ͼ
	void on_actionRight_triggered();
	//������ɫ����
	void PressBtn_backgroundSelect();
	//�̸߳�ɫ
	void PressBtn_zcolor();
	//���ɫ
	void PressBtn_classcolor();
	//������ɫ
	void PressBtn_ptcolor();
	//���ƴ�С����
	void PressBtn_slider();
	void slider_set(QString data);
	//�����˲�����
	void PressBtn_voxelfilter();
	void voxel_set(QString voxelnum, QString voxelsize);
	//ֱͨ�˲�����
	void PressBtn_zhitongfilter();
	void zhitong_set(int data1,QString data2, QString data3);
	//�̶������������
	void PressBtn_staticrandownsamp();
	void staticrandownsamp_set(QString nums);
	//�����˲�����
	void PressBtn_cursamp();
	void cursamp_set(QString K_search, QString nums);
	//ͳ���˲�����
	void PressBtn_statisticalnoise();
	void statisticalnoise_set(QString d_max, QString nums,int id);
	//�뾶�˲�����
	void PressBtn_radiusdenoise();
	void radiusdenoise_set(QString radius, QString nums);
	//��˹�˲�����
	void PressBtn_gausismooth();
	void gausismooth_set(QString sigma, QString sigmarela, QString threshold, QString radius);
	//˫���˲�����
	void PressBtn_doublelinesmooth();
	void doublelinesmooth_set(QString nums, QString sigmadis, QString sigmanormal);
	//�ƶ���С���˷��ϲ�������
	void PressBtn_mlsupsamp();
	void mlsupsamp_set(QString radius, QString step);
	//��Զ�������
	void PressBtn_farthestdownsamp();
	void farthestdownsamp_set(QString nums);
	void onfarthestdownsampFinished(int mid, std::vector<int>sam_indice);
	//������̬ѧ�˲�
	void PressBtn_morphfilter();
	void morphfilter_set(QString grid_max, QString slope, QString height_init, QString height_max);
	void onmorphfilterFinished(int mid, std::vector<int>sam_indice);
	//pointnet++����ָ�
	void PressBtn_pointnet_sem_seg();
	void dplsemseg_set(QString blocksize, QString blockstride, QString model_path);
	void onpointnet2semseg_Finished(int mid, std::vector<int>label_indice);
};

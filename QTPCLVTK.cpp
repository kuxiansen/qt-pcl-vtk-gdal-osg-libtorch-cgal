#include "QTPCLVTK.h"

int QTPCLVTK::color_set = 0;
QStandardItemModel* model;
QStandardItem* itemFolder; 
QStatusBar* stbar;
std::vector<QString>cloud_path;
osg::ref_ptr<osg::Group>m_pRootGroup = new osg::Group;
osg::ref_ptr<osgViewer::Viewer>pViewer; 
QAction* actptpick;
int actptpick_num = 0;
QTPCLVTK::QTPCLVTK(QWidget *parent)
    : QMainWindow(parent)
{
    ui->setupUi(this);

	//点云初始化
	m_Cloud.reset(new PointCloudT);
	//可视化对象初始化
	viewer.reset(new pcl::visualization::PCLVisualizer("viewer", false));
	//设置VTK可视化窗口指针
	ui->qvtkWidget->SetRenderWindow(viewer->getRenderWindow());
	//设置窗口交互，窗口可接受键盘等事件
	viewer->setupInteractor(ui->qvtkWidget->GetInteractor(), ui->qvtkWidget->GetRenderWindow());
	
	//设置坐标系
	//viewer->addCoordinateSystem(1);
	//工具栏中创建打开文件
	actptpick = new QAction(QIcon(":/img/image/open.png"), tr("单选点云"), this);
	actptpick->setCheckable(true);
	actptpick->setChecked(FALSE);
	ui->toolBar->addAction(actptpick);
	ui->toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);//设置文字图标下
	
	//***********************************connect设置*******************************************//
	ui->behind_view->setIcon(QIcon(":/img/image/back.png"));
	ui->up_view->setIcon(QIcon(":/img/image/up.png"));
	ui->down_view->setIcon(QIcon(":/img/image/bottom.png"));
	ui->front_view->setIcon(QIcon(":/img/image/front.png"));
	ui->left_view->setIcon(QIcon(":/img/image/left.png"));
	ui->right_view->setIcon(QIcon(":/img/image/right.png"));
	ui->open_action->setIcon(QIcon(":/img/image/open.png"));
	stbar = this->statusBar();
	this->setStatusBar(stbar);
	//设置树形控件
	model = new QStandardItemModel(this);
	//设置表头隐藏
	//ui->treeView->setHeaderHidden(true);
	model->setHorizontalHeaderLabels(QStringList() << "点云列表" );
	//设置model 
	ui->treeView->setModel(model);
	//设置展开
	ui->treeView->expandAll();
	ui->progressBar->setVisible(false);
	//初始化osgviewer
	//osgWidget = new osgQOpenGLWidget(ui->qosgWidget);
	//osgWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	connect(ui->opendem_action, SIGNAL(triggered()), this, SLOT(on_actionOpendem_triggered()));//读取dem文件
	connect(ui->open_action, SIGNAL(triggered()), this, SLOT(on_actionOpen_triggered()));//读取文件
	connect(ui->quit_action, SIGNAL(triggered()), this, SLOT(on_actionClose_triggered()));//关闭文件
	connect(ui->save_action, SIGNAL(triggered()), this, SLOT(on_actionSave_triggered()));//保存文件
	connect(model, SIGNAL(clicked(const QModelIndex & index)), this, SLOT(on_treeView_clicked(const QModelIndex & index)));
	//connect(model,SIGNAL(itemChanged(QStandardItem*)),this,SLOT(treeItemChanged(QStandardItem*)));
	connect(ui->qosgWidget, SIGNAL(initialized()), this, SLOT(initOSG()));//osg初始化
	connect(actptpick, SIGNAL(triggered()), this, SLOT(on_point_select_clicked()));//单点选取交互
	connect(ui->up_view, SIGNAL(triggered()), this, SLOT(on_actionUp_triggered()));//俯视图
	connect(ui->down_view, SIGNAL(triggered()), this, SLOT(on_actionBottom_triggered()));//底视图
	connect(ui->front_view, SIGNAL(triggered()), this, SLOT(on_actionFront_triggered()));//前视图
	connect(ui->behind_view, SIGNAL(triggered()), this, SLOT(on_actionBack_triggered()));//后视图
	connect(ui->left_view, SIGNAL(triggered()), this, SLOT(on_actionLeft_triggered()));//左视图
	connect(ui->right_view, SIGNAL(triggered()), this, SLOT(on_actionRight_triggered()));//右视图
	connect(ui->view_background, SIGNAL(triggered(bool)), this, SLOT(PressBtn_backgroundSelect()));//设置背景颜色
	connect(ui->z_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_zcolor()));//高程赋色
	connect(ui->class_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_classcolor()));//类别赋色
	connect(ui->pt_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_ptcolor()));//点云赋色
	connect(ui->pt_size, SIGNAL(triggered(bool)), this, SLOT(PressBtn_slider()));//点云大小
	connect(ui->filtervoxel_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_voxelfilter()));//体素滤波
	connect(ui->filterzhitong_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_zhitongfilter()));//直通滤波
	connect(ui->curdownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_cursamp()));//曲率滤波采样
	connect(ui->staticrandownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_staticrandownsamp()));//固定点采样
	connect(ui->statisticalnoising_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_statisticalnoise()));//统计滤波去噪
	connect(ui->radiusnoising_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_radiusdenoise()));//半径滤波去噪
	connect(ui->gaussfilter_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_gausismooth()));//高斯滤波平滑
	connect(ui->doubleline_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_doublelinesmooth()));//双边滤波平滑
	connect(ui->mlsupsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_mlsupsamp()));//双边滤波平滑
	connect(ui->morphologicalfilter_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_morphfilter()));	//渐进形态学滤波
	connect(ui->farthestdownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_farthestdownsamp()));	//最远点采样法
	connect(ui->s3dis_semsegaction, SIGNAL(triggered(bool)), this, SLOT(PressBtn_pointnet_sem_seg()));	//最远点采样法
	ui->qvtkWidget->update();
}

QTPCLVTK::~QTPCLVTK()
{
	//delete cloud_path;
	/*if (osgWidget != NULL)
	{
		delete osgWidget;
	}*/
	if (dialog_colorselect != NULL)
	{
		delete dialog_colorselect;
	}
	if (dialog_slider != NULL)
	{
		delete dialog_slider;
	}
	delete ui;
}


void QTPCLVTK::on_point_select_clicked()
{
	HKL hCurKL = NULL;
	//强制英文输入法
	hCurKL = GetKeyboardLayout(0);
	LoadKeyboardLayout((LPCSTR)QString("0x0409").utf16(), KLF_ACTIVATE);
	QApplication::processEvents();
	
	//当界面上未选取时，实际传进来的是复选之后的即true
	//if (actptpick->isChecked())
	//{
	//	//设置鼠标交互
	//	
	//	actptpick_num++;
	//}
	//else
	//{
	//	actptpick_num++;
	//	return;
	//}
	if (actptpick_num==0)
	{
		viewer->registerPointPickingCallback(pp_callback_PointsSelect, this);
		actptpick_num++;
	}
	else
	{
		//模拟键盘输入“x”
		INPUT input[2] = {};
		// 设置按下 'X' 键的输入事件
		input[0].type = INPUT_KEYBOARD;
		input[0].ki.wVk = 0;  // 没有虚拟键代码
		input[0].ki.wScan = 0x2D; // 扫描码，'X' 键的扫描码
		input[0].ki.dwFlags = KEYEVENTF_SCANCODE; // 扫描码模式
		input[0].ki.time = 0;
		input[0].ki.dwExtraInfo = 0;

		// 设置松开 'X' 键的输入事件
		input[1].type = INPUT_KEYBOARD;
		input[1].ki.wVk = 0;  // 没有虚拟键代码
		input[1].ki.wScan = 0x2D; // 扫描码，'X' 键的扫描码
		input[1].ki.dwFlags = KEYEVENTF_KEYUP | KEYEVENTF_SCANCODE; // 按键释放
		input[1].ki.time = 0;
		input[1].ki.dwExtraInfo = 0;
		SendInput(2, input, sizeof(INPUT));
	}
}

void QTPCLVTK::pp_callback_PointsSelect(const pcl::visualization::PointPickingEvent& event, void* args)
{
	QTPCLVTK* data = (QTPCLVTK*)args;
	// 如果点击无效，返回
	if (event.getPointIndex() == -1)
		return;
	PointT current_point;
	event.getPoint(current_point.x, current_point.y, current_point.z);
	std::string info = "X: "+std::to_string(current_point.x) +", "+ "Y: " +std::to_string(current_point.y) +", " + "Z: " +std::to_string(current_point.z);
	stbar->clearMessage();
	stbar->showMessage(QString::fromStdString(info));
	//data->ui->qvtkWidget->update();
	//QMessageBox::information(data, "提示", std::to_string(current_point.x).c_str());
	 // 更新点云中选中的点
	//data->clicked_points_3d->points.push_back(current_point);
	// 更新点云显示，绘制为红色
	//pcl::visualization::PointCloudColorHandlerCustom<PointT> red(data->clicked_points_3d, 255, 0, 0);
	//data->viewerPtr->removePointCloud("clicked_points");
	//data->viewerPtr->addPointCloud(data->clicked_points_3d, red, "clicked_points");
	//data->viewerPtr->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "clicked_points");
}

void QTPCLVTK::initOSG()
{
	pViewer = ui->qosgWidget->getOsgViewer();
	pViewer->getCamera()->setClearColor(osg::Vec4(0.0f, 0.0f, 0.0f, 1.0f));  // 设置清除颜色为黑色
	osg::ref_ptr<osgGA::CameraManipulator>pManipulator = new osgGA::TrackballManipulator();
	pViewer->setCameraManipulator(pManipulator);
}

void QTPCLVTK::updataosgview()
{
	// 计算新的场景外包
	osg::BoundingSphere boundingSphere;
	osg::ComputeBoundsVisitor cbVisitor;
	m_pRootGroup->accept(cbVisitor);
	osg::BoundingBox& bb = cbVisitor.getBoundingBox();
	// 更新相机位置,俯视角度观察模型
	double radius = osg::maximum(double(boundingSphere.radius()), 1e-6);
	pViewer->setCameraManipulator(NULL);
	//osgGA::CameraManipulator* pManipulator = new osgGA::TrackballManipulator();
	pViewer->setCameraManipulator(new osgGA::TrackballManipulator);
	//根据分辨率确定合适的投影来保证显示的图形不变形
	double fovy, aspectRatio, zNear, zFar;
	pViewer->getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, zFar);
	double newAspectRatio = double(ui->qosgWidget->width()) / double(ui->qosgWidget->height());
	double aspectRatioChange = newAspectRatio / aspectRatio;
	if (aspectRatioChange != 1.0)
	{
		//设置投影矩阵
		pViewer->getCamera()->getProjectionMatrix() *= osg::Matrix::scale(1.0 / aspectRatioChange, 1.0, 1.0);
	}
}

void QTPCLVTK::CreateCloudFromTxt(const std::string& file_path, pcl::PointCloud<pcl::PointXYZL>::Ptr cloud)
{
	std::ifstream file(file_path.c_str());
	std::string line;
	pcl::PointXYZL point;
	//float nx, ny, nz;
	while (getline(file, line)) {
		std::stringstream ss(line);
		ss >> point.x;
		ss >> point.y;
		ss >> point.z;
		ss >> point.label;
		//ss >> nx;
		//ss >> ny;
		//ss >> nz;
		cloud->push_back(point);
	}
	file.close();
}

void QTPCLVTK::CreateCloudFromlas(const std::string& file_path, pcl::PointCloud<pcl::PointXYZL>::Ptr cloud)
{
	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(file_path.c_str());
	LASreader* lasreader = lasreadopener.open();
	pcl::PointXYZL point;
	while (lasreader->read_point())
	{
		point.x = lasreader->point.get_x();
		point.y = lasreader->point.get_y();
		point.z = lasreader->point.get_z();
		point.label = lasreader->point.classification;
		cloud->push_back(point);
	}
	lasreader->close();
	delete lasreader;
}
//打开dem
void QTPCLVTK::on_actionOpendem_triggered()
{
	// 获取dem路径
	QString path = QFileDialog::getOpenFileName(this, "选择dem", ".//", "dem文件(*.osg *.tif );;所有文件(*.*)");
	if (path.isEmpty())
		return;
	else
	{
		// 初始化 GDAL
		GDALAllRegister();
		// 打开 DEM 文件
		GDALDataset* poDataset = (GDALDataset*)GDALOpen(path.toStdString().c_str(), GA_ReadOnly);
		if (!poDataset) 
		{
			qDebug() << "读取dem文件失败" << endl;
			return;
		}
		// 获取地理变换信息
		double gdalGeoTransform[6];
		if (poDataset->GetGeoTransform(gdalGeoTransform) != CE_None) 
		{
			qDebug() << "读取dem文件GeoTransform信息失败" << endl;
			GDALClose(poDataset);
			return;
		}
		// 创建 osg::HeightField 对象
		osg::ref_ptr<osg::HeightField> hf = new osg::HeightField();
		hf->allocate(poDataset->GetRasterXSize(), poDataset->GetRasterYSize());
		hf->setOrigin(osg::Vec3(gdalGeoTransform[0], gdalGeoTransform[3], 0));
		hf->setXInterval(gdalGeoTransform[1]);
		hf->setYInterval(gdalGeoTransform[5]);
		// 读取 DEM 高程数据
		int rasterXSize = poDataset->GetRasterXSize();
		int rasterYSize = poDataset->GetRasterYSize();
		float* heightData = new float[rasterXSize * rasterYSize];
		// 获取栅格数据
		if (poDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, rasterXSize, rasterYSize, heightData, rasterXSize, rasterYSize, GDT_Float32, 0, 0) != CE_None) 
		{
			qDebug() << "读取dem文件栅格数据失败" << endl;
			delete[] heightData;
			GDALClose(poDataset);
			return ;
		}
		// 设置 NoData 值
		float noDataValue = poDataset->GetRasterBand(1)->GetNoDataValue();
		float noDataValueFill = 0.0f;
		// 将数据填充到 HeightField
		float* heightPtr = heightData;
		for (int r = rasterYSize - 1; r >= 0; --r) 
		{
			for (int c = 0; c < rasterXSize; ++c) 
			{
				float h = *heightPtr++;
				if (h != noDataValue)
					hf->setHeight(c, r, h);
				else
					hf->setHeight(c, r, noDataValueFill); // 填充 NoData 区域
			}
		}
		// 创建 Locator 和 TerrainTile
		osg::ref_ptr<osgTerrain::Locator> locator = new osgTerrain::Locator();
		double minX = std::min(gdalGeoTransform[0], gdalGeoTransform[0] + rasterXSize * gdalGeoTransform[1]);
		double minY = std::min(gdalGeoTransform[3], gdalGeoTransform[3] + rasterYSize * gdalGeoTransform[5]);
		double maxX = std::max(gdalGeoTransform[0], gdalGeoTransform[0] + rasterXSize * gdalGeoTransform[1]);
		double maxY = std::max(gdalGeoTransform[3], gdalGeoTransform[3] + rasterYSize * gdalGeoTransform[5]);

		locator->setTransformAsExtents(minX, minY, maxX, maxY);

		// 创建 HeightFieldLayer 并与 TerrainTile 绑定
		osg::ref_ptr<osgTerrain::HeightFieldLayer> hfl = new osgTerrain::HeightFieldLayer();
		hfl->setHeightField(hf);
		hfl->setLocator(locator.get());

		// 创建 ImageLayer 并计算颜色
		osg::ref_ptr<osgTerrain::ImageLayer> colorLayer = new osgTerrain::ImageLayer();
		osg::ref_ptr<osg::Image> colorImage = new osg::Image();
		colorImage->allocateImage(rasterXSize, rasterYSize, 1, GL_RGBA, GL_UNSIGNED_BYTE);

		// 获取最小和最大高度值用于颜色映射
		float minHeight = std::numeric_limits<float>::max();
		float maxHeight = std::numeric_limits<float>::lowest();
		for (int r = 0; r < rasterYSize; ++r) {
			for (int c = 0; c < rasterXSize; ++c) {
				float h = hf->getHeight(c, r);
				if (h != noDataValueFill) {
					if (h < minHeight) minHeight = h;
					if (h > maxHeight) maxHeight = h;
				}
			}
		}

		// 填充颜色图像
		for (int r = 0; r < rasterYSize; ++r) {
			for (int c = 0; c < rasterXSize; ++c) {
				float h = hf->getHeight(c, r);
				if (h != noDataValueFill) {
					float normalizedHeight = (h - minHeight) / (maxHeight - minHeight);
					unsigned char rColor = static_cast<unsigned char>((1.0f - normalizedHeight) * 255);
					unsigned char gColor = static_cast<unsigned char>(normalizedHeight * 255);
					unsigned char bColor = 0;
					unsigned char aColor = 255; // alpha通道

					// 获取图像的像素数据
					unsigned char* pixel = colorImage->data() + (r * rasterXSize + c) * 4;  // 每个像素占4个字节（RGBA）
					pixel[0] = rColor;  // 红色
					pixel[1] = gColor;  // 绿色
					pixel[2] = bColor;  // 蓝色
					pixel[3] = aColor;  // Alpha通道
				}
				else {
					// 设置NoData区域为白色
					unsigned char* pixel = colorImage->data() + (r * rasterXSize + c) * 4;
					pixel[0] = 0;  // 红色
					pixel[1] = 0;  // 绿色
					pixel[2] = 0;  // 蓝色
					pixel[3] = 0;  // Alpha通道
				}
			}
		}


		colorLayer->setImage(colorImage.get());

		osg::ref_ptr<osgTerrain::TerrainTile> terrainTile = new osgTerrain::TerrainTile();
		terrainTile->setElevationLayer(hfl);
		terrainTile->setColorLayer(1, colorLayer.get());
		// 创建 StateSet，用于设置渲染状态（材质）
		osg::ref_ptr<osg::StateSet> stateSet = terrainTile->getOrCreateStateSet();
		// 创建材质对象并将其应用到 StateSet
		osg::ref_ptr<osg::Material> material = new osg::Material();
		material->setDiffuse(osg::Material::FRONT_AND_BACK, osg::Vec4(0.5f, 0.5f, 0.5f, 1.0f));  // 设置漫反射颜色
		material->setSpecular(osg::Material::FRONT_AND_BACK, osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));  // 设置镜面反射颜色
		material->setShininess(osg::Material::FRONT_AND_BACK, 50.0f);  // 设置光泽度
		stateSet->setAttributeAndModes(material.get(), osg::StateAttribute::ON);  // 将材质应用到 StateSet
		m_pRootGroup->addChild(terrainTile.get());
		pViewer->setSceneData(m_pRootGroup);
		pViewer->addEventHandler(new osgViewer::WindowSizeHandler());
		// 释放内存
		delete[] heightData;
		GDALClose(poDataset);
		updataosgview();
	}
}
//打开点云
void QTPCLVTK::on_actionOpen_triggered()
{

	//获取点云路径
	QString path = QFileDialog::getOpenFileName(this, "选择点云", ".//", "点云文件(*.ply *.pcd *.txt *.las);;所有文件(*.*)");
	
	pcl::PointCloud<pcl::PointXYZL>::Ptr cloud_tmp(new pcl::PointCloud<pcl::PointXYZL>);

	if (path.isEmpty())
		return;

	if (path.endsWith(".pcd", Qt::CaseInsensitive))
	{
		qDebug() << path;
		if (pcl::io::loadPCDFile(path.toStdString(), *cloud_tmp) == -1)
		{
			qDebug()<< "读取pcd点云失败" << endl;
			return;
		}
		else
		{
			QMessageBox::information(this, "提示", "pcd点云读取完毕");
		}
	}
	else if (path.endsWith(".ply", Qt::CaseInsensitive))
	{
		qDebug() << path;
		if (pcl::io::loadPLYFile(path.toStdString(), *cloud_tmp) == -1)
		{
			qDebug() << "读取ply点云失败" << endl;
			return;
		}
		else
		{
			QMessageBox::information(this, "提示", "ply点云读取完毕");
		}
	}
	else if (path.endsWith(".txt", Qt::CaseInsensitive))
	{
		qDebug() << path;
		CreateCloudFromTxt(path.toStdString(), cloud_tmp);
		QMessageBox::information(this, "提示", "txt点云读取完毕");
	}
	else if (path.endsWith(".las", Qt::CaseInsensitive))
	{
		qDebug() << path;
		CreateCloudFromlas(path.toStdString(), cloud_tmp);
		QMessageBox::information(this, "提示", "las点云读取完毕");
	}
	else
	{
		QMessageBox::warning(this, "Warning", "点云读取格式错误！");
	}


	//清空点云
	//m_Cloud->clear();
	//viewer->removeAllPointClouds();
	//viewer->removeAllCoordinateSystems();

	if (cloud_tmp->is_dense)
		pcl::copyPointCloud(*cloud_tmp, *m_Cloud);
	else
	{
		PCL_WARN("Cloud is not dense! Non finite points will be removed\n");
		std::vector<int> vec;
		pcl::removeNaNFromPointCloud(*cloud_tmp, *m_Cloud, vec);
	}
	
	cloud_vec.push_back(m_Cloud->makeShared());
	cloud_path.push_back(path);
	cloud_index.push_back(1);
	//cloud_tmp->clear();
	//m_Cloud->clear();
	if (cloud_vec.size()>=1)
	{
		//QStandardItem* itemFolder = new QStandardItem(QStringLiteral("cloud%1").arg(cloud_vec.size() - 1));
		itemFolder = new QStandardItem(path);
		itemFolder->setCheckable(true);
		itemFolder->setCheckState(Qt::Checked);//获取选中状态
		model->appendRow(itemFolder);
	}
	view_updata(cloud_vec, cloud_index);
	//PointT max_p;
	//PointT min_p;
	//int nums = 0;
	//for (size_t i = 0; i < cloud_vec.size(); i++)
	//{
	//	//添加到窗口
	//	
	//	switch (color_set)
	//	{
	//	case 0:
	//	{
	//		pcl::visualization::PointCloudColorHandlerGenericField<PointT>fildcolor(cloud_vec[i], "label");
	//		viewer->addPointCloud(cloud_vec[i], fildcolor, "Number" + std::to_string(i));
	//		break;
	//	}
	//		
	//	case 1:
	//	{
	//		pcl::visualization::PointCloudColorHandlerGenericField<PointT>fildcolor(cloud_vec[i], "z");
	//		viewer->addPointCloud(cloud_vec[i], fildcolor, "Number" + std::to_string(i));
	//		break;
	//	}
	//	}
	//	
	//	pcl::getMinMax3D(*cloud_vec[i], min_p, max_p);
	//	if (i==0)
	//	{
	//		p_max = max_p;
	//		p_min = min_p;
	//	}
	//	else
	//	{
	//		if (min_p.x< p_min.x)
	//		{
	//			p_min.x = min_p.x;
	//		}
	//		if (min_p.y < p_min.y)
	//		{
	//			p_min.y = min_p.y;
	//		}
	//		if (min_p.z < p_min.z)
	//		{
	//			p_min.z = min_p.z;
	//		}
	//		if (max_p.x > p_max.x)
	//		{
	//			p_max.x = max_p.x;
	//		}
	//		if (max_p.y > p_max.y)
	//		{
	//			p_max.y = max_p.y;
	//		}
	//		if (max_p.z > p_max.z)
	//		{
	//			p_max.z = max_p.z;
	//		}
	//	}
	//	nums += cloud_vec[i]->size();
	//}
	//maxLen = getMaxValue(p_max, p_min);
	////重设视角
	//viewer->resetCamera();
	//std::string num = std::to_string(nums);
	//ui->canshu->setText("点云数量为:"+ QString::fromStdString(num));
	//ui->canshu->append("点云x轴范围:" + QString::fromStdString(std::to_string(p_min.x))+" -" + QString::fromStdString(std::to_string(p_max.x)));
	//ui->canshu->append("点云y轴范围:" + QString::fromStdString(std::to_string(p_min.y)) + " -" + QString::fromStdString(std::to_string(p_max.y)));
	//ui->canshu->append("点云z轴范围:" + QString::fromStdString(std::to_string(p_min.z)) + " -" + QString::fromStdString(std::to_string(p_max.z)));
	////刷新窗口
	//ui->qvtkWidget->update();
}
//关闭点云
void QTPCLVTK::on_actionClose_triggered()
{
	if (cloud_vec.size() != 0)
	{
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 0)
			{
				cloud_vec[i]->clear();
				cloud_path[i].clear();
			}
		}
		auto it_vec = cloud_vec.begin();
		auto it_index = cloud_index.begin();
		auto it_path = cloud_path.begin();
		while (it_vec!= cloud_vec.end())
		{
			if ((*it_vec)->empty())
			{
				it_vec = cloud_vec.erase(it_vec);
			}
			else
			{
				it_vec++;
			}
		}
		while (it_index != cloud_index.end())
		{
			if ((*it_index)==0)
			{
				it_index = cloud_index.erase(it_index);
			}
			else
			{
				it_index++;
			}
		}
		while (it_path != cloud_path.end())
		{
			if ((*it_path).isEmpty())
			{
				it_path = cloud_path.erase(it_path);
			}
			else
			{
				it_path++;
			}
		}
		model->clear();
		model->setHorizontalHeaderLabels(QStringList() << "点云列表");
		if (cloud_vec.size() >= 1)
		{
			for (size_t i = 0; i < cloud_vec.size(); i++)
			{
				itemFolder = new QStandardItem(cloud_path[i]);
				itemFolder->setCheckable(true);
				itemFolder->setCheckState(Qt::Checked);//获取选中状态
				model->appendRow(itemFolder);
			}
			//QStandardItem* itemFolder = new QStandardItem(QStringLiteral("cloud%1").arg(cloud_vec.size() - 1)
		}
		view_updata(cloud_vec, cloud_index);
	}
	else
	{
		QMessageBox::warning(this, "Warning", "视图中没有点云！");
		return;
	}
}
//保存点云
void QTPCLVTK::on_actionSave_triggered()
{
	if (cloud_vec.size() != 0)
	{
		QString path = QFileDialog::getSaveFileName(this, "选择点云", ".//", "点云文件(*.ply *.pcd);;所有文件(*.*)");
		if (path.isEmpty()) return;
		if (path.endsWith(".pcd", Qt::CaseInsensitive))
		{
			PointCloudT::Ptr cloud_merged(new PointCloudT);
			for (size_t i = 0; i < cloud_vec.size(); i++)
			{
				if (cloud_index[i] == 1)
				{
					*cloud_merged += *cloud_vec[i];
				}
			}
			pcl::io::savePCDFileBinary(path.toStdString(), *cloud_merged);
			QMessageBox::information(this, "提示", "pcd点云保存成功");
		}
		else if (path.endsWith(".ply", Qt::CaseInsensitive))
		{
			PointCloudT::Ptr cloud_merged(new PointCloudT);
			for (size_t i = 0; i < cloud_vec.size(); i++)
			{
				if (cloud_index[i] == 1)
				{
					*cloud_merged += *cloud_vec[i];
				}
			}
			pcl::io::savePLYFileBinary(path.toStdString(), *cloud_merged);
			QMessageBox::information(this, "提示", "ply点云保存成功");
		}
		else
		{
			PointCloudT::Ptr cloud_merged(new PointCloudT);
			for (size_t i = 0; i < cloud_vec.size(); i++)
			{
				if (cloud_index[i] == 1)
				{
					*cloud_merged += *cloud_vec[i];
				}
			}
			path.append(".ply");
			pcl::io::savePLYFileBinary(path.toStdString(), *cloud_merged);
			QMessageBox::information(this, "提示", "ply点云保存成功");
		}
	}
	else
	{
		QMessageBox::warning(this, "Warning", "视图中没有点云！");
	}
}



void QTPCLVTK::on_treeView_clicked(const QModelIndex& index)
{
	QStandardItem* item = model->itemFromIndex(index);
	if (item ==nullptr)
	{
		return;
	}
	if (item->isCheckable())
	{
		Qt::CheckState state = item->checkState();//获取当前选择状态
		if (state==Qt::Checked)
		{
			cloud_index[index.row()] = 1;
		}
		if (state == Qt::Unchecked)
		{
			cloud_index[index.row()] = 0;
		}
		view_updata(cloud_vec, cloud_index);
	}
}

void QTPCLVTK::view_updata(std::vector<PointCloudT::Ptr>& cloud_vec, std::vector<int>& cloud_index)
{
	this->viewer->removeAllPointClouds();
	int nums = 0;
	PointT max_p;
	PointT min_p;
	for (size_t i = 0; i < cloud_vec.size(); i++)
	{
		//添加到窗口
		if (cloud_index[i] == 1)
		{
			switch (color_set)
			{
			case 0:
			{
				pcl::visualization::PointCloudColorHandlerGenericField<PointT>fildcolor(cloud_vec[i], "label");
				viewer->addPointCloud(cloud_vec[i], fildcolor, "Number" + std::to_string(i));
				break;
			}

			case 1:
			{
				pcl::visualization::PointCloudColorHandlerGenericField<PointT>fildcolor(cloud_vec[i], "z");
				viewer->addPointCloud(cloud_vec[i], fildcolor, "Number" + std::to_string(i));
				break;
			}
			}
			pcl::getMinMax3D(*cloud_vec[i], min_p, max_p);
			if (i == 0)
			{
				p_max = max_p;
				p_min = min_p;
			}
			else
			{
				if (min_p.x < p_min.x)
				{
					p_min.x = min_p.x;
				}
				if (min_p.y < p_min.y)
				{
					p_min.y = min_p.y;
				}
				if (min_p.z < p_min.z)
				{
					p_min.z = min_p.z;
				}
				if (max_p.x > p_max.x)
				{
					p_max.x = max_p.x;
				}
				if (max_p.y > p_max.y)
				{
					p_max.y = max_p.y;
				}
				if (max_p.z > p_max.z)
				{
					p_max.z = max_p.z;
				}
			}
			nums += cloud_vec[i]->size();
		}
	}
	std::string num = std::to_string(nums);
	maxLen = getMaxValue(p_max, p_min);
	viewer->resetCamera();
	ui->canshu->setText("点云数量为:" + QString::fromStdString(num));
	ui->canshu->append("点云x轴范围:" + QString::fromStdString(std::to_string(p_min.x)) + " -" + QString::fromStdString(std::to_string(p_max.x)));
	ui->canshu->append("点云y轴范围:" + QString::fromStdString(std::to_string(p_min.y)) + " -" + QString::fromStdString(std::to_string(p_max.y)));
	ui->canshu->append("点云z轴范围:" + QString::fromStdString(std::to_string(p_min.z)) + " -" + QString::fromStdString(std::to_string(p_max.z)));
	//ui->canshu->append("点云类别范围:" + QString::fromStdString(std::to_string(p_min.label)) + " -" + QString::fromStdString(std::to_string(p_max.label)));
	ui->qvtkWidget->update();
}




//俯视图
void QTPCLVTK::on_actionUp_triggered()
{
	if (cloud_vec.size()!=0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_max.z + 2 * maxLen, 0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_max.z, 0, 1, 0);
		ui->qvtkWidget->update();
	}
}
//底视图
void QTPCLVTK::on_actionBottom_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_min.z - 2 * maxLen, 0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_min.z, 0, 1, 0);
		ui->qvtkWidget->update();
	}
}
//前视图
void QTPCLVTK::on_actionFront_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), p_min.y - 2 * maxLen, 0.5 * (p_min.z + p_max.z), 0.5 * (p_min.x + p_max.x), p_min.y, 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//后视图
void QTPCLVTK::on_actionBack_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), p_max.y + 2 * maxLen, 0.5 * (p_min.z + p_max.z), 0.5 * (p_min.x + p_max.x), p_min.y, 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//左视图
void QTPCLVTK::on_actionLeft_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(p_min.x - 2 * maxLen, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), p_max.x, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//右视图
void QTPCLVTK::on_actionRight_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(p_max.x + 2 * maxLen, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), p_max.x, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}

double QTPCLVTK::getMinValue(PointT p1, PointT p2)
{
	double min = 0;

	if (p1.x - p2.x > p1.y - p2.y)
	{
		min = p1.y - p2.y;
	}
	else
	{
		min = p1.x - p2.x;
	}

	if (min > p1.z - p2.z)
	{
		min = p1.z - p2.z;
	}

	return min;
}
double QTPCLVTK::getMaxValue(PointT p1, PointT p2)
{
	double max = 0;

	if (p1.x - p2.x > p1.y - p2.y)
	{
		max = p1.x - p2.x;

	}
	else
	{
		max = p1.y - p2.y;
	}

	if (max < p1.z - p2.z)
	{
		max = p1.z - p2.z;
	}

	return max;
}

//设置背景颜色
void QTPCLVTK::PressBtn_backgroundSelect()
{
	dialog_colorselect = new pcl_view_select_color();
	QColor color = dialog_colorselect->getColor();
	viewer->setBackgroundColor(color.redF(), color.greenF(), color.blueF());
	return;
}
//高程赋色
void QTPCLVTK::PressBtn_zcolor()
{
	if (cloud_vec.size() != 0)
	{
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i]==1)
			{
				pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZL>render(cloud_vec[i], "z");
				viewer->updatePointCloud(cloud_vec[i], render, "Number" + std::to_string(i));
			}
		}
		color_set = 1;
		ui->qvtkWidget->update();
	}
	else
	{
		QMessageBox::warning(this, "Warning", "视图中没有点云！");
	}
}
void QTPCLVTK::PressBtn_classcolor()
{
	if (cloud_vec.size() != 0)
	{
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZL>render(cloud_vec[i], "label");
				viewer->updatePointCloud(cloud_vec[i], render, "Number" + std::to_string(i));
			}
		}
		color_set = 0;
		ui->qvtkWidget->update();
	}
	else
	{
		QMessageBox::warning(this, "Warning", "视图中没有点云！");
	}
}
//设置点云颜色
void QTPCLVTK::PressBtn_ptcolor()
{
	if (cloud_vec.size() != 0)
	{
		dialog_colorselect = new pcl_view_select_color;
		QColor color = dialog_colorselect->getColor();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZL>selct_color(cloud_vec[i], color.redF() * 255, color.greenF() * 255, color.blueF() * 255);
				viewer->updatePointCloud(cloud_vec[i], selct_color,"Number" + std::to_string(i));
			}
		}
		
		ui->qvtkWidget->update();
	}
	else 
	{
		QMessageBox::warning(this, "Warning", "视图中没有点云！");
	}
}

//点的大小设置
void QTPCLVTK::PressBtn_slider()
{
	dialog_slider = new View_slider(this);
	connect(dialog_slider, SIGNAL(senddata(QString)), this, SLOT(slider_set(QString)));
	dialog_slider->exec();
	/*if (dialog_slider->exec() == QDialog::Accepted) 
	{
	}
	else 
	{
		
	}*/
	delete dialog_slider;
}
void QTPCLVTK::slider_set(QString data)
{

	int p_size = data.toInt();
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, p_size);


}
//体素滤波窗口
void QTPCLVTK::PressBtn_voxelfilter()
{
	voxel_filter* dialog_voxelfilter = new voxel_filter(this);
	connect(dialog_voxelfilter, SIGNAL(send_voxel(QString, QString)), this, SLOT(voxel_set(QString, QString)));
	dialog_voxelfilter->exec();
}
void QTPCLVTK::voxel_set(QString voxelnum, QString voxelsize)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		float size = voxelsize.toFloat();
		int minnum = voxelnum.toUInt();
		//QMessageBox::warning(this, "warning", voxelnum);
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_filter_voxel(cloud_vec[i]->makeShared(), size, minnum);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//直通滤波窗口
void QTPCLVTK::PressBtn_zhitongfilter()
{
	zhitong_filter* dialog_zhitongfilter = new zhitong_filter(this);
	connect(dialog_zhitongfilter, SIGNAL(sendcanshu(int, QString, QString)), this, SLOT(zhitong_set(int, QString, QString)));
	dialog_zhitongfilter->exec();
}
void QTPCLVTK::zhitong_set(int data1, QString data2, QString data3)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		double minsize = data2.toFloat();
		double maxsize = data3.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_filter_zhitong(cloud_vec[i]->makeShared(), data1, minsize, maxsize);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//固定随机采样窗口
void QTPCLVTK::PressBtn_staticrandownsamp()
{
	staticrandownsamp* dialog_staticrandownsamp = new staticrandownsamp(this);
	connect(dialog_staticrandownsamp, SIGNAL(send_nums(QString)), this, SLOT(staticrandownsamp_set(QString)));
	dialog_staticrandownsamp->show();
}
void QTPCLVTK::staticrandownsamp_set(QString nums)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int num = nums.toInt();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_static_randomsamp(cloud_vec[i]->makeShared(), num);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//曲率采样窗口
void QTPCLVTK::PressBtn_cursamp()
{
	cur_downsamp* dialog_curdownsamp = new cur_downsamp(this);
	connect(dialog_curdownsamp,SIGNAL(send_cursamp(QString, QString)),this,SLOT(cursamp_set(QString, QString)));
	dialog_curdownsamp->exec();
}
void QTPCLVTK::cursamp_set(QString K_search, QString nums)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int K = K_search.toInt();
		float num = nums.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_cursamp(cloud_vec[i]->makeShared(), K, num);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
				
		}
		view_updata(cloud_vec, cloud_index);
	}
}


//直通滤波去噪窗口
void QTPCLVTK::PressBtn_statisticalnoise()
{
	statiscal_denoise* dialog_statiscaldenoise = new statiscal_denoise(this);
	connect(dialog_statiscaldenoise, SIGNAL(send_std(QString, QString,int)), this, SLOT(statisticalnoise_set(QString, QString,int)));
	dialog_statiscaldenoise->exec();
}
void QTPCLVTK::statisticalnoise_set(QString d_max, QString nums,int id)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int num = nums.toInt();
		double stddev = d_max.toDouble();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				if (id == 1)
				{
					auto cloudout = pcl_statiscal_denoise(cloud_vec[i]->makeShared(), stddev, num);
					if (cloudout->empty())
					{
						QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
					}
					else {
						cloud_vec[i] = std::move(cloudout);
					}
				}
				if (id==2)
				{
					auto cloudout = my_statiscal_denoise(cloud_vec[i]->makeShared(), stddev, num);
					if (cloudout->empty())
					{
						QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
					}
					else {
						cloud_vec[i] = std::move(cloudout);
					}
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//半径滤波窗口
void QTPCLVTK::PressBtn_radiusdenoise()
{
	radius_denoise* dialog_radiusdenoise = new radius_denoise(this);
	connect(dialog_radiusdenoise, SIGNAL(send_radius_nois(QString, QString)), this, SLOT(radiusdenoise_set(QString, QString)));
	dialog_radiusdenoise->exec();
}
void QTPCLVTK::radiusdenoise_set(QString radius, QString nums)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		float s_radius = radius.toFloat();
		int num = nums.toInt();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_radiussamp(cloud_vec[i]->makeShared(), s_radius, num);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//高斯滤波窗口
void QTPCLVTK::PressBtn_gausismooth()
{
	gaussian_filter* dialog_gaussianfilter = new gaussian_filter(this);
	connect(dialog_gaussianfilter, SIGNAL(send_gaussian(QString, QString, QString, QString)), this, SLOT(gausismooth_set(QString, QString, QString, QString)));
	dialog_gaussianfilter->exec();
}
void QTPCLVTK::gausismooth_set(QString sigma, QString sigmarela, QString threshold, QString radius)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		float s_sigma = sigma.toFloat();
		float s_sigmarela = sigmarela.toFloat();
		float s_threshold = threshold.toFloat();
		float s_radius = radius.toFloat();
		
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_gausifilter(cloud_vec[i]->makeShared(), s_sigma, s_sigmarela, s_threshold,s_radius);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//双边滤波窗口
void QTPCLVTK::PressBtn_doublelinesmooth()
{
	doubleline_filter* dialog_doublelinefilter = new doubleline_filter(this);
	connect(dialog_doublelinefilter, SIGNAL(send_doubleline(QString, QString, QString)), this, SLOT(doublelinesmooth_set(QString , QString , QString)));
	dialog_doublelinefilter->exec();
}
void QTPCLVTK::doublelinesmooth_set(QString nums, QString sigmadis, QString sigmanormal)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int s_num = nums.toInt();
		float s_sigmadis = sigmadis.toFloat();
		float s_sigmanormal = sigmanormal.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = my_doublelinefilter(cloud_vec[i]->makeShared(), s_num, s_sigmadis, s_sigmanormal);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}
//mls上采样窗口
void QTPCLVTK::PressBtn_mlsupsamp()
{
	mls_upsamp* dialog_mlsupsamp = new mls_upsamp(this);
	connect(dialog_mlsupsamp, SIGNAL(send_mlsup(QString, QString)), this, SLOT(mlsupsamp_set(QString, QString)));
	dialog_mlsupsamp->exec();
}

void QTPCLVTK::mlsupsamp_set(QString radius, QString step)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		float s_radius = radius.toFloat();
		float s_step = step.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				auto cloudout = pcl_mlsupsamp(cloud_vec[i]->makeShared(), s_radius, s_step);
				if (cloudout->empty())
				{
					QMessageBox::warning(this, "warning", "点云输出为空或有误,请重新输入参数");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}
//最远点采样法
void QTPCLVTK::PressBtn_farthestdownsamp()
{
	farthest_downsamp* dialog_farthestdownsamp = new farthest_downsamp(this);
	connect(dialog_farthestdownsamp, SIGNAL(send_farthestsamp(QString)), this, SLOT(farthestdownsamp_set(QString)));
	dialog_farthestdownsamp->exec();
}

void QTPCLVTK::farthestdownsamp_set(QString nums)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int s_nums = nums.toInt();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			
			if (cloud_index[i] == 1)
			{
				farest_sampling* funcer = new farest_sampling(i,cloud_vec[i]->makeShared(), s_nums);
				QThread* thread = new QThread;
				ui->progressBar->setRange(0, 100);  // 设置进度条范围
				ui->progressBar->setValue(0);  // 设置初始值
				ui->progressBar->setVisible(true);
				funcer->moveToThread(thread);
				QMessageBox::about(this, "about", "任务转移至子线程");
				connect(thread, &QThread::started, funcer, &farest_sampling::sampling);
				connect(funcer, &farest_sampling::finished, thread, &QThread::quit);
				connect(funcer, &farest_sampling::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(funcer, &farest_sampling::finished, this, &QTPCLVTK::onfarthestdownsampFinished,Qt::QueuedConnection); // 连接带返回值的完成信号	
				connect(funcer, &farest_sampling::finished, funcer, &QObject::deleteLater);
				connect(thread, &QThread::finished, thread, &QThread::deleteLater);
				thread->start();
			}
		}
	}
}


void QTPCLVTK::onfarthestdownsampFinished(int mid, std::vector<int>sam_indice)
{
	ui->progressBar->setVisible(false);
	if (sam_indice.empty())
	{
		QMessageBox::warning(this, "warning", "数值过大，返回原始点云");
	}
	else
	{
		PointCloudT::Ptr cloud_filter(new PointCloudT);
		pcl::copyPointCloud(*cloud_vec[mid], sam_indice, *cloud_filter);
		cloud_vec[mid] = cloud_filter;
		view_updata(cloud_vec, cloud_index);
	}
	//auto cloudout = cloud_final;
	//cloud_vec[mid] = cloudout;
	//view_updata(cloud_vec, cloud_index);
}
//渐进形态学滤波
void QTPCLVTK::PressBtn_morphfilter()
{
	morphological_filter* dialog_morphfilter = new morphological_filter(this);
	connect(dialog_morphfilter, SIGNAL(send_morph(QString, QString, QString, QString)), this, SLOT(morphfilter_set(QString, QString, QString, QString)));
	dialog_morphfilter->exec();
}

void QTPCLVTK::morphfilter_set(QString grid_max, QString slope, QString height_init, QString height_max)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		int max_grid = grid_max.toInt();
		float my_slope = slope.toFloat();
		float init_height = height_init.toFloat();
		float max_height = height_max.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				morph_filter* morph = new morph_filter(i, cloud_vec[i]->makeShared(), max_grid, my_slope, init_height, max_height);
				QThread* thread_m = new QThread;
				ui->progressBar->setRange(0, 100);  // 设置进度条范围
				ui->progressBar->setValue(0);  // 设置初始值
				ui->progressBar->setVisible(true);
				morph->moveToThread(thread_m);
				QMessageBox::about(this, "about", "任务转移至子线程");
				connect(thread_m, &QThread::started, morph, &morph_filter::my_filtering);
				connect(morph, &morph_filter::output, thread_m, &QThread::quit);
				connect(morph, &morph_filter::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(morph, &morph_filter::output, this, &QTPCLVTK::onmorphfilterFinished, Qt::QueuedConnection); // 连接带返回值的完成信号	
				connect(morph, &morph_filter::output, morph, &QObject::deleteLater);
				connect(thread_m, &QThread::finished, thread_m, &QThread::deleteLater);
				thread_m->start();
			}
		}
	}
}

void QTPCLVTK::onmorphfilterFinished(int mid, std::vector<int> sam_indice)
{
	ui->progressBar->setVisible(false);
	if (sam_indice.empty())
	{
		QMessageBox::warning(this, "warning", "数值过大，返回原始点云");
	}
	else
	{
		//QMessageBox::warning(this, "warning", "点云输入");
		PointCloudT::Ptr cloud_filter(new PointCloudT);
		pcl::copyPointCloud(*cloud_vec[mid], sam_indice, *cloud_filter);
		cloud_vec[mid] = cloud_filter;
		view_updata(cloud_vec, cloud_index);
	}
}

////pointnet++室内语义分割
void QTPCLVTK::PressBtn_pointnet_sem_seg()
{
	dpl_sem_seg* dialog_dplsemseg = new dpl_sem_seg(this);
	connect(dialog_dplsemseg, SIGNAL(send_block(QString, QString, QString)), this, SLOT(dplsemseg_set(QString, QString, QString)));
	dialog_dplsemseg->exec();
}

void QTPCLVTK::dplsemseg_set(QString blocksize, QString blockstride, QString model_path)
{
	if (cloud_vec.empty())
	{
		QMessageBox::warning(this, "warning", "无点云输入");
		return;
	}
	else
	{
		float s_blocksize = blocksize.toFloat();
		float s_blockstride = blockstride.toFloat();
		for (size_t i = 0; i < cloud_vec.size(); i++)
		{
			if (cloud_index[i] == 1)
			{
				pointnet2_semseg* pointnet2 = new pointnet2_semseg(i, cloud_vec[i]->makeShared(), s_blocksize, s_blockstride, model_path);
				QThread* thread_m = new QThread;
				ui->progressBar->setRange(0, 100);  // 设置进度条范围
				ui->progressBar->setValue(0);  // 设置初始值
				ui->progressBar->setVisible(true);
				pointnet2->moveToThread(thread_m);
				QMessageBox::about(this, "about", "任务转移至子线程");
				connect(thread_m, &QThread::started, pointnet2, &pointnet2_semseg::semseg);
				connect(pointnet2, &pointnet2_semseg::finished, thread_m, &QThread::quit);
				connect(pointnet2, &pointnet2_semseg::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(pointnet2, &pointnet2_semseg::finished, this, &QTPCLVTK::onpointnet2semseg_Finished, Qt::QueuedConnection); // 连接带返回值的完成信号	
				connect(pointnet2, &pointnet2_semseg::finished, pointnet2, &QObject::deleteLater);
				connect(thread_m, &QThread::finished, thread_m, &QThread::deleteLater);
				thread_m->start();
			}
		}
	}
}
void QTPCLVTK::onpointnet2semseg_Finished(int mid, std::vector<int> label_indice)
{
	if (label_indice.empty())
	{
		QMessageBox::warning(this, "warning", "预测发生错误，返回原始点云");
	}
	else
	{
		//QMessageBox::warning(this, "warning", "点云输入");
		PointCloudT::Ptr cloud_filter(new PointCloudT);
		for (size_t i = 0; i < label_indice.size(); i++)
		{
			cloud_vec[mid]->points[i].label = label_indice[i];
		}
		view_updata(cloud_vec, cloud_index);
	}
	ui->progressBar->setVisible(false);
}


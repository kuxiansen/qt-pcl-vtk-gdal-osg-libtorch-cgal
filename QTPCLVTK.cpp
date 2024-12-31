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

	//���Ƴ�ʼ��
	m_Cloud.reset(new PointCloudT);
	//���ӻ������ʼ��
	viewer.reset(new pcl::visualization::PCLVisualizer("viewer", false));
	//����VTK���ӻ�����ָ��
	ui->qvtkWidget->SetRenderWindow(viewer->getRenderWindow());
	//���ô��ڽ��������ڿɽ��ܼ��̵��¼�
	viewer->setupInteractor(ui->qvtkWidget->GetInteractor(), ui->qvtkWidget->GetRenderWindow());
	
	//��������ϵ
	//viewer->addCoordinateSystem(1);
	//�������д������ļ�
	actptpick = new QAction(QIcon(":/img/image/open.png"), tr("��ѡ����"), this);
	actptpick->setCheckable(true);
	actptpick->setChecked(FALSE);
	ui->toolBar->addAction(actptpick);
	ui->toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);//��������ͼ����
	
	//***********************************connect����*******************************************//
	ui->behind_view->setIcon(QIcon(":/img/image/back.png"));
	ui->up_view->setIcon(QIcon(":/img/image/up.png"));
	ui->down_view->setIcon(QIcon(":/img/image/bottom.png"));
	ui->front_view->setIcon(QIcon(":/img/image/front.png"));
	ui->left_view->setIcon(QIcon(":/img/image/left.png"));
	ui->right_view->setIcon(QIcon(":/img/image/right.png"));
	ui->open_action->setIcon(QIcon(":/img/image/open.png"));
	stbar = this->statusBar();
	this->setStatusBar(stbar);
	//�������οؼ�
	model = new QStandardItemModel(this);
	//���ñ�ͷ����
	//ui->treeView->setHeaderHidden(true);
	model->setHorizontalHeaderLabels(QStringList() << "�����б�" );
	//����model 
	ui->treeView->setModel(model);
	//����չ��
	ui->treeView->expandAll();
	ui->progressBar->setVisible(false);
	//��ʼ��osgviewer
	//osgWidget = new osgQOpenGLWidget(ui->qosgWidget);
	//osgWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	connect(ui->opendem_action, SIGNAL(triggered()), this, SLOT(on_actionOpendem_triggered()));//��ȡdem�ļ�
	connect(ui->open_action, SIGNAL(triggered()), this, SLOT(on_actionOpen_triggered()));//��ȡ�ļ�
	connect(ui->quit_action, SIGNAL(triggered()), this, SLOT(on_actionClose_triggered()));//�ر��ļ�
	connect(ui->save_action, SIGNAL(triggered()), this, SLOT(on_actionSave_triggered()));//�����ļ�
	connect(model, SIGNAL(clicked(const QModelIndex & index)), this, SLOT(on_treeView_clicked(const QModelIndex & index)));
	//connect(model,SIGNAL(itemChanged(QStandardItem*)),this,SLOT(treeItemChanged(QStandardItem*)));
	connect(ui->qosgWidget, SIGNAL(initialized()), this, SLOT(initOSG()));//osg��ʼ��
	connect(actptpick, SIGNAL(triggered()), this, SLOT(on_point_select_clicked()));//����ѡȡ����
	connect(ui->up_view, SIGNAL(triggered()), this, SLOT(on_actionUp_triggered()));//����ͼ
	connect(ui->down_view, SIGNAL(triggered()), this, SLOT(on_actionBottom_triggered()));//����ͼ
	connect(ui->front_view, SIGNAL(triggered()), this, SLOT(on_actionFront_triggered()));//ǰ��ͼ
	connect(ui->behind_view, SIGNAL(triggered()), this, SLOT(on_actionBack_triggered()));//����ͼ
	connect(ui->left_view, SIGNAL(triggered()), this, SLOT(on_actionLeft_triggered()));//����ͼ
	connect(ui->right_view, SIGNAL(triggered()), this, SLOT(on_actionRight_triggered()));//����ͼ
	connect(ui->view_background, SIGNAL(triggered(bool)), this, SLOT(PressBtn_backgroundSelect()));//���ñ�����ɫ
	connect(ui->z_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_zcolor()));//�̸߳�ɫ
	connect(ui->class_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_classcolor()));//���ɫ
	connect(ui->pt_color, SIGNAL(triggered(bool)), this, SLOT(PressBtn_ptcolor()));//���Ƹ�ɫ
	connect(ui->pt_size, SIGNAL(triggered(bool)), this, SLOT(PressBtn_slider()));//���ƴ�С
	connect(ui->filtervoxel_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_voxelfilter()));//�����˲�
	connect(ui->filterzhitong_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_zhitongfilter()));//ֱͨ�˲�
	connect(ui->curdownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_cursamp()));//�����˲�����
	connect(ui->staticrandownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_staticrandownsamp()));//�̶������
	connect(ui->statisticalnoising_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_statisticalnoise()));//ͳ���˲�ȥ��
	connect(ui->radiusnoising_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_radiusdenoise()));//�뾶�˲�ȥ��
	connect(ui->gaussfilter_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_gausismooth()));//��˹�˲�ƽ��
	connect(ui->doubleline_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_doublelinesmooth()));//˫���˲�ƽ��
	connect(ui->mlsupsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_mlsupsamp()));//˫���˲�ƽ��
	connect(ui->morphologicalfilter_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_morphfilter()));	//������̬ѧ�˲�
	connect(ui->farthestdownsamp_action, SIGNAL(triggered(bool)), this, SLOT(PressBtn_farthestdownsamp()));	//��Զ�������
	connect(ui->s3dis_semsegaction, SIGNAL(triggered(bool)), this, SLOT(PressBtn_pointnet_sem_seg()));	//��Զ�������
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
	//ǿ��Ӣ�����뷨
	hCurKL = GetKeyboardLayout(0);
	LoadKeyboardLayout((LPCSTR)QString("0x0409").utf16(), KLF_ACTIVATE);
	QApplication::processEvents();
	
	//��������δѡȡʱ��ʵ�ʴ��������Ǹ�ѡ֮��ļ�true
	//if (actptpick->isChecked())
	//{
	//	//������꽻��
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
		//ģ��������롰x��
		INPUT input[2] = {};
		// ���ð��� 'X' ���������¼�
		input[0].type = INPUT_KEYBOARD;
		input[0].ki.wVk = 0;  // û�����������
		input[0].ki.wScan = 0x2D; // ɨ���룬'X' ����ɨ����
		input[0].ki.dwFlags = KEYEVENTF_SCANCODE; // ɨ����ģʽ
		input[0].ki.time = 0;
		input[0].ki.dwExtraInfo = 0;

		// �����ɿ� 'X' ���������¼�
		input[1].type = INPUT_KEYBOARD;
		input[1].ki.wVk = 0;  // û�����������
		input[1].ki.wScan = 0x2D; // ɨ���룬'X' ����ɨ����
		input[1].ki.dwFlags = KEYEVENTF_KEYUP | KEYEVENTF_SCANCODE; // �����ͷ�
		input[1].ki.time = 0;
		input[1].ki.dwExtraInfo = 0;
		SendInput(2, input, sizeof(INPUT));
	}
}

void QTPCLVTK::pp_callback_PointsSelect(const pcl::visualization::PointPickingEvent& event, void* args)
{
	QTPCLVTK* data = (QTPCLVTK*)args;
	// ��������Ч������
	if (event.getPointIndex() == -1)
		return;
	PointT current_point;
	event.getPoint(current_point.x, current_point.y, current_point.z);
	std::string info = "X: "+std::to_string(current_point.x) +", "+ "Y: " +std::to_string(current_point.y) +", " + "Z: " +std::to_string(current_point.z);
	stbar->clearMessage();
	stbar->showMessage(QString::fromStdString(info));
	//data->ui->qvtkWidget->update();
	//QMessageBox::information(data, "��ʾ", std::to_string(current_point.x).c_str());
	 // ���µ�����ѡ�еĵ�
	//data->clicked_points_3d->points.push_back(current_point);
	// ���µ�����ʾ������Ϊ��ɫ
	//pcl::visualization::PointCloudColorHandlerCustom<PointT> red(data->clicked_points_3d, 255, 0, 0);
	//data->viewerPtr->removePointCloud("clicked_points");
	//data->viewerPtr->addPointCloud(data->clicked_points_3d, red, "clicked_points");
	//data->viewerPtr->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "clicked_points");
}

void QTPCLVTK::initOSG()
{
	pViewer = ui->qosgWidget->getOsgViewer();
	pViewer->getCamera()->setClearColor(osg::Vec4(0.0f, 0.0f, 0.0f, 1.0f));  // ���������ɫΪ��ɫ
	osg::ref_ptr<osgGA::CameraManipulator>pManipulator = new osgGA::TrackballManipulator();
	pViewer->setCameraManipulator(pManipulator);
}

void QTPCLVTK::updataosgview()
{
	// �����µĳ������
	osg::BoundingSphere boundingSphere;
	osg::ComputeBoundsVisitor cbVisitor;
	m_pRootGroup->accept(cbVisitor);
	osg::BoundingBox& bb = cbVisitor.getBoundingBox();
	// �������λ��,���ӽǶȹ۲�ģ��
	double radius = osg::maximum(double(boundingSphere.radius()), 1e-6);
	pViewer->setCameraManipulator(NULL);
	//osgGA::CameraManipulator* pManipulator = new osgGA::TrackballManipulator();
	pViewer->setCameraManipulator(new osgGA::TrackballManipulator);
	//���ݷֱ���ȷ�����ʵ�ͶӰ����֤��ʾ��ͼ�β�����
	double fovy, aspectRatio, zNear, zFar;
	pViewer->getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, zFar);
	double newAspectRatio = double(ui->qosgWidget->width()) / double(ui->qosgWidget->height());
	double aspectRatioChange = newAspectRatio / aspectRatio;
	if (aspectRatioChange != 1.0)
	{
		//����ͶӰ����
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
//��dem
void QTPCLVTK::on_actionOpendem_triggered()
{
	// ��ȡdem·��
	QString path = QFileDialog::getOpenFileName(this, "ѡ��dem", ".//", "dem�ļ�(*.osg *.tif );;�����ļ�(*.*)");
	if (path.isEmpty())
		return;
	else
	{
		// ��ʼ�� GDAL
		GDALAllRegister();
		// �� DEM �ļ�
		GDALDataset* poDataset = (GDALDataset*)GDALOpen(path.toStdString().c_str(), GA_ReadOnly);
		if (!poDataset) 
		{
			qDebug() << "��ȡdem�ļ�ʧ��" << endl;
			return;
		}
		// ��ȡ����任��Ϣ
		double gdalGeoTransform[6];
		if (poDataset->GetGeoTransform(gdalGeoTransform) != CE_None) 
		{
			qDebug() << "��ȡdem�ļ�GeoTransform��Ϣʧ��" << endl;
			GDALClose(poDataset);
			return;
		}
		// ���� osg::HeightField ����
		osg::ref_ptr<osg::HeightField> hf = new osg::HeightField();
		hf->allocate(poDataset->GetRasterXSize(), poDataset->GetRasterYSize());
		hf->setOrigin(osg::Vec3(gdalGeoTransform[0], gdalGeoTransform[3], 0));
		hf->setXInterval(gdalGeoTransform[1]);
		hf->setYInterval(gdalGeoTransform[5]);
		// ��ȡ DEM �߳�����
		int rasterXSize = poDataset->GetRasterXSize();
		int rasterYSize = poDataset->GetRasterYSize();
		float* heightData = new float[rasterXSize * rasterYSize];
		// ��ȡդ������
		if (poDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, rasterXSize, rasterYSize, heightData, rasterXSize, rasterYSize, GDT_Float32, 0, 0) != CE_None) 
		{
			qDebug() << "��ȡdem�ļ�դ������ʧ��" << endl;
			delete[] heightData;
			GDALClose(poDataset);
			return ;
		}
		// ���� NoData ֵ
		float noDataValue = poDataset->GetRasterBand(1)->GetNoDataValue();
		float noDataValueFill = 0.0f;
		// ��������䵽 HeightField
		float* heightPtr = heightData;
		for (int r = rasterYSize - 1; r >= 0; --r) 
		{
			for (int c = 0; c < rasterXSize; ++c) 
			{
				float h = *heightPtr++;
				if (h != noDataValue)
					hf->setHeight(c, r, h);
				else
					hf->setHeight(c, r, noDataValueFill); // ��� NoData ����
			}
		}
		// ���� Locator �� TerrainTile
		osg::ref_ptr<osgTerrain::Locator> locator = new osgTerrain::Locator();
		double minX = std::min(gdalGeoTransform[0], gdalGeoTransform[0] + rasterXSize * gdalGeoTransform[1]);
		double minY = std::min(gdalGeoTransform[3], gdalGeoTransform[3] + rasterYSize * gdalGeoTransform[5]);
		double maxX = std::max(gdalGeoTransform[0], gdalGeoTransform[0] + rasterXSize * gdalGeoTransform[1]);
		double maxY = std::max(gdalGeoTransform[3], gdalGeoTransform[3] + rasterYSize * gdalGeoTransform[5]);

		locator->setTransformAsExtents(minX, minY, maxX, maxY);

		// ���� HeightFieldLayer ���� TerrainTile ��
		osg::ref_ptr<osgTerrain::HeightFieldLayer> hfl = new osgTerrain::HeightFieldLayer();
		hfl->setHeightField(hf);
		hfl->setLocator(locator.get());

		// ���� ImageLayer ��������ɫ
		osg::ref_ptr<osgTerrain::ImageLayer> colorLayer = new osgTerrain::ImageLayer();
		osg::ref_ptr<osg::Image> colorImage = new osg::Image();
		colorImage->allocateImage(rasterXSize, rasterYSize, 1, GL_RGBA, GL_UNSIGNED_BYTE);

		// ��ȡ��С�����߶�ֵ������ɫӳ��
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

		// �����ɫͼ��
		for (int r = 0; r < rasterYSize; ++r) {
			for (int c = 0; c < rasterXSize; ++c) {
				float h = hf->getHeight(c, r);
				if (h != noDataValueFill) {
					float normalizedHeight = (h - minHeight) / (maxHeight - minHeight);
					unsigned char rColor = static_cast<unsigned char>((1.0f - normalizedHeight) * 255);
					unsigned char gColor = static_cast<unsigned char>(normalizedHeight * 255);
					unsigned char bColor = 0;
					unsigned char aColor = 255; // alphaͨ��

					// ��ȡͼ�����������
					unsigned char* pixel = colorImage->data() + (r * rasterXSize + c) * 4;  // ÿ������ռ4���ֽڣ�RGBA��
					pixel[0] = rColor;  // ��ɫ
					pixel[1] = gColor;  // ��ɫ
					pixel[2] = bColor;  // ��ɫ
					pixel[3] = aColor;  // Alphaͨ��
				}
				else {
					// ����NoData����Ϊ��ɫ
					unsigned char* pixel = colorImage->data() + (r * rasterXSize + c) * 4;
					pixel[0] = 0;  // ��ɫ
					pixel[1] = 0;  // ��ɫ
					pixel[2] = 0;  // ��ɫ
					pixel[3] = 0;  // Alphaͨ��
				}
			}
		}


		colorLayer->setImage(colorImage.get());

		osg::ref_ptr<osgTerrain::TerrainTile> terrainTile = new osgTerrain::TerrainTile();
		terrainTile->setElevationLayer(hfl);
		terrainTile->setColorLayer(1, colorLayer.get());
		// ���� StateSet������������Ⱦ״̬�����ʣ�
		osg::ref_ptr<osg::StateSet> stateSet = terrainTile->getOrCreateStateSet();
		// �������ʶ��󲢽���Ӧ�õ� StateSet
		osg::ref_ptr<osg::Material> material = new osg::Material();
		material->setDiffuse(osg::Material::FRONT_AND_BACK, osg::Vec4(0.5f, 0.5f, 0.5f, 1.0f));  // ������������ɫ
		material->setSpecular(osg::Material::FRONT_AND_BACK, osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));  // ���þ��淴����ɫ
		material->setShininess(osg::Material::FRONT_AND_BACK, 50.0f);  // ���ù����
		stateSet->setAttributeAndModes(material.get(), osg::StateAttribute::ON);  // ������Ӧ�õ� StateSet
		m_pRootGroup->addChild(terrainTile.get());
		pViewer->setSceneData(m_pRootGroup);
		pViewer->addEventHandler(new osgViewer::WindowSizeHandler());
		// �ͷ��ڴ�
		delete[] heightData;
		GDALClose(poDataset);
		updataosgview();
	}
}
//�򿪵���
void QTPCLVTK::on_actionOpen_triggered()
{

	//��ȡ����·��
	QString path = QFileDialog::getOpenFileName(this, "ѡ�����", ".//", "�����ļ�(*.ply *.pcd *.txt *.las);;�����ļ�(*.*)");
	
	pcl::PointCloud<pcl::PointXYZL>::Ptr cloud_tmp(new pcl::PointCloud<pcl::PointXYZL>);

	if (path.isEmpty())
		return;

	if (path.endsWith(".pcd", Qt::CaseInsensitive))
	{
		qDebug() << path;
		if (pcl::io::loadPCDFile(path.toStdString(), *cloud_tmp) == -1)
		{
			qDebug()<< "��ȡpcd����ʧ��" << endl;
			return;
		}
		else
		{
			QMessageBox::information(this, "��ʾ", "pcd���ƶ�ȡ���");
		}
	}
	else if (path.endsWith(".ply", Qt::CaseInsensitive))
	{
		qDebug() << path;
		if (pcl::io::loadPLYFile(path.toStdString(), *cloud_tmp) == -1)
		{
			qDebug() << "��ȡply����ʧ��" << endl;
			return;
		}
		else
		{
			QMessageBox::information(this, "��ʾ", "ply���ƶ�ȡ���");
		}
	}
	else if (path.endsWith(".txt", Qt::CaseInsensitive))
	{
		qDebug() << path;
		CreateCloudFromTxt(path.toStdString(), cloud_tmp);
		QMessageBox::information(this, "��ʾ", "txt���ƶ�ȡ���");
	}
	else if (path.endsWith(".las", Qt::CaseInsensitive))
	{
		qDebug() << path;
		CreateCloudFromlas(path.toStdString(), cloud_tmp);
		QMessageBox::information(this, "��ʾ", "las���ƶ�ȡ���");
	}
	else
	{
		QMessageBox::warning(this, "Warning", "���ƶ�ȡ��ʽ����");
	}


	//��յ���
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
		itemFolder->setCheckState(Qt::Checked);//��ȡѡ��״̬
		model->appendRow(itemFolder);
	}
	view_updata(cloud_vec, cloud_index);
	//PointT max_p;
	//PointT min_p;
	//int nums = 0;
	//for (size_t i = 0; i < cloud_vec.size(); i++)
	//{
	//	//��ӵ�����
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
	////�����ӽ�
	//viewer->resetCamera();
	//std::string num = std::to_string(nums);
	//ui->canshu->setText("��������Ϊ:"+ QString::fromStdString(num));
	//ui->canshu->append("����x�᷶Χ:" + QString::fromStdString(std::to_string(p_min.x))+" -" + QString::fromStdString(std::to_string(p_max.x)));
	//ui->canshu->append("����y�᷶Χ:" + QString::fromStdString(std::to_string(p_min.y)) + " -" + QString::fromStdString(std::to_string(p_max.y)));
	//ui->canshu->append("����z�᷶Χ:" + QString::fromStdString(std::to_string(p_min.z)) + " -" + QString::fromStdString(std::to_string(p_max.z)));
	////ˢ�´���
	//ui->qvtkWidget->update();
}
//�رյ���
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
		model->setHorizontalHeaderLabels(QStringList() << "�����б�");
		if (cloud_vec.size() >= 1)
		{
			for (size_t i = 0; i < cloud_vec.size(); i++)
			{
				itemFolder = new QStandardItem(cloud_path[i]);
				itemFolder->setCheckable(true);
				itemFolder->setCheckState(Qt::Checked);//��ȡѡ��״̬
				model->appendRow(itemFolder);
			}
			//QStandardItem* itemFolder = new QStandardItem(QStringLiteral("cloud%1").arg(cloud_vec.size() - 1)
		}
		view_updata(cloud_vec, cloud_index);
	}
	else
	{
		QMessageBox::warning(this, "Warning", "��ͼ��û�е��ƣ�");
		return;
	}
}
//�������
void QTPCLVTK::on_actionSave_triggered()
{
	if (cloud_vec.size() != 0)
	{
		QString path = QFileDialog::getSaveFileName(this, "ѡ�����", ".//", "�����ļ�(*.ply *.pcd);;�����ļ�(*.*)");
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
			QMessageBox::information(this, "��ʾ", "pcd���Ʊ���ɹ�");
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
			QMessageBox::information(this, "��ʾ", "ply���Ʊ���ɹ�");
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
			QMessageBox::information(this, "��ʾ", "ply���Ʊ���ɹ�");
		}
	}
	else
	{
		QMessageBox::warning(this, "Warning", "��ͼ��û�е��ƣ�");
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
		Qt::CheckState state = item->checkState();//��ȡ��ǰѡ��״̬
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
		//��ӵ�����
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
	ui->canshu->setText("��������Ϊ:" + QString::fromStdString(num));
	ui->canshu->append("����x�᷶Χ:" + QString::fromStdString(std::to_string(p_min.x)) + " -" + QString::fromStdString(std::to_string(p_max.x)));
	ui->canshu->append("����y�᷶Χ:" + QString::fromStdString(std::to_string(p_min.y)) + " -" + QString::fromStdString(std::to_string(p_max.y)));
	ui->canshu->append("����z�᷶Χ:" + QString::fromStdString(std::to_string(p_min.z)) + " -" + QString::fromStdString(std::to_string(p_max.z)));
	//ui->canshu->append("�������Χ:" + QString::fromStdString(std::to_string(p_min.label)) + " -" + QString::fromStdString(std::to_string(p_max.label)));
	ui->qvtkWidget->update();
}




//����ͼ
void QTPCLVTK::on_actionUp_triggered()
{
	if (cloud_vec.size()!=0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_max.z + 2 * maxLen, 0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_max.z, 0, 1, 0);
		ui->qvtkWidget->update();
	}
}
//����ͼ
void QTPCLVTK::on_actionBottom_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_min.z - 2 * maxLen, 0.5 * (p_min.x + p_max.x), 0.5 * (p_min.y + p_max.y), p_min.z, 0, 1, 0);
		ui->qvtkWidget->update();
	}
}
//ǰ��ͼ
void QTPCLVTK::on_actionFront_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), p_min.y - 2 * maxLen, 0.5 * (p_min.z + p_max.z), 0.5 * (p_min.x + p_max.x), p_min.y, 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//����ͼ
void QTPCLVTK::on_actionBack_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(0.5 * (p_min.x + p_max.x), p_max.y + 2 * maxLen, 0.5 * (p_min.z + p_max.z), 0.5 * (p_min.x + p_max.x), p_min.y, 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//����ͼ
void QTPCLVTK::on_actionLeft_triggered()
{
	if (cloud_vec.size() != 0)
	{
		viewer->setCameraPosition(p_min.x - 2 * maxLen, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), p_max.x, 0.5 * (p_min.y + p_max.y), 0.5 * (p_min.z + p_max.z), 0, 0, 1);
		ui->qvtkWidget->update();
	}
}
//����ͼ
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

//���ñ�����ɫ
void QTPCLVTK::PressBtn_backgroundSelect()
{
	dialog_colorselect = new pcl_view_select_color();
	QColor color = dialog_colorselect->getColor();
	viewer->setBackgroundColor(color.redF(), color.greenF(), color.blueF());
	return;
}
//�̸߳�ɫ
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
		QMessageBox::warning(this, "Warning", "��ͼ��û�е��ƣ�");
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
		QMessageBox::warning(this, "Warning", "��ͼ��û�е��ƣ�");
	}
}
//���õ�����ɫ
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
		QMessageBox::warning(this, "Warning", "��ͼ��û�е��ƣ�");
	}
}

//��Ĵ�С����
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
//�����˲�����
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//ֱͨ�˲�����
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//�̶������������
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//���ʲ�������
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
				
		}
		view_updata(cloud_vec, cloud_index);
	}
}


//ֱͨ�˲�ȥ�봰��
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
						QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
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
						QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
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

//�뾶�˲�����
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//��˹�˲�����
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}

//˫���˲�����
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}
//mls�ϲ�������
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
					QMessageBox::warning(this, "warning", "�������Ϊ�ջ�����,�������������");
				}
				else {
					cloud_vec[i] = std::move(cloudout);
				}
			}
		}
		view_updata(cloud_vec, cloud_index);
	}
}
//��Զ�������
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
				ui->progressBar->setRange(0, 100);  // ���ý�������Χ
				ui->progressBar->setValue(0);  // ���ó�ʼֵ
				ui->progressBar->setVisible(true);
				funcer->moveToThread(thread);
				QMessageBox::about(this, "about", "����ת�������߳�");
				connect(thread, &QThread::started, funcer, &farest_sampling::sampling);
				connect(funcer, &farest_sampling::finished, thread, &QThread::quit);
				connect(funcer, &farest_sampling::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(funcer, &farest_sampling::finished, this, &QTPCLVTK::onfarthestdownsampFinished,Qt::QueuedConnection); // ���Ӵ�����ֵ������ź�	
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
		QMessageBox::warning(this, "warning", "��ֵ���󣬷���ԭʼ����");
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
//������̬ѧ�˲�
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
				ui->progressBar->setRange(0, 100);  // ���ý�������Χ
				ui->progressBar->setValue(0);  // ���ó�ʼֵ
				ui->progressBar->setVisible(true);
				morph->moveToThread(thread_m);
				QMessageBox::about(this, "about", "����ת�������߳�");
				connect(thread_m, &QThread::started, morph, &morph_filter::my_filtering);
				connect(morph, &morph_filter::output, thread_m, &QThread::quit);
				connect(morph, &morph_filter::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(morph, &morph_filter::output, this, &QTPCLVTK::onmorphfilterFinished, Qt::QueuedConnection); // ���Ӵ�����ֵ������ź�	
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
		QMessageBox::warning(this, "warning", "��ֵ���󣬷���ԭʼ����");
	}
	else
	{
		//QMessageBox::warning(this, "warning", "��������");
		PointCloudT::Ptr cloud_filter(new PointCloudT);
		pcl::copyPointCloud(*cloud_vec[mid], sam_indice, *cloud_filter);
		cloud_vec[mid] = cloud_filter;
		view_updata(cloud_vec, cloud_index);
	}
}

////pointnet++��������ָ�
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
		QMessageBox::warning(this, "warning", "�޵�������");
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
				ui->progressBar->setRange(0, 100);  // ���ý�������Χ
				ui->progressBar->setValue(0);  // ���ó�ʼֵ
				ui->progressBar->setVisible(true);
				pointnet2->moveToThread(thread_m);
				QMessageBox::about(this, "about", "����ת�������߳�");
				connect(thread_m, &QThread::started, pointnet2, &pointnet2_semseg::semseg);
				connect(pointnet2, &pointnet2_semseg::finished, thread_m, &QThread::quit);
				connect(pointnet2, &pointnet2_semseg::progress, this, [=](int progressValue) {ui->progressBar->setValue(progressValue); });
				connect(pointnet2, &pointnet2_semseg::finished, this, &QTPCLVTK::onpointnet2semseg_Finished, Qt::QueuedConnection); // ���Ӵ�����ֵ������ź�	
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
		QMessageBox::warning(this, "warning", "Ԥ�ⷢ�����󣬷���ԭʼ����");
	}
	else
	{
		//QMessageBox::warning(this, "warning", "��������");
		PointCloudT::Ptr cloud_filter(new PointCloudT);
		for (size_t i = 0; i < label_indice.size(); i++)
		{
			cloud_vec[mid]->points[i].label = label_indice[i];
		}
		view_updata(cloud_vec, cloud_index);
	}
	ui->progressBar->setVisible(false);
}


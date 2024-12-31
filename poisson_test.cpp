//#include <iostream>
//#include <string>
//#include <pcl/memory.h>
//#include <pcl/pcl_macros.h>
//#include <pcl/surface/reconstruction.h>
//#include <pcl\io\pcd_io.h>
//#include <pcl\io\ply_io.h>
//#include <pcl\point_types.h>
//#include <pcl\kdtree\kdtree_flann.h>
//#include <pcl\features\normal_3d.h>
//#include <pcl\features\normal_3d_omp.h>
//#include <pcl/common/common.h>
//#include <pcl/common/vector_average.h>
//#include <pcl/Vertices.h>
//#include <pcl\surface\poisson.h>
//#include <pcl\surface\gp3.h>
//#include <pcl\visualization\cloud_viewer.h>
//#include <pcl\visualization\pcl_visualizer.h>
//#include <pcl/surface/3rdparty/poisson4/octree_poisson.h>
//#include <pcl/surface/3rdparty/poisson4/sparse_matrix.h>
//#include <pcl/surface/3rdparty/poisson4/function_data.h>
//#include <pcl/surface/3rdparty/poisson4/ppolynomial.h>
//#include <pcl/surface/3rdparty/poisson4/multi_grid_octree_data.h>
//#include <pcl/surface/3rdparty/poisson4/geometry.h>
//
//#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
//#include <cstdarg>
////���߳�
//#include <boost\thread\thread.hpp>
//
//#include <vector>
//using namespace std;
//
//
//template <int Degree> void ipsr(const pcl::PCLPointCloud2::ConstPtr& input, pcl::PolygonMesh& output_mesh,
//	float samples_per_node, int min_depth, float scale, int depth, int solver_divide, int iso_divide, int num_iteration,
//	int k_neighbors, float point_weight, double resolution)
//{
//
//}
//
//int main()
//{
//	pcl::PointCloud<pcl::PointXYZ> ::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::io::loadPLYFile("K:\\scans\\train\\scene0706_00_vh_clean_2.labels.ply",*cloud);
//	cout << cloud->points.size() << endl;
//
//	// ���㷨����
//	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>); //���������ƶ���ָ��
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;//���߹��ƶ���
//	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);//�洢���Ƶķ��ߵ�ָ��
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	tree->setInputCloud(cloud);
//	n.setInputCloud(cloud);
//	n.setSearchMethod(tree);
//	n.setKSearch(20);
//	n.compute(*normals); //���㷨�ߣ�����洢��normals��
//	//�����ƺͷ��߷ŵ�һ��
//	pcl::concatenateFields(*cloud, *normals, *cloud_with_normals);
//	//����Poisson���󣬲����ò���
//	pcl::Poisson<pcl::PointNormal> pn;
//	pcl::poisson::Point3D<float> center;
//	const int degree = 2;//���ò���degree[1,5],ֵԽ��Խ��ϸ����ʱԽ�á�
//	int depth = 8;//���������ȣ����2^d x 2^d x 2^d������Ԫ�����ڰ˲�������Ӧ�����ܶȣ�ָ��ֵ��Ϊ�����ȡ�
//	int min_depth = 5;//������С���
//	int pnthreads = 1;//�߳���
//	bool confidense = false;//�Ƿ�ʹ�÷������Ĵ�С��Ϊ������Ϣ�����false�����з���������һ����
//	int SolverDivide = 8;//����������Է������Gauss-Seidel�������������
//	bool non_adaptive_weights = false;
//	int iso_divide = 8;//������ȡISO��ֵ����㷨�����
//	bool show_residual = false;
//	int min_iterations = 8;
//	float solver_accuracy = 1e-3f;
//	bool manifold = true;//�Ƿ���Ӷ���ε����ģ�����������ǻ�ʱ�� �������б�־���������Ϊtrue����Զ���ν���ϸ�����ǻ�ʱ������ģ�����false�����
//	bool output_polygons = false;
//	pcl::poisson::Real iso_value = 0;
//	pcl::poisson::Real scale = 1.1; //���������ع���������ֱ���������߽�������ֱ���ı��ʡ�
//	pcl::poisson::Real scalenum = 2;//���������ع���������ֱ���������߽�������ֱ���ı�ֵ��
//	pcl::poisson::Real samplenode = 1;//��������һ���˲�������е����������С������
//	pcl::poisson::Real pointweight = 4;//���õ��Ȩ��
//	pcl::poisson::TreeNodeData::UseIndex = 1;
//	pn.setConfidence(confidense); 
//	pn.setDegree(degree); 
//	pn.setDepth(depth); 
//	pn.setIsoDivide(iso_divide); 
//	pn.setManifold(manifold); 
//	pn.setOutputPolygons(false); //�Ƿ������������񣨶��������ǻ��ƶ�������Ľ����
//	pn.setSamplesPerNode(samplenode); 
//	pn.setScale(scale); 
//	pn.setSolverDivide(SolverDivide); 
//	pn.setInputCloud(cloud_with_normals);
//	//����������������ڴ洢���
//	pcl::PolygonMesh mesh;
//
//
//
//	// ����һ��Octreeʵ��
//	pcl::poisson::Octree<Degree> octree;
//	octree.threads = pnthreads;
//	if (SolverDivide < min_depth)
//	{
//		SolverDivide = min_depth;
//	}
//	if (iso_divide < min_depth)
//	{
//		iso_divide = min_depth;
//	}
//	pcl::poisson::TreeOctNode::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);// ���ð˲����ڵ������
//	int kernel_depth_ = depth - 2;
//	octree.setBSplineData(depth, pcl::poisson::Real(1.0 / (1 << depth)), true);// ����BSpline����
//	octree.maxMemoryUsage = 0;
//	int point_count = octree.setTree<pcl::PointNormal>
//		(cloud_with_normals,
//			depth, min_depth, kernel_depth_,
//			samplenode, scale, center, scalenum,
//			confidense, pointweight, !non_adaptive_weights);//���ð˲���
//	octree.ClipTree();// �ü��˲���
//	octree.finalize();// ��ɰ˲�������
//	octree.RefineBoundary(iso_divide);// ��ϸ���߽磬octree�������
//	PCL_DEBUG("Input Points: %d\n", point_count);
//	PCL_DEBUG("Leaves/Nodes: %d/%d\n", octree.tree.leaves(), octree.tree.nodes());
//	//��ȡoctree���
//	int octree_depth = octree.tree.maxDepth();
//
//
//
//	//ִ���ع�
//	pn.performReconstruction(mesh);
//	//��������ͼ
//	pcl::io::savePLYFile("FanBlade.ply", mesh);
//	// ��ʾ���ͼ
//	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D viewer"));
//	viewer->setBackgroundColor(0.5, 0.5, 1);
//	viewer->addPolygonMesh(mesh, "Blade");
//	viewer->addCoordinateSystem(50.0);
//	viewer->initCameraParameters();
//	while (!viewer->wasStopped()) {
//		viewer->spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//	system("pause");
//	return 0;
//}
//
//

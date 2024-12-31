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
////多线程
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
//	// 计算法向量
//	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>); //法向量点云对象指针
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;//法线估计对象
//	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);//存储估计的法线的指针
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	tree->setInputCloud(cloud);
//	n.setInputCloud(cloud);
//	n.setSearchMethod(tree);
//	n.setKSearch(20);
//	n.compute(*normals); //计算法线，结果存储在normals中
//	//将点云和法线放到一起
//	pcl::concatenateFields(*cloud, *normals, *cloud_with_normals);
//	//创建Poisson对象，并设置参数
//	pcl::Poisson<pcl::PointNormal> pn;
//	pcl::poisson::Point3D<float> center;
//	const int degree = 2;//设置参数degree[1,5],值越大越精细，耗时越久。
//	int depth = 8;//树的最大深度，求解2^d x 2^d x 2^d立方体元。由于八叉树自适应采样密度，指定值仅为最大深度。
//	int min_depth = 5;//树的最小深度
//	int pnthreads = 1;//线程数
//	bool confidense = false;//是否使用法向量的大小作为置信信息。如果false，所有法向量均归一化。
//	int SolverDivide = 8;//设置求解线性方程组的Gauss-Seidel迭代方法的深度
//	bool non_adaptive_weights = false;
//	int iso_divide = 8;//用于提取ISO等值面的算法的深度
//	bool show_residual = false;
//	int min_iterations = 8;
//	float solver_accuracy = 1e-3f;
//	bool manifold = true;//是否添加多边形的重心，当多边形三角化时。 设置流行标志，如果设置为true，则对多边形进行细分三角话时添加重心，设置false则不添加
//	bool output_polygons = false;
//	pcl::poisson::Real iso_value = 0;
//	pcl::poisson::Real scale = 1.1; //设置用于重构的立方体直径和样本边界立方体直径的比率。
//	pcl::poisson::Real scalenum = 2;//设置用于重构的立方体直径和样本边界立方体直径的比值。
//	pcl::poisson::Real samplenode = 1;//设置落入一个八叉树结点中的样本点的最小数量。
//	pcl::poisson::Real pointweight = 4;//设置点的权重
//	pcl::poisson::TreeNodeData::UseIndex = 1;
//	pn.setConfidence(confidense); 
//	pn.setDegree(degree); 
//	pn.setDepth(depth); 
//	pn.setIsoDivide(iso_divide); 
//	pn.setManifold(manifold); 
//	pn.setOutputPolygons(false); //是否输出多边形网格（而不是三角化移动立方体的结果）
//	pn.setSamplesPerNode(samplenode); 
//	pn.setScale(scale); 
//	pn.setSolverDivide(SolverDivide); 
//	pn.setInputCloud(cloud_with_normals);
//	//创建多变形网格，用于存储结果
//	pcl::PolygonMesh mesh;
//
//
//
//	// 创建一个Octree实例
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
//	pcl::poisson::TreeOctNode::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);// 设置八叉树节点分配器
//	int kernel_depth_ = depth - 2;
//	octree.setBSplineData(depth, pcl::poisson::Real(1.0 / (1 << depth)), true);// 设置BSpline数据
//	octree.maxMemoryUsage = 0;
//	int point_count = octree.setTree<pcl::PointNormal>
//		(cloud_with_normals,
//			depth, min_depth, kernel_depth_,
//			samplenode, scale, center, scalenum,
//			confidense, pointweight, !non_adaptive_weights);//设置八叉树
//	octree.ClipTree();// 裁剪八叉树
//	octree.finalize();// 完成八叉树构建
//	octree.RefineBoundary(iso_divide);// 精细化边界，octree构建完成
//	PCL_DEBUG("Input Points: %d\n", point_count);
//	PCL_DEBUG("Leaves/Nodes: %d/%d\n", octree.tree.leaves(), octree.tree.nodes());
//	//读取octree深度
//	int octree_depth = octree.tree.maxDepth();
//
//
//
//	//执行重构
//	pn.performReconstruction(mesh);
//	//保存网格图
//	pcl::io::savePLYFile("FanBlade.ply", mesh);
//	// 显示结果图
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

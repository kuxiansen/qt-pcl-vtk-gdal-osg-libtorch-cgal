#pragma once
#include <torch/script.h>
#include <pcl/common/geometry.h>
#include <pcl/segmentation/progressive_morphological_filter.h>
#include<pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include<pcl/kdtree/kdtree_flann.h>
#include <pcl/octree/octree_search.h>
#include<pcl/filters/random_sample.h>
#include<pcl/filters/statistical_outlier_removal.h>
#include<pcl/filters/radius_outlier_removal.h>
#include<pcl/features/normal_3d_omp.h>
#include<pcl/filters/conditional_removal.h>
#include<pcl/filters/convolution_3d.h>
#include<pcl/surface/mls.h>
#include<vector>
#include<thread>
#include<random>
#include <qmessagebox.h>
#include<time.h>
#include <limits>
#include<omp.h>
#include <algorithm>
#include <QObject>
#include <QThread>
#include <QString>
typedef pcl::PointXYZL PointT;
typedef pcl::PointCloud<PointT> PointCloudT;


//最远点采样----------------------------------------------------------------
class farest_sampling : public QObject {
    Q_OBJECT
public:
    farest_sampling(int id,PointCloudT::Ptr cloud, int samples):myid(id),cloud_in(cloud),sample_nums(samples){}
public slots:
    //最远点采样----------------------------------------------------------------
    void sampling();
signals:
    void progress(int value); // 进度更新信号
    void finished(int mid, std::vector<int>sam_indice);
private:
    PointCloudT::Ptr cloud_in;
    int sample_nums;
    int myid;
};

//渐进形态学滤波----------------------------------------------------------------
class morph_filter : public QObject {
    Q_OBJECT
public:
    morph_filter(int id, PointCloudT::Ptr cloud, int max_grid, float m_slope, float m_heightin, float m_heightmax) 
        :myid(id), cloud_in(cloud), grid(max_grid), slope(m_slope),height_in(m_heightin), height_max(m_heightmax){}
public slots:
    void my_filtering();
    void filtering();
signals:
    void progress(int value); // 进度更新信号
    void output(int mid, std::vector<int>sam_indice);
private:
    PointCloudT::Ptr cloud_in;
    int myid;
    int grid;
    float slope;
    float height_in;
    float height_max;
};

//pointnet++s3dis_semseg----------------------------------------------------------------
class pointnet2_semseg :public QObject
{
    Q_OBJECT
public:
    pointnet2_semseg(int id,PointCloudT::Ptr cloud, float block_size, float block_stride, QString model_path,int point_num = 4096, int class_num = 13)
        :myid(id),my_blocksize(block_size), cloud_in(cloud), my_blockstride(block_stride), mymodel_path(model_path), mypointnum(point_num), myclassnum(class_num) {}
public slots:
    void semseg();
signals:
    void progress(int value); // 进度更新信号
    void finished(int mid, std::vector<int>sam_indice);
private:
    PointCloudT::Ptr cloud_in;
    QString mymodel_path;
    int myid;
    float my_blocksize;
    float my_blockstride;
    int mypointnum;
    int myclassnum;
};

//法向量计算----------------------------------------------------------------
pcl::PointCloud<pcl::Normal>::Ptr computeNormal(pcl::PointCloud<pcl::PointXYZ>::ConstPtr target_cloud,const int& k_search);
//voxel_filter----------------------------------------------------------------
PointCloudT::Ptr pcl_filter_voxel(PointCloudT::Ptr cloud_in, const float& leaf_size, const unsigned int& minnum);
PointCloudT::Ptr pcl_kdtree(PointCloudT::Ptr cloud_in, PointCloudT::Ptr cloud_filter,
    const int& k, std::vector<int>& Idxsearch, std::vector<float>& IdxDistance);
//zhitong_filter----------------------------------------------------------------
PointCloudT::Ptr pcl_filter_zhitong(PointCloudT::Ptr cloud_in, const int& canshu, const float& min_size, const float& max_size);
//zhitong_filter_label----------------------------------------------------------------
PointCloudT::Ptr pcl_filter_zhitong_label(PointCloudT::Ptr cloud_in, const float& min_size, const float& max_size);
//固定点随机采样----------------------------------------------------------------
PointCloudT::Ptr pcl_static_randomsamp(PointCloudT::Ptr cloud_in, const unsigned int& num);
//条件曲率滤波采样----------------------------------------------------------------
PointCloudT::Ptr pcl_cursamp(PointCloudT::Ptr cloud_in, const int& K_Search, const float& num);
//移动最小二乘法上采样----------------------------------------------------------------
PointCloudT::Ptr pcl_mlsupsamp(PointCloudT::Ptr cloud_in, const float& radius, const float& step);
//直通滤波去噪----------------------------------------------------------------pcl版
PointCloudT::Ptr pcl_statiscal_denoise(PointCloudT::Ptr cloud_in, const double& stddev, const int& num);
//直通滤波去噪----------------------------------------------------------------openmp+手写版
PointCloudT::Ptr my_statiscal_denoise(PointCloudT::Ptr cloud_in, const double& stddev, const int& num);
//半径滤波去噪----------------------------------------------------------------
PointCloudT::Ptr pcl_radiussamp(PointCloudT::Ptr cloud_in, const float& radius, const int& num);
//高斯滤波平滑----------------------------------------------------------------
PointCloudT::Ptr pcl_gausifilter(PointCloudT::Ptr cloud_in, const float& sigma, const float& sigmarela, const float& threshold, const float& radius);
//双边滤波平滑----------------------------------------------------------------
PointCloudT::Ptr my_doublelinefilter(PointCloudT::Ptr cloud_in, const int& nums, const float& sigmadis, const float& sigmanormal);

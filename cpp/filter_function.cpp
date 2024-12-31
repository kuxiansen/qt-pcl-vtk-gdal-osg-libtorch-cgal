#include "filter_function.h"

//法向量计算----------------------------------------------------------------
pcl::PointCloud<pcl::Normal>::Ptr computeNormal(pcl::PointCloud<pcl::PointXYZ>::ConstPtr target_cloud, const int& k_search)
{
    pcl::search::KdTree <pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree <pcl::PointXYZ>);
    pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal>n;
    kdtree->setInputCloud(target_cloud);
    n.setInputCloud(target_cloud);
    n.setSearchMethod(kdtree);
    n.setKSearch(k_search);
    n.setNumberOfThreads(std::thread::hardware_concurrency());//设置omp线程数
    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    n.compute(*normals);
    return normals;
}

//voxel_filter----------------------------------------------------------------
PointCloudT::Ptr pcl_filter_voxel(PointCloudT::Ptr cloud_in, const float& leaf_size, const unsigned int& minnum)
{
    pcl::VoxelGrid<pcl::PointXYZL> voxel_grid;
    voxel_grid.setInputCloud(cloud_in);
    voxel_grid.setLeafSize(leaf_size, leaf_size, leaf_size);
    voxel_grid.setMinimumPointsNumberPerVoxel(minnum);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    voxel_grid.filter(*cloud_filter);
    int K = 1;
    std::vector<int>pointIdxKNNsearch(K);
    std::vector<float>pointIdxSquaredDistance(K);
    return pcl_kdtree(cloud_in, cloud_filter, K, pointIdxKNNsearch, pointIdxSquaredDistance);
}
PointCloudT::Ptr pcl_kdtree(PointCloudT::Ptr cloud_in, PointCloudT::Ptr cloud_filter,
    const int&k, std::vector<int>&Idxsearch, std::vector<float>&IdxDistance)
{
    pcl::KdTreeFLANN<pcl::PointXYZL>kdtree;
    kdtree.setInputCloud(cloud_in);
    std::vector<int>dices;//采样后根据最邻近点提取的样本点下标索引
    for (size_t i = 0; i < cloud_filter->points.size(); i++)
    {
        kdtree.nearestKSearch(cloud_filter->points[i], k, Idxsearch, IdxDistance);
        dices.push_back(Idxsearch[0]);
    }
    PointCloudT::Ptr cloud_final(new PointCloudT);
    pcl::copyPointCloud(*cloud_in, dices, *cloud_final);
    return cloud_final;
}
//zhitong_filter----------------------------------------------------------------
PointCloudT::Ptr pcl_filter_zhitong(PointCloudT::Ptr cloud_in, const int& canshu, const float& min_size, const float& max_size)
{
    pcl::PassThrough<pcl::PointXYZL>pass(true);
    pass.setInputCloud(cloud_in);
    switch (canshu)
    {
    case 0:
        pass.setFilterFieldName("x");
        break;
    case 1:
        pass.setFilterFieldName("y");
        break;
    case 2:
        pass.setFilterFieldName("z");
        break;
    case 3:
        return pcl_filter_zhitong_label(cloud_in, min_size, max_size);
        //pass.setFilterFieldName("label");
        //break;
    }
    pass.setFilterLimits(min_size, max_size);
    pass.setNegative(false);//获取给定范围点
    std::vector<int>indices;
    pass.filter(indices);
    PointCloudT::Ptr cloud_out(new PointCloudT);
    pcl::copyPointCloud(*cloud_in, indices, *cloud_out);
    return cloud_out;
}
PointCloudT::Ptr pcl_filter_zhitong_label(PointCloudT::Ptr cloud_in, const float& min_size, const float& max_size)
{
    PointCloudT::Ptr cloud_out(new PointCloudT);

    for (const auto& point : cloud_in->points)
    {
        if (point.label<=max_size&& point.label>=min_size)
        {
            cloud_out->points.push_back(point);
        }
    }
    return cloud_out;
}
//固定点随机采样----------------------------------------------------------------
PointCloudT::Ptr pcl_static_randomsamp(PointCloudT::Ptr cloud_in, const unsigned int& num)
{
    pcl::RandomSample<pcl::PointXYZL>rs;
    rs.setInputCloud(cloud_in);
    rs.setSample(num);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    rs.filter(*cloud_filter);
    return cloud_filter;
}
//条件曲率滤波采样----------------------------------------------------------------
PointCloudT::Ptr pcl_cursamp(PointCloudT::Ptr cloud_in, const int& K_Search, const float& num)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr XYZ(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::copyPointCloud(*cloud_in, *XYZ);
    pcl::PointCloud<pcl::Normal>::Ptr normals(computeNormal(XYZ,K_Search));
    pcl::PointCloud<pcl::PointXYZLNormal>::Ptr cnormals(new pcl::PointCloud<pcl::PointXYZLNormal>);
    float threshold = num;
    pcl::concatenateFields(*cloud_in, *normals, *cnormals);
    pcl::ConditionOr<pcl::PointXYZLNormal>::Ptr range_cond(new pcl::ConditionOr<pcl::PointXYZLNormal>);
    range_cond->addComparison(pcl::FieldComparison<pcl::PointXYZLNormal>::ConstPtr
    (new pcl::FieldComparison<pcl::PointXYZLNormal>("curvature", pcl::ComparisonOps::GT, threshold)));
    pcl::ConditionalRemoval<pcl::PointXYZLNormal>condrem;
    condrem.setCondition(range_cond);
    condrem.setInputCloud(cnormals);
    pcl::PointCloud<pcl::PointXYZLNormal>::Ptr filter_normals(new pcl::PointCloud<pcl::PointXYZLNormal>);
    condrem.filter(*filter_normals);
    PointCloudT::Ptr cloud_final(new PointCloudT);
    //pcl::copyPointCloud(*cnormals, *cloud_final);
    pcl::copyPointCloud(*filter_normals,  *cloud_final);
    return cloud_final;
    
}
//移动最小二乘法上采样----------------------------------------------------------------
PointCloudT::Ptr pcl_mlsupsamp(PointCloudT::Ptr cloud_in, const float& radius, const float& step)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr XYZ(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::copyPointCloud(*cloud_in, *XYZ);
    pcl::MovingLeastSquaresOMP<pcl::PointXYZ, pcl::PointXYZ>up;
    up.setNumberOfThreads(std::thread::hardware_concurrency());
    //up.setNumberOfThreads(5);
    up.setInputCloud(XYZ);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree;
    up.setSearchMethod(tree);
    up.setSearchRadius(radius * 2.5);
    up.setUpsamplingMethod(pcl::MovingLeastSquaresOMP<pcl::PointXYZ, pcl::PointXYZ>::SAMPLE_LOCAL_PLANE);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filter(new pcl::PointCloud<pcl::PointXYZ>);
    up.setUpsamplingRadius(radius);
    up.setUpsamplingStepSize(step);
    up.process(*cloud_filter);
    PointCloudT::Ptr cloud_final(new PointCloudT);
    pcl::copyPointCloud(*cloud_filter, *cloud_final);
    return cloud_final;
}

//直通滤波去噪----------------------------------------------------------------
PointCloudT::Ptr pcl_statiscal_denoise(PointCloudT::Ptr cloud_in, const double& stddev, const int& num)
{
    pcl::StatisticalOutlierRemoval<PointT>sor;
    sor.setInputCloud(cloud_in);
    sor.setMeanK(num);
    sor.setStddevMulThresh(stddev);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    sor.filter(*cloud_filter);
    return cloud_filter;
}
//直通滤波去噪----------------------------------------------------------------openmp+手写版
PointCloudT::Ptr my_statiscal_denoise(PointCloudT::Ptr cloud_in, const double& stddev, const int& num)
{
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    std::vector<int>indices;
    pcl::KdTreeFLANN<PointT>kdtree;
    kdtree.setInputCloud(cloud_in);
    std::vector<float> distances(cloud_in->size());
    int mythread_nums = std::thread::hardware_concurrency();
    int mean_k = num + 1;
    int valid_distances = 0;
    omp_set_num_threads(mythread_nums);
    #pragma omp parallel reduction(+:valid_distances)
    {
        #pragma omp for
        for (int i = 0; i < cloud_in->size(); i++)
        {
            std::vector<int>pointIdxKNNsearch(num);
            std::vector<float>pointIdxSquaredDistance(num);
            if (kdtree.nearestKSearch(cloud_in->points[i], mean_k, pointIdxKNNsearch, pointIdxSquaredDistance) > 0)
            {
                double dist_sum = 0.0;
                for (int k = 1; k < mean_k; ++k)  // k = 0 is the query point
                    dist_sum += sqrt(pointIdxSquaredDistance[k]);
                #pragma omp critical
                {
                    distances[i] = static_cast<float> (dist_sum / num);
                }
            }
            else
            {
                #pragma omp critical
                {
                    distances[i] = 0.0;
                }
            }
            #pragma omp critical
            {
                valid_distances++;
            }
        }
    }
    double sum = 0, sq_sum = 0;
    for (const auto& dist : distances) {
        sum += dist;
        sq_sum += dist * dist;
    }
    double mean_dist = sum / static_cast<double>(valid_distances);//求出平均距离
    double variance = (sq_sum - sum * sum / static_cast<double>(valid_distances)) / (static_cast<double>(valid_distances) - 1);
    double std_dev = sqrt(variance);
    double distance_threshold = mean_dist + std_dev * stddev;
    for (int i = 0; i < cloud_in->size(); i++)
    {
        if (distances[i]<= distance_threshold)
        {
            indices.push_back(i);
        }
    }     
    pcl::copyPointCloud(*cloud_in, indices, *cloud_filter);
    return cloud_filter;
}
//半径滤波去噪----------------------------------------------------------------
PointCloudT::Ptr pcl_radiussamp(PointCloudT::Ptr cloud_in, const float& radius, const int& num)
{
    pcl::RadiusOutlierRemoval<PointT>ror;
    ror.setInputCloud(cloud_in);
    ror.setRadiusSearch(radius);
    ror.setMinNeighborsInRadius(num);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    ror.filter(*cloud_filter);
    return cloud_filter;
}

PointCloudT::Ptr pcl_gausifilter(PointCloudT::Ptr cloud_in, const float& sigma, const float& sigmarela, const float& threshold, const float& radius)
{
    pcl::filters::GaussianKernel<PointT, PointT>kernel;
    kernel.setSigma(sigma);
    kernel.setThresholdRelativeToSigma(sigmarela);
    kernel.setThreshold(threshold);
    pcl::search::KdTree<PointT>::Ptr kdtree(new pcl::search::KdTree<PointT>);
    kdtree->setInputCloud(cloud_in);
    pcl::filters::Convolution3D< PointT, PointT, pcl::filters::GaussianKernel<PointT, PointT>>convolution;
    convolution.setKernel(kernel);
    convolution.setInputCloud(cloud_in);
    convolution.setNumberOfThreads(std::thread::hardware_concurrency());
    convolution.setSearchMethod(kdtree);
    convolution.setRadiusSearch(radius);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    convolution.convolve(*cloud_filter);
    return cloud_filter;
}

PointCloudT::Ptr my_doublelinefilter(PointCloudT::Ptr cloud_in, const int& nums, const float& sigmadis, const float& sigmanormal)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr XYZ(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::copyPointCloud(*cloud_in, *XYZ);
    PointCloudT::Ptr cloud_filter(new PointCloudT);
    pcl::copyPointCloud(*cloud_in, *cloud_filter);
    pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>);
    tree->setInputCloud(cloud_in);
    std::vector<int>k_indices;//存储索引下标
    std::vector<float>k_distances;//存储距离值
    pcl::PointCloud<pcl::Normal>::Ptr normals_input(computeNormal(XYZ, nums));//计算法向量

    for (size_t i = 0; i < cloud_in->size(); i++)
    {
        float BF = 0;
        float W = 0;
        tree->radiusSearch(i, sigmadis * 2, k_indices, k_distances);
        Eigen::Vector3f normal = normals_input->at(i).getNormalVector3fMap();//获取3维vector存放法向量
        for (size_t j = 0; j < k_indices.size(); j++)
        {
            int id = k_indices[j];
            float dist = sqrt(k_distances[j]);//计算欧氏距离
            Eigen::Vector3f point_p = cloud_in->points[i].getVector3fMap();//获取当前点坐标
            Eigen::Vector3f point_q = cloud_in->points[k_indices[j]].getVector3fMap();//获取当前邻域点坐标
            float normal_dist = normal.dot(point_q - point_p);//点积计算法向量距离
            //计算核函数
            float w_a = (std::exp(-(dist * dist) / (2 * sigmadis * sigmadis)));
            float w_b = (std::exp(-(normal_dist * normal_dist) / (2 * sigmanormal * sigmanormal)));
            float weight = w_a * w_b;
            BF += weight * normal_dist;
            W += weight;
        }
        Eigen::Vector3f point_filter = cloud_in->points[i].getVector3fMap() + (BF / W) * normal;
        cloud_filter->points[i].x = point_filter[0];
        cloud_filter->points[i].y = point_filter[1];
        cloud_filter->points[i].z = point_filter[2];
    }
    return cloud_filter;
}
//最远点采样----------------------------------------------------------------
//PointCloudT::Ptr filterfunc::my_farthestdownsamp(PointCloudT::Ptr cloud_in, const int& sample_nums)
//{
//    if (cloud_in->size() <= sample_nums)
//    {
//        return cloud_in;
//    }
//    else
//    {
//        int num_points = cloud_in->points.size();
//        std::vector<int>sampled_indices;
//        sampled_indices.reserve(sample_nums);
//        std::vector<float> distances_to_selected_points(num_points, std::numeric_limits<float>::max());//用于存储每个点到采样集中最近点的距离，初始值设为float的最大值。
//        // 使用时间生成随机种子
//        std::mt19937 rng(static_cast<unsigned int>(std::time(0)));
//        // 设置随机数生成范围
//        std::uniform_int_distribution<int> dist(0, num_points - 1);
//        int max_index = dist(rng);
//        distances_to_selected_points[max_index] = -1.0;//将第一个采样点的距离设置为-1.0，以表示该点已被选为采样点，不会再被选择。
//        sampled_indices.push_back(max_index);//将该索引加入采样点索引集合
//        int mythread_nums = std::thread::hardware_concurrency();
//        omp_set_num_threads(mythread_nums);
//        for (int i = 1; i < sample_nums; i++)
//        {
//            int farthest_index = -1;
//            float max_distance = -std::numeric_limits<float>::max();
//#pragma omp parallel for
//            for (int j = 0; j < num_points; j++)
//            {
//                if (distances_to_selected_points[j] == -1.0)
//                    continue;
//                float to_mindistance = distances_to_selected_points[j];//计算当前点到采样集中最近点的距离
//                for (int k = 0; k < sampled_indices.size(); k++)
//                {
//                    float dx = cloud_in->points[j].x - cloud_in->points[sampled_indices[k]].x;
//                    float dy = cloud_in->points[j].y - cloud_in->points[sampled_indices[k]].y;
//                    float dz = cloud_in->points[j].z - cloud_in->points[sampled_indices[k]].z;
//                    float distance = dx * dx + dy * dy + dz * dz;
//                    if (distance < to_mindistance)
//                    {
//                        to_mindistance = distance;
//                    }
//                }
//#pragma omp critical
//                {
//                    distances_to_selected_points[j] = to_mindistance;
//                    if (to_mindistance > max_distance)//如果当前距离大于此前所得最大距离，更新最大距离与最远点
//                    {
//                        max_distance = to_mindistance;
//                        farthest_index = j;
//                    }
//                }
//
//            }
//            max_index = farthest_index;
//            distances_to_selected_points[max_index] = -1.0;
//            sampled_indices.push_back(max_index);
//        }
//        PointCloudT::Ptr cloud_final(new PointCloudT);
//        pcl::copyPointCloud(*cloud_in, sampled_indices, *cloud_final);
//        return cloud_final;
//    }
//}


//最远点采样----------------------------------------------------------------
void farest_sampling::sampling()
{
    std::vector<int>sampled_indices;
    if (cloud_in->size() <= sample_nums)
    {
        //emit finished(myid);
        emit finished(myid, sampled_indices);
    }
    else
    {
        int num_points = cloud_in->points.size();
        sampled_indices.reserve(sample_nums);
        std::vector<float> distances_to_selected_points(num_points, std::numeric_limits<float>::max());//用于存储每个点到采样集中最近点的距离，初始值设为float的最大值。
        // 使用时间生成随机种子
        std::mt19937 rng(static_cast<unsigned int>(std::time(0)));
        // 设置随机数生成范围
        std::uniform_int_distribution<int> dist(0, num_points - 1);
        int max_index = dist(rng);
        distances_to_selected_points[max_index] = -1.0;//将第一个采样点的距离设置为-1.0，以表示该点已被选为采样点，不会再被选择。
        sampled_indices.push_back(max_index);//将该索引加入采样点索引集合
        /*int mythread_nums = std::thread::hardware_concurrency();
        omp_set_num_threads(mythread_nums);*/
        int progress_count = 0; // 记录总进度
#pragma omp parallel for schedule(dynamic)
        for (int i = 1; i < sample_nums; i++)
        {
            int farthest_index = -1;
            float max_distance = -std::numeric_limits<float>::max();
            for (int j = 0; j < num_points; j++)
            {
                if (distances_to_selected_points[j] == -1.0)
                    continue;
                float to_mindistance = distances_to_selected_points[j];//计算当前点到采样集中最近点的距离
                for (int k = 0; k < sampled_indices.size(); k++)
                {
                    float dx = cloud_in->points[j].x - cloud_in->points[sampled_indices[k]].x;
                    float dy = cloud_in->points[j].y - cloud_in->points[sampled_indices[k]].y;
                    float dz = cloud_in->points[j].z - cloud_in->points[sampled_indices[k]].z;
                    float distance = dx * dx + dy * dy + dz * dz;
                    if (distance < to_mindistance)
                    {
                        to_mindistance = distance;
                    }
                }
#pragma omp critical
                {
                    distances_to_selected_points[j] = to_mindistance;
                    if (to_mindistance > max_distance)//如果当前距离大于此前所得最大距离，更新最大距离与最远点
                    {
                        max_distance = to_mindistance;
                        farthest_index = j;
                    }
                }

            }
#pragma omp atomic
            progress_count++;
            max_index = farthest_index;
            distances_to_selected_points[max_index] = -1.0;
            sampled_indices.push_back(max_index);
            if (progress_count % (sample_nums / 20) == 0) {
                emit progress((progress_count * 100) / sample_nums); // 发射进度信号
            }
        }
        //PointCloudT::Ptr cloud_filter(new PointCloudT);
        //pcl::copyPointCloud(*cloud_in, sampled_indices, *cloud_filter);
        emit finished(myid, sampled_indices);
        //emit finished(myid);
    }
}

void morph_filter::filtering()
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr XYZ(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::copyPointCloud(*cloud_in, *XYZ);
    pcl::ProgressiveMorphologicalFilter<pcl::PointXYZ> pmf;
    std::vector<int>sampled_indices;
    pmf.setInputCloud(XYZ);
    pmf.setMaxWindowSize(grid);
    pmf.setSlope(slope);
    pmf.setInitialDistance(height_in);
    pmf.setMaxDistance(height_max);
    pmf.extract(sampled_indices);
    emit output(myid, sampled_indices);
}

void morph_filter::my_filtering()
{
    std::vector<float> height_thresholds;
    std::vector<float> window_sizes;
    std::vector<int>groud_indices;
    groud_indices.reserve(cloud_in->points.size());
    int iteration = 0;
    float window_size = 0.0f;
    float height_threshold = 0.0f;
    float cell_size_ = 2.0f;
    for (size_t i = 0; i < cloud_in->points.size(); i++)
    {
        groud_indices.push_back(i);
    }
    while (window_size < grid)//判断窗口大小是否小于阈值
    {
        window_size = cell_size_*(4.0f * (iteration + 1) + 1.0f);
        if (iteration == 0)
            height_threshold = height_in;
        else
            height_threshold = slope * (window_size - window_sizes[iteration - 1]) * cell_size_ + height_in;
        if (height_threshold > height_max)
            height_threshold = height_max;
        window_sizes.push_back(window_size);
        height_thresholds.push_back(height_threshold);
        iteration++;
    } 
    int progress_count = 0; // 记录总进度
    for (std::size_t i = 0; i < window_sizes.size(); ++i)
    {
        PointCloudT::Ptr cloud(new PointCloudT);
        pcl::copyPointCloud<PointT>(*cloud_in, groud_indices, *cloud);
        PointCloudT::Ptr cloud_filter(new PointCloudT);
        pcl::copyPointCloud(*cloud, *cloud_filter);
        pcl::octree::OctreePointCloudSearch<PointT> tree(window_sizes[i]);
        tree.setInputCloud(cloud);
        tree.addPointsFromInputCloud();
        float half_res = window_sizes[i] / 2.0f;//外包大小
        {
            pcl::PointCloud<PointT> cloud_temp;
            pcl::copyPointCloud(*cloud, cloud_temp);
            int mythread_nums = std::thread::hardware_concurrency() / 2 + 1;
            omp_set_num_threads(mythread_nums);
#pragma omp parallel for
            for (int p_idx = 0; p_idx < cloud_temp.points.size(); ++p_idx)//腐蚀
            {
                Eigen::Vector3f bbox_min, bbox_max;
                std::vector<int> pt_indices;
                float minx = cloud_temp.points[p_idx].x - half_res;
                float miny = cloud_temp.points[p_idx].y - half_res;
                float minz = -std::numeric_limits<float>::max();
                float maxx = cloud_temp.points[p_idx].x + half_res;
                float maxy = cloud_temp.points[p_idx].y + half_res;
                float maxz = std::numeric_limits<float>::max();
                bbox_min = Eigen::Vector3f(minx, miny, minz);
                bbox_max = Eigen::Vector3f(maxx, maxy, maxz);
                tree.boxSearch(bbox_min, bbox_max, pt_indices);
                if (!pt_indices.empty())
                {
                    Eigen::Vector4f min_pt, max_pt;
#pragma omp critical
                    {
                        pcl::getMinMax3D<PointT>(cloud_temp, pt_indices, min_pt, max_pt);
                        cloud_filter->points[p_idx].z = min_pt.z();
                    }  
                }
            }
            cloud_temp.swap(*cloud_filter);
            progress_count++;
            emit progress((progress_count * 50) / window_sizes.size()); // 发射进度信号
#pragma omp parallel for
            for (int p_idx = 0; p_idx < cloud_temp.points.size(); ++p_idx)//膨胀
            {
                Eigen::Vector3f bbox_min, bbox_max;
                std::vector<int> pt_indices;
                float minx = cloud_temp.points[p_idx].x - half_res;
                float miny = cloud_temp.points[p_idx].y - half_res;
                float minz = -std::numeric_limits<float>::max();
                float maxx = cloud_temp.points[p_idx].x + half_res;
                float maxy = cloud_temp.points[p_idx].y + half_res;
                float maxz = std::numeric_limits<float>::max();
                bbox_min = Eigen::Vector3f(minx, miny, minz);
                bbox_max = Eigen::Vector3f(maxx, maxy, maxz);
                tree.boxSearch(bbox_min, bbox_max, pt_indices);
                if (!pt_indices.empty())
                {
                    Eigen::Vector4f min_pt, max_pt;
#pragma omp critical
                    {
                        pcl::getMinMax3D<PointT>(cloud_temp, pt_indices, min_pt, max_pt);
                        cloud_filter->points[p_idx].z = max_pt.z();
                    }
                }
            }
        }
        std::vector<int> p_indices;
        for (std::size_t p_idx = 0; p_idx < groud_indices.size(); ++p_idx)
        {
            float diff = cloud->points[p_idx].z - cloud_filter->points[p_idx].z;
            if (diff < height_thresholds[i])
                p_indices.push_back(groud_indices[p_idx]);
        }
        groud_indices.swap(p_indices);
        progress_count++;
        emit progress((progress_count * 50) / window_sizes.size()); // 发射进度信号
        }
    emit output(myid, groud_indices);
}

void pointnet2_semseg::semseg()
{
    int points_num =  cloud_in->points.size();
    
    PointT min_p, max_p;
    pcl::getMinMax3D(*cloud_in, min_p, max_p);
    // 根据点云的范围计算网格数量
    int grid_x = ceil((max_p.x - min_p.x - my_blocksize) / my_blockstride) + 1;
    int grid_y = ceil((max_p.y - min_p.y - my_blocksize) / my_blockstride) + 1;
    std::vector<int> index_room;
    PointCloudT::Ptr cloud_room(new PointCloudT);
    srand(time(0));
    // 遍历所有网格，进行块划分
#pragma omp parallel for
    for (int index_y = 0; index_y < grid_y; index_y++)
    {
        for (int index_x = 0; index_x < grid_x; index_x++)
        {
            float s_x = min_p.x + index_x * my_blockstride;
            float e_x = std::min(s_x + my_blocksize, max_p.x);
            s_x = e_x - my_blocksize;
            float s_y = min_p.y + index_y * my_blockstride;
            float e_y = std::min(s_y + my_blocksize, max_p.y);
            s_y = e_y - my_blocksize;
            std::vector<int> point_idxs;
            for (int i = 0; i < points_num; i++)
            {
                if (cloud_in->points[i].x >= s_x && cloud_in->points[i].x <= e_x && cloud_in->points[i].y >= s_y && cloud_in->points[i].y <= e_y)
                    point_idxs.push_back(i);
            }
            if (point_idxs.empty()) continue;
            int num_batch = ceil(point_idxs.size() * 1.0 / mypointnum);
            int point_size = num_batch * mypointnum;
            std::vector<int> point_idxs_repeat;
            for (int i = 0; i < point_size - point_idxs.size(); i++)
            {
                int id = rand() % point_idxs.size();
                point_idxs_repeat.push_back(point_idxs[id]);
            }
            point_idxs.insert(point_idxs.end(), point_idxs_repeat.begin(), point_idxs_repeat.end());
            // 打乱点索引
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(point_idxs.begin(), point_idxs.end(), g);
            PointCloudT::Ptr cloud_batch(new PointCloudT);
            for (size_t i = 0; i < point_idxs.size(); i++)
            {
                cloud_batch->points.push_back(cloud_in->points[point_idxs[i]]);
            }
            // 数据标准化
#pragma omp critical
            {
                for (int i = 0; i < point_size; i++)
                {
                    //cloud_batch->points[i].m_normal_x = cloud_batch->points[i].x / max_p.x;
                    //cloud_batch->points[i].m_normal_y = cloud_batch->points[i].y / max_p.y;
                    //cloud_batch->points[i].m_normal_z = cloud_batch->points[i].z / max_p.z;
                    cloud_batch->points[i].x -= (s_x + my_blocksize / 2.0);
                    cloud_batch->points[i].y -= (s_y + my_blocksize / 2.0);
                    //cloud_batch->points[i].r /= 255.0;
                    //cloud_batch->points[i].g /= 255.0;
                    //cloud_batch->points[i].b /= 255.0;
                    cloud_room->points.push_back(cloud_batch->points[i]);
                    index_room.push_back(point_idxs[i]);
                }
            }
        }
    }
    
    int n = mypointnum, m = index_room.size() / n;
    std::vector<PointCloudT::Ptr> data_rooms(m);
    for (int i = 0; i < m; ++i) {
        data_rooms[i] = PointCloudT::Ptr(new PointCloudT); // 或者使用 std::make_shared<PointCloudT>();
        data_rooms[i]->resize(n);
    }
    std::vector<std::vector<int>> index_rooms(m, std::vector<int>(n));

    // 将数据划分为批次
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            data_rooms[i]->points[j] = cloud_room->points[i * n + j];
            index_rooms[i][j] = index_room[i * n + j];
        }
    }
    std::vector<std::vector<int>> vote_label_pool(points_num, std::vector<int>(myclassnum, 0));
    int num_blocks = data_rooms.size();
    // 加载PyTorch模型
    torch::jit::script::Module module = torch::jit::load(mymodel_path.toStdString());
    //module.to(torch::kCUDA);
    module.to(torch::kCPU);
    // 对每个批次进行推理
    int progress_count = 0;
#pragma omp parallel for
    for (int sbatch = 0; sbatch < num_blocks; sbatch++)
    {
        PointCloudT::Ptr batch_data= data_rooms[sbatch];
        std::vector<int> point_idx = index_rooms[sbatch];
        std::vector<float> batch(mypointnum * 9);
        for (size_t i = 0; i < mypointnum; i++)
        {
            batch[9 * i + 0] = batch_data->points[i].x;
            batch[9 * i + 1] = batch_data->points[i].y;
            batch[9 * i + 2] = batch_data->points[i].z;
           /* batch[9 * i + 3] = batch_data->points[i].r;
            batch[9 * i + 4] = batch_data->points[i].g;
            batch[9 * i + 5] = batch_data->points[i].b;
            batch[9 * i + 6] = batch_data->points[i].normal_x;
            batch[9 * i + 7] = batch_data->points[i].normal_y;
            batch[9 * i + 8] = batch_data->points[i].normal_z;*/
            batch[9 * i + 3] = 0;
            batch[9 * i + 4] = 0;
            batch[9 * i + 5] = 0;
            batch[9 * i + 6] = batch_data->points[i].x/ max_p.x;
            batch[9 * i + 7] = batch_data->points[i].y / max_p.y;
            batch[9 * i + 8] = batch_data->points[i].z / max_p.z;
        }
        torch::Tensor inputs = torch::from_blob(batch.data(), { 1, mypointnum, 9 }, torch::kFloat);
        inputs = inputs.permute({ 0, 2, 1 });
        inputs = inputs.to(torch::kCPU);
        //inputs = inputs.to(torch::kCUDA);
        auto outputs = module.forward({ inputs }).toTuple();
        torch::Tensor out0 = outputs->elements()[0].toTensor();

        auto max_index = std::get<1>(torch::max(out0, 2));
        max_index = torch::squeeze(max_index).to(torch::kCPU).to(torch::kInt);
        std::vector<int> pred_label(max_index.data_ptr<int>(), max_index.data_ptr<int>() + max_index.numel());
#pragma omp critical
        {
            for (size_t i = 0; i < pred_label.size(); i++)
            {
                int idx = point_idx[i];
                vote_label_pool[idx][pred_label[i]] += 1;
            }
            progress_count++;
            if (progress_count % (num_blocks / 20) == 0) {
                emit progress((progress_count * 100) / num_blocks); // 发射进度信号
            }
        }
    }
  
    std::vector<int>label_indices;
    for (size_t i = 0; i < points_num; i++)
    {
        int max_index = std::max_element(vote_label_pool[i].begin(), vote_label_pool[i].end()) - vote_label_pool[i].begin();
        label_indices.push_back(max_index);
    }
    emit finished(myid, label_indices);
}
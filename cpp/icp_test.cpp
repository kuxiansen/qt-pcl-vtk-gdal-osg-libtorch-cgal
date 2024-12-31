//#include <iostream>
//#include <string>
//#include <pcl/io/ply_io.h>
//#include <pcl/point_types.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <pcl/console/time.h>   // TicToc
//
///**
// * The bool will help us know when the user asks for the next iteration of ICP
// */
//bool next_iteration = false;
//
///**
// * This function takes the reference of a 4x4 matrix and prints the rigid transformation (刚体变换)in an human readable (可读的) way.
// * %6.3f 是指：要输出的浮点数总位数（包括小数点）大于6位的话，按全宽输出，小于 6 位时，小数点后输出3位小数，右对齐，左边不足的位用空格填充
// */
//void print4x4Matrix(const Eigen::Matrix4d& matrix) {
//    printf("Rotation matrix :\n");
//    printf("    | %6.3f %6.3f %6.3f | \n", matrix(0, 0), matrix(0, 1), matrix(0, 2));
//    printf("R = | %6.3f %6.3f %6.3f | \n", matrix(1, 0), matrix(1, 1), matrix(1, 2));
//    printf("    | %6.3f %6.3f %6.3f | \n", matrix(2, 0), matrix(2, 1), matrix(2, 2));
//    printf("Translation vector :\n");
//    printf("t = < %6.3f, %6.3f, %6.3f >\n\n", matrix(0, 3), matrix(1, 3), matrix(2, 3));
//}
//
///**
// * This function is the callback for the viewer.
// * This function will be called whenever a key is pressed when the viewer window is on top.
// * If “space” is hit; set the bool(next_iteration) to true.
// */
//void keyboardEventOccurred(const pcl::visualization::KeyboardEvent& event, void*) {
//    if (event.getKeySym() == "space" && event.keyDown())
//        next_iteration = true;
//}
//
//
//
//int main(int argc, char* argv[]) {
//    // The point clouds we will be using. The 3 point clouds we will use to store the data.
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);  // Original point cloud
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tr(new pcl::PointCloud<pcl::PointXYZ>);  // Transformed point cloud
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_icp(new pcl::PointCloud<pcl::PointXYZ>);  // ICP output point cloud
//
//    // Checking program arguments
//    if (argc < 2) {
//        printf("Usage : ");
//        printf("%s file.ply number_of_ICP_iterations\n", argv[0]);
//        PCL_ERROR("Provide one ply file.\n");
//        return (-1);
//    }
//
//    int iterations = 1;  // Default number of ICP iterations
//    if (argc > 2) {
//        // If the user passed the number of iteration as an argument
//        iterations = atoi(argv[2]);
//        if (iterations < 1) {
//            PCL_ERROR("Number of initial iterations must be >= 1\n");
//            return (-1);
//        }
//    }
//
//    pcl::console::TicToc time;
//    time.tic();
//    if (pcl::io::loadPLYFile("bimba.ply", *cloud_in) < 0) {
//        PCL_ERROR("Error loading cloud %s.\n", argv[1]);
//        return (-1);
//    }
//    std::cout << "\nLoaded file " << argv[1] << " (" << cloud_in->size() << " points) in " << time.toc() << " ms\n" << std::endl;
//
//    // Defining a rotation matrix and translation vector.
//    // We check the arguments of the program, set the number of initial ICP iterations and try to load the PLY file.
//    Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
//
//    // A rotation matrix (see https://en.wikipedia.org/wiki/Rotation_matrix)
//    double theta = M_PI / 8;  // The angle of rotation in radians
//    transformation_matrix(0, 0) = std::cos(theta);
//    transformation_matrix(0, 1) = -sin(theta);
//    transformation_matrix(1, 0) = sin(theta);
//    transformation_matrix(1, 1) = std::cos(theta);
//
//    // A translation on Z axis (0.4 meters) Z 方向平移 0.4 米
//    transformation_matrix(2, 3) = 0.4;
//
//    // Display in terminal the transformation matrix
//    std::cout << "Applying this rigid transformation to: cloud_in -> cloud_icp" << std::endl;
//    print4x4Matrix(transformation_matrix);
//
//    cout << "transformation_matrix\n" << transformation_matrix << endl;
//
//    // Executing the transformation
//    pcl::transformPointCloud(*cloud_in, *cloud_icp, transformation_matrix);
//
//    // We backup cloud_icp into cloud_tr for later use
//    *cloud_tr = *cloud_icp;
//    /**
//     * The Icp(Iterative Closest Point algorithm)
//     * We transform the original point cloud using a rigid matrix transformation.
//     * cloud_in contains the original point cloud.
//     * cloud_tr and cloud_icp contains the translated/rotated point cloud.
//     * cloud_tr is a backup we will use for display (green point cloud).
//     * This is the creation of the ICP object. We set the parameters of the ICP algorithm.
//     * setMaximumIterations(iterations) sets the number of initial iterations to do (1 is the default value).
//     * We then transform the point cloud into cloud_icp.
//     * After the first alignment we set ICP max iterations to 1 for all the next times this ICP object will be used (when the user presses “space”).
//     */
//    time.tic();
//    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//    icp.setMaximumIterations(iterations);
//    icp.setInputSource(cloud_icp);
//    icp.setInputTarget(cloud_in);
//    // 输出配准后点云
//    icp.align(*cloud_icp);
//    // We set this variable to 1 for the next time we will call .align () function
//    icp.setMaximumIterations(1);
//    std::cout << "Applied " << iterations << " ICP iteration(s) in " << time.toc() << " ms" << std::endl;
//
//    /**
//     * Check if the ICP algorithm converged(收敛 ); otherwise exit the program.
//     * In case of success we store the transformation matrix in a 4x4 matrix and then print the rigid matrix transformation.
//     * The reason why we store this matrix is explained later.
//     */
//    if (icp.hasConverged()) {
//        std::cout << "\nICP has converged, score is " << icp.getFitnessScore() << std::endl;
//        std::cout << "\nICP transformation " << iterations << " : cloud_icp -> cloud_in" << std::endl;
//        transformation_matrix = icp.getFinalTransformation().cast<double>();
//        print4x4Matrix(transformation_matrix);
//        cout << "transformation_matrix\n" << transformation_matrix << endl;
//    }
//    else {
//        PCL_ERROR("\nICP has not converged.\n");
//        return (-1);
//    }
//
//    /**
//     * For the visualization we create two viewports in the visualizer vertically separated(垂直分离的可视化器 ).
//     * bckgr_gray_level and txt_gray_lvl are variables to easily switch
//     * from white background & black text/point cloud to black background & white text/point cloud.
//     */
//    pcl::visualization::PCLVisualizer viewer("ICP demo");
//    // Create two vertically separated viewports
//    int v1(0);
//    int v2(1);
//    // 零点在左下角，x轴向右，ｙ轴向上。
//    viewer.createViewPort(0.0, 0.0, 0.5, 1.0, v1);
//    viewer.createViewPort(0.5, 0.0, 1.0, 1.0, v2);
//
//    // The color we will be using -> 背景颜色
//    float bckgr_gray_level = 0.0;  // Black
//    float txt_gray_lvl = 1.0 - bckgr_gray_level;
//
//    /**
//     * We add the original point cloud in the 2 viewports and display it the same color as txt_gray_lvl.
//     * We add the point cloud we transformed using the matrix in the left viewport in green and the point cloud aligned with ICP in red (right viewport).
//     */
//     // Original point cloud is white　-> 原始点云是白色的，旋转45°的绿色的，红色的是icp每一步运行运行的结果。
//     // 这里的txt_gray_lvl，是上面定义的显示界面的文字的颜色，这个颜色和背景颜色相加的和是“１“
//    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_in_color_h(
//        cloud_in,
//        (int)255 * txt_gray_lvl,
//        (int)255 * txt_gray_lvl,
//        (int)255 * txt_gray_lvl);
//    viewer.addPointCloud(cloud_in, cloud_in_color_h, "cloud_in_v1", v1);
//    viewer.addPointCloud(cloud_in, cloud_in_color_h, "cloud_in_v2", v2);
//
//    // Transformed point cloud is green
//    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_tr_color_h(cloud_tr, 20, 180, 20);
//    viewer.addPointCloud(cloud_tr, cloud_tr_color_h, "cloud_tr_v1", v1);
//
//    // ICP aligned point cloud is red
//    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_icp_color_h(cloud_icp, 180, 20, 20);
//    viewer.addPointCloud(cloud_icp, cloud_icp_color_h, "cloud_icp_v2", v2);
//
//    /**
//     * Adding text descriptions in each viewport
//     * We add descriptions for the point clouds in each viewport so the user knows what is what.
//     * The string stream ss is needed to transform the integer iterations into a string.
//     */
//    viewer.addText("White: Original point cloud\nGreen: Matrix transformed point cloud", 10, 15, 16, txt_gray_lvl,
//        txt_gray_lvl, txt_gray_lvl, "icp_info_1", v1);
//    viewer.addText("White: Original point cloud\nRed: ICP aligned point cloud", 10, 15, 16, txt_gray_lvl, txt_gray_lvl,
//        txt_gray_lvl, "icp_info_2", v2);
//
//    std::stringstream ss;
//    ss << iterations;
//    std::string iterations_cnt = "ICP iterations = " + ss.str();
//    viewer.addText(iterations_cnt, 10, 60, 16, txt_gray_lvl, txt_gray_lvl, txt_gray_lvl, "iterations_cnt", v2);
//
//    /**
//     * We set the two viewports background color according to bckgr_gray_level.
//     * To get the camera parameters I simply pressed “C” in the viewer.
//     * Then I copied the parameters into this function to save the camera position / orientation / focal point.
//     * The function registerKeyboardCallback allows us to call a function whenever the users pressed a keyboard key when viewer windows is on top.
//     */
//     // Set background color
//    viewer.setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, v1);
//    viewer.setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, v2);
//
//    // Set camera position and orientation
//    viewer.setCameraPosition(-3.68332, 2.94092, 5.71266, 0.289847, 0.921947, -0.256907, 0);
//    viewer.setSize(1280, 1024);  // Visualiser window size
//
//    // Register keyboard callback :
//    viewer.registerKeyboardCallback(&keyboardEventOccurred, (void*)NULL);
//
//    // Display the visualiser
//    // This is the normal behaviour if no key is pressed. The viewer waits to exit.
//    while (!viewer.wasStopped()) {
//        viewer.spinOnce();
//        /**
//         * If the user press any key of the keyboard, the function keyboardEventOccurred is called;
//         * this function checks if the key is “space” or not. If yes the global bool next_iteration is set to true,
//         * allowing the viewer loop to enter the next part of the code: the ICP object is called to align the meshes.
//         * Remember we already configured this object input/output clouds and we set max iterations to 1 in lines 90-93.
//         *
//         * As before we check if ICP as converged, if not we exit the program. printf(“033[11A”);
//         * is a little trick to go up 11 lines in the terminal to write over the last matrix displayed.
//         * In short it allows to replace text instead of writing new lines; making the output more readable.
//         * We increment iterations to update the text value in the visualizer.
//
//         * Now we want to display the rigid transformation from the original transformed point cloud to the current alignment made by ICP.
//         * The function getFinalTransformation() returns the rigid matrix transformation done during the iterations (here: 1 iteration).
//         * This means that if you have already done 10 iterations this function returns the matrix to transform the point cloud from the iteration 10 to 11.
//
//         * This is not what we want.
//         * If we multiply the last matrix with the new one the result is the transformation matrix from the start to the current iteration.
//         * This is basically how it works
//         *
//         * While this is mathematically true, you will easily notice that this is not true in this program due to roundings.
//         * This is why I introduced the initial ICP iteration parameters.
//         * Try to launch the program with 20 initial iterations and save the matrix in a text file.
//         * Launch the same program with 1 initial iteration and press space till you go to 20 iterations.
//         * You will a notice a slight difference.
//         * The matrix with 20 initial iterations is much more accurate than the one multiplied 19 times.
//         */
//         // The user pressed "space" :
//        if (next_iteration) {
//            // The Iterative Closest Point algorithm
//            time.tic();
//            icp.align(*cloud_icp);
//            std::cout << "Applied 1 ICP iteration in " << time.toc() << " ms" << std::endl;
//
//            if (icp.hasConverged()) {
//                printf("\033[11A");  // Go up 11 lines in terminal output.
//                printf("\nICP has converged, score is %+.0e\n", icp.getFitnessScore());
//                std::cout << "\nICP transformation " << ++iterations << " : cloud_icp -> cloud_in" << std::endl;
//                transformation_matrix *= icp.getFinalTransformation().cast<double>();  // WARNING /!\ This is not accurate! For "educational" purpose only!
//                print4x4Matrix(transformation_matrix);  // Print the transformation between original pose and current pose
//
//                ss.str("");
//                ss << iterations;
//                std::string iterations_cnt = "ICP iterations = " + ss.str();
//                viewer.updateText(iterations_cnt, 10, 60, 16, txt_gray_lvl, txt_gray_lvl, txt_gray_lvl, "iterations_cnt");
//                viewer.updatePointCloud(cloud_icp, cloud_icp_color_h, "cloud_icp_v2");
//            }
//            else {
//                PCL_ERROR("\nICP has not converged.\n");
//                return (-1);
//            }
//        }
//        // We set the bool to false and the rest is the ending of the program.
//        next_iteration = false;
//    }
//    return (0);
//}

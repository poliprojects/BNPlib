#include <iostream>
#include <fstream>
#include <chrono>

#include "math.h"

int main(int argc, char *argv[]){
    // PROVA TIME //////////////////////

    Eigen::MatrixXd data_test = Eigen::MatrixXd::Identity(1000,2);
    // Try different cluster sizes
    //std::vector<int> clust_idxs = {0,30,20,5,6,7,22,34,91,4,35,76,90};
    //std::vector<int> clust_idxs = {1};
    std::vector<int> clust_idxs(1000); 

    std::iota(std::begin(clust_idxs), std::end(clust_idxs), 0); // fill with
        //0, 1, ..., 999

    std::chrono::time_point<std::chrono::system_clock> start_1, end_1,
        start_2, end_2;

    // FIRST WAY 
    start_1 = std::chrono::system_clock::now();    
    Eigen::MatrixXd curr_data(data_test.rows(), data_test.cols());
    int k=0;

    for(auto &idx : clust_idxs){
        curr_data.row(k) = data_test.row(idx);  
        k += 1;
    }
    
    curr_data.conservativeResize(k,Eigen::NoChange); 
    end_1 = std::chrono::system_clock::now();
    

    // SECOND WAY
    start_2 = std::chrono::system_clock::now();    
    std::vector<double> curr_data_2;

        for(auto &idx : clust_idxs){
            for (int i = 0; i < data_test.cols(); i++){
                curr_data_2.push_back(data_test.row(idx)(i));
            }
  }
  

            
    Eigen::MatrixXd curr_data_map = Eigen::Map<Eigen::MatrixXd>(
        curr_data_2.data(), curr_data_2.size()/data_test.cols(),
        data_test.cols()); 
    end_2 = std::chrono::system_clock::now();

    typedef std::chrono::duration<int, std::ratio<1, 100000000>> shakes;

    int elapsed_seconds_1 = std::chrono::duration_cast<shakes>(
        end_1-start_1).count();
    int elapsed_seconds_2 = std::chrono::duration_cast<shakes>(
        end_2-start_2).count();
    std::cout << "first  time:" << elapsed_seconds_1 << std::endl;
    std::cout << "second time:" << elapsed_seconds_2 << std::endl;
    return 0;
}

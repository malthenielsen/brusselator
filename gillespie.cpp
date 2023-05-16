#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <random>
#include <algorithm>
#include <omp.h>
#include <fstream>
#include <cmath>
#include <set>
#include <utility>
#include <cstdio>
using namespace std;

void delete_file(const char* filename) {
    if (std::remove(filename) != 0) {
        std::perror("Error deleting file");
    }
}

std::vector<int> findPeaks(const std::vector<double>& data, int windowSize, int distance) {
    std::vector<double> smoothedData;
    for (int i = 0; i < data.size(); i++) {
        int start = std::max(0, i - windowSize);
        int end = std::min(static_cast<int>(data.size()), i + windowSize + 1);
        double sum = 0.0;
        for (int j = start; j < end; j++) {
            sum += data[j];
        }
        double average = sum / (end - start);
        smoothedData.push_back(average);
    }

    std::vector<int> peaks;
    int prev_idx = 0;
    for (int i = 1; i < smoothedData.size() - 1; i++) {
        if (smoothedData[i] > smoothedData[i - 1] && smoothedData[i] > smoothedData[i + 1] && i - 0 > distance) {
            peaks.push_back(i);
            prev_idx = i;
        }
    }

    return peaks;
}


std::pair<double, double> f1(double x1,double x2, double a, double b, double c, double d){
  double resp = a + c*x1*x1*x2;
  double resn =  (b + d)*x1;
  return std::make_pair(resp, resn);
}

std::pair<double, double> f2(double x1,double x2, double b, double c){
  double resp = b * x1;
  double resn = c*x1*x1*x2;
  return std::make_pair(resp, resn);
}

int single_run(){
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> distrib(0, 1);

  double tau0 = .1;
  double a = 3;
  double c = distrib(gen)*5;
  double d = distrib(gen)*10;
  double b = (d + c*a*a/(d*d) - tau0);
  double RealTime = 0;
  double vol = 1.0*pow(10,-15);
  double NA = 6.03*pow(10,23);
  double cal = NA*vol*pow(10,-6);

  double x1 = a/d *cal ;
  double x2 = b*d/(c*a)*cal ;

  a = a * cal;
  c = c / (cal * cal);
  b = b;
  d = d;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;

  int max_step = 1000000;
  vector<double> X1 = vector<double>(max_step,0);
  vector<double> X2 = vector<double>(max_step,0);
  vector<double> RT = vector<double>(max_step,0);


  std::cout << x1 << "  " << x2 << std::endl;

  int idx = 0;
  int idx_vec = 0;

  std::pair<double, double> res1 = f1(x1, x2, a, b, c, d);
  std::pair<double, double> res2 = f2(x1, x2, b, c);

  vector<double> dt = vector<double>(4,0);
  double DT = 0;
  float prob = 0;
  float sum = 0;

  for (int i =0; i < max_step*10; i++){
     res1 = f1(x1, x2, a, b, c, d);
     res2 = f2(x1, x2, b, c);
     dt[0] = res1.first;
     dt[1] = res1.second;
     dt[2] = res2.first;
     dt[3] = res2.second;
     double total_rate = dt[0] + dt[1] + dt[2] + dt[3];
     DT = 1/total_rate * log(1/distrib(gen));
     prob = distrib(gen);
     dt[0] = dt[0]/total_rate;
     dt[1] = dt[1]/total_rate;
     dt[2] = dt[2]/total_rate;
     dt[3] = dt[3]/total_rate;
     int max_iter = 0;
     sum = 0;
     for (int j  = 0; j < 4; j++){
       sum = sum + dt[j];
       if (sum > prob){
         max_iter = j;
         break;
       }
     }
     if (max_iter == 0){
       x1 += 1;
     }else if (max_iter == 1){
       x1 -= 1;
     }else if (max_iter == 2){
       x2 += 1;
     }else if (max_iter == 3){
       x2 -= 1;
     }
     RealTime += DT;

     if (idx%10 == 0){ 
       X1[idx_vec] = x1;
       X2[idx_vec] = x2;
       RT[idx_vec] = RealTime;
       idx_vec += 1;
     }
    idx += 1;
  
  }
  delete_file("peaks.bin");
  std::ofstream outPeak("peaks.bin", std::ios::binary | std::ios::app);
  std::vector<int> peakInds = findPeaks(X1, 200, 1000);
  outPeak.write(reinterpret_cast<const char*>(peakInds.data()), peakInds.size() * sizeof(int));
  outPeak.close();

  delete_file("results.bin");
  std::ofstream outFile("results.bin", std::ios::binary | std::ios::app);
  outFile.write(reinterpret_cast<const char*>(X1.data()), X1.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(X2.data()), X2.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(RT.data()), RT.size() * sizeof(double));
  outFile.close();

  return 1;

}




int main(){

  int res = single_run();

  return 0;

}

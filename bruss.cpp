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
using namespace std;


long double f1(long double x1,long double x2, long double a, long double b, long double c, long double d, long double osc){
  long double res = a *(1 + osc) - (b + d) * x1 + c*x1*x1*x2;
  return res;
}

long double f2(long double x1,long double x2, long double a, long double b, long double c, long double d, long double osc){
  long double res = b * x1 - c*x1*x1*x2;
  return res;
}

vector<int> find_peaks(const vector<double long>& input_vec){
  vector<int> output_vec;
  int vector_size = input_vec.size();
  for (int i = 1; i < vector_size-1; i++){
    if (input_vec[i] > input_vec[i-1]){
      if (input_vec[i] > input_vec[i+1]){
        output_vec.push_back(i);
      }
    }
  }
  return output_vec;
}

float mean(const vector<int>& input_vec){
  float sum;
  int vector_size = input_vec.size();
  for (int i = 0; i < vector_size; i++){
    sum = sum + (float)input_vec[i];
  }
  sum = sum / static_cast<float>(vector_size);
  return sum;
}

float STD(const vector<int>& input_vec){
  float mean_val = mean(input_vec);
  int vector_size = input_vec.size();
  float diff = 0;
  float sum = 0;
  for (int i = 0; i < vector_size; i++){
    diff = (float)input_vec[i] - mean_val;
    sum += diff*diff;
  }
  sum = sum / input_vec.size();
  sum = std::sqrt(sum);
  return sum;
}

vector<int> diff_func(const vector<int>& input){
  int vector_size = input.size();
  vector<int> output_vec;
  output_vec.resize(vector_size-2, 0);
  for (int i = 0; i < vector_size-1; i++){
    output_vec[i] = input[i + 1] - input[i];
  }
  return output_vec;
}


int main() {

  int N = 200;
  vector<float> D;
  D.resize(N,0);
  vector<float> amp;
  amp.resize(N,0);

  // #pragma omp parallel for num_threads(7)
  for (int idx = 0; idx < N; idx++) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(0, 1);
    long double tau0 = .1;
    long double a = 3;
    long double c = distrib(gen)*5;
    long double d = distrib(gen)*10;
    long double b = (d + c*a*a/(d*d) - tau0);

    long double A1 = 0.01;
    long double omega = 0.86;
    long double h  = .01;
    long double RT = 0;
    std::cout << b - d << std::endl;
    float diff = b - d;
    D[idx] = diff;
    
    long double Tmax = 500;
    int steps = Tmax/h;
    vector<long double> X1;
    X1.resize(steps,0);
    vector<long double> X2;
    X2.resize(steps,0);
    vector<long double> osc;
    osc.resize(steps,0);

    long double x1 = a/d;
    long double x2 = b*d/(c*a);
    // std::cout << x1 << "  " << x2 << std::endl;
    
    // std::ofstream outfile("results.txt");

    for (int i = 0; i < steps; i++){

      long double osc = A1*sin(RT*omega);

      long double l1 = f1(x1, x2, a, b, c, d, osc)*h;
      long double k1 = f2(x1, x2, a, b, c, d, osc)*h;

      long double l2 = f1(x1 + .5*l1, x2 + .5*k1, a, b, c, d, osc)*h;
      long double k2 = f2(x1 + .5*l1, x2 + .5*k1, a, b, c, d, osc)*h;

      long double l3 = f1(x1 + .5*l2, x2 + .5*k2, a, b, c, d, osc)*h;
      long double k3 = f2(x1 + .5*l2, x2 + .5*k2, a, b, c, d, osc)*h;

      long double l4 = f1(x1 + l3, x2 + k3, a, b, c, d, osc)*h;
      long double k4 = f2(x1 + l3, x2 + k3, a, b, c, d, osc)*h;

      x1 = x1 + 1./6*(l1 + 2*l2 + 2*l3 + l4);
      x2 = x2 + 1./6*(k1 + 2*k2 + 2*k3 + k4);
      RT = RT + h;


      X1[i] = x1;
      X2[i] = x2;
      // outfile << x1 << "," << x2 << "\n";

    }
    // outfile.close();

    auto max_it = std::max_element(X2.begin(), X2.end(), [](const float& a, const float& b) { return a < b; });
    float max_val = *max_it;
    amp[idx] = max_val;
    vector<int> arr = find_peaks(X1);
    vector<int> diff_arr = diff_func(arr);
    // std::set<int> freq_set(diff_arr.begin(), diff_arr.end());
    std::cout << mean(diff_arr) << "  " << STD(diff_arr) << std::endl;

  }

  // std::ofstream outfile("scatter.txt");

  // for (int h = 0; h < N; h++) {
    // outfile << D[h] << "," << amp[h] << "\n";
  // }
  // outfile.close();

  

  // std::cout << max_val << std::endl;
}






  // std::cout << steps << std::endl;








  // vector<int> W;
  // W.resize(N,0);
  // for (int i = 0; i < steps; i++){
  // }



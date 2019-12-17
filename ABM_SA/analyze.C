#include <stdio.h>
#include <new>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sstream>

using namespace std;

int main ( int argc, char *argv[] ){
  int N,K;
  double x;
  if ( 1 < argc ){N = atoi( argv[1] );}else{cout << "Enter the number of sample points:" << endl;cin >> N;}
  if ( 2 < argc ){K = atoi( argv[2] );}else{cout << "Enter the number of parameters to be tested by the sensitivity:" << endl;cin >> K;}
  double *YA = new double[N];
  double *YB = new double[N];
  double **YC = new double*[K];
  for(int i = 0; i < K; i++)
    YC[i] = new double[N];
  ifstream input;
  input.open ("results-A.txt");
  if (input.is_open()) {
    for(int i = 0; i < N; i++){
      for(int j = 0; j < K; j++){
	input >> x;
      }
      input >> YA[i];
    }
  }
  input.close();
  input.open ("results-B.txt");
  if (input.is_open()) {
    for(int i = 0; i < N; i++){
      for(int j = 0; j < K; j++){
	input >> x;
      }
      input >> YB[i];
    }
  }
  input.close();
  std::string base_name = "results-C";
  for(int p = 0; p < K; p++){
    std::stringstream oss;
    oss << base_name << p << ".txt";
    input.open(std::string(oss.str()).c_str());
    if (input.is_open()) {
      for(int i = 0; i < N; i++){
	for(int j = 0; j < K; j++){
	  input >> x;
	}
	input >> YC[p][i];
      }
    }
    input.close();
  }
  ofstream output;
  std::string base_n = "nft-N";
  std::stringstream os;
  os << base_n << N << ".txt";
  output.open(std::string(os.str()).c_str());
  for(int n = 1; n <= N; n++){
    output << scientific << n;
    double yaa = 0., ya2 = 0.;
    for(int i = 0; i < n; i++){
      ya2 += YA[i];
      yaa += YA[i]*YA[i];
    }
    ya2 = pow(ya2/n,2);
    yaa = yaa/n;
    double vy = yaa-ya2;
    for(int k = 0; k < K; k++){
      double si = 0.,st = 0.;
      for(int i = 0; i < n; i++){
	si += YB[i]*(YC[k][i]-YA[i]);
	st += YA[i]*(YA[i]-YC[k][i]);
      }
      si = si/n;
      st = st/n;
      output << scientific << " " << si/vy << " " << st/vy;
    }
    output << endl;
  }
  output.close();
  for(int k = 0; k < K; k++)
    delete[] YC[k];
  delete[] YC;
  delete[] YA;
  delete[] YB;
  return 0;
}

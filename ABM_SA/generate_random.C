#include <iostream>
#include <new>
#include "nr3.h"
#include "ran.h"
#include "matriz.h"
#include <sstream>

using namespace std;

int main ( int argc, char *argv[] ){
  Ran ran(7);
  int N,K;
  if ( 1 < argc ){N = atoi( argv[1] );}else{cout << "Enter the number of sample points:" << endl;cin >> N;}
  if ( 2 < argc ){K = atoi( argv[2] );}else{cout << "Enter the number of parameters to be tested by the sensitivity:" << endl;cin >> K;}
  Matrix A(N,K),B(N,K);
  for(int i = 0; i < A.LC(0); i++)
    for(int j = 0; j < A.LC(1); j++)
      A[i][j] = ran.doub();
  for(int i = 0; i < B.LC(0); i++)
    for(int j = 0; j < B.LC(1); j++)
      B[i][j] = ran.doub();
  ofstream output;
  output.open ("files-A.txt");
  for(int i = 0; i < A.LC(0); i++){
    for(int j = 0; j < A.LC(1); j++)
      output << scientific << A[i][j] << " ";
    output << endl;
  }
  output.close();
  output.open ("files-B.txt");
  for(int i = 0; i < B.LC(0); i++){
    for(int j = 0; j < B.LC(1); j++)
      output << scientific << B[i][j] << " ";
    output << endl;
  }
  output.close();
  std::string base_name = "files-C";
  for(int p = 0; p < K; p++){
    std::stringstream oss;
    oss << base_name << p << ".txt";
    output.open(std::string(oss.str()).c_str());
    for(int i = 0; i < B.LC(0); i++){
      for(int j = 0; j < B.LC(1); j++){
	if(j==p)
	  output << scientific << B[i][j] << " ";
	else
	  output << scientific << A[i][j] << " ";
      }
      output << endl;
    }
    output.close();
  }
  return 0;
}

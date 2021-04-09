// Testing Jabc for random states
#include "../include/Jabc.hpp"
#include <armadillo>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
  int a = atoi(argv[1]);
  int b = atoi(argv[2]);
  int c = atoi(argv[3]);
  int n_q = atoi(argv[4]);

  cx_mat psi_r, psi_i, psi;
  psi_r.randn(1<<n_q, 1);
  psi_i.randn(1<<n_q, 1);
  psi = psi_r + cx_double(0.0, 1.0) * psi_i;
  
  
  psi /= norm(psi);

  cout << Jabc(a,b,c,psi)<< endl;
}

#include "../include/Jabc.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>


using namespace std;
using namespace arma;


void transform_xlogx_mat(cx_mat &X)
{
  cx_mat U, V;
  vec s, sp;

  svd_econ(U, s, V, X);

  sp = s;

  for (int i=0; i<s.size(); i++)
    {
      sp(i) = -2.0*s(i) * log(s(i));
    }

  transform_xlogx_vec(s);
  X = U * diagmat(s) * V.t();
}

void transform_xlogx_vec(vec &s)
{
  for (int i=0; i<s.size(); i++)
    {
      s(i) = -2.0*s(i) * log(s(i));
    }
}

int log2_int(int n)
{
  int k=-1;
  while (n>0)
    {
      n/=2;
      k++;
    }
  return k;
}

// |\psi\rangle \to K_X|\psi\rangle, where X is the first k qubits.
// (Note that this is equal to K_{\bar{X}})|\psi\rangle).
// n: Number of qubits
// Input format: psi should be a column vector.
void apply_modular_op(int k, int n, cx_mat &psi)
{
  psi.reshape(1<<k, 1<<(n-k));
  transform_xlogx_mat(psi);
  psi.reshape(1<<n,1);
}


// Given |\psi>_{ABCD}, return \log \rho_{AB} |\psi>_{ABCD}
// a, b, c: Number of bits in A, B, C.
// Format: First a bits are A. Next b bits are B. The next c bits are C.
void transform_ab(int a, int b, int c, cx_mat &psi)
{
  int dim = psi.size();
  int n = log2_int(dim);
  apply_modular_op(a+b, n, psi);
}

// Given |\psi>_{ABCD}, return \log \rho_{BC} |\psi>_{ABCD}
// a, b, c: Number of bits in A, B, C.
// Format: First a bits are A. Next b bits are B. The next c bits are C.
void transform_bc(int a, int b, int c, cx_mat &psi)
{
  int dim = psi.size();
  int n = log2_int(dim);
  psi.reshape(1<<a, 1<<(n-a));
  psi = psi.st();
  psi.reshape(1<<n,1);
  apply_modular_op(b+c, n, psi);
  psi.reshape(1<<(n-a), 1<<a);
  psi = psi.st();
  psi.reshape(1<<n,1);
}


double Jabc(int a, int b, int c, cx_mat &psi)
{
  cx_mat psi_cpy = psi;

  transform_ab(a,b,c,psi);  
  transform_bc(a,b,c,psi_cpy);
  
  cx_mat result = psi_cpy.t() * psi;

  cout << "result" << endl;
  cout << result << endl;
  
  cx_double inner_product = result(0,0);
  return inner_product.imag() / 2;  
}

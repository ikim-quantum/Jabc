#ifndef Jabc_HPP_
#define Jabc_HPP_

#include <armadillo>
using namespace arma;

void transform_xlogx_mat(cx_mat &X);
void transform_xlogx_vec(vec &s);
int log2_int(int n);
void apply_modular_op(int k, int n, cx_mat &psi);
void transform_ab(int a, int b, int c, cx_mat &psi);
void transform_bc(int a, int b, int c, cx_mat &psi);
double Jabc(int a, int b, int c, cx_mat &psi);

#endif // Jabc_HPP_

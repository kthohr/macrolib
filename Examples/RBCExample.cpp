/*
 * Solving a basic RBC model using gensys++
 *
 * Keith O'Hara
 * 06/29/15
 *
 * This version:
 * 06/30/15
 *
 * For example, on a Mac you can compile an executable using
 *
 * g++ -I/opt/local/include/ gensys.cpp -c -o gensys.o -O2
 * g++ -I/opt/local/include/ RBCExample.cpp -c -o RBCExample.o -O2
 * g++ -o rbc RBCExample.o gensys.o -lm -framework Accelerate
 *
 * Replace '-framework Accelerate' with '-lblas -llapack' on a Linux-based OS
 * (or, if available, you can substitute '-lblas' and '-llapack' with high-speed
 * replacements, such as OpenBLAS: '-lopenblas -llapack').
 *
 */
#include <iostream>
#include "armadillo"
#include "armadillo_qz"
#include "gensys.h"
//
using namespace arma;
using namespace std;
//
int main(int argc, char** argv)
{
    //
    double beta      = 0.99;
    double alpha     = .33;
    double delta     = .015;
    double eta       = 1.0;
    double rho       = 0.95;
    double sigmaTech = 1;
    //
    double RSS  = 1/beta;
    double YKSS = (RSS + delta - 1)/alpha;
    double IKSS = delta;
    double IYSS = ((alpha*delta)/(RSS + delta - 1));
    double CYSS = 1 - IYSS;
    //
    double Gam62 = alpha*YKSS/(RSS);
    double Gam63 = alpha*YKSS/(RSS);
    //
    arma::mat Gamma0(9,9);
    arma::mat Gamma1(9,9);
    //
    Gamma0 <<       1 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<      -1 <<   1/eta << endr
           <<   -CYSS <<       1 <<       0 <<       0 <<       0 <<   -IYSS <<       0 <<       0 <<       0 << endr
           <<       0 <<       1 <<       0 <<-1+alpha <<       0 <<       0 <<      -1 <<       0 <<       0 << endr
           <<       0 <<       0 <<       1 <<       0 <<       0 <<   -IKSS <<       0 <<       0 <<       0 << endr
           <<    -eta <<       1 <<       0 <<      -1 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<   Gam62 <<       0 <<       0 <<      -1 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       1 <<       0 <<       0 << endr
           <<       1 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       1 <<       0 <<       0 <<       0 <<       0 << endr;
    //
    Gamma1 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<   alpha <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 << 1-delta <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<   Gam63 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<     rho <<       0 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       1 <<       0 << endr
           <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       0 <<       1 << endr;
    //
    arma::mat C = arma::zeros<arma::mat>(9,1);
    //
    arma::mat Psi = arma::zeros<arma::mat>(9,1);
    Psi(6,0) = 1;
    //
    arma::mat Pi = arma::zeros<arma::mat>(9,2);
    Pi(7,0) = 1; Pi(8,1) = 1;
    //
    arma::mat G1; arma::mat Cons; arma::mat impact;
    arma::cx_mat fmat; arma::cx_mat fwt; arma::cx_mat ywt; arma::cx_mat gev;
    arma::vec eu; arma::mat loose;
    //
    gensys(Gamma0, Gamma1, C, Psi, Pi, G1, Cons, impact, fmat, fwt, ywt, gev, eu, loose);
    //
    arma::cout << G1 << arma::endl;
    arma::cout << Cons << arma::endl;
    arma::cout << impact << arma::endl;
    arma::cout << fmat << arma::endl;
    arma::cout << fwt << arma::endl;
    arma::cout << ywt << arma::endl;
    arma::cout << gev << arma::endl;
    arma::cout << eu << arma::endl;
    arma::cout << loose << arma::endl;
    //
    return 0;
}
//
//
//END
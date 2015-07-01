/*
 * Solving a simple 3-equation New-Keynesian Model using gensys++
 *
 * Keith O'Hara
 * 06/29/15
 *
 * This version:
 * 06/30/15
 *
 * For example, on a Mac you can compile an executable form using
 *
 * g++ -I/opt/local/include/ gensys.cpp -c -o gensys.o -O2
 * g++ -I/opt/local/include/ NKExample.cpp -c -o NKExample.o -O2
 * g++ -o nkm NKExample.o gensys.o -lm -framework Accelerate
 *
 * Replace '-framework Accelerate' with '-lblas -llapack' on a Linux-based OS
 * (or, if available, you can substitute '-lblas' and '-llapack' with high-speed
 * replacements, such as OpenBLAS: '-lopenblas -llapack').
 *
 */
#include <iostream>
#include "armadillo_qz"
#include "gensys.h"
//
using namespace arma;
using namespace std;
//
int main(int argc, char** argv)
{
    //
    double sigma     = 1;
    double kappa     = 0.2;
    double beta      = 0.99;
    double phipi     = 1.5;
    double phix      = 0.5;
    double rhog      = 0.5;
    double rhou      = 0.3;
    double sigmai    = 0.0025;
    double sigmag    = 0.0025;
    double sigmau    = 0.0025;
    //
    arma::mat Gamma0(9,9);
    arma::mat Gamma1(9,9);
    //
    Gamma0 <<       1  <<   sigma  <<       0  <<      -1  <<       0  <<      -1  <<  -sigma  << endr
           <<  -kappa  <<       0  <<       1  <<       0  <<      -1  <<       0  <<   -beta  << endr
           <<   -phix  <<       1  <<  -phipi  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       1  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       1  <<       0  <<       0  << endr
           <<       1  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       1  <<       0  <<       0  <<       0  <<       0  << endr;
    //
    Gamma1 <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<    rhog  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<    rhou  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  << endr;
    //
    arma::mat C = arma::zeros<arma::mat>(7,1);
    //
    arma::mat Psi = arma::zeros<arma::mat>(7,3);
    Psi(2,2) = sigmai;
    Psi(3,0) = sigmag;
    Psi(4,1) = sigmau;
    //
    arma::mat Pi = arma::zeros<arma::mat>(7,2);
    Pi(5,0) = 1;
    Pi(6,1) = 1;
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
    arma::cout << eu << arma::endl;
    /*
    arma::cout << fmat << arma::endl;
    arma::cout << fwt << arma::endl;
    arma::cout << ywt << arma::endl;
    arma::cout << gev << arma::endl;
    arma::cout << loose << arma::endl;*/
    //
    return 0;
}
//
//
//END
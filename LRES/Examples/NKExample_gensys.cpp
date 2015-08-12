/*
 * Solving a simple 3-equation New-Keynesian Model using gensys++
 *
 * Keith O'Hara
 * 06/29/15
 *
 * This version:
 * 07/15/15
 *
 * For example, on a Mac (depending on where you installed the Armadillo header files) 
 * you can compile an executable using
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
#include "armadillo"
#include "gensys.h"
//
using namespace arma;
using namespace std;
//
int main(int argc, char** argv)
{
    //
    double alpha    = 0.33;
    double beta     = 0.99;
    double vartheta = 6;
    double theta    = 0.6667;
    //
    double eta      = 1;
    double phi      = 1;
    double phi_pi   = 1.5;
    double phi_y    = 0.125;
    double rho_a    = 0.90;
    double rho_v    = 0.5;
    //
    double BigTheta = (1-alpha)/(1-alpha+alpha*vartheta);
    double kappa = (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)));
    double psi = (eta*(1+phi))/(1-alpha+eta*(phi + alpha));
    //
    double G0_47 = (1/eta)*psi*(rho_a - 1);
    //
    arma::mat Gamma0(10,10);
    arma::mat Gamma1(10,10);
    //
    Gamma0 <<      -1  <<       0  <<       0  <<     eta  <<  -eta/4  <<       0  <<       0  <<       0  <<       1  <<   eta/4  << endr
           <<   kappa  <<       0  <<   -0.25  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<  beta/4  << endr
           <<   phi_y  <<       0  <<phi_pi/4  <<       0  <<   -0.25  <<       0  <<       0  <<       1  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<      -1  <<       0  <<       0  <<   G0_47  <<       0  <<       0  <<       0  << endr
           <<       0  <<      -1  <<       0  <<       0  <<       0  << 1-alpha  <<       1  <<       0  <<       0  <<       0  << endr
           <<      -1  <<       1  <<       0  <<       0  <<       0  <<       0  <<    -psi  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  <<       0  <<       0  << endr
           <<       1  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       1  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr;
    //
    Gamma1 <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<   rho_a  <<       0  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<   rho_v  <<       0  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  <<       0  << endr
           <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       0  <<       1  << endr;
    //
    arma::mat C = arma::zeros<arma::mat>(10,1);
    //
    arma::mat Psi = arma::zeros<arma::mat>(10,2);
    Psi(6,0) = 1;
    Psi(7,1) = 1;
    //
    arma::mat Pi = arma::zeros<arma::mat>(10,2);
    Pi(8,0) = 1;
    Pi(9,1) = 1;
    /*
    arma::cout << Gamma0 << arma::endl;
    arma::cout << Gamma1 << arma::endl;
    arma::cout << Psi << arma::endl;
    arma::cout << Pi << arma::endl;
    */
    arma::mat G1; arma::mat Cons; arma::mat impact;
    arma::cx_mat fmat; arma::cx_mat fwt; arma::cx_mat ywt; arma::cx_mat gev;
    arma::vec eu; arma::mat loose;
    //
    gensys(G1, Cons, impact, fmat, fwt, ywt, gev, eu, loose, Gamma0, Gamma1, C, Psi, Pi);
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
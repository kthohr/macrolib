#ifndef gensys_H
#define gensys_H

#include <iostream>
#include <stdlib.h>
#include "armadillo"

int gensys(arma::mat Gamma0, arma::mat Gamma1, arma::mat C, arma::mat Psi, arma::mat Pi, arma::mat & G1, arma::mat & Cons, arma::mat & impact, arma::cx_mat & fmat, arma::cx_mat & fwt, arma::cx_mat & ywt, arma::cx_mat & gev, arma::vec & eu, arma::mat & loose);

#endif



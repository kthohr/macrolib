/*
 * Uhlig's method for solving linear rational
 * expectations models using a generalized eigen
 * decomposition.
 *
 * Written by
 * Keith O'Hara
 * 07/01/12
 *
 * This version:
 * 07/03/15
 *
 * Requires the 'Armadillo' linear algebra library.
 */
#ifndef uhlig_H
#define uhlig_H

#include <stdlib.h>
#include "armadillo"

int uhlig(arma::mat &P, arma::mat &Q, arma::mat &R, arma::mat &S,
          arma::cx_vec &eigenvals, arma::cx_mat &eigenvecs,
          arma::mat A, arma::mat B, arma::mat C, arma::mat D, arma::mat F, arma::mat G, arma::mat H,
          arma::mat J, arma::mat K, arma::mat L, arma::mat M, arma::mat N,
          arma::vec whichEig);

#endif



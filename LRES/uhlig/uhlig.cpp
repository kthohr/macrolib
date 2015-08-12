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

#include <stdlib.h>
#include "armadillo"

int uhlig(arma::mat &P, arma::mat &Q, arma::mat &R, arma::mat &S,
          arma::cx_vec &eigenvals, arma::cx_mat &eigenvecs,
          arma::mat A, arma::mat B, arma::mat C, arma::mat D,
          arma::mat F, arma::mat G, arma::mat H,
          arma::mat J, arma::mat K, arma::mat L, arma::mat M, arma::mat N,
          arma::vec which_eig)
{
    //
    int l = C.n_rows;
    int n = C.n_cols;
    int k = std::min(N.n_rows,N.n_cols);
    //
    int m = 0;
    if( l == 0 ){
        m = F.n_cols;
    }else{ // l > 0
        m = A.n_cols;
    }
    //
    int pick_eig = 0;
    if(which_eig.n_elem > 0){
        pick_eig = 1;
    }
    //
    int i;
    double bignum = 1e+07;
    //
    arma::mat Xi(2*m,2*m);
    arma::mat Delta(2*m,2*m);
    //
    arma::mat Psi = F;
    arma::mat Gamma = -G;
    arma::mat Theta = -H;
    //
    arma::mat Cpinv = C;
    if(l > 0){
        Cpinv = arma::pinv(C);
        Psi = F - J*Cpinv*A;
        Gamma = J*Cpinv*B - G + K*Cpinv*A;
        Theta = K*Cpinv*B - H;
    }
    //
    Xi(arma::span(0,m-1),arma::span(0,m-1)) = Gamma;
    Xi(arma::span(0,m-1),arma::span(m,2*m-1)) = Theta;
    Xi(arma::span(m,2*m-1),arma::span(0,m-1)) = arma::eye(m,m);
    Xi(arma::span(m,2*m-1),arma::span(m,2*m-1)).zeros();
    //
    Delta(arma::span(0,m-1),arma::span(0,m-1)) = Psi;
    Delta(arma::span(0,m-1),arma::span(m,2*m-1)).zeros();
    Delta(arma::span(m,2*m-1),arma::span(0,m-1)).zeros();
    Delta(arma::span(m,2*m-1),arma::span(m,2*m-1)) = arma::eye(m,m);
    /*
     * Next perform a generalized eigen decomp:
     */
    arma::cx_vec EigValue;
    arma::cx_mat EigVec;
    // generalized eigen decomp:
    arma::eig_pair(EigValue, EigVec, Xi, Delta);
    //
    arma::vec EigValueAbs = abs(EigValue);
    arma::vec EigValueReal = real(EigValue);
    /*
     * Deal with infinite values (otherwise sort will return an error):
     */
    arma::uvec infindices = find_nonfinite(EigValueAbs);
    arma::uvec neginfind = arma::find(EigValueReal < -bignum);
    //
    int infin2 = infindices.n_elem;
    int infin3 = neginfind.n_elem;
    //
    if(infin2 > 0){
        arma::cx_vec BigNum(infin2);
        BigNum.ones();
        //
        EigValueAbs.elem(infindices) = arma::ones(infin2)*bignum;
        EigValue.elem(infindices) = BigNum*bignum;
    }
    if(infin3 > 0){
        arma::cx_vec BigNum2(infin3);
        BigNum2.ones();
        //
        EigValue.elem(neginfind) = - BigNum2*bignum;
    }
    /*
     * Now sort the eigenvalues and eigenvectors...
     */
    arma::uvec indices = arma::sort_index(EigValueAbs);
    // ... from smallest to largest in absolute value
    arma::vec EigValueAbsSorted = EigValueAbs.elem(indices);
    arma::cx_vec EigValueSorted = EigValue.elem(indices);
    arma::cx_mat EigVecSorted = EigVec.cols(indices);
    /*
     * If the user prefers to 'choose' which eigenvalues to use...
     */
    if(pick_eig > 0){
        /*
         * by 'egvecind-1' this implies that the elements of 'which_eig' begin at 1 and not 0
         */
        indices = arma::conv_to<arma::uvec>::from(which_eig - 1);
        //
        EigValueSorted = EigValueSorted.elem(indices);
        EigVecSorted = EigVecSorted.cols(indices);
    }
    //
    //
    //
    //arma::cx_mat EigVecSortedRet = EigVecSorted;
    //arma::cx_vec EigValueSortedRet = EigValueSorted;
    eigenvals = EigValueSorted;
    eigenvecs = EigVecSorted;
    //
    arma::cx_vec LambdaVec = EigValueSorted.rows(0,m-1);
    arma::cx_mat Lambda = arma::diagmat(LambdaVec);
    //
    arma::cx_mat Omega = EigVecSorted(arma::span(m,2*m-1),arma::span(0,m-1));
    //
    P = arma::real( Omega*Lambda*arma::pinv(Omega) );
    /*
     * Now calculate Q, R, and S
     */
    Cpinv = C;
    //
    arma::mat LNM = L*N + M;
    arma::colvec LNMStacked = arma::vectorise(LNM);
    //
    arma::mat V = arma::kron(arma::trans(N),F) + arma::kron(arma::eye(k,k),(F*P+G));
    //
    arma::vec QS;
    R.set_size(0,P.n_cols); R.zeros();
    if(l == 0){
        QS = -arma::inv(V)*LNMStacked;
    }else{
        /*
         * First R, ...
         */
        Cpinv = pinv(C);
        R = - Cpinv*(A*P + B);
        /*
         * then Q
         */
        arma::mat V2 = arma::zeros<arma::mat>((k*A.n_rows) + V.n_rows, (k*A.n_rows) + V.n_rows);
        //
        V2(arma::span(0,(k*A.n_rows)-1),arma::span(0,(k*A.n_cols)-1)) = arma::kron(arma::eye(k,k),A);
        V2(arma::span(0,(k*A.n_rows)-1),arma::span(k*A.n_cols,(k*A.n_cols)+(k*C.n_cols)-1)) = arma::kron(arma::eye(k,k),C);
        V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(0,(k*F.n_cols)-1)) = arma::kron(trans(N),F) + arma::kron(arma::eye(k,k),(F*P + J*R + G));
        V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(k*F.n_cols,(k*F.n_cols)+(k*J.n_cols)-1)) = arma::kron(trans(N),J) + arma::kron(arma::eye(k,k),K);
        //
        arma::colvec DStacked = arma::vectorise(D);
        //
        arma::mat DLNM(DStacked.n_rows + LNMStacked.n_rows,1);
        DLNM.rows(0,DStacked.n_rows-1) = DStacked;
        DLNM.rows(DStacked.n_rows,DStacked.n_rows + LNMStacked.n_rows - 1) = LNMStacked;
        //
        QS = - arma::inv(V2)*DLNM;
    }
    arma::vec QVec = QS.rows(0,(m*k)-1);
    Q = arma::reshape(QVec, m, k);
    /*
     * Finally, S...
     */
    S.set_size(0,Q.n_cols); S.zeros();
    if(l > 0){
        arma::vec Sc = QS.rows(m*k,(m+n)*k-1);
        S = arma::reshape(Sc, n, k);
    }
    //
    return 0;
}
//
//
//END
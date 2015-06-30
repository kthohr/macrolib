/* 
 * Chris Sims' gensys program for solving
 * linear rational expectations models using
 * a complex generalized Schur (QZ) decomposition.
 *
 * Ported and adapted to C++ by
 * Keith O'Hara
 * 06/29/15
 *
 * Requires the 'Armadillo' linear algebra library, and
 * additional files to call LAPACK's ZGGES routine.
 *
 */
#include <iostream>
#include "armadillo"
#include "armadillo_qz"
//
using namespace std;
//
int qzswitch(int i, arma::cx_mat &A, arma::cx_mat &B, arma::cx_mat &Q, arma::cx_mat &Z){
    //
    double realsmall = 1e-08;
    //
    std::complex<double> a = A(i,i);         std::complex<double> d = B(i,i);
    std::complex<double> b = A(i,i+1);     std::complex<double> e = B(i,i+1);
    std::complex<double> c = A(i+1,i+1); std::complex<double> f = B(i+1,i+1);
    //
    arma::cx_mat wz(2,1); arma::cx_mat wz1(2,1); arma::cx_mat wz2(2,2); arma::cx_mat wzt;
    arma::cx_mat xy(1,2); arma::cx_mat xy1; arma::cx_mat xy2(2,2); arma::cx_mat xyt;
    //
    std::complex<double> n,m;
    //
    //
    //
    if( std::abs(c) < realsmall && std::abs(f) < realsmall ){
        if( std::abs(a) < realsmall ){
            // l.r. coincident 0's with u.l. of A=0; do nothing
            return 0;
        }
        else{
            // l.r. coincident zeros; put 0 in u.l. of a
            wz.set_size(2,1);
            //
            wz(0,0) = b;
            wz(1,0) = -a;
            //
            wz = wz / arma::as_scalar(arma::sqrt(arma::trans(wz)*wz));
            wz1 = wz;
            wzt = arma::trans(wz);
            //
            wz2(arma::span(),0) = wz1;
            //
            wz2(0,1) = wzt(0,1);
            wz2(1,1) = -wzt(0,0);
            //
            xy2 = arma::eye<arma::cx_mat>(2,2);
        }
    }
    else if( std::abs(a) < realsmall && std::abs(d) < realsmall ){
        if(std::abs(c) < realsmall){
            // u.l. coincident zeros with l.r. of A=0; do nothing
            return 0;
        }
        else{
            // u.l. coincident zeros; put 0 in l.r. of A
            wz2 = arma::eye<arma::cx_mat>(2,2);
            //
            xy.set_size(1,2);
            xy(0,0) = c; xy(0,1) = - b;
            //
            xy = xy / arma::as_scalar(arma::sqrt(xy*arma::trans(xy)));
            xy1 = xy;
            xyt = arma::trans(xy);
            //
            xy2(0,0) = xyt(1,0);
            xy2(0,1) = -xyt(0,0);
            //
            xy2(1,arma::span()) = xy1;
        }
    }
    else{
        // usual case
        wz.set_size(1,2);
        wz(0,0) = c*e - f*b; wz(0,1) = std::conj(c*d - f*a);
        //
        xy.set_size(1,2);
        xy(0,0) = std::conj(b*d - e*a); xy(0,1) = std::conj(c*d - f*a);
        //
        n = arma::as_scalar(arma::sqrt(wz*arma::trans(wz)));
        m = arma::as_scalar(arma::sqrt(xy*arma::trans(xy)));
        //
        if( std::abs(m) < 1e-14){
            // all elements of A and B proportional
            return 0;
        }
        //
        wz = wz / n;
        xy = xy / m;
        //
        wz1 = wz; xy1 = xy;
        wzt = arma::trans(wz); xyt = arma::trans(xy);
        //
        wz2(0,arma::span()) = wz1;
        wz2(1,0) = - wzt(1,0); wz2(1,1) = wzt(0,0);
        //
        xy2(0,arma::span()) = xy1;
        xy2(1,0) = - xyt(1,0);
        xy2(1,1) = xyt(0,0);
        //
    }
    //
    wz = wz2;
    xy = xy2;
    //
    A(arma::span(i,i+1),arma::span()) = xy*A(arma::span(i,i+1),arma::span());
    B(arma::span(i,i+1),arma::span()) = xy*B(arma::span(i,i+1),arma::span());
    A(arma::span(),arma::span(i,i+1)) = A(arma::span(),arma::span(i,i+1))*wz;
    B(arma::span(),arma::span(i,i+1)) = B(arma::span(),arma::span(i,i+1))*wz;
    Z(arma::span(),arma::span(i,i+1)) = Z(arma::span(),arma::span(i,i+1))*wz;
    Q(arma::span(i,i+1),arma::span()) = xy*Q(arma::span(i,i+1),arma::span());
    //
    return 0;
    //
}
//
int qzdiv(double stake, arma::cx_mat &A, arma::cx_mat &B, arma::cx_mat &Q, arma::cx_mat &Z){
    //
    int n = A.n_rows;
    //
    arma::cx_colvec a = A.diag();
    arma::cx_colvec b = B.diag();
    //
    arma::mat roots(n,2);
    roots(arma::span(),0) = arma::abs(a);
    roots(arma::span(),1) = arma::abs(b);
    //
    int i,j,k,m;
    for(i = 0; i < n; ++i){
        if( roots(i,0) < 1.e-13 ){
            roots(i,0) += - (roots(i,0) + roots(i,1));
        }
    }
    //
    roots(arma::span(),1) = roots(arma::span(),1) / roots(arma::span(),0);
    //
    //
    //
    double tmp = 0;
    //
    for(i = n; i >= 1; --i){
        //
        m=0;
        //
        for(j = i; j >= 1; --j){
            if (roots(j-1,1) > stake || roots(j-1,1) < -0.1){
                m=j;
                break;
            }
        }
        if(m==0){
            return 0;
        }
        //
        for(k=m; k <= i-1; ++k){
            qzswitch(k-1,A,B,Q,Z); // notice that we have k - 1 here, not k as in Matlab code
            tmp = roots(k-1,1);
            roots(k-1,1) = roots(k,1);
            roots(k,1) = tmp;
        }
    }
    //
    return 0;
    
}
/*
 *  Main gensys function, based on a discrete-time dynamic system of the form
 *
 *         Gamma0*y(t) = C + Gamma1*y(t-1) + Psi*z(t) + Pi*eta(t)
 *
 *  where:         z(t) is a vector of exogenous shocks,
 *               eta(t) is a vector of 'one-step-ahead' expectational errors
 *
 *  Output is of the form:
 *
 *         y(t) = Cons + G1*y(t-1) + impact*z(t) + ywt*inv(I - fmat*inv(L))*fwt*z(t+1).
 *
 *  In a lot of cases z is iid, and so the final term drops out.
 */
int gensys(arma::mat Gamma0, arma::mat Gamma1, arma::mat C, arma::mat Psi, arma::mat Pi, arma::mat & G1, arma::mat & Cons, arma::mat & impact, arma::cx_mat & fmat, arma::cx_mat & fwt, arma::cx_mat & ywt, arma::cx_mat & gev, arma::vec & eu, arma::mat & loose)
{
    //
    int n = Gamma0.n_rows;
    //
    eu.set_size(2);
    eu.zeros();
    //
    arma::cx_mat xGamma0 = arma::zeros<arma::cx_mat>(n,n);
    arma::cx_mat xGamma1 = arma::zeros<arma::cx_mat>(n,n);
    //
    xGamma0.set_real(Gamma0);
    xGamma1.set_real(Gamma1);
    /*
     *
     * QZ decomposition: Gamma0 = Q S Z**H,  Gamma1 = Q T Z**H
     * Note: this is not equivalent to the output produced by Matlab's 'qz' function.
     *
     */
    arma::cx_mat Q; arma::cx_mat Z; arma::cx_mat S; arma::cx_mat T;
    //
    arma::qz(Q, Z, S, T, xGamma0, xGamma1);
    /* 
     *
     * Next, transpose Q as our lapack output is of the form:
     *                Gamma0 = Q S Z**H,    Gamma1 = Q T Z**H,
     * whereas Matlab output is of the form:
     *                Gamma0 = Q**H S Z**H, Gamma1 = Q**H T Z**H
     *
     */
    Q = arma::trans(Q);
    //
    arma::cx_vec alpha_mat = S.diag();
    arma::cx_vec beta_mat = T.diag();
    //
    int nstable = 0; // this is 'nunstab' in Sims' code
    double zxz = 0;
    double divhat;
    double div = 1.01;
    double small = 1e-06;
    //
    int i,j;
    for(i = 0; i < n; ++i){
        //
        if( std::abs(alpha_mat(i)) > 0 ){
            divhat = std::abs(beta_mat(i))/std::abs(alpha_mat(i));
            //
            if( 1 + small < divhat && divhat <= div){
                div = 0.5*(1 + divhat);
            }
            //
        }
        //
        if( std::abs(beta_mat(i)) > div*std::abs(alpha_mat(i)) ){
            nstable = nstable + 1;
        }
        //
        if( std::abs(alpha_mat(i)) < small && std::abs(beta_mat(i)) < small ){
            zxz = 1;
        }
    }
    //
    if(zxz == 0){
        qzdiv(div,S,T,Q,Z);
    }
    else{ // zxz == 1
        std::cout << "Coincident zeros. Indeterminacy and/or nonexistence." << std::endl;
        eu(0) = -2; eu(1) = -2;
        //
        return 0;
    }
    //
    gev.set_size(n,2);
    gev(arma::span(),0) = S.diag();
    gev(arma::span(),1) = T.diag();
    //
    arma::cx_mat Q1 = Q(arma::span(0,n-nstable-1),arma::span());
    arma::cx_mat Q2 = Q(arma::span(n-nstable,n-1),arma::span());
    arma::cx_mat Z1 = arma::trans(Z(arma::span(),arma::span(0,n-nstable-1)));
    arma::cx_mat Z2 = arma::trans(Z(arma::span(),arma::span(n-nstable,n-1)));
    arma::cx_mat S2 = S(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1));
    arma::cx_mat T2 = T(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1));
    //
    arma::cx_mat etawt = Q2*Pi;
    int neta = Pi.n_cols;
    //
    arma::cx_mat ueta = arma::zeros<arma::cx_mat>(0,0);
    arma::vec deta = arma::zeros<arma::vec>(0);
    arma::cx_mat veta = arma::zeros<arma::cx_mat>(neta,0);
    //
    if(nstable > 0){
        arma::cx_mat tueta; arma::vec tdeta; arma::cx_mat tveta;
        arma::uvec bigev;
        //
        arma::svd(tueta,tdeta,tveta,etawt);
        //
        bigev = arma::find(tdeta > small);
        //
        ueta = tueta.cols(bigev);
        veta = tveta.cols(bigev);
        deta = tdeta(bigev);
        //
        if(bigev.n_elem >= nstable){ // existence
            eu(0) = 1;
        }
    }
    //
    //
    //
    arma::cx_mat etawt1 = arma::zeros<arma::cx_mat>(0,neta);
    arma::cx_mat ueta1 = arma::zeros<arma::cx_mat>(0,0);
    arma::vec deta1 = arma::zeros<arma::vec>(0);
    arma::cx_mat veta1 = arma::zeros<arma::cx_mat>(neta,0);
    //
    if(nstable != n){
        arma::cx_mat tueta1; arma::vec tdeta1; arma::cx_mat tveta1;
        //
        etawt1 = Q1*Pi;
        //
        arma::svd(tueta1,tdeta1,tveta1,etawt1);
        //
        arma::uvec bigev2 = arma::find(tdeta1 > small);
        //
        ueta1 = tueta1.cols(bigev2);
        veta1 = tveta1.cols(bigev2);
        deta1 = tdeta1(bigev2);
    }
    //
    //
    //
    arma::cx_mat loose_temp;
    int uniq = 0;
    int nloose = 0;
    if(veta.n_rows==0){
        uniq = 1;
    }
    else{
        loose_temp = veta1 - veta*arma::trans(veta)*veta1;
        //
        arma::cx_mat ul; arma::vec dl; arma::cx_mat vl;
        arma::svd(ul,dl,vl,loose_temp);
        //
        arma::uvec kfind = arma::find( arma::abs(dl) > small*n);
        //
        if(kfind.n_elem == 0){
            uniq = 1;
        }
        else{
            uniq = 0;
        }
        nloose = kfind.n_elem;
    }
    //
    if(uniq == 1){ // uniqueness
        eu(1) = 1;
    }else{
        std::cout << "Indeterminacy. " << nloose << " loose endog errors." << std::endl;
    }
    /*
     *
     * Now put it all together
     *
     */
    arma::mat detamat = arma::eye(nstable,nstable); // put the singular values in diagonal matrix form
    detamat.diag() = deta;
    //
    arma::mat deta1mat = arma::eye(nstable,nstable);
    deta1mat.diag() = deta1;
    //
    arma::cx_mat tmat(n - nstable,n);
    tmat(arma::span(),arma::span(0,n-nstable-1)) = arma::eye<arma::cx_mat>(n-nstable,n-nstable);
    tmat(arma::span(),arma::span(n-nstable,n-1)) = - arma::trans((ueta*(arma::inv(detamat)*arma::trans(veta))*veta1*deta1mat*arma::trans(ueta1)));
    //
    //
    //
    arma::cx_mat G0(n,n);
    G0.zeros();
    G0(arma::span(0,n-nstable-1),arma::span()) = tmat * S;
    G0(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1)) = arma::eye<arma::cx_mat>(nstable,nstable);
    //
    arma::cx_mat G0I = arma::inv(G0);
    //
    arma::cx_mat G1_temp(n,n);
    G1_temp.zeros();
    G1_temp(arma::span(0,n-nstable-1),arma::span()) = tmat * T;
    //
    G1_temp = G0I*G1_temp;
    //
    //
    //
    int usix = n - nstable + 1;
    //
    arma::cx_mat Cons_temp;
    bool Cstatus = arma::any(arma::vectorise(C));
    if(Cstatus==true){ // if any of the elements of 'C' are non-zero...
        arma::cx_mat C2(n,C.n_cols);
        C2.zeros();
        //
        C2(arma::span(0,n-nstable-1),arma::span()) = tmat * Q * C;
        C2(arma::span(n-nstable,n-1),arma::span()) = arma::inv( S(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) - T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * Q2 * C;
        //
        Cons_temp = G0I*C2;
    }else{ // ... otherwise set to zero. This avoids unnecessary calculations.
        Cons_temp = arma::zeros<arma::cx_mat>(n,C.n_cols);
    }
    //
    //
    //
    arma::cx_mat impact_temp(n,Psi.n_cols);
    impact_temp.zeros();
    impact_temp(arma::span(0,n-nstable-1),arma::span()) = tmat * Q * Psi;
    impact_temp = G0I * impact_temp;
    //
    //
    //
    fmat = arma::inv( T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * S(arma::span(usix-1,n-1),arma::span(usix-1,n-1));
    //
    fwt = - arma::inv( T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * Q2 * Psi;
    //
    ywt = G0I(arma::span(),arma::span(usix-1,n-1));
    //
    loose_temp.set_size(n,neta);
    loose_temp.zeros();
    loose_temp(arma::span(0,n-nstable-1),arma::span()) = etawt1 * (arma::eye<arma::cx_mat>(neta,neta) - veta * arma::trans(veta));
    loose_temp = G0I*loose_temp;
    //
    //
    //
    G1 = arma::real(Z * G1_temp * arma::trans(Z));
    Cons = arma::real(Z * Cons_temp);
    impact = arma::real(Z * impact_temp);
    //
    loose = arma::real(Z * loose_temp);
    ywt = Z * ywt;
    //
    return 0;
}
//
//
//END
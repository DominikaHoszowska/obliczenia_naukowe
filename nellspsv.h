//
// Created by dominika on 03.06.19.
//
#include <iostream>
#include <cblas.h>
#include "lapacke.h"

#ifndef OBLICZENIA_NAUKOWE_NELLSPSV_H
#define OBLICZENIA_NAUKOWE_NELLSPSV_H
void nellspsv(double* A, double* B, int M, int N)
{
    double* C;
    C=  (double *)malloc(sizeof(double)*N*N);
    for(int i=0;i<N*N;i++)
        C[i]=0;
    const double alpha=1;
    const double beta=0;
    CBLAS_LAYOUT Layout;
    Layout = CblasColMajor;
    CBLAS_TRANSPOSE A1,A2;
    A1=CblasTrans;
    A2=CblasNoTrans;
    int lda=M, ldb=M, ldc=N;
    int m=N,n=N,k=M;
    cblas_dgemm(Layout,A1,A2,m,n,k,alpha,A,lda,A,ldb,beta,C,ldc);
    std::cout<<"Macierz A'*A: ";
    for(int i=0;i<N*N;i++)
        std::cout<<C[i]<<" ";
    std::cout<<std::endl;
    double* g;
    g= (double *)malloc(sizeof(double)*N*1);
    std::cout<<"Macierz g: ";

    for(int i=0;i<N;i++)
        g[i]=0;
    m=N;
    k=M;
    n=1;
    cblas_dgemm(Layout,A1,A2,m,n,k,alpha,A,lda,B,ldb,beta,g,ldc);
    for(int i=0;i<N;i++)
        std::cout<<g[i]<<" ";
    std::cout<<std::endl;
    int* X;
    X= (int *)malloc(sizeof(int)*N*1);
    for(int i=0;i<N;i++)
        X[i]=0;
    //BX=g
    int a=LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,C,N,X,g,1);
    std::cout<<"Macierz X: ";

    for(int i=0;i<N;i++)
        std::cout<<g[i]<<" ";
    std::cout<<std::endl;
    free(g);
    free(C);
}
#endif //OBLICZENIA_NAUKOWE_NELLSPSV_H

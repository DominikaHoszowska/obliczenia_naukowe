#include <iostream>
#include <cblas.h>
#include "nellspsv.cpp"
#include "lapacke.h"
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
    for(int i=0;i<N*N;i++)
        std::cout<<C[i]<<" ";
    std::cout<<std::endl;
    double* g;
    g= (double *)malloc(sizeof(double)*N*1);
    for(int i=0;i<N;i++)
        g[i]=0;
    m=N;
    k=M;
    n=1;
    cblas_dgemm(Layout,A1,A2,m,n,k,alpha,A,lda,B,ldb,beta,g,ldc);
    for(int i=0;i<N;i++)
        std::cout<<g[i]<<" ";
    std::cout<<std::endl;
    double* X;
    X= (double *)malloc(sizeof(double)*N*1);
    //BX=g
    n=N;
    int NRHS=1;
    ldb=N;

    free(g);
    free(C);
}
int main() {
    int M=3,N=2;
    double* A;
    A = (double *)malloc(sizeof(double)*M*N);
    double* B, *g;
    B = (double *)malloc(sizeof(double)*M);
    A[0]=1;
    A[1]=1;
    A[2]=1;
    A[3]=1;
    A[4]=2;
    A[5]=3;
    B[0]=1;
    B[1]=1;
    B[2]=1;

    nellspsv(A,B,M,N);
    free(A);
    free(B);
    free(g);
    return(0);
}
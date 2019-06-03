#include <iostream>
#include <cblas.h>
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
int main() {
    std::cout<<"Przykład 1"<<std::endl;
    int M=3,N=2;
    double* A;
    A = (double *)malloc(sizeof(double)*M*N);
    double* B;
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
    std::cout<<"Przykład 2"<<std::endl;
    M=2;
    N=2;
    A = (double *)malloc(sizeof(double)*M*N);
    B = (double *)malloc(sizeof(double)*M);
    A[0]=1;
    A[1]=1;
    A[2]=2;
    A[3]=1;
    B[0]=1;
    B[1]=0;
    nellspsv(A,B,M,N);
    nellspsv(A,B,M,N);
    A[0]=2;
    A[1]=1;
    A[2]=2;
    A[3]=0;
    B[0]=2;
    B[1]=0;
    nellspsv(A,B,M,N);

    return(0);
}
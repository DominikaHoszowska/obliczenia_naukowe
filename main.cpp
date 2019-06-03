#include <iostream>
#include <cblas.h>
#include "lapacke.h"
#include "nellspsv.h"


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
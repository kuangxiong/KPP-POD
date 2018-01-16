/*************************************************************************
	> File Name: KPP_SVD_C.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2017年12月19日 星期二 20时33分26秒
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"mpi.h"
#include"complex.h" 


void KPP_Build_FFTMatrix(int N, int ncup, int localN, complex *mpiUP_hat, 
                    int myrank, int nprocs, complex *tmp_hat1, MPI_Comm comm)
{
    int i, j, k, dest1, dest2, N1, N2, N3;
    complex tmp_v[localN *(N/2-1)*ncup];
    MPI_Status status;

	N1 = localN*N;
	N2 = localN*(N/2+1);
	N3 = localN*(N/2-1);

    //build matrix tmp_hat1  for SVD decomposition
    for(k = 0; k < ncup; k++){
        for(i=0; i< localN; i++){
            for(j=0; j<= N/2; j++){
                tmp_hat1[k*N1+ j + i * N] = mpiUP_hat[k*N2 + j + i * (N/2+1)];   
            }
        }
    }
    //build matrix tmp_v prepare for tmp_hat1 
    for(k=0; k< ncup; k++){
        for(i=0; i< localN; i++){
            for(j=0; j< N/2-1; j++){
                tmp_v[k*N3 + j + i * (N/2-1)] = conj(mpiUP_hat[k*N2 + N/2 - 1 - j + i * (N/2+1)]);   
            }
        }
    }
    // get target process for echo rank 
    dest1 = nprocs - myrank -1;
    dest2 = nprocs - myrank;
    if(myrank==0||myrank==nprocs/2)
        dest2 = myrank; 

//    printf("testrank:%d\t%d\n", myrank, dist2);
    
    for(k=0; k<ncup; k++){
        for(i = 1; i<localN; i++){
          MPI_Sendrecv(&tmp_v[k * N3 + 0+i * (N/2-1)], N/2-1, MPI_COMPLEX, dest1, 991,
                       &tmp_hat1[k *N1 + N/2 + 1 + (localN - i) * N], N/2-1, MPI_COMPLEX,
                       dest1, 991, comm, &status);
        }
    }
    
    for(k=0; k<ncup; k++){
         MPI_Sendrecv(&tmp_v[k*N3 + 0], N/2-1, MPI_COMPLEX, dest2, 992,
                      &tmp_hat1[k*N1 + 0], N/2-1, MPI_COMPLEX,
                     dest2, 992, comm, &status);
    }

    if(myrank==0 || myrank==nprocs/2){
        for(k=0; k< ncup; k++){
            for(i=1; i< N/2; i++){
                tmp_hat1[k*N1 + N/2+i] = conj(tmp_hat1[k *N1 + N/2-i]);
            }
        }
    }
}



void KPP_Build_FFTMatrix1(int N, int ncup, int localN, complex *mpiUP_hat, 
       complex *tmp_hat1)
{
    int i, j, k, N1, N2;
    MPI_Status status;

//	N1 = localN*N;
//	N2 = localN*(N/2+1);

    //build matrix tmp_hat1  for SVD decomposition
    for(k = 0; k < ncup; k++){
        for(i=0; i< localN; i++){
            for(j=0; j<= N/2; j++){
                tmp_hat1[k*(localN*N)+i*N+j]=mpiUP_hat[k*localN*(N/2+1)+i*(N/2+1)+j];   
            }
        }
    }
    //build matrix tmp_v prepare for tmp_hat1 
    for(k=0; k< ncup; k++){
        for(i=0; i< localN; i++){
            for(j=N/2+1; j< N; j++){
                tmp_hat1[k*(localN*N)+j+i*N] = conj(mpiUP_hat[k*localN*(N/2+1)+N-j+i*(N/2+1)]);   
            }
        }
    }
}


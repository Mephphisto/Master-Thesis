//
// Created by jakob on 20.01.20.
//
#include  <time.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  <cuda_runtime.h>
#include  <cusolverDn.h>

#define  BILLION  1000000000L;
#define RND_MX rand()/ (double) RAND_MAX

int main(int argc, char *argv[]) {
    struct timespec start, stop;            //  variables  for  timing
    double accum;                            //  elapsed  time  variable
    cusolverDnHandle_t cusolverH;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat = cudaSuccess;
    const int m = 8192;         //  number  of rows  and  columns  of A
    const int lda = m;                     //  leading  dimension  of A
    cuDoubleComplex *A;                                             // mxm  matrix
    cuDoubleComplex *V;                         // mxm  matrix  of  eigenvectors
    double *W;                            //   m-vector  of  eigenvalues//  prepare  memory  on the  host
    A = (cuDoubleComplex *) malloc(lda * m * sizeof(cuDoubleComplex));
    V = (cuDoubleComplex *) malloc(lda * m * sizeof(cuDoubleComplex));
    W = (double *) malloc(m * sizeof(double));//  define  random A
    for (int i = 0; i < lda * m; i++) {
        double x = rand() / RND_MX;
        double y = rand() / RND_MX;
        A[i] = make_cuDoubleComplex(x, y);
    }//  declare  arrays  on the  device

    cuDoubleComplex *d_A;                       // mxm  matrix A on the  device
    double *d_W;          // m-vector  of  eigenvalues  on the  device-
    int *devInfo;                               // info on the  device
    cuDoubleComplex *d_work;                       //  workspace  on the  device
    int lwork = 0;                                  //  workspace  size
    int info_gpu = 0;           // info  copied  from  device  to host//  create  cusolver  handle
    cusolver_status = cusolverDnCreate(&cusolverH);//  prepare  memory  on the  device
    cudaStat = cudaMalloc((void **) &d_A, sizeof(cuDoubleComplex) * lda * m);
    cudaStat = cudaMalloc((void **) &d_W, sizeof(double) * m);
    cudaStat = cudaMalloc((void **) &devInfo, sizeof(int));
    cudaStat = cudaMemcpy(d_A, A, sizeof(cuDoubleComplex) * lda * m, cudaMemcpyHostToDevice);
    // copy A->d_A//  compute  eigenvalues  and  eigenvectors
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;// use  lower  left  triangle  of the  matrix
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;//  compute  buffer  size  and  prepare  workspace
    cusolver_status = cusolverDnZheevd_bufferSize(cusolverH, jobz, uplo, m, d_A, lda, d_W, &lwork);
    cudaStat = cudaMalloc((void **) &d_work, sizeof(cuDoubleComplex) * lwork);
    clock_gettime(CLOCK_REALTIME, &start);           // start  timer
    //  compute  the  eigenvalues  and  eigenvectors  for a symmetric ,
    // real  mxm  matrix (only  the  lower  left  triangle  af A is used)
    cusolver_status = cusolverDnZheevd(cusolverH, jobz, uplo, m, d_A, lda, d_W, d_work, lwork, devInfo);
    cudaStat = cudaDeviceSynchronize();
    clock_gettime(CLOCK_REALTIME, &stop);             // stop
    accum = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) / (double) BILLION;
    printf("Dsyevd  time: %lf sec.\n", accum);// print  elapsed  time

    cudaStat = cudaMemcpy(V, d_A, sizeof(cuDoubleComplex) * lda * m, cudaMemcpyDeviceToHost);// copy d_A ->V
    cudaStat = cudaMemcpy(W, d_W, sizeof(double) * m, cudaMemcpyDeviceToHost);             // copy d_W ->W
    cudaStat = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost); // copy  devInfo ->info_gpu
    printf("after  syevd: info_gpu = %d\n", info_gpu);
    printf("eigenvalues :\n");           //  print  first  eigenvalues
    for (int i = 0; i < 3; i++) {
        printf("W[%d] = %E\n", i + 1, W[i]);
    }// free  memory
    cudaFree(d_A);
    cudaFree(d_W);
    cudaFree(devInfo);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaDeviceReset();
    return 0;
}
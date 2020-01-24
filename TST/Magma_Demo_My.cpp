#include  <cuda.h>
#include <iostream>
#include <complex>
#include "magma_v2.h"
#include "magma_lapack.h"

int main(int argc, char **argv) {
    magma_init(); //  initialize  Magma
    magma_queue_t queue = NULL;
    magma_int_t dev = 0; // Devce iD
    magma_queue_create(dev, &queue);
    //for (magma_int_t n = 4; n <= 15000; n = sqrt(2)*n ) {
        //std::cout << "n = " << n << std::endl;
        double gpu_time, cpu_time; // variable for GPU Calculaton
        magma_int_t n = 8192, n2 = n * n;
        //magma_int_t n2 = n * n;
        magmaDoubleComplex *Matrix_A_host, *Matrix_R_host;           // Matrix_A_host, r - nxn  matrices  on the  host
        magmaDoubleComplex *Matrix_A_device;            // nxn  matrix  on the  device
        magmaDoubleComplex *h_work;                  //  workspace
        double *rwork;
        magma_int_t lwork, lrwork;           //  h_work  size
        magma_int_t *iwork;         //  workspace
        magma_int_t liwork;        // iwork  size
        double *w1, *w2;         // w1,w2 - vectors  of   eigenvalues
        double error, work[1];
        magma_int_t ione = 1, info;
        double mione = -1.0;
        magma_int_t incr = 1;
        magma_int_t ISEED[4] = {0, 0, 0, 1};                    // seed
        magma_dmalloc_cpu(&w1, n);               // host  memory  for  real
        magma_dmalloc_cpu(&w2, n);                         //  eigenvalues
        magma_zmalloc_cpu(&Matrix_A_host, n2);                  // host  memory  for Matrix_A_host
        magma_zmalloc_cpu(&Matrix_R_host, n2);                  // host  memory  for r
        magma_zmalloc(&Matrix_A_device, n2);                //  device  memory  for Matrix_A_device//  Query  for  workspace  sizes
        magmaDoubleComplex aux_work[1];
        double aux_rwork[1];
        magma_int_t aux_iwork[1];
        magma_zheevd_gpu(
                MagmaVec, // jobz
                MagmaLower, // uplo
                n,
                Matrix_A_device, //dA
                n, //dlla
                w1, //w
                Matrix_R_host,
                n,
                aux_work,
                -1,
                aux_rwork,
                -1,
                aux_iwork,
                -1,
                &info);
        lwork = static_cast<magma_int_t>(aux_work[0].x);
        lrwork = static_cast<magma_int_t>(aux_rwork[0]);
        liwork = aux_iwork[0];
        iwork = (magma_int_t *) malloc(liwork * sizeof(magma_int_t));
        magma_zmalloc_cpu(&h_work, lwork); //  memory  for  workspace
        magma_dmalloc_cpu(&rwork, lrwork);
        //  Randomize  the  matrix Matrix_A_host and  copy  Matrix_A_host -> r
        lapackf77_zlarnv(&ione, ISEED, &n2, Matrix_A_host);
        lapackf77_zlacpy(MagmaFullStr, &n, &n, Matrix_A_host, &n, Matrix_R_host, &n);
        magma_zsetmatrix(n, n, Matrix_A_host, n, Matrix_A_device, n, queue);
        // copy Matrix_A_host -> Matrix_A_device
        // compute  the  eigenvalues   and  eigenvectors  for Matrix_A_host symmetric , real  nxn  matrix
        // Magma  version
        gpu_time = magma_sync_wtime(NULL);
        magma_zheevd_gpu(
                MagmaVec,
                MagmaLower,
                n,
                Matrix_A_device,
                n,
                w1,
                Matrix_R_host,
                n,
                h_work,
                lwork,
                rwork,
                lrwork,
                iwork,
                liwork,
                &info);
        gpu_time = magma_sync_wtime(NULL) - gpu_time;
        printf("ssyevd  gpu  time: %7.5f sec.\n", gpu_time);

        cpu_time = magma_wtime();
        lapackf77_zheevd(
                "V",
                "L",
                &n,
                Matrix_A_host,
                &n,
                w2,
                h_work,
                &lwork,
                rwork,
                &lrwork,
                iwork,
                &liwork,
                &info);
        cpu_time = magma_wtime() - cpu_time;
        printf("ssyevd  cpu  time: %7.5f sec.\n", cpu_time); //  Lapack
//  difference  in   eigenvalues                                // time
        blasf77_daxpy(&n, &mione, w1, &incr, w2, &incr);
        error = lapackf77_dlange("M", &n, &ione, w2, &n, work);
        printf("difference  in  eigenvalues: %e\n", error);
        free(w1);                       // free  host
        free(w2);            // free  host
        free(Matrix_A_host);             // free  host
        free(Matrix_R_host);                // free  host
        free(h_work);       // free  host
        free(rwork);
        magma_free(Matrix_A_device); // free  devic
    //}
    magma_queue_destroy(queue);         //  destroy  queue
    magma_finalize();            //  finalize  Magma
    return EXIT_SUCCESS;
}
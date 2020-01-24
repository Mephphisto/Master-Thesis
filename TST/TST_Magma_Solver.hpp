//
// Created by jakob on 21.01.20.
//

#pragma once

#include <Eigen/src/Eigenvalues/Tridiagonalization.h>
#include <cuda.h>
#include "magma_v2.h"
#include "magma_lapack.h"

namespace magma {
    template<typename _MatrixType>
    class SelfAdjointEigenSolver {
        typedef _MatrixType MatrixType;
        enum {
            Size = MatrixType::RowsAtCompileTime,
            ColsAtCompileTime = MatrixType::ColsAtCompileTime,
            Options = MatrixType::Options,
            MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
        };
        typedef Eigen::Index Index;
        typedef typename MatrixType::Scalar Scalar;
        typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;
        typedef typename Eigen::internal::plain_col_type<MatrixType, RealScalar>::type RealVectorType;
        typedef Eigen::Tridiagonalization<MatrixType> TridiagonalizationType;
        typedef typename TridiagonalizationType::SubDiagonalType SubDiagonalType;
        typedef Eigen::Matrix<Scalar, Size, Size, Eigen::ColMajor, MaxColsAtCompileTime, MaxColsAtCompileTime> EigenvectorsType;
    private:
        magma_queue_t queue = NULL;
        magma_int_t dev = 0;
        magma_int_t n, n2;
        magmaDoubleComplex *d_r;                         // nxn  matrix  on the  device
        magmaDoubleComplex *h_work;
        double *rwork; //  workspace
        magma_int_t lrwork, lwork;                                  //  h_work  size
        magma_int_t *iwork;                                  //  workspace
        magma_int_t liwork;                                 // iwork  size

    protected:
        EigenvectorsType m_eivec;
        RealVectorType m_eivalues;
        typename TridiagonalizationType::SubDiagonalType m_subdiag;
        magma_int_t m_info;
        bool m_isInitialized;
        bool m_eigenvectorsOk;

    public:
        SelfAdjointEigenSolver(Index size) {
            magma_init();                                  //  initialize  Magma
            magma_queue_create(dev, &queue);
            n = static_cast<magma_int_t>(size);
            n2 = n * n;
            magma_zmalloc(&d_r, n2);                //  device  memory  for d_r
            //  Query  for  workspace  sizes
            magmaDoubleComplex aux_work[1];
            double aux_rwork[1];
            magma_int_t aux_iwork[1];
            double *aux_w;
            magmaDoubleComplex *aux_r;
            magma_zheevd_gpu(
                    MagmaVec, MagmaLower,
                    n, d_r, n, aux_w, aux_r, n,
                    aux_work, -1,
                    aux_rwork, -1,
                    aux_iwork, -1,
                    &m_info);
            lwork = (magma_int_t) aux_work[0].x;
            lrwork = static_cast<magma_int_t>(aux_rwork[0]);
            liwork = aux_iwork[0];
            iwork = (magma_int_t *) malloc(liwork * sizeof(magma_int_t));
            magma_zmalloc_cpu(&h_work, lwork);
            //  memory  for  workspace//  Randomize  the  matrix a and  copy a -> r
            magma_dmalloc_cpu(&rwork, lrwork);
        }

        void compute(MatrixType M) {

            eigen_assert(matrix.cols() == matrix.rows());
            eigen_assert((options & ~(EigVecMask | GenEigMask)) == 0
                         && (options & EigVecMask) != EigVecMask
                         && "invalid option parameter");

            bool computeEigenvectors = true;
            magma_int_t n = Eigen::internal::convert_index<magma_int_t>(M.cols());
            m_eivalues.resize(n, 1);
            m_subdiag.resize(n - 1);
            m_eivec = M;

            if (n == 1) {
                m_eivalues.coeffRef(0, 0) = Eigen::numext::real(m_eivec.coeff(0, 0));
                if (computeEigenvectors) m_eivec.setOnes(n, n);
                m_info = Eigen::Success;
                m_isInitialized = true;
                m_eigenvectorsOk = computeEigenvectors;
                return;
            }

            lapackf77_zlacpy(MagmaFullStr, &n, &n, (magmaDoubleComplex *) m_eivec.data(), &n,
                             (magmaDoubleComplex *) m_eivec.data(), &n);
            magma_zsetmatrix(n, n, (magmaDoubleComplex *) m_eivec.data(), n, d_r, n, queue);

            //  compute  the  eigenvalues   and  eigenvectors  for a symmetric ,// real  nxn  matrix; Magma  version
            magma_zheevd_gpu(
                    MagmaVec, MagmaLower,
                    n, d_r, n,
                    m_eivalues.data(),
                    (magmaDoubleComplex *) m_eivec.data(), n,
                    h_work, lwork,
                    rwork, lrwork,
                    iwork, liwork,
                    &m_info);
            m_isInitialized = true;
            m_eigenvectorsOk = computeEigenvectors;
            return;
        }

        const EigenvectorsType eigenvectors() const {
            eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
            eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
            return m_eivec;
        }

        const RealVectorType &eigenvalues() const {
            eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
            return m_eivalues;
        }

        ~SelfAdjointEigenSolver() {
            free(h_work);                                  // free  host  memory
            magma_free(d_r);                            // free  device  memory
            magma_queue_destroy(queue);                     //  destroy  queue
            magma_finalize();
        }
    };
}
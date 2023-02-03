!===============================================================================
! Copyright 2021-2022 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!      Intel(R) oneAPI Math Kernel Library (oneMKL)
!      FORTRAN OpenMP offload examples for SGETRS_BATCH_STRIDED
!*******************************************************************************

include "mkl_omp_offload.f90"

program solve_batched_linear_systems

! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
    use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
    use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif

    integer,   parameter :: n = 50, batch_size = 20, nrhs = 5
    integer              :: lda, stride_a, stride_ipiv
    integer              :: ldb, stride_b
    real,    allocatable :: a(:,:), b(:,:)
    integer, allocatable :: ipiv(:,:), info_rf(:), info_rs(:)

    integer i, j, imat, allocstat
    character (len = 132) :: allocmsg

    lda         = n
    stride_a    = n * lda
    stride_ipiv = n
    ldb         = n
    stride_b    = n * nrhs

    ! Allocate required memory
    allocate (a(stride_a, batch_size), b(n, batch_size*nrhs), ipiv(stride_ipiv, batch_size), &
              info_rf(batch_size), info_rs(batch_size), stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

    ! Initialize the matrices with a random number in the interval (-0.5, 0.5)
    call random_seed()
    call random_number(a)
    a = 0.5 - a
    ! Make diagonal entries larger to ensure well-conditioned matrices
    do i = 1, n
        a(i + (i-1)*lda, :) = a(i + (i-1)*lda, :) + 5.0
    enddo
    ! Initialize the RHS with a random number in the interval (-2.5, 2.5)
    call random_number(b)
    b = 2.5 - (5.0 * b)

    ! Compute the LU factorizations and solve the linear systems using OpenMP offload.
    ! On entry, "a" contains the input matrices. On exit, it contains the factored matrices.
    !$omp target data map(to:a, ipiv) map(tofrom: b) map(from:info_rf, info_rs)
        !$omp dispatch
        call sgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info_rf)

        !$omp dispatch
        call sgetrs_batch_strided('N', n, nrhs, a, lda, stride_a, ipiv, stride_ipiv, &
                                   b, ldb, stride_b, batch_size, info_rs)
    !$omp end target data

    if (any(info_rf .ne. 0)) then
        print *, 'Error: getrf_batch_strided returned with errors.'
    elseif (any(info_rs .ne. 0)) then
        print *, 'Error: getrs_batch_strided returned with errors.'
    else
        print '("Computation executed successfully")'
    endif

    ! Clean up
    deallocate (a, b, ipiv, info_rf, info_rs)
end program solve_batched_linear_systems

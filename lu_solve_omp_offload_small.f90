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
!
!  Content:
!    Intel(R) oneAPI Math Kernel Library (oneMKL)
!    FORTRAN OpenMP offload examples for solving batched linear systems.
!
! Compile for CPU:
!   ifx -i8 -qmkl -free lu_solve_omp_offload_small.f90 -o lu_solve_offload_small
!
! Compile for GPU:
!   ifx -i8 -DMKL_ILP64 -qopenmp -fopenmp-targets=spir64 -fsycl -free \
!       lu_solve_omp_offload_small.f90 -o lu_solve_ex3_omp_offload_small \
!       -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 \
!       -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
!
!*******************************************************************************

include "mkl_omp_offload.f90"

program solve_batched_linear_systems

! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
    use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
    use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif

    integer,   parameter :: n = 3, batch_size = 2, nrhs = 1
    integer              :: lda, stride_a, stride_ipiv
    integer              :: ldb, stride_b
    real,    allocatable :: a(:,:), b(:,:)
    integer, allocatable :: ipiv(:,:), info(:)

    integer i, j, imat, allocstat
    character (len = 132) :: allocmsg

    lda         = n
    stride_a    = n * lda
    stride_ipiv = n
    ldb         = n
    stride_b    = n * nrhs

    ! Allocate required memory
    allocate (a(stride_a, batch_size), b(n, batch_size*nrhs), ipiv(stride_ipiv, batch_size), &
              info(batch_size), stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

    ! Initialize the matrices. Remember that Fortran is a column-major language.
    a(:,1) = (/2.0, 1.0, -6.0, 4.0, -4.0, -9.0, -4.0,  3.0, 10.0/)
    a(:,2) = (/0.0, 2.0,  6.0, 0.0,  4.0,  5.0,  3.0, -1.0,  5.0/)
    b(:,1) = (/12.0, -21.0, -24.0/)        ! x = (-4,  2, -3)
    b(:,2) = (/ 6.0, -10.0,  -7.0/)        ! x = (-2, -1,  2)

    ! Compute the LU factorizations and solve the linear systems using OpenMP offload.
    ! On entry, "a" contains the input matrices. On exit, it contains the factored matrices.
    !$omp target data map(tofrom:a) map(from:ipiv) map(from:info)
        !$omp dispatch
        call sgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info)
    !$omp end target data
    print '("Finished call to sgetrf_batch_strided")'

    if (any(info .ne. 0)) then
        print *, 'Error: getrf_batch_strided returned with errors.'
    else
        ! Solving the linear systems. On exit, the solutions are stored in b.
        !$omp target data map(to:a) map(to:ipiv) map(tofrom: b) map(from:info)
            !$omp dispatch
            call sgetrs_batch_strided('N', n, nrhs, a, lda, stride_a, ipiv, stride_ipiv, &
                                                    b, ldb, stride_b, batch_size, info)
        !$omp end target data
        print '("Finished call to sgetrs_batch_strided")'
        if (any(info .ne. 0)) then
            print *, 'Error: getrs_batch_strided returned with errors.'
        else
            print '("Computation executed successfully")'
            print *, b(:,1)
            print *, b(:,2)
        endif
    endif

    ! Clean up
    deallocate (a, b, ipiv, info)
end program solve_batched_linear_systems

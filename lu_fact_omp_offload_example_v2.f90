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
!      FORTRAN OpenMP offload examples for DGETRI_OOP_BATCH_STRIDED
!*******************************************************************************

include "mkl_omp_offload.f90"

program sgetri_oop_batch_strided_example

#if defined(MKL_ILP64)
    use onemkl_lapack_omp_offload_ilp64
#else
    use onemkl_lapack_omp_offload_lp64
#endif

    integer,   parameter :: n = 10, batch_size = 4
    integer              :: lda, ldainv, stride_a, stride_ainv, stride_ipiv
    real,    allocatable :: a(:,:), ainv(:,:)
    integer, allocatable :: ipiv(:,:), info_rf(:), info_ri(:)

    integer i, j, imat, allocstat
    character (len = 132) :: allocmsg
    real temp

    lda         = n
    ldainv      = n
    stride_a    = n * lda
    stride_ainv = n * ldainv
    stride_ipiv = n

    print '("========================================================================================")'
    print '("Compute the matrix inverse of the matrices in A. Store the results in Ainv.")'
    print '("The paramters for the batch are as follows:")'
    print '("========================================================================================")'
    print '("   Batch size:                          " (I8) )', batch_size
    print '("   Matrix order:                        " (I8) )', n
    print '("   Leading dimension for A matrices:    " (I8) )', lda
    print '("   Leading dimension for Ainv matrices: " (I8) )', ldainv
    print '("   Stride for A matrices:               " (I8) )', stride_a
    print '("   Stride for Ainv matrices:            " (I8) )', stride_ainv
    print '("   Stride for pivot arrays:             " (I8) )', stride_ipiv
    print '("========================================================================================"/)'

    ! Allocate required memory
    allocate(a(stride_a, batch_size), ainv(stride_ainv, batch_size),                  &
             ipiv(stride_ipiv, batch_size), info_rf(batch_size), info_ri(batch_size), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

    ! Initialize the matrices with a random number in the interval (-0.5, 0.5)
    call random_seed()
    call random_number(a)
    a = 0.5 - a
    ! Make diagonal entries larger to ensure well-conditioned matrices
    do i = 1, n
        a(i + (i-1)*lda, :) = a(i + (i-1)*lda, :) + 5.0
    enddo

    print '("====================================================================")'
    print '("                         Input matrices                             ")'
    print '("====================================================================")'
    do imat = 1, batch_size
        do j = 1, n
            do i = 1, n
                write(6, advance='no', fmt="(f8.4)") a(i + (j - 1) * lda, imat)
            enddo
            print *
        enddo
        print *
    enddo

    ! Offload LU factorization and matrix inverse using OpenMP
    print '(//"====================================================================")'
    print '(" Compute LU factorizations and inverses of the LU factored matrices ")'
    print '("====================================================================")'

    !$omp target data map(to:a) map(from:ainv) map(from:info_rf) map(from:info_ri)

        ! On entry, "a" contains the input matrix. On exit, it contains the LU factorization.
        !$omp dispatch
        call sgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info_rf)

        ! Compute the inverse and store it in ainv.
        !$omp dispatch
        call sgetri_oop_batch_strided(n, a, lda, stride_a, ipiv, stride_ipiv, ainv, ldainv, stride_ainv, batch_size, info_ri)

    !$omp end target data

    if (any(info_rf .ne. 0)) then
        print *, 'Error: getrf_batch_strided returned with errors.'
    elseif (any(info_ri .ne. 0)) then
        print *, 'Error: getri_oop_batch_strided returned with errors.'
    else
        print '("Computation executed successfully")'
    endif

    print '(//"====================================================================")'
    print '("                         Inverted matrices                          ")'
    print '("====================================================================")'
    do imat = 1, batch_size
        do j = 1, n
            do i = 1, n
                write(6, advance='no', fmt="(f8.4)") ainv(i + (j - 1) * lda, imat)
            enddo
            print *
        enddo
        print *
    enddo

    ! Clean up
    deallocate(a, ainv, ipiv, info_rf, info_ri)
end program sgetri_oop_batch_strided_example

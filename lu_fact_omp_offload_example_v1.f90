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
!      FORTRAN OpenMP offload examples for SGETRI_OOP_BATCH_STRIDED
!*******************************************************************************

include "mkl_omp_offload.f90"

program sgetri_oop_batch_strided_example

! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
    use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
    use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif

    integer,   parameter :: n = 10, batch_size = 4
    integer              :: lda, ldainv, stride_a, stride_ainv, stride_ipiv
    real,    allocatable :: a(:,:), ainv(:,:)
    integer, allocatable :: ipiv(:,:), info(:)

    integer i, j, imat
    real r_num

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
    allocate (a(stride_a, batch_size), ainv(stride_ainv, batch_size), &
              ipiv(stride_ipiv, batch_size), info(batch_size))

    if ( .not.allocated(a)    .OR. &
         .not.allocated(ainv) .OR. &
         .not.allocated(ipiv) .OR. &
         .not.allocated(info) ) then
         print '("[FAILED] Failed allocation")'
         stop 1
    endif

    ! Initialize the matrices with a random number in the interval (-0.5, 0.5)
    call random_seed()
    call random_number(a)
    a = 0.5 - a
    ! Make diagonal entries larger to ensure well-conditioned matrices
    do imat = 1, batch_size
        do i = 1, n
            a(i + (i-1)*lda, imat) = a(i + (i-1)*lda, imat) + 5.0
        enddo
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

    print '(//"====================================================================")'
    print '(" Compute LU factorizations and inverses of the LU factored matrices ")'
    print '("====================================================================")'

    ! Compute the LU factorization using OpenMP offload
    ! On entry, "a" contains the input matrix. On exit, it contains the LU factorization.
    !$omp target data map(tofrom:a) map(from:ipiv) map(from:info)
        !$omp dispatch
        call sgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info)
    !$omp end target data
    print '("Finished call to sgetrf_batch_strided")'

    print '(//"====================================================================")'
    print '("                         Factored matrices                          ")'
    print '("====================================================================")'
    do imat = 1, batch_size
       do j = 1, n
          do i = 1, n
             write(6, advance='NO', fmt="(f8.4)") a(i + (j - 1) * lda, imat)
          enddo
          print *
       enddo
       print *
    enddo

    if (any(info .ne. 0)) then
        print *, 'Error: getrf_batch_strided returned with errors.'
    else
        ! Compute the matrix inverse using OpenMP offload. On exit, the inverse is stored in Ainv.
        !$omp target data map(to:a) map(to:ipiv) map(from:ainv) map(from:info)
            !$omp dispatch
            call sgetri_oop_batch_strided(n, a, lda, stride_a, ipiv, stride_ipiv, ainv, ldainv, &
                                          stride_ainv, batch_size, info)
        !$omp end target data
        print '("Finished call to sgetri_oop_batch_strided")'
    endif

    if (any(info .ne. 0)) then
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
    deallocate (a, ainv, ipiv, info)
end program sgetri_oop_batch_strided_example

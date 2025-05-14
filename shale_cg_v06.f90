! SHALE - spherical harmonic analysis with least squares
! conjugate gradient version
!
! V0.6 - parallelised with OpenMP and MPI
!
! Tobias Wittwer, 2006
program shale
    implicit none
    ! settings
    integer*4, parameter :: nmax = 300              ! maximum degree/order
    integer*4, parameter :: nobs = 259200                 ! number of observations
    character*80, parameter :: obsfile = 'pot_glob_0_300.txt'   ! observation file
    character*80, parameter :: outfile = 'cnm_300.txt'   ! output file
    integer*4, parameter :: u = nmax**2+2*nmax+1     ! number of unknowns
    integer*4, parameter :: size_nbd = 2*(nmax*(nmax+1)*(2*nmax+1)/6)+(nmax+1)**2

    ! constants
    real*8, parameter :: pi = 3.141592653589793
    real*8, parameter :: rho = 180.0d0/pi
    real*8, parameter :: GM = 3986005.0d8
    real*8, parameter :: rearth = 6371000.0d0
    
    ! arrays
    real*8, allocatable :: a_row(:)                    ! design matrix line
    real*8, allocatable :: N(:)                   ! normal matrix (preconditioner)
    real*8, allocatable :: y(:)                      ! observation vector
    real*8, allocatable :: b(:)                      ! rhs vector
    real*8, allocatable :: bold(:)                      ! rhs vector
    real*8, allocatable :: x(:)                      ! unknown parameters
    real*8, allocatable :: p(:)                      ! preconditioning vector
    real*8, allocatable :: a(:)
    real*8, allocatable :: pnm(:,:)                    ! legendre functions
    real*8, allocatable :: long(:),lat(:),rad(:),obs(:)   ! positions
    real*8, allocatable :: bdiff(:),xdiff(:),xinc(:)          ! differences
    real*8, allocatable :: bt(:),xt(:),e(:),xold(:),pt(:),ptold(:)                ! more vectors
    real*8, allocatable :: b2(:,:)

    ! other
    integer*4 :: i, it, j, k
    real*8 :: value, theta, latold, factor
    real*8 :: alpha, beta, maxx,br, norm0, gamma, help
    real*8, parameter :: eps1 = 1e-6
    real*8, parameter :: eps2 = 1e-14
    real*8 :: omp_get_wtime, timetot, timeprecon, timeit
    integer*4 :: omp_get_num_threads,omp_get_thread_num,numthreads,threadnum

    ! MPI
    integer*4 :: ierr, nprocs, myrank

    interface
        real*8 function DDOT(n,dx,incx,dy,incy)
            integer*4 :: n,incx,incy
            real*8 :: dx(*),dy(*)
        end function ddot
        real*8 FUNCTION DNRM2 ( N, X, INCX )
            INTEGER*4 :: INCX, N
            real*8 :: X( * )
        end function dnrm2
    end interface

    include 'mpif.h'

    timetot = omp_get_wtime()

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

    if(myrank.eq.0) then
        write(*,'(A)') '* SHALE CG V0.6 *'
        write(*,'(A,I4)') 'spherical harmonic analysis up to degree ', nmax
        write(*,'(A,I6,A)') 'using ', nobs, ' observations'
        write(*,*) ' '
    end if

    ! allocate arrays
    allocate(N(size_nbd),y(nobs),b(u))
    allocate(long(nobs),lat(nobs),rad(nobs),obs(nobs))
    allocate(x(u),p(u),a(u),bold(u),xold(u))
    allocate(bdiff(u),xdiff(u),bt(u),xt(u),e(u))
    allocate(ptold(u),pt(u),xinc(u))

    ! data input
    open(1,file=obsfile)
    do i=1,nobs
        read(1,*) long(i), lat(i), rad(i), y(i)
        long(i) = long(i)/rho
        lat(i) = lat(i)/rho
    end do
    close(1)

    ! build preconditioner
    if(myrank.eq.0) write(*,'(A)') 'building preconditioner'
    timeprecon = omp_get_wtime()
    call build_precon(long,lat,rad,N,nmax,size_nbd,nobs,u,rearth)
    if(myrank.eq.0) then
        write(*,'(A,F10.3)') 'time spent [s]:', omp_get_wtime()-timeprecon
        write(*,*) ' '
    end if

    !$OMP PARALLEL
    numthreads = omp_get_num_threads()
    !$OMP END PARALLEL
    allocate(b2(u,numthreads))
    
    ! initialisation
    b(:) = 0.0d0
    b2(:,:)=0.0d0
    x(:) = 0.0d0
    xold=x
    xt=x
    !$OMP PARALLEL PRIVATE(a_row,latold,pnm,threadnum)
    latold = 1.6d0
    allocate(pnm(0:nmax,0:nmax))
    allocate(a_row(u))
    threadnum = omp_get_thread_num()+1
    !$OMP DO
    do i=myrank+1,nobs,nprocs
        call build_a_line(a_row,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
        latold = lat(i)
        call daxpy(u,y(i),a_row,1,b2(1,threadnum),1)
    end do
    !$OMP END DO
    deallocate(pnm,a_row)
    !$OMP END PARALLEL
    do j=2,numthreads
        b2(:,1) = b2(:,1)+b2(:,j)
    end do
    call mpi_allreduce(b2,b,u,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    bt = b
    call mult_nbd(p,N,b,nmax,u,size_nbd)
    norm0 = dnrm2(u, bt, 1)
    pt = p

    timeit = omp_get_wtime()
    do it=1,u
        if(myrank.eq.0) write(*,'(A,I4)') "iteration", it
        a(:) = 0.0d0
        b2(:,:) = 0.0d0
        !$OMP PARALLEL PRIVATE(a_row,latold,pnm,help,threadnum)
        latold = 1.6d0
        allocate(pnm(0:nmax,0:nmax))
        allocate(a_row(u))
        threadnum = omp_get_thread_num()+1
        !$OMP DO
        do i=myrank+1,nobs,nprocs
            call build_a_line(a_row,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
            latold = lat(i)
            help = ddot(u,a_row,1,p,1)
            !call daxpy(u,help,a_row,1,a,1)
            call daxpy(u,help,a_row,1,b2(1,threadnum),1)
        end do
        !$OMP END DO
        deallocate(pnm,a_row)
        !$OMP END PARALLEL
        do j=2,numthreads
            b2(:,1) = b2(:,1)+b2(:,j)
        end do
        call mpi_allreduce(b2,a,u,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        alpha = ddot(u, b, 1, p, 1)/ddot(u, a, 1, p, 1)
        call daxpy(u, alpha, p, 1, x, 1)
        bold = b
        call daxpy(u, -alpha, a, 1, b, 1)
        bdiff = bt-b
        xdiff = xt-x
        call mult_nbd(e,N,bdiff,nmax,u,size_nbd)
        gamma = -ddot(u,b,1,e,1)/ddot(u,bdiff,1,e,1)
        xold = xt
        do i=1,u
            xt(i)=x(i)+gamma*xdiff(i)
            bt(i)=b(i)+gamma*bdiff(i)
        end do
        xinc = xt-xold

        br = dnrm2(u, bt, 1)/norm0
        maxx = maxval(dabs(xinc))
        if((br.lt.eps1).and.(maxx.lt.eps2)) then
            if (myrank.eq.0) then
                write(*,*) ' '
                write(*,'(A,I4)') "solution reached in iteration", it
            end if
            exit
        end if

        ptold = pt
        call mult_nbd(pt,N,b,nmax,u,size_nbd)
        beta = ddot(u,b,1,pt,1)/ddot(u, bold, 1, ptold, 1);
        do i=1,u
            p(i) = beta*p(i)+pt(i)
        end do
    end do

    if(myrank.eq.0) then
        write(*,'(A,F10.3)') 'time in solving [s]:', omp_get_wtime()-timeit
    
        ! write output
        open(1,file=outfile)
        factor = rearth/GM
        k=1
        do j=0,nmax
            write(1,'(2I4,2E20.12)') j,0,factor*xt(k),0.0d0
            k = k+1
        end do
        do i=1,nmax
            do j=i,nmax
                write(1,'(2I4,2E20.12)') j,i,factor*xt(k),factor*xt(k+nmax+1-i)
                k = k+1
            end do
            k = k+nmax+1-i
        end do
        close(1)
        write(*,*) ' '
        write(*,'(A,F10.3)') 'total runtime [s]:', omp_get_wtime()-timetot
    end if

    deallocate(N,y,b)
    deallocate(b2)
    deallocate(a,x,p,bold)
    deallocate(xold)
    deallocate(bdiff,xdiff,bt,xt,e)
    deallocate(ptold,pt,xinc)
    deallocate(long,lat,rad,obs)

    call mpi_finalize(ierr)
end program shale

subroutine build_a_line(a,nmax,u,long,lat,rad,latold,p,rearth)
    ! build one line of design matrix
    ! input:
    !   a               line of design matrix
    !   nmax            maximum degree/order
    !   u               number of unknowns
    !   long            longitude of point []
    !   lat             latitude of point []
    !   rad             radius of point [m]
    !   latold          latitude of previous point []
    !   p               array of legendre functions
    !   rearth          earth radius
    implicit none
    real*8, parameter :: pi = 3.141592653589793
    integer*4 :: nmax, u, m, n, j,i,blocksize
    real*8 :: a(u), long, lat, rad, latold
    real*8 :: p(0:nmax,0:nmax)
    real*8 :: cosml, sinml, rearth
    real*8, allocatable :: rrm(:)

    allocate(rrm(0:nmax))

    if(lat.ne.latold) call pnm(nmax,lat,p)
    do n=0,nmax
        rrm(n) = (rearth/rad)**(n+1)
    end do
    j = 1
    do m=0,nmax
        cosml = dcos(dble(m)*long)
        sinml = dsin(dble(m)*long)
        ! cosine terms
        do n=m,nmax
            a(j) = rrm(n)*cosml*p(n,m)
            j = j+1
        end do
        if(m.gt.0) then
            ! sine terms
            do n=m,nmax
                a(j) = rrm(n)*sinml*p(n,m)
                j = j+1
            end do
        end if
    end do
    deallocate(rrm)
end subroutine build_a_line

subroutine build_precon(long,lat,rad,N,nmax,size_nbd,nobs,u,rearth)
    ! build preconditioner
    ! input:
    !   long            longitude of point []
    !   lat             latitude of point []
    !   rad             radius of point [m]
    !   N               preconditioner
    !   nmax            maximum degree/order
    !   size_nbd        size of preconditioner
    !   nobs            number of observations
    !   u               number of unknowns
    !   pnm             array of legendre functions
    !   rearth          earth radius
    implicit none
    integer*4 :: nmax, size_nbd, nobs,u,idx,offset,i,j,k,m,info,lwork,start
    real*8 :: long(nobs),lat(nobs),rad(nobs)
    real*8 :: N(size_nbd),latold, rearth
    real*8, allocatable :: a_block(:,:), pnm(:,:), N2(:,:)
    integer*4 :: omp_get_num_threads,omp_get_thread_num,numthreads,threadnum
    integer*4, parameter :: blocksize = 256
    integer*4 :: numblocks
    integer*4 :: myrank, nprocs, ierr
    include 'mpif.h'

    N(:) = 0.0d0

    numblocks = nobs/blocksize+1
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

    !$OMP PARALLEL
    numthreads = omp_get_num_threads()
    !$OMP END PARALLEL
    allocate(N2(size_nbd,numthreads))
    N2(:,:) = 0.0d0
    
    ! build of preconditioner
    !$OMP PARALLEL PRIVATE(latold,a_block,pnm,offset,idx,threadnum,start,i)
    allocate(pnm(0:nmax,0:nmax))
    allocate(a_block(u,blocksize))
    latold=1.6d0
    threadnum = omp_get_thread_num()+1
    !$OMP DO
    do j=myrank+1,numblocks,nprocs
        do i=1,blocksize
            idx = (j-1)*blocksize+i
            if(idx.gt.nobs) then
                a_block(:,i) = 0.0d0
            else
                call build_a_line(a_block(1,i),nmax,u,long(idx),lat(idx),rad(idx),latold,pnm,rearth)
                latold = lat(idx)
            end if
        end do

        call dsyrk('L','N',nmax+1,blocksize,1.0d0,a_block,u,1.0d0,N2(1,threadnum),nmax+1)
        start = 1+(nmax+1)**2
        offset=nmax+2
        ! order 1 to nmax
        do i=1,nmax
            ! cosine terms
            call dsyrk('L','N',nmax-i+1,blocksize,1.0d0,a_block(offset,1),u,1.0d0,N2(start,threadnum),nmax-i+1)
            start = start+(nmax-i+1)**2
            offset=offset+nmax-i+1
            ! sine terms
            call dsyrk('L','N',nmax-i+1,blocksize,1.0d0,a_block(offset,1),u,1.0d0,N2(start,threadnum),nmax-i+1)
            start = start+(nmax-i+1)**2
            offset=offset+nmax-i+1
        end do
    end do
    !$OMP END DO
    deallocate(a_block)
    deallocate(pnm)
    !$OMP END PARALLEL

    do j=2,numthreads
        N2(:,1) = N2(:,1)+N2(:,j)
    end do
    call mpi_allreduce(N2,N,size_nbd,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    ! blockwise inversion of preconditioner
    ! order 0
    call dpotrf('L',nmax+1,N,nmax+1,info)
    call dpotri('L',nmax+1,N,nmax+1,info)
    call sym2full(nmax+1,N)

    start = 1+(nmax+1)**2
    ! order 1 to nmax
    do i=1,nmax
        ! cosine terms
        call dpotrf('L',nmax-i+1,N(start),nmax-i+1,info)
        call dpotri('L',nmax-i+1,N(start),nmax-i+1,info)
        call sym2full(nmax-i+1,N(start))
        start = start+(nmax-i+1)**2
        ! sine terms
        call dpotrf('L',nmax-i+1,N(start),nmax-i+1,info)
        call dpotri('L',nmax-i+1,N(start),nmax-i+1,info)
        call sym2full(nmax-i+1,N(start))
        start = start+(nmax-i+1)**2
    end do
end subroutine build_precon

subroutine mult_nbd(p,N,b,nmax,u,size_nbd)
    ! multiply vector with preconditioner
    ! p         result vector
    ! N         preconditioner
    ! b         vector
    ! nmax      maximum degree/order
    ! u         number of unknown parameters
    ! size_nbd  number of elements in preconditioner
    implicit none
    integer*4 :: nmax,u,size_nbd
    integer*4 :: i,j,start,offset
    real*8 :: p(u),b(u),N(size_nbd)
    interface
        real*8 function DDOT(n,dx,incx,dy,incy)
            integer*4 :: n,incx,incy
            real*8 :: dx(*),dy(*)
        end function ddot
    end interface

    ! order 0
    do i=1,nmax+1
        p(i) = ddot(nmax+1,N((i-1)*(nmax+1)+1),1,b,1)
    end do

    ! order 1 to nmax
    offset = nmax+1
    start = (nmax+1)**2
    do i=1,nmax
        ! cosine terms
        do j=1,nmax-i+1
            p(offset+j) = ddot(nmax-i+1,N(start+1+(j-1)*(nmax-i+1)),1,b(offset+1),1)
        end do
        start = start+(nmax-i+1)**2
        offset = offset+nmax-i+1
        ! sine terms
        do j=1,nmax-i+1
            p(offset+j) = ddot(nmax-i+1,N(start+1+(j-1)*(nmax-i+1)),1,b(offset+1),1)
        end do
        offset = offset+nmax-i+1
        start = start+(nmax-i+1)**2
    end do
end subroutine mult_nbd

subroutine sym2full(n,A)
    implicit none
    integer*4 :: n,i,j
    real*8 :: A(n**2)
    do i=1,n-1
        do j=1,i
            A(i*n+j) = A((j-1)*n+i+1)
        end do
    end do
end subroutine sym2full


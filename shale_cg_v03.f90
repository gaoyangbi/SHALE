! SHALE - spherical harmonic analysis with least squares
! conjugate gradient version
!
! V0.2 - with time measurement
!
! Tobias Wittwer, 2006
program shale
    implicit none
    ! settings
    integer*4, parameter :: nmax = 50                ! maximum degree/order
    integer*4, parameter :: nobs = 16200                 ! number of observations
    character*80, parameter :: obsfile = 'pot_glob_0_50.txt'   ! observation file
    character*80, parameter :: outfile = 'cnm_50.txt'   ! output file
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

    ! other
    integer*4 :: i, it, j, k
    real*8 :: value, theta, latold, factor
    real*8 :: alpha, beta, maxx,br, norm0, gamma, help
    real*8, parameter :: eps1 = 1e-6
    real*8, parameter :: eps2 = 1e-14
    real*8 :: omp_get_wtime, timetot, timeprecon, timeit

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

    timetot = omp_get_wtime()
    write(*,'(A)') '* SHALE CG V0.3 *'
    write(*,'(A,I4)') 'spherical harmonic analysis up to degree ', nmax
    write(*,'(A,I6,A)') 'using ', nobs, ' observations'
    write(*,*) ' '

    ! allocate arrays
    allocate(a_row(u),N(size_nbd),y(nobs),b(u))
    allocate(pnm(0:nmax,0:nmax))
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
    write(*,'(A)') 'building preconditioner'
    timeprecon = omp_get_wtime()
    call build_precon(long,lat,rad,N,nmax,size_nbd,nobs,u,pnm,rearth)
    write(*,'(A,F10.3)') 'time spent [s]:', omp_get_wtime()-timeprecon
    write(*,*) ' '

    ! initialisation
    b(:) = 0.0d0
    x(:) = 0.0d0
    xold=x
    xt=x
    latold = 1.6d0
    do i=1,nobs
        call build_a_line(a_row,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
        latold = lat(i)
        call daxpy(u,y(i),a_row,1,b,1)
    end do
    bt = b
    call mult_nbd(p,N,b,nmax,u,size_nbd)
    norm0 = dnrm2(u, bt, 1)
    pt = p

    timeit = omp_get_wtime()  ! Add by ztwang.
    do it=1,u
        write(*,'(A,I4)') "iteration", it
        latold = 1.6d0
        a(:) = 0.0d0
        do i=1,nobs
            call build_a_line(a_row,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
            latold = lat(i)
            help = ddot(u,a_row,1,p,1)
            call daxpy(u,help,a_row,1,a,1)
        end do

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
            write(*,'(A,I4)') "solution reached in iteration", it
            exit
        end if

        ptold = pt
        call mult_nbd(pt,N,b,nmax,u,size_nbd)
        beta = ddot(u,b,1,pt,1)/ddot(u, bold, 1, ptold, 1);
        do i=1,u
            p(i) = beta*p(i)+pt(i)
        end do
    end do
    write(*,'(A,F15.3)') 'time in solving [s]:', omp_get_wtime()-timeit

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

    deallocate(a_row,N,y,b)
    deallocate(pnm)
    deallocate(a,x,p,bold)
    deallocate(xold)
    deallocate(bdiff,xdiff,bt,xt,e)
    deallocate(ptold,pt,xinc)
    deallocate(long,lat,rad,obs)
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
    integer*4 :: nmax, u, m, n, j
    real*8 :: a(u), long, lat, rad, latold
    real*8 :: p(0:nmax,0:nmax)
    real*8 :: cosnl, sinnl, rearth
    real*8, allocatable :: rrm(:)

    allocate(rrm(0:nmax))

    if(lat.ne.latold) call pnm(nmax,lat,p)
    do n=0,nmax
        rrm(n) = (rearth/rad)**(n+1)
    end do
    j = 1
    do m=0,nmax
        cosnl = dcos(dble(m)*long)
        sinnl = dsin(dble(m)*long)
        ! cosine terms
        do n=m,nmax
            A(j) = rrm(n)*cosnl*p(n,m)
            j = j+1
        end do
        if(m.gt.0) then
            ! sine terms
            do n=m,nmax
                A(j) = rrm(n)*sinnl*p(n,m)
                j = j+1
            end do
        end if
    end do
    deallocate(rrm)
end subroutine build_a_line

subroutine build_precon(long,lat,rad,N,nmax,size_nbd,nobs,u,pnm,rearth)
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
    real*8, allocatable :: a_row(:), work(:)
    real*8 :: pnm(0:nmax,0:nmax)
    integer*4, allocatable :: ipiv(:)

    allocate(ipiv(nmax+1))
    allocate(work((nmax+1)*(nmax+1)))
    lwork = (nmax+1)*(nmax+1)

    ! build of preconditioner
    allocate(a_row(u))
    latold=1.6d0
    do i=1,nobs
        call build_a_line(a_row,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
        latold = lat(i)

        idx = 1
        ! elements of order 0
        do j=1,nmax+1
            do k=1,nmax+1
                N(idx) = N(idx) + a_row(j)*a_row(k)
                idx = idx+1
            end do
        end do
        
        offset=nmax+1
        ! elements of order 1 to nmax
        do m=1,nmax
            ! cosine
            do j=1,nmax-m+1
                do k=1,nmax-m+1
                    N(idx) = N(idx) + a_row(j+offset)*a_row(k+offset)
                    idx = idx+1
                end do
            end do
            offset=offset+nmax-m+1
            ! sine
            do j=1,nmax-m+1
                do k=1,nmax-m+1
                    N(idx) = N(idx) + a_row(j+offset)*a_row(k+offset)
                    idx = idx+1
                end do
            end do
            offset=offset+nmax-m+1
        end do
    end do
    deallocate(a_row)

    ! blockwise inversion of preconditioner

    ! order 0
    call dgetrf(nmax+1,nmax+1,N,nmax+1,ipiv,info)
    call dgetri(nmax+1,N,nmax+1,ipiv,work,lwork,info)

    start = 1+(nmax+1)**2
    ! order 1 to nmax
    do i=1,nmax
        ! cosine terms
        call dgetrf(nmax-i+1,nmax-i+1,N(start),nmax-i+1,ipiv,info)
        call dgetri(nmax-i+1,N(start),nmax-i+1,ipiv,work,lwork,info)
        start = start+(nmax-i+1)**2
        ! sine terms
        call dgetrf(nmax-i+1,nmax-i+1,N(start),nmax-i+1,ipiv,info)
        call dgetri(nmax-i+1,N(start),nmax-i+1,ipiv,work,lwork,info)
        start = start+(nmax-i+1)**2
    end do
    write(2,*) N

    deallocate(ipiv)
    deallocate(work)
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

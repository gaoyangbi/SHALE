! SHALE - spherical harmonic analysis with least squares
!
! V0.4 - parallelised with OpenMP
!
! Tobias Wittwer, 2006
program shale
    implicit none
    ! settings
    integer*4, parameter :: nmax = 100                ! maximum degree/order
    integer*4, parameter :: nobs = 259200                  ! number of observations
    character*80, parameter :: obsfile = 'pot_glob_0_300.txt'   ! observation file
    character*80, parameter :: outfile = 'cnm_100_ls.txt'   ! output file
    integer*4, parameter :: u = nmax**2+2*nmax+1     ! number of unknowns

    ! constants
    real*8, parameter :: pi = 3.141592653589793
    real*8, parameter :: rho = 180.0d0/pi
    real*8, parameter :: GM = 3986005.0d8
    real*8, parameter :: rearth = 6371000.0d0
    
    ! arrays
    real*8, allocatable :: A(:,:)                    ! design matrix
    real*8, allocatable :: N(:,:)                    ! normal matrix
    real*8, allocatable :: y(:)                      ! observation vector
    real*8, allocatable :: b(:)                      ! rhs vector
    real*8, allocatable :: pnm(:,:)                  ! legendre functions
    real*8, allocatable :: long(:), lat(:), rad(:)   ! positions

    ! other
    integer*4 :: i, info, j, k
    real*8 :: latold, factor
    real*8 :: omp_get_wtime, timetot, time

    timetot = omp_get_wtime()

    write(*,'(A)') '* SHALE V0.4 *'
    write(*,'(A,I4)') 'spherical harmonic analysis up to degree ', nmax
    write(*,'(A,I6,A)') 'using ', nobs, ' observations'
    write(*,*) ' '

    ! allocate arrays
    allocate(A(nobs,u),N(u,u),y(nobs),b(u))
    allocate(long(nobs),lat(nobs),rad(nobs))

    ! data input
    open(1,file=obsfile)
    do i=1,nobs
        read(1,*) long(i), lat(i), rad(i), y(i)
        long(i) = long(i)/rho
        lat(i) = lat(i)/rho
    end do
    close(1)

    ! build design matrix
    time = omp_get_wtime()
    !$OMP PARALLEL PRIVATE(latold,pnm)
    allocate(pnm(0:nmax,0:nmax))
    latold = 1.6d0

    !$OMP DO
    do i=1,nobs
        call build_a_line(A,i,nobs,nmax,u,long(i),lat(i),rad(i),latold,pnm,rearth)
        latold = lat(i)
    end do
    !$OMP END DO
    deallocate(pnm)
    !$OMP END PARALLEL
    write(*,'(A,F10.3)') 'build A [s]      :', omp_get_wtime()-time

    ! build normal equation matrix and rhs vector
    time = omp_get_wtime()
    call DSYRK('L','T', u, nobs, 1.0d0, A, nobs, 0.0d0, N, u)
    write(*,'(A,F10.3)') 'build N [s]      :', omp_get_wtime()-time
    time = omp_get_wtime()
    call DGEMV('T', nobs, u, 1.0d0, A, nobs, y, 1, 0.0d0, b, 1)
    write(*,'(A,F10.3)') 'build b [s]      :', omp_get_wtime()-time

    ! solve
    time = omp_get_wtime()
    call DPOSV('L', u, 1, N, u, b, u, info)
    write(*,'(A,F10.3)') 'solving [s]      :', omp_get_wtime()-time
    write(*,'(A,I10)') 'info             :', info

    ! write output
    open(1,file=outfile)
    factor = rearth/GM
    k=1
    do j=0,nmax
        write(1,'(2I4,2E20.12)') j,0,factor*b(k),0.0d0
        k = k+1
    end do
    do i=1,nmax
        do j=i,nmax
            write(1,'(2I4,2E20.12)') j,i,factor*b(k),factor*b(k+nmax+1-i)
            k = k+1
        end do
        k = k+nmax+1-i
    end do
    close(1)

    deallocate(A,N,y,b)
    deallocate(long,lat,rad)
    write(*,*) ' '
    write(*,'(A,F10.3)') 'total runtime [s]:', omp_get_wtime()-timetot
end program shale

subroutine build_a_line(a,i,nobs,nmax,u,long,lat,rad,latold,p,rearth)
    ! build one line of design matrix
    ! input:
    !   a               design matrix
    !   i               line of design matrix
    !   nobs            number of observations
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
    integer*4 :: i, nobs, nmax, u, m, n, j
    real*8 :: a(nobs,u), long, lat, rad, latold
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
            A(i,j) = rrm(n)*cosml*p(n,m)
            j = j+1
        end do
        if(m.gt.0) then
            ! sine terms
            do n=m,nmax
                A(i,j) = rrm(n)*sinml*p(n,m)
                j = j+1
            end do
        end if
    end do
    deallocate(rrm)
end subroutine build_a_line


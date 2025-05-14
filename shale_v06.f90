! SHALE - spherical harmonic analysis with least squares
!
! V0.6 - parallelised with OpenMP and ScaLAPACK, row and column distribution
!
! Tobias Wittwer, 2006
program shale
    implicit none
    ! settings
    integer*4, parameter :: nmax = 100                ! maximum degree/order
    integer*4, parameter :: nobs = 64800                  ! number of observations
    character*80, parameter :: obsfile = 'pot_glob_0_100.txt'   ! observation file
    character*80, parameter :: outfile = 'cnm_100.txt'   ! output file
    integer*4, parameter :: u = nmax**2+2*nmax+1     ! number of unknowns

    ! constants
    real*8, parameter :: pi = 3.141592653589793
    real*8, parameter :: rho = 180.0d0/pi
    real*8, parameter :: GM = 3986005.0d8
    real*8, parameter :: rearth = 6371000.0d0
    
    ! arrays
    real*8, allocatable :: A(:,:),a_row(:)           ! design matrix
    real*8, allocatable :: N(:,:)                    ! normal matrix
    real*8, allocatable :: y(:),obs(:)               ! observation vector
    real*8, allocatable :: b(:)                      ! rhs vector
    real*8, allocatable :: pnm(:,:)                  ! legendre functions
    real*8, allocatable :: long(:), lat(:), rad(:)   ! positions

    ! other
    integer*4 :: i, info, j, k
    real*8 :: latold, factor
    real*8 :: omp_get_wtime, timetot, time

    ! ScaLAPACK
    integer*4 :: iam, nprocs, ictxt, idx, idx2
    integer*4 :: nprow, npcol, myrow, mycol
    integer*4 :: desca(9), descb(9), descn(9), descy(9)
    integer*4 :: ra,ca,rb,cb,rn,cn,ry,cy
    real*8, allocatable :: x(:)
    integer*4, parameter :: blocksize = 64

    interface
        INTEGER*4 FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
            INTEGER*4 :: INDXLOC, IPROC, ISRCPROC, NB, NPROCS
        end FUNCTION INDXL2G
    end interface


    timetot = omp_get_wtime()

    ! initialise BLACS process grid
    call blacs_pinfo(iam,nprocs)
    call blacs_get( -1, 0, ictxt )
    call blacs_gridinit(ictxt,'R',nprocs/4,4)

    if(iam.eq.0) then
        write(*,'(A)') '* SHALE V0.6 *'
        write(*,'(A,I4)') 'spherical harmonic analysis up to degree ', nmax
        write(*,'(A,I6,A)') 'using ', nobs, ' observations'
        write(*,*) ' '
    end if

    ! calculate matrix sizes
    call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
    !write(*,*) iam, myrow,mycol
    call calc_sizes_descinit(blocksize,ictxt,myrow,mycol,nprow,npcol,nobs,u,desca,ra,ca)
    call calc_sizes_descinit(blocksize,ictxt,myrow,mycol,nprow,npcol,u,1,descb,rb,cb)
    call calc_sizes_descinit(blocksize,ictxt,myrow,mycol,nprow,npcol,u,u,descn,rn,cn)
    call calc_sizes_descinit(blocksize,ictxt,myrow,mycol,nprow,npcol,nobs,1,descy,ry,cy)

    ! allocate arrays
    allocate(A(ra,ca),N(rn,cn),y(ry*cy),b(rb*cb))
    allocate(long(nobs),lat(nobs),rad(nobs),obs(nobs))
    allocate(x(u))

    ! data input
    open(1,file=obsfile)
    do i=1,nobs
        read(1,*) long(i), lat(i), rad(i), obs(i)
        long(i) = long(i)/rho
        lat(i) = lat(i)/rho
    end do
    close(1)

    !write(*,*) iam,ry,cy,rb,cb

    ! build design matrix
    time = omp_get_wtime()
    !$OMP PARALLEL PRIVATE(latold,pnm,idx,a_row,idx2,j)
    allocate(pnm(0:nmax,0:nmax),a_row(u))
    latold = 1.6d0
    !$OMP DO
    do i=1,ra
        idx=indxl2g(i, blocksize, myrow, 0, nprow)
        !write(*,*) iam,i,idx
        call build_a_line(a_row,i,nobs,ra,nmax,u,long(idx),lat(idx),rad(idx),latold,pnm,rearth)
        do j=1,ca
            idx2=indxl2g(j,blocksize,mycol,0,npcol)
            !write(*,*) i,idx,j,idx2
            A(i,j) = a_row(idx2)
        end do
        if(cy.gt.0) y(i)=obs(idx)
        latold = lat(idx)
    end do
    !$OMP END DO
    deallocate(pnm,a_row)
    !$OMP END PARALLEL
    if(iam.eq.0) write(*,'(A,F10.3)') 'build A [s]      :', omp_get_wtime()-time
    ! build normal equation matrix and rhs vector
    time = omp_get_wtime()
    call PDSYRK('L','T', u, nobs, 1.0d0, A, 1,1, desca, 0.0d0, N,1,1, descn)
    if(iam.eq.0) write(*,'(A,F10.3)') 'build N [s]      :', omp_get_wtime()-time
    time = omp_get_wtime()
    call PDGEMV('T', nobs, u, 1.0d0, A,1,1,desca, y, 1,1,descy,1, 0.0d0, b,1,1,descb, 1)
    if(iam.eq.0) write(*,'(A,F10.3)') 'build b [s]      :', omp_get_wtime()-time

    ! solve
    time = omp_get_wtime()
    call PDPOSV('L', u, 1, N, 1,1,descn, b,1,1,descb, info)
    if(iam.eq.0) then
        write(*,'(A,F10.3)') 'solving [s]      :', omp_get_wtime()-time
        write(*,'(A,I10)') 'info             :', info
    end if

    ! gather estimated coefficients
    x(:) = 0.0d0
    if(cb.gt.0) then
        do i=1,rb
            idx=indxl2g(i, blocksize, myrow, 0, nprow)
            x(idx) = b(i)
        end do
    end if
    call dgsum2d(ictxt,'A',' ',1,u,x,1,-1,-1)

    if(iam.eq.0) then
        ! write output
        open(1,file=outfile)
        factor = rearth/GM
        k=1
        do j=0,nmax
            write(1,'(2I4,2E20.12)') j,0,factor*x(k),0.0d0
            k = k+1
        end do
        do i=1,nmax
            do j=i,nmax
                write(1,'(2I4,2E20.12)') j,i,factor*x(k),factor*x(k+nmax+1-i)
                k = k+1
            end do
            k = k+nmax+1-i
        end do
        close(1)
        write(*,*) ' '
        write(*,'(A,F10.3)') 'total runtime [s]:', omp_get_wtime()-timetot
    end if

    deallocate(x,A,N,y,b)
    deallocate(long,lat,rad,obs)
    call blacs_gridexit(ictxt)
    call blacs_exit(0)
end program shale

subroutine build_a_line(a,i,nobs,ra,nmax,u,long,lat,rad,latold,p,rearth)
    ! build one line of design matrix
    ! input:
    !   a               design matrix
    !   i               line of design matrix
    !   nobs            number of observations
    !   ra              number of rows in local A
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
    integer*4 :: i, nobs, nmax, u, m, n, j, ra
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
            !A(i,j) = rrm(n)*cosml*p(n,m)
            a(j) = rrm(n)*cosml*p(n,m)
            j = j+1
        end do
        if(m.gt.0) then
            ! sine terms
            do n=m,nmax
                !A(i,j) = rrm(n)*sinml*p(n,m)
                a(j) = rrm(n)*sinml*p(n,m)
                j = j+1
            end do
        end if
    end do
    deallocate(rrm)
end subroutine build_a_line

subroutine calc_sizes_descinit(blocksize,ictxt,myrow,mycol,nprow,npcol,m,n,desc,r,c)
    ! calculate matrix sizes and initialise desc vector
    ! input:
    !  blocksize      blocksize for matrix distribution
    !  ictxt          BLACS communication context
    !  myrow          row of the process
    !  mycol          column of the process
    !  nprow          total number of process rows
    !  npcol          total number of process columns
    !  m              number of rows of the global matrix
    !  n              number of columns of the global matrix
    !  desc           description vector
    !  r              number of rows of the local matrix
    !  c              number of columns of the local matrix
    implicit none
    integer*4 :: m,n,desc(9),info,r,c,ictxt,ld
    integer*4 :: blocksize, myrow, mycol, nprow, npcol
    interface
        INTEGER*4 FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
            INTEGER*4 :: IPROC, ISRCPROC, N, NB, NPROCS
        end FUNCTION NUMROC
      
        INTEGER*4 FUNCTION NPREROC( N, NB, IPROC, ISRCPROC, NPROCS )
            INTEGER*4 :: IPROC, ISRCPROC, N, NB, NPROCS
        end FUNCTION NPREROC
    end interface
    
    r=numroc(m,blocksize,myrow,0,nprow)
    c=numroc(n,blocksize,mycol,0,npcol)
    if(r.eq.0) then
        ld=1
    else
        ld=r
    end if
    call descinit(desc,m,n,blocksize,blocksize,0,0,ictxt,ld,info)
end subroutine calc_sizes_descinit


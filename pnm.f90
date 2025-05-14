subroutine pnm(nmax, lat, p)
    implicit none
    integer*4 :: nmax
    real*8 :: lat,p(0:nmax,0:nmax)
    integer*4 :: l,m,k
    real*8 :: ff1,ff2,fk
    real*8 :: f1(0:1), f2(0:1)
    real*8 :: eta, sigma, tau, ny
    real*8 :: sinlat, coslat, sinpilat, cospilat

    sinlat = dsin(lat)
    coslat = dcos(lat)

    f1(0)=dsqrt(3.0d0)
    f2(0)=dsqrt(5.0d0)
    f1(1)=dsqrt(15.0d0)
    f2(1)=f1(1)/2.0d0

    p(0,0)=1.0d0
    p(1,0)=f1(0)*sinlat
    p(1,1)=f1(0)*coslat
    p(2,0)=f2(0)*(1.5d0*sinlat*sinlat-0.5d0)
    p(2,1)=f1(1)*sinlat*coslat
    P(2,2)=f2(1)*coslat*coslat

    do l=3, nmax
        do m=0,nmax-2
            ff1=dsqrt(dble(4*l*l-1))/dsqrt(dble(l*l-m*m))
            ff2=dsqrt(dble((2*l+1)*(l+m-1)*(l-m-1)))/dsqrt(dble((l*l-m*m)*(2*l-3)))

            p(l,m) = ff1*sinlat*p(l-1,m)-ff2*p(l-2,m)
        end do

        do k=1,2
            fk=dsqrt(dble(2*l+1))/dsqrt(dble(2*(l-k+1)))

            p(l,l+1-k)=fk*coslat*p(l-1,l-k)
        end do
    end do
end subroutine pnm
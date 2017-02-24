program main
!implicit none
!
! test GLL library in SPECFEM2D
! GLL points, weights
! define_derivation_matrices
! written by Yangkang Chen (2017)
! gfortran test_gll_library.f90 && ./a.out
! GLL points are also saved in cykdata/earthquake_packages/SEM_matlab/SEM_matlab/xwh_gll

        double precision, parameter :: alpha = 0.d0, beta = 0.d0
         integer,parameter :: np=5   !number of GLL points, polynomial degree + 1 
        double precision,dimension(np) :: z
        real(kind=4), dimension(np,np) :: H

! function for calculating derivatives of Lagrange polynomials
        double precision, external :: lagrange_deriv_GLL
              
        real(kind=4)  :: w(np)   
!       real(kind=CUSTOM_REAL,dimension(np) :: w      
        call zwgljd(z,w,np,alpha,beta)

        write(*,*) 'Hello'
        write(*,*)"GLL points",z
        write(*,*)"GLL weights",w
        write(*,*)"GLL alpha",alpha
        write(*,*)"GLL beta",beta
        
        ! calculate derivatives of the Lagrange polynomials
        ! and precalculate some products in double precision
        ! hprime(i,j) = h'_j(xgll_i) by definition of the derivation matrix
          do i1=1,np
            do i2=1,np
              H(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,z,np)
            enddo
          enddo        

        write(*,*)"derivation matrix",H
                               
        call timestamp()

             
!        call timestamp()

!        double precision, parameter :: alpha = 0.d0, beta = 0.d0        
!
!         integer,parameter :: np=5   
!        double precision,dimension(np) :: z
!        double precision,dimension(np) :: w        
!        zwgljd(z,w,np,alpha,beta)
               
end program main

subroutine timestamp()
implicit none

        character(8) date
        character(5) time
        real,parameter::d=1

        call date_and_time(date,time)
        write(*,'(a8,2x,a10)') date,time

        
end subroutine



!=======================================================================
!
!  Library to compute the Gauss-Lobatto-Legendre points and weights
!  Based on Gauss-Lobatto routines from M.I.T.
!  Department of Mechanical Engineering
!  https://svn.mcs.anl.gov/repos/nek5/branches/adapt_mesh/speclib.f
!
!=======================================================================

  double precision function endw1(n,alpha,beta)

  implicit none

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  double precision, external :: gammaf
  integer i

  f3 = zero
  apb = alpha+beta
  if(n == 0) then
    endw1 = zero
    return
  endif
  f1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1 = f1*(apb+two)*two**(apb+two)/two
  if(n == 1) then
   endw1 = f1
   return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if(n == 2) then
   endw1 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw1  = f3

  end function endw1

!
!=======================================================================
!

  double precision function endw2(n,alpha,beta)

  implicit none

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  double precision, external :: gammaf
  integer i

  apb   = alpha+beta
  f3 = zero
  if (n == 0) then
   endw2 = zero
   return
  endif
  f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw2 = f1
   return
  endif
  fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw2 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw2  = f3

  end function endw2

!
!=======================================================================
!

  double precision function gammaf (x)

  implicit none

  double precision, parameter :: pi = 3.141592653589793d0

  double precision x

  double precision, parameter :: half=0.5d0,one=1.d0,two=2.d0

  gammaf = one

  if (x == -half) gammaf = -two*sqrt(pi)
  if (x ==  half) gammaf =  sqrt(pi)
  if (x ==  one ) gammaf =  one
  if (x ==  two ) gammaf =  one
  if (x ==  1.5d0) gammaf =  sqrt(pi)/2.d0
  if (x ==  2.5d0) gammaf =  1.5d0*sqrt(pi)/2.d0
  if (x ==  3.5d0) gammaf =  2.5d0*1.5d0*sqrt(pi)/2.d0
  if (x ==  3.d0 ) gammaf =  2.d0
  if (x ==  4.d0 ) gammaf = 6.d0
  if (x ==  5.d0 ) gammaf = 24.d0
  if (x ==  6.d0 ) gammaf = 120.d0

  end function gammaf

!
!=====================================================================
!

  subroutine jacg (xjac,np,alpha,beta)

!=======================================================================
!
! computes np Gauss points, which are the zeros of the
! Jacobi polynomial with parameters alpha and beta
!
!                  .alpha = beta =  0.0  ->  Legendre points
!                  .alpha = beta = -0.5  ->  Chebyshev points
!
!=======================================================================

  implicit none

  integer np
  double precision alpha,beta
  double precision xjac(np)

  integer k,j,i,jmin,jm,n
  double precision xlast,dth,x,x1,x2,recsum,delx,xmin,swap
  double precision p,pd,pm1,pdm1,pm2,pdm2

  integer, parameter :: K_MAX_ITER = 10
  double precision, parameter :: zero = 0.d0, eps = 1.0d-12

  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  xlast = 0.d0
  n   = np-1
  dth = 4.d0*atan(1.d0)/(2.d0*dble(n)+2.d0)
  p = 0.d0
  pd = 0.d0
  jmin = 0
  do j=1,np
   if(j == 1) then
      x = cos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
   else
      x1 = cos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
      x2 = xlast
      x  = (x1+x2)/2.d0
   endif
   do k=1,K_MAX_ITER
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
      recsum = 0.d0
      jm = j-1
      do i=1,jm
         recsum = recsum+1.d0/(x-xjac(np-i+1))
      enddo
      delx = -p/(pd-recsum*p)
      x    = x+delx
      if(abs(delx) < eps) goto 31
   enddo
 31      continue
   xjac(np-j+1) = x
   xlast        = x
  enddo
  do i=1,np
   xmin = 2.d0
   do j=i,np
      if(xjac(j) < xmin) then
         xmin = xjac(j)
         jmin = j
      endif
   enddo
   if(jmin /= i) then
      swap = xjac(i)
      xjac(i) = xjac(jmin)
      xjac(jmin) = swap
   endif
  enddo

  end subroutine jacg

!
!=====================================================================
!

  subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,bet,x)

!=======================================================================
!
! Computes the Jacobi polynomial of degree n and its derivative at x
!
!=======================================================================

  implicit none

  double precision poly,pder,polym1,pderm1,polym2,pderm2,alp,bet,x
  integer n

  double precision apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave
  integer k

  apb  = alp+bet
  poly = 1.d0
  pder = 0.d0
  psave = 0.d0
  pdsave = 0.d0

  if (n == 0) return

  polyl = poly
  pderl = pder
  poly  = (alp-bet+(apb+2.d0)*x)/2.d0
  pder  = (apb+2.d0)/2.d0
  if (n == 1) return

  do k=2,n
    dk = dble(k)
    a1 = 2.d0*dk*(dk+apb)*(2.d0*dk+apb-2.d0)
    a2 = (2.d0*dk+apb-1.d0)*(alp**2-bet**2)
    b3 = (2.d0*dk+apb-2.d0)
    a3 = b3*(b3+1.d0)*(b3+2.d0)
    a4 = 2.d0*(dk+alp-1.d0)*(dk+bet-1.d0)*(2.d0*dk+apb)
    polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
    pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
    psave  = polyl
    pdsave = pderl
    polyl  = poly
    poly   = polyn
    pderl  = pder
    pder   = pdern
  enddo

  polym1 = polyl
  pderm1 = pderl
  polym2 = psave
  pderm2 = pdsave

  end subroutine jacobf

!
!------------------------------------------------------------------------
!

  double precision FUNCTION PNDLEG (Z,N)

!------------------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
  implicit none

  double precision z
  integer n

  double precision P1,P2,P1D,P2D,P3D,FK,P3
  integer k

  P1   = 1.d0
  P2   = Z
  P1D  = 0.d0
  P2D  = 1.d0
  P3D  = 1.d0

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
    P3D = ((2.d0*FK+1.d0)*P2 + (2.d0*FK+1.d0)*Z*P2D - FK*P1D) / (FK+1.d0)
    P1  = P2
    P2  = P3
    P1D = P2D
    P2D = P3D
  enddo

  PNDLEG = P3D

  end function pndleg

!
!------------------------------------------------------------------------
!

  double precision FUNCTION PNLEG (Z,N)

!------------------------------------------------------------------------
!
!     Compute the value of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
  implicit none

  double precision z
  integer n

  double precision P1,P2,P3,FK
  integer k

  P1   = 1.d0
  P2   = Z
  P3   = P2

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
    P1  = P2
    P2  = P3
  enddo

  PNLEG = P3

  end function pnleg

!
!------------------------------------------------------------------------
!

  double precision function pnglj(z,n)

!------------------------------------------------------------------------
!
!     Compute the value of the Nth order polynomial of the
!     Gauss-Lobatto-Jacobi (0,1) at Z. from Legendre polynomials.
!
!------------------------------------------------------------------------

  implicit none
  include 'constants.h'

  double precision z
  integer n
  double precision, external :: pnleg

  if (abs(z+1.d0) > TINYVAL) then  ! if (z /= -1.d0)
    pnglj = (pnleg(z,n)+pnleg(z,n+1))/(ONE+z)
  else
    pnglj = (dble(n)+ONE)*(-1)**n
  endif

  end function pnglj

!
!------------------------------------------------------------------------
!

  double precision function pndglj(z,n)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of Nth order polynomial of the
!     Gauss-Lobatto-Jacobi (0,1) at Z. from Legendre polynomials.
!
!------------------------------------------------------------------------

  implicit none
  include 'constants.h'

  double precision z
  integer n
  double precision, external :: pnleg, pndleg

  if (abs(z+1.d0) > TINYVAL) then  ! if (z /= -1.d0)
    pndglj = (pndleg(z,n)+pndleg(z,n+1))/(ONE+z) - (pnleg(z,n)+pnleg(z,n+1))/((ONE+z)**2)
  else
    pndglj = pnleg(-1.d0,n)+pnleg(-1.d0,n+1)
  endif

  end function pndglj

!
!------------------------------------------------------------------------
!

  double precision function pnormj (n,alpha,beta)

  implicit none

  double precision alpha,beta
  integer n

  double precision one,two,dn,const,prod,dindx,frac
  double precision, external :: gammaf
  integer i

  one   = 1.d0
  two   = 2.d0
  dn    = dble(n)
  const = alpha+beta+one

  if (n <= 1) then
    prod   = gammaf(dn+alpha)*gammaf(dn+beta)
    prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
    pnormj = prod * two**const/(two*dn+const)
    return
  endif

  prod  = gammaf(alpha+one)*gammaf(beta+one)
  prod  = prod/(two*(one+const)*gammaf(const+one))
  prod  = prod*(one+alpha)*(two+alpha)
  prod  = prod*(one+beta)*(two+beta)

  do i=3,n
    dindx = dble(i)
    frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
    prod  = prod*frac
  enddo

  pnormj = prod * two**const/(two*dn+const)

  end function pnormj

!
!------------------------------------------------------------------------
!

  subroutine zwgjd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g j d : Generate np Gauss-Jacobi points and weights
!                 associated with Jacobi polynomial of degree n = np-1
!
!     Note : Coefficients alpha and beta must be greater than -1.
!     ----
!=======================================================================

  implicit none
  include 'constants.h'

  !double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision z(np)
  real(kind=CUSTOM_REAL)  :: w(np)
  double precision alpha,beta

  integer n,np1,np2,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef
  double precision, external :: gammaf,pnormj

  pd = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n    = np-1
  apb  = alpha+beta
  p    = zero
  pdm1 = zero

  if (np == 1) then
   z(1) = (beta-alpha)/(apb+two)
   w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
   return
  endif

  call jacg(z,np,alpha,beta)

  np1   = n+1
  np2   = n+2
  dnp1  = dble(np1)
  dnp2  = dble(np2)
  fac1  = dnp1+alpha+beta+one
  fac2  = fac1+dnp1
  fac3  = fac2+one
  fnorm = pnormj(np1,alpha,beta)
  rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
  do i=1,np
    call jacobf(p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
    w(i) = -rcoef/(p*pdm1)
  enddo

  end subroutine zwgjd

!
!------------------------------------------------------------------------
!

  subroutine zwgljd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!     Note : alpha and beta coefficients must be greater than -1.
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================

  implicit none
  include 'constants.h'


  !double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision alpha,beta
  double precision z(np)
  real(kind=CUSTOM_REAL)  :: w(np)

  integer n,nm1,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision alpg,betg
  double precision, external :: endw1,endw2

  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero


  if (nm1 > 0) then
    alpg  = alpha+one
    betg  = beta+one
    call zwgjd(z(2),w(2),nm1,alpg,betg)
  endif

  z(1)  = - one
  z(np) =  one

  do i=2,np-1
   w(i) = w(i)/(one-z(i)**2)
  enddo

  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1(n,alpha,beta)/(two*pd)
  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2(n,alpha,beta)/(two*pd)

  end subroutine zwgljd



! The following is from src/specfem2d/lagrange_poly.f90
!========================================================================

  double precision function hgll(I,Z,ZGLL,NZ)

!-------------------------------------------------------------
!
!  Compute the value of the Lagrangian interpolant L through
!  the NZ Gauss-Lobatto Legendre points ZGLL at point Z
!
!-------------------------------------------------------------

  implicit none

  integer i,nz
  double precision z
  double precision ZGLL(0:nz-1)

  integer n
  double precision EPS,DZ,ALFAN
  double precision, external :: PNLEG,PNDLEG

  EPS = 1.d-5
  DZ = Z - ZGLL(I)
  if(abs(DZ) < EPS) then
   HGLL = 1.d0
   return
  endif
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+1.d0)
  HGLL = - (1.d0-Z*Z)*PNDLEG(Z,N)/ (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))

  end function hgll

!
!=====================================================================
!

  double precision function hglj(I,Z,ZGLJ,NZ)

!-------------------------------------------------------------
!
!  Compute the value of the Lagrangian interpolant L through
!  the NZ Gauss-Lobatto Jacobi points ZGLJ at point Z
!
!-------------------------------------------------------------

  implicit none

  integer i,nz
  double precision z
  double precision ZGLJ(0:nz-1)

  integer n
  double precision EPS,DZ,ALFAN
  double precision, external :: PNGLJ,PNDGLJ

  EPS = 1.d-5
  DZ = Z - ZGLJ(I)
  if(abs(DZ) < EPS) then
   HGLJ = 1.d0
   return
  endif
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+2.d0)
  HGLJ = - (1.d0-Z*Z)*PNDGLJ(Z,N)/ (ALFAN*PNGLJ(ZGLJ(I),N)*(Z-ZGLJ(I)))

  end function hglj

!
!=====================================================================
!

  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  integer NGLL
  double precision xi,xigll(NGLL),h(NGLL),hprime(NGLL)

  integer dgr,i,j
  double precision prod1,prod2

  do dgr=1,NGLL

    prod1 = 1.0d0
    prod2 = 1.0d0
    do i=1,NGLL
      if(i /= dgr) then
        prod1 = prod1*(xi-xigll(i))
        prod2 = prod2*(xigll(dgr)-xigll(i))
      endif
    enddo
    h(dgr)=prod1/prod2

    hprime(dgr)=0.0d0
    do i=1,NGLL
      if(i /= dgr) then
        prod1=1.0d0
        do j=1,NGLL
          if(j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
        enddo
        hprime(dgr) = hprime(dgr)+prod1
      endif
    enddo
    hprime(dgr) = hprime(dgr)/prod2

  enddo

  end subroutine lagrange_any

!
!=====================================================================
!

! subroutine to compute the derivative of the Lagrange interpolants
! at the GLL points at any given GLL point

  double precision function lagrange_deriv_GLL(I,j,ZGLL,NZ)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrange interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)
!
!------------------------------------------------------------------------

  implicit none

  integer i,j,nz
  double precision zgll(0:nz-1)

  integer degpoly

  double precision, external :: pnleg,pndleg

  degpoly = nz - 1
  if (i == 0 .and. j == 0) then
    lagrange_deriv_GLL = - dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == degpoly .and. j == degpoly) then
    lagrange_deriv_GLL = dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == j) then
    lagrange_deriv_GLL = 0.d0
  else
    lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
      (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
      + (1.d0-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) / (dble(degpoly)* &
      (dble(degpoly)+1.d0)*pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  end function lagrange_deriv_GLL

!
!=======================================================================
!

! subroutine to compute the derivative of the interpolants of the GLJ
! quadrature at the GLJ points at any given GLJ point

!  double precision function poly_deriv_GLJ(I,j,ZGLJ,NZ)
!
!!------------------------------------------------------------------------
!!
!!     Compute the value of the derivative of the I-th
!!     polynomial interpolant of the GLJ quadrature through the
!!     NZ Gauss-Lobatto-Jacobi (0,1) points ZGLJ at point ZGLJ(j)
!!
!!------------------------------------------------------------------------
!
!  implicit none
!
!  integer i,j,nz
!  double precision zglj(0:nz-1)
!
!  integer degpoly
!
!  double precision, external :: pnglj
!
!  degpoly = nz - 1
!
!  if (i == 0 .and. j == 0) then ! Case 1
!    poly_deriv_GLJ = -dble(degpoly)*(dble(degpoly)+2.d0)/6.d0
!  else if (i == 0 .and. 0 < j .and. j < degpoly) then ! Case 2
!    poly_deriv_GLJ = 2.d0*(-1)**degpoly*pnglj(zglj(j),degpoly)/((1.d0+zglj(j))*(dble(degpoly)+1.d0))
!  else if (i == 0 .and. j == degpoly) then ! Case 3
!    poly_deriv_GLJ = (-1)**degpoly/(dble(degpoly)+1.d0)
!  else if (0 < i .and. i < degpoly .and. j == 0) then ! Case 4
!    poly_deriv_GLJ = (-1)**(degpoly+1)*(dble(degpoly)+1.d0)/(2.d0*pnglj(zglj(i),degpoly)*(1.d0+zglj(i)))
!  else if (0 < i .and. i < degpoly .and. 0 < j .and. j < degpoly .and. i /= j) then ! Case 5
!    poly_deriv_GLJ = 1.d0/(zglj(j)-zglj(i))*pnglj(zglj(j),degpoly)/pnglj(zglj(i),degpoly)
!  else if (0 < i .and. i < degpoly .and. i == j) then  ! Case 6
!    poly_deriv_GLJ = -1.d0/(2.d0*(1.d0+zglj(i)))
!  else if (0 < i .and. i < degpoly .and. j == degpoly) then ! Case 7
!    poly_deriv_GLJ = 1.d0/(pnglj(zglj(i),degpoly)*(1.d0-zglj(i)))
!  else if (i == degpoly .and. j == 0) then ! Case 8
!    poly_deriv_GLJ = (-1)**(degpoly+1)*(dble(degpoly)+1.d0)/4.d0
!  else if (i == degpoly .and. 0 < j .and. j < degpoly) then ! Case 9
!    poly_deriv_GLJ = -1.d0/(1.d0-zglj(j))*pnglj(zglj(j),degpoly)
!  else if (i == degpoly .and. j == degpoly) then ! Case 10
!    poly_deriv_GLJ = (dble(degpoly)*(dble(degpoly)+2.d0)-1.d0)/4.d0
!  else
!    call exit_MPI('Problem in poly_deriv_GLJ: in a perfect world this would NEVER appear')
!  endif
!
!  end function poly_deriv_GLJ

































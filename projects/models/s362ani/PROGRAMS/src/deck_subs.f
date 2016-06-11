	subroutine read_deck (n,nic,noc,nmoho,nsl,r,rho,
     #      vph,vpv,vsh,vsv,eta,
     #		qro,qshear,qkappa,vphspl,
     #      vpvspl,vshspl,vsvspl,etaspl,
     #		qmiuspl,qkapspl,refstr,iin,iout)

	implicit none

	integer maxdeck
	parameter(maxdeck=1000)
	character*(*) refstr
	integer*4 ititle(20)
	real*8 vpv(maxdeck),vph(maxdeck),vsv(maxdeck),qkappa(maxdeck)
	real*8 vsh(maxdeck),eta(maxdeck),qshear(maxdeck)
	real*8 r(maxdeck),rho(maxdeck),qro(3,maxdeck)
	real*8 vphspl(3,maxdeck),vpvspl(3,maxdeck),vshspl(3,maxdeck)
	real*8 vsvspl(3,maxdeck),etaspl(3,maxdeck)
	real*8 qmiuspl(3,maxdeck),qkapspl(3,maxdeck)
	real*8 pi,tref,rx,r1,r2,dr,rt,wrk(10000),val
	integer ifanis
	integer nic,noc,n,nsl
	integer iin,iout,i,ifdeck,nmoho,nreg,knt
	integer jj,nn,nlay,j,ind,k
	logical exists
c
	integer nlevout
	real*4 depmin,depmax
c
	real*8 tau,bigg
      	data bigg,tau/6.6723d-11,1.d3/

        inquire(file=refstr,exist=exists)
	if(.not.exists) then
	  write(6,"('reference model file does not exist')")
	  stop
	endif

	open(iin,file=refstr,status='old')

c --- this program should also work for polynomial prem
c	open(iin,file='aniprmc',status='old')

	open(iout,file='premout.txt')

      pi=3.14159265358979d0
      read(iin,100) (ititle(i),i=1,20)
  100 format(20a4)
      read(iin,*) ifanis,tref,ifdeck
      if(ifdeck.eq.0) go to 1000
c*** card deck model ***
      read(iin,*) n,nic,noc,nmoho,nsl
      read(iin,*) (r(i),rho(i),vpv(i),vsv(i),
     +     qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
	if(i.gt.maxdeck) stop 'you need larger maxdeck in deck_subs.f'
  105 format(f8.0,3f9.2,2f9.1,2f9.2,f9.5)
      go to 2000
c*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      rx=rx*tau
      n=0
      knt=0
      jj=5
      if(ifanis.ne.0) jj=8
      do 10 nn=1,nreg
      read(iin,*) nlay,r1,r2
      r1=r1*tau
      r2=r2*tau
      dr=(r2-r1)/float(nlay-1)
      do 15 i=1,nlay
      n=n+1
   15 r(n)=r1+dr*float(i-1)
      do 20 j=1,jj
      read(iin,*) (wrk(i),i=1,5)
      do 20 i=1,nlay
      ind=knt+i
      rt=r(ind)/rx
      val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
      if(j.eq.1) rho(ind)=val*tau
      if(j.eq.2) vpv(ind)=val*tau
      if(j.eq.3) vsv(ind)=val*tau
c      if(j.eq.4) qshear(ind)=val
c      if(j.eq.5) qkappa(ind)=val
      if(j.eq.4) qshear(ind)=1.0/val
      if(j.eq.5) qkappa(ind)=1.0/val
c
      if(ifanis.eq.0) goto 20
      if(j.eq.6) vph(ind)=val*tau
      if(j.eq.7) vsh(ind)=val*tau
      if(j.eq.8) eta(ind)=val
   20 continue
   10 knt=knt+nlay
 2000 if(ifanis.ne.0) go to 3000
      do 25 i=1,n
      vph(i)=vpv(i)
      vsh(i)=vsv(i)
   25 eta(i)=1.d0
 3000 continue
c-----------------------------------------------------------------
         nlevout=0
      do i=2,n
        if(r(n)-r(i).gt.dble(depmin*1000.-500.).and.
     #		r(n)-r(i).lt.dble(depmax*1000.+500.)) then
         nlevout=nlevout+1
        endif
      enddo
      if(nlevout.gt.0) nlevout=nlevout+1
c      write(6,"('depmin,depmax,nlevout:',2g15.5,i10)") depmin,depmax,nlevout
c-----------------------------------------------------------------
      if(iout.lt.0) goto 30
c*** write out model ***
       write(iout,900) (ititle(k),k=1,20),tref
900    format(1x,20a4,' ref per =',f6.1,' secs',///,2x,'level',
     1 4x,'radius',8x,'rho',9x,'vpv',9x,'vph',9x,'vsv',
     2 9x,'vsh',9x,'eta',9x,'qmu ',8x,'qkap',/)
       write(iout,905) (i,r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
     1 eta(i),qshear(i),qkappa(i),i=1,n)
  905 format(3x,i3,f12.1,5f12.2,f12.5,2f12.2)

c --- spline

30  	do 45 i=1,n
	if(i.gt.1.and.dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
	if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
45	if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)

	call drspln(1,n,r,vph,vphspl,wrk)
	call drspln(1,n,r,vpv,vpvspl,wrk)
	call drspln(1,n,r,vsh,vshspl,wrk)
	call drspln(1,n,r,vsv,vsvspl,wrk)
	call drspln(1,n,r,eta,etaspl,wrk)
	call drspln(1,n,r,rho,qro,wrk)
	call drspln(1,n,r,qshear,qmiuspl,wrk)
	call drspln(1,n,r,qkappa,qkapspl,wrk)

	close(iin)
	close(iout)
      end


c ------------
c --- return one parameter at dep[km] in model from read_deck

c --- ityp=1-rho,2-vpv,3-vsv,4-qk,5-qm,
c --- 6-vph,7-vsh,8-eta,9-isoP,10-isoS
c --- in the output velocities will be in km/s

	subroutine eval_deck (ityp,dep,n,nic,noc,r,rho,
     #      vph,vpv,vsh,vsv,eta,
     #		qro,qmiu,qkap,vphspl,
     #      vpvspl,vshspl,vsvspl,etaspl,
     # 		qmiuspl,qkapspl,velout)

	implicit none
	integer maxdeck
	parameter(maxdeck=1000)	
	real*8 r(maxdeck),rho(maxdeck),qro(3,maxdeck)
	real*8 qmiu(maxdeck),qkap(maxdeck)
	real*8 vphspl(3,maxdeck),vpvspl(3,maxdeck),vshspl(3,maxdeck)
	real*8 vsvspl(3,maxdeck),etaspl(3,maxdeck)
	real*8 qmiuspl(maxdeck),qkapspl(maxdeck)
	real*8 xa,xc,xf,xn,xl,dep,rr,velout
	integer ityp,n,nic,noc
	real*8 vpv(maxdeck),vph(maxdeck),vsv(maxdeck)
	real*8 vsh(maxdeck),eta(maxdeck)
	real*8 xkappa,xmiu
	real*8 xrho,vvsh,vvsv,vvph,vvpv,vveta
	real*8 drsple

	if(n.gt.maxdeck) stop 'you need larger maxdeck in deck_subs.f'
c	if(ityp.eq.4.or.ityp.eq.5) stop 'elast_ref is not for q'

	rr=6371000.d0-dep*1000.d0

	if(ityp.eq.1) then
	 velout=drsple(1,n,r,rho,qro,rr) 
	else if(ityp.eq.2) then
	 velout=drsple(1,n,r,vpv,vpvspl,rr) 
	else if(ityp.eq.3) then
	 velout=drsple(1,n,r,vsv,vsvspl,rr) 
	else if(ityp.eq.4) then
	 velout=drsple(1,n,r,qkap,qkapspl,rr)*1000.d0 
	else if(ityp.eq.5) then
	 velout=drsple(1,n,r,qmiu,qmiuspl,rr)*1000.d0 
	else if(ityp.eq.6) then
	 velout=drsple(1,n,r,vph,vphspl,rr) 
	else if(ityp.eq.7) then
	 velout=drsple(1,n,r,vsh,vshspl,rr) 
	else if(ityp.eq.8) then
	 velout=drsple(1,n,r,eta,etaspl,rr)*1000.d0
	else if(ityp.eq.9) then
	 xrho=drsple(1,n,r,rho,qro,rr) 
	 vvsv=drsple(1,n,r,vsv,vsvspl,rr) 
	 vvsh=drsple(1,n,r,vsh,vshspl,rr) 
	 vvph=drsple(1,n,r,vph,vphspl,rr) 
	 vvpv=drsple(1,n,r,vpv,vpvspl,rr) 
	 vveta=drsple(1,n,r,eta,etaspl,rr) 
	 xa=xrho*vvph*vvph
	 xc=xrho*vvpv*vvpv
	 xn=xrho*vvsh*vvsh
	 xl=xrho*vvsv*vvsv
	 xf=vveta*(xa-2.d0*xl)
	 xkappa=(4.d0*(xa+xf-xn)+xc)/9.d0
	 xmiu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
	 velout=dsqrt((xkappa+4.d0*xmiu/3.d0)/xrho)
	else if(ityp.eq.10) then
	 xrho=drsple(1,n,r,rho,qro,rr) 
	 vvsv=drsple(1,n,r,vsv,vsvspl,rr) 
	 vvsh=drsple(1,n,r,vsh,vshspl,rr) 
	 vvph=drsple(1,n,r,vph,vphspl,rr) 
	 vvpv=drsple(1,n,r,vpv,vpvspl,rr) 
	 vveta=drsple(1,n,r,eta,etaspl,rr) 
	 xa=xrho*vvph*vvph
	 xc=xrho*vvpv*vvpv
	 xn=xrho*vvsh*vvsh
	 xl=xrho*vvsv*vvsv
	 xf=vveta*(xa-2.d0*xl)
	 xmiu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
 	 velout=dsqrt(xmiu/xrho)
	else 
	 stop 'this ityp is not allowed in elast_ref'
	endif

	velout=velout/1000.d0
	end

c ------------
c ---  return all parameters at dep[km] in model from read_deck
c --- in the output: velocities in [km/s]

	subroutine eval_all (dep,n,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #		qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,
     # 		qmiuspl,qkapspl,
     #		xrho,vvpv,vvsv,vvqk,vvqm,vvph,vvsh,vveta,vvp,vvs)

	implicit none
	integer maxdeck
	parameter(maxdeck=1000)	
	real*8 r(maxdeck),rho(maxdeck),qro(3,maxdeck)
	real*8 qmiu(maxdeck),qkap(maxdeck)
	real*8 vphspl(3,maxdeck),vpvspl(3,maxdeck),vshspl(3,maxdeck)
	real*8 vsvspl(3,maxdeck),etaspl(3,maxdeck)
	real*8 qmiuspl(maxdeck),qkapspl(maxdeck)
	real*8 xa,xc,xf,xn,xl,dep,rr
	integer n,nic,noc
	real*8 vpv(maxdeck),vph(maxdeck),vsv(maxdeck)
	real*8 vsh(maxdeck),eta(maxdeck)
	real*8 xkappa,xmiu
	real*8 xrho,vvsh,vvsv,vvph,vvpv,vveta,vvp,vvs,vvqk,vvqm
	real*8 drsple

	if(n.gt.maxdeck) stop 'you need larger maxdeck in deck_subs.f'
	rr=6371000.d0-dep*1000.d0

	xrho=drsple(1,n,r,rho,qro,rr)/1000.d0
	vvpv=drsple(1,n,r,vpv,vpvspl,rr)/1000.d0 
	vvsv=drsple(1,n,r,vsv,vsvspl,rr)/1000.d0 
	vvqk=drsple(1,n,r,qkap,qkapspl,rr) 
	vvqm=drsple(1,n,r,qmiu,qmiuspl,rr)
	vvph=drsple(1,n,r,vph,vphspl,rr)/1000.d0
	vvsh=drsple(1,n,r,vsh,vshspl,rr)/1000.d0 
	vveta=drsple(1,n,r,eta,etaspl,rr)

	xa=xrho*vvph*vvph
	xc=xrho*vvpv*vvpv
	xn=xrho*vvsh*vvsh
	xl=xrho*vvsv*vvsv
	xf=vveta*(xa-2.d0*xl)
	xkappa=(4.d0*(xa+xf-xn)+xc)/9.d0
	xmiu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
	vvp=dsqrt((xkappa+4.d0*xmiu/3.d0)/xrho)
        vvs=dsqrt(xmiu/xrho)
	end




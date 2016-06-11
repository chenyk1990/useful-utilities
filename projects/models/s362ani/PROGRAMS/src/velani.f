c Evaluate vsh, vsv, Voigt average, and vsh-vsv 
c in km/s (if absolute) , or in per cent (if relative).
c For vsh and vsv evaluations, relative to vsh and vsv
c in the reference model, respectively.
c For both Voigt-average and vsh-vsv evaluations, relative
c to the Voigt average in the reference model.
c Evaluate velocity and anisotropy only within the selected
c rectangular region of interest.  
c Remove the average of all evaluations, if requested.

	implicit none

	integer maxker,maxl,maxcoe,maxver,maxhpa,maxdep
	real twopi

	parameter (maxker=100)		! --- radial kernels
	parameter (maxl=72)		! --- sph harmonics
	parameter (maxcoe=2000)		! --- sph splines
	parameter (maxver=1500)		! --- contibuting splines

	parameter (maxhpa=2)		
        parameter (twopi=6.2831853)
	parameter (maxdep=1000)		! --- max depth points

c --- radial kernels

	character*80 kerstr
	integer iker,numker,ierror
	real*4 vercof(maxker)		! --- radial splines 
	real*4 vercofd(maxker)
	character*80 hsplfl(maxhpa)
	character*40 dskker(maxker)
	character*40 varstr(maxker)	! --- description of the variable
	character*40 vstr
	integer ivarkern(maxker)

c --- horizontal kernels

	integer numcof,numhpa,ihpa,i
	integer itypehpa(maxhpa)
	integer ihpakern(maxker)
	integer numcoe(maxhpa)
	real*4 xlat,xlon

c --- array of pixels

	integer mxpx
	parameter(mxpx=264800)
	real*4 xlatarr(mxpx)
	real*4 xlonarr(mxpx)
	real*4 areaarr(mxpx)
	integer ipx,npx

c --- spherical splines

	integer nconpt(maxhpa),iver
	integer iconpt(maxver,maxhpa)
	real*4 conpt(maxver,maxhpa)	
	real*4 xlaspl(maxcoe,maxhpa)		! --- lat of vertices
	real*4 xlospl(maxcoe,maxhpa)		! --- lon of vertices
	real*4 radspl(maxcoe,maxhpa)		! --- rad of vertices

c --- spherical harmonics

	real*4 ylmcof((maxl+1)**2,maxhpa)
	real*4 wk1(maxl+1)
	real*4 wk2(maxl+1)
	real*4 wk3(maxl+1)
	integer lmax,nylm
	integer lmxhpa(maxhpa)

c --- model coefficients

	real*4 coe(maxcoe,maxker)

c --- areas used to weight averages

	real*4 area,totarea,arealatzone

c --- card deck reference model 

	integer maxdeck
	parameter (maxdeck=1000)	! --- layers in card deck ref model
	real*8 r(maxdeck),rho(maxdeck),qro(3,maxdeck)
	real*8 vphspl(3,maxdeck),vpvspl(3,maxdeck),vshspl(3,maxdeck)
	real*8 vsvspl(3,maxdeck),etaspl(3,maxdeck)
	real*8 qmiuspl(3,maxdeck),qkapspl(3,maxdeck)
	real*8 vpv(maxdeck),vph(maxdeck),vsv(maxdeck)
	real*8 vsh(maxdeck),eta(maxdeck)
	real*8 qmiu(maxdeck),qkap(maxdeck)
	integer ndeck,nic,noc,nmoho,nsl
	character*80 refmdl		! --- reference model name
	integer iin,iout		! --- lu for read_deck 

c --- model evaluation

	integer ish	! --- 0 if SV, 1 if SH
	integer ieval	! --- 1 for velocity, 2 for anisotropy
	real*4 valu(2)	! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
	real*4 value	! --- used in single evaluation of perturbation
	integer isel	! --- if variable should be included
	integer irel	! --- 0 if relative pert, 1 if absolute velocity
	real*4 depth	! --- depth
	real*8 dep8	! --- double precision depth
	real*4 x,y	! --- lat lon
	integer ityp	! --- defines variable in the reference model
	real*4 vsh3drel ! --- relative perturbation
	real*4 vsv3drel ! --- relative perturbation

c --- velocity, anisotropy

	real*4 vsh3d(mxpx)	
	real*4 vsv3d(mxpx)	
	real*4 aniso(mxpx)
	real*4 avervsh,avervsv,averani
	real*4 vsv1d	! --- vsv in a reference model
	real*4 vsh1d	! --- vsh in a reference model
	real*4 vph1d,vpv1d,eta1d,vpterm,avervg

c --- Voight average 

	real*4 vg1d
	real*4 vg3d(mxpx)	

c --- region of interest and grid spacing

	real*4 size	! --- size of a pixel 
	real*4 xminla	! --- definition of region of interest
	real*4 xmaxla
	real*4 xminlo
	real*4 xmaxlo

c --- 

	character*128 string
	character*128 model
	character*128 modeldef
	data modeldef/' '/
	logical exists
	logical gotmodel
	integer itpspl(maxcoe,maxhpa)
	integer lu,imap,numvar,lstr
	integer lnblnk
	integer luout
	character*8 remav		! --- remove average
	integer ifremav	     		! --- remove average
	real*8 vel8
        real*8 xrad

c -------------------------------------

	lu=1 			! --- input : 3-D model 
	luout=2			! --- output: a map.vg.gmt
	iin=10			! -- read_deck
	iout=11			! -- read_deck
        xrad=3.14159265358979/180.d0

c --- read the model 

	imap=1
	gotmodel=.false.
	modeldef=''
	do while(.not.gotmodel)
	  if(lnblnk(modeldef).eq.0) call getcurrmodels(modeldef,imap)
            write(6,"('change earth model [',a,']:')")
     #        modeldef(1:lnblnk(modeldef))
          read(5,"(a)") string
          if(lnblnk(string).ne.0.or.modeldef.ne.string) then
	    if(lnblnk(string).eq.0) string=modeldef
            inquire(file=string,exist=exists)
            if(exists) then
              call gt3dmodl(lu,string,
     #        maxhpa,maxker,maxcoe,
     #        numhpa,numker,numcoe,lmxhpa,
     #        ihpakern,itypehpa,coe,
     #        itpspl,xlaspl,xlospl,radspl,
     #        numvar,ivarkern,varstr,
     #        refmdl,kerstr,hsplfl,dskker,ierror) 
c   	      write(6,"(i2,' variables')") numvar
c              do i=1,numvar
c                write(6,"(i2,1x,a)") i,varstr(i)
c              enddo
c              do i=1,numker
c                write(6,"(i2,i3)") i,ivarkern(i)
c              enddo
c              do i=1,numhpa
c                write(6,"(i2,i5,i3)") i,numcoe(i),lmxhpa(i)
c              enddo
            if(ierror.eq.0) then
	      call putcurrmodels(string,imap)
              model=string
              modeldef=model
              gotmodel=.true.
            endif
          endif
	  else
            gotmodel=.true.
	  endif
	enddo

c --- check arrays

	if(numker.gt.maxker) stop 'numker.gt.maxker'
        do ihpa=1,numhpa
         if(itypehpa(ihpa).eq.1) then
	  if(lmxhpa(ihpa).gt.maxl) stop 'lmxhpa(ihpa).gt.maxl'
	 else if(itypehpa(ihpa).eq.2) then
 	  if(numcoe(ihpa).gt.maxcoe) stop 'numcoe(ihpa).gt.maxcoe'
	 else 
	  stop 'plot.f: problem with itypehpa'
	 endif
	enddo

c --- select if absolute or relative vsh and vsv should be plotted

	write(6,"('evaluate absolute (=0), or relative (=1) ',$)")
	read(5,*),irel
	if(irel.eq.0) then
	  print*,'selected absolute velocity v'
	else if(irel.eq.1) then
	  print*,'selected relative velocity dv/v'
	else
	  stop 'irel must be =0 or 1'
	endif

c --- read reference model

	write(6,"('select reference model:',$)")
	read(5,"(a)") refmdl
	write(6,"('reading reference model: ',a)")
     #refmdl(1:lnblnk(refmdl))
	call read_deck (ndeck,nic,noc,nmoho,nsl,r,rho,vph,vpv,vsh,vsv,eta,
     #		qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,
     #		qmiuspl,qkapspl,refmdl,iin,iout)

c --- define a rectangle; a region of interest

	write(6,"('pixel size (grid spacing): ',$)")
	read(5,*)size
	if(size.lt.0.5.or.size.gt.100) stop 'size.lt.0.5.or.size.gt.100'
	write(6,"('region of interest: xminla,xmaxla,xminlo,xmaxlo ')")
	read(5,*)xminla,xmaxla,xminlo,xmaxlo
c	xminla=-90.0
c	xmaxla=90.0
c	xminlo=0.0
c	xmaxlo=360.0
	if(xminla.ge.xmaxla) stop 'xminla.ge.xmaxla'
	if(xminlo.ge.xmaxlo) stop 'xminlo.ge.xmaxlo'

c --- remove the average?

	write(6,"('remove average in region of interest? (y or n) ',$)")
	read(5,"(a)") remav
	if(remav(1:1).eq.'y'.or.remav(1:3).eq.'yes'
     #.or.remav(1:1).eq.'Y') then
	  ifremav=1
c	  print*,'will remove the average'
	else if(remav(1:1).eq.'n'.or.remav(1:3).eq.
     #'no'.or.remav(1:1).eq.'N') then
	  ifremav=0
c	  print*,'will not remove the average'
	else
	  stop 'do not know if average should be removed'
	endif

c --- find pixels for evaluation

      ipx=0
      xlat=xminla-size
      do while(xlat.lt.xmaxla)
        xlat=xlat+size 	
        xlon=xminlo-size
        do while(xlon.lt.xmaxlo)
          xlon=xlon+size	
          ipx=ipx+1
          if(ipx.gt.mxpx) stop 'ipx.gt.mxpx'
          xlatarr(ipx)=xlat
          xlonarr(ipx)=xlon	
        enddo
      enddo
      npx=ipx
c      write(6,"('number of pixels= ',i)") npx

c --- find areas used to calculate the spherical average	

      if(ifremav.eq.1) then
	  totarea=0.0	 
	  do ipx=1,npx
	    xlat=xlatarr(ipx)
	    xlon=xlonarr(ipx)
	    arealatzone=sin(xrad*(xlat+size/2.0))-
     #sin(xrad*(xlat-size/2.0))
   	    area=twopi*arealatzone*size/360.0
	    areaarr(ipx)=area
	    totarea=totarea+area
	  enddo
      endif

c --- select depth, evaluate reference model 
c --- ityp=1-rho,2-vpv,3-vsv,4-qk,5-qm,6-vph,7-vsh,8-eta,9-isoP,10-isoS

	write(6,"('select depth ')")
	read(5,*) depth

	dep8=depth
	ityp=7
	call eval_deck(ityp,dep8,ndeck,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #	       qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,qmiuspl,
     #         qkapspl,vel8)
	vsh1d=sngl(vel8)

	ityp=3
	call eval_deck(ityp,dep8,ndeck,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #	       qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,qmiuspl,
     #         qkapspl,vel8)
	vsv1d=sngl(vel8)

	ityp=6
	call eval_deck(ityp,dep8,ndeck,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #	       qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,qmiuspl,
     #         qkapspl,vel8)
	vph1d=sngl(vel8)

	ityp=2
	call eval_deck(ityp,dep8,ndeck,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #	       qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,qmiuspl,
     #         qkapspl,vel8)
	vpv1d=sngl(vel8)

	ityp=8
	call eval_deck(ityp,dep8,ndeck,nic,noc,r,rho,vph,vpv,vsh,vsv,eta,
     #	       qro,qmiu,qkap,vphspl,vpvspl,vshspl,vsvspl,etaspl,qmiuspl,
     #         qkapspl,vel8)
	eta1d=sngl(vel8)

	write(6,"('evaluating the model ...')")

c --- evaluate radial basis functions 

	call evradker (depth,kerstr,numker,vercof,vercofd,ierror)
	if(ierror.ne.0) stop 'ierror 3423'

c --- loop over pixels to evaluate the model

	do ipx=1,npx 

	    xlat=xlatarr(ipx)
	    xlon=xlonarr(ipx)
	    area=areaarr(ipx)

c --- loop over evaluations (sv=0,sh=1)

            do ish=0,1			

c             --- contributing horizontal basis functions at xlat,xlon 
	  
	      y=xlat
	      x=xlon
              do ihpa=1,numhpa
                if(itypehpa(ihpa).eq.1) then
                  lmax=lmxhpa(ihpa)
                  call ylm(y,x,lmax,ylmcof(1,ihpa),wk1,wk2,wk3)
                else if(itypehpa(ihpa).eq.2) then
                  numcof=numcoe(ihpa)
                  call splcon(y,x,numcof,xlaspl(1,ihpa),
     #            xlospl(1,ihpa),radspl(1,ihpa),
     #            nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
                else
                  write(6,"('problem 1')") 
                endif
              enddo

c             --- evaluate 3-D perturbations in velocity and anisotropy 
          
	      valu(1)=0. ! --- velocity 
	      valu(2)=0. ! --- anisotropy

	      do ieval=1,2
	        value=0.
	        do iker=1,numker
	          isel=0
	          lstr=lnblnk(varstr(ivarkern(iker)))
	          vstr=(varstr(ivarkern(iker)))
	          if(ieval.eq.1) then
	            if(vstr(1:lstr).eq.'UM (SH+SV)*0.5,'.or.
     #		       vstr(1:lstr).eq.'LM (SH+SV)*0.5,'.or.
     #		       vstr(1:lstr).eq.'EA (SH+SV)*0.5,') then
	              isel=1
		    endif
	          else if(ieval.eq.2) then
	            if(vstr(1:lstr).eq.'UM SH-SV,'.or.
     #	       	       vstr(1:lstr).eq.'LM SH-SV,'.or.
     #	       	       vstr(1:lstr).eq.'EA SH-SV,') then
	              isel=1
	            endif
	          endif

	          if(isel.eq.1) then
	            if(vercof(iker).ne.0.) then 
                      if(itypehpa(ihpakern(iker)).eq.1) then
		        ihpa=ihpakern(iker)
                        nylm=(lmxhpa(ihpakern(iker))+1)**2
                        do i=1,nylm
                          value=value+vercof(iker)*ylmcof(i,ihpa)
     #                                *coe(i,iker)
                        enddo
                      else if(itypehpa(ihpakern(iker)).eq.2) then
		        ihpa=ihpakern(iker)
                        do i=1,nconpt(ihpa)
                          iver=iconpt(i,ihpa)
                          value=value+vercof(iker)*conpt(i,ihpa)
     #                                *coe(iver,iker)
                        enddo
                      else
                        write(6,"('problem 2')")
                        stop
                      endif ! --- itypehpa
	            endif ! --- vercof(iker).ne.0.
	          endif ! --- isel.eq.1
	        enddo ! --- end of do iker=1,numker
	   
	        valu(ieval)=value
	      enddo ! --- ieval

c 	      --- evaluate perturbations in vsh and vsv 

              if(ish.eq.1) then 
                vsh3drel=valu(1)+0.5*valu(2) 	! --- dvSH
              else if(ish.eq.0) then
                vsv3drel=valu(1)-0.5*valu(2) 	! --- dvSV
              else
                stop 'something wrong'
              endif

            enddo ! --- by isht

c --- save absolute velocities in 1D and 3D models

	    vsh3d(ipx)=(vsh3drel+1.0)*vsh1d
	    vsv3d(ipx)=(vsv3drel+1.0)*vsv1d

	enddo ! --- loop on ipx

c --- calculate the Voigt average

	vpterm=(vpv1d*vpv1d+vph1d*vph1d*(1.0-2.0*eta1d))
	vg1d=5.0*vsh1d*vsh1d+(6.0+4.0*eta1d)*vsv1d*vsv1d
	vg1d=sqrt((vpterm+vg1d)/15.0)
	do ipx=1,npx
	  vg3d(ipx)=5.0*vsh3d(ipx)*vsh3d(ipx)+
     #              (6.0+4.0*eta1d)*vsv3d(ipx)*vsv3d(ipx)
	  vg3d(ipx)=sqrt((vpterm+vg3d(ipx))/15.0)
	enddo ! --- loop on ipx

c --- calculate the anisotropy

	do ipx=1,npx
	  aniso(ipx)=vsh3d(ipx)-vsv3d(ipx)
	enddo

c --- convert into relative perturbations, if necessary  

	if(irel.eq.1) then
	  do ipx=1,npx
	    vg3d(ipx)=100.0*(vg3d(ipx)-vg1d)/vg1d
	    aniso(ipx)=100.0*aniso(ipx)/vg1d
	    vsh3d(ipx)=100.0*(vsh3d(ipx)-vsh1d)/vsh1d
	    vsv3d(ipx)=100.0*(vsv3d(ipx)-vsv1d)/vsv1d
	  enddo
	endif

c --- calculate and remove the average if necessary

	if(ifremav.eq.1) then
	  avervg=0.0	
	  averani=0.0	
	  avervsh=0.0	
	  avervsv=0.0	
	  do ipx=1,npx
	    avervg=avervg+vg3d(ipx)*areaarr(ipx) 
	    averani=averani+aniso(ipx)*areaarr(ipx) 
	    avervsh=avervsh+vsh3d(ipx)*areaarr(ipx) 
	    avervsv=avervsv+vsv3d(ipx)*areaarr(ipx) 
	  enddo ! --- loop on ipx
	  avervg=avervg/totarea
	  averani=averani/totarea
	  avervsh=avervsh/totarea
	  avervsv=avervsv/totarea
	  do ipx=1,npx
	    vg3d(ipx)=vg3d(ipx)-avervg
	    aniso(ipx)=aniso(ipx)-averani
	    vsh3d(ipx)=vsh3d(ipx)-avervsh
	    vsv3d(ipx)=vsv3d(ipx)-avervsv
	  enddo
	endif 

c --- adjust xlon 

	if(xmaxlo.gt.180.0) then
          do ipx=1,npx
      	    if(xlonarr(ipx).lt.0.0) xlonarr(ipx)=xlonarr(ipx)+360.0
          enddo
	endif
	
c --- write out the Voigt average  

	open(luout,file='voigt.txt')
	do ipx=1,npx
	  write(luout,"(2f7.2,f9.4)")xlonarr(ipx),xlatarr(ipx),vg3d(ipx)	
	enddo
	close(luout)

c --- write out vsh  

	open(luout,file='vsh.txt')
	do ipx=1,npx
	  write(luout,"(2f7.2,f9.4)")xlonarr(ipx),xlatarr(ipx),vsh3d(ipx)	
	enddo
	close(luout)

c --- write out vsv

	open(luout,file='vsv.txt')
	do ipx=1,npx
	  write(luout,"(2f7.2,f9.4)")xlonarr(ipx),xlatarr(ipx),vsv3d(ipx)	
	enddo
	close(luout)

c --- write out anisotropy 

	open(luout,file='ani.txt')
	do ipx=1,npx
	  write(luout,"(2f7.2,f9.4)")xlonarr(ipx),xlatarr(ipx),aniso(ipx)	
	enddo
	close(luout)
	end


	


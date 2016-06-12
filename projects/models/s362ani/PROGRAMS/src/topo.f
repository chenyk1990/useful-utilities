c --- evaluate topography (depression) in km

	implicit none

	real twopi
	integer maxker,maxl,maxcoe,maxver,maxhpa,mxpx
	parameter (maxker=100)		! --- radial kernels
	parameter (maxl=72)		! --- sph harmonics
	parameter (maxcoe=1500)		! --- sph splines
	parameter (maxver=1000)		! --- contibuting splines

	parameter (maxhpa=2)		
c	parameter (mxpx=64800)		! --- max pts on the map
	parameter (mxpx=266000)		! --- max pts on the map
        parameter (twopi=6.2831853)

c --- array of pixels

	integer ipx,npx
	real*4 xlatarr(mxpx)
	real*4 xlonarr(mxpx)
	real*4 areaarr(mxpx)
	real*4 xvalarr410(mxpx)
	real*4 xvalarr650(mxpx)

c --- radial kernels

	character*80 kerstr
	integer iker,numker,ierror

c --- horizontal kernels

	integer numcof,numhpa,ihpa,i
	integer itypehpa(maxhpa)
	integer ihpakern(maxker)
	integer numcoe(maxhpa)
	real*4 xlat,xlon
	real*4 xminla,xmaxla,xminlo,xmaxlo

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

c --- area used to find an average

	real*4 aver,area,totarea,arealatzone
	integer ifremav
	character*64 avstr

c --- other

	integer filmap
	real*4 size,value410,value650,x,y
	
	character*128 string
	character*128 model
	character*128 modeldef
	data modeldef/' '/

	logical exists
	logical gotmodel
c
	integer itpspl(maxcoe,maxhpa)
	character*80 refmdl
	character*80 hsplfl(maxhpa)
	character*40 dskker(maxker)
c
	character*40 varstr(maxker)
	integer ivarkern(maxker)
c
	integer lu,imap,numvar
	integer lnblnk  
        real*8 xrad
c
c -------------------------------------

      xrad=3.14159265358979/180.d0

c --- read the model file

      lu=19 !--- model 
      imap=1

      gotmodel=.false.
      modeldef=''
      do while(.not.gotmodel)
	if(lnblnk(modeldef).eq.0) call getcurrmodels(modeldef,imap)
        write(6,"('change earth model [',a,']:')") 
     #		modeldef(1:lnblnk(modeldef))
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
   	    write(6,"(i2,' variables')") numvar
              do i=1,numvar
                write(6,"(i2,1x,a)") i,varstr(i)
              enddo
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

      if(numker.eq.0) stop 'numker.eq.0'
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

c -- remove average ?

      write(6,"('remove average? (y or n)? ')")
      read(5,*) avstr
      if(avstr(1:1).eq.'y'.or.avstr(1:3).eq.
     #		'yes'.or.avstr(1:1).eq.'Y') then
        ifremav=1
      else if(avstr(1:1).eq.'n'.or.avstr(1:2).eq.
     #		'no'.or.avstr(1:1).eq.'N') then
        ifremav=0
      else
        write(6,"('do not know if average should be removed, error')")
	pause
      endif

c --- specify block size and the region of interest

      write(6,"('block size: ')")
      read(5,*)size
c      write(6,"('region of interest: xminla,xmaxla,xminlo,xmaxlo ')")
c      read(5,*) xminla,xmaxla,xminlo,xmaxlo
      xminla=-90.0
      xmaxla=90.0
      xminlo=-180.0
      xmaxlo=180.0

c --- find pixels

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
c      write(6,"('number of pixels= ',*)")npx

c --- find areas used to calculate the spherical average	

      totarea=0.0	 
      do ipx=1,npx
        xlat=xlatarr(ipx)
	xlon=xlonarr(ipx)
        arealatzone=sin(xrad*(xlat+size/2.0))-sin(xrad*(xlat-size/2.0))
        area=twopi*arealatzone*size/360.0
        areaarr(ipx)=area
        totarea=totarea+area
      enddo

c --- loop by pixels and evaluate the model 

      do ipx=1,npx
        xlat=xlatarr(ipx)
        xlon=xlonarr(ipx)

c --- evaluate contributing horizontal basis functions 

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
     #		  nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
          else
            write(6,"('problem 1')") 
          endif
        enddo

c --- evaluate the model at lat and lon

        value410=0.
        value650=0.
        do iker=1,numker
          if(varstr(ivarkern(iker)).eq.'Topo 400,') then
              if(itypehpa(ihpakern(iker)).eq.1) then
                ihpa=ihpakern(iker)
                nylm=(lmxhpa(ihpakern(iker))+1)**2
                do i=1,nylm
                  value410=value410+ylmcof(i,ihpa)*coe(i,iker)
                enddo
              else if(itypehpa(ihpakern(iker)).eq.2) then
                ihpa=ihpakern(iker)
                do i=1,nconpt(ihpa)
                  iver=iconpt(i,ihpa)
                  value410=value410+conpt(i,ihpa)*coe(iver,iker)
                enddo
              else
                write(6,"('problem 2')")
                stop
              endif
          else if(varstr(ivarkern(iker)).eq.'Topo 670,') then
              if(itypehpa(ihpakern(iker)).eq.1) then
                ihpa=ihpakern(iker)
                nylm=(lmxhpa(ihpakern(iker))+1)**2
                do i=1,nylm
                  value650=value650+ylmcof(i,ihpa)*coe(i,iker)
                enddo
              else if(itypehpa(ihpakern(iker)).eq.2) then
                ihpa=ihpakern(iker)
                do i=1,nconpt(ihpa)
                  iver=iconpt(i,ihpa)
                  value650=value650+conpt(i,ihpa)*coe(iver,iker)
                enddo
              else
                write(6,"('problem 2')")
                stop
              endif
          endif ! --- ivarsel(ivarkern(iker)).eq.1
        enddo ! --- end of do iker=1,numker

	xvalarr410(ipx)=value410
	xvalarr650(ipx)=value650
      enddo	! --- end of loop by ipx

c --- remove the average

      if(ifremav.eq.1) then
        aver=0.0
        do ipx=1,npx
      	  aver=aver+xvalarr410(ipx)*areaarr(ipx) 
        enddo
        aver=aver/totarea
	do ipx=1,npx
	  xvalarr410(ipx)=xvalarr410(ipx)-aver	  
	enddo
        aver=0.0
        do ipx=1,npx
      	  aver=aver+xvalarr650(ipx)*areaarr(ipx) 
        enddo
        aver=aver/totarea
	do ipx=1,npx
	  xvalarr650(ipx)=xvalarr650(ipx)-aver	  
	enddo
      endif

c --- write out everything

      filmap=19
      open(filmap,file='topo410.txt')
      do ipx=1,npx
	write(filmap,"(2f9.3,f11.4)")xlonarr(ipx),
     #xlatarr(ipx),xvalarr410(ipx)
      enddo
      close(filmap)
      open(filmap,file='topo650.txt')
      do ipx=1,npx
	write(filmap,"(2f9.3,f11.4)")xlonarr(ipx),
     #xlatarr(ipx),xvalarr650(ipx)
      enddo
      close(filmap)
      end




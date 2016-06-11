      subroutine splcon(xlat,xlon,nver,verlat,verlon,verrad,
     #                  ncon,icon,con)
      dimension verlat(1)
      dimension verlon(1)
      dimension verrad(1)
      dimension icon(1)
      dimension con(1)
c
c---- all calculations done in double precision
c---- but results are single
c
      real*8 dd
      real*8 rn
      real*8 dr
c
      ncon=0
      do iver=1,nver
        if(xlat.gt.verlat(iver)-2.*verrad(iver)) then
          if(xlat.lt.verlat(iver)+2.*verrad(iver)) then
            dd=dsind(dble(verlat(iver)))*dsind(dble(xlat))
            dd=dd+dcosd(dble(verlat(iver)))*dcosd(dble(xlat))*
     #         dcosd(dble(xlon-verlon(iver)))
            dd=dacosd(dd)
            if(dd.gt.dble(verrad(iver))*2.d0) then
            else
              ncon=ncon+1
              icon(ncon)=iver
              rn=dd/dble(verrad(iver))
              dr=rn-1.d0
              if(rn.le.1.d0) then
                con(ncon)=(0.75d0*rn-1.5d0)*(rn**2)+1.d0
              else if(rn.gt.1.d0) then
                con(ncon)=((-0.25d0*dr+0.75d0)*dr-0.75d0)*dr+0.25d0
              else
                con(ncon)=0.
              endif
            endif
          endif
        endif
      enddo
      return
      end

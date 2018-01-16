!     ------------------------------------------------------------------
!
!     CapCalc5.f
!
!
!     Created on Sep 25 13:19:10 2008
!     Copyright (c) 2008 R.Levicky. All rights reserved.
!
!
!     ------------------------------------------------------------------
      program capcalc5

C units: program uses kg, m, C, s, K, V     
C version 3 implements finite thickness of the Morpholino layer
C still, only MCP, Morpholino, and double layer in solution are 
C included
C MCP is assumed to be a perfect dielectric
C cdna is concentration of DNA charge 'in the Morpholino layer'
      
      implicit none
      
      integer lrs
      character*1 ch
      real*8 kb,eo,e,na

      common/consts/kb,eo,e,na,lrs

      kb=1.38065d-23;eo=8.85419d-12;e=1.60218d-19;na=6.02d23
      print*,'number of layers (1 to 10)'
      read*,lrs

      call inputs

      print*
      print*,'1). Calculate a V(x) curve at a single V'
      print*,'2). Calculate a C(V) curve over entire V range'
      read*,ch

      if (ch.eq.'1') then
        call getv
      elseif (ch.eq.'2') then
        call getcap
      endif

      pause 'Results saved'

      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine inputs
      
      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      integer i
      character*40 infile
      character*10 dummy
      
      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol

      data lt/70*0.d0/
      
      print*,'inputs file'
      read*,infile
      open(unit=11,file=infile,status='old')
C      pause 'opened file'
      read(11,30) dummy,v1l
      read(11,30) dummy,v1h
      read(11,30) dummy,vstep
      read(11,30) dummy,t
      read(11,30) dummy,cptot
      read(11,30) dummy,ka1
      read(11,30) dummy,ka2
      read(11,30) dummy,ka3
      read(11,30) dummy,ph
      read(11,30) dummy,esol
30    format(a10,3x,d23.16)
      do 1 i=1,lrs
        read(11,3) (lt(i,j),j=1,7)
1     enddo
3     format(d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16,
     &3x,d23.16)
      close(11)
      
C col 1 of lt is dielectric constant (dimensionless)
C col 2 of lt is layer thickness in m
C col 3 of lt is beta na (dimensionless)
C col 4 of lt is beta pm
C col 5 of lt is beta pd
C col 6 of lt is beta pt
C col 7 of lt is charge concentration in mol/m^3

      cpm=cptot/(10.d0**(ka1-ph)+1.d0+10.d0**(ph-ka2)+10.d0**
     &(2.d0*ph-ka3-ka2))
      cpd=cpm*10.d0**(ph-ka2)
      cpt=cpm*10.d0**(2.d0*ph-ka2-ka3)
      cp=cpm*10.d0**(ka1-ph)
      cna=cpm+2.d0*cpd+3.0d0*cpt
     
      print*,'V1L=',v1l,'  V1H=',v1h,'  Vstep=',vstep
      print*,'t=',t,'  cptot=',cp+cpm+cpd+cpt,'  ka1=',ka1
      print*,'ka2=',ka2,'  ka3=',ka3,'  ph=',ph
      print*,'cna=',cna,'  cpm=',cpm,'  cpd=',cpd
      print*,'cpt=',cpt,'  cp=',cp,'  esol=',esol
      do 5 i=1,lrs
         print*,i,' ',lt(i,1),' ',lt(i,2),' ',lt(i,3),' ',lt(i,4),
     &' ',lt(i,5),' ',lt(i,6),' ',lt(i,7)
5     enddo
      print*

      points=idint((v1h-v1l)/vstep + 2)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getcap
      
      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      real*8 v1,slope1,x,stepx,xend,y(2)
      real*8 param(50)
      integer pos,ido,lastpt,lpos
      real*8 results(1:points,1:lrs+4),width,dummy1,dummy2

      external eqnslr

      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol
      common/ldesig/lpos,width

      data param/50*0.d0/

      pos=1
      v1=v1l
      width=0.d0
      do 16 i=1,lrs
        width=width+lt(i,2)
16    enddo        

8     if (v1.gt.0.d0) then
        slope1=-gradsqrdl(v1)**0.5d0
      elseif (v1.le.0.d0) then
        slope1=gradsqrdl(v1)**0.5d0
      endif
      lpos=lrs
      ido=1
      param(3)=0.1
      param(4)=1.d6
      param(10)=1.0

      x=1.d0
      stepx=-x/200.d0
      results(pos,lrs+4)=v1
      
18    if (lpos.eq.lrs) then
        y(1)=esol/lt(lpos,1)*slope1*e/(kb*t)*width
        y(2)=v1*e/(kb*t)
      else
        y(1)=lt(lpos+1,1)/lt(lpos,1)*y(1)
      endif
        
      xend=x+stepx
      xswitch=x-lt(lpos,2)/width

4     if (xend.le.xswitch) goto 5
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      xend=xend+stepx
      goto 4
 
5     xend=xswitch
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      
C     if (lpos.eq.2) dummy1=(y(1)*kb*t/(e*width))**2.d0
      
      if (lpos.eq.1) then
        results(pos,1)=y(2)*kb*t/e
      else
        results(pos,lpos+3)=y(2)*kb*t/e
      endif

      lpos=lpos-1
      if (lpos.ge.1) then
        x=xswitch
        goto 18
      endif
     
C results array has following columns: 1. electrode potential, 2. capacitance, 3. surface charge, 4. electric field at surface
C 5 to lrs+4: potentials at interfaces between layers (i.e. column 5 is potential at layer 1/2 interface, column lrs+4 is 
C potential at lrs/electrolyte interface

      results(pos,4)=y(1)*kb*t/(e*width)
      results(pos,3)=-y(1)*lt(1,1)*eo*kb*t/(e*width)

      ido=3
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      
C check accuracy of RK numerical integration

C      dummy2=gradsqrmorph(results(pos,5),results(pos,6))
C      print*,'dv/dx^2|layer1/2 from RK slope=',dummy1
C      print*,'dv/dx^2|layer1/2 from once integrated expression=',dummy2
C      print*,'% difference=', dabs((dummy1-dummy2)/dummy2)*100.d0
C      pause

      v1=v1+vstep
      pos=pos+1
      if (v1.le.v1h) goto 8
      
C calculate capacitance from the derivative dq/dv(x=0); this misses the first and last points to avoid extrapolation

      lastpt=pos-1
      pos=2
12    results(pos,2)=0.5d0*((results(pos+1,3)-results(pos,3))/
     &(results(pos+1,1)-results(pos,1))+
     &(results(pos,3)-results(pos-1,3))/
     &(results(pos,1)-results(pos-1,1)))
      pos=pos+1
      if (pos.lt.lastpt) goto 12

C save results array and other parameters - exclude 1st and last points where cap was not calculated

      open(unit=11,file='capresults.txt',status='unknown')
      do 20 i=2,lastpt-1
        write(11,10) (results(i,j),j=1,lrs+3)
        write(11,22) results(i,lrs+4)
20    enddo
10    format(d23.16,',',\)
22    format(d23.16)
      close(11)

      open(unit=11,file='parameters.txt',status='unknown')
      write(11,14) 'v1l       ',v1l
      write(11,14) 'v1h       ',v1h
      write(11,14) 'vstep     ',vstep
      write(11,14) 'temperatur',t
      write(11,14) 'phosph tot',cptot
      write(11,14) 'ka1       ',ka1
      write(11,14) 'ka2       ',ka2
      write(11,14) 'ka3       ',ka3
      write(11,14) 'pH        ',ph
      write(11,14) 'esol      ',esol
14    format(a10,3x,d23.16)
      do 34 i=1,lrs
        write(11,32) (lt(i,j),j=1,7)
34    enddo
32    format(d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x
     &,d23.16)

      close(11)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getv
      
      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      real*8 v1,slope1,x,stepx,xend,y(2)
      real*8 param(50)
      integer ido,lpos
      real*8 width,varr(200,6)

      external eqnslr,eqnslre

      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol
      common/ldesig/lpos,width

      data param/50*0.d0/
      data varr/1200*0.d0/

      print*
      print*,'potential (V)'
      read*,v1
      print*

      open(unit=11,file='v.txt',status='unknown')

      if (v1.gt.0.d0) then
        slope1=-gradsqrdl(v1)**0.5d0
      elseif (v1.le.0.d0) then
        slope1=gradsqrdl(v1)**0.5d0
      endif

      width=0.d0
      do 16 i=1,lrs
        width=width+lt(i,2)
16    enddo        

      ido=1
      param(3)=0.1
      param(4)=1.d6
      param(10)=1.0

      x=1.d0
      stepx=x/200.d0
      xend=x+stepx
      
      varr(1,1)=x*width
      varr(1,2)=v1
      varr(1,3)=cna*dexp(-e*varr(1,2)/(kb*t))
      varr(1,4)=cpm*dexp(e*varr(1,2)/(kb*t))
      varr(1,5)=cpd*dexp(2.d0*e*varr(1,2)/(kb*t))
      varr(1,6)=cpt*dexp(3.d0*e*varr(1,2)/(kb*t))
      
      y(1)=slope1*e/(kb*t)*width
      y(2)=v1*e/(kb*t)

      do 42 i=2,200
        call divprk(ido,2,eqnslre,x,xend,1.d-5,param,y)
        varr(i,1)=xend*width
        varr(i,2)=y(2)*kb*t/e
        varr(i,3)=cna*dexp(-e*varr(i,2)/(kb*t))
        varr(i,4)=cpm*dexp(e*varr(i,2)/(kb*t))
        varr(i,5)=cpd*dexp(2.d0*e*varr(i,2)/(kb*t))
        varr(i,6)=cpt*dexp(3.d0*e*varr(i,2)/(kb*t))
        xend=xend+stepx
42    enddo

      ido=3
      call divprk(ido,2,eqnslre,x,xend,1.d-5,param,y)

      do 44 i=200,1,-1
        write(11,40) (varr(i,j),j=1,6)
44    enddo

      if (v1.gt.0.d0) then
        slope1=-gradsqrdl(v1)**0.5d0
      elseif (v1.le.0.d0) then
        slope1=gradsqrdl(v1)**0.5d0
      endif
      lpos=lrs
      ido=1
      param(3)=0.1
      param(4)=1.d6
      param(10)=1.0

      x=1.d0
      stepx=-x/200.d0

C      write(11,40) width,v1
40    format(d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16,3x,d23.16)
      
18    if (lpos.eq.lrs) then
        y(1)=esol/lt(lpos,1)*slope1*e/(kb*t)*width
        y(2)=v1*e/(kb*t)
      else
        y(1)=lt(lpos+1,1)/lt(lpos,1)*y(1)
      endif
        
      xend=x+stepx
      xswitch=x-lt(lpos,2)/width

4     if (xend.le.xswitch) goto 5
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      write(11,40) xend*width,y(2)*kb*t/e,cna*dexp(-y(2)-lt(lpos,3)),
     &cpm*dexp(y(2)-lt(lpos,4)),cpd*dexp(2.d0*y(2)-lt(lpos,5)),
     &cpt*dexp(3.d0*y(2)-lt(lpos,6))
      xend=xend+stepx
      goto 4
 
5     xend=xswitch
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      write(11,40) xend*width,y(2)*kb*t/e,cna*dexp(-y(2)-lt(lpos,3)),
     &cpm*dexp(y(2)-lt(lpos,4)),cpd*dexp(2.d0*y(2)-lt(lpos,5)),
     &cpt*dexp(3.d0*y(2)-lt(lpos,6))
            
      lpos=lpos-1
      if (lpos.ge.1) then
        x=xswitch
        goto 18
      endif
     
      ido=3
      call divprk(ido,2,eqnslr,x,xend,1.d-5,param,y)
      
      close(11)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine eqnslre(n,x,y,yprime)

      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      integer n,lpos
      real*8 x,y(n),yprime(n),width

      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol
      common/ldesig/lpos,width

      rho=cna*dexp(-y(2))
      rho=rho-cpm*dexp(y(2))
      rho=rho-2.d0*cpd*dexp(2.d0*y(2))
      rho=rho-3.d0*cpt*dexp(3.d0*y(2))
      rho=na*e*rho
      yprime(1)=-rho*e*width**2.d0/(esol*eo*kb*t)
      yprime(2)=y(1)

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine eqnslr(n,x,y,yprime)

      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      integer n,lpos
      real*8 x,y(n),yprime(n),width

      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol
      common/ldesig/lpos,width

      rho=cna*dexp(-y(2)-lt(lpos,3))
      rho=rho-cpm*dexp(y(2)-lt(lpos,4))
      rho=rho-2.d0*cpd*dexp(2.d0*y(2)-lt(lpos,5))
      rho=rho-3.d0*cpt*dexp(3.d0*y(2)-lt(lpos,6))
      rho=rho+lt(lpos,7)
      rho=na*e*rho
      yprime(1)=-rho*e*width**2.d0/(lt(lpos,1)*eo*kb*t)
      yprime(2)=y(1)

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function gradsqrdl(v1)

      real*8 lt(1:10,1:7)
      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
      real*8 cpm,cpd,cpt,cp,cna,esol
      integer lrs,points
      real*8 kb,eo,e,na
      real*8 v1
      real*8 dummy

      common/consts/kb,eo,e,na,lrs
      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
     &cpt,cp,lt,points,esol

C double layer is in solution, so there are no partitioning beta coeffs or immobilized charge

      dummy=cna*dexp(-e*v1/(kb*t))*
     &(dexp(e*v1/(kb*t))-1.d0)
      dummy=dummy+cpm*dexp(e*v1/(kb*t))*
     &(dexp(-e*v1/(kb*t))-1.d0)
      dummy=dummy+cpd*dexp(2.d0*e*v1/(kb*t))*
     &(dexp(-2.d0*e*v1/(kb*t))-1.d0)
      dummy=dummy+cpt*dexp(3.d0*e*v1/(kb*t))*
     &(dexp(-3.d0*e*v1/(kb*t))-1.d0)
      dummy=(kb*t/e)*dummy
      dummy=-2.d0*na*e/(esol*eo)*dummy
      
      gradsqrdl=dummy
            
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      real*8 function gradsqrmorph(v2,v1)

C      real*8 lt(1:10,1:7)
C      real*8 v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph
C      real*8 cpm,cpd,cpt,cp,cna,esol
C      integer lrs,points
C      real*8 kb,eo,e,na
C      real*8 v1,v2
C      real*8 part1,part2

C      common/consts/kb,eo,e,na,lrs
C      common/params/v1l,v1h,vstep,t,cptot,ka1,ka2,ka3,ph,cna,cpm,cpd,
C     &cpt,cp,lt,points,esol

C      part1=(esol**2.d0/emorph**2.d0)*gradsqrdl(v1)

C      part2=cna*dexp(-betana)*dexp(-e*v2/(kb*t))*
C     &(dexp(e*(v2-v1)/(kb*t))-1.d0)
C      part2=part2+cpm*dexp(-betapm)*dexp(e*v2/(kb*t))*
C     &(dexp(-e*(v2-v1)/(kb*t))-1.d0)
C      part2=part2+cpd*dexp(-betapd)*dexp(2.d0*e*v2/(kb*t))*
C     &(dexp(-2.d0*e*(v2-v1)/(kb*t))-1.d0)
C      part2=part2+cpt*dexp(-betapt)*dexp(3.d0*e*v2/(kb*t))*
C     &(dexp(-3.d0*e*(v2-v1)/(kb*t))-1.d0)
        
C      part2=(kb*t/e)*part2+qprob*(v2-v1)
C      part2=-2.d0*na*e/(emorph*eo)*part2
      
C      gradsqrmorph=part1+part2
            
C      return
C      end
      

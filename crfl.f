c     crfl.f
c     this program was rewritten from the original reflectivity
c     version (fuchs and mueller) and vectorized by
c     k.-j. sandmeier
c     institute of geophysics
c     hertzstr. 16
c     d-7500 karlsruhe, f.r.g.
c
c     program p-sv und sh
c
c     maximum number of frequencies:    nfr=1024            see
c     maximum number of ray-parameters: npa=2500         parameter
c     maximum number of layers:         nla=1000         statement
c     maximum number of distancies:     ndis=200          statement
c
c     Altered to these values: NS 1/2/06
c     maximum number of frequencies:    nfr=8192            see
c     maximum number of ray-parameters: npa=5000         parameter
c     maximum number of layers:         nla=1000         statement
c     maximum number of distancies:     ndis=200          statement
c
c     Altered to these values: NS 1/31/06
c     maximum number of frequencies:    nfr=131072            see
c     maximum number of ray-parameters: npa=5000         parameter
c     maximum number of layers:         nla=1000         statement
c     maximum number of distancies:     ndis=200          statement
c
c     Moment tensor read in degrees
c
c     maximum number of sources = 100  for double couple
c
c     program main(input,output,a2,a3,a21,a22,a23,
c    &unit5=input,unit6=output,unit2=a2,unit3=a3,
c    &unit21=a21,unit22=a22,unit23=a23)
c      parameter(nfr=1024,npa=2500,nla=1000,ndis=200)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000,ndis=200)
      common/bsh1/xs(100),ys(100),es(100),rr(100),dz(100),vv(100),
     &ts(100),zs(100),az(1000),alpr(nla),betr(nla),ztrue(100)
      common/blo1/z(0:nla),d(nla),alp(nla),alpi(nla),bet(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,
     &nrho,dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),
     &phasc(nfr),sphasc(nfr),cphasc(nfr)
      common x1r(nfr),x1i(nfr),x2r(nfr),x2i(nfr),x3r(nfr),x3i(nfr),
     &x4r(nfr),x4i(nfr),x5r(nfr),x5i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),r1r(nfr),r1i(nfr),r2r(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr),xsum1r(nfr,ndis),
     &xsum1i(nfr,ndis),xsum2r(nfr,ndis),xsum2i(nfr,ndis)
      common/integ1/tshr(npa),tshi(npa),rshd1r(npa),rshd1i(npa),
     &rshd2r(npa),rshd2i(npa),rshu1r(npa),rshu1i(npa),rshu2r(npa),
     &rshu2i(npa)
      common/integ2/xx1r(npa),xx1i(npa),yy1r(npa),yy1i(npa),xx2r(npa),
     &xx2i(npa),xx6i(npa)
      common/integ3/ersin(npa),ersin1(npa),xas(npa),xj0r(npa),xj0i(npa),
     &xj1r(npa),xj1i(npa)
      common/bsh2/omega2(nfr),phase(0:nfr),
c     &s(8193),wavel(4096),balla(8192),iballa(8192),window(npa),
c     &phase1(nfr),sc1(4096)
c     &s(32769),wavel(16384),balla(32768),iballa(32768),window(npa),
c     &phase1(nfr),sc1(16384)
c     &s(65537),wavel(32768),balla(65536),iballa(65536),window(npa),
c     &phase1(nfr),sc1(32768)
     &s(131073),wavel(65536),balla(131072),iballa(131072),window(npa),
     &phase1(nfr),sc1(65536)
      common/bsh3/cwil,cwir,c1,c2,dnue,dp,ein,einh,fo,fu,fwil,fwir,
     &inpt,io,iu,kom,na,nent,nh,tsigma
      common/bsh4/kappa1,kappa2,kappa3,lamda1,lamda2,m11,m12,m13,m22,
     &m23,m33,lamda
      complex ein,einh,wavel,s,am,bm,al(nla),be(nla),bmiso,bmiso2,sc1,
     &cv
      real kappa1,kappa2,kappa3,lamda1,lamda2,m11,m12,m13,m22,m23,m33,
     &lamda
      dimension rppd1r(nfr),rppd1i(nfr),rspd1r(nfr),rspd1i(nfr),
     &rpsd1r(nfr),rpsd1i(nfr),rssd1r(nfr),rssd1i(nfr),
     &rppu1r(nfr),rppu1i(nfr),rspu1r(nfr),rspu1i(nfr),rpsu1r(nfr),
     &rpsu1i(nfr),rssu1r(nfr),rssu1i(nfr),
     &rppd2r(nfr),rppd2i(nfr),rspd2r(nfr),rspd2i(nfr),rpsd2r(nfr),
     &rpsd2i(nfr),rssd2r(nfr),rssd2i(nfr),
     &rppu2r(nfr),rppu2i(nfr),rspu2r(nfr),rspu2i(nfr),rpsu2r(nfr),
     &rpsu2i(nfr),rssu2r(nfr),rssu2i(nfr),
     &tppu1r(nfr),tppu1i(nfr),tspu1r(nfr),tspu1i(nfr),tpsu1r(nfr),
     &tpsu1i(nfr),tssu1r(nfr),tssu1i(nfr),
     &tppd2r(nfr),tppd2i(nfr),tspd2r(nfr),tspd2i(nfr),tpsd2r(nfr),
     &tpsd2i(nfr),tssd2r(nfr),tssd2i(nfr),
     &tppu2r(nfr),tppu2i(nfr),tspu2r(nfr),tspu2i(nfr),tpsu2r(nfr),
     &tpsu2i(nfr),tssu2r(nfr),tssu2i(nfr),
     &a11r(nfr),a11i(nfr),a12r(nfr),a12i(nfr),a21r(nfr),a21i(nfr),
     &a22r(nfr),a22i(nfr),b11r(nfr),b11i(nfr),b12r(nfr),b12i(nfr),
     &b21r(nfr),b21i(nfr),b22r(nfr),b22i(nfr),c11r(nfr),c11i(nfr),
     &c12r(nfr),c12i(nfr),c21r(nfr),c21i(nfr),c22r(nfr),c22i(nfr),
     &d11r(nfr),d11i(nfr),d12r(nfr),d12i(nfr),d21r(nfr),d21i(nfr),
     &d22r(nfr),d22i(nfr)
      dimension s11dr(nfr),s11di(nfr),s12dr(nfr),s12di(nfr),
     &s21dr(nfr),s21di(nfr),s22dr(nfr),s22di(nfr),
     &s11ur(nfr),s11ui(nfr),s12ur(nfr),s12ui(nfr),s21ur(nfr),s21ui(nfr),
     &s22ur(nfr),s22ui(nfr),
     &ua1r(nfr),ua1i(nfr),ua2r(nfr),ua2i(nfr),ua3r(nfr),ua3i(nfr),
     &ua4r(nfr),ua4i(nfr),
     &ub1r(nfr),ub1i(nfr),ub2r(nfr),ub2i(nfr),ub3r(nfr),ub3i(nfr),
     &ub4r(nfr),ub4i(nfr),
     &r11r(nfr),r11i(nfr),r12r(nfr),r12i(nfr),r21r(nfr),r21i(nfr),
     &r22r(nfr),r22i(nfr),
     &x11r(nfr),x11i(nfr),x12r(nfr),x12i(nfr),x21r(nfr),x21i(nfr),
     &x22r(nfr),x22i(nfr),
     &y11r(nfr),y11i(nfr),y12r(nfr),y12i(nfr),y21r(nfr),y21i(nfr),
     &y22r(nfr),y22i(nfr)
      character*80 dentif
      character*7 filename

c     call link("unit5=(crfl.dat,open,text),unit6=(crfl.out,create,
c    &           text),unit2=(crfl.psv,create,text),
c    &           unit3=(crfl.sh,create,text),
c    &           unit21=(a21,create,text),unit22=(a22,create,text),
c    &           unit23=(a23,create,text)//")
      open(5,file='crfl.dat',access='sequential',form='formatted')
      open(6,file='crfl.out',access='sequential',form='formatted')
cc    open(2,file='crfl.psv',access='sequential',form='formatted')
      rewind 5
      rewind 6
cc    rewind 2
      z(0)=0.
      iouts=2
      ioutsh=3
      iins=5
      pi=3.14159265
      pi2=2.*pi
      pi4=pi/4.
      ein=cmplx(0.,1.)
      einh=cmplx(0.,0.5)
      read(iins,501) dentif
 501  format(a80)
      write(6,501) dentif
      read(iins,502) (iss(i),i=1,25)
 502  format(5(5i2,2x))
      write(6,699) (iss(i),i=1,25)
 699  format(5(5i2,2x)/)
c     bedeutung der schalter
      if(iss(1).eq.0) write(6,671)
 671  format(' the full wavefield is calculated')
      if(iss(1).eq.1) write(6,672)
 672  format(' only pp-reflections are considered')
      if(iss(1).eq.2) write(6,673)
 673  format(' only ss-reflections are considered')
      if(iss(1).eq.3) write(6,674)
 674  format(' pp- and ss-reflections are considered')
      if(iss(2).eq.0) write(6,675)
 675  format(' for transmitting layers - conversion is considered')
      if(iss(2).eq.1) write(6,676)
 676  format(' for transmitting layers - no conversion (tps=tsp=0)')
      if(iss(3).eq.0) write(6,677)
 677  format(' seismogram for displacement')
      if(iss(3).eq.1) write(6,678)
 678  format(' seismogram for velocity')
      if(iss(3).eq.2) write(6,679)
 679  format(' seismogram for acceleration')
cc    if(iss(5).eq.0) then
cc    open(3,file='crfl.sh',access='sequential',form='formatted')
cc    rewind 3
cc    endif
      if(iss(6).eq.1) then
      open(21,file='s21',access='sequential',form='formatted')
      open(22,file='s22',access='sequential',form='formatted')
      open(23,file='s23',access='sequential',form='formatted')
      rewind 21
      rewind 22
      rewind 23
      endif
      if(iss(7).eq.1) write(6,680)
 680  format(' no upgoing waves from the source are considered')
c     if(iss(13).eq.0) write(6,681)
      write(6,681)
 681  format(' bessel-function: only far-field approximation')
c     if(iss(13).eq.1) write(6,682)
c682  format(' bessel-function: near-field approximation is applied for
c    &small arguments')
      if(iss(16).eq.1) write(6,685)
 685  format(' no direct wave is calculated if mdeck.ge.iso')
      if(iss(16).ne.1) write(6,686)
 686  format(' the direct wave is calculated even if mdeck.ge.iso')
      if(iss(17).eq.0) write(6,687)
 687  format(' explosive source')
      if(iss(17).eq.1) write(6,688)
 688  format(' double couple source')
      if(iss(17).eq.2) write(6,788)
 788  format(' point source radiating p-,sv-and sh-energy,p-independent'
     &)
      if(iss(17).eq.4) write(6,689)
 689  format(' line source (2-dimensional)')
      if(iss(18).eq.0) write(6,690)
 690  format(' fuchs-mueller signal')
      if(iss(18).eq.1) write(6,691)
 691  format(' delta-impulse')
      if(iss(18).eq.2) write(6,692)
 692  format(' step-impulse')
      if(iss(18).eq.3) write(6,693)
 693  format(' momentfunction after bruestle')
      if(iss(18).eq.4) write(6,694)
 694  format(' no source-function = spike')
      if(iss(18).eq.5) write(6,695)
 695  format(' digitized source signal')
      if(iss(18).eq.6) write(6,796)
 796  format(' ricker wavelet')
      if(iss(21).eq.0) write(6,696)
 696  format(' no multiple reflections between source and receiver')
      if(iss(21).eq.1) write(6,697)
 697  format(' multiple reflections between source and receiver are cons
     &idered')
      npa1=npa
      nfr1=nfr
      nla1=nla
      ndis1=ndis
      write(6,698) npa1,nfr1,nla1,ndis1
 698  format(//' max. number of ray parameters: ',i6/' max. number of fr
     &equencies: ',i6/' max. number of layers: ',i6/' max. number of dis
     &tancies: ',i6)
      call erdmod(iins,nh)
      if(iss(17).eq.1) then
c
c     empfaenger- und herdkoordinaten und double couple-orientierung
c
         read(iins,504) dre
         do 104 i=1,nh
         read(iins,504) xs(i),ys(i),zs(i),ts(i),es(i)
c
c
c
         ztrue(i)=zs(i)
 104     if(iss(12).eq.1) zs(i)=3390.*alog(3390./(3390.-zs(i)))
c104     if(iss(12).eq.1) zs(i)=3390.*(1.-exp(-zs(i)/3390.))
c
c
c
c
         dso=zs(1)
         nso=1
         nre=1
         do 105 i=1,nlay
         diffs=dso-z(i)
         if(diffs.lt.0) iso=nso
         if(diffs.ge.0) nso=nso+1
         diffr=dre-z(i)
         if(diffr.lt.0) ire=nre
         if(diffr.ge.0) nre=nre+1
 105     continue
         if(dre.eq.0.) ire=0
         write(6,612) iso
         write(6,613) (i,xs(i),ys(i),zs(i),ts(i),es(i),i=1,nh)
 612     format(//' index of the layer with the sources',i5//
     &   ' coordinates of the sources (i,xs,ys,zs,ts,es): '/)
 613     format(i10,5f10.4)
         if(ire.eq.0) write(6,617)
         if(ire.ne.0) write(6,618) ire,dre
 617     format(//,1x,'the receivers are situated at the free surface')
 618     format(//,1x,'index of the layer with the receiver = ',i5/1x,
     1'        depth of the receiver = ',f8.2,'  km',//)
         if(iss(19)-1) 101,102,103
 101     read(iins,504) m11,m12,m13,m22,m23,m33
         write(6,615) m11,m12,m13,m12,m22,m23,m13,m23,m33
         go to 51
 102     read(iins,504) phis,delta,lamda
	 delta=delta/57.29578
	 lamda=lamda/57.29578
	 phis=phis/57.29578
         write(6,616) phis,delta,lamda
         sind=sin(delta)
         cosd=cos(delta)
         sinl=sin(lamda)
         cosl=cos(lamda)
         sinp=sin(phis)
         cosp=cos(phis)
         sin2d=sin(2.*delta)
         sin2p=sin(2.*phis)
         m11=-sind*cosl*sin2p-sin2d*sinl*sinp*sinp
         m12=sind*cosl*cos(2.*phis)+0.5*sin2d*sinl*sin2p
         m13=-cosd*cosl*cosp-cos(2.*delta)*sinl*sinp
         m22=sind*cosl*sin2p-sin2d*sinl*cosp*cosp
         m23=-cosd*cosl*sinp+cos(2.*delta)*sinl*cosp
         m33=sin2d*sinl
         write(6,615) m11,m12,m13,m12,m22,m23,m13,m23,m33
         go to 51
 103     read(iins,504) f1,f2,f3,an1,an2,an3
         write(6,614) f1,f2,f3,an1,an2,an3
         m11=2.*f1*an1
         m22=2.*f2*an2
         m33=2.*f3*an3
         m12=f1*an2+f2*an1
         m13=f1*an3+f3*an1
         m23=f2*an3+f3*an2
         write(6,615) m11,m12,m13,m12,m22,m23,m13,m23,m33
         go to 51
 614     format(//' double couple (f1,f2,f3,an1,an2,an3):'//,6f10.4)
 615     format(//' moment tensor(m11,....):'//,3(5x,3f10.4/))
 616     format(//' strike: ',f10.4/' dip:    ',f10.4/'slip:   ',f10.4)
         end if
c
c     read distances
c
 51   read(iins,503) x1,x2,dx,azi,nent
 503  format(4f10.4,i10)
      if(nent.eq.0) then
         nent=(x2-x1)/dx+1
         x(1)=x1
         az(1)=azi
         do 100 i=2,nent
         az(i)=azi
 100     x(i)=x(i-1)+dx
      else
         read(iins,504) (x(i),i=1,nent)
         read(iins,504) (az(i),i=1,nent)
 504     format(8f10.4)
      end if
      if(nent.gt.ndis) nent=ndis
c
c     write out distances
c
      write(6,601)
 601  format(1x,'list of the distances :'/)
      write(6,602) (x(i),i=1,nent)
      write(6,661)
 661  format(1x,'list of the azimuths :'/)
      write(6,602) (az(i),i=1,nent)
 602  format(1x,7f10.3)
c
c     reduction velocity vred and minimal reduced time tmin
c
      read(iins,505) vred,tmin
 505  format(2f10.3)
      if(vred.eq.0.) vred=100000.
      write(6,603) vred,tmin
 603  format(///,1x,'reduction velocity vred =',f10.3,' km/sec'//,1x,
     1'the seismograms begin at the reduced time tmin =',f10.3,' sec'//)
c
c     read integration parameter
c     calculate slowness increment
c
      read(iins,506) c2,cwir,cwil,c1,np
 506  format(4f10.4,i10)
      if(np.gt.npa) np=npa
      if(c1) 400,400,401
 400  c1=9999.999
 401  c2r=1./c2
      c1r=1./c1
      cwilr=1./cwil
      cwirr=1./cwir
      np1=np-1
      dp=(c2r-c1r)/float(np1)
      ifenr=ifix((c2r-cwirr)/dp)
      ifenl=ifix((cwilr-c1r)/dp)
      if(ifenl.le.1) then
         ifenl=2
         cwilr=float(ifenl)*dp+c1r
         cwil=1./cwilr
      end if
      if(ifenr.le.1) then
         ifenr=2
         cwirr=-float(ifenr)*dp+c2r
         cwir=1./cwirr
      end if
      diss=float(ifenl)/float(ifenr)
      nfensr=np-ifenr
      nfensl=ifenl+1
      dfen=pi/float(ifenl)
      do 402 i=1,np
      p(i)=float(i-1)*dp+c1r
      p2(i)=p(i)*p(i)
      if(i.lt.nfensr.and.i.gt.nfensl) then
         window(i)=1.
      else if(i.le.nfensl) then
         fif=i-1
         window(i)=(0.5-cos(dfen*fif)/2.)**iss(10)
      else
         fif=i-nfensr
         window(i)=(cos(dfen*diss*fif)/2.+0.5)**iss(10)
      end if
c
c     *******      line-source      ********
c
      if(iss(17).eq.4) window(i)=2.*window(i)/p(i)
c
 402  continue
      write(6,604)c2,c1,np,dp
 604  format(//1x,'integration parameter :'//1x,'minimum phase ',
     &'velocity c2 =',f10.3,' km/sec'/1x,
     &'maximum phase velocity c1 =',f10.3,' km/sec  ',
     &/1x,'number of slownessses used for integration np =',i5/1x,
     &'slowness increment dp =',f10.6,' sec/km'/)
      write(6,605)iss(10),c2,cwir
      write(6,605)iss(10),cwil,c1
 605  format(1x,'cos**',i2,' window from ',f10.4,' km/sec to ',f10.4,' k
     1m/sec'/)
c
c     read frequency passband
c
      read(iins,507) fu,fwil,fwir,fo,fr,gamma
 507  format(6f10.4)
      if(iss(17).eq.1.and.fr.ne.0.) then
         fr=0.
         write(6,650)
 650     format(' in the case of a double couble source'/' the causal co
     &mplex velocity law is not implemented')
      end if
      if(fwil.eq.0.)  fwil=fu
      if(fwir.eq.0.)  fwir=fo
      write(6,606)fu,fo,fr,gamma
 606  format(//1x,'frequency window :'//1x,'low frequency cut at ',
     1f5.2,' hz'/1x,'high frequency cut at ',f5.2,' hz'/1x,'reference fr
     2equency =',f5.2,' hz'/1x,'exponent gamma for power law =',f5.2)
      if(iss(11).ne.0)write(6,607)iss(11),fu,fwil,iss(11),fwir,fo
 607  format(1x,'cos**',i2,' taper from ',f5.2,' hz to ',f5.2,' hz'/
     11x,'cos**',i2,' taper from ',f5.2,' hz to ',f5.2,' hz')
c
c     read signal parameters
c

      read(iins,508) dt,npts,na,n,t,tsigma
 508  format(f10.4,3i10,2f10.4)
      write(6,608) dt,npts,na,n,t,tsigma
 608  format(//1x,'signal parameter :'//,12x,'dt',8x,'npts',8x,'na',9x,
     &'n',4x,'tsignal'/1x,f15.3,3i10,2f10.4)
      omegai=0.
      if(tsigma.ne.0.)omegai=-1./tsigma
      fnpt=npts
      dauer=dt*fnpt
      dnue=1./dauer
      do 599 i=1,npts
 599  s(i)=(0.,0.)
      if(iss(18).eq.0) then
c     mueller-signal
         call signa0(t,dt,n,npts,na,s,tsigma)
      else if(iss(18).eq.1.or.iss(18).eq.2) then
c     delta-impuls oder stufen-impuls
         call signa1(t,dt,s,tsigma)
      else if(iss(18).eq.3) then
c     bruestle-momentensignal
         call signa3(t,dt,s,tsigma)
      else if(iss(18).eq.4) then
c     no source-signal - spike
         s(1)=cmplx(1.0,0.)
      else if(iss(18).eq.5) then
c     digitized source signal read in from cards
         call signa5(t,dt,s,tsigma)
      end if
      do 600 i=1,npts
 600  s(i)=s(i)/fnpt
c
c     read other meta-data
c
!       read(iins,509) elat,elon
! 509   format(2f10.4)
!       read(iins,510) (slat(i),i=1,nent)
! 510   format(8f10.4)
!       write(6,511)
! 511   format(1x,'list of the lat :'/)
!       write(6,512) (slat(i),i=1,nent)
! 512   format(1x,7f10.3)
!       read(iins,513) (slon(i),i=1,nent)
! 513   format(8f10.4)
!       write(6,514)
! 514   format(1x,'list of the lon :'/)
!       write(6,515) (slon(i),i=1,nent)
! 515   format(1x,7f10.3)
!       read(iins,516) (baz(i),i=1,nent)
! 516   format(8f10.4)
!       write(6,517)
! 517   format(1x,'list of the baz :'/)
!       write(6,518) (baz(i),i=1,nent)
! 518   format(1x,7f10.3)
!       read(iins,519) (tp1(i),i=1,nent)
! 519   format(8f10.4)
!       write(6,520)
! 520   format(1x,'list of the SS :'/)
!       write(6,521) (tp1(i),i=1,nent)
! 521   format(1x,7f10.3)
!       read(iins,522) (tp2(i),i=1,nent)
! 522   format(8f10.4)
!       write(6,523)
! 523   format(1x,'list of the S410S :'/)
!       write(6,524) (tp2(i),i=1,nent)
! 524   format(1x,7f10.3)
!       read(iins,525) (tp3(i),i=1,nent)
! 525   format(8f10.4)
!       write(6,526)
! 526   format(1x,'list of the S660S :'/)
!       write(6,527) (tp3(i),i=1,nent)
! 527   format(1x,7f10.3)
!       read(iins,528) (tp4(i),i=1,nent)
! 528   format(8f10.4)
!       write(6,529)
! 529   format(1x,'list of the ScSScS :'/)
!       write(6,530) (tp4(i),i=1,nent)
! 530   format(1x,7f10.3)
!       read(iins,531) (tp5(i),i=1,nent)
! 531   format(8f10.4)
!       write(6,532)
! 532   format(1x,'list of the Sv660sS :'/)
!       write(6,533) (tp5(i),i=1,nent)
! 533   format(1x,7f10.3)
c
c     calculation of the complex spectrum
c
      call zzlog(npts,mga)
      if(iss(18).ne.6) call cool(mga,s,-1.0)
      inpt = npts/2
      iu = fu*dauer+1.
      if(iu.eq.1)  iu=2
      io = fo*dauer + 1.
      ileft = fwil*dauer + 1.
      if(ileft.eq.1)  ileft=2
      iright = fwir*dauer + 1.
      nleft = ileft - iu
      nright = io - iright
      if(nleft.gt.0) dleft = pi/float(nleft)
      if(nright.gt.0) dright = pi/float(nright)
      if(iss(18).eq.6) then
         d1=1./dauer
         do 720 i=1,inpt
         om=float(i)*d1*pi2
         s(i)=om*om*exp(-om*om*t*t/(pi2*pi2))
 720     continue
      end if
      do 700 i=1,inpt
      fif = i - iright
      faf = i - iu
      if(i.ge.ileft.and.i.le.iright) wavel(i) = s(i)
      if(i.le.iu.or.i.ge.io) wavel(i) = (0.0,0.0)
      if(i.gt.iu.and.i.lt.ileft)
     1wavel(i) = s(i)*(0.5-cos(dleft*faf)/2.)**iss(11)
      if(i.gt.iright.and.i.lt.io)
     1wavel(i) = s(i)*(cos(dright*fif)/2.+0.5)**iss(11)
 700  continue
      nfreq=io-iu+1
      if(nfreq.gt.nfr) go to 6591
c
c     write on output file
c
cc      call output(dentif,nent,iouts)
cc      call output(dentif,nent,ioutsh)
      if(iss(6).eq.1) then
         call output(dentif,nent,21)
         call output(dentif,nent,22)
         call output(dentif,nent,23)
      end if
      dauer1=pi2/dauer
      iu1=iu-1
      do 800 i=1,nfreq
      omega(i)=float(iu1)*dauer1
      omegac(i)=sqrt(omega(i)*omega(i)+omegai*omegai)
      phasc(i)=atan(omegai/omega(i))
      sphasc(i)=sin(phasc(i))
      cphasc(i)=cos(phasc(i))
      phasc(i)=phasc(i)*0.5-pi4
      if(iss(17).eq.0.or.iss(17).eq.2) then
         omega2(i)=(omega(i)*omega(i)+omegai*omegai)
         phase1(i)=2.*atan(omegai/omega(i))
      else if(iss(17).eq.4) then
         omega2(i)=sqrt(omega(i)*omega(i)+omegai*omegai)
         phase1(i)=atan(omegai/omega(i))
      else if(iss(18).ne.3.and.iss(18).ne.5) then
         omega2(i)=(omega(i)*omega(i)+omegai*omegai)**0.75
         phase1(i)=1.5*atan(omegai/omega(i))
      else
         omega2(i)=-1./(omega(i)*omega(i)+omegai*omegai)**0.25
         phase1(i)=-0.5*atan(omegai/omega(i))
      end if
 800  iu1=iu1+1
      omegar=pi2*fr
      write(*,*) 'fede', fr, omegar, gamma
      if(omegar.ne.0.and.gamma.eq.0.) then
         do 810 i=1,nfreq
         bomega=sqrt(omega(i)**2+omegai**2)
         pomega=atan(omegai/omega(i))
         cv(i)=(alog(bomega/omegar)+ein*(pomega+1.570796))/pi
         cvr(i)=real(cv(i))
         cvi(i)=aimag(cv(i))
         x1r(i)=real(cv(i))
         x1i(i)=aimag(cv(i))
 810     continue
         ii=dauer/t+1
         do 812 j=1,nlay
         alpr(j)=alp(j)*(1.+x1r(ii)/qp(j))
         alpi(j)=alp(j)*x1i(ii)/qp(j)
         betr(j)=bet(j)*(1.+x1r(ii)/qs(j))
         beti(j)=bet(j)*x1i(ii)/qs(j)
 812     continue
      else if(omegar.ne.0.and.gamma.ne.0.) then
         cgam=1./tan(gamma*pi*0.5)
         do 820 i=1,nfreq
         x1r(i)=0.5*(omegar/omega(i))**gamma
         cvr(i)=cgam*(0.5-x1r(i))
         cvi(i)=x1r(i)
         cv(i)=cmplx(cvr(i),cvi(i))
 820     continue
         do 822 j=1,nlay
         alpr(j)=alp(j)
         betr(j)=bet(j)
 822     continue
      else
         do 814 j=1,nlay
         alpr(j)=alp(j)
         betr(j)=bet(j)
 814     continue
      end if
      write(*,*) 'Here comes the dragon!'
      do j =1,nlay
           if (z(j)<=1500.) then
               write(*,*) z(j), j, betr(j), beti(j), alpr(j), alpi(j)
           endif
      enddo
!      STOP
      jf1=0
      al(iso)=cmplx(alpr(iso),alpi(iso))
      do 805 jf=1,inpt
      if(jf.lt.iu.or.jf.gt.io) go to 805
      jf1=jf1+1
      wavel(jf)=wavel(jf)*omega2(jf1)*cexp(ein*phase1(jf1))
      if(iss(3).eq.1) wavel(jf)=wavel(jf)*cmplx(-omegai,omega(jf1))
      if(iss(3).eq.2) wavel(jf)=-wavel(jf)*cmplx(omega(jf1),omegai)*
     &cmplx(omega(jf1),omegai)
      if(iss(18).eq.0) then
         wavel(jf)=wavel(jf)*al(iso)**2
      else if(iss(18).eq.2) then
         wavel(jf)=wavel(jf)/cmplx(-omegai,omega(jf1))
      else if(iss(18).eq.5.or.iss(18).eq.6) then
         wavel(jf)=wavel(jf)/cmplx(-omegai,omega(jf1))
         wavel(jf)=wavel(jf)/cmplx(-omegai,omega(jf1))
      end if
 805  continue
      if(iss(21).eq.1) then
         if(mdeck.gt.ire.or.mdeck.gt.iso) iss(21)=0
         if(mdeck.gt.ire.or.mdeck.gt.iso) write(6,619)
      end if
 619  format(' switch iss(21) has been changed to 0 because mdeck is gre
     &ater then'/' iso(source layer) or ire(receiver layer)')
      if(iss(5).ne.1) then
c
c     ber. der sh-seismogramme
c
         call shseis(ioutsh)
         if(iss(4).eq.3) go to 6602
      end if
      if(iss(17).eq.1) then
         d11x=-1./(pi2*rho(iso)*sqrt(pi2))
      else
         d11x=1.
      end if
      do 1000 i=1,np
      id=-1
      ih=1
      xx1=1.
      u=p(i)

c
c     bestimmung von ua und ub, ire=receiver
c
      if(ire.ne.0) then
         diffr=dre-z(ire)
         al(ire)=cmplx(alpr(ire),alpi(ire))
         be(ire)=cmplx(betr(ire),beti(ire))
         am=csqrt(1./(al(ire)*al(ire))-p2(i))
         if(real(am).eq.0.) am=-am
         amr=real(am)
         ami=aimag(am)
c        amdr=amr*diffr
c        amdi=ami*diffr
         bm=csqrt(1./(be(ire)*be(ire))-p2(i))
         if(real(bm).eq.0.) bm=-bm
         bmr=real(bm)
         bmi=aimag(bm)
c        bmdr=bmr*diffr
c        bmdi=bmi*diffr
         if(omegar.eq.0.) then
            call phas3(diffr,omegai,am,bm,x1r,x1i,x2r,x2i,x3r,x3i,x4r,
     &                 x4i,omega,nfreq)
         else
            p2i=p2(i)
            call phas4(diffr,omegai,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,
     &                 p2i,omega,nfreq,ire)
         end if
         do 1001 j=1,nfreq
         a11r(j)=-x1i(j)-x3r(j)
         a11i(j)=x1r(j)-x3i(j)
         a22r(j)=-x2i(j)-x4r(j)
         a22i(j)=x2r(j)-x4i(j)
c        a11r(j)=-omega(j)*amdi-omegai*amdr
c        a11i(j)=omega(j)*amdr-omegai*amdi
c        a22r(j)=-omega(j)*bmdi-omegai*bmdr
c        a22i(j)=omega(j)*bmdr-omegai*bmdi
         b11r(j)=exp(a11r(j))*cos(a11i(j))
         b11i(j)=exp(a11r(j))*sin(a11i(j))
         b12r(j)=exp(a22r(j))*cos(a22i(j))
         b12i(j)=exp(a22r(j))*sin(a22i(j))
         b22r(j)=exp(-a11r(j))*cos(a11i(j))
         b22i(j)=-exp(-a11r(j))*sin(a11i(j))
         b21r(j)=exp(-a22r(j))*cos(a22i(j))
         b21i(j)=-exp(-a22r(j))*sin(a22i(j))
         ua1r(j)=p(i)*b11r(j)
         ua1i(j)=p(i)*b11i(j)
         ua2r(j)=bmr*b12r(j)-bmi*b12i(j)
         ua2i(j)=bmi*b12r(j)+bmr*b12i(j)
         ua3r(j)=amr*b11r(j)-ami*b11i(j)
         ua3i(j)=ami*b11r(j)+amr*b11i(j)
         ua4r(j)=-p(i)*b12r(j)
         ua4i(j)=-p(i)*b12i(j)
         ub1r(j)=p(i)*b22r(j)
         ub1i(j)=p(i)*b22i(j)
         ub2r(j)=-bmr*b21r(j)+bmi*b21i(j)
         ub2i(j)=-bmi*b21r(j)-bmr*b21i(j)
         ub3r(j)=-amr*b22r(j)+ami*b22i(j)
         ub3i(j)=-ami*b22r(j)-amr*b22i(j)
         ub4r(j)=-p(i)*b21r(j)
         ub4i(j)=-p(i)*b21i(j)
 1001    continue
      else
         call set(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,1.,nfreq)
         call set(ub1r,ub1i,ub2r,ub2i,ub3r,ub3i,ub4r,ub4i,0.,nfreq)
      end if
      if(iss(17).eq.1) go to 1003
c
c******ber.des quellensignals mit zeitl.verschiebung ea und eb*********
c
      p3(i)=p2(i)*p(i)
      diffs=dso-z(iso-1)
      al(iso)=cmplx(alpr(iso),alpi(iso))
      be(iso)=cmplx(betr(iso),beti(iso))
      beal=betr(iso)**3/alpr(iso)**3
      am=csqrt(1./(al(iso)*al(iso))-p2(i))
      if(real(am).eq.0.) am=-am
      amr=real(am)
      ami=aimag(am)
c     amdr=amr*diffs
c     amdi=ami*diffs
      amb=1./(amr*amr+ami*ami)
      bm=csqrt(1./(be(iso)*be(iso))-p2(i))
      if(real(bm).eq.0.) bm=-bm
      bmr=real(bm)
      bmi=aimag(bm)
      bmb=1./(bmr*bmr+bmi*bmi)
c     bmdr=bmr*diffs
c     bmdi=bmi*diffs
         if(omegar.eq.0.) then
            call phas3(diffs,omegai,am,bm,x1r,x1i,x2r,x2i,x3r,x3i,x4r,
     &                 x4i,omega,nfreq)
         else
            p2i=p2(i)
            call phas4(diffs,omegai,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,
     &                 p2i,omega,nfreq,iso)
         end if
      do 1002 j=1,nfreq
      a11r(j)=-x1i(j)-x3r(j)
      a11i(j)=x1r(j)-x3i(j)
      a22r(j)=-x2i(j)-x4r(j)
      a22i(j)=x2r(j)-x4i(j)
c     a11r(j)=-omega(j)*amdi-omegai*amdr
c     a11i(j)=omega(j)*amdr-omegai*amdi
c     a22r(j)=-omega(j)*bmdi-omegai*bmdr
c     a22i(j)=omega(j)*bmdr-omegai*bmdi
c     b11r(j)=exp(2.*a11r(j))*cos(2.*a11i(j))
c     b11i(j)=exp(2.*a11r(j))*sin(2.*a11i(j))
c     b12r(j)=exp(2.*a22r(j))*cos(2.*a22i(j))
c     b12i(j)=exp(2.*a22r(j))*sin(2.*a22i(j))
      b11r(j)=exp(a11r(j))*cos(a11i(j))
      b11i(j)=exp(a11r(j))*sin(a11i(j))
      b12r(j)=exp(a22r(j))*cos(a22i(j))
      b12i(j)=exp(a22r(j))*sin(a22i(j))
      b22r(j)=exp(-a11r(j))*cos(a11i(j))
      b22i(j)=-exp(-a11r(j))*sin(a11i(j))
      b21r(j)=exp(-a22r(j))*cos(a22i(j))
      b21i(j)=-exp(-a22r(j))*sin(a22i(j))
      s11ur(j)=p3(i)*amb*(ami*b22r(j)-amr*b22i(j))
      s11ui(j)=p3(i)*amb*(amr*b22r(j)+ami*b22i(j))
      s12ur(j)=-p2(i)*b21i(j)
      s12ui(j)=p2(i)*b21r(j)
      s21ur(j)=p(i)*(-ami*b22r(j)-amr*b22i(j))
      s21ui(j)=p(i)*(amr*b22r(j)-ami*b22i(j))
      s22ur(j)=-s12ur(j)
      s22ui(j)=-s12ui(j)
      s11dr(j)=p3(i)*amb*(ami*b11r(j)-amr*b11i(j))
      s11di(j)=p3(i)*amb*(amr*b11r(j)+ami*b11i(j))
      s12dr(j)=p2(i)*b12i(j)
      s12di(j)=-p2(i)*b12r(j)
      s21dr(j)=p(i)*(-ami*b11r(j)-amr*b11i(j))
      s21di(j)=p(i)*(amr*b11r(j)-ami*b11i(j))
      s22dr(j)=-s12dr(j)
      s22di(j)=-s12di(j)
c     s11dr(j)=s11ur(j)*b11r(j)-s11ui(j)*b11i(j)
c     s11di(j)=s11ui(j)*b11r(j)+s11ur(j)*b11i(j)
c     s12dr(j)=-s12ur(j)*b12r(j)+s12ui(j)*b12i(j)
c     s12di(j)=-s12ui(j)*b12r(j)-s12ur(j)*b12i(j)
c     s21dr(j)=s21ur(j)*b11r(j)-s21ui(j)*b11i(j)
c     s21di(j)=s21ui(j)*b11r(j)+s21ur(j)*b11i(j)
c     s22dr(j)=-s12dr(j)
c     s22di(j)=-s12di(j)
 1002 continue
      if(iss(17).eq.2.and.iss(1).ne.1) then
c
c     point source for p-and sv-waves
c
         do 1102 j=1,nfreq
         s12ur(j)=p3(i)*bmb*(bmi*b21r(j)-bmr*b21i(j))*beal
         s12ui(j)=p3(i)*bmb*(bmr*b21r(j)+bmi*b21i(j))*beal
         s22ur(j)=p(i)*(-bmi*b21r(j)-bmr*b21i(j))*beal
         s22ui(j)=p(i)*(bmr*b21r(j)-bmi*b21i(j))*beal
         s12dr(j)=p3(i)*bmb*(bmi*b12r(j)-bmr*b12i(j))*beal
         s12di(j)=p3(i)*bmb*(bmr*b12r(j)+bmi*b12i(j))*beal
         s22dr(j)=p(i)*(-bmi*b12r(j)-bmr*b12i(j))*beal
         s22di(j)=p(i)*(bmr*b12r(j)-bmi*b12i(j))*beal
 1102    continue
         if(iss(1).eq.2) then
            do 1103 j=1,nfreq
            s11ur(j)=0.
            s11ui(j)=0.
            s21ur(j)=0.
            s21ui(j)=0.
            s11dr(j)=0.
            s11di(j)=0.
            s21dr(j)=0.
 1103       s21di(j)=0.
         end if
      end if
      if(iss(7).eq.1) then
c     no upgoing waves from the source are considered
         do 1012 j=1,nfreq
         s11ur(j)=0.
         s11ui(j)=0.
         s12ur(j)=0.
         s12ui(j)=0.
         s21ur(j)=0.
         s21ui(j)=0.
         s22ur(j)=0.
         s22ui(j)=0.
 1012    continue
      end if
 1003 iso1=iso-1
      if(iso1.lt.ire) go to 1200
      if(iso1.gt.ire) go to 1100
      if(iso1.eq.ire.and.ire.ne.0) go to 1100
c***********************************************************************
c     source and receiver are situated at the free surface
c***********************************************************************
c     turs=1.    und   rdrs=0.
c
      call freesu(alpr,alpi,betr,beti,u,tppu1r,tppu1i,tspu1r,
     &tspu1i,tssu1r,tssu1i,tpsu1r,tpsu1i,nfreq)
      if(mdeck.lt.iso) then
c
c     ber. der reflmatrix rdsl
c
      call rekurs(iso,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,1)
c
c     ber. der reflmatrix rusm (rusf,falls mdeck=0)
c
      call rekurs(iso,mdeck,ih,alpr,alpi,betr,beti,rho,d,u,rppu2r,
     &rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,0)
c
c     ber. von a1=tusr*(1-rdsl*rusf)**-1
c
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,rppu2r,rppu2i,rspu2r,rspu2i,rpsu2r,rpsu2i,rssu2r,rssu2i,
     &b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,nfreq)
      call matinv(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,d11r,d11i,
     &d12r,d12i,d21r,d21i,d22r,d22i,nfreq)
      call mmtml(tppu1r,tppu1i,tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,
     &tssu1i,d11r,d11i,d12r,d12i,d21r,d21i,d22r,d22i,a11r,a11i,a12r,
     &a12i,a21r,a21i,a22r,a22i,nfreq)
      go to 1330
      else
c
c     ber. der reflmatrix rdsl=tums*rdml*tdsm
c
      call rekurs(mdeck,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,1)
      if(mdeck.eq.iso) go to 1010
      call transm(iso,mdeck,id,alpr,alpi,betr,beti,rho,d,u,tppd2r,
     &tppd2i,tspd2r,tspd2i,tssd2r,tssd2i,tpsd2r,tpsd2i,1)
      call transm(mdeck,iso,ih,alpr,alpi,betr,beti,rho,d,u,tppu2r,
     &tppu2i,tspu2r,tspu2i,tssu2r,tssu2i,tpsu2r,tpsu2i,1)
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,tppd2r,tppd2i,tspd2r,tspd2i,tpsd2r,tpsd2i,tssd2r,tssd2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call mmtml(tppu2r,tppu2i,tspu2r,tspu2i,tpsu2r,tpsu2i,tssu2r,
     &tssu2i,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,rppd2r,rppd2i,
     &rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,nfreq)
c
c     ber. von y1=ua*tusr*rdsl, xx1=0.
c
 1010 call mmtml(tppu1r,tppu1i,tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,
     &tssu1i,rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call mmtml(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,a11r,a11i,
     &a12r,a12i,a21r,a21i,a22r,a22i,y11r,y11i,y12r,y12i,y21r,y21i,y22r
     &,y22i,nfreq)
      xx1=0.
      go to 1400
      end if
c***********************************************************************
c     the source is situated under the receiver
c     case - receiver at the free surface - is included
c***********************************************************************
c
c     ber. der transmissivitaetsmatrix tusr
c
 1100 if(ire.ne.0.and.iss(21).eq.1) then
         call transn(iso,ire,ih,alpr,alpi,betr,beti,rho,d,u,tppu1r,
     &   tppu1i,tspu1r,tspu1i,tssu1r,tssu1i,tpsu1r,tpsu1i,
     &   rppu2r,rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,0)
      else if(ire.eq.0.and.iss(21).eq.1) then
         call transn(iso,ire,ih,alpr,alpi,betr,beti,rho,d,u,c11r,
     &   c11i,c12r,c12i,c22r,c22i,c21r,c21i,
     &   rppu2r,rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,0)
      else if(ire.ne.0) then
         call transm(iso,ire,ih,alpr,alpi,betr,beti,rho,d,u,tppu1r,
     &   tppu1i,tspu1r,tspu1i,tssu1r,tssu1i,tpsu1r,tpsu1i,0)
      else
         call transm(iso,1,ih,alpr,alpi,betr,beti,rho,d,u,c11r,c11i,c12r
     &   ,c12i,c22r,c22i,c21r,c21i,1)
      end if
      if(ire.eq.0) then
         call freesu(alpr,alpi,betr,beti,u,b11r,b11i,b12r,b12i,
     &   b22r,b22i,b21r,b21i,nfreq)
         call mmtml(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,c11r,c11i,
     &   c12r,c12i,c21r,c21i,c22r,c22i,tppu1r,tppu1i,tspu1r,tspu1i,
     &   tpsu1r,tpsu1i,tssu1r,tssu1i,nfreq)
      end if
      if(mdeck.lt.ire) then
c
c     ber. der reflmatrix rdrs
c
      call rekurs(ire,iso,id,alpr,alpi,betr,beti,rho,d,u,rppd1r,rppd1i,
     &rspd1r,rspd1i,rssd1r,rssd1i,rpsd1r,rpsd1i,0)
c
c     ber. der reflmatrix rurm (rurf,falls mdeck=0)
c
      call rekurs(ire,mdeck,ih,alpr,alpi,betr,beti,rho,d,u,rppu1r,
     &rppu1i,rspu1r,rspu1i,rssu1r,rssu1i,rpsu1r,rpsu1i,1)
c
c     ber. der reflmatrix rdsl
c
      call rekurs(iso,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,1)
c
c     ber. der reflmatrix rusm (rusf,falls mdeck=0)
c
      call rekurs(iso,mdeck,ih,alpr,alpi,betr,beti,rho,d,u,rppu2r,
     &rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,0)
c
c     ber. von (1-rdrs*rurf)**-1
c
      call mmtml(rppd1r,rppd1i,rspd1r,rspd1i,rpsd1r,rpsd1i,rssd1r,
     &rssd1i,rppu1r,rppu1i,rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,b11r,b11i,
     &b12r,b12i,b21r,b21i,b22r,b22i,nfreq)
c
c     ber. von (1-rdsl*rusf)**-1
c
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,rppu2r,rppu2i,rspu2r,rspu2i,rpsu2r,rpsu2i,rssu2r,rssu2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,c11r,c11i,
     &c12r,c12i,c21r,c21i,c22r,c22i,nfreq)
c
c     ber. von a1=(1-rdrs*rurf)**-1*tusr*(1-rdsl*rusf)**-1
c
      call mmtml(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,tppu1r,tppu1i,
     &tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,tssu1i,d11r,d11i,d12r,d12i,
     &d21r,d21i,d22r,d22i,nfreq)
      call mmtml(d11r,d11i,d12r,d12i,d21r,d21i,d22r,d22i,c11r,c11i,c12r
     &,c12i,c21r,c21i,c22r,c22i,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
     &nfreq)
      go to 1320
      else if(mdeck.lt.iso) then

c*****rurf=0  ---  rdrs braucht nicht berechnet zu werden
c*****mdeck=ire ist eingeschlossen
c
c     ber. der reflmatrix rdsl
c
      call rekurs(iso,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,1)
c
c     ber. der reflmatrix rusm
c
      if(iss(21).ne.1) call rekurs(iso,mdeck,ih,alpr,alpi,betr,beti,rho,
     &d,u,rppu2r,rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,0)
c
c     ber. von a1=tusr*(1-rdsl*rusm)**-1
c
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,rppu2r,rppu2i,rspu2r,rspu2i,rpsu2r,rpsu2i,rssu2r,rssu2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,c11r,c11i,
     &c12r,c12i,c21r,c21i,c22r,c22i,nfreq)
      call mmtml(tppu1r,tppu1i,tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,
     &tssu1i,c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,a11r,a11i,a12r,
     &a12i,a21r,a21i,a22r,a22i,nfreq)
      go to 1330
      else

c*****rurf=0   rdrs=0   rusf=0
c
c     ber. von rdsl=tums*rdml*tdsm
c
      call rekurs(mdeck,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,1)
      if(mdeck.eq.iso) go to 1110
      call transm(iso,mdeck,id,alpr,alpi,betr,beti,rho,d,u,tppd2r,
     &tppd2i,tspd2r,tspd2i,tssd2r,tssd2i,tpsd2r,tpsd2i,1)
      call transm(mdeck,iso,ih,alpr,alpi,betr,beti,rho,d,u,tppu2r,
     &tppu2i,tspu2r,tspu2i,tssu2r,tssu2i,tpsu2r,tpsu2i,1)
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,tppd2r,tppd2i,tspd2r,tspd2i,tpsd2r,tpsd2i,tssd2r,tssd2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call mmtml(tppu2r,tppu2i,tspu2r,tspu2i,tpsu2r,tpsu2i,tssu2r,
     &tssu2i,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,rppd2r,rppd2i,
     &rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,nfreq)
 1110 do 1150 j=1,nfreq
      a11r(j)=tppu1r(j)
      a11i(j)=tppu1i(j)
      a12r(j)=tspu1r(j)
      a12i(j)=tspu1i(j)
      a21r(j)=tpsu1r(j)
      a21i(j)=tpsu1i(j)
      a22r(j)=tssu1r(j)
 1150 a22i(j)=tssu1i(j)
      go to 1330
      end if

c***********************************************************************
c     the source is situated above the receiver
c     case - source at the free surface - is included
c***********************************************************************
c
c     ber. der transmissivtaetsmatrix tdsr
c

 1200 if(iss(21).eq.1) then
         call transn(iso,ire,id,alpr,alpi,betr,beti,rho,d,u,tppu1r,
     &   tppu1i,tspu1r,tspu1i,tssu1r,tssu1i,tpsu1r,tpsu1i,
     &   rppu2r,rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,2)
      else
         call transm(iso,ire,id,alpr,alpi,betr,beti,rho,d,u,tppu1r,
     &   tppu1i,tspu1r,tspu1i,tssu1r,tssu1i,tpsu1r,tpsu1i,2)
      end if
      if(mdeck.lt.iso) then
c
c     ber. der reflmatrix rurs
c
      call rekurs(ire,iso,ih,alpr,alpi,betr,beti,rho,d,u,rppd1r,rppd1i,
     &rspd1r,rspd1i,rssd1r,rssd1i,rpsd1r,rpsd1i,1)
c
c     ber. der reflmatrix rdrl
c
      call rekurs(ire,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppu1r,
     &rppu1i,rspu1r,rspu1i,rssu1r,rssu1i,rpsu1r,rpsu1i,0)
c
c     ber. der reflmatrix rusm (rusf,falls mdeck=0)
c
      call rekurs(iso,mdeck,ih,alpr,alpi,betr,beti,rho,d,u,rppd2r,
     &rppd2i,rspd2r,rspd2i,rssd2r,rssd2i,rpsd2r,rpsd2i,0)
c
c     ber. der reflmatrix rdsl
c
      call rekurs(iso,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppu2r,
     &rppu2i,rspu2r,rspu2i,rssu2r,rssu2i,rpsu2r,rpsu2i,1)
c
c     ber. von a1=(1-rurs*rdrl)**-1*tdsr*(1-rusf*rdsl)**-1
c
      call mmtml(rppd1r,rppd1i,rspd1r,rspd1i,rpsd1r,rpsd1i,rssd1r,
     &rssd1i,rppu1r,rppu1i,rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,b11r,b11i,
     &b12r,b12i,b21r,b21i,b22r,b22i,nfreq)
      call mmtml(rppd2r,rppd2i,rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,
     &rssd2i,rppu2r,rppu2i,rspu2r,rspu2i,rpsu2r,rpsu2i,rssu2r,rssu2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,c11r,c11i,
     &c12r,c12i,c21r,c21i,c22r,c22i,nfreq)
      call mmtml(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,tppu1r,tppu1i,
     &tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,tssu1i,d11r,d11i,d12r,d12i,
     &d21r,d21i,d22r,d22i,nfreq)
      call mmtml(d11r,d11i,d12r,d12i,d21r,d21i,d22r,d22i,c11r,c11i,c12r
     &,c12i,c21r,c21i,c22r,c22i,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
     &nfreq)
      go to 1310
      else if(mdeck.lt.ire) then
c*****rusf=0 -- rdsl ist unnoetig geworden
c*****mdeck=iso ist eingeschlossen, auch iso=1
c
c     ber. der reflmatrix rurm
c
      call rekurs(ire,mdeck,ih,alpr,alpi,betr,beti,rho,d,u,rppd1r,
     &rppd1i,rspd1r,rspd1i,rssd1r,rssd1i,rpsd1r,rpsd1i,1)
c
c     ber. der reflmatrix rdrl
c
      call rekurs(ire,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppu1r,
     &rppu1i,rspu1r,rspu1i,rssu1r,rssu1i,rpsu1r,rpsu1i,0)
c
c     ber. von a1=(1-rurm*rdrl)**-1*tdsr
c
      call mmtml(rppd1r,rppd1i,rspd1r,rspd1i,rpsd1r,rpsd1i,rssd1r,
     &rssd1i,rppu1r,rppu1i,rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,b11r,b11i,
     &b12r,b12i,b21r,b21i,b22r,b22i,nfreq)
      call mmtml(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,tppu1r,tppu1i,
     &tspu1r,tspu1i,tpsu1r,tpsu1i,tssu1r,tssu1i,a11r,a11i,a12r,a12i,
     &a21r,a21i,a22r,a22i,nfreq)
      xx1=0.
      go to 1340
      else
c*****rusf=0    rurs=0
c
c     ber. von rdrl=tumr*rdml*tdrm
c
      it=0
      if(mdeck.ne.ire) it=1
      call rekurs(mdeck,nlay,id,alpr,alpi,betr,beti,rho,d,u,rppu1r,
     &rppu1i,rspu1r,rspu1i,rssu1r,rssu1i,rpsu1r,rpsu1i,it)
      if(mdeck.eq.ire) go to 1210
      call transm(ire,mdeck,id,alpr,alpi,betr,beti,rho,d,u,tppd2r,
     &tppd2i,tspd2r,tspd2i,tssd2r,tssd2i,tpsd2r,tpsd2i,0)
      call transm(mdeck,ire,ih,alpr,alpi,betr,beti,rho,d,u,tppu2r,
     &tppu2i,tspu2r,tspu2i,tssu2r,tssu2i,tpsu2r,tpsu2i,0)
      call mmtml(rppu1r,rppu1i,rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,
     &rssu1i,tppd2r,tppd2i,tspd2r,tspd2i,tpsd2r,tpsd2i,tssd2r,tssd2i,
     &a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,nfreq)
      call mmtml(tppu2r,tppu2i,tspu2r,tspu2i,tpsu2r,tpsu2i,tssu2r,
     &tssu2i,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,rppu1r,rppu1i,
     &rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,nfreq)
 1210 do 1250 j=1,nfreq
      a11r(j)=tppu1r(j)
      a11i(j)=tppu1i(j)
      a12r(j)=tspu1r(j)
      a12i(j)=tspu1i(j)
      a21r(j)=tpsu1r(j)
      a21i(j)=tpsu1i(j)
      a22r(j)=tssu1r(j)
 1250 a22i(j)=tssu1i(j)
      xx1=0.
      go to 1340
      end if
c***********************************************************************
c     ende der berechnung der einzelnen matrizen
c***********************************************************************
c
c     ber. der vollstaendigen seismogramme
c
c     multiplikation von ua mit rdrl - anteil der von unten einf. (auf-
c     tauchenden welle) - upgoing - b1
c
 1310 call mmtml(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,rppu1r,rppu1i,
     &rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,b11r,b11i,b12r,b12i,
     &b21r,b21i,b22r,b22i,nfreq)
c
c     addition von ub und ua*rdrl - c1
c
      call matadd(ub1r,ub1i,ub2r,ub2i,ub3r,ub3i,ub4r,ub4i,b11r,b11i,b12r
     &,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,
     &nfreq)
c
c     multiplikation von (ub+ua*rdrl) und a1 - y1
c
      call mmtml(c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,a11r,a11i,
     &a12r,a12i,a21r,a21i,a22r,a22i,y11r,y11i,y12r,y12i,y21r,y21i,y22r,
     &y22i,nfreq)
c
c     multiplikation von y1 mit rusf fuer upgoing wave-field
c                            - x1
c
      call mmtml(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,rppd2r,rppd2i,
     &rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,x11r,x11i,x12r,x12i,
     &x21r,x21i,x22r,x22i,nfreq)
      go to 1400
c
c     multiplikation von ub mit rufr - anteil der von oben einf. (ab-
c     tauchenden welle) - downgoing - b1
c
 1320 call mmtml(ub1r,ub1i,ub2r,ub2i,ub3r,ub3i,ub4r,ub4i,rppu1r,rppu1i,
     &rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,b11r,b11i,b12r,b12i,
     &b21r,b21i,b22r,b22i,nfreq)
c
c     addition von ua und ub*rufr - c1
c
      call matadd(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,b11r,b11i,b12r
     &,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,
     &nfreq)
c
c     multiplikation von (ua+ub*rufr) und a1 - x1
c
      call mmtml(c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,a11r,a11i,
     &a12r,a12i,a21r,a21i,a22r,a22i,x11r,x11i,x12r,x12i,x21r,x21i,x22r,
     &x22i,nfreq)
c
c     multiplikation von x1 mit rdsl  fuer downgoing wave-field
c                            - y1
c
      call mmtml(x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,rppd2r,rppd2i,
     &rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,y11r,y11i,y12r,y12i,
     &y21r,y21i,y22r,y22i,nfreq)
      go to 1400
 1330 if(xx1.eq.0.) then
c
c     ber. von (ua+ub*rurf)*a1 = y1
c
      call mmtml(ub1r,ub1i,ub2r,ub2i,ub3r,ub3i,ub4r,ub4i,rppu1r,rppu1i,
     &rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,b11r,b11i,b12r,b12i,
     &b21r,b21i,b22r,b22i,nfreq)
      call matadd(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,b11r,b11i,b12r
     &,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,
     &nfreq)
      call mmtml(c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,a11r,a11i,
     &a12r,a12i,a21r,a21i,a22r,a22i,y11r,y11i,y12r,y12i,y21r,y21i,y22r,
     &y22i,nfreq)
      go to 1400
      else
c
c     multiplikation von ua mit a1 -- x1
c
      call mmtml(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,a11r,a11i,a12r
     &,a12i,a21r,a21i,a22r,a22i,x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,
     &nfreq)
c
c     multiplikation von ua*a1 mit rdsl -- y1
c
      call mmtml(x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,rppd2r,rppd2i,
     &rspd2r,rspd2i,rpsd2r,rpsd2i,rssd2r,rssd2i,y11r,y11i,y12r,y12i,
     &y21r,y21i,y22r,y22i,nfreq)
      if(iss(16).eq.1.and.mdeck.ge.iso) xx1=0.
      go to 1400
      end if

c
c     ber. von (ub+ua*rufr)*a1 = y1
c
 1340 call mmtml(ua1r,ua1i,ua2r,ua2i,ua3r,ua3i,ua4r,ua4i,rppu1r,rppu1i,
     &rspu1r,rspu1i,rpsu1r,rpsu1i,rssu1r,rssu1i,b11r,b11i,b12r,b12i,
     &b21r,b21i,b22r,b22i,nfreq)
      if(iss(16).eq.1.and.mdeck.ge.ire) then
         call mmtml(b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,a11r,a11i,
     &   a12r,a12i,a21r,a21i,a22r,a22i,y11r,y11i,y12r,y12i,y21r,y21i,
     &   y22r,y22i,nfreq)
      else
         call matadd(ub1r,ub1i,ub2r,ub2i,ub3r,ub3i,ub4r,ub4i,b11r,b11i,
     &   b12r,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,c21r,c21i,
     &   c22r,c22i,nfreq)
         call mmtml(c11r,c11i,c12r,c12i,c21r,c21i,c22r,c22i,a11r,a11i,
     &   a12r,a12i,a21r,a21i,a22r,a22i,y11r,y11i,y12r,y12i,y21r,y21i,
     &   y22r,y22i,nfreq)
      end if
c
c     multiplikation mit quellen-amplituden - r1,r2,r3
c
 1400 if(iss(17).eq.1) go to 1950
      if(xx1.eq.0.) then
      call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s11dr,s11di,
     &s12dr,s12di,r11r,r11i,r12r,r12i,nfreq)
      call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s21dr,s21di,
     &s22dr,s22di,r21r,r21i,r22r,r22i,nfreq)
      else
      call matmus(x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,s11ur,s11ui,
     &s12ur,s12ui,a11r,a11i,a12r,a12i,nfreq)
      call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s11dr,s11di,
     &s12dr,s12di,b11r,b11i,b12r,b12i,nfreq)
      call matadt(a11r,a11i,a12r,a12i,b11r,b11i,b12r,b12i,r11r,r11i,r12r
     &,r12i,nfreq)
      call matmus(x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,s21ur,s21ui,
     &s22ur,s22ui,a11r,a11i,a12r,a12i,nfreq)
      call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s21dr,s21di,
     &s22dr,s22di,b11r,b11i,b12r,b12i,nfreq)
      call matadt(a11r,a11i,a12r,a12i,b11r,b11i,b12r,b12i,r21r,r21i,r22r
     &,r22i,nfreq)

      end if

c
c     ber. von kappa1....kappa2  fuer explosionsquelle-entfernungsunab-
c     haengig
c
      kappa1=-window(i)
      kappa2=kappa1
c
c     multiplikation mit kappa1 - kappa2
c
      do 1500 j=1,nfreq
      r11r(j)=r11r(j)*kappa1
      r11i(j)=r11i(j)*kappa1
      r12r(j)=r12r(j)*kappa1
      r12i(j)=r12i(j)*kappa1
      r21r(j)=r21r(j)*kappa2
      r21i(j)=r21i(j)*kappa2
      r22r(j)=r22r(j)*kappa2
      r22i(j)=r22i(j)*kappa2
 1500 continue
c
c     addition von r1*kappa1 und r2*kappa2
c
      call matadt(r11r,r11i,r12r,r12i,r21r,r21i,r22r,r22i,z11r,z11i,
     &z12r,z12i,nfreq)
c
c     beginn entfernungsschleife
c
 1950 do 3000 ke=1,nent
      if(i.eq.1) then
         do 2090 j=1,nfreq
         xsum1r(j,ke)=0.
         xsum1i(j,ke)=0.
         xsum2r(j,ke)=0.
 2090    xsum2i(j,ke)=0.
      end if
      if(iss(17).eq.1) then
c
c     ber. von kappa1....kappa3
c
         phi=az(ke)/57.29578
         cphi=cos(phi)
         sphi=sin(phi)
         kappa1=0.5*(m11*cphi*cphi+m22*sphi*sphi+m12*sin(2.*phi))
         kappa2=0.5*m33
         kappa3=m13*cphi+m23*sphi
         xr=x(ke)*cphi
         yr=x(ke)*sphi
c     schleife ueber anzahl der punktquellen-ber.d.staerke und d.lage
         do 2002 js=1,nh
         dx=xr-xs(js)
         dy=yr-ys(js)
         rr(js)=sqrt(dx*dx+dy*dy)
         dz(js)=zs(js)-z(iso-1)
 2002    vv(js)=es(js)/sqrt(rr(js))
c     ber. von kappai*su(bzw.sd) ohne zeitl. verschiebung ea bzw. eb
c      *1./sqrt(u) * window(i)
         sp=window(i)/sqrt(p(i))
         p3(i)=p2(i)*p(i)
         al(iso)=cmplx(alpr(iso),alpi(iso))
         be(iso)=cmplx(betr(iso),beti(iso))
         am=csqrt(1./(al(iso)*al(iso))-p2(i))
         if(real(am).eq.0.) am=-am
         amr=real(am)
         ami=aimag(am)
         amb=1./(amr*amr+ami*ami)
         bm=csqrt(1./(be(iso)*be(iso))-p2(i))
         if(real(bm).eq.0.) bm=-bm
         bmiso2=be(iso)*be(iso)
         bmr=real(bm)
         bmi=aimag(bm)
         s1ar=p3(i)*amb*amr*sp
         s1ai=-p3(i)*amb*ami*sp
         s1br=p2(i)*sp
c        s1br=-p2(i)*sp
         s1bi=0.
         s2ar=p(i)*amr*sp
         s2ai=p(i)*ami*sp
         s2br=-s1br
         s2bi=0.
         s3ar=-p2(i)*sp
         s3ai=0.
         bmiso=p(i)/(2.*bm)*(2.*p2(i)-1./bmiso2)*sp
c        bmiso=-p(i)/(2.*bm)*(2.*p2(i)-1./bmiso2)*sp
         s3br=real(bmiso)
         s3bi=aimag(bmiso)
c     upgoing
         xu1r=kappa1*s1ar+kappa2*s2ar+kappa3*s3ar
         xu1i=kappa1*s1ai+kappa2*s2ai+kappa3*s3ai
         xu2r=kappa1*s1br+kappa2*s2br+kappa3*s3br
         xu2i=kappa1*s1bi+kappa2*s2bi+kappa3*s3bi
c     downgoing
         xd1r=kappa1*s1ar+kappa2*s2ar-kappa3*s3ar
         xd1i=kappa1*s1ai+kappa2*s2ai-kappa3*s3ai
         xd2r=-kappa1*s1br-kappa2*s2br+kappa3*s3br
         xd2i=-kappa1*s1bi-kappa2*s2bi+kappa3*s3bi
c     ber. von ea bzw. eb sowie v.bessel-funktion(fernfeldnaeherung)
         oamr=omegai*amr
         oami=omegai*ami
         obmr=omegai*bmr
         obmi=omegai*bmi
         do 2008 jf=1,nfreq
         b11r(jf)=0.
         b11i(jf)=0.
         b12r(jf)=0.
         b12i(jf)=0.
         b21r(jf)=0.
         b21i(jf)=0.
         b22r(jf)=0.
         b22i(jf)=0.
2008     continue
         do 2010 js=1,nh
         do 2012 jf=1,nfreq
         uur=(p(i)*rr(js)+ts(js))*omega(jf)
         uui=(p(i)*rr(js)+ts(js))*omegai
         xar=(-omega(jf)*ami-oamr)*dz(js)
         xai=(omega(jf)*amr-oami)*dz(js)
         xbr=(-omega(jf)*bmi-obmr)*dz(js)
         xbi=(omega(jf)*bmr-obmi)*dz(js)
         eu1r=exp(-xar+uui)
         eu2r=exp(-xbr+uui)
         ed1r=exp(xar+uui)
         ed2r=exp(xbr+uui)
         a11r(jf)=eu1r*cos(xai+uur)*vv(js)
         a11i(jf)=-eu1r*sin(xai+uur)*vv(js)
         a12r(jf)=eu2r*cos(xbi+uur)*vv(js)
         a12i(jf)=-eu2r*sin(xbi+uur)*vv(js)
         a21r(jf)=ed1r*cos(xai-uur)*vv(js)
         a21i(jf)=ed1r*sin(xai-uur)*vv(js)
         a22r(jf)=ed2r*cos(xbi-uur)*vv(js)
 2012    a22i(jf)=ed2r*sin(xbi-uur)*vv(js)
         do 2014 jf=1,nfreq
         b11r(jf)=b11r(jf)+a11r(jf)
         b11i(jf)=b11i(jf)+a11i(jf)
         b12r(jf)=b12r(jf)+a12r(jf)
         b12i(jf)=b12i(jf)+a12i(jf)
         b21r(jf)=b21r(jf)+a21r(jf)
         b21i(jf)=b21i(jf)+a21i(jf)
         b22r(jf)=b22r(jf)+a22r(jf)
 2014    b22i(jf)=b22i(jf)+a22i(jf)
 2010    continue
c     mult. mit kappai*s1u(s1d)
         do 2016 jf=1,nfreq
         s11ur(jf)=xu1r*b11r(jf)-xu1i*b11i(jf)
         s11ui(jf)=xu1i*b11r(jf)+xu1r*b11i(jf)
         s12ur(jf)=xu2r*b12r(jf)-xu2i*b12i(jf)
         s12ui(jf)=xu2i*b12r(jf)+xu2r*b12i(jf)
         s11dr(jf)=xd1r*b21r(jf)-xd1i*b21i(jf)
         s11di(jf)=xd1i*b21r(jf)+xd1r*b21i(jf)
         s12dr(jf)=xd2r*b22r(jf)-xd2i*b22i(jf)
         s12di(jf)=xd2i*b22r(jf)+xd2r*b22i(jf)
 2016    continue
         if(iss(7).eq.1) xx1=0.
         if(xx1.eq.0.) go to 2020
        call matmus(x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,s11ur,s11ui,
     &   s12ur,s12ui,a11r,a11i,a12r,a12i,nfreq)
        call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s11dr,s11di,
     &   s12dr,s12di,b11r,b11i,b12r,b12i,nfreq)
         call matadt(a11r,a11i,a12r,a12i,b11r,b11i,b12r,b12i,z11r,z11i,
     &   z12r,z12i,nfreq)
         go to 2004
 2020   call matmus(y11r,y11i,y12r,y12i,y21r,y21i,y22r,y22i,s11dr,s11di,
     &   s12dr,s12di,z11r,z11i,z12r,z12i,nfreq)
c
c     lesen der daten, dann integration
c
 2004    if(iss(4).eq.1) go to 2050
c     horizontalkomponente
         do 2041 jf=1,nfreq
         xsum1r(jf,ke)=xsum1r(jf,ke)+z11r(jf)
         xsum1i(jf,ke)=xsum1i(jf,ke)+z11i(jf)
 2041    continue
         if(i.eq.1.or.i.eq.np) then
         do 2042 jf=1,nfreq
         xsum1r(jf,ke)=xsum1r(jf,ke)-z11r(jf)*0.5
         xsum1i(jf,ke)=xsum1i(jf,ke)-z11i(jf)*0.5
 2042    continue
         end if
 2050    if(iss(4).eq.2) go to 3000
c     vertikalkomponente
         do 2052 jf=1,nfreq
         xsum2r(jf,ke)=xsum2r(jf,ke)+z12r(jf)
         xsum2i(jf,ke)=xsum2i(jf,ke)+z12i(jf)
 2052    continue
         if(i.eq.1.or.i.eq.np) then
         do 2053 jf=1,nfreq
         xsum2r(jf,ke)=xsum2r(jf,ke)-z12r(jf)*0.5
         xsum2i(jf,ke)=xsum2i(jf,ke)-z12i(jf)*0.5
 2053    continue
         end if
         go to 3000
      end if
c
c     set parameter for the approximation of the bessel-function
c
      if(iss(4).eq.1) go to 2500
c
c     ber. der horizontalkomponente
c
      ij=1
      call bessel(ij,i,ke)
      do 2200 jf=1,nfreq
      x1r(jf)=-x3r(jf)*z11r(jf)+x3i(jf)*z11i(jf)
      x1i(jf)=-x3i(jf)*z11r(jf)-x3r(jf)*z11i(jf)
      xsum1r(jf,ke)=xsum1r(jf,ke)+x1r(jf)
      xsum1i(jf,ke)=xsum1i(jf,ke)+x1i(jf)
 2200 continue
      if(i.eq.1.or.i.eq.np) then
         do 2250 jf=1,nfreq
         xsum1r(jf,ke)=xsum1r(jf,ke)-x1r(jf)*0.5
         xsum1i(jf,ke)=xsum1i(jf,ke)-x1i(jf)*0.5
 2250    continue
      end if
 2500 if(iss(4).eq.2) go to 3000
c
c     ber. der vertikalkomponente
c
      ij=0
      call bessel(ij,i,ke)
      do 2700 jf=1,nfreq
      x1r(jf)=-x2i(jf)*z12r(jf)-x2r(jf)*z12i(jf)
      x1i(jf)=x2r(jf)*z12r(jf)-x2i(jf)*z12i(jf)
      xsum2r(jf,ke)=xsum2r(jf,ke)+x1r(jf)
      xsum2i(jf,ke)=xsum2i(jf,ke)+x1i(jf)
 2700 continue
      if(i.eq.1.or.i.eq.np) then
         do 2750 jf=1,nfreq
         xsum2r(jf,ke)=xsum2r(jf,ke)-x1r(jf)*0.5
         xsum2i(jf,ke)=xsum2i(jf,ke)-x1i(jf)*0.5
 2750    continue
      end if
c
c     multiply with wavelet, after inverse transform vertical
c     seismograms are stored on real part of s, horizontal
c     ones on the imaginary part
c
 3000 continue
 1000 continue
c
c     ende winkelschleife
c
      do 2600 ke=1,nent
      jf1=0
      t0=x(ke)/vred+tmin
      do 2091 j=1,nfreq
 2091 phase(j)=omega(j)*t0
      if(iss(17).eq.1) then
         do 2092 j=1,nfreq
 2092    phase(j)=phase(j)+pi4
      end if
      do 2650 jf=1,inpt
      jz=npts-jf+1
      s(jf)=0.
      s(jz)=0.
      if(jf.lt.iu.or.jf.gt.io) go to 2650
      jf1=jf1+1
      if(jf.eq.iu.or.jf.eq.io) go to 2650
      jz=npts-jf+1
      xsum1r(jf1,ke)=xsum1r(jf1,ke)*d11x*dp
      xsum1i(jf1,ke)=xsum1i(jf1,ke)*d11x*dp
      s(jz)=cmplx(xsum1r(jf1,ke),xsum1i(jf1,ke))*cexp(ein*phase(jf1))
      xsum2r(jf1,ke)=xsum2r(jf1,ke)*d11x*dp
      xsum2i(jf1,ke)=xsum2i(jf1,ke)*d11x*dp
      s(jf)=cmplx(xsum2r(jf1,ke),xsum2i(jf1,ke))*cexp(ein*phase(jf1))
 2650 continue
      if(iss(6).eq.1) then
         do 3001 jf=1,inpt
 3001    sc1(jf)=s(jf)*wavel(jf)
         smax=abs(real(sc1(1)))
         if(abs(aimag(sc1(1))).gt.smax) smax=abs(aimag(sc1(1)))
         do 3002 jf=2,inpt
         if(abs(real(sc1(jf))).gt.smax) smax=abs(real(sc1(jf)))
 3002    if(abs(aimag(sc1(jf))).gt.smax) smax=abs(aimag(sc1(jf)))
         do 3003 jf=1,inpt
 3003    sc1(jf)=sc1(jf)/smax
         write(21,5004) fu,fwil,fwir,fo,dnue,c2,cwir,cwil,c1,np
         write(21,5005) x(ke),npts,smax
         write(21,5006) (sc1(jf),jf=1,inpt)
         inpts=inpt+1
         do 3004 jf=npts,inpts,-1
 3004    sc1(jf)=s(jf)*wavel(npts-jf+1)
         smax=abs(real(sc1(npts)))
         if(abs(aimag(sc1(npts))).gt.smax) smax=abs(aimag(sc1(npts)))
         do 3005 jf=npts-1,inpts,-1
         if(abs(real(sc1(jf))).gt.smax) smax=abs(real(sc1(jf)))
 3005    if(abs(aimag(sc1(jf))).gt.smax) smax=abs(aimag(sc1(jf)))
         do 3006 jf=npts,inpts,-1
 3006    sc1(jf)=sc1(jf)/smax
         write(22,5004) fu,fwil,fwir,fo,dnue,c2,cwir,cwil,c1,np
         write(22,5005) x(ke),npts,smax
         write(22,5006) (sc1(jf),jf=npts,inpts,-1)
      end if

      s(1)=(s(1)+ein*s(npts))*wavel(1)
      do 3100 kf=2,inpt
      kf1 = npts - kf + 1
      kf2 = kf1 + 1
      s(kf2)=(conjg(s(kf))+ein*conjg(s(kf1)))*conjg(wavel(kf))
 3100 s(kf) = (s(kf) + ein*s(kf1))*wavel(kf)
      inp1=inpt+1
      s(inp1)=(0.0,0.0)
      call cool(mga,s,1.)
      ddd=1.
      ccc=1.

      CALL FLUSH()
      if(tsigma.ne.0.) then
         ddd=exp(t0/tsigma)
         ccc=exp(dt/tsigma)
      end if
      dd=1./ccc
      do 3105 i=1,npts
      dd=dd*ccc
 3105 s(i)=s(i)*dd*ddd
c
c     write out seismograms on units 6 and iouts
c

      if(iss(4) - 1) 3200,3300,3400
 3200 l1 = 1
      l2 = 2
      go to 3500
 3300 l1 = 1
      l2 = 1
      go to 3500
 3400 l1 = 2
      l2 = 2
 3500 sphred=1.
      if(iss(12).eq.1)  sphred=sqrt((x(ke)/3390.)/abs(sin(x(ke)/3390.)))
     1*((3390.-dso)/3390.)**(float(nrho+1)/2.)
      do 3600 kom = l1,l2
      go to(3610,3630),kom
 3610 do 3620 j=1,npts
 3620 balla(j)=sphred*real(s(j))
      go to 3650
 3630 do 3640 j=1,npts
 3640 balla(j)=sphred*aimag(s(j))
 3650 balmax = amax(balla,npts)
      bal = abs(amin(balla,npts))
      if(bal.gt.balmax) balmax = bal

c
c create sac file, changed 19.11.2009 from H. Paulsen
c B. Tauzin, Universiteit Utrecht
c
      call newhdr
      call setfhv('B', t0, nerr)
      call setfhv('DELTA', dt, nerr)
      call setfhv('EVDP',ztrue(1),nerr)
      call setnhv('NPTS', npts, nerr)
      call setlhv('LEVEN', .true., nerr)
      call setfhv('DIST', x(ke), nerr)
      call setfhv('GCARC', x(ke)/59.16, nerr)
      call setfhv('AZ', az(ke), nerr)
      call base(filename,ke)
      if(kom.eq.1) filename(6:7)=".z"
      if(kom.eq.2) filename(6:7)=".r"
      call wsac0(filename,balla(1),balla(1),nerr)
c
      do 3700 i = 1,npts
 3700 iballa(i) = ifix(999.10*balla(i)/balmax)
c

      if(iss(15).eq.1) go to 3800
      epi=x(ke)*0.0169
      go to(3901,3902),kom
 3901 write(6,4001) x(ke),epi
 4001 format(///,1x,'seismogram (vertical displacement) for the distance
     1',f10.3,' km   =',f8.3,' degree'/)
      go to 3910
 3902 write(6,4002) x(ke),epi
 4002 format(///,1x,'seismogram (horizontal displacement) for the distan
     1ce',f10.3,' km   =',f8.3,' degree'/)
 3910 write(6,4003) npts,balmax
 4003 format(i10,e15.4)
      write(6,4004) (iballa(i), i = 1,npts)
 4004 format(20i4)
 3800 do 3920 i=1,npts
      balla(i) = 4.99*balla(i)/balmax + 5.001
 3920 iballa(i) = 10000.*balla(i)
      abal = tmin - na*dt + x(ke)/vred
cc    write(iouts,5001)x(ke),abal,kom,npts,balmax
      write(6,4005) x(ke),abal,kom,npts,balmax
 5001 format(2f15.5,i5/,i10,5x,e15.4)
 4005 format(2f10.4,i5,i10,5x,e15.4)
cc    write(iouts,5002)(iballa(i),i=1,npts)
 5002 format(16i5)
cc    write(iouts,5003)x(ke)
 5003 format(f15.5)
 5004 format(/' frequenzband:'/4f10.4/' frequenzintervall:'/f10.6/
     &' strahlparameterfenster:'/4f10.4/' anzahl der strahlparameter:'/
     &i10)
 5005 format(/f10.4,i10,e15.5)
 5006 format(10f8.5)
 3600 continue
 2600 continue
      go to 6602
c
c     calculation of seismograms is finished
c
 6591 nfrequ=nfr
      write(6,9301) nfrequ,nfreq
 9301 format(//1x,'the maximum number of frequencies is',i10/1x,'the cal
     &culated number nfreq = dt*npts*(fo-fu)'/29x,'=',i10)
 6601 write(6,9302)
 9302 format(/1x,'the job was cancelled because of the input error menti
     &oned above')
 6602 stop
      end
c==================================================================
      subroutine bessel(ij,i,ke)
c==================================================================
c
c        the subroutine bessel computes the bessel functions of the 0.
c                    and 1.order for complex arguments
c
c      parameter(nfr=1024,npa=2500,nla=1000,ndis=200)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000,ndis=200)
      common xa(nfr),xas(nfr),xj0r(nfr),xj0i(nfr),xj1r(nfr),xj1i(nfr),
     &arg2(nfr),ersin(nfr),x5r(nfr),x5i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),rrr(nfr),rri(nfr),rsr(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr),xsum1r(nfr,ndis),
     &xsum1i(nfr,ndis),xsum2r(nfr,ndis),xsum2i(nfr,ndis)
      common/blo1/z(0:nla),d(nla),alpr(nla),alpi(nla),betr(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      data sc/0.39894228/
      complex cv
      x1=x(ke)*p(i)
      do 100 jf=1,nfreq
      xa(jf)=x1*omegac(jf)
      arg2(jf)=xa(jf)*cphasc(jf)+phasc(jf)
      ersin(jf)=exp(xa(jf)*sphasc(jf))
 100  continue
      if(ij.eq.1) go to 400
c
c   the argument for the bessel function of the 0.order is calculated
c
      if(iss(17).eq.4) then
c
c               line source
c
         do 244 jf=1,nfreq
         xj0r(jf)=cos(arg2(jf))
 244     xj0i(jf)=0.
      else
c
c               point source
c
          do 240 jf=1,nfreq
          xas(jf)=1./sqrt(xa(jf))*sc
          xj0r(jf)=xas(jf)*cos(arg2(jf))*ersin(jf)
 240      xj0i(jf)=-xas(jf)*sin(arg2(jf))*ersin(jf)
      end if
 411  return
c
c   the argument for the bessel function of the 1.order is calculated
c
 400  if(iss(17).eq.4) then
c
c               line source
c
         do 344 jf=1,nfreq
         xj1r(jf)=sin(arg2(jf))
 344     xj1i(jf)=0.
      else
c
c               point source
c
          do 340 jf=1,nfreq
          xas(jf)=1./sqrt(xa(jf))*sc
          xj1r(jf)=xas(jf)*sin(arg2(jf))*ersin(jf)
 340      xj1i(jf)=xas(jf)*cos(arg2(jf))*ersin(jf)
      end if
 511  return
      end
c==================================================================
      subroutine rekurs(n1,n2,i1,alpr,alpi,betr,beti,rho,d,p,rppr,rpp
c==================================================================
     &i,rspr,rspi,rssr,rssi,rpsr,rpsi,it)
c
c     rekurs computes the plane wave reflection coefficients for one
c     slowness (p) and all required frequencies for n layers with a
c             recursive algorthm described by kennett:
c        theoretical reflection seismograms for elastic media,
c                    geoph. prosp.,v.27,p.301-321,1979
c
c     by complex velocities damping is considered
c
c     i1 =  1  -  reflecitivity for an incident upgoing wavefield
c        = -1  -  reflecitivity for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer, if n2=0 - free surface
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch oberste bzw. unterste schicht
c        = 1   -  transmission durch oberste bzw. unterste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex c1,al(nla),be(nla),aq1,bq1,cosi,cosj
      dimension alpr(nla),betr(nla),alpi(nla),beti(nla),rho(nla),
     &d(nla)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/refk/rppdr,rppdi,rspdr,rspdi,rpsdr,rpsdi,rssdr,rssdi,
     &rppur,rppui,rspur,rspui,rpsur,rpsui,rssur,rssui/transk/tppdr,
     &tppdi,tspdr,tspdi,tpsdr,tpsdi,tssdr,tssdi,tppur,tppui,tspur,
     &tspui,tpsur,tpsui,tssur,tssui
      common x1r(nfr),x1i(nfr),x2r(nfr),x2i(nfr),x3r(nfr),x3i(nfr),
     &x4r(nfr),x4i(nfr),x5r(nfr),x5i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),r1r(nfr),r1i(nfr),r2r(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr)
      dimension rppr(nfr),rppi(nfr),rspr(nfr),rspi(nfr),rssr(nfr),rssi(n
     &fr),rpsr(nfr),rpsi(nfr)
c
c     computation of the reflection and transmission coefficients at
c                             1 interface
c
      if(n1.eq.n2) then
          call set(rppr,rppi,rspr,rspi,rpsr,rpsi,rssr,rssi,0.,nfreq)
          return
      end if
      if(n2.ne.0) then
        n3=n2+i1
        al(n2)=cmplx(alpr(n2),alpi(n2))
        be(n2)=cmplx(betr(n2),beti(n2))
        al(n3)=cmplx(alpr(n3),alpi(n3))
        be(n3)=cmplx(betr(n3),beti(n3))
        call reflex(p,rho(n3),rho(n2),al(n3),al(n2),be(n3),be(n2))
      else
        n3=n2+1
        al(n3)=cmplx(alpr(n3),alpi(n3))
        be(n3)=cmplx(betr(n3),beti(n3))
        call refls(p,al(n3),be(n3))
      end if
      if(i1.eq.1.and.n2.ne.0) then
         rpsdr=-rpsdr
         rpsdi=-rpsdi
         rspdr=-rspdr
         rspdi=-rspdi
         rpsur=-rpsur
         rpsui=-rpsui
         rspur=-rspur
         rspui=-rspui
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
         tpsur=-tpsur
         tpsui=-tpsui
         tspur=-tspur
         tspui=-tspui
      end if
      do 217 i=1,nfreq
      r1r(i)=rppdr
      r1i(i)=rppdi
      r2r(i)=rspdr
      r2i(i)=rspdi
      r3r(i)=rssdr
 217  r3i(i)=rssdi
      nb=n2+i1
      ne=n1-i1
      nend=nb-ne
      if(i1.gt.0) then
         nend=-nend
      end if
      pi=3.14159265*2.
      if(nend.lt.0) go to 11
c
c     begin of the loop over the layers
c     the values which are independent of the frequency are calculated
c
      do 10 j=nb,ne,i1
      al(j)=cmplx(alpr(j),alpi(j))
      al(j+i1)=cmplx(alpr(j+i1),alpi(j+i1))
      be(j)=cmplx(betr(j),beti(j))
      be(j+i1)=cmplx(betr(j+i1),beti(j+i1))
      aq1=csqrt(1./(al(j)*al(j))-p*p)
      if(real(aq1).eq.0.) aq1=-aq1
      bq1=csqrt(1./(be(j)*be(j))-p*p)
      if(real(bq1).eq.0.) bq1=-bq1
      c1=aq1/bq1
      c1r=real(c1)
      c1i=aimag(c1)
      call reflex(p,rho(j+i1),rho(j),al(j+i1),al(j),be(j+i1),be(j))
      if(i1.eq.1) then
         rpsdr=-rpsdr
         rpsdi=-rpsdi
         rspdr=-rspdr
         rspdi=-rspdi
         rpsur=-rpsur
         rpsui=-rpsui
         rspur=-rspur
         rspui=-rspui
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
         tpsur=-tpsur
         tpsui=-tpsui
         tspur=-tspur
         tspui=-tspui
      end if
      pp1r=tppur*tppdr-tppui*tppdi
      pp1i=tppui*tppdr+tppur*tppdi
      pp2r=tspur*tpsdr-tspui*tpsdi
      pp2i=tspui*tpsdr+tspur*tpsdi
      pp3r=tppur*tpsdr-tppui*tpsdi
      pp3i=tppui*tpsdr+tppur*tpsdi
      pp4r=tspur*tppdr-tspui*tppdi
      pp4i=tspui*tppdr+tspur*tppdi
      pp5r=pp1r*rssur-pp1i*rssui
      pp5i=pp1i*rssur+pp1r*rssui
      pp6r=pp2r*rppur-pp2i*rppui
      pp6i=pp2i*rppur+pp2r*rppui
      pp7r=pp3r*rspur-pp3i*rspui
      pp7i=pp3i*rspur+pp3r*rspui
      pp8r=pp4r*rpsur-pp4i*rpsui
      pp8i=pp4i*rpsur+pp4r*rpsui
      pp3r=-pp4r*c1r+pp4i*c1i+pp3r
      pp3i=-pp4i*c1r-pp4r*c1i+pp3i
      pp4r=-pp5r-pp6r+pp7r+pp8r
      pp4i=-pp5i-pp6i+pp7i+pp8i
      pp5r=pp4r*c1r-pp4i*c1i
      pp5i=pp4i*c1r+pp4r*c1i
      sp1r=tppur*tspdr-tppui*tspdi
      sp1i=tppui*tspdr+tppur*tspdi
      sp2r=tspur*tssdr-tspui*tssdi
      sp2i=tspui*tssdr+tspur*tssdi
      sp3r=tppur*tssdr-tppui*tssdi
      sp3i=tppui*tssdr+tppur*tssdi
      sp4r=tspur*tspdr-tspui*tspdi
      sp4i=tspui*tspdr+tspur*tspdi
      sp5r=sp1r*rssur-sp1i*rssui
      sp5i=sp1i*rssur+sp1r*rssui
      sp6r=sp2r*rppur-sp2i*rppui
      sp6i=sp2i*rppur+sp2r*rppui
      sp7r=sp3r*rspur-sp3i*rspui
      sp7i=sp3i*rspur+sp3r*rspui
      sp8r=sp4r*rpsur-sp4i*rpsui
      sp8i=sp4i*rpsur+sp4r*rpsui
      sp3r=-sp4r*c1r+sp4i*c1i+sp3r
      sp3i=-sp4i*c1r-sp4r*c1i+sp3i
      sp4r=-sp5r-sp6r+sp7r+sp8r
      sp4i=-sp5i-sp6i+sp7i+sp8i
      sp5r=sp4r*c1r-sp4i*c1i
      sp5i=sp4i*c1r+sp4r*c1i
      ss1r=tpsur*tspdr-tpsui*tspdi
      ss1i=tpsui*tspdr+tpsur*tspdi
      ss2r=tssur*tssdr-tssui*tssdi
      ss2i=tssui*tssdr+tssur*tssdi
      ss3r=tpsur*tssdr-tpsui*tssdi
      ss3i=tpsui*tssdr+tpsur*tssdi
      ss4r=tssur*tspdr-tssui*tspdi
      ss4i=tssui*tspdr+tssur*tspdi
      ss5r=ss1r*rssur-ss1i*rssui
      ss5i=ss1i*rssur+ss1r*rssui
      ss6r=ss2r*rppur-ss2i*rppui
      ss6i=ss2i*rppur+ss2r*rppui
      ss7r=ss3r*rspur-ss3i*rspui
      ss7i=ss3i*rspur+ss3r*rspui
      ss8r=ss4r*rpsur-ss4i*rpsui
      ss8i=ss4i*rpsur+ss4r*rpsui
      ss3r=-ss4r*c1r+ss4i*c1i+ss3r
      ss3i=-ss4i*c1r-ss4r*c1i+ss3i
      ss4r=-ss5r-ss6r+ss7r+ss8r
      ss4i=-ss5i-ss6i+ss7i+ss8i
      ss5r=ss4r*c1r-ss4i*c1i
      ss5i=ss4i*c1r+ss4r*c1i
      det1r=-rppur
      det1i=-rppui
      det2r=-rssur
      det2i=-rssui
      det3r=rspur*c1r-rspui*c1i-rpsur
      det3i=rspui*c1r+rspur*c1i-rpsui
      det4r=rssur*rppur-rssui*rppui-rspur*rpsur+rspui*rpsui
      det4i=rssui*rppur+rssur*rppui-rspui*rpsur-rspur*rpsui
      det5r=det4r*c1r-det4i*c1i
      det5i=det4i*c1r+det4r*c1i
c
c     begin of the vectorizable loop over the frequencies
c
c     cosi=csqrt(1./(al(j)*al(j))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(j)*be(j))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     tpr=-d(j)*cosir
c     tpi=-d(j)*cosii
c     tsr=-d(j)*cosjr
c     tsi=-d(j)*cosji
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(j)
      if(omegar.eq.0.) then
         cosi=csqrt(1./(al(j)*al(j))-p*p)
         cosj=csqrt(1./(be(j)*be(j))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,j)
      end if
      do 216 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      z11r(i)=z1r(i)*z1r(i)-z1i(i)*z1i(i)
      z11i(i)=2.*z1r(i)*z1i(i)
      z22r(i)=z2r(i)*z2r(i)-z2i(i)*z2i(i)
      z22i(i)=2.*z2r(i)*z2i(i)
      z12r(i)=z1r(i)*z2r(i)-z1i(i)*z2i(i)
      z12i(i)=z1i(i)*z2r(i)+z1r(i)*z2i(i)
      x1r(i)=z11r(i)*r1r(i)-z11i(i)*r1i(i)
      x1i(i)=z11i(i)*r1r(i)+z11r(i)*r1i(i)
      x2r(i)=z22r(i)*r3r(i)-z22i(i)*r3i(i)
      x2i(i)=z22i(i)*r3r(i)+z22r(i)*r3i(i)
      x3r(i)=z12r(i)*r2r(i)-z12i(i)*r2i(i)
      x3i(i)=z12i(i)*r2r(i)+z12r(i)*r2i(i)
      x4r(i)=x1r(i)*x2r(i)-x1i(i)*x2i(i)
      x4i(i)=x1i(i)*x2r(i)+x1r(i)*x2i(i)
      x5r(i)=x3r(i)*x3r(i)-x3i(i)*x3i(i)
      x5i(i)=2.*x3r(i)*x3i(i)
      detr(i)=1.+det1r*x1r(i)-det1i*x1i(i)+det2r*x2r(i)-det2i*x2i(i)+
     &det3r*x3r(i)-det3i*x3i(i)+det4r*x4r(i)-det4i*x4i(i)+det5r*x5r(i)-
     &det5i*x5i(i)
      deti(i)=det1i*x1r(i)+det1r*x1i(i)+det2i*x2r(i)+det2r*x2i(i)+
     &det3i*x3r(i)+det3r*x3i(i)+det4i*x4r(i)+det4r*x4i(i)+det5i*x5r(i)+
     &det5r*x5i(i)
      dett(i)=detr(i)*detr(i)+deti(i)*deti(i)
      dett(i)=1./dett(i)
      w1r(i)=pp1r*x1r(i)-pp1i*x1i(i)+pp2r*x2r(i)-pp2i*x2i(i)+pp3r*x3r(i)
     &-pp3i*x3i(i)+pp4r*x4r(i)-pp4i*x4i(i)+pp5r*x5r(i)-pp5i*x5i(i)
      w1i(i)=pp1i*x1r(i)+pp1r*x1i(i)+pp2i*x2r(i)+pp2r*x2i(i)+pp3i*x3r(i)
     &+pp3r*x3i(i)+pp4i*x4r(i)+pp4r*x4i(i)+pp5i*x5r(i)+pp5r*x5i(i)
      w2r(i)=sp1r*x1r(i)-sp1i*x1i(i)+sp2r*x2r(i)-sp2i*x2i(i)+sp3r*x3r(i)
     &-sp3i*x3i(i)+sp4r*x4r(i)-sp4i*x4i(i)+sp5r*x5r(i)-sp5i*x5i(i)
      w2i(i)=sp1i*x1r(i)+sp1r*x1i(i)+sp2i*x2r(i)+sp2r*x2i(i)+sp3i*x3r(i)
     &+sp3r*x3i(i)+sp4i*x4r(i)+sp4r*x4i(i)+sp5i*x5r(i)+sp5r*x5i(i)
      w3r(i)=ss1r*x1r(i)-ss1i*x1i(i)+ss2r*x2r(i)-ss2i*x2i(i)+ss3r*x3r(i)
     &-ss3i*x3i(i)+ss4r*x4r(i)-ss4i*x4i(i)+ss5r*x5r(i)-ss5i*x5i(i)
      w3i(i)=ss1i*x1r(i)+ss1r*x1i(i)+ss2i*x2r(i)+ss2r*x2i(i)+ss3i*x3r(i)
     &+ss3r*x3i(i)+ss4i*x4r(i)+ss4r*x4i(i)+ss5i*x5r(i)+ss5r*x5i(i)
      r1r(i)=rppdr+(w1r(i)*detr(i)+w1i(i)*deti(i))*dett(i)
      r1i(i)=rppdi+(w1i(i)*detr(i)-w1r(i)*deti(i))*dett(i)
      r2r(i)=rspdr+(w2r(i)*detr(i)+w2i(i)*deti(i))*dett(i)
      r2i(i)=rspdi+(w2i(i)*detr(i)-w2r(i)*deti(i))*dett(i)
      r3r(i)=rssdr+(w3r(i)*detr(i)+w3i(i)*deti(i))*dett(i)
      r3i(i)=rssdi+(w3i(i)*detr(i)-w3r(i)*deti(i))*dett(i)
216   continue
10    continue
11    al(n1)=cmplx(alpr(n1),alpi(n1))
      be(n1)=cmplx(betr(n1),beti(n1))
      aq1=csqrt(1./(al(n1)*al(n1))-p*p)
      if(real(aq1).eq.0.) aq1=-aq1
      bq1=csqrt(1./(be(n1)*be(n1))-p*p)
      if(real(bq1).eq.0.) bq1=-bq1
      c1=aq1/bq1
      c1r=real(c1)
      c1i=aimag(c1)
      if(it.eq.1) then
c
c     transmission durch oberste bzw. unterste schicht
c
c     cosi=csqrt(1./(al(n1)*al(n1))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(n1)*be(n1))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     tpr=-d(n1)*cosir
c     tpi=-d(n1)*cosii
c     tsr=-d(n1)*cosjr
c     tsi=-d(n1)*cosji
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(n1)
      if(omegar.eq.0.) then
         cosi=csqrt(1./(al(n1)*al(n1))-p*p)
         cosj=csqrt(1./(be(n1)*be(n1))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,n1)
      end if
      do 219 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      z11r(i)=z1r(i)*z1r(i)-z1i(i)*z1i(i)
      z11i(i)=2.*z1r(i)*z1i(i)
      z22r(i)=z2r(i)*z2r(i)-z2i(i)*z2i(i)
      z22i(i)=2.*z2r(i)*z2i(i)
      z12r(i)=z1r(i)*z2r(i)-z1i(i)*z2i(i)
      z12i(i)=z1i(i)*z2r(i)+z1r(i)*z2i(i)
      rppr(i)=z11r(i)*r1r(i)-z11i(i)*r1i(i)
      rppi(i)=z11i(i)*r1r(i)+z11r(i)*r1i(i)
      rspr(i)=z12r(i)*r2r(i)-z12i(i)*r2i(i)
      rspi(i)=z12i(i)*r2r(i)+z12r(i)*r2i(i)
      rssr(i)=z22r(i)*r3r(i)-z22i(i)*r3i(i)
      rssi(i)=z22i(i)*r3r(i)+z22r(i)*r3i(i)
      rpsr(i)=-rspr(i)*c1r+rspi(i)*c1i
      rpsi(i)=-rspi(i)*c1r-rspr(i)*c1i
219   continue
      else
      do 220 i=1,nfreq
      rppr(i)=r1r(i)
      rppi(i)=r1i(i)
      rspr(i)=r2r(i)
      rspi(i)=r2i(i)
      rssr(i)=r3r(i)
      rssi(i)=r3i(i)
      rpsr(i)=-rspr(i)*c1r+rspi(i)*c1i
 220  rpsi(i)=-rspi(i)*c1r-rspr(i)*c1i
      end if
      if(iss(1).eq.1) then
         do 230 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
         rspi(i)=0.
         rssr(i)=0.
 230     rssi(i)=0.
      else if(iss(1).eq.2) then
         do 240 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
         rspi(i)=0.
         rppr(i)=0.
 240     rppi(i)=0.
      else if(iss(1).eq.3) then
         do 250 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
 250     rspi(i)=0.
      end if
      return
      end
c==================================================================
      subroutine transm(n1,n2,in,alpr,alpi,betr,beti,rho,d,p,tppr,tppi,
     &tspr,tspi,tssr,tssi,tpsr,tpsi,it)
c==================================================================
c
c     ber. des transmissionsvermoegens eines schichtenpaketes
c
c     by complex velocities damping is considered
c
c     if iss(2).ne.0 no converted waves are considered (tps and tsp = 0)
c
c     in =  1  -  transmission for an incident upgoing wavefield
c        = -1  -  transmission for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch unterste und oberste schicht
c        = 1   -  transmission durch unterste bzw. oberste schicht
c        = 2   -  transmission durch unterste und oberste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex al(nla),be(nla),cosi,cosj
      dimension d(nla),alpr(nla),alpi(nla),betr(nla),beti(nla),
     &rho(nla)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/transc/tppdr,tppdi,tspdr,tspdi,tpsdr,tpsdi,tssdr,tssdi
      common x1r(nfr),x1i(nfr),x2r(nfr),x2i(nfr),x3r(nfr),x3i(nfr),
     &x4r(nfr),x4i(nfr),x5r(nfr),x5i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),r1r(nfr),r1i(nfr),r2r(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr)
      dimension tppr(nfr),tppi(nfr),tspr(nfr),tspi(nfr),tssr(nfr),tssi(n
     &fr),tpsr(nfr),tpsi(nfr),x11r(nfr),x11i(nfr),x12r(nfr),x12i(nfr),
     &x21r(nfr),x21i(nfr),x22r(nfr),x22i(nfr),xppr(nfr),xppi(nfr),xspr(n
     &fr),xspi(nfr),xpsr(nfr),xpsi(nfr),xssr(nfr),xssi(nfr)
      i1=-in
      n11=n1
      n22=n2
      tpr=0.
      tpi=0.
      tsr=0.
      tsi=0.
      yppr=1.
      yppi=0.
      yssr=1.
      yssi=0.
      if(iss(2).eq.0) then
         call set(tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,1.,nfreq)
         call set(xppr,xppi,xspr,xspi,xpsr,xpsi,xssr,xssi,1.,nfreq)
      end if
      if(n11.eq.n22.and.it.ne.2) then
         call set(tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,1.,nfreq)
         return
      end if
      if(in.eq.1) then
         nx=n11
         ny=n22+1
         ni1=0
         ni2=1
      else
         nx=n22
         ny=n11
         ni1=1
         ni2=0
      end if
      if(it-1) 2,1,3
c
c     transmission through the lowermost layer
c
 3    al(nx)=cmplx(alpr(nx),alpi(nx))
      be(nx)=cmplx(betr(nx),beti(nx))
      cosi=csqrt(1./(al(nx)*al(nx))-p*p)
      cosir=real(cosi)
      cosii=aimag(cosi)
      cosii=abs(cosii)
      cosj=csqrt(1./(be(nx)*be(nx))-p*p)
      cosjr=real(cosj)
      cosji=aimag(cosj)
      cosji=abs(cosji)
      tpr=-d(nx)*cosir
      tpi=-d(nx)*cosii
      tsr=-d(nx)*cosjr
      tsi=-d(nx)*cosji
      if(iss(2).ne.0.and.n11.ne.n22) go to 1
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(nx)
      if(omegar.eq.0.) then
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,nx)
      end if
      do 10 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      x11r(i)=x1r(i)*cos(x1i(i))
      x11i(i)=x1r(i)*sin(x1i(i))
      x22r(i)=x2r(i)*cos(x2i(i))
      x22i(i)=x2r(i)*sin(x2i(i))
 10   continue
      if(in.eq.1.or.n11.eq.n22) then
         do 11 i=1,nfreq
         tppr(i)=x11r(i)
         tppi(i)=x11i(i)
         tssr(i)=x22r(i)
         tssi(i)=x22i(i)
         tpsr(i)=0.
         tpsi(i)=0.
         tspr(i)=0.
         tspi(i)=0.
 11      continue
      else
         do 12 i=1,nfreq
         xppr(i)=x11r(i)
         xppi(i)=x11i(i)
         xssr(i)=x22r(i)
         xssi(i)=x22i(i)
 12      continue
      end if
      go to 1
 2    n3=ny+i1
      al(ny)=cmplx(alpr(ny),alpi(ny))
      al(n3)=cmplx(alpr(n3),alpi(n3))
      be(ny)=cmplx(betr(ny),beti(ny))
      be(n3)=cmplx(betr(n3),beti(n3))
      call tkoeff(p,rho(ny),rho(n3),al(ny),al(n3),be(ny),be(n3))
      if(in.eq.1) then
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
      end if
      if(iss(1).ne.0.or.iss(2).ne.0) then
         tspdr=0.
         tspdi=0.
         tpsdr=0.
         tpsdi=0.
      end if
      yppr=tppdr
      yppi=tppdi
      yssr=tssdr
      yssi=tssdi
      n11=n11+ni1
      n22=n22+ni2
      if(iss(2).ne.0.and.n11.ne.n22) go to 1
      if(in.eq.-1.or.n11.eq.n22) then
         do 20 i=1,nfreq
         tppr(i)=tppdr
         tppi(i)=tppdi
         tspr(i)=tspdr
         tspi(i)=tspdi
         tssr(i)=tssdr
         tssi(i)=tssdi
         tpsr(i)=tpsdr
         tpsi(i)=tpsdi
 20      continue
      else
         do 21 i=1,nfreq
         xppr(i)=tppdr
         xppi(i)=tppdi
         xspr(i)=tspdr
         xspi(i)=tspdi
         xssr(i)=tssdr
         xssi(i)=tssdi
         xpsr(i)=tpsdr
         xpsi(i)=tpsdi
 21      continue
      end if
 1    if(n11.eq.n22) go to 2000
      if(iss(2).eq.0) go to 1000
c
c     no converted waves are considered (tps and tsp = 0)
c
      ne=n22+in
      do 100 ns=n11,ne,i1
      if(in.eq.1) then
         nx=ns+i1
      else
         nx=ns
      end if
      al(ns)=cmplx(alpr(ns),alpi(ns))
      be(ns)=cmplx(betr(ns),beti(ns))
      n3=ns+i1
      al(n3)=cmplx(alpr(n3),alpi(n3))
      be(n3)=cmplx(betr(n3),beti(n3))
      call tkoeff(p,rho(ns),rho(n3),al(ns),al(n3),be(ns),be(n3))
      w11r=tppdr*yppr-tppdi*yppi
      w11i=tppdr*yppi+tppdi*yppr
      w22r=tssdr*yssr-tssdi*yssi
      w22i=tssdr*yssi+tssdi*yssr
      yppr=w11r
      yppi=w11i
      yssr=w22r
      yssi=w22i
      al(nx)=cmplx(alpr(nx),alpi(nx))
      be(nx)=cmplx(betr(nx),beti(nx))
      cosi=csqrt(1./(al(nx)*al(nx))-p*p)
c     if(real(cosi).eq.0.) cosi=-cosi
      cosir=real(cosi)
      cosii=aimag(cosi)
      cosii=abs(cosii)
      cosi=csqrt(1./(be(nx)*be(nx))-p*p)
c     if(real(cosi).eq.0.) cosi=-cosi
      cosjr=real(cosi)
      cosji=aimag(cosi)
      cosji=abs(cosji)
      tpr=-d(nx)*cosir+tpr
      tpi=-d(nx)*cosii+tpi
      tsr=-d(nx)*cosjr+tsr
      tsi=-d(nx)*cosji+tsi
 100  continue
      tpri=omegai*tpr
      tpii=omegai*tpi
      tsri=omegai*tsr
      tsii=omegai*tsi
      do 200 i=1,nfreq
      x1r(i)=omega(i)*tpi-tpri
      x1i(i)=omega(i)*tpr+tpii
      x2r(i)=omega(i)*tsi-tsri
      x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      tppr(i)=yppr*z1r(i)-yppi*z1i(i)
      tppi(i)=yppi*z1r(i)+yppr*z1i(i)
      tssr(i)=yssr*z2r(i)-yssi*z2i(i)
      tssi(i)=yssi*z2r(i)+yssr*z2i(i)
      tspr(i)=0.
      tspi(i)=0.
      tpsr(i)=0.
      tpsi(i)=0.
200   continue
      go to 2000
c
c     converted waves are considered (tps and tsp ne 0)
c
1000  ne=n22+in
      do 300 ns=n11,ne,i1
      if(in.eq.1) then
         nx=ns+i1
      else
         nx=ns
      end if
      al(ns)=cmplx(alpr(ns),alpi(ns))
      be(ns)=cmplx(betr(ns),beti(ns))
      n3=ns+i1
      al(n3)=cmplx(alpr(n3),alpi(n3))
      be(n3)=cmplx(betr(n3),beti(n3))
      call tkoeff(p,rho(ns),rho(n3),al(ns),al(n3),be(ns),be(n3))
      if(in.eq.1) then
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
      end if
c     al(nx)=cmplx(alpr(nx),alpi(nx))
c     be(nx)=cmplx(betr(nx),beti(nx))
c     cosi=csqrt(1./(al(nx)*al(nx))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(nx)*be(nx))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     tpr=-d(nx)*cosir
c     tpi=-d(nx)*cosii
c     tsr=-d(nx)*cosjr
c     tsi=-d(nx)*cosji
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(nx)
      if(omegar.eq.0.) then
         al(nx)=cmplx(alpr(nx),alpi(nx))
         be(nx)=cmplx(betr(nx),beti(nx))
         cosi=csqrt(1./(al(nx)*al(nx))-p*p)
         cosj=csqrt(1./(be(nx)*be(nx))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,nx)
      end if
      do 310 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
310   z2i(i)=x2r(i)*sin(x2i(i))
      if(in.eq.1) then
         do 320 i=1,nfreq
         x1r(i)=tppdr*z1r(i)-tppdi*z1i(i)
         x1i(i)=tppdi*z1r(i)+tppdr*z1i(i)
         x2r(i)=tspdr*z1r(i)-tspdi*z1i(i)
         x2i(i)=tspdi*z1r(i)+tspdr*z1i(i)
         x3r(i)=tssdr*z2r(i)-tssdi*z2i(i)
         x3i(i)=tssdi*z2r(i)+tssdr*z2i(i)
         x4r(i)=tpsdr*z2r(i)-tpsdi*z2i(i)
320      x4i(i)=tpsdi*z2r(i)+tpsdr*z2i(i)
      else
         do 330 i=1,nfreq
         x1r(i)=tppdr*z1r(i)-tppdi*z1i(i)
         x1i(i)=tppdi*z1r(i)+tppdr*z1i(i)
         x2r(i)=tspdr*z2r(i)-tspdi*z2i(i)
         x2i(i)=tspdi*z2r(i)+tspdr*z2i(i)
         x3r(i)=tssdr*z2r(i)-tssdi*z2i(i)
         x3i(i)=tssdi*z2r(i)+tssdr*z2i(i)
         x4r(i)=tpsdr*z1r(i)-tpsdi*z1i(i)
330      x4i(i)=tpsdi*z1r(i)+tpsdr*z1i(i)
      end if
      call mmtml(x1r,x1i,x2r,x2i,x4r,x4i,x3r,x3i,tppr,tppi,tspr,tspi,
     &tpsr,tpsi,tssr,tssi,x11r,x11i,x12r,x12i,x21r,x21i,x22r,x22i,nfreq)
      do 340 i=1,nfreq
      tppr(i)=x11r(i)
      tppi(i)=x11i(i)
      tspr(i)=x12r(i)
      tspi(i)=x12i(i)
      tpsr(i)=x21r(i)
      tpsi(i)=x21i(i)
      tssr(i)=x22r(i)
      tssi(i)=x22i(i)
340   continue
300   continue
      call mmtml(xppr,xppi,xspr,xspi,xpsr,xpsi,xssr,xssi,x11r,x11i,x12r
     &,x12i,x21r,x21i,x22r,x22i,tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,
     &nfreq)
      return
 2000 if(iss(1).eq.0) then
         return
      else if(iss(1).eq.1) then
         do 2100 i=1,nfreq
         tssr(i)=tppr(i)
 2100    tssi(i)=tppi(i)
      else
         do 2200 i=1,nfreq
         tppr(i)=tssr(i)
 2200    tppi(i)=tssi(i)
      end if
      return
      end
c==================================================================
      subroutine transn(n1,n2,i1,alpr,alpi,betr,beti,rho,d,p,tppr,tpp
     &i,tspr,tspi,tssr,tssi,tpsr,tpsi,rppr,rppi,rspr,rspi,rssr,rssi,rpsr
     &,rpsi,it)
c==================================================================
c
c     transn computes the plane wave transmissivity for one
c     slowness (p) and all required frequencies for n layers with a
c             recursive algorthm described by kennett:
c        theoretical reflection seismograms for elastic media,
c                    geoph. prosp.,v.27,p.301-321,1979
c
c     by complex velocities damping is considered
c
c     i1 =  1  -  transmissivity for an incident upgoing wavefield
c        = -1  -  transmissivity for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer, if n2=0 - free surface
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch oberste bzw. unterste schicht
c        = 1   -  transmission durch oberste bzw. unterste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex c1,al(nla),be(nla),aq1,bq1,cosi,cosj
      dimension alpr(nla),betr(nla),alpi(nla),beti(nla),rho(nla),
     &d(nla)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/refk/rppdr,rppdi,rspdr,rspdi,rpsdr,rpsdi,rssdr,rssdi,
     &rppur,rppui,rspur,rspui,rpsur,rpsui,rssur,rssui/transk/tppdr,
     &tppdi,tspdr,tspdi,tpsdr,tpsdi,tssdr,tssdi,tppur,tppui,tspur,
     &tspui,tpsur,tpsui,tssur,tssui
      common x1r(nfr),x1i(nfr),x2r(nfr),x2i(nfr),x3r(nfr),x3i(nfr),
     &x4r(nfr),x4i(nfr),w4r(nfr),w4i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),r1r(nfr),r1i(nfr),r2r(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr)
      dimension rppr(nfr),rppi(nfr),rspr(nfr),rspi(nfr),rssr(nfr),rssi(n
     &fr),rpsr(nfr),rpsi(nfr),
     &tppr(nfr),tppi(nfr),tspr(nfr),tspi(nfr),tssr(nfr),tssi(nfr),
     &tpsr(nfr),tpsi(nfr)
c
c     computation of the reflection and transmission coefficients at
c                             1 interface
c
      if(n1.eq.n2) then
         call set(tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,1.,nfreq)
         go to 11
      end if
      if(n2.ne.0) then
        n3=n2+i1
        al(n2)=cmplx(alpr(n2),alpi(n2))
        be(n2)=cmplx(betr(n2),beti(n2))
        al(n3)=cmplx(alpr(n3),alpi(n3))
        be(n3)=cmplx(betr(n3),beti(n3))
        call reflex(p,rho(n3),rho(n2),al(n3),al(n2),be(n3),be(n2))
      else
        n3=n2+1
        al(n3)=cmplx(alpr(n3),alpi(n3))
        be(n3)=cmplx(betr(n3),beti(n3))
        call refls(p,al(n3),be(n3))
      end if
      if(i1.eq.1.and.n2.ne.0) then
         rpsdr=-rpsdr
         rpsdi=-rpsdi
         rspdr=-rspdr
         rspdi=-rspdi
         rpsur=-rpsur
         rpsui=-rpsui
         rspur=-rspur
         rspui=-rspui
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
         tpsur=-tpsur
         tpsui=-tpsui
         tspur=-tspur
         tspui=-tspui
      end if
      if(n2.eq.0) then
         call set(tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,1.,nfreq)
      else
         do 301 i=1,nfreq
         tppr(i)=tppdr
         tppi(i)=tppdi
         tspr(i)=tspdr
         tspi(i)=tspdi
         tpsr(i)=tpsdr
         tpsi(i)=tpsdi
         tssr(i)=tssdr
 301     tssi(i)=tssdi
      end if
      do 217 i=1,nfreq
      r1r(i)=rppdr
      r1i(i)=rppdi
      r2r(i)=rspdr
      r2i(i)=rspdi
      r3r(i)=rssdr
      r3i(i)=rssdi
      rpsr(i)=rpsdr
 217  rpsi(i)=rpsdi
      nb=n2+i1
      ne=n1-i1
      nend=nb-ne
      if(i1.gt.0) then
         nend=-nend
      end if
      pi=3.14159265*2.
      if(nend.lt.0) go to 11
c
c     begin of the loop over the layers
c     the values which are independent of the frequency are calculated
c
      do 10 j=nb,ne,i1
      al(j)=cmplx(alpr(j),alpi(j))
      al(j+i1)=cmplx(alpr(j+i1),alpi(j+i1))
      be(j)=cmplx(betr(j),beti(j))
      be(j+i1)=cmplx(betr(j+i1),beti(j+i1))
      aq1=csqrt(1./(al(j)*al(j))-p*p)
      if(real(aq1).eq.0.) aq1=-aq1
      bq1=csqrt(1./(be(j)*be(j))-p*p)
      if(real(bq1).eq.0.) bq1=-bq1
      c1=aq1/bq1
      c1r=real(c1)
      c1i=aimag(c1)
      call reflex(p,rho(j+i1),rho(j),al(j+i1),al(j),be(j+i1),be(j))
      if(i1.eq.1) then
         rpsdr=-rpsdr
         rpsdi=-rpsdi
         rspdr=-rspdr
         rspdi=-rspdi
         rpsur=-rpsur
         rpsui=-rpsui
         rspur=-rspur
         rspui=-rspui
         tpsdr=-tpsdr
         tpsdi=-tpsdi
         tspdr=-tspdr
         tspdi=-tspdi
         tpsur=-tpsur
         tpsui=-tpsui
         tspur=-tspur
         tspui=-tspui
      end if
      pp1r=tppur*tppdr-tppui*tppdi
      pp1i=tppui*tppdr+tppur*tppdi
      pp2r=tspur*tpsdr-tspui*tpsdi
      pp2i=tspui*tpsdr+tspur*tpsdi
      pp3r=tppur*tpsdr-tppui*tpsdi
      pp3i=tppui*tpsdr+tppur*tpsdi
      pp4r=tspur*tppdr-tspui*tppdi
      pp4i=tspui*tppdr+tspur*tppdi
      pp5r=pp1r*rssur-pp1i*rssui
      pp5i=pp1i*rssur+pp1r*rssui
      pp6r=pp2r*rppur-pp2i*rppui
      pp6i=pp2i*rppur+pp2r*rppui
      pp7r=pp3r*rspur-pp3i*rspui
      pp7i=pp3i*rspur+pp3r*rspui
      pp8r=pp4r*rpsur-pp4i*rpsui
      pp8i=pp4i*rpsur+pp4r*rpsui
      pp3r=-pp4r*c1r+pp4i*c1i+pp3r
      pp3i=-pp4i*c1r-pp4r*c1i+pp3i
      pp4r=-pp5r-pp6r+pp7r+pp8r
      pp4i=-pp5i-pp6i+pp7i+pp8i
      pp5r=pp4r*c1r-pp4i*c1i
      pp5i=pp4i*c1r+pp4r*c1i
      sp1r=tppur*tspdr-tppui*tspdi
      sp1i=tppui*tspdr+tppur*tspdi
      sp2r=tspur*tssdr-tspui*tssdi
      sp2i=tspui*tssdr+tspur*tssdi
      sp3r=tppur*tssdr-tppui*tssdi
      sp3i=tppui*tssdr+tppur*tssdi
      sp4r=tspur*tspdr-tspui*tspdi
      sp4i=tspui*tspdr+tspur*tspdi
      sp5r=sp1r*rssur-sp1i*rssui
      sp5i=sp1i*rssur+sp1r*rssui
      sp6r=sp2r*rppur-sp2i*rppui
      sp6i=sp2i*rppur+sp2r*rppui
      sp7r=sp3r*rspur-sp3i*rspui
      sp7i=sp3i*rspur+sp3r*rspui
      sp8r=sp4r*rpsur-sp4i*rpsui
      sp8i=sp4i*rpsur+sp4r*rpsui
      sp3r=-sp4r*c1r+sp4i*c1i+sp3r
      sp3i=-sp4i*c1r-sp4r*c1i+sp3i
      sp4r=-sp5r-sp6r+sp7r+sp8r
      sp4i=-sp5i-sp6i+sp7i+sp8i
      sp5r=sp4r*c1r-sp4i*c1i
      sp5i=sp4i*c1r+sp4r*c1i
      ss1r=tpsur*tspdr-tpsui*tspdi
      ss1i=tpsui*tspdr+tpsur*tspdi
      ss2r=tssur*tssdr-tssui*tssdi
      ss2i=tssui*tssdr+tssur*tssdi
      ss3r=tpsur*tssdr-tpsui*tssdi
      ss3i=tpsui*tssdr+tpsur*tssdi
      ss4r=tssur*tspdr-tssui*tspdi
      ss4i=tssui*tspdr+tssur*tspdi
      ss5r=ss1r*rssur-ss1i*rssui
      ss5i=ss1i*rssur+ss1r*rssui
      ss6r=ss2r*rppur-ss2i*rppui
      ss6i=ss2i*rppur+ss2r*rppui
      ss7r=ss3r*rspur-ss3i*rspui
      ss7i=ss3i*rspur+ss3r*rspui
      ss8r=ss4r*rpsur-ss4i*rpsui
      ss8i=ss4i*rpsur+ss4r*rpsui
      ss3r=-ss4r*c1r+ss4i*c1i+ss3r
      ss3i=-ss4i*c1r-ss4r*c1i+ss3i
      ss4r=-ss5r-ss6r+ss7r+ss8r
      ss4i=-ss5i-ss6i+ss7i+ss8i
      ss5r=ss4r*c1r-ss4i*c1i
      ss5i=ss4i*c1r+ss4r*c1i
      det1r=-rppur
      det1i=-rppui
      det2r=-rssur
      det2i=-rssui
      det3r=rspur*c1r-rspui*c1i-rpsur
      det3i=rspui*c1r+rspur*c1i-rpsui
      det4r=rssur*rppur-rssui*rppui-rspur*rpsur+rspui*rpsui
      det4i=rssui*rppur+rssur*rppui-rspui*rpsur-rspur*rpsui
      det5r=det4r*c1r-det4i*c1i
      det5i=det4i*c1r+det4r*c1i
c
c     begin of the vectorizable loop over the frequencies
c
      al(j+i1)=cmplx(alpr(j+i1),alpi(j+i1))
      be(j+i1)=cmplx(betr(j+i1),beti(j+i1))
      aq1=csqrt(1./(al(j+i1)*al(j+i1))-p*p)
      if(real(aq1).eq.0.) aq1=-aq1
      bq1=csqrt(1./(be(j+i1)*be(j+i1))-p*p)
      if(real(bq1).eq.0.) bq1=-bq1
      c1=aq1/bq1
      c1r=real(c1)
      c1i=aimag(c1)
c     cosi=csqrt(1./(al(j)*al(j))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(j)*be(j))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     tpr=-d(j)*cosir
c     tpi=-d(j)*cosii
c     tsr=-d(j)*cosjr
c     tsi=-d(j)*cosji
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(j)
      if(omegar.eq.0.) then
         cosi=csqrt(1./(al(j)*al(j))-p*p)
         cosj=csqrt(1./(be(j)*be(j))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,j)
      end if
      do 216 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      z11r(i)=z1r(i)*z1r(i)-z1i(i)*z1i(i)
      z11i(i)=2.*z1r(i)*z1i(i)
      z22r(i)=z2r(i)*z2r(i)-z2i(i)*z2i(i)
      z22i(i)=2.*z2r(i)*z2i(i)
      z12r(i)=z1r(i)*z2r(i)-z1i(i)*z2i(i)
      z12i(i)=z1i(i)*z2r(i)+z1r(i)*z2i(i)
c     calculation of transmissivity
      w1r(i)=tppr(i)
      w1i(i)=tppi(i)
      w2r(i)=tspr(i)
      w2i(i)=tspi(i)
      w3r(i)=tpsr(i)
      w3i(i)=tpsi(i)
      w4r(i)=tssr(i)
      w4i(i)=tssi(i)
      tppr(i)=z1r(i)*w1r(i)-z1i(i)*w1i(i)
      tppi(i)=z1i(i)*w1r(i)+z1r(i)*w1i(i)
      tspr(i)=z2r(i)*w2r(i)-z2i(i)*w2i(i)
      tspi(i)=z2i(i)*w2r(i)+z2r(i)*w2i(i)
      tpsr(i)=z1r(i)*w3r(i)-z1i(i)*w3i(i)
      tpsi(i)=z1i(i)*w3r(i)+z1r(i)*w3i(i)
      tssr(i)=z2r(i)*w4r(i)-z2i(i)*w4i(i)
      tssi(i)=z2i(i)*w4r(i)+z2r(i)*w4i(i)
      x1r(i)=z11r(i)*r1r(i)-z11i(i)*r1i(i)
      x1i(i)=z11i(i)*r1r(i)+z11r(i)*r1i(i)
      x2r(i)=z22r(i)*r3r(i)-z22i(i)*r3i(i)
      x2i(i)=z22i(i)*r3r(i)+z22r(i)*r3i(i)
      x3r(i)=z12r(i)*r2r(i)-z12i(i)*r2i(i)
      x3i(i)=z12i(i)*r2r(i)+z12r(i)*r2i(i)
      rppr(i)=x1r(i)
      rppi(i)=x1i(i)
      rspr(i)=x3r(i)
      rspi(i)=x3i(i)
      rssr(i)=x2r(i)
      rssi(i)=x2i(i)
      x4r(i)=z12r(i)*rpsr(i)-z12i(i)*rpsi(i)
      z11i(i)=z12i(i)*rpsr(i)+z12r(i)*rpsi(i)
      z11r(i)=x4r(i)
      x4r(i)=x1r(i)*x2r(i)-x1i(i)*x2i(i)
      x4i(i)=x1i(i)*x2r(i)+x1r(i)*x2i(i)
      w4r(i)=x3r(i)*x3r(i)-x3i(i)*x3i(i)
      w4i(i)=2.*x3r(i)*x3i(i)
      detr(i)=1.+det1r*x1r(i)-det1i*x1i(i)+det2r*x2r(i)-det2i*x2i(i)+
     &det3r*x3r(i)-det3i*x3i(i)+det4r*x4r(i)-det4i*x4i(i)+det5r*w4r(i)-
     &det5i*w4i(i)
      deti(i)=det1i*x1r(i)+det1r*x1i(i)+det2i*x2r(i)+det2r*x2i(i)+
     &det3i*x3r(i)+det3r*x3i(i)+det4i*x4r(i)+det4r*x4i(i)+det5i*w4r(i)+
     &det5r*w4i(i)
      dett(i)=detr(i)*detr(i)+deti(i)*deti(i)
      dett(i)=1./dett(i)
      w1r(i)=pp1r*x1r(i)-pp1i*x1i(i)+pp2r*x2r(i)-pp2i*x2i(i)+pp3r*x3r(i)
     &-pp3i*x3i(i)+pp4r*x4r(i)-pp4i*x4i(i)+pp5r*w4r(i)-pp5i*w4i(i)
      w1i(i)=pp1i*x1r(i)+pp1r*x1i(i)+pp2i*x2r(i)+pp2r*x2i(i)+pp3i*x3r(i)
     &+pp3r*x3i(i)+pp4i*x4r(i)+pp4r*x4i(i)+pp5i*w4r(i)+pp5r*w4i(i)
      w2r(i)=sp1r*x1r(i)-sp1i*x1i(i)+sp2r*x2r(i)-sp2i*x2i(i)+sp3r*x3r(i)
     &-sp3i*x3i(i)+sp4r*x4r(i)-sp4i*x4i(i)+sp5r*w4r(i)-sp5i*w4i(i)
      w2i(i)=sp1i*x1r(i)+sp1r*x1i(i)+sp2i*x2r(i)+sp2r*x2i(i)+sp3i*x3r(i)
     &+sp3r*x3i(i)+sp4i*x4r(i)+sp4r*x4i(i)+sp5i*w4r(i)+sp5r*w4i(i)
      w3r(i)=ss1r*x1r(i)-ss1i*x1i(i)+ss2r*x2r(i)-ss2i*x2i(i)+ss3r*x3r(i)
     &-ss3i*x3i(i)+ss4r*x4r(i)-ss4i*x4i(i)+ss5r*w4r(i)-ss5i*w4i(i)
      w3i(i)=ss1i*x1r(i)+ss1r*x1i(i)+ss2i*x2r(i)+ss2r*x2i(i)+ss3i*x3r(i)
     &+ss3r*x3i(i)+ss4i*x4r(i)+ss4r*x4i(i)+ss5i*w4r(i)+ss5r*w4i(i)
      r1r(i)=rppdr+(w1r(i)*detr(i)+w1i(i)*deti(i))*dett(i)
      r1i(i)=rppdi+(w1i(i)*detr(i)-w1r(i)*deti(i))*dett(i)
      r2r(i)=rspdr+(w2r(i)*detr(i)+w2i(i)*deti(i))*dett(i)
      r2i(i)=rspdi+(w2i(i)*detr(i)-w2r(i)*deti(i))*dett(i)
      r3r(i)=rssdr+(w3r(i)*detr(i)+w3i(i)*deti(i))*dett(i)
      r3i(i)=rssdi+(w3i(i)*detr(i)-w3r(i)*deti(i))*dett(i)
      rpsr(i)=-r2r(i)*c1r+r2i(i)*c1i
      rpsi(i)=-r2i(i)*c1r-r2r(i)*c1i
216   continue
      do 211 i=1,nfreq
      x1r(i)=rppur
      x1i(i)=rppui
      x2r(i)=rspur
      x2i(i)=rspui
      x3r(i)=rpsur
      x3i(i)=rpsui
      x4r(i)=rssur
 211  x4i(i)=rssui
      call mmtml(x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,rppr,rppi,rspr,rspi,
     &z11r,z11i,rssr,rssi,w1r,w1i,w2r,w2i,w3r,w3i,w4r,w4i,nfreq)
      call matinv(w1r,w1i,w2r,w2i,w3r,w3i,w4r,w4i,x1r,x1i,x2r,x2i,x3r,
     &x3i,x4r,x4i,nfreq)
      do 212 i=1,nfreq
      w1r(i)=tppdr
      w1i(i)=tppdi
      w2r(i)=tspdr
      w2i(i)=tspdi
      w3r(i)=tpsdr
      w3i(i)=tpsdi
      w4r(i)=tssdr
 212  w4i(i)=tssdi
      call mmtml(x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,w1r,w1i,w2r,w2i,w3r,
     &w3i,w4r,w4i,z1r,z1i,z2r,z2i,z11r,z11i,z22r,z22i,nfreq)
      call mmtml(tppr,tppi,tspr,tspi,tpsr,tpsi,tssr,tssi,z1r,z1i,z2r,
     &z2i,z11r,z11i,z22r,z22i,w1r,w1i,w2r,w2i,w3r,w3i,w4r,w4i,nfreq)
      do 213 i=1,nfreq
      tppr(i)=w1r(i)
      tppi(i)=w1i(i)
      tspr(i)=w2r(i)
      tspi(i)=w2i(i)
      tpsr(i)=w3r(i)
      tpsi(i)=w3i(i)
      tssr(i)=w4r(i)
 213  tssi(i)=w4i(i)
10    continue
11    al(n1)=cmplx(alpr(n1),alpi(n1))
      be(n1)=cmplx(betr(n1),beti(n1))
      aq1=csqrt(1./(al(n1)*al(n1))-p*p)
      if(real(aq1).eq.0.) aq1=-aq1
      bq1=csqrt(1./(be(n1)*be(n1))-p*p)
      if(real(bq1).eq.0.) bq1=-bq1
      c1=aq1/bq1
      c1r=real(c1)
      c1i=aimag(c1)
      if(it.eq.2) then
c
c     transmission durch oberste und unterste schicht
c
c     cosi=csqrt(1./(al(n1)*al(n1))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(n1)*be(n1))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     tpr=-d(n1)*cosir
c     tpi=-d(n1)*cosii
c     tsr=-d(n1)*cosjr
c     tsi=-d(n1)*cosji
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(n1)
      if(omegar.eq.0.) then
         cosi=csqrt(1./(al(n1)*al(n1))-p*p)
         cosj=csqrt(1./(be(n1)*be(n1))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,n1)
      end if
      do 219 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      w1r(i)=z1r(i)*tppr(i)-z1i(i)*tppi(i)
      w1i(i)=z1i(i)*tppr(i)+z1r(i)*tppi(i)
      w2r(i)=z2r(i)*tspr(i)-z2i(i)*tspi(i)
      w2i(i)=z2i(i)*tspr(i)+z2r(i)*tspi(i)
      w3r(i)=z1r(i)*tpsr(i)-z1i(i)*tpsi(i)
      w3i(i)=z1i(i)*tpsr(i)+z1r(i)*tpsi(i)
      w4r(i)=z2r(i)*tssr(i)-z2i(i)*tssi(i)
      w4i(i)=z2i(i)*tssr(i)+z2r(i)*tssi(i)
219   continue
      if(n1.eq.n2) go to 222
c     cosi=csqrt(1./(al(n2)*al(n2))-p*p)
c     cosir=real(cosi)
c     cosii=aimag(cosi)
c     cosii=abs(cosii)
c     cosi=csqrt(1./(be(n2)*be(n2))-p*p)
c     cosjr=real(cosi)
c     cosji=aimag(cosi)
c     cosji=abs(cosji)
c     if(n1.ne.n2) then
c        tpr=-d(n2)*cosir
c        tpi=-d(n2)*cosii
c        tsr=-d(n2)*cosjr
c        tsi=-d(n2)*cosji
c     else
c        tpr=0.
c        tpi=0.
c        tsr=0.
c        tsi=0.
c     end if
c     tpri=omegai*tpr
c     tpii=omegai*tpi
c     tsri=omegai*tsr
c     tsii=omegai*tsi
      dj=-d(n2)
      if(omegar.eq.0.) then
         cosi=csqrt(1./(al(n2)*al(n2))-p*p)
         cosj=csqrt(1./(be(n2)*be(n2))-p*p)
         call phas1(dj,omegai,cosi,cosj,z1r,z1i,z2r,z2i,x3r,x3i,x4r,
     &              x4i,omega,nfreq)
      else
         p2i=p*p
         call phas2(dj,omegai,z1r,z1i,z2r,z2i,x3r,x3i,x4r,x4i,
     &              p2i,omega,nfreq,n2)
      end if
      do 220 i=1,nfreq
      x1r(i)=z1i(i)-x3r(i)
      x1i(i)=z1r(i)+x3i(i)
      x2r(i)=z2i(i)-x4r(i)
      x2i(i)=z2r(i)+x4i(i)
c     x1r(i)=omega(i)*tpi-tpri
c     x1i(i)=omega(i)*tpr+tpii
c     x2r(i)=omega(i)*tsi-tsri
c     x2i(i)=omega(i)*tsr+tsii
      x1r(i)=exp(x1r(i))
      x2r(i)=exp(x2r(i))
      z1r(i)=x1r(i)*cos(x1i(i))
      z1i(i)=x1r(i)*sin(x1i(i))
      z2r(i)=x2r(i)*cos(x2i(i))
      z2i(i)=x2r(i)*sin(x2i(i))
      tppr(i)=z1r(i)*w1r(i)-z1i(i)*w1i(i)
      tppi(i)=z1i(i)*w1r(i)+z1r(i)*w1i(i)
      tspr(i)=z1r(i)*w2r(i)-z1i(i)*w2i(i)
      tspi(i)=z1i(i)*w2r(i)+z1r(i)*w2i(i)
      tpsr(i)=z2r(i)*w3r(i)-z2i(i)*w3i(i)
      tpsi(i)=z2i(i)*w3r(i)+z2r(i)*w3i(i)
      tssr(i)=z2r(i)*w4r(i)-z2i(i)*w4i(i)
      tssi(i)=z2i(i)*w4r(i)+z2r(i)*w4i(i)
220   continue
      go to 224
222   do 223 i=1,nfreq
      tppr(i)=w1r(i)
      tppi(i)=w1i(i)
      tspr(i)=w2r(i)
      tspi(i)=w2i(i)
      tpsr(i)=w3r(i)
      tpsi(i)=w3i(i)
      tssr(i)=w4r(i)
      tssi(i)=w4i(i)
223   continue
224   end if
      do 225 i=1,nfreq
      rppr(i)=r1r(i)
      rppi(i)=r1i(i)
      rspr(i)=r2r(i)
      rspi(i)=r2i(i)
      rssr(i)=r3r(i)
      rssi(i)=r3i(i)
      rpsr(i)=-rspr(i)*c1r+rspi(i)*c1i
 225  rpsi(i)=-rspi(i)*c1r-rspr(i)*c1i
      if(iss(1).eq.1) then
         do 230 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
         rspi(i)=0.
         rssr(i)=0.
         rssi(i)=0.
         tpsr(i)=0.
         tpsi(i)=0.
         tspr(i)=0.
         tspi(i)=0.
         tssr(i)=0.
 230     tssi(i)=0.
      else if(iss(1).eq.2) then
         do 240 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
         rspi(i)=0.
         rppr(i)=0.
         rppi(i)=0.
         tpsr(i)=0.
         tpsi(i)=0.
         tspr(i)=0.
         tspi(i)=0.
         tppr(i)=0.
 240     tppi(i)=0.
      else if(iss(1).eq.3) then
         do 250 i=1,nfreq
         rpsr(i)=0.
         rpsi(i)=0.
         rspr(i)=0.
         rspi(i)=0.
         tpsr(i)=0.
         tpsi(i)=0.
         tspr(i)=0.
 250     tspi(i)=0.
      end if
      return
      end
c==================================================================
      subroutine tkoeff(u,rho1,rho2,al1,al2,be1,be2)
c==================================================================
      common/transc/tppdr,tppdi,tspdr,tspdi,tpsdr,tpsdi,tssdr,tssdi
c
c       matrizen der transmissionskoeffizienten einer trennflaeche
c       zwischen zwei halbraeumen fuer welleneinfall von oben
c
        complex tppd,tspd,tpsd,tssd,
     *  a1,b1,a2,b2,d1d,d2d,
     *  d,a1b1,a2b2,a1b2,a2b1,t1,t2,t3,t4,t5,t6,t7,t10,
     *  c,c0,c1,c2,c3,aq1,aq2,bq1,bq2,al1,al2,be1,be2
        aq1=1./(al1*al1)
        aq2=1./(al2*al2)
        uq=u*u
        a1=csqrt(aq1-uq)
        if(real(a1).eq.0.) a1=-a1
        a2=csqrt(aq2-uq)
        if(real(a2).eq.0.) a2=-a2
        if(be1.eq.0.) then
           bq2=1./(be2*be2)
           b2=csqrt(bq2-uq)
           if(real(b2).eq.0.) b2=-b2
           t3=2.*be2**2*rho2*rho1*a1
           c2=1./be2**2-2.*uq
           t1=(be2**2*rho2)**2*a1*(c2**2+4.*uq*a2*b2)
           c0=rho1*rho2*a2
           d=1./(c0+t1)
           tppd=t3*c2*d
           tpsd=cmplx(0.,1.)*2.*t3*u*a2*d
           tssd=0.
           tspd=0.
           go to 100
        else if(be2.eq.0.) then
           bq1=1./(be1*be1)
           b1=csqrt(bq1-uq)
           if(real(b1).eq.0.) b1=-b1
           t3=2.*be1**2*rho1*rho2*a2
           c2=1./be1**2-2.*uq
           t1=(be1**2*rho1)**2*a2*(c2**2+4.*uq*a1*b1)
           c0=rho2*rho1*a1
           d=1./(c0+t1)
           tppd=t3*c2*d*(rho1*a1)/(rho2*a2)
           tspd=cmplx(0.,1.)*2.*t3*u*a1*d*(rho1*b1)/(rho2*a2)
           tssd=0.
           tpsd=0.
           go to 100
        end if
        bq1=1./(be1*be1)
        bq2=1./(be2*be2)
        c=2.*(rho1/bq1-rho2/bq2)
        b1=csqrt(bq1-uq)
        if(real(b1).eq.0.) b1=-b1
        b2=csqrt(bq2-uq)
        if(real(b2).eq.0.) b2=-b2
        c0=c*uq
        c1=c0-rho1
        c2=c0+rho2
        c3=c1+rho2
        a1b1=a1*b1
        a2b2=a2*b2
        a1b2=a1*b2
        a2b1=a2*b1
        rho12=rho1*rho2
        t1=c1*c1*a2b2+rho12*a2b1
        t2=c2*c2*a1b1+rho12*a1b2
        t3=c3*c3*uq
        t4=c*c0*a1b1*a2b2
        d1d=t3+t1
        d2d=t4+t2
        d=d1d+d2d
        t5=2./d
        t6=a1*t5
        t7=b1*t5
        t10=c2*b1-c1*b2
        tppd=rho1*t6*t10
        t10=u*(c3+c*a2b1)
        tpsd=-rho1*t6*t10
        t10=c2*a1-c1*a2
        tssd=rho1*t7*t10
        t10=u*(c3+c*a1b2)
        tspd=rho1*t7*t10
100   tppdr=real(tppd)
      tppdi=aimag(tppd)
      tspdr=real(tspd)
      tspdi=aimag(tspd)
      tpsdr=real(tpsd)
      tpsdi=aimag(tpsd)
      tssdr=real(tssd)
      tssdi=aimag(tssd)
        return
        end
c==================================================================
        subroutine reflex(u,rho1,rho2,al1,al2,be1,be2)
c==================================================================
      common/refk/rppdr,rppdi,rspdr,rspdi,rpsdr,rpsdi,rssdr,rssdi,
     &rppur,rppui,rspur,rspui,rpsur,rpsui,rssur,rssui/transk/tppdr,
     &tppdi,tspdr,tspdi,tpsdr,tpsdi,tssdr,tssdi,tppur,tppui,tspur,
     &tspui,tpsur,tpsui,tssur,tssui
c
c       matrizen der reflexions- und transmissionskoeffizienten
c       einer trennflaeche zwischen zwei halbraeumen
c
        complex rppd,rspd,rpsd,rssd,tppd,tspd,tpsd,tssd,rppu,rspu,
     *  rpsu,rssu,tppu,tspu,tpsu,tssu,a1,b1,a2,b2,d1d,d2d,d1u,d2u,
     *  d,a1b1,a2b2,a1b2,a2b1,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,
     *  c,c0,c1,c2,c3,aq1,aq2,bq1,bq2,al1,al2,be1,be2
        aq1=1./(al1*al1)
        aq2=1./(al2*al2)
        uq=u*u
        a1=csqrt(aq1-uq)
        if(real(a1).eq.0.) a1=-a1
        a2=csqrt(aq2-uq)
        if(real(a2).eq.0.) a2=-a2
        if(be1.eq.0.) then
           bq2=1./(be2*be2)
           b2=csqrt(bq2-uq)
           if(real(b2).eq.0.) b2=-b2
           c2=1./be2**2-2.*uq
           t1=(be2**2*rho2)**2*a1*(c2**2+4.*uq*a2*b2)
           t2=(be2**2*rho2)**2*a1*(c2**2-4.*uq*a2*b2)
           t3=2.*be2**2*rho2*rho1*a1
           t4=4.*(be2**2*rho2)**2*u
           c0=rho1*rho2*a2
           d=1./(c0+t1)
           rppd=(-c0+t1)*d
           rppu=(c0-t2)*d
           rssu=(-c0-t2)*d
           rssd=0.
           rspu=(-t4*cmplx(0.,1.)*a1*b2*c2)*d
           rspd=0.
           rpsu=rspu*a2/b2
           rpsd=0.
           tppd=t3*c2*d
           tppu=tppd*(rho2*a2)/(rho1*a1)
           tpsd=cmplx(0.,1.)*2.*t3*u*a2*d
           tpsu=0.
           tssu=0.
           tssd=0.
           tspd=0.
           tspu=tpsd*(rho2*b2)/(rho1*a1)
           go to 100
        else if(be2.eq.0.) then
           bq1=1./(be1*be1)
           b1=csqrt(bq1-uq)
           if(real(b1).eq.0.) b1=-b1
           c2=1./be1**2-2.*uq
           t1=(be1**2*rho1)**2*a2*(c2**2+4.*uq*a1*b1)
           t2=(be1**2*rho1)**2*a2*(c2**2-4.*uq*a1*b1)
           t3=2.*be1**2*rho1*rho2*a2
           t4=4.*(be1**2*rho1)**2*u
           c0=rho2*rho1*a1
           d=1./(c0+t1)
           rppu=(-c0+t1)*d
           rppd=(c0-t2)*d
           rssd=(-c0-t2)*d
           rssu=0.
           rspd=(-t4*cmplx(0.,1.)*a2*b1*c2)*d
           rspu=0.
           rpsd=rspd*a1/b1
           rspu=0.
           tppu=t3*c2*d
           tppd=tppu*(rho1*a1)/(rho2*a2)
           tpsu=cmplx(0.,1.)*2.*t3*u*a1*d
           tpsd=0.
           tssu=0.
           tssd=0.
           tspu=0.
           tspd=tpsu*(rho1*b1)/(rho2*a2)
           go to 100
        end if
        bq1=1./(be1*be1)
        bq2=1./(be2*be2)
        c=2.*(rho1/bq1-rho2/bq2)
        b1=csqrt(bq1-uq)
        if(real(b1).eq.0.) b1=-b1
        b2=csqrt(bq2-uq)
        if(real(b2).eq.0.) b2=-b2
        c0=c*uq
        c1=c0-rho1
        c2=c0+rho2
        c3=c1+rho2
        a1b1=a1*b1
        a2b2=a2*b2
        a1b2=a1*b2
        a2b1=a2*b1
        rho12=rho1*rho2
        t1=c1*c1*a2b2+rho12*a2b1
        t2=c2*c2*a1b1+rho12*a1b2
        t3=c3*c3*uq
        t4=c*c0*a1b1*a2b2
        d1d=t3+t1
        d2d=t4+t2
        d=d1d+d2d
        d1u=t3+t2
        d2u=t4+t1
        t5=2./d
        t6=a1*t5
        t7=b1*t5
        t8=a2*t5
        t9=b2*t5
        rppd=(d2d-d1d)/d
        rppu=(d2u-d1u)/d
        t10=u*(c3*c2+c*c1*a2b2)
        rpsd=-t6*t10
        rspd=t7*t10
        t10=rho12*(a1b2-a2b1)*t5
        rssd=rppd-t10
        rssu=rppu+t10
        t10=u*(c3*c1+c*c2*a1b1)
        rpsu=t8*t10
        rspu=-t9*t10
        t10=c2*b1-c1*b2
        tppd=rho1*t6*t10
        tppu=rho2*t8*t10
        t10=u*(c3+c*a2b1)
        tpsd=-rho1*t6*t10
        tspu=rho2*t9*t10
        t10=c2*a1-c1*a2
        tssd=rho1*t7*t10
        tssu=rho2*t9*t10
        t10=u*(c3+c*a1b2)
        tspd=rho1*t7*t10
        tpsu=-rho2*t8*t10
 100  tppdr=real(tppd)
      tppdi=aimag(tppd)
      tspdr=real(tspd)
      tspdi=aimag(tspd)
      tpsdr=real(tpsd)
      tpsdi=aimag(tpsd)
      tssdr=real(tssd)
      tssdi=aimag(tssd)
      tppur=real(tppu)
      tppui=aimag(tppu)
      tspur=real(tspu)
      tspui=aimag(tspu)
      tpsur=real(tpsu)
      tpsui=aimag(tpsu)
      tssur=real(tssu)
      tssui=aimag(tssu)
      rppdr=real(rppd)
      rppdi=aimag(rppd)
      rspdr=real(rspd)
      rspdi=aimag(rspd)
      rpsdr=real(rpsd)
      rpsdi=aimag(rpsd)
      rssdr=real(rssd)
      rssdi=aimag(rssd)
      rppur=real(rppu)
      rppui=aimag(rppu)
      rspur=real(rspu)
      rspui=aimag(rspu)
      rpsur=real(rpsu)
      rpsui=aimag(rpsu)
      rssur=real(rssu)
      rssui=aimag(rssu)
        return
        end
c==================================================================
        subroutine refls(u,al1,be1)
c==================================================================
c
c       reflexionskoeffizienten einer freien oberflaeche fuer wellen-
c       einfall von unten
c
      common/refk/rppdr,rppdi,rspdr,rspdi,rpsdr,rpsdi,rssdr,rssdi,
     &rppur,rppui,rspur,rspui,rpsur,rpsui,rssur,rssui
        complex rppu,rspu,rpsu,rssu,a1,b1,d1,d2,d,t1,t2,t3,aq1,bq1,al1,
     &  be1
        aq1=1./(al1*al1)
        uq=u*u
        a1=csqrt(aq1-uq)
        if(real(a1).eq.0.) a1=-a1
        if(be1.eq.0.) then
            rppu=-1.
            rpsu=0.
            rspu=0.
            rssu=0.
            go to 100
        end if
        bq1=1./(be1*be1)
        b1=csqrt(bq1-uq)
        if(real(b1).eq.0.) b1=-b1
        t1=2./bq1
        t2=t1*uq-1.
        d1=t2*t2
        d2=t1*t1*uq*a1*b1
        d=d1+d2
        t3=2.*t1*u*t2/d
        rppu=(d2-d1)/d
        rspu=-b1*t3
        rpsu=a1*t3
        rssu=rppu
100   rppdr=real(rppu)
      rppdi=aimag(rppu)
      rspdr=real(rspu)
      rspdi=aimag(rspu)
      rpsdr=real(rpsu)
      rpsdi=aimag(rpsu)
      rssdr=real(rssu)
      rssdi=aimag(rssu)
        return
        end
c==================================================================
        subroutine freesu(alpr,alpi,betr,beti,u,h11r,h11i,h12r,
     *  h12i,h22r,h22i,h21r,h21i,nfreq)
c==================================================================
c
c       matrix 2.*h
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
        dimension alpr(nla),alpi(nla),betr(nla),beti(nla),
     *  h11r(nfr),h11i(nfr),h12r(nfr),h12i(nfr),h21r(nfr),
     *  h21i(nfr),h22r(nfr),h22i(nfr)
        real h11r,h11i,h12r,h12i,h21r,h21i,h22r,h22i
        complex l11,l12,l21,l22,a1,b1,d1,d2,d,t1,t2,t3,aq1,bq1,al1,be1
        al1=cmplx(alpr(1),alpi(1))
        be1=cmplx(betr(1),beti(1))
        aq1=1./(al1*al1)
        bq1=1./(be1*be1)
        uq=u*u
        a1=csqrt(aq1-uq)
        if(real(a1).eq.0.) a1=-a1
        b1=csqrt(bq1-uq)
        if(real(b1).eq.0.) b1=-b1
        t1=2./bq1
        t2=t1*uq-1.
        d1=t2*t2
        d2=t1*t1*uq*a1*b1
        d=d1+d2
        t3=2./d
        l11=t1*u*a1*b1*t3
        l12=-b1*t2*t3
        l21=-a1*t2*t3
        l22=-l11
        do 100 i=1,nfreq
        h11r(i)=real(l11)
        h11i(i)=aimag(l11)
        h12r(i)=real(l12)
        h12i(i)=aimag(l12)
        h21r(i)=real(l21)
        h21i(i)=aimag(l21)
        h22r(i)=real(l22)
 100    h22i(i)=aimag(l22)
        return
        end
c==================================================================
        subroutine mmtml(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
     *  b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,
     *  c21r,c21i,c22r,c22i,n)
c==================================================================
        dimension a11r(n),a11i(n),a12r(n),a12i(n),a21r(n),a21i(n),
     *  a22r(n),a22i(n),
     *  b11r(n),b11i(n),b12r(n),b12i(n),b21r(n),b21i(n),b22r(n),b22i(n),
     *  c11r(n),c11i(n),c12r(n),c12i(n),c21r(n),c21i(n),c22r(n),c22i(n)
        do 100 i=1,n
        c11r(i)=a11r(i)*b11r(i)-a11i(i)*b11i(i)+a12r(i)*b21r(i)-
     *  a12i(i)*b21i(i)
        c11i(i)=a11i(i)*b11r(i)+a11r(i)*b11i(i)+a12i(i)*b21r(i)+
     *  a12r(i)*b21i(i)
        c12r(i)=a11r(i)*b12r(i)-a11i(i)*b12i(i)+a12r(i)*b22r(i)-
     *  a12i(i)*b22i(i)
        c12i(i)=a11i(i)*b12r(i)+a11r(i)*b12i(i)+a12i(i)*b22r(i)+
     *  a12r(i)*b22i(i)
        c21r(i)=a21r(i)*b11r(i)-a21i(i)*b11i(i)+a22r(i)*b21r(i)-
     *  a22i(i)*b21i(i)
        c21i(i)=a21i(i)*b11r(i)+a21r(i)*b11i(i)+a22i(i)*b21r(i)+
     *  a22r(i)*b21i(i)
        c22r(i)=a21r(i)*b12r(i)-a21i(i)*b12i(i)+a22r(i)*b22r(i)-
     *  a22i(i)*b22i(i)
        c22i(i)=a21i(i)*b12r(i)+a21r(i)*b12i(i)+a22i(i)*b22r(i)+
     *  a22r(i)*b22i(i)
 100    continue
        return
        end
c==================================================================
        subroutine matmus(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
c==================================================================
     *  b11r,b11i,b21r,b21i,c11r,c11i,c21r,c21i,n)
        dimension a11r(n),a11i(n),a12r(n),a12i(n),a21r(n),a21i(n),
     *  a22r(n),a22i(n),
     *  b11r(n),b11i(n),b21r(n),b21i(n),
     *  c11r(n),c11i(n),c21r(n),c21i(n)
        do 100 i=1,n
        c11r(i)=a11r(i)*b11r(i)-a11i(i)*b11i(i)+a12r(i)*b21r(i)-
     *  a12i(i)*b21i(i)
        c11i(i)=a11i(i)*b11r(i)+a11r(i)*b11i(i)+a12i(i)*b21r(i)+
     *  a12r(i)*b21i(i)
        c21r(i)=a21r(i)*b11r(i)-a21i(i)*b11i(i)+a22r(i)*b21r(i)-
     *  a22i(i)*b21i(i)
        c21i(i)=a21i(i)*b11r(i)+a21r(i)*b11i(i)+a22i(i)*b21r(i)+
     *  a22r(i)*b21i(i)
 100    continue
        return
        end
c==================================================================
        subroutine matinv(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
c==================================================================
     *  b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,n)
c        parameter(npa=2500,nfr=1024,nla=1000)
c NS_CHANGE
        parameter(npa=5000,nfr=131072,nla=1000)
        dimension a11r(n),a11i(n),a12r(n),a12i(n),a21r(n),a21i(n),
     *  a22r(n),a22i(n),
     *  b11r(n),b11i(n),b12r(n),b12i(n),b21r(n),b21i(n),b22r(n),b22i(n),
     *  d1r(nfr),d1i(nfr),d2r(nfr),d2i(nfr)
        do 50 i=1,n
        a11r(i)=1.-a11r(i)
        a11i(i)=-a11i(i)
        a22r(i)=1.-a22r(i)
        a22i(i)=-a22i(i)
        a12r(i)=-a12r(i)
        a12i(i)=-a12i(i)
        a21r(i)=-a21r(i)
 50     a21i(i)=-a21i(i)
        do 100 i=1,n
        d1r(i)=a11r(i)*a22r(i)-a11i(i)*a22i(i)
        d2r(i)=a12r(i)*a21r(i)-a12i(i)*a21i(i)
        d1i(i)=a11i(i)*a22r(i)+a11r(i)*a22i(i)
        d2i(i)=a12i(i)*a21r(i)+a12r(i)*a21i(i)
        d1r(i)=d1r(i)-d2r(i)
        d1i(i)=d1i(i)-d2i(i)
        d2r(i)=d1r(i)*d1r(i)+d1i(i)*d1i(i)
        d1r(i)=d1r(i)/d2r(i)
        d1i(i)=-d1i(i)/d2r(i)
        b11r(i)=a22r(i)*d1r(i)-a22i(i)*d1i(i)
        b11i(i)=a22i(i)*d1r(i)+a22r(i)*d1i(i)
        b12r(i)=-a12r(i)*d1r(i)+a12i(i)*d1i(i)
        b12i(i)=-a12i(i)*d1r(i)-a12r(i)*d1i(i)
        b21r(i)=-a21r(i)*d1r(i)+a21i(i)*d1i(i)
        b21i(i)=-a21i(i)*d1r(i)-a21r(i)*d1i(i)
        b22r(i)=a11r(i)*d1r(i)-a11i(i)*d1i(i)
        b22i(i)=a11i(i)*d1r(i)+a11r(i)*d1i(i)
 100    continue
        return
        end
c==================================================================
        subroutine matadd(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
c==================================================================
     *  b11r,b11i,b12r,b12i,b21r,b21i,b22r,b22i,c11r,c11i,c12r,c12i,
     *  c21r,c21i,c22r,c22i,n)
        dimension a11r(n),a11i(n),a12r(n),a12i(n),a21r(n),a21i(n),
     *  a22r(n),a22i(n),
     *  b11r(n),b11i(n),b12r(n),b12i(n),b21r(n),b21i(n),b22r(n),b22i(n),
     *  c11r(n),c11i(n),c12r(n),c12i(n),c21r(n),c21i(n),c22r(n),c22i(n)
        do 100 i=1,n
        c11r(i)=a11r(i)+b11r(i)
        c11i(i)=a11i(i)+b11i(i)
        c12r(i)=a12r(i)+b12r(i)
        c12i(i)=a12i(i)+b12i(i)
        c21r(i)=a21r(i)+b21r(i)
        c21i(i)=a21i(i)+b21i(i)
        c22r(i)=a22r(i)+b22r(i)
        c22i(i)=a22i(i)+b22i(i)
 100    continue
        return
        end
c==================================================================
        subroutine matadt(a11r,a11i,a12r,a12i,b11r,b11i,b12r,b12i,
c==================================================================
     *  c11r,c11i,c12r,c12i,n)
        dimension a11r(n),a11i(n),a12r(n),a12i(n),
     *  b11r(n),b11i(n),b12r(n),b12i(n),c11r(n),c11i(n),c12r(n),c12i(n)
        do 100 i=1,n
        c11r(i)=a11r(i)+b11r(i)
        c11i(i)=a11i(i)+b11i(i)
        c12r(i)=a12r(i)+b12r(i)
        c12i(i)=a12i(i)+b12i(i)
 100    continue
        return
        end
c==================================================================
        subroutine set(a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,xr,n)
c==================================================================
        dimension a11r(n),a11i(n),a12r(n),a12i(n),a21r(n),a21i(n),
     *  a22r(n),a22i(n)
        do 100 i=1,n
        a11r(i)=xr
        a11i(i)=0.
        a12r(i)=0.
        a12i(i)=0.
        a21r(i)=0.
        a21i(i)=0.
        a22r(i)=xr
        a22i(i)=0.
 100    continue
        return
        end
c==================================================================
      subroutine erdmod(iins,nh)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      common/blo1/z(0:nla),d(nla),ar(nla),ai(nla),br(nla),bi(nla),
     &rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      complex cv
      read(iins,51) iso,ire,mdeck,nrho,nh
 51   format(5i5)
      if(iss(20).eq.1.and.iss(17).ne.1) then
         read(iins,52) dso,dre
      end if
      if(iss(8).eq.1) go to 1000
      qsum=0.
      i=0
 100  i=i+1
      read(iins,52) d(i),ar(i),qp(i),br(i),qs(i),rho(i)
 52   format(6f10.4)
      if(qs(i).eq.0.) qs(i)=qp(i)
      if(qp(i).eq.0.) ai(i)=0.
      qsum=qsum+qp(i)+qs(i)
      if(qs(i).eq.0.) bi(i)=0.
      if(qp(i).ne.0.) ai(i)=ar(i)/(2.*qp(i))
      if(br(i).eq.0.)br(i)=ar(i)/1.732
      if(br(i).eq.-1.)br(i)=0.
      if(qs(i).ne.0.) bi(i)=br(i)/(2.*qs(i))
      if(rho(i).eq.0.)  rho(i)=0.252+0.3788*ar(i)
      if(d(i)+ar(i)+br(i))100,200,100
 200  nlay=i-1
      n1=nlay-1
      if(iss(9).eq.1) then
         do 300 i=1,nlay
 300     z(i)=d(i)
         d(1)=z(1)
         if(n1.lt.2) go to 210
         do 400 i=2,n1
 400     d(i)=z(i)-z(i-1)
      else
         z(1)=d(1)
         if(n1.lt.2) go to 210
         do 500 i=2,n1
 500     z(i)=z(i-1)+d(i)
      end if
c
c     write out layer parameter
c
 210  write(6,61)
 61   format(/1x,'layer parameter :'/)
 1200 write(6,62)
 62   format(7x, 'depth   thickness  p-velocity   qp     s-velocity   qs
     1      density      '/8x,2hkm,8x,2hkm,9x,4hkm/s,16x,4hkm/s,14x,  7h
     2g/cm**3/)
      do 600 i=1,n1
      if(i.gt.mdeck) go to 700
      write(6,63) i,z(i),d(i),ar(i),qp(i),br(i),qs(i),rho(i)
 63   format(1x,i4,7f12.4)
      go to 600
 700  write(6,64) i,z(i),d(i),ar(i),qp(i),br(i),qs(i),rho(i)
 64   format(1x,i4,7f12.4,5x,'refl. zone')
 600  continue
      z(nlay)=9999.999
      d(nlay)=9999.999
      write(6,65) ar(nlay),qp(nlay),br(nlay),qs(nlay),rho(nlay)
 65   format(24x,5f12.4,3x,'refl. zone')
      if(iss(14).eq.1) then
c **********************************************************************
c     calculation of the critical distancies for all layers
c
      write(6,8510)
 8510 format('      depth         p-veloci. crit.dist.   time    s-veloc
     &i. crit.dist.   time')
      do 8512 l1=2,nlay
      tia=0.
      tib=0.
      xcra=0.
      xcrb=0.
      do 8514 l2=1,l1-1
      if(ar(l2).ge.ar(l1).or.br(l2).ge.br(l1)) then
         write(6,8520) l1
 8520    format(' layer with index',i5,' not critical')
         go to 8512
      end if
      gama=asin(ar(l2)/ar(l1))
      gamb=asin(br(l2)/br(l1))
      ff=2.
      if(l2.lt.iso) ff=1.
      xcra=xcra+ff*d(l2)*tan(gama)
      xcrb=xcrb+ff*d(l2)*tan(gamb)
      tia=tia+ff*cos(gama)*d(l2)/ar(l2)
      tib=tib+ff*cos(gamb)*d(l2)/br(l2)
 8514 continue
      gama=asin(ar(1)/ar(l1))
      gamb=asin(br(1)/br(l1))
      tia=tia+xcra*sin(gama)/ar(1)
      tib=tib+xcrb*sin(gamb)/br(1)
      write(6,8516) l1,z(l1-1),ar(l1),xcra,tia,br(l1),xcrb,tib
 8516 format(i4,f10.4,' / ',3f10.4,' / ',3f10.4)
 8512 continue
c
c     calculation is finished
c **********************************************************************
      end if
      if(iss(17).eq.1) go to 1100
      if(iss(20).ne.1) then
         dso=0.
         ii=iso-1
         if(ii.eq.0)  go to 900
         do  910  i=1,ii
 910     dso=dso+d(i)
c
c
c
         if(iss(12).eq.1) dso=3390.*(1.-exp(-dso/3390.))
 900     dre=0.
         ii=ire
         if(ii.eq.0)  go to 920
         do  930  i=1,ii
 930     dre=dre+d(i)
      else
         nso=1
         nre=1
         do 980 i=1,nlay
         diffs=dso-z(i)
         if(diffs.lt.0) iso=nso
         if(diffs.ge.0) nso=nso+1
         diffr=dre-z(i)
         if(diffr.lt.0) ire=nre
         if(diffr.ge.0) nre=nre+1
 980     continue
         if(dre.eq.0) ire=0
      end if
 920  write(6,66) iso,dso
 66   format(///,1x,'index of the layer with the source =',i5,'    depth
     1 of the source = ',f8.2,'  km',/)
      if(ire.eq.0) then
         write(6,67)
 67      format(1x,'the receivers are situated at the free surface')
         go to 1100
      else
         write(6,68) ire,dre
 68   format(///,1x,'index of the layer with the receiver = ',i5/1x,'
     1         depth of the receiver = ',f8.2,'  km',//)
         go to 1100
      end if
c
c     automatic generation of the depth-velocity distribution (iss(8)=1)
c
 1000 call inhom(iins)
      n1=nlay-1
      write(6,69)
 69   format(//1x,44hlayer parameter after automatic generation :/)
      go to 1200
 1100 if(nlay.lt.2.or.nlay.gt.1000) then
         write(6,70)
 70      format(/1x,'check number of layers')
      end if
      return
      end
c==================================================================
      subroutine output(dentif,nent,iouts)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      common/blo1/z(0:nla),d(nla),ar(nla),ai(nla),br(nla),bi(nla),
     &rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      complex cv
      character*80 dentif
cc      write(iouts,100) dentif
  100 format(10x,a42)
      iss4=iss(4)
      if(iouts.eq.3) iss4=3
      if(iouts.eq.21) iss4=1
      if(iouts.eq.22) iss4=2
      if(iouts.eq.23) iss4=3
cc      write(iouts,110) nlay,mdeck,iso,iss4
  110 format(4i5,5x,'nlay, mdeck,  iso,  iss(4)')
      nl=nlay-1
      md=mdeck+1
      if(mdeck.gt.0) then
cc      write(iouts,120) (z(i),d(i),ar(i),br(i),rho(i), i=1,mdeck)
  120 format(5f10.4)
      end if
      if(md.gt.nl) go to 1
cc      write(iouts,130) (z(i),d(i),ar(i),br(i),rho(i), i=md,nl)
  130 format(5f10.4,5x,'refl. zone')
cc    1 write(iouts,140) ar(nlay),br(nlay),rho(nlay)
    1 continue
  140 format(20x,3f10.4,5x,'refl. zone')
cc      write(iouts,150) nent
  150 format(i5,' zahl der seismogramme = zahl der distanzen')
cc      write(iouts,160) (x(i), i=1,nent)
  160 format(7f10.3)
cc      write(iouts,170) vred,tmin,dt
  170 format(2f10.3,f10.4,5x,'v-red - initial time - increment')
      return
      end
c==================================================================
      subroutine inhom(iins)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
c
c     inhom approximates stepwise linear velocity depth function
c     with thin homogeneous layers.applied only if iss(8)=1
c     if iss(12).eq.1 an earth flattenimg approximation is applied
c
      dimension zz(nla),aa(nla),bb(nla),rrho(nla),nhs(nla),qpp(nla),
     1qss(nla),aai(nla),bbi(nla)
      common/blo1/z(0:nla),d(nla),ar(nla),ai(nla),br(nla),bi(nla),
     &rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      complex cv
      qsum=0.
      i=0
  100 i=i+1
      read(iins,9000) zz(i),aa(i),qpp(i),bb(i),qss(i),rrho(i),nhs(i)
 9000 format(6f10.4,i10)
      if(qpp(i).eq.0.) aai(i)=0.
      if(qss(i).eq.0.) qss(i)=qpp(i)
      qsum=qsum+qpp(i)+qss(i)
      if(qss(i).eq.0.) bbi(i)=0.
      if(qpp(i).ne.0.) aai(i)=aa(i)/(2.*qpp(i))
      if(rrho(i).eq.0.)  rrho(i)=0.252+0.3788*aa(i)
      if(bb(i).eq.0.)bb(i)=aa(i)/1.732
      if(bb(i).eq.-1.)bb(i)=0.
      if(qss(i).ne.0.) bbi(i)=bb(i)/(2.*qss(i))
      if(zz(i)+aa(i)+bb(i))100,200,100
  200 npkt = i-1
      write(6,9010)
9010  format(///,1x,'inhomogeneous velocity-depth distribution :'//,4x,
     &'depth              p-velocity   qp     s-velocity   qs      densi
     &ty     number of'/,5x,2hkm,19x,4hkm/s,16x,4hkm/s,14x,7hg/cm**3,6x,
     &6hlayers/)
      write(6,9020) (zz(i),aa(i),qpp(i),bb(i),qss(i),rrho(i),nhs(i),i=1,
     1npkt)
 9020 format(1x,f12.4,10x,5f12.4,i10)
      if(iss(12).ne.1)  go to 700
      write(6,9030)
 9030 format(//,1x,'modification after earth-flattening approximation :'
     &/)
      do  300  i=1,npkt
      y0=3390./(3390.-zz(i))
      zz(i)=3390.*alog(3390./(3390.-zz(i)))
      aa(i)=y0*aa(i)
      aai(i)=y0*aai(i)
      bb(i)=y0*bb(i)
      bbi(i)=y0*bbi(i)
  300 rrho(i)=y0**nrho*rrho(i)
      write(6,9020) (zz(i),aa(i),qpp(i),bb(i),qss(i),rrho(i),nhs(i),i=1,
     1npkt)
  700 j=0
      do 1000 i=2,npkt
      if(nhs(i).eq.0) go to 1000
      k = j+1
      j = j+nhs(i)
      l = i-1
      y0 = (zz(i)-zz(l))/nhs(i)
      y1 = (aa(i)-aa(l))/nhs(i)
      y1i= (aai(i)-aai(l))/nhs(i)
      y2 = (bb(i)-bb(l))/nhs(i)
      y2i= (bbi(i)-bbi(l))/nhs(i)
      y3 = (rrho(i)-rrho(l))/nhs(i)
      y4 = (qpp(i)-qpp(l))/nhs(i)
      y5 = (qss(i)-qss(l))/nhs(i)
      if(nhs(l).eq.0)  go to 701
      ar(k) = aa(l)+y1/2.
      ai(k) = aai(l)+y1i/2.
      br(k) = bb(l)+y2/2.
      bi(k) = bbi(l)+y2i/2.
      rho(k) = rrho(l)+y3/2.
      qp(k) = qpp(l)+y4/2.
      qs(k) = qss(l)+y5/2.
      go to 702
  701 ar(k)=aa(l)
      ai(k)=aai(l)
      br(k)=bb(l)
      bi(k)=bbi(l)
      rho(k)=rrho(l)
      qp(k)=qpp(l)
      qs(k)=qss(l)
  702 z(k) = zz(l)+y0
      d(k) = y0
      if(nhs(i).eq.1) go to 1000
      if(i.lt.npkt.and.nhs(i+1).eq.0)  go to 703
      ar(j) = aa(i)-y1/2.
      ai(j) = aai(i)-y1i/2.
      br(j) = bb(i)-y2/2.
      bi(j) = bbi(i)-y2i/2.
      rho(j) = rrho(i)-y3/2.
      qp(j) = qpp(i)-y4/2.
      qs(j) = qss(i)-y5/2.
      go to 704
  703 ar(j)=aa(i)
      ai(j)=aai(i)
      br(j)=bb(i)
      bi(j)=bbi(i)
      rho(j)=rrho(i)
      qp(j)=qpp(i)
      qs(j)=qss(i)
  704 z(j) = zz(i)
      d(j) = y0
      if(nhs(i).eq.2) go to 1000
      n = k+1
      nn = j-1
      do 800 m=n,nn
      l = m-1
      ar(m) = ar(l)+y1
      ai(m) = ai(l)+y1i
      br(m) = br(l)+y2
      bi(m) = bi(l)+y2i
      rho(m) = rho(l)+y3
      qp(m) = qp(l)+y4
      qs(m) = qs(l)+y5
      z(m) = z(l)+y0
  800 d(m) = y0
 1000 continue
      nlay = j+1
      ar(nlay) = aa(npkt)
      ai(nlay) = aai(npkt)
      br(nlay) = bb(npkt)
      bi(nlay) = bbi(npkt)
      rho(nlay) = rrho(npkt)
      qp(nlay) = qpp(npkt)
      qs(nlay) = qss(npkt)
      return
      end
c==================================================================
      subroutine signa0(t,dt,n,npts,na,s,tsigma)
c==================================================================
c
c     signa0 calculates the integrated form of the source signal
c         sin(n*pi*time/t)-n/(n+2)*sin((n+2)/n*pi*time/t)
c
      dimension s(1)
      complex s
      do 100 i = 1,na
 100  s(i) = (0.0,0.0)
      ns = t/dt + 1.
      fn = n
      fm = (fn + 2.)/fn
      fmq = 1./(fm*fm)
      d = fn*3.14159265/t
      tt=-dt
      a=0.0
      do 200 i=1,ns
      tt=tt+dt
      t1=d*tt
      t2=fm*t1
      t3=3.36*t1
      c=(1.-cos(t1)+fmq*(cos(t2)-1.))/d
      s(na+i)=cmplx(c,0.0)
      if(tsigma.ne.0) s(na+i)=s(na+i)*exp(-tt/tsigma)
      c=abs(sin(t1)-sin(t2)/fm)
      if(c.gt.a) a = c
 200  continue
      a = 1./a
      do 300 i = 1,ns
      j = na + i
 300  s(j) = a*s(j)
      j = na + ns + 1
      a=3.14159265/float(npts-j+1)
      do 400 i = j,npts
  400 s(i)=0.5*(1.+cos(float(i-j+1)*a))*s(ns+na)
      return
      end
c==================================================================
      subroutine signa1(t,dt,s,tsigma)
c==================================================================
c
c     signa1 calculates the signal of a delta impulse
c
      dimension s(1)
      complex s
      pi=3.14159265
      ns=t/dt+1.
      dx=dt/t
      do 100 i=1,ns
      t1=float(i-1)*dt
      s(i)=float(i-1)*dx
      if(tsigma.ne.0) s(i)=s(i)*exp(-t1/tsigma)
 100  continue
      ns1=2*ns-1
      do 200 i=ns+1,ns1
      t1=float(i-1)*dt
      s(i)=float(ns1-i)*dx
 200  if(tsigma.ne.0) s(i)=s(i)*exp(-t1/tsigma)
      write(6,300)
      write(6,301) (s(i),i=1,ns1)
 300  format(/'delta-impuls: '/)
 301  format(10f10.4)
      return
      end
c==================================================================
      subroutine signa3(t,dt,s,tsigma)
c==================================================================
c
c     signa3 calculates the second derivation of the source signal
c         9*m0/16*(1.-cos(pi*t1/t)+1./9*(cos(3*pi*t1/t)-1))
c                      see bruestle and mueller
c
      dimension s(1),ss(8192)
      complex s
      pi=3.14159265
      ns=t/dt+1.
      f1=pi/t
      f2=3.*f1
      f3=f2*f2/16.*10.**5
      do 200 i=1,ns
      t1=float(i-1)*dt
      s(i)=f3*(cos(f1*t1)-cos(f2*t1))
      if(tsigma.ne.0) s(i)=s(i)*exp(-t1/tsigma)
 200  ss(i)=(sin(f1*t1)-sin(f2*t1)/3.)/1.33333
      write(6,300)
      write(6,301) (ss(i),i=1,ns)
 300  format(/'normierte fernfeldverschiebung des signals: '/)
 301  format(10f10.4)
      return
      end
c==================================================================
      subroutine signa5(t,dt,s,tsigma)
c==================================================================
c
c     signa5 - digitized source signal, read in from cards
c
      dimension s(1),x(10)
      complex s
      endata=-9999.
      do 200 j=1,119,8
      jdc=j+7
      read(5,100) (x(i),i=1,8)
      do 210 i=j,jdc
      t1=float(i-1)*dt
      ijk=i-j+1
      if(x(ijk)-endata) 211,220,211
 211  s(i)=cmplx(x(ijk),0.)
      if(tsigma.ne.0) s(i)=s(i)*exp(-t1/tsigma)
 210  continue
 200  continue
 220  write(6,110)
      write(6,100) (s(j),j=1,i)
 110  format(/'digitized source signal: '/)
 100  format(8f10.4)
      return
      end
c
c     sub10
c
c==================================================================
      subroutine cool(n,x,sign)
c==================================================================
c
c     cool executes the fast fourier transform
c
      dimension m(25)
      dimension x(2)
      complex x,hold,wk,q,v,di
      di=(0.0,1.0)
      lx=2**n
      do 1 i=1,n
1     m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=float(k)
      flx=float(lx)
      v=6.2831853*sign*fk/flx*di
      wk=cexp(v)
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q
2     continue
      do 3 i=2,n
      ii=i
      if ( k.lt.m(i)) go to 4
3     k=k-m(i)
4     k=k+m(ii)
      k=0
      do 7 j=1,lx
      if (k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
5     do 6 i=1,n
      ii=i
      if (k.lt.m(i)) go to 7
6     k=k-m(i)
7     k=k+m(ii)
      return
      end
c
c     sub11
c
c==================================================================
      subroutine zzlog(n,m)
c==================================================================
c     do 10 i=1,14 max 16384; now 262144 NS 1/31/06
      do 10 i=1,18
      if(n-2**i) 10,5,10
 5    m=i
      go to 20
 10   continue
 20   return
      end
c==================================================================
      function amax(a,n)
c==================================================================
      dimension a(1)
      amax=a(1)
      do 10 i=2,n
      if(a(i)-amax) 10,10,9
9     amax=a(i)
10    continue
      return
      end
c
c     sub15
c
c==================================================================
      function amin(a,n)
c==================================================================
      dimension a(1)
      amin=a(1)
      do 10 i=2,n
      if(a(i)-amin) 9,10,10
9     amin=a(i)
10    continue
      return
      end
c==================================================================
      subroutine matsh1(x1r,x1i,x2r,x2i,x3r,x3i,y1r,y1i,n)
c==================================================================
c *** ber. von x1/(1-x2*x3) ****
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      dimension x1r(n),x1i(n),x2r(n),x2i(n),x3r(n),x3i(n),y1r(n),y1i(n)
      common/integ3/det(npa),z1r(npa),z1i(npa),z2r(npa),z2i(npa),
     &z3r(npa),z3i(npa)
      do 100 i=1,n
      z1r(i)=1.-x2r(i)*x3r(i)+x2i(i)*x3i(i)
      z1i(i)=-x2i(i)*x3r(i)-x2r(i)*x3i(i)
      det(i)=1./(z1r(i)*z1r(i)+z1i(i)*z1i(i))
      y1r(i)=(x1r(i)*z1r(i)+x1i(i)*z1i(i))*det(i)
 100  y1i(i)=(x1i(i)*z1r(i)-x1r(i)*z1i(i))*det(i)
      return
      end
c==================================================================
      subroutine matsh2(x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,y1r,y1i,
     &n)
c==================================================================
c *** ber. von x1/((1-x2*x3)*(1-x4*x5)) ****
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      dimension x1r(n),x1i(n),x2r(n),x2i(n),x3r(n),x3i(n),x4r(n),x4i(n),
     &x5r(n),x5i(n),y1r(n),y1i(n)
      common/integ3/det(npa),z1r(npa),z1i(npa),z2r(npa),z2i(npa),
     &z3r(npa),z3i(npa)
      do 100 i=1,n
      z1r(i)=1.-x2r(i)*x3r(i)+x2i(i)*x3i(i)
      z1i(i)=-x2i(i)*x3r(i)-x2r(i)*x3i(i)
      z2r(i)=1.-x4r(i)*x5r(i)+x4i(i)*x5i(i)
      z2i(i)=-x4i(i)*x5r(i)-x4r(i)*x5i(i)
      z3r(i)=z1r(i)*z2r(i)-z1i(i)*z2i(i)
      z3i(i)=z1i(i)*z2r(i)+z1r(i)*z2i(i)
      det(i)=1./(z3r(i)*z3r(i)+z3i(i)*z3i(i))
      y1r(i)=(x1r(i)*z3r(i)+x1i(i)*z3i(i))*det(i)
 100  y1i(i)=(x1i(i)*z3r(i)-x1r(i)*z3i(i))*det(i)
      return
      end
c==================================================================
      subroutine rekush(n1,n2,i1,omega,rssr,rssi,it)
c==================================================================
c
c     rekush computes the plane wave reflection coefficients for an
c     incident sh wave for all required slownesses for n layers
c
c     by complex velocities damping is considered
c
c     i1 =  1  -  reflecitivity for an incident upgoing wavefield
c        = -1  -  reflecitivity for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer, if n2=0 - free surface
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch oberste bzw. unterste schicht
c        = 1   -  transmission durch oberste bzw. unterste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex bc(nla),mue(nla),bq1
      real muer,muei
      dimension domr(nla),domi(nla),muer(nla),muei(nla),rssr(npa),
     &rssi(npa),rdsr(npa),rdsi(npa)
      common/blo1/z(0:nla),d(nla),alpr(nla),alpi(nla),betr(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/integ2/bc1r(npa),bc1i(npa),bc2r(npa),bc2i(npa),x1r(npa),
     &x1i(npa),x2r(npa)
      common/integ3/x2i(npa),x3r(npa),x3i(npa),x4r(npa),x4i(npa),
     &x5r(npa),x5i(npa)
      complex cv
      if(n1.eq.n2) then
         do 85 i=1,np
         rssr(i)=0.
 85      rssi(i)=0.
         return
      end if
      n11=n1-i1
      do 100 j=n2,n1,i1
      if(j.eq.0) go to 100
      bc(j)=cmplx(betr(j),beti(j))
      bc(j)=1./(bc(j)*bc(j))
      mue(j)=rho(j)/bc(j)
      muer(j)=real(mue(j))
      muei(j)=aimag(mue(j))
      domr(j)=-2.*d(j)*omega
      domi(j)=-2.*d(j)*omegai
 100  continue
      nn=n1-mdeck
      if(nn.lt.0.and.i1.eq.-1) n11=mdeck-i1
      do 150 j=n2,n11,i1
      if(j.eq.0) then
         do 175 i=1,np
         x5r(i)=0.
         x5i(i)=0.
         rssr(i)=1.
         rssi(i)=0.
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 175     bc2i(i)=aimag(bq1)
         go to 150
      else if(j.eq.n2) then
         do 200 i=1,np
         x5r(i)=0.
         x5i(i)=0.
         bq1=csqrt(bc(j)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc1r(i)=real(bq1)
         bc1i(i)=aimag(bq1)
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 200     bc2i(i)=aimag(bq1)
         do 250 i=1,np
         x1r(i)=bc2r(i)*muer(j+i1)-bc2i(i)*muei(j+i1)
         x1i(i)=bc2i(i)*muer(j+i1)+bc2r(i)*muei(j+i1)
         x2r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x2i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x3r(i)=x1r(i)-x2r(i)
         x3i(i)=x1i(i)-x2i(i)
         x4r(i)=x1r(i)+x2r(i)
         x4i(i)=x1i(i)+x2i(i)
         x1r(i)=1./(x4r(i)*x4r(i)+x4i(i)*x4i(i))
         rdsr(i)=(x3r(i)*x4r(i)+x3i(i)*x4i(i))*x1r(i)
         rdsi(i)=(x3i(i)*x4r(i)-x3r(i)*x4i(i))*x1r(i)
         rssr(i)=rdsr(i)
250      rssi(i)=rdsi(i)
      else
         do 300 i=1,np
         bc1r(i)=bc2r(i)
         bc1i(i)=bc2i(i)
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 300     bc2i(i)=aimag(bq1)
         do 350 i=1,np
         x1r(i)=bc2r(i)*muer(j+i1)-bc2i(i)*muei(j+i1)
         x1i(i)=bc2i(i)*muer(j+i1)+bc2r(i)*muei(j+i1)
         x2r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x2i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x3r(i)=x1r(i)-x2r(i)
         x3i(i)=x1i(i)-x2i(i)
         x4r(i)=x1r(i)+x2r(i)
         x4i(i)=x1i(i)+x2i(i)
         x1r(i)=1./(x4r(i)*x4r(i)+x4i(i)*x4i(i))
         rdsr(i)=(x3r(i)*x4r(i)+x3i(i)*x4i(i))*x1r(i)
         rdsi(i)=(x3i(i)*x4r(i)-x3r(i)*x4i(i))*x1r(i)
         x1i(i)=domr(j)*bc1r(i)-domi(j)*bc1i(i)
         x1r(i)=-domr(j)*bc1i(i)-domi(j)*bc1r(j)
         x2r(i)=exp(x1r(i))*cos(x1i(i))
         x2i(i)=exp(x1r(i))*sin(x1i(i))
         x1r(i)=x2r(i)*rssr(i)-x2i(i)*rssi(i)
         x1i(i)=x2i(i)*rssr(i)+x2r(i)*rssi(i)
         x2r(i)=1.+rdsr(i)*x1r(i)-rdsi(i)*x1i(i)
         x2i(i)=rdsi(i)*x1r(i)+rdsr(i)*x1i(i)
         x3r(i)=1./(x2r(i)*x2r(i)+x2i(i)*x2i(i))
         x4r(i)=rdsr(i)+x1r(i)
         x4i(i)=rdsi(i)+x1i(i)
         rssr(i)=(x4r(i)*x2r(i)+x4i(i)*x2i(i))*x3r(i)
350      rssi(i)=(x4i(i)*x2r(i)-x4r(i)*x2i(i))*x3r(i)
      end if
150   continue
      if(nn.lt.0.and.i1.eq.-1) then
         n11=n1-i1
         do 500 j=mdeck,n11,i1
         do 525 i=1,np
         bc1r(i)=bc2r(i)
         bc1i(i)=bc2i(i)
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 525     bc2i(i)=aimag(bq1)
         do 550 i=1,np
         x5i(i)=domr(j)*bc1r(i)-domi(j)*bc1i(i)+x5i(i)
         x5r(i)=-domr(j)*bc1i(i)-domi(j)*bc1r(j)+x5r(i)
         x1r(i)=bc2r(i)*muer(j+i1)-bc2i(i)*muei(j+i1)
         x1i(i)=bc2i(i)*muer(j+i1)+bc2r(i)*muei(j+i1)
         x2r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x2i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x3r(i)=x1r(i)+x2r(i)
         x3i(i)=x1i(i)+x2i(i)
         x4r(i)=x3r(i)*x3r(i)-x3i(i)*x3i(i)
         x4i(i)=2.*x3r(i)*x3i(i)
         x3r(i)=x1r(i)*x2r(i)-x1i(i)*x2i(i)
         x3i(i)=x1i(i)*x2r(i)+x1r(i)*x2i(i)
         x2r(i)=1./(x4r(i)*x4r(i)+x4i(i)*x4i(i))
         rdsr(i)=4.*(x4r(i)*x3r(i)+x4i(i)*x3i(i))*x2r(i)
         rdsi(i)=4.*(x4i(i)*x3r(i)-x4r(i)*x3i(i))*x2r(i)
         x1r(i)=rdsr(i)*rssr(i)-rdsi(i)*rssi(i)
         x1i(i)=rdsi(i)*rssr(i)+rdsr(i)*rssi(i)
         rssr(i)=x1r(i)
 550     rssi(i)=x1i(i)
 500     continue
      end if
      if(it.eq.1) then
c
c     transmission durch oberste bzw. unterste schicht
c
         do 400 i=1,np
         x5i(i)=domr(n1)*bc2r(i)-domi(n1)*bc2i(i)+x5i(i)
 400     x5r(i)=-domr(n1)*bc2i(i)-domi(n1)*bc2r(j)+x5r(i)
      end if
      do 405 i=1,np
      x2r(i)=exp(x5r(i))*cos(x5i(i))
      x2i(i)=exp(x5r(i))*sin(x5i(i))
      x1r(i)=x2r(i)*rssr(i)-x2i(i)*rssi(i)
      x1i(i)=x2i(i)*rssr(i)+x2r(i)*rssi(i)
      rssr(i)=x1r(i)
 405  rssi(i)=x1i(i)
      return
      end
c==================================================================
      subroutine tramsh(n1,n2,i1,omega,tssr,tssi,rssr,rssi,it)
c==================================================================
c
c         tramsh computes the transmissivity for an icnident
c     sh wave for all required slownesses and one frequency omega
c
c     by complex velocities damping is considered
c
c     i1 =  1  -  reflecitivity for an incident upgoing wavefield
c        = -1  -  reflecitivity for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer, if n2=0 - free surface
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch oberste bzw. unterste schicht
c        = 2   -  transmission durch oberste und unterste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex bc(nla),mue(nla),bq1
      real muer,muei
      dimension domr(nla),domi(nla),muer(nla),muei(nla),rssr(npa),
     &rssi(npa),rdsr(npa),rdsi(npa),tssr(npa),tssi(npa),x6r(npa)
      common/blo1/z(0:nla),d(nla),alpr(nla),alpi(nla),betr(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/integ2/bc1r(npa),bc1i(npa),bc2r(npa),bc2i(npa),x1r(npa),
     &x1i(npa),x2r(npa)
      common/integ3/x2i(npa),x3r(npa),x3i(npa),x4r(npa),x4i(npa),
     &x5r(npa),x5i(npa)
      complex cv
      if(n1.eq.n2) then
         do 90 i=1,np
         tssr(i)=1.
 90      tssi(i)=0.
         go to 11
      end if
      n11=n1-i1
      do 100 j=n2,n1,i1
      if(j.eq.0) go to 100
      bc(j)=cmplx(betr(j),beti(j))
      bc(j)=1./(bc(j)*bc(j))
      mue(j)=rho(j)/bc(j)
      muer(j)=real(mue(j))
      muei(j)=aimag(mue(j))
      domr(j)=-2.*d(j)*omega
      domi(j)=-2.*d(j)*omegai
 100  continue
      do 150 j=n2,n11,i1
      if(j.eq.0) then
         do 175 i=1,np
         rssr(i)=1.
         rssi(i)=0.
         tssr(i)=2.
         tssi(i)=0.
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 175     bc2i(i)=aimag(bq1)
         go to 150
      else if(j.eq.n2) then
         do 200 i=1,np
         bq1=csqrt(bc(j)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc1r(i)=real(bq1)
         bc1i(i)=aimag(bq1)
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 200     bc2i(i)=aimag(bq1)
         do 250 i=1,np
         x1r(i)=bc2r(i)*muer(j+i1)-bc2i(i)*muei(j+i1)
         x1i(i)=bc2i(i)*muer(j+i1)+bc2r(i)*muei(j+i1)
         x2r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x2i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x3r(i)=x1r(i)-x2r(i)
         x3i(i)=x1i(i)-x2i(i)
         x4r(i)=x1r(i)+x2r(i)
         x4i(i)=x1i(i)+x2i(i)
         x6r(i)=1./(x4r(i)*x4r(i)+x4i(i)*x4i(i))
         rdsr(i)=(x3r(i)*x4r(i)+x3i(i)*x4i(i))*x6r(i)
         rdsi(i)=(x3i(i)*x4r(i)-x3r(i)*x4i(i))*x6r(i)
         tssr(i)=2.*(x1r(i)*x4r(i)+x1i(i)*x4i(i))*x6r(i)
         tssi(i)=2.*(x1i(i)*x4r(i)-x1r(i)*x4i(i))*x6r(i)
         rssr(i)=rdsr(i)
250      rssi(i)=rdsi(i)
      else
         do 300 i=1,np
         bc1r(i)=bc2r(i)
         bc1i(i)=bc2i(i)
         bq1=csqrt(bc(j+i1)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 300     bc2i(i)=aimag(bq1)
         do 350 i=1,np
         x1r(i)=bc2r(i)*muer(j+i1)-bc2i(i)*muei(j+i1)
         x1i(i)=bc2i(i)*muer(j+i1)+bc2r(i)*muei(j+i1)
         x2r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x2i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x3r(i)=x1r(i)-x2r(i)
         x3i(i)=x1i(i)-x2i(i)
         x4r(i)=x1r(i)+x2r(i)
         x4i(i)=x1i(i)+x2i(i)
         x6r(i)=1./(x4r(i)*x4r(i)+x4i(i)*x4i(i))
         rdsr(i)=(x3r(i)*x4r(i)+x3i(i)*x4i(i))*x6r(i)
         rdsi(i)=(x3i(i)*x4r(i)-x3r(i)*x4i(i))*x6r(i)
         x5r(i)=2.*(x1r(i)*x4r(i)+x1i(i)*x4i(i))*x6r(i)
         x5i(i)=2.*(x1i(i)*x4r(i)-x1r(i)*x4i(i))*x6r(i)
         x1i(i)=domr(j)*bc1r(i)-domi(j)*bc1i(i)
         x1r(i)=-domr(j)*bc1i(i)-domi(j)*bc1r(j)
         x2r(i)=exp(x1r(i))*cos(x1i(i))
         x2i(i)=exp(x1r(i))*sin(x1i(i))
         x3r(i)=exp(x1r(i)*0.5)*cos(x1i(i)*0.5)
         x3i(i)=exp(x1r(i)*0.5)*sin(x1i(i)*0.5)
         x1r(i)=x2r(i)*rssr(i)-x2i(i)*rssi(i)
         x1i(i)=x2i(i)*rssr(i)+x2r(i)*rssi(i)
         x4r(i)=tssr(i)
         x4i(i)=tssi(i)
         tssr(i)=x4r(i)*x3r(i)-x4i(i)*x3i(i)
         tssi(i)=x4i(i)*x3r(i)+x4r(i)*x3i(i)
         x2r(i)=1.+rdsr(i)*x1r(i)-rdsi(i)*x1i(i)
         x2i(i)=rdsi(i)*x1r(i)+rdsr(i)*x1i(i)
         x3r(i)=1./(x2r(i)*x2r(i)+x2i(i)*x2i(i))
         x4r(i)=rdsr(i)+x1r(i)
         x4i(i)=rdsi(i)+x1i(i)
         rssr(i)=(x4r(i)*x2r(i)+x4i(i)*x2i(i))*x3r(i)
         rssi(i)=(x4i(i)*x2r(i)-x4r(i)*x2i(i))*x3r(i)
         x1r(i)=(tssr(i)*x2r(i)+tssi(i)*x2i(i))*x3r(i)
         x1i(i)=(tssi(i)*x2r(i)-tssr(i)*x2i(i))*x3r(i)
         tssr(i)=x5r(i)*x1r(i)-x5i(i)*x1i(i)
350      tssi(i)=x5i(i)*x1r(i)+x5r(i)*x1i(i)
      end if
150   continue
 11   if(it.eq.2) then
c
c     transmission durch oberste und unterste schicht
c
         do 400 i=1,np
         bq1=csqrt(bc(n2)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc1r(i)=real(bq1)
         bc1i(i)=aimag(bq1)
         x5i(i)=domr(n2)*bc1r(i)-domi(n2)*bc1i(i)
         x5r(i)=-domr(n2)*bc1i(i)-domi(n2)*bc1r(j)
         x5i(i)=(domr(n1)*bc2r(i)-domi(n1)*bc2i(i)+x5i(i))*0.5
         x5r(i)=(-domr(n1)*bc2i(i)-domi(n1)*bc2r(j)+x5r(i))*0.5
         x2r(i)=exp(x5r(i))*cos(x5i(i))
         x2i(i)=exp(x5r(i))*sin(x5i(i))
         x1r(i)=x2r(i)*tssr(i)-x2i(i)*tssi(i)
         x1i(i)=x2i(i)*tssr(i)+x2r(i)*tssi(i)
         tssr(i)=x1r(i)
 400     tssi(i)=x1i(i)
      end if
      return
      end
c==================================================================
      subroutine transh(n1,n2,i1,omega,tssr,tssi,it)
c==================================================================
c
c     transh computes the transmission for an
c     incident sh wave for all required slownesses for n layers
c
c     by complex velocities damping is considered
c
c     i1 =  1  -  reflecitivity for an incident upgoing wavefield
c        = -1  -  reflecitivity for an incident downgoing wavefield
c
c     n1       -  first layer
c                 if i1=1  n1=lowermost layer
c                 if i1=-1 n1=uppermost layer
c
c     n2       -  last layer
c                 if i1=1  n2=uppermost layer, if n2=0 - free surface
c                 if i1=-1 n2=lowermost layer
c
c     it = 0   -  keine transmission durch oberste bzw. unterste schicht
c        = 1   -  transmission durch oberste bzw. unterste schicht
c        = 2   -  transmission durch oberste und unterste schicht
c
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      complex bc(nla),mue(nla),bq1
      real muer,muei
      dimension domr(nla),domi(nla),muer(nla),muei(nla),tssr(npa),
     &tssi(npa),rdsr(npa),rdsi(npa)
      common/blo1/z(0:nla),d(nla),alpr(nla),alpi(nla),betr(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common/integ2/bc1r(npa),bc1i(npa),bc2r(npa),bc2i(npa),x1r(npa),
     &x1i(npa),x2r(npa)
      common/integ3/x2i(npa),x3r(npa),x3i(npa),x4r(npa),x4i(npa),
     &x5r(npa),x5i(npa)
      complex cv
      if(n1.eq.n2.and.it.ne.2) then
         do 85 i=1,np
         tssr(i)=0.
 85      tssi(i)=0.
         return
      end if
      n22=n2+i1
      in=-i1
      do 100 j=n1,n2,in
      bc(j)=cmplx(betr(j),beti(j))
      bc(j)=1./(bc(j)*bc(j))
      mue(j)=rho(j)/bc(j)
      muer(j)=real(mue(j))
      muei(j)=aimag(mue(j))
      domr(j)=-d(j)*omega
 100  domi(j)=-d(j)*omegai
      if(n1.eq.n2) go to 500
      do 150 j=n1,n22,in
      if(j.eq.n1) then
         do 200 i=1,np
         bq1=csqrt(bc(j)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc1r(i)=real(bq1)
         bc1i(i)=aimag(bq1)
         bq1=csqrt(bc(j+in)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 200     bc2i(i)=aimag(bq1)
         do 250 i=1,np
         x1r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x1i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x2r(i)=bc2r(i)*muer(j+in)-bc2i(i)*muei(j+in)
         x2i(i)=bc2i(i)*muer(j+in)+bc2r(i)*muei(j+in)
         x3r(i)=x1r(i)+x2r(i)
         x3i(i)=x1i(i)+x2i(i)
         x4r(i)=1./(x3r(i)*x3r(i)+x3i(i)*x3i(i))
         rdsr(i)=2.*(x1r(i)*x3r(i)+x1i(i)*x3i(i))*x4r(i)
         rdsi(i)=2.*(x1i(i)*x3r(i)-x1r(i)*x3i(i))*x4r(i)
         tssr(i)=rdsr(i)
         tssi(i)=rdsi(i)
         x5r(i)=0.
 250     x5i(i)=0.
         if(it.ne.2.and.in.eq.-1) go to 150
         if(it.eq.0) go to 150
         do 260 i=1,np
         x5i(i)=domr(j)*bc1r(i)-domi(j)*bc1i(i)
 260     x5r(i)=-domr(j)*bc1i(i)-domi(j)*bc1r(j)
         go to 150
      else
         do 300 i=1,np
         bc1r(i)=bc2r(i)
         bc1i(i)=bc2i(i)
         bq1=csqrt(bc(j+in)-p2(i))
         if(real(bq1).eq.0.) bq1=-bq1
         bc2r(i)=real(bq1)
 300     bc2i(i)=aimag(bq1)
         do 350 i=1,np
         x5i(i)=domr(j)*bc1r(i)-domi(j)*bc1i(i)+x5i(i)
         x5r(i)=-domr(j)*bc1i(i)-domi(j)*bc1r(j)+x5r(i)
         x1r(i)=bc1r(i)*muer(j)-bc1i(i)*muei(j)
         x1i(i)=bc1i(i)*muer(j)+bc1r(i)*muei(j)
         x2r(i)=bc2r(i)*muer(j+in)-bc2i(i)*muei(j+in)
         x2i(i)=bc2i(i)*muer(j+in)+bc2r(i)*muei(j+in)
         x3r(i)=x1r(i)+x2r(i)
         x3i(i)=x1i(i)+x2i(i)
         x4r(i)=1./(x3r(i)*x3r(i)+x3i(i)*x3i(i))
         rdsr(i)=2.*(x1r(i)*x3r(i)+x1i(i)*x3i(i))*x4r(i)
         rdsi(i)=2.*(x1i(i)*x3r(i)-x1r(i)*x3i(i))*x4r(i)
         x1r(i)=rdsr(i)*tssr(i)-rdsi(i)*tssi(i)
         x1i(i)=rdsi(i)*tssr(i)+rdsr(i)*tssi(i)
         tssr(i)=x1r(i)
 350     tssi(i)=x1i(i)
      end if
 150  continue
      if(it.ne.2.and.in.eq.1) go to 400
      if(it.eq.0) go to 400
      do 450 i=1,np
      x5i(i)=domr(n2)*bc2r(i)-domi(n2)*bc2i(i)+x5i(i)
 450  x5r(i)=-domr(n2)*bc2i(i)-domi(n2)*bc2r(j)+x5r(i)
 400  do 405 i=1,np
      x2r(i)=exp(x5r(i))*cos(x5i(i))
      x2i(i)=exp(x5r(i))*sin(x5i(i))
      x1r(i)=x2r(i)*tssr(i)-x2i(i)*tssi(i)
      x1i(i)=x2i(i)*tssr(i)+x2r(i)*tssi(i)
      tssr(i)=x1r(i)
 405  tssi(i)=x1i(i)
      return
 500  do 550 i=1,np
      bq1=csqrt(bc(n1)-p2(i))
      if(real(bq1).eq.0.) bq1=-bq1
      bc1r(i)=real(bq1)
 550  bc1i(i)=aimag(bq1)
      do 560 i=1,np
      x5i(i)=domr(n1)*bc1r(i)-domi(n1)*bc1i(i)
      x5r(i)=-domr(n1)*bc1i(i)-domi(n1)*bc1r(j)
      tssr(i)=exp(x5r(i))*cos(x5i(i))
 560  tssi(i)=exp(x5r(i))*sin(x5i(i))
      return
      end
c==================================================================
      subroutine phas1(a1,a2,c1,c2,e1,e2,e3,e4,f1,f2,f3,f4,omega,n)
c==================================================================
      dimension e1(n),e2(n),e3(n),e4(n),f1(n),f2(n),f3(n),f4(n),omega(n)
      complex c1,c2
      xr=real(c1)
      xi=aimag(c1)
      xi=abs(xi)
      yr=real(c2)
      yi=aimag(c2)
      yi=abs(yi)
      x1=a1*xr
      x2=a1*xi
      x3=a1*yr
      x4=a1*yi
      y1=a2*x1
      y2=a2*x2
      y3=a2*x3
      y4=a2*x4
      do 100 i=1,n
      e1(i)=x1*omega(i)
      e2(i)=x2*omega(i)
      e3(i)=x3*omega(i)
      e4(i)=x4*omega(i)
      f1(i)=y1
      f2(i)=y2
      f3(i)=y3
      f4(i)=y4
 100  continue
      return
      end
c==================================================================
      subroutine phas3(a1,a2,c1,c2,e1,e2,e3,e4,f1,f2,f3,f4,omega,n)
c==================================================================
      dimension e1(n),e2(n),e3(n),e4(n),f1(n),f2(n),f3(n),f4(n),omega(n)
      complex c1,c2
      xr=real(c1)
      xi=aimag(c1)
      yr=real(c2)
      yi=aimag(c2)
      x1=a1*xr
      x2=a1*xi
      x3=a1*yr
      x4=a1*yi
      y1=a2*x1
      y2=a2*x2
      y3=a2*x3
      y4=a2*x4
      do 100 i=1,n
      e1(i)=x1*omega(i)
      e2(i)=x2*omega(i)
      e3(i)=x3*omega(i)
      e4(i)=x4*omega(i)
      f1(i)=y1
      f2(i)=y2
      f3(i)=y3
      f4(i)=y4
 100  continue
      return
      end
c==================================================================
      subroutine phas2(a1,a2,e1,e2,e3,e4,f1,f2,f3,f4,p2i,omega,n,ilay)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      dimension e1(n),e2(n),e3(n),e4(n),f1(n),f2(n),f3(n),f4(n),omega(n)
      common/blo1/z(0:nla),d(nla),alp(nla),alpi(nla),bet(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      complex cv
      dimension alr(nfr),ali(nfr),alr1(nfr),ali1(nfr)
      dimension ber(nfr),bei(nfr),ber1(nfr),bei1(nfr)
      dimension alb(nfr),alh(nfr)
      dimension beb(nfr),beh(nfr)
      qp1=1./qp(ilay)
      qs1=1./qs(ilay)
      do 100 i=1,n
      alr(i)=alp(ilay)*(1.+cvr(i)*qp1)
      ali(i)=alp(ilay)*cvi(i)*qp1
      ber(i)=bet(ilay)*(1.+cvr(i)*qs1)
      bei(i)=bet(ilay)*cvi(i)*qs1
      alr1(i)=alr(i)*alr(i)-ali(i)*ali(i)
      ali1(i)=2.*alr(i)*ali(i)
      ber1(i)=ber(i)*ber(i)-bei(i)*bei(i)
      bei1(i)=2.*ber(i)*bei(i)
      alr(i)=1./(alr1(i)*alr1(i)+ali1(i)*ali1(i))
      ber(i)=1./(ber1(i)*ber1(i)+bei1(i)*bei1(i))
      alr1(i)=alr(i)*alr1(i)-p2i
      ali1(i)=-alr(i)*ali1(i)
      ber1(i)=ber(i)*ber1(i)-p2i
      bei1(i)=-ber(i)*bei1(i)
      alb(i)=sqrt(alr1(i)*alr1(i)+ali1(i)*ali1(i))
      alh(i)=atan(ali1(i)/alr1(i))*0.5
      beb(i)=sqrt(ber1(i)*ber1(i)+bei1(i)*bei1(i))
      beh(i)=atan(bei1(i)/ber1(i))*0.5
      alr(i)=sqrt(alb(i))
      ber(i)=sqrt(beb(i))
 100  continue
      do 110 i=1,n
      if(alr1(i).lt.0.) alh(i)=alh(i)-0.5*pi
      if(ber1(i).lt.0.) beh(i)=beh(i)-0.5*pi
 110  continue
      do 120 i=1,n
      f1(i)=alr(i)*cos(alh(i))
      f2(i)=alr(i)*sin(alh(i))
      f2(i)=abs(f2(i))
      f3(i)=ber(i)*cos(beh(i))
      f4(i)=ber(i)*sin(beh(i))
      f4(i)=abs(f4(i))
 120  continue
      a12=a1*a2
      do 200 i=1,n
      e1(i)=a1*f1(i)*omega(i)
      e2(i)=a1*f2(i)*omega(i)
      e3(i)=a1*f3(i)*omega(i)
      e4(i)=a1*f4(i)*omega(i)
      f1(i)=a12*f1(i)
      f2(i)=a12*f2(i)
      f3(i)=a12*f3(i)
      f4(i)=a12*f4(i)
200   continue
c     complex al(nfr),be(nfr),am(nfr),bm(nfr),cv
c     do 100 i=1,n
c     al(i)=alp(ilay)*(1.+cv(i)/qp(ilay))
c     be(i)=bet(ilay)*(1.+cv(i)/qs(ilay))
c     am(i)=csqrt(1./(al(i)*al(i))-p2i)
c     bm(i)=csqrt(1./(be(i)*be(i))-p2i)
c     if(real(am(i)).eq.0.) am(i)=-am(i)
c     if(real(bm(i)).eq.0.) bm(i)=-bm(i)
c00   continue
      return
      end
c==================================================================
      subroutine phas4(a1,a2,e1,e2,e3,e4,f1,f2,f3,f4,p2i,omega,n,ilay)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000)
      dimension e1(n),e2(n),e3(n),e4(n),f1(n),f2(n),f3(n),f4(n),omega(n)
      common/blo1/z(0:nla),d(nla),alp(nla),alpi(nla),bet(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      complex al(nfr),be(nfr),am(nfr),bm(nfr),cv
      do 100 i=1,n
      al(i)=alp(ilay)*(1.+cv(i)/qp(ilay))
      be(i)=bet(ilay)*(1.+cv(i)/qs(ilay))
      am(i)=csqrt(1./(al(i)*al(i))-p2i)
      bm(i)=csqrt(1./(be(i)*be(i))-p2i)
      if(real(am(i)).eq.0.) am(i)=-am(i)
      if(real(bm(i)).eq.0.) bm(i)=-bm(i)
      f1(i)=real(am(i))
      f2(i)=aimag(am(i))
      f3(i)=real(bm(i))
      f4(i)=aimag(bm(i))
100   continue
      a12=a1*a2
      do 200 i=1,n
      e1(i)=a1*f1(i)*omega(i)
      e2(i)=a1*f2(i)*omega(i)
      e3(i)=a1*f3(i)*omega(i)
      e4(i)=a1*f4(i)*omega(i)
      f1(i)=a12*f1(i)
      f2(i)=a12*f2(i)
      f3(i)=a12*f3(i)
      f4(i)=a12*f4(i)
200   continue
      return
      end
c==================================================================
      subroutine shseis(ioutsh)
c==================================================================
c      parameter(nfr=1024,npa=2500,nla=1000,ndis=200)
c NS_CHANGE
      parameter(nfr=131072,npa=5000,nla=1000,ndis=200)
      common/bsh1/xs(100),ys(100),es(100),rr(100),dz(100),vv(100),
     &ts(100),zs(100),az(1000),alpr(nla),betr(nla),ztrue(100)
      common/blo1/z(0:nla),d(nla),alp(nla),alpi(nla),bet(nla),
     &beti(nla),rho(nla),qp(nla),qs(nla),qsum,nlay,mdeck,iso,ire,
     &p(npa),p2(npa),p3(npa),pi,pi2,x(1000),vred,tmin,dt,npts,dso,nrho,
     &dre,np,pi4,cv(nfr),cvr(nfr),cvi(nfr)
      common/blo2/nfreq,iss(25),omegai,omegar,omegac(nfr),phasc(nfr),
     &sphasc(nfr),cphasc(nfr)
      common x1r(nfr),x1i(nfr),x2r(nfr),x2i(nfr),x3r(nfr),x3i(nfr),
     &x4r(nfr),x4i(nfr),x5r(nfr),x5i(nfr),z1r(nfr),z1i(nfr),z2r(nfr),
     &z2i(nfr),z11r(nfr),z11i(nfr),z12r(nfr),z12i(nfr),z22r(nfr),
     &z22i(nfr),detr(nfr),deti(nfr),dett(nfr),w1r(nfr),w1i(nfr),
     &w2r(nfr),w2i(nfr),w3r(nfr),w3i(nfr),r1r(nfr),r1i(nfr),r2r(nfr),
     &r2i(nfr),r3r(nfr),r3i(nfr),omega(nfr),xsum1r(nfr,ndis),
     &xsum1i(nfr,ndis),xsum2r(nfr,ndis),xsum2i(nfr,ndis)
      common/integ1/tshr(npa),tshi(npa),rshd1r(npa),rshd1i(npa),
     &rshd2r(npa),rshd2i(npa),rshu1r(npa),rshu1i(npa),rshu2r(npa),
     &rshu2i(npa)
      common/integ2/xx1r(npa),xx1i(npa),yy1r(npa),yy1i(npa),xx2r(npa),
     &xx2i(npa),xx6i(npa)
      common/integ3/bmr(npa),bmi(npa),bmb(npa),uashr(npa),
     &uashi(npa),ubshr(npa),ubshi(npa)
      common/bsh2/omega2(nfr),phase(0:nfr),
c     &s(8193),wavel(4096),balla(8192),iballa(8192),window(npa),
c     &phase1(nfr),sc1(4096)
c     &s(32769),wavel(16384),balla(32768),iballa(32768),window(npa),
c     &phase1(nfr),sc1(16384)
c     &s(65537),wavel(32768),balla(65536),iballa(65536),window(npa),
c     &phase1(nfr),sc1(32768)
     &s(131073),wavel(65536),balla(131072),iballa(131072),window(npa),
     &phase1(nfr),sc1(65536)
      common/bsh3/cwil,cwir,c1,c2,dnue,dp,ein,einh,fo,fu,fwil,fwir,
     &inpt,io,iu,kom,na,nent,nh,tsigma
      common/bsh4/kappa1,kappa2,kappa3,lamda1,lamda2,m11,m12,m13,m22,
     &m23,m33,lamda
      complex ein,einh,wavel,s,am,bm,al(nla),be(nla),bmiso,bmiso2,sc1,
     &cv
      real kappa1,kappa2,kappa3,lamda1,lamda2,m11,m12,m13,m22,m23,m33,
     &lamda
      dimension su1r(npa),su1i(npa),sd1r(npa),sd1i(npa),sp(npa)
      character*7 filename
      d11x=-1./(pi2*rho(iso)*sqrt(pi2))
      be(iso)=cmplx(betr(iso),beti(iso))
      bmiso=1./(be(iso)*be(iso))
      bmiso2=bmiso/2.
      sp1=window(1)*sqrt(p(1))
      do 4000 jf=1,nfreq
      omeg=omega(jf)
      if(ire.ne.0) then
         diffr=dre-z(ire)
         be(ire)=cmplx(betr(ire),beti(ire))
         do 1001 i=1,np
         bm=csqrt(1./(be(ire)*be(ire))-p2(i))
         if(real(bm).eq.0.) bm=-bm
         bm1r=real(bm)
         bm1i=aimag(bm)
         bmdr=bm1r*diffr
         bmdi=bm1i*diffr
         bomr=-omega(jf)*bmdi-omegai*bmdr
         bomi=omega(jf)*bmdr-omegai*bmdi
         ebomp=exp(bomr)
         ebomm=exp(-bomr)
         sd1r(i)=ebomp*cos(bomi)
         sd1i(i)=ebomp*sin(bomi)
         su1r(i)=ebomm*cos(bomi)
         su1i(i)=-ebomm*sin(bomi)
 1001    continue
      else
         do 1002 i=1,np
         sd1r(i)=1.
         sd1i(i)=0.
         su1r(i)=1.
         su1i(i)=0.
 1002    continue
      end if
      id=-1
      ih=1
      xx1=1.
      iso1=iso-1
      if(iso1.lt.ire) go to 4200
      if(iso1.gt.ire) go to 4100
      if(iso1.eq.ire.and.ire.ne.0) go to 4100
c***********************************************************************
c     source and receiver are situated at the free surface
c***********************************************************************
c     turs=1.    und   rdrs=0.
c
      do 4015 i=1,np
      tshr(i)=2.
 4015 tshi(i)=0.
      if(mdeck.lt.iso) then
c  ber. der reflmatrix rdsl
      call rekush(iso,nlay,id,omeg,rshd2r,rshd2i,1)
c  ber. der reflmatrix rusm (rusf,falls mdeck=0)
      call rekush(iso,mdeck,ih,omeg,rshu2r,rshu2i,0)
c  ber. von a1=tusr*(1-rdsl*rusf)**-1
      call matsh1(tshr,tshi,rshd2r,rshd2i,rshu2r,rshu2i,xx1r,xx1i,
     &np)
      go to 4330
      else
c  ber. der reflmatrix rdsl=tums*rdml*tdsm
      call rekush(iso,nlay,id,omeg,rshd2r,rshd2i,1)
c  ber. von y1=tusr*rdsl, xx1=0.
      do 4050 i=1,np
      xx2r(i)=rshd2r(i)*sd1r(i)-rshd2i(i)*sd1i(i)
      xx2i(i)=rshd2i(i)*sd1r(i)+rshd2r(i)*sd1i(i)
      yy1r(i)=tshr(i)*xx2r(i)-tshi(i)*xx2i(i)
4050  yy1i(i)=tshi(i)*xx2r(i)+tshr(i)*xx2i(i)
      xx1=0.
      go to 4400
      end if
c***********************************************************************
c     the source is situated under the receiver
c     case - receiver at the free surface - is included
c***********************************************************************
c  ber. der transmissivitaetsmatrix tusr
 4100 if(iss(21).eq.1) then
         call tramsh(iso,ire,ih,omeg,tshr,tshi,rshu2r,rshu2i,0)
      else if(ire.ne.0) then
         call transh(iso,ire,ih,omeg,tshr,tshi,0)
      else
         call transh(iso,1,ih,omeg,tshr,tshi,1)
         do 4115 i=1,np
         tshr(i)=2.*tshr(i)
 4115    tshi(i)=2.*tshi(i)
      end if
      if(mdeck.lt.ire) then
c  ber. der reflmatrix rdrs
      call rekush(ire,iso,id,omeg,rshd1r,rshd1i,0)
c  ber. der reflmatrix rurm (rurf,falls mdeck=0)
      call rekush(ire,mdeck,ih,omeg,rshu1r,rshu1i,1)
c  ber. der reflmatrix rdsl
      call rekush(iso,nlay,id,omeg,rshd2r,rshd2i,1)
c   ber. der reflmatrix rusm (rusf,falls mdeck=0)
      call rekush(iso,mdeck,ih,omeg,rshu2r,rshu2i,0)
c     ber. von (1-rdrs*rurf)**-1*tusr*(1-rdsl*rusm)**-1
      call matsh2(tshr,tshi,rshd1r,rshd1i,rshu1r,rshu1i,rshd2r,rshd2i,
     &rshu2r,rshu2i,yy1r,yy1i,np)
      go to 4320
      else if(mdeck.lt.iso) then
c*****rurf=0  ---  rdrs braucht nicht berechnet zu werden
c*****mdeck=ire ist eingeschlossen
c  ber. der reflmatrix rdsl
      call rekush(iso,nlay,id,omeg,rshd2r,rshd2i,1)
c     ber. der reflmatrix rusm
      if(iss(21).ne.1)call rekush(iso,mdeck,ih,omeg,rshu2r,rshu2i,0)
c  ber. von a1=tusr*(1-rdsl*rusm)**-1
      call matsh1(tshr,tshi,rshd2r,rshd2i,rshu2r,rshu2i,xx1r,xx1i,np)
      go to 4330
      else
c*****rurf=0   rdrs=0   rusf=0
c  ber. von rdsl=tums*rdml*tdsm
      call rekush(iso,nlay,id,omeg,rshd2r,rshd2i,1)
      do 4150 i=1,np
      xx1r(i)=tshr(i)
 4150 xx1i(i)=tshi(i)
      go to 4330
      end if
c***********************************************************************
c     the source is situated above the receiver
c     case - source at the free surface - is included
c***********************************************************************
c  ber. der transmissivtaetsmatrix tdsr
 4200 if(iss(21).eq.1) then
         call tramsh(iso,ire,id,omeg,tshr,tshi,rshu2r,rshu2i,2)
      else
         call transh(iso,ire,id,omeg,tshr,tshi,2)
      end if
      if(mdeck.lt.iso) then
c  ber. der reflmatrix rurs
      call rekush(ire,iso,ih,omeg,rshd1r,rshd1i,1)
c  ber. der reflmatrix rdrl
      call rekush(ire,nlay,id,omeg,rshu1r,rshu1i,0)
c  ber. der reflmatrix rusm (rusf,falls mdeck=0)
      call rekush(iso,mdeck,ih,omeg,rshd2r,rshd2i,0)
c  ber. der reflmatrix rdsl
      call rekush(iso,nlay,id,omeg,rshu2r,rshu2i,1)
c  ber. von a1=(1-rurs*rdrl)**-1*tdsr*(1-rusf*rdsl)**-1
      call matsh2(tshr,tshi,rshd1r,rshd1i,rshu1r,rshu1i,rshd2r,rshd2i,
     &rshu2r,rshu2i,xx1r,xx1i,np)
      go to 4310
      else if(mdeck.lt.ire) then
c*****rusf=0 -- rdsl ist unnoetig geworden
c*****mdeck=iso ist eingeschlossen, auch iso=1
c  ber. der reflmatrix rurm
      call rekush(ire,iso,ih,omeg,rshd1r,rshd1i,1)
c  ber. der reflmatrix rdrl
      call rekush(ire,nlay,id,omeg,rshu1r,rshu1i,0)
c     ber. von a1=(1-rurm*rdrl)**-1*tdsr
      call matsh1(tshr,tshi,rshd1r,rshd1i,rshu1r,rshu1i,xx1r,xx1i,np)
      xx1=0.
      go to 4340
      else
c*****rusf=0    rurs=0
c     ber. von rdrl=tumr*rdml*tdrm
      it=0
      if(mdeck.ne.ire) it=1
      call rekush(ire,nlay,id,omeg,rshu1r,rshu1i,it)
      do 4250 i=1,np
      xx1r(i)=tshr(i)
 4250 xx1i(i)=tshi(i)
      xx1=0.
      go to 4340
      end if
c***********************************************************************
c     ende der berechnung der einzelnen matrizen
c***********************************************************************
 4310 do 4315 i=1,np
      xx2r(i)=su1r(i)+rshu1r(i)*sd1r(i)-rshu1i(i)*sd1i(i)
      xx2i(i)=su1i(i)+rshu1i(i)*sd1r(i)+rshu1r(i)*sd1i(i)
      yy1r(i)=xx1r(i)*xx2r(i)-xx1i(i)*xx2i(i)
      yy1i(i)=xx1i(i)*xx2r(i)+xx1r(i)*xx2i(i)
      xx1r(i)=yy1r(i)*rshd2r(i)-yy1i(i)*rshd2i(i)
 4315 xx1i(i)=yy1i(i)*rshd2r(i)+yy1r(i)*rshd2i(i)
      go to 4400
c*****************************************************
 4320 do 4325 i=1,np
      xx2r(i)=sd1r(i)+rshu1r(i)*su1r(i)-rshu1i(i)*su1i(i)
      xx2i(i)=sd1i(i)+rshu1i(i)*su1r(i)+rshu1r(i)*su1i(i)
      xx1r(i)=yy1r(i)*xx2r(i)-yy1i(i)*xx2i(i)
      xx1i(i)=yy1i(i)*xx2r(i)+yy1r(i)*xx2i(i)
      yy1r(i)=xx1r(i)*rshd2r(i)-xx1i(i)*rshd2i(i)
 4325 yy1i(i)=xx1i(i)*rshd2r(i)+xx1r(i)*rshd2i(i)
      go to 4400
c*****************************************************
 4330 do 4335 i=1,np
      xx2r(i)=rshd2r(i)*sd1r(i)-rshd2i(i)*sd1i(i)
      xx2i(i)=rshd2i(i)*sd1r(i)+rshd2r(i)*sd1i(i)
      yy1r(i)=xx1r(i)*xx2r(i)-xx1i(i)*xx2i(i)
 4335 yy1i(i)=xx1i(i)*xx2r(i)+xx1r(i)*xx2i(i)
      if(iss(16).eq.1.and.mdeck.ge.iso) xx1=0.
      go to 4400
c*****************************************************
 4340 xyr=1.
      if(iss(16).eq.1.and.mdeck.ge.ire) xyr=0.
      do 4345 i=1,np
      xx2r(i)=xyr*su1r(i)+rshu1r(i)*sd1r(i)-rshu1i(i)*sd1i(i)
      xx2i(i)=xyr*su1i(i)+rshu1i(i)*sd1r(i)+rshu1r(i)*sd1i(i)
      yy1r(i)=xx1r(i)*xx2r(i)-xx1i(i)*xx2i(i)
 4345 yy1i(i)=xx1i(i)*xx2r(i)+xx1r(i)*xx2i(i)
c*****************************************************
c
c     ende der ber. des versch.potentials
c
c     beginn entfernungsschleife
c
 4400 do 4500 ke=1,nent
      xsum1r(jf,ke)=0.
      xsum1i(jf,ke)=0.
      if(iss(17).eq.2) then
c
c     sh-point source, p-independent
c
         xr=x(ke)
         nh=1
         rr(1)=xr
         dz(1)=dso-z(iso-1)
         ts(1)=0.
         vv(1)=1./sqrt(xr*omeg)
         lamda1=2.
         lamda2=0.
         go to 4520
      end if
c
c     ber. von lamda1,lamda2
c
      phi=az(ke)/57.29578
      cphi=cos(phi)
      sphi=sin(phi)
      c2phi=cos(2.*phi)
      s2phi=sin(2.*phi)
      lamda1=0.5*(m11-m22)*s2phi-m12*c2phi
      lamda2=m13*sphi-m23*cphi
      xr=x(ke)*cphi
      yr=x(ke)*sphi
c     schleife ueber anzahl der punktquellen-ber.d.staerke und d.lage
      do 4510 js=1,nh
      dx=xr-xs(js)
      dy=yr-ys(js)
      rr(js)=sqrt(dx*dx+dy*dy)
      dz(js)=zs(js)-z(iso-1)
 4510 vv(js)=es(js)/sqrt(rr(js))
c     ber. von kappai*su(bzw.sd) ohne zeitl. verschiebung ea bzw. eb
c      *1./sqrt(u) * window(i)
 4520 do 4650 i=1,np
      sp(i)=window(i)*sqrt(p(i))
      bm=csqrt(bmiso-p2(i))
      if(real(bm).eq.0.) bm=-bm
      bmr(i)=real(bm)
      bmi(i)=aimag(bm)
      bmb(i)=1./(bmr(i)*bmr(i)+bmi(i)*bmi(i))
      rshu1r(i)=0.
      rshu1i(i)=0.
      rshd1r(i)=0.
      rshd1i(i)=0.
4650  continue
      do 4610 js=1,nh
      do 4612 i=1,np
      obmr=omegai*bmr(i)
      obmi=omegai*bmi(i)
      uur=(p(i)*rr(js)+ts(js))*omega(jf)
      uui=(p(i)*rr(js)+ts(js))*omegai
      xbr=(-omega(jf)*bmi(i)-obmr)*dz(js)
      xbi=(omega(jf)*bmr(i)-obmi)*dz(js)
      exb1=exp(-xbr+uui)
      exb2=exp(xbr+uui)
      rshu2r(i)=exb1*cos(xbi+uur)*vv(js)
      rshu2i(i)=-exb1*sin(xbi+uur)*vv(js)
      rshd2r(i)=exb2*cos(xbi-uur)*vv(js)
 4612 rshd2i(i)=exb2*sin(xbi-uur)*vv(js)
      do 4614 i=1,np
      rshu1r(i)=rshu2r(i)+rshu1r(i)
      rshu1i(i)=rshu2i(i)+rshu1i(i)
      rshd1r(i)=rshd2r(i)+rshd1r(i)
      rshd1i(i)=rshd2i(i)+rshd1i(i)
 4614 continue
 4610 continue
      do 4616 i=1,np
      rshu2r(i)=-bmr(i)*p(i)*lamda1*bmb(i)+lamda2
      rshu2i(i)=bmi(i)*p(i)*lamda1*bmb(i)
      rshd2r(i)=rshu2r(i)-2.*lamda2
      rshd2i(i)=rshu2i(i)
      su1r(i)=rshu1r(i)*rshu2r(i)-rshu1i(i)*rshu2i(i)
      su1i(i)=rshu1i(i)*rshu2r(i)+rshu1r(i)*rshu2i(i)
      sd1r(i)=rshd1r(i)*rshd2r(i)-rshd1i(i)*rshd2i(i)
      sd1i(i)=rshd1i(i)*rshd2r(i)+rshd1r(i)*rshd2i(i)
      xx2r(i)=yy1r(i)*sd1r(i)-yy1i(i)*sd1i(i)
      xx2i(i)=yy1i(i)*sd1r(i)+yy1r(i)*sd1i(i)
 4616 continue
      do 4618 i=1,np
      if(xx1.eq.0.) go to 4620
      xx2r(i)=xx2r(i)+xx1r(i)*su1r(i)-xx1i(i)*su1i(i)
      xx2i(i)=xx2i(i)+xx1i(i)*su1r(i)+xx1r(i)*su1i(i)
 4620 xsum1r(jf,ke)=xsum1r(jf,ke)+xx2r(i)*sp(i)
      xsum1i(jf,ke)=xsum1i(jf,ke)+xx2i(i)*sp(i)
 4618 continue
      xsum1r(jf,ke)=xsum1r(jf,ke)-(xx2r(1)*sp1+xx2r(np)*sp(np))*0.5
      xsum1i(jf,ke)=xsum1i(jf,ke)-(xx2i(1)*sp1+xx2i(np)*sp(np))*0.5
      xsum1r(jf,ke)=xsum1r(jf,ke)*d11x*dp
      xsum1i(jf,ke)=xsum1i(jf,ke)*d11x*dp
 4500 continue
 4000 continue
      do 4550 ke=1,nent
      t0=x(ke)/vred+tmin
      do 4505 j=1,nfreq
 4505 phase(j)=omega(j)*t0+pi4
      do 4507 j=1,npts
 4507 s(j)=(0.,0.)
      jf1=1
      do 4600 jf=1,inpt
      if(jf.le.iu.or.jf.ge.io) go to 4600
      jf1=jf1+1
      s(jf)=cmplx(xsum1r(jf1,ke),xsum1i(jf1,ke))*cexp(ein*phase(jf1))*
     &bmiso2*wavel(jf)
      s(npts-jf+2)=conjg(s(jf))
 4600 continue
      if(iss(6).eq.1) then
         smax=abs(real(s(1)))
         if(abs(aimag(s(1))).gt.smax) smax=abs(aimag(s(1)))
         do 4701 jf=2,inpt
         if(abs(real(s(jf))).gt.smax) smax=abs(real(s(jf)))
 4701    if(abs(aimag(s(jf))).gt.smax) smax=abs(aimag(s(jf)))
         do 4702 jf=1,inpt
 4702    sc1(jf)=s(jf)/smax
         write(23,5004) fu,fwil,fwir,fo,dnue,c2,cwir,cwil,c1,np
         write(23,5005) x(ke),npts,smax
         write(23,5006) (sc1(jf),jf=1,inpt)
 5004 format(/' frequenzband:'/4f10.4/' frequenzintervall:'/f10.6/
     &' strahlparameterfenster:'/4f10.4/' anzahl der strahlparameter:'/
     &i10)
 5005 format(/f10.4,i10,e15.5)
 5006 format(10f8.5)
      end if
      inp1=inpt+1
      s(inp1)=(0.0,0.0)
      call zzlog(npts,mga)
      call cool(mga,s,1.)
      ddd=1.
      ccc=1.
      if(tsigma.ne.0.) then
         ddd=exp(t0/tsigma)
         ccc=exp(dt/tsigma)
      end if
      dd=1./ccc
      do 4700 i=1,npts
      dd=dd*ccc
 4700 s(i)=s(i)*dd*ddd
c
c     write out seismograms on units 6 and ioutsh
c
      sphred=1.
      if(iss(12).eq.1)  sphred=sqrt((x(ke)/3390.)/abs(sin(x(ke)/3390.)))
     1*((3390.-dso)/3390.)**(float(nrho+1)/2.)
      do 4710 j=1,npts
 4710 balla(j)=sphred*real(s(j))
      balmax = amax(balla,npts)
      bal = abs(amin(balla,npts))
      if(bal.gt.balmax) balmax = bal
c
c create sac file
c
      call newhdr
      call setfhv('B', t0, nerr)
      call setfhv('DELTA', dt, nerr)
      call setfhv('EVDP', zs(1),nerr)
      call setnhv('NPTS', npts, nerr)
      call setlhv('LEVEN', .true., nerr)
      call setfhv('DIST', x(ke), nerr)
      call setfhv('GCARC', x(ke)/59.16, nerr)
      call setfhv('AZ', az(ke), nerr)
      call base(filename,ke)
      filename(6:7)=".t"
      call wsac0(filename,balla(1),balla(1),nerr)
c
      do 4720 i = 1,npts
 4720 iballa(i) = ifix(999.10*balla(i)/balmax)
      if(iss(15).eq.1) go to 4750
      epi=x(ke)*0.0169
      write(6,4730) x(ke),epi
 4730 format(///,1x,'sh-seismogram  for the distance',f10.3,' km  =',
     1f8.3,' degree'/)
      write(6,4013) npts,balmax
      write(6,4004) (iballa(i), i = 1,npts)
 4004 format(20i4)
 4013 format(i10,e15.4)
 4750 do 4740 i=1,npts
      balla(i) = 4.99*balla(i)/balmax + 5.001
 4740 iballa(i) = 10000.*balla(i)
      abal = tmin - na*dt + x(ke)/vred
      ikom=3
cc      write(ioutsh,5001)x(ke),abal,ikom,npts,balmax
 5001 format(2f15.5,i5/,i10,5x,e15.4) 
      write(6,4005) x(ke),abal,kom,npts,balmax
 4005 format(2f10.4,i5,i10,5x,e15.4)
cc      write(ioutsh,5002)(iballa(i),i=1,npts)
cc      write(ioutsh,5003)x(ke)
 5002 format(16i5)
 5003 format(f15.5)
 4550 continue
      return
      end
c==================================================================
      subroutine base(filename,num)
c==================================================================
      integer num
      character*7 filename
 
      filename(1:5) = 'st000'
      IF (NUM.le.9.) THEN
         filename(5:5) = CHAR(NUM+48)
      ELSE IF (NUM.le.99) THEN
         filename(4:4) = CHAR( NUM/10 +48 )
         filename(5:5) = CHAR( MOD(NUM,10) +48 )
      ELSE
         filename(3:3) = CHAR( NUM/100 +48 )
         filename(4:4) = CHAR( MOD(NUM,100)/10 +48 )
         filename(5:5) = CHAR( MOD( MOD(NUM,100),10 ) +48 )
      ENDIF
      RETURN
      END

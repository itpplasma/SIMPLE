ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 08.03.2014      subroutine gbstfp(scon,tcon,phi,bscon,btcon,bfcon
      subroutine gbstfp(xcon,zcon,phic,bxcon,bzcon,bfcon
c 08.03.2014     *,gsr,gsf,gsz,gtr,gtf,gtz,gfr,gff,gfz
     *,gxcr,gxcf,gxcz,gzcr,gzcf,gzcz,gfr,gff,gfz
     *,r,z,br,bf,bz
     *,brr,brf,brz,bfr,bff,bfz,bzr,bzf,bzz)
      implicit real*8 (a-h,o-z),integer(i-n)

c scon, tcon, phi are the input parameters, these are the contravariant
c coordinates of the background coordinate system; scon is the modified
c radial coordinate, tcon (theta) is the poloidal angle, phi is the
c toroidal angle in the cylindrical coordinates;

c the rest parameters are the output parameters, all they are calculated
c trough scon,tcon and phi;
c bscon, btcon and bfcon   are the contravariant components of the magnetic
c field in the coordinates scon,tcon and phi;
c r and z are the cylindrical coordinates;
c br, bf and bz are the magnetic field components in the cylindrical
c coordinates;
c brr, ..., bzz are the derivatives of the magnetic field components
c in the cylindrical coordinates with respect to the cylindrical coordinates;

c the scon value changes from 0 to s_max where s_max coincides with the value
c of the maximum radial index;

c the maximum radial index, the maxima poloidal and toroidal harmonic numbers
c as  well as the number of field periods are presented in the file param.f
c as parameters, these quantities are also at the beginning of the file
c pies_extrac, which contains the data;

c the program reads the necessary data from the file pies_extrac with the
c operator read(71,*).

c 08.03.2014
c 28.09.2016 
      common/csoi/nsoi0
c 28.09.2016 end
c 27.03.2016
      common/cismax/ismax
      logical prop
      data prop/.false./
      ismax=0
c 27.03.2016 end
      scon=dsqrt(xcon**2+zcon**2)
      tcon=datan2(zcon,xcon)
      cost=dcos(tcon)
      sint=dsin(tcon)
      dxcds=cost
      dzcds=sint
      dxcdt=-scon*sint
      dzcdt=scon*cost
c 08.03.2014 end

c 25.03.2016
c 09.07.2016      facr=0.01d0
      facr=1d0
c 25.03.2016 end

      soi=scon
      fii=phic
      teti=tcon
c 28.09.2016
      nsoi=soi
      nsoi0=nsoi+1
c 28.09.2016 end

      call gbpi(soi,fii,teti
     *,roe,ze,brc,bfc,bzc
c 01.06.2016    *,dbrods,dbrodt,dbrodf,dbteds,dbtedt,dbtedf
c 01.06.2016     *,dbfids,dbfidt,dbfidf,bfco
c 22.08.2016     *,bfco
c 06.07.2016     *,e1r,e1z,e1f,e2r,e2z,e2f,e3r,e3z,e3f
     *,grsr,grsz,grsf,grtr,grtz,grtf,grfr,grfz,grff
     *,brs,brt,brfi,bfs,bft,bffi,bzs,bzt,bzfi)

c 27.03.2016
c 11.07.2016      if(ismax.eq.1) return
c      if(ismax.eq.1) then
c      print*,'ismax=1, gbstfp, stop'
c      stop
c      end if
c 27.03.2016 end
c 02.10.2016
      if(ismax.eq.1) return
c 02.10.2016 end
      grxr=dxcds*grsr+dxcdt*grtr
      grxf=dxcds*grsf+dxcdt*grtf
      grxz=dxcds*grsz+dxcdt*grtz

      grzr=dzcds*grsr+dzcdt*grtr
      grzf=dzcds*grsf+dzcdt*grtf
      grzz=dzcds*grsz+dzcdt*grtz

c 25.03.2016 (facr)
      gxcr=grxr*facr
      gxcf=grxf*facr
      gxcz=grxz*facr

      gzcr=grzr*facr
      gzcf=grzf*facr
      gzcz=grzz*facr
c 25.03.2016 (facr) end

c 08.03.2014 end

c 25.03.2016 (facr) 
      gfr=grfr*facr
      gff=grff*facr
      gfz=grfz*facr

      r=roe/facr
      z=ze/facr
c 25.03.2016 (facr) end

      br=brc
      bf=bfc
      bz=bzc
c 10.07.2016
      bs=br*grsr+bf*grsf+bz*grsz      
      bt=br*grtr+bf*grtf+bz*grtz
c 10.07.2016 end

c 08.03.2014
      bxcon=br*grxr+bf*grxf+bz*grxz
      bzcon=br*grzr+bf*grzf+bz*grzz
      bfcon=br*grfr+bf*grff+bz*grfz
c 08.03.2014 end

c 25.03.2016 (facr)
      brr=(brs*grsr+brt*grtr+brfi*grfr)*facr
      brf=(brs*grsf+brt*grtf+brfi*grff)*roe
      brz=(brs*grsz+brt*grtz+brfi*grfz)*facr

      bfr=(bfs*grsr+bft*grtr+bffi*grfr)*facr
      bff=(bfs*grsf+bft*grtf+bffi*grff)*roe
      bfz=(bfs*grsz+bft*grtz+bffi*grfz)*facr

      bzr=(bzs*grsr+bzt*grtr+bzfi*grfr)*facr
      bzf=(bzs*grsf+bzt*grtf+bzfi*grff)*roe
      bzz=(bzs*grsz+bzt*grtz+bzfi*grfz)*facr
c 25.03.2016 (facr) end
      if(prop) return
      prop=.true.
   61 format(3e14.5,2x,a16)
      print 61,xcon,zcon,phic,'xcon,zcon,phic'
      print 61,soi,fii,teti,'soi,fii,teti'
      print 61,bxcon,bzcon,bfcon,'bxcon,bzcon,bfcon'
      print 61,bs,bt,bfcon,'bs,bt,bfcon'
      print 61,gxcr,gxcf,gxcz,'gxcr,gxcf,gxcz'
      print 61,gzcr,gzcf,gzcz,'gzcr,gzcf,gzcz'
      print 61,gfr,gff,gfz,'gfr,gff,gfz'
      print*,r,z,'r,z'
      print 61,br,bf,bz,'br,bf,bz'
      print 61,brr,brf,brz,'brr,brf,brz'
      print 61,bfr,bff,bfz,'bfr,bff,bfz'
      print 61,bzr,bzf,bzz,'bzr,bzf,bzz'
      print*,'1st end of gbstfp'
c      write(99,61)xcon,zcon,phic,'xcon,zcon,phic'
c      write(99,61)soi,fii,teti,'soi,fii,teti'
c      write(99,*)nsoi0,'nsoi0'
c      write(99,61)bxcon,bzcon,bfcon,'bxcon,bzcon,bfcon'
c      write(99,61)bs,bt,bfcon,'bs,bt,bfcon'
c      write(99,61)gxcr,gxcf,gxcz,'gxcr,gxcf,gxcz'
c      write(99,61)gzcr,gzcf,gzcz,'gzcr,gzcf,gzcz'
c      write(99,61)gfr,gff,gfz,'gfr,gff,gfz'
c      write(99,*)r,z,'r,z'
c      write(99,61)br,bf,bz,'br,bf,bz'
c      write(99,61)brr,brf,brz,'brr,brf,brz'
c      write(99,61)bfr,bff,bfz,'bfr,bff,bfz'
c      write(99,61)bzr,bzf,bzz,'bzr,bzf,bzz'
c      write(99,61)br/bf,bf/bf,bz/bf,'bro,bfo,bzo'
c      write(99,61)brr/bf,brf/bf,brz/bf,'brro,brfo,brzo'
c      write(99,61)bfr/bf,bff/bf,bfz/bf,'bfro,bffo,bfzo'
c      write(99,61)bzr/bf,bzf/bf,bzz/bf,'bzro,bzfo,bzzo'
c      write(99,*)'1st end of gbstfp'
      close(99)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 01.06.2016      subroutine gbpi(soi,fii,teti,brodbf,btdbf,ros
      subroutine gbpi(soi,fii,teti
     *,roe,ze,br,bf,bz
     *,grsr,grsz,grsf,grtr,grtz,grtf,grfr,grfz,grff
     *,brs,brt,brf,bfs,bft,bff,bzs,bzt,bzf)
!
      use vmec_stuff_mod                                                !<=2017 NEW
!
      implicit real*8 (a-h,o-z),integer(i-n)
c 28.69.2016

      parameter(pi=3.14159265358979d0)
!<=2017 NEW     include 'parvmec.f'
c      parameter(nsurfm=64,nstrm=61,nper=5,kpar=63)
!<=2017 NEW     save prop,kpai,kstri,kspl,twopdn,ss,
      save prop,kpai,kstri,kspl,twopdn,
!<=2017 NEW     *smin,smax,srmn,szmn,slmn,axm,axn,saiota,ssps,sphi,psi1
     *smin,smax,psi1
      common/cflux/flux1
c 27.03.2016
      save/cismax/
c 27.03.2016 end
!<=2017 NEW      dimension soa(0:kpar)
c 27.03.2016
!<2017 NEW      dimension srmn(4,nstrm,kpar)
!<2017 NEW      dimension szmn(4,nstrm,kpar)
!<2017 NEW      dimension slmn(4,nstrm,kpar)
!<2017 NEW      dimension saiota(4,kpar)
!<2017 NEW      dimension ssps(4,kpar)      
!<2017 NEW      dimension sphi(4,kpar)
!<2017 NEW      dimension ss(4,kpar)
!<2017 NEW      dimension axm(nstrm),axn(nstrm)
c 01.06.2016 end
      common/cismax/ismax
c 27.03.2016 end
c 12.07.2016
      common/csps/spsc
c 12.07.2016 end
c 27.09.2016
      common/cisc/iscc
      save/cisc/
c      save iscc
c 27.09.2016 end
      logical prop
      data prop/.false./
      if(prop) go to 100
c 10.07.2016      prop=.true.
      kpai=kpar
      kstri=nstrm
      kspl=4
c 28.06.2016      
      twopdn=2d0*pi/nper
c      psi1=flux1/(2d0*pi)
c 28.06.2016 end

      do is=0,kpai
      soa(is)=is
      end do
      smin=soa(0)
      smax=soa(kpai)

c 22.08.2016      call vvobooz(ssps,saiota,scurt,scurp,sbmn,srmn,szmn,sgmn
      call vvovmec(ssps,sphi,saiota,srmn,szmn,slmn,axm,axn,ss
     *,kspl,kstri,kpai)

      psi1=flux1/(2d0*pi)
  100 continue
!
c 27.03.2016      ierrfield=0
      ismax=0
!
      ier=0
      so=soi
c 19.09.2016      fin=fii*nper
      fin=fii
      tet=teti
      if(so.lt.smin) then
      ier=1
      print*,'  so < smin'
      print*,so,'=so',smin,'=smin'
      end if
      if(so.ge.smax) then
      ier=1
      print*,'  so > smax'
      print*,so,'=so',smax,'=smax'
      end if

      if(ier.eq.1) then
      ismax=1
      return
      end if

      so0=smin
      do is=1,kpai
      if(is.gt.1) so0=soa(is-1)
      if (so.lt.soa(is)) then
      ds=so-so0
      isc=is
      go to 110
      end if
      end do
  110 continue

      rmnb=0d0
      rmnbs=0d0
      rmnbt=0d0
      rmnbf=0d0
      rmnbss=0d0

      zmnb=0d0
      zmnbs=0d0
      zmnbt=0d0
      zmnbf=0d0
      zmnbss=0d0

      almnb=0d0
      almnbs=0d0
      almnbt=0d0
      almnbf=0d0
      almnbss=0d0
      almnbts=0d0
      almnbfs=0d0
      almnbtt=0d0
      almnbtf=0d0
      almnbff=0d0

      aiotb=0d0
      aiotbs=0d0
      aiotbss=0d0

      ssb=0d0      

      spsb=0d0
      spsbs=0d0
      spsbss=0d0

      phib=0d0
      phibs=0d0
      phibss=0d0

      e1rsb=0d0
      e1rtb=0d0
      e1rfb=0d0
      e2rtb=0d0
      e2rfb=0d0
      e3rfb=0d0
      e1zsb=0d0
      e1ztb=0d0
      e1zfb=0d0
      e2ztb=0d0
      e2zfb=0d0
      e3zfb=0d0

      rob=0d0
      rosb=0d0
      drdsb=0d0
      drdtb=0d0
      drdfb=0d0
      zb=0d0
      dzdsb=0d0
      dzdtb=0d0
      dzdfb=0d0

      ssb=((ss(4,isc)*ds+ss(3,isc))*ds+
     *ss(2,isc))*ds+ss(1,isc)

      aiotb=((saiota(4,isc)*ds+saiota(3,isc))*ds+
     *saiota(2,isc))*ds+saiota(1,isc)
      aiotbs=(3d0*saiota(4,isc)*ds+2d0*saiota(3,isc))*ds+
     *saiota(2,isc)
      aiotbss=6d0*saiota(4,isc)*ds+2d0*saiota(3,isc)

      spsb=((ssps(4,isc)*ds+ssps(3,isc))*ds+
     *ssps(2,isc))*ds+ssps(1,isc)
      spsbs=(3d0*ssps(4,isc)*ds+2d0*ssps(3,isc))*ds+
     *ssps(2,isc)
      spsbss=6d0*ssps(4,isc)*ds+2d0*ssps(3,isc)
      spsc=spsb

      phib=((sphi(4,isc)*ds+sphi(3,isc))*ds+
     *sphi(2,isc))*ds+sphi(1,isc)
      phibs=(3d0*sphi(4,isc)*ds+2d0*sphi(3,isc))*ds+
     *sphi(2,isc)
      phibss=6d0*sphi(4,isc)*ds+2d0*sphi(3,isc)

c      do km=0,mpai
      do kst=1,kstri
c      do kn=-npai,npai

         alpha=0.5d0*axm(kst)                               !<=NEW
         salp=so**alpha                                     !<=NEW
         dsalp=alpha*salp/so                                !<=NEW
         ddsalp=(alpha-1.d0)*dsalp/so                       !<=NEW

      rbh=((srmn(4,kst,isc)*ds+srmn(3,kst,isc))*ds+
     *srmn(2,kst,isc))*ds+srmn(1,kst,isc)
      rbhs=(3d0*srmn(4,kst,isc)*ds+2d0*srmn(3,kst,isc))*ds+
     *srmn(2,kst,isc)
      rbhss=6d0*srmn(4,kst,isc)*ds+2d0*srmn(3,kst,isc)

         rbhss=ddsalp*rbh+2.d0*dsalp*rbhs+salp*rbhss        !<=NEW
         rbhs=dsalp*rbh+salp*rbhs                           !<=NEW
         rbh=salp*rbh                                       !<=NEW

      zbh=((szmn(4,kst,isc)*ds+szmn(3,kst,isc))*ds+
     *szmn(2,kst,isc))*ds+szmn(1,kst,isc)
      zbhs=(3d0*szmn(4,kst,isc)*ds+2d0*szmn(3,kst,isc))*ds+
     *szmn(2,kst,isc)
      zbhss=6d0*szmn(4,kst,isc)*ds+2d0*szmn(3,kst,isc)

         zbhss=ddsalp*zbh+2.d0*dsalp*zbhs+salp*zbhss        !<=NEW
         zbhs=dsalp*zbh+salp*zbhs                           !<=NEW
         zbh=salp*zbh                                       !<=NEW


      albh=((slmn(4,kst,isc)*ds+slmn(3,kst,isc))*ds+
     *slmn(2,kst,isc))*ds+slmn(1,kst,isc)
      albhs=(3d0*slmn(4,kst,isc)*ds+2d0*slmn(3,kst,isc))*ds+
     *slmn(2,kst,isc)
      albhss=6d0*slmn(4,kst,isc)*ds+2d0*slmn(3,kst,isc)

         albhss=ddsalp*albh+2.d0*dsalp*albhs+salp*albhss    !<=NEW
         albhs=dsalp*albh+salp*albhs                        !<=NEW
         albh=salp*albh                                     !<=NEW

c      knnper=kn*nper
c 01.06.2016      cosf=dcos(fin*kn-tet*km)
c 01.06.2016      sinf=dsin(fin*kn-tet*km)
c 17.09.2016      cosf=dcos(fin*kn+tet*km)
c 17.09.2016      sinf=dsin(fin*kn+tet*km)
      xm=axm(kst)
      xn=axn(kst)
      cosf=dcos(xm*tet-xn*fin)
      sinf=dsin(xm*tet-xn*fin)
c 01.06.2016
      rmnb=rmnb+rbh*cosf
      rmnbs=rmnbs+rbhs*cosf
      rmnbt=rmnbt-xm*rbh*sinf
      rmnbf=rmnbf+xn*rbh*sinf
      rmnbss=rmnbss+rbhss*cosf

      zmnb=zmnb+zbh*sinf
      zmnbs=zmnbs+zbhs*sinf
      zmnbt=zmnbt+xm*zbh*cosf
      zmnbf=zmnbf-xn*zbh*cosf
      zmnbss=zmnbss+zbhss*sinf

      almnb=almnb+albh*sinf
      almnbs=almnbs+albhs*sinf
      almnbt=almnbt+xm*albh*cosf
c      gmnbf=gmnbf+xn*gbh*cosf
      almnbf=almnbf-xn*albh*cosf
      almnbss=almnbss+albhss*sinf
      almnbts=almnbts+xm*albhs*cosf
c      gmnbfs=gmnbfs+xn*gbhs*cosf
      almnbfs=almnbfs-xn*albhs*cosf
      almnbtt=almnbtt-xm*xm*albh*sinf
c      gmnbtf=gmnbtf-xm*xn*gbh*sinf
      almnbtf=almnbtf+xm*xn*albh*sinf
c      gmnbff=gmnbff-xn*xn*gbh*sinf
      almnbff=almnbff-xn*xn*albh*sinf

      rob=rmnb
c?      rosb=rosb+x1bhs

      drdsb=rmnbs
      drdtb=rmnbt
      drdfb=rmnbf

      zb=zmnb
      dzdsb=zmnbs
      dzdtb=zmnbt
      dzdfb=zmnbf

      e1rsb=rmnbss
      e1rtb=e1rtb-xm*rbhs*sinf
      e1rfb=e1rfb+xn*rbhs*sinf

      e2rtb=e2rtb-xm*xm*rbh*cosf
      e2rfb=e2rfb+xn*xm*rbh*cosf
      e3rfb=e3rfb-xn*xn*rbh*cosf

      e1zsb=zmnbss
      e1ztb=e1ztb+xm*zbhs*cosf
      e1zfb=e1zfb-xn*zbhs*cosf

      e2ztb=e2ztb-xm*xm*zbh*sinf
      e2zfb=e2zfb+xn*xm*zbh*sinf
      e3zfb=e3zfb-xn*xn*zbh*sinf

      end do

      almnbft=almnbtf

      e1r=drdsb
      e1z=dzdsb
c 20.09.2016      e1f=fcylsb*rob
      e1f=0d0

      e2r=drdtb
      e2z=dzdtb
c 20.09.2016      e2f=fcyltb*rob
      e2f=0d0

      e3r=drdfb
      e3z=dzdfb
c 20.09.2016      e3f=fcylfb*rob
      e3f=rob

      e1rs=e1rsb
      e1rt=e1rtb
      e1rf=e1rfb

      e2rs=e1rtb
      e2rt=e2rtb
      e2rf=e2rfb

      e3rs=e1rfb
      e3rt=e2rfb
      e3rf=e3rfb

      e1fs=0d0
      e1ft=0d0
      e1ff=0d0

      e2fs=0d0
      e2ft=0d0
      e2ff=0d0
     
      e3fs=e1r
      e3ft=e2r
      e3ff=e3r

      e1zs=e1zsb
      e1zt=e1ztb
      e1zf=e1zfb

      e2zs=e1ztb
      e2zt=e2ztb
      e2zf=e2zfb

      e3zs=e1zfb
      e3zt=e2zfb
      e3zf=e3zfb

      roe=rob
c      ros=rosb
      ros=1d0/smax
      ze=zb
c 06.07.2016
      ve1e2r=0d0
      ve1e2f=e1z*e2r-e1r*e2z
      ve1e2z=0d0
c      sqg=ve1e2r*e3r+ve1e2f*e3f+ve1e2z*e3z
      sqg=ve1e2f*e3f
      osqg=1d0/sqg

      aiota=aiotb
      aiotas=aiotbs
      aiotat=0d0
      aiotap=0d0
      aiotbt=0d0
      aiotbf=0d0
      phibst=0d0
      phibsf=0d0
c 07.07.2016
c   61 format(3e14.5,2x,a16)
c      print 61,e3r,e3f,e3z,'e3r,e3f,e3z'
c      print 61,e2r,e2f,e2z,'e2r,e2f,e2z'
c      print 61,aiota,osqg,psi1,'aiota,osqg,psi1'
c      print*,spshs,'spshs',isc,'=isc'

      btet=phibs*(aiotb-almnbf)*osqg
      bfii=phibs*(1d0+almnbt)*osqg

      br=btet*e2r+bfii*e3r
      bf=bfii*e3f
      bz=btet*e2z+bfii*e3z



c 23.09.2016      bro=(e3r+aiota*e2r)*osqg
c 23.09.2016      bfo=(e3f+aiota*e2f)*osqg
c 23.09.2016      bzo=(e3z+aiota*e2z)*osqg

c 23.09.2016      br=bro*spshs*psi1
c 23.09.2016      bf=bfo*spshs*psi1
c 23.09.2016      bz=bzo*spshs*psi1
c 07.07.2016 end
      ve1e2fs=e1zs*e2r-e1rs*e2z+e1z*e2rs-e1r*e2zs
      ve1e2ft=e1zt*e2r-e1rt*e2z+e1z*e2rt-e1r*e2zt
      ve1e2ff=e1zf*e2r-e1rf*e2z+e1z*e2rf-e1r*e2zf

      sqgs=ve1e2fs*e3f+ve1e2f*e3fs
      sqgt=ve1e2ft*e3f+ve1e2f*e3ft
      sqgf=ve1e2ff*e3f+ve1e2f*e3ff

c 16.10.2016
      phibss=0d0
c 16.10.2016 end

      btets=(phibss*(aiotb-almnbf)+phibs*(aiotbs-almnbfs)
     *-btet*sqgs)*osqg
      btett=(phibst*(aiotb-almnbf)+phibs*(aiotbt-almnbft)
     *-btet*sqgt)*osqg
      btetf=(phibsf*(aiotb-almnbf)+phibs*(aiotbf-almnbff)
     *-btet*sqgf)*osqg

      bfiis=(phibss*(1d0+almnbt)+phibs*almnbts-bfii*sqgs)*osqg
      bfiit=(phibst*(1d0+almnbt)+phibs*almnbtt-bfii*sqgt)*osqg
      bfiif=(phibsf*(1d0+almnbt)+phibs*almnbtf-bfii*sqgf)*osqg

      brs=btets*e2r+bfiis*e3r+btet*e2rs+bfii*e3rs
      brt=btett*e2r+bfiit*e3r+btet*e2rt+bfii*e3rt
      brf=btetf*e2r+bfiif*e3r+btet*e2rf+bfii*e3rf

      bfs=bfiis*e3f+bfii*e3fs
      bft=bfiit*e3f+bfii*e3ft
      bff=bfiif*e3f+bfii*e3ff

      bzs=btets*e2z+bfiis*e3z+btet*e2zs+bfii*e3zs
      bzt=btett*e2z+bfiit*e3z+btet*e2zt+bfii*e3zt
      bzf=btetf*e2z+bfiif*e3z+btet*e2zf+bfii*e3zf

c 23.09.2016     brs=(e3rs+aiota*e2rs+aiotas*e2r-bro*sqgs)*osqg*spshs+spshss*bro
c 23.09.2016      brs=brs*psi1
c 23.09.2016      brt=(e3rt+aiota*e2rt+aiotat*e2r-bro*sqgt)*osqg*spshs*psi1
c 23.09.2016      brf=(e3rf+aiota*e2rf+aiotaf*e2r-bro*sqgp)*osqg*spshs*psi1

c 23.09.2016     bfs=(e3fs+aiota*e2fs+aiotas*e2f-bfo*sqgs)*osqg*spshs+spshss*bfo
c 23.09.2016      bfs=bfs*psi1
c 23.09.2016      bft=(e3ft+aiota*e2ft+aiotat*e2f-bfo*sqgt)*osqg*spshs*psi1
c 23.09.2016      bff=(e3ff+aiota*e2ff+aiotaf*e2f-bfo*sqgp)*osqg*spshs*psi1

c 23.09.2016     bzs=(e3zs+aiota*e2zs+aiotas*e2z-bzo*sqgs)*osqg*spshs+spshss*bzo
c 23.09.2016      bzs=bzs*psi1
c 23.09.2016      bzt=(e3zt+aiota*e2zt+aiotat*e2z-bzo*sqgt)*osqg*spshs*psi1
c 23.09.2016      bzp=(e3zp+aiota*e2zp+aiotap*e2z-bzo*sqgp)*osqg*spshs*psi1

      grfr=0d0
      grff=ve1e2f*osqg
      grfz=0d0

      ve2e3r=-e2z*e3f
      ve2e3f=e2z*e3r-e2r*e3z
      ve2e3z=e2r*e3f

      grsr=ve2e3r*osqg
      grsf=ve2e3f*osqg
      grsz=ve2e3z*osqg

      ve3e1r=e3f*e1z-e3z*e1f
      ve3e1f=e3z*e1r-e3r*e1z
      ve3e1z=e3r*e1f-e3f*e1r

      grtr=ve3e1r*osqg
      grtf=ve3e1f*osqg
      grtz=ve3e1z*osqg
c 06.07.2016 end
c 10.07.2016
      bs=br*grsr+bf*grsf+bz*grsz
      bt=br*grtr+bf*grtf+bz*grtz
      bfco=br*grfr+bf*grff+bz*grfz
c      print 61,bs,bt,bp,'bs,bt,bp'
c 10.07.2016 end
ccc      br=brb
c      bf=(rt+roe)*bfb
ccc      bf=roe*bup3b
ccc      bz=bzb

c 06.07.2016      bro=bup1b*bup3b
c 06.07.2016      bte=bup2b*bup3b
c      bfco=bup3b

c 06.07.2016      br=bro*e1r*smax+bte*e2r+bup3b*e3r
c 06.07.2016      bf=bro*e1f*smax+bte*e2f+bup3b*e3f
c 06.07.2016      bz=bro*e1z*smax+bte*e2z+bup3b*e3z

c 06.07.2016      brs=(dbrods*e1r+bro*e1rs)*smax+dbteds*e2r+bte*e2rs
c 06.07.2016     *+dbfids*e3r+bup3b*e3rs
c 06.07.2016      brt=(dbrodt*e1r+bro*e1rt)*smax+dbtedt*e2r+bte*e2rt
c 06.07.2016     *+dbfidt*e3r+bup3b*e3rt
c 06.07.2016      brfi=(dbrodf*e1r+bro*e1rf)*smax+dbtedf*e2r+bte*e2rf
c 06.07.2016     *+dbfidf*e3r+bup3b*e3rf

c 06.07.2016      bfs=(dbrods*e1f+bro*e1fs)*smax+dbteds*e2f+bte*e2fs
c 06.07.2016     *+dbfids*e3f+bup3b*e3fs
c 06.07.2016      bft=(dbrodt*e1f+bro*e1ft)*smax+dbtedt*e2f+bte*e2ft
c 06.07.2016     *+dbfidt*e3f+bup3b*e3ft
c 06.07.2016      bffi=(dbrodf*e1f+bro*e1ff)*smax+dbtedf*e2f+bte*e2ff
c 06.07.2016     *+dbfidf*e3f+bup3b*e3ff

c 06.07.2016      bzs=(dbrods*e1z+bro*e1zs)*smax+dbteds*e2z+bte*e2zs
c 06.07.2016     *+dbfids*e3z+bup3b*e3zs
c 06.07.2016      bzt=(dbrodt*e1z+bro*e1zt)*smax+dbtedt*e2z+bte*e2zt
c 06.07.2016     *+dbfidt*e3z+bup3b*e3zt
c 06.07.2016      bzfi=(dbrodf*e1z+bro*e1zf)*smax+dbtedf*e2z+bte*e2zf
c 06.07.2016     *+dbfidf*e3z+bup3b*e3zf

c 02.06.2016      broso=dbrods
c 02.06.2016      brote=dbrodt
c 02.06.2016      brofi=dbrodf

c 02.06.2016      bteso=dbteds
c 02.06.2016      btete=dbtedt
c 02.06.2016      btefi=dbtedf

c 02.06.2016      bfiso=dbfids
c 02.06.2016      bfite=dbfidt
c 02.06.2016      bfifi=dbfidf
      if(prop) return
      prop=.true.
   61 format(3e14.5,2x,a16)
      print*,fin,'=fii, gbpi'
      print 61,ssb,roe,ze,'=s, r, z'      
c      print 61,fcylsb,fcyltb,fcylfb,'fcylsb,fcyltb,fcylfb'
      print 61,e1r,e1f,e1z,'e3r,e3f,e3z'
      print 61,e2r,e2f,e2z,'e2r,e2f,e2z'
      print 61,e3r,e3f,e3z,'e3r,e3f,e3z'
      print 61,grfr,grff,grfz,'grpr,grpf,grpz'
      print 61,aiotb,osqg,psi1,'aiota,osqg,psi1'
      print 61,almnb,almnbt,almnbf,'lmnb,lmnbt,lmnbf'
      print 61,btet,bfii,phibs,'btet,bfii,phibs'
      print 61,bs,bt,bfco,'bs,bt,bfco'
      print 61,br,bf,bz,'br,bf,bz'
      print 61,almnb,almnbt,almnbf,'lmnb,lmnbt,lmnbf'    
      print*,twopdn,'=twopdn, gbpi 1st end'
c      write(99,*)fin,'=fii, gbpi'
c      write(99,61)ssb,roe,ze,'=s, r, z'
cc      write(99,61)fcylsb,fcyltb,fcylfb,'fcylsb,fcyltb,fcylfb'
c      write(99,61)e1r,e1f,e1z,'e3r,e3f,e3z'
c      write(99,61)e2r,e2f,e2z,'e2r,e2f,e2z'
c      write(99,61)e3r,e3f,e3z,'e3r,e3f,e3z'
c      write(99,61)grfr,grff,grfz,'grfr,grff,grfz'
c      write(99,61)aiotb,osqg,psi1,'aiota,osqg,psi1'
c      write(99,61)almnb,almnbt,almnbf,'lmnb,lmnbt,lmnbf'
c      write(99,61)btet,bfii,phibs,'btet,bfii,phibs'
c      write(99,61)bs,bt,bfco,'bs,bt,bfco'
c      write(99,61)br,bf,bz,'br,bf,bz'
      iscc=isc
c      write(99,*)spsbs,'spsbs',isc,'=isc'
c      write(99,*)phibs,'phibs',phibss,'=phibss'
c      write(99,*)twopdn,'=twopdn, gbpi 1st end'

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 02.06.2016ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine vvopie(sbup1,sbup2,sbup3,sbxby1,sbxby2,sx1,sx2
c 22.08.2016      subroutine vvobooz(ssps,saiota,scurt,scurp,sbmn,srmn,szmn,sgmn
      subroutine vvovmec(ssps,sphi,saiota,srmn,szmn,slmn,axm,axn,ss
     *,kspl,nstri,kpai)
!
      use vmec_stuff_mod, only : nsurfm,nstrm,nper,kpar                 !<=2017 NEW
!
      implicit real*8 (a-h,o-z),integer(i-n)
c      include 'parncsx.f'
!<=2017 NEW      include 'parvmec.f'
c      parameter(nsurfm=64,nstrm=61,nper=5,kpar=63)

      dimension srmn(kspl,nstri,kpai)
      dimension szmn(kspl,nstri,kpai)
      dimension slmn(kspl,nstri,kpai)
      dimension axm(nstri),axn(nstri)
      dimension saiota(kspl,kpai),sphi(kspl,kpai),ssps(kspl,kpai)
      dimension ss(kspl,kpai)
!<=2017 NEW      dimension soa(0:kpar)
!<=2017 NEW      dimension xs(kpar),ys(kpar),as(kpar),bs(kpar),cs(kpar),ds(kpar)

!<=2017 NEW      dimension rmn(nstrm,0:kpar)
!<=2017 NEW      dimension zmn(nstrm,0:kpar)
!<=2017 NEW      dimension almn(nstrm,0:kpar)
!<=2017 NEW      dimension arm(nstrm),arn(nstrm)
!<=2017 NEW      dimension aiota(0:kpar),phi(0:kpar),sps(0:kpar),s(0:kpar)
      double precision, dimension(:),   allocatable :: soa,xs,ys,       !<=2017 NEW
     *                        as,bs,cs,ds,arm,arn,aiota,phi,sps,s       !<=2017 NEW
      double precision, dimension(:,:), allocatable :: rmn,zmn,almn     !<=2017 NEW

      common/cflux/flux
c 28.09.2016
      common/csoi/nsoi
!
      allocate(soa(0:kpar))                                             !<=2017 NEW
      allocate(xs(kpar),ys(kpar),as(kpar),bs(kpar),cs(kpar),ds(kpar))   !<=2017 NEW
      allocate(rmn(nstrm,0:kpar),zmn(nstrm,0:kpar),almn(nstrm,0:kpar))  !<=2017 NEW
      allocate(arm(nstrm),arn(nstrm))                                   !<=2017 NEW
      allocate(aiota(0:kpar),phi(0:kpar),sps(0:kpar),s(0:kpar))         !<=2017 NEW
!
      isc=nsoi
c 28.09.2016 end

c      mpas=mpai
c      npas=npai
      npas=nstri
      kpas=kpai

c      mpav=mpar
c      npav=npar
c      kpav=kpar
      nsurfb=nsurfm
      nstrb=nstrm
      kparb=kpar

      do is=0,kpas
      soa(is)=is
      end do

c      call boozin(bmn,rmn,zmn,gmn,aiota,curt,curp,sps,mpav,npav,kpav
c     *,flux1)
      call vmecin(rmn,zmn,almn,aiota,phi,sps,arm,arn,s
     *,nsurfb,nstrb,kparb,flux1)

      flux=flux1

c      do km=0,mpas
            do ks=0,kpas
c            write(94,*)s(ks)
            end do

      do kn=1,npas
      
      axm(kn)=arm(kn)
      axn(kn)=arn(kn)

                alpha=0.5d0*axm(kn)                             !<=NEW

            do ks=1,kpas
            ys(ks)=rmn(kn,ks)
            xs(ks)=soa(ks)
                ys(ks)=ys(ks)/xs(ks)**alpha                     !<=NEW
            end do
            y0=rmn(kn,0)
            x0=soa(0)
                y0=ys(1)+(ys(2)-ys(1))*(x0-xs(1))/(xs(2)-xs(1)) !<=NEW
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            srmn(1,kn,ks)=as(ks)
            srmn(2,kn,ks)=bs(ks)
            srmn(3,kn,ks)=cs(ks)
            srmn(4,kn,ks)=ds(ks)
            end do

            do ks=1,kpas
            ys(ks)=zmn(kn,ks)
            xs(ks)=soa(ks)
                ys(ks)=ys(ks)/xs(ks)**alpha                     !<=NEW
            end do
            y0=zmn(kn,0)
            x0=soa(0)
                y0=ys(1)+(ys(2)-ys(1))*(x0-xs(1))/(xs(2)-xs(1)) !<=NEW
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            szmn(1,kn,ks)=as(ks)
            szmn(2,kn,ks)=bs(ks)
            szmn(3,kn,ks)=cs(ks)
            szmn(4,kn,ks)=ds(ks)
            end do

   91 format(i3,2f7.1,1pe25.16,i3,'  kn,xm,xn,lmn,ks vvovmec')
            do ks=1,kpas
            ys(ks)=almn(kn,ks)
            xs(ks)=soa(ks)
                ys(ks)=ys(ks)/xs(ks)**alpha                     !<=NEW
c 27.09.2016            if(ks.eq.16) then
            if(ks.eq.isc) then
c            write(93,91)kn,axm(kn),axn(kn),almn(kn,ks),ks
            end if
            end do
            y0=almn(kn,0)
            x0=soa(0)
                y0=ys(1)+(ys(2)-ys(1))*(x0-xs(1))/(xs(2)-xs(1)) !<=NEW
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            slmn(1,kn,ks)=as(ks)
            slmn(2,kn,ks)=bs(ks)
            slmn(3,kn,ks)=cs(ks)
            slmn(4,kn,ks)=ds(ks)
            end do

      end do
c      end do
            do ks=1,kpas
            ys(ks)=aiota(ks)
            xs(ks)=soa(ks)
            end do
            y0=aiota(0)
            x0=soa(0)
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            saiota(1,ks)=as(ks)
            saiota(2,ks)=bs(ks)
            saiota(3,ks)=cs(ks)
            saiota(4,ks)=ds(ks)
            end do

            do ks=1,kpas
            ys(ks)=phi(ks)
            xs(ks)=soa(ks)
            end do
            y0=phi(0)
            x0=soa(0)
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            sphi(1,ks)=as(ks)
            sphi(2,ks)=bs(ks)
            sphi(3,ks)=cs(ks)
            sphi(4,ks)=ds(ks)
            end do

            do ks=1,kpas
            ys(ks)=sps(ks)
            xs(ks)=soa(ks)
            end do
            y0=sps(0)
            x0=soa(0)
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            ssps(1,ks)=as(ks)
            ssps(2,ks)=bs(ks)
            ssps(3,ks)=cs(ks)
            ssps(4,ks)=ds(ks)
            end do

            do ks=1,kpas
            ys(ks)=s(ks)
            xs(ks)=soa(ks)
            end do
            y0=s(0)
            x0=soa(0)
            call spltr(kpas,x0,y0,xs,ys,as,bs,cs,ds)
            do ks=1,kpas
            ss(1,ks)=as(ks)
            ss(2,ks)=bs(ks)
            ss(3,ks)=cs(ks)
            ss(4,ks)=ds(ks)
            end do

      print*,flux,'=flux',isc,'=isc,  vvobooz is performed, vvobooz'
      close(71)
      deallocate(soa,xs,ys,as,bs,cs,ds,rmn,zmn,almn,arm,arn,aiota,phi,  !<=2017 NEW
     *                                                           sps,s) !<=2017 NEW
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spltr(nxi,x0i,y0i,xi,yi,ase,bse,cse,dse)
!<=2017 NEW      implicit double precision(a-h,o-z),integer(i-n)
      implicit none                                                     !<=2017 NEW
      integer :: nxi,nx,nx1,nx2,ki                                      !<=2017 NEW
      double precision :: x0i,y0i,det,hi,him1,x0,y0,yim1,yim2           !<=2017 NEW
      double precision, dimension(nxi) :: xi,yi,ase,bse,cse,dse         !<=2017 NEW
      double precision, dimension(:), allocatable :: x,y,a,b,c,d        !<=2017 NEW
      double precision, dimension(:), allocatable :: as,bs,cs,ds        !<=2017 NEW
!<=2017 NEW      parameter(n=253)
!<=2017 NEW      dimension x(n),y(n),a(n),b(n),c(n),d(n),as(n),bs(n),cs(n),ds(n)
!<=2017 NEW      dimension xi(1),yi(1),ase(1),bse(1),cse(1),dse(1)
c
ccc   y=y(x) is an initial graph, ase,bse,cse,dse are the coefficients
ccc   of splain for the initial graph
c
      nx=nxi
      allocate(x(nx),y(nx),a(nx),b(nx),c(nx),d(nx))                     !<=2017 NEW
      allocate(as(nx),bs(nx),cs(nx),ds(nx))                             !<=2017 NEW
cc      if(nx.ge.n) nx=n-1
c add
!<=2017 NEW      if(nx.gt.n) then
!<=2017 NEW      print*,n,'=n',nx,'=nx','  nx>n'
!<=2017 NEW      stop
!<=2017 NEW      end if
c add
      nx1=nx-1
      nx2=nx-2
      x0=x0i
      y0=y0i
ccc      print*,'  nx=',nx ,'  x0=',x0,'  y0=',y0
      do 10 ki=1,nx
      x(ki)=xi(ki)
      y(ki)=yi(ki)
ccc      print*,'  ki=',ki ,'   x=',x(ki),'   y=',y(ki)
   10 continue
      a(1)=0d0
      b(1)=1d0
      c(1)=0d0
      d(1)=0d0
      him1=x(1)-x0
      yim2=y0
      as(1)=y0
      do 20 ki=2,nx
      as(ki)=y(ki-1)
      hi=x(ki)-x(ki-1)
      a(ki)=him1
      b(ki)=-2d0*(him1+hi)
      c(ki)=hi
c      print*,' ki=',ki,' hi=',hi,' him1=',him1
      d(ki)=3d0*((y(ki)-y(ki-1))/hi-(y(ki-1)-yim2)/him1)
      him1=x(ki)-x(ki-1)
      yim2=y(ki-1)
   20 continue
ccc      print*,'  obr k prog'
      call prog(nx,cs,a,b,c,d,det)

ccc      print*,'  prog prorab'
      yim1=y0
      hi=x(1)-x0
      do 30 ki=1,nx1
      ds(ki)=(cs(ki+1)-cs(ki))/(3d0*hi)
      bs(ki)=(y(ki)-yim1)/hi-(cs(ki+1)+2d0*cs(ki))*hi/3d0
      hi=c(ki+1)
      yim1=y(ki)
   30 continue
      ds(nx)=-cs(nx)/(3d0*hi)
      bs(nx)=(y(nx)-y(nx-1))/hi-2d0*cs(nx)*hi/3d0

      do 40 ki=1,nx
      ase(ki)=as(ki)
      bse(ki)=bs(ki)
      cse(ki)=cs(ki)
      dse(ki)=ds(ki)
   40 continue
      deallocate(x,y,a,b,c,d,as,bs,cs,ds)                               !<=2017 NEW
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prog(n,x,a,b,c,d,det)
!      implicit double precision(a-h,o-z),integer(i-n)
      implicit none                                                     !<=2017 NEW
      integer :: n,ni,ni1,ni2,ki                                        !<=2017 NEW
      double precision :: det,aki,akski,bki,cki,dki,detki,etaki,xki,zn  !<=2017 NEW
      double precision, dimension(n) :: x,a,b,c,d                       !<=2017 NEW
      double precision, dimension(:), allocatable :: eta,aks            !<=2017 NEW
!<=2017 NEW      parameter(nsp=253,nmax=nsp+1)
!<=2017 NEW      dimension x(n),a(n),b(n),c(n),d(n),eta(nmax),aks(nmax)
      allocate(eta(n+1),aks(n+1))                                       !<=2017 NEW
      ni=n
cc      if(ni.ge.nmax) ni=nmax-1
c add
!<=2017 NEW      if(ni.gt.nsp) then
!<=2017 NEW      print*,nsp,'=nsp',ni,'=ni','  ni>nsp'
!<=2017 NEW      stop
!<=2017 NEW      end if
c add
      ni1=ni+1
      ni2=ni+2
      akski=0d0
      etaki=0d0
      aks(1)=0d0
      eta(1)=0d0
      aki=0d0
      detki=1d0
      do 10 ki=1,ni
      bki=b(ki)
      cki=c(ki)
      dki=d(ki)
      zn=bki-aki*akski
      akski=cki/zn
      etaki=(aki*etaki-dki)/zn
      detki=-detki*zn
      aks(ki+1)=akski
      eta(ki+1)=etaki
      if(ki.lt.ni) aki=a(ki+1)
   10 continue
      det=detki
      xki=0d0
      do 20 ki=1,ni
      xki=aks(ni2-ki)*xki+eta(ni2-ki)
      x(ni1-ki)=xki
   20 continue
      deallocate(eta,aks)                                               !<=2917 NEW
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE tj2vvo(RT0,R0i,L1i,cbfi,BY0i,bf0)
      double precision rt0,r0i,cbfi,by0i,bf0
c 22.09.2010
c      real*8 facr
c      facr=0.01d0      
c      facr=1.0d0
c 22.09.2010 end
cc
      INTEGER BMMA,BNC
      LOGICAL NOPRI1,NOPRI2
      double precision RC(6),KBZ,rt,r0,by0
c 08.01.2011
      real*8 facr
      facr=0.01d0
c 08.07.2016      facr=1.0d0
c 08.01.2011 end
      BNC=1
      NOPRI1=.false.
      NOPRI2=.FALSE.

      open(17,form='FORMATTED',file='pchsm.d')
      READ(17,340)nopri1,nopri2
  340 FORMAT(7X,L5,6X,7X,L5)

      READ(17,120)NMAX,MMA,NMA,RT,R0,L1,NC
c  120 FORMAT(5X,I2,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
  120 FORMAT(4X,I3,2X,4X,I3,2X,5X,I2,2X,2(3X,E10.4,2X),2(3X,I3,2X))
      IF(NC.GT.BNC) NC=BNC
!<=2017      PRINT 130,NMAX,BMMA,MMA,NMA,L1,RT,R0,NC
  130 FORMAT(//5HNMAX=,I3,4X,4X,5HBMMA=,I3,4X//
     *4HMMA=,I3,4X,4HNMA=,I3,4X,3HL1=,I3,4X,
     *3HRT=,E14.7,4X,3HR0=,E14.7,4X,3HNC=,I3/)
      READ(17,420)KBZ,RC
  420 FORMAT(5X,E13.7/3(3X,E13.7)/3(3X,E13.7))
      READ(17,180)BY0
  180 FORMAT(4X,E14.7)
      PRINT 250,BY0,KBZ
  250 FORMAT(/4HBY0=,E14.7,3X,4HKBZ=,E14.7/)
      PRINT 480,RC
  480 FORMAT(/3HRC=,6(E14.7,2X))
c 20.09.2010      CALL COINhs(RT,R0,RC,NC,BMMA,MMA,NMA,NMAX,
c 20.09.2010     *KBZ,L1,NOPRI1,NOPRI2)

      close(17)
      RT0=rt/facr
      R0i=r0/facr
      L1i=l1
      cbf=kbz/facr
      cbfi=cbf
      BY0i=by0
      bf0=-cbf*L1/rt0

      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MD_CHX
        integer npar,nhist
        parameter(npar=454,nhist=10000)
C
        character*256 fpth,fdeg
        integer nmax,nfile
        double precision T
C
        character*2048 line
        integer path(npar),pmin(npar)
        integer i,j,npath,i1,i2,i0,ngo,nmin,npmin,igo
        double precision enow,lmin,lmin1,lmin2,lmin0,lnow
        double precision q0(npar),q1(npar),qtot
        double precision dl
        double precision deg(nhist),prob(nhist)
        integer imin,il
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,'(a)') fdeg
        read(20,*) nmax
        read(20,*) T
        read(20,*) dl
        read(20,*) nfile
C
        i0=58
C
        do ngo=1,nfile
          read(20,'(a)') fpth
          open(30,FILE=fpth,STATUS='old')
C
          read(30,"(a)") line
          read(line(1:18),*) lmin
          read(line(i0:i0+2),*) npmin
          do j=1,npmin
            read(line(i0+j*4:i0+2+j*4),*) pmin(j)
          enddo
          lmin2=2.d0*(lmin+1.d0)
          nmin=1
          do i=2,nmax
            read(30,"(a)") line
            read(line(1:18),*) lnow
            if (lnow.lt.lmin2) then
              read(line(i0:i0+2),*) npath
              do j=1,npath
                read(line(i0+j*4:i0+2+j*4),*) path(j)
              enddo
              if(npmin.eq.npath) then
                igo=1
                do j=1,npath
                  if(pmin(j).ne.path(j)) igo=0
                enddo
                if(igo.eq.1) then
                  nmin=nmin+1
                else
                  if(lnow.lt.lmin) then
                    npmin=npath
                    do j=1,npath
                      pmin(j)=path(j)
                    enddo
                    lmin2=lmin
                    lmin=lnow
                    nmin=1
                  else
                    lmin2=lnow
                  endif
                endif
              else
                if(lnow.lt.lmin) then
                  npmin=npath
                  do j=1,npath
                    pmin(j)=path(j)
                  enddo
                  lmin2=lmin
                  lmin=lnow
                  nmin=1
                else
                  lmin2=lnow
                endif
              endif
            endif
          enddo
          close(30)
          if(ngo.eq.1) then
            lmin1=lmin
            lmin0=lmin
          elseif(lmin.lt.lmin0) then
            lmin0=lmin
          endif
          q0(ngo)=dble(nmax)/dble(nmin)*dexp((lmin1-lmin)/T)
        enddo
        close(20)
C
        qtot=0.d0
        do i=1,nfile
          qtot=qtot+q0(i)
        enddo
        qtot=qtot*dble(nmax)*dl
        do i=1,nfile
          q1(i)=q0(i)/dble(nmax)
          q0(i)=q0(i)/qtot
        enddo
        imin=int(lmin0/dl)-1
C
        do i=1,nhist
          prob(i)=0.d0
          deg(i)=0.d0
        enddo
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,*) 
        read(20,*) 
        read(20,*) 
        read(20,*) 
        read(20,*) 
        do ngo=1,nfile
          read(20,'(a)') fpth
          open(30,FILE=fpth,STATUS='old')
          do i=1,nmax
            read(30,*) lnow
            il=int(lnow/dl)-imin
            if(il.lt.nhist) then
              if (il.le.0) write(0,*) 'shit',il
              prob(il)=prob(il)+q0(ngo)
              deg(il)=deg(il)+q1(ngo)*dexp((lnow-lmin1)/T)
            else
              prob(nhist)=prob(nhist)+q0(ngo)
              deg(nhist)=deg(nhist)+q1(ngo)*dexp((lnow-lmin1)/T)
            endif
          enddo
          close(30)
        enddo
        close(20)
C
        open(20,FILE=fdeg,STATUS='unknown')
        do i=1,nhist
          write(20,*) (dble(i+imin)-0.5d0)*dl,prob(i),deg(i)
        enddo
        close(20)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

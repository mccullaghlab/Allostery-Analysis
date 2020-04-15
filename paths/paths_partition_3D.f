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
        integer i,j,npath,i1,i2,i0,ngo,nmin,npmin,igo,n
        double precision enow,lmin,lmin1,lmin2,lmin0,lnow
        double precision q(npar),u(npar),s(npar),a(npar)
        double precision dl
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
          u(ngo)=0.d0
          read(20,'(a)') fpth
          open(30,FILE=fpth,STATUS='old')
C
          read(30,"(a)") line
         read(line(1:18),*) lmin
          u(ngo)=u(ngo)+lmin
          read(line(i0:i0+2),*) npmin
          do j=1,npmin
            read(line(i0+j*4:i0+2+j*4),*) pmin(j)
          enddo
          lmin2=2.d0*(lmin+1.d0)
          nmin=1
          do i=2,nmax
            read(30,"(a)") line
            read(line(1:18),*) lnow
            u(ngo)=u(ngo)+lnow
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
          q(ngo)=dble(nmax)/dble(nmin)*dexp((lmin1-lmin)/T)
          u(ngo)=u(ngo)/dble(nmax)
          a(ngo)=-T*dlog(q(ngo))+lmin1
          s(ngo)=u(ngo)-a(ngo)
        enddo
        close(20)
C
        open(20,FILE=fdeg,STATUS='unknown')
C
        n=0
        write(20,*) 'partition:'
        do i=1,4
          do j=1,4
            n=n+1
            write(20,999) q(n)*dexp(-lmin1/T)
          enddo
          write(20,*) ''
        enddo
        write(20,*) ''
C
        n=0
        write(20,*) 'Hemholtz free energy:'
        do i=1,4
          do j=1,4
            n=n+1
            write(20,999) a(n)
          enddo
          write(20,*) ''
        enddo
        write(20,*) ''
C
        n=0
        write(20,*) 'Internal energy:'
        do i=1,4
          do j=1,4
            n=n+1
            write(20,999) u(n)
          enddo
          write(20,*) ''
        enddo
        write(20,*) ''
C
        n=0
        write(20,*) 'Entropy (TS):'
        do i=1,4
          do j=1,4
            n=n+1
            write(20,999) s(n)
          enddo
          write(20,*) ''
        enddo
        write(20,*) ''
        close(20)
C
999     format(e25.18,' ',$)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

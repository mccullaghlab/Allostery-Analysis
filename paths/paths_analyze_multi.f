CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MD_CHX
        integer npar
        parameter(npar=454)
C
        character*256 fpth,fnod,fedg
        integer nmax,nfile
        double precision T
C
        character*2048 line
        integer path(npar),pmin(npar)
        integer i,j,npath,i1,i2,i0,ngo,nmin,npmin,igo
        double precision enow,lmin,lmin2,lmin0,lnow
        double precision node(npar)
        double precision edge(npar,npar)
        double precision q0(npar),qtot
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,'(a)') fnod
        read(20,'(a)') fedg
        read(20,*) nmax
        read(20,*) T
        read(20,*) nfile
C
        i0=20
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
          if(ngo.eq.1) lmin0=lmin
          q0(ngo)=dble(nmax)/dble(nmin)*dexp((lmin0-lmin)/T)
        enddo
        close(20)
C
        qtot=0.d0
        do i=1,nfile
          qtot=qtot+q0(i)
        enddo
        qtot=qtot*dble(nmax)
        do i=1,nfile
          q0(i)=q0(i)/qtot
        enddo
C
        do i=1,npar
          node(i)=0.d0
          do j=1,npar
            edge(j,i)=0.d0
          enddo
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
            read(30,"(a)") line
            read(line(i0:i0+2),*) npath
            do j=1,npath
              read(line(i0+j*4:i0+2+j*4),*) path(j)
            enddo
            i1=path(1)
            node(i1)=node(i1)+q0(ngo)
            do j=2,npath
              i2=path(j)
              node(i2)=node(i2)+q0(ngo)
              edge(i2,i1)=edge(i2,i1)+q0(ngo)
              i1=i2
            enddo
          enddo
          close(30)
        enddo
        close(20)
C
        open(20,FILE=fnod,STATUS='unknown')
        do i=1,npar
          write(20,*) i,node(i)
        enddo
        close(20)
C
        open(20,FILE=fedg,STATUS='unknown')
        do i=1,npar
          do j=i+1,npar
            enow=edge(j,i)+edge(i,j)
            if(enow.gt.1.d-16) then
             write(20,*) i,j,enow
            endif
          enddo
        enddo
        close(20)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

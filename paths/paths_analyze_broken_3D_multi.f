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
        integer i,j,npath,i1,i2,j1,j2,i0,ngo,nmin,npmin,igo
        double precision enow,lmin,lmin2,lmin0,lnow,enow1,enow2,enow3
        double precision node(2,npar)
        double precision edge(3,npar,npar)
        double precision q0(npar),qtot
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,'(a)') fnod
        read(20,'(a)') fedg
        read(20,*) nmax
        read(20,*) T
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
          read(line(i0:i0+3),*) npmin
          do j=1,npmin
            read(line(i0+j*5:i0+3+j*5),*) pmin(j)
          enddo
          lmin2=2.d0*(lmin+1.d0)
          nmin=1
          do i=2,nmax
            read(30,"(a)") line
            read(line(1:18),*) lnow
            if (lnow.lt.lmin2) then
              read(line(i0:i0+3),*) npath
              do j=1,npath
                read(line(i0+j*5:i0+3+j*5),*) path(j)
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
          node(1,i)=0.d0
          node(2,i)=0.d0
          do j=1,npar
            edge(1,j,i)=0.d0
            edge(2,j,i)=0.d0
            edge(3,j,i)=0.d0
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
            read(line(i0:i0+3),*) npath
            do j=1,npath
              read(line(i0+j*5:i0+3+j*5),*) path(j)
            enddo
            i1=iabs(path(1))
            j1=(1-path(1)/i1)/2+1
            node(j1,i1)=node(j1,i1)+q0(ngo)
            do j=2,npath
              i2=iabs(path(j))
              j2=(1-path(j)/i2)/2+1
              node(j2,i2)=node(j2,i2)+q0(ngo)
              if(j1.eq.j2) then
                edge(j2,i2,i1)=edge(j2,i2,i1)+q0(ngo)
              else
                edge(3,i2,i1)=edge(3,i2,i1)+q0(ngo)
              endif
              i1=i2
              j1=j2
            enddo
          enddo
          close(30)
        enddo
        close(20)
C
        open(20,FILE=fnod,STATUS='unknown')
        do i=1,npar
          write(20,*) i,node(1,i),node(2,i)
        enddo
        close(20)
C
        open(20,FILE=fedg,STATUS='unknown')
        do i=1,npar
          do j=i+1,npar
            enow1=edge(1,j,i)+edge(1,i,j)
            enow2=edge(2,j,i)+edge(2,i,j)
            enow3=edge(3,j,i)+edge(3,i,j)
            enow=enow1+enow2+enow3
            if(enow.gt.1.d-16) then
             write(20,*) i,j,enow1,enow2,enow3
            endif
          enddo
        enddo
        close(20)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

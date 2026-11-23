CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MD_CHX
        integer npar
        parameter(npar=454)
C
        character*256 fpth,fout
        integer nmax,nfile
C
        character*1024 line
        integer path(npar)
        integer plist(npar,1000),nlist(1000)
        double precision dlist(1000)
        double precision dmax,dlen,dmin
        integer i,j,k,npath,i1,i2,imax,inew,jnew,imin,ngo
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,'(a)') fout
        read(20,*) nmax
        read(20,*) nfile
C
        do i=1,1000
          dlist(i)=1.d99
          plist(1,i)=0
          nlist(i)=0
        enddo
        dmax=dlist(1)
        imax=1
C
        do ngo=1,nfile
          read(20,'(a)') fpth
          open(30,FILE=fpth,STATUS='old')
          do i=1,nmax
            read(30,"(a)") line
            read(line(1:18),*) dlen
            read(line(20:23),*) npath
            do j=1,npath
              read(line(20+j*5:23+j*5),*) path(j)
            enddo
C
            if(dlen.le.dmax) then
              inew=1
              j=1
              do while(inew.eq.1)
                if(j.eq.1000) inew=2
                if (plist(1,j).eq.npath) then
                  jnew=1
                  k=1
                  do while(jnew.eq.1)
                    if (k.eq.npath) jnew=0
                    if(path(k).ne.plist(k+1,j)) jnew=2
                    k=k+1
                  enddo
                  if(jnew.eq.0) then
                    inew=0
                    nlist(j)=nlist(j)+1
                  endif
                endif
                j=j+1
              enddo
              if (inew.eq.2) then
                plist(1,imax)=npath
                do j=1,npath
                  plist(j+1,imax)=path(j)
                enddo
                dlist(imax)=dlen
                nlist(imax)=1
C
                dmax=dlist(1)
                imax=1
                do j=2,1000
                  if (dlist(j).gt.dmax) then
                    dmax=dlist(j)
                    imax=j
                  endif
                enddo
              endif
            endif
          enddo
          close(30)
        enddo
        close(20)
C
        open(20,FILE=fout,STATUS='unknown')
        do i=1,1000
          dmin=dlist(1)
          imin=1
          do j=2,1000
            if(dlist(j).lt.dmin) then
              dmin=dlist(j)
              imin=j
            endif
          enddo
C
          write(20,888) dlist(imin),nlist(imin),plist(1,imin)
          do k=2,plist(1,imin)+1
            write(20,887) plist(k,imin)
          enddo
          write(20,*) ''
888       format(e18.12,1X,i6,1X,i4,$)
887       format(1X,i4,$)
C
          dlist(imin)=1.d99
        enddo
        close(20)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

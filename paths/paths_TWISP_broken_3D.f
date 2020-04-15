CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MD_CHX
        integer npar
        parameter(npar=454)
C
        integer path(npar),ipath(npar)
C        integer minp(npar),minn
        double precision x(3,npar)
        double precision rhat(3,npar,npar)
        double precision minl
        integer ntry(npar)
        double precision wtry(npar)
        double precision adj(npar,npar)
        double precision anow,ran3,try,prob,wgt,dwgt,T,bet,wgt2
        integer i,j,k,isrc,isnk,i1,j1,n,npath,nswp,ngo,mgo
        integer ncon,i2,igo,jgo,nmax,i3,i4
        integer j2,j3,j4,mswp
        double precision wtot,wmax,x2,dx
        character*11349 line
        character*256 fout,fadj,fxyz
C random
        integer ma(55),inext,inextp
C
        common /random/ ma,inext,inextp
C  Get random number seed
        open(20,FILE='ran3.dat',STATUS='old')
        read(20,*) inext,inextp
        do i=1,55
          read(20,*) ma(i)
        enddo
        close(20)
C
        open(20,FILE='input.dat',STATUS='old')
        read(20,'(a)') fadj
        read(20,'(a)') fxyz
        read(20,'(a)') fout
        read(20,*) isrc
        read(20,*) isnk
        read(20,*) nmax
        read(20,*) T
        close(20)
        bet=1.d0/T
C
        open(20,FILE=fadj,STATUS='old')
        do i=1,npar
          read(20,'(a11349)') line
          do j=1,npar
            read(line(25*j-24:25*j),'(e24.18)') anow
            if(anow.lt.1.d-10) then
              adj(j,i)=-1.d0
            else
              adj(j,i)=-dlog(anow)
            endif
          enddo
        enddo
        close(20)
C  define starting path
        npath=2
        path(1)=isrc
        path(2)=-isnk
        do i=1,npar
          ipath(i)=0
        enddo
        ipath(isrc)=1
        ipath(isnk)=-2
C     
        open(20,FILE=fxyz,STATUS='old')
        do i=1,npar
          read(20,*)  x(1,i),x(2,i),x(3,i)
        enddo
        close(20)
C
        do i=1,npar
          do j=1,npar
            if(i.ne.j) then
              x2=0.d0
              do k=1,3
                dx=x(k,j)-x(k,i)
                rhat(k,j,i)=dx
                x2=x2+dx*dx
              enddo
              dx=x2**0.5d0
              rhat(1,j,i)=rhat(1,j,i)/dx
              rhat(2,j,i)=rhat(2,j,i)/dx
              rhat(3,j,i)=rhat(3,j,i)/dx
            endif
          enddo
        enddo
C
        open(20,FILE=fout,STATUS='unknown')
        do ngo=1,nmax
          do mgo=1,npar
            n=int(ran3()*dble(npar))+1
            do while (n.eq.isrc.or.n.eq.isnk)
              n=int(ran3()*dble(npar))+1
            enddo
C Calculate ncon
            ncon=0
            i1=iabs(path(1))
            j1=path(1)/i1
            jgo=1
            wmax=-1.d12
            do i=2,npath
              i2=iabs(path(i))
              j2=path(i)/i2
              if(j1.eq.j2) then
                if(i2.ne.n) then
                  if(adj(n,i1).gt.-0.1d0.and.adj(n,i2).gt.-0.1d0) then
                    ncon=ncon+1
                    ntry(ncon)=j1*i
                    wtry(ncon)=adj(i1,i2)-adj(n,i1)-adj(n,i2)
                    if(wtry(ncon).gt.wmax) then
                      wmax=wtry(ncon)
                    endif
                  endif
                  i1=i2
                  j1=j2
                else
                  i2=iabs(path(i+1))
                  j2=path(i+1)/i2
                  if(j1.eq.j2.and.adj(i1,i2).lt.-0.1d0) jgo=0
                endif
              else
                if(i2.ne.n) then
                  if(adj(n,i1).gt.-0.1d0) then
                    ncon=ncon+1
                    ntry(ncon)=j1*i
                    wtry(ncon)=-adj(n,i1)
                    if(wtry(ncon).gt.wmax) then
                      wmax=wtry(ncon)
                    endif
                  endif
                  if(adj(n,i2).gt.-0.1d0) then
                    ncon=ncon+1
                    ntry(ncon)=j2*i
                    wtry(ncon)=-adj(n,i2)
                    if(wtry(ncon).gt.wmax) then
                      wmax=wtry(ncon)
                    endif
                  endif
                  i1=i2
                  j1=j2
                endif
              endif
            enddo
            wtot=0.d0
            do i=1,ncon
              wtot=wtot+dexp((wtry(i)-wmax)*bet)
              wtry(i)=wtot
            enddo
C add to path
            if(ipath(n).eq.0) then
              if(ncon.ne.0) then
                prob=dexp(wmax*bet)*wtot
                if(prob.lt.1.d0) then
                  try=ran3()
                  if (try.lt.prob) then
                    igo=1
                  else
                    igo=0
                  endif
                else
                  igo=1
                endif
                if(igo.eq.1) then
                  try=ran3()*wtot
                  nswp=1
                  do while (try.gt.wtry(nswp))
                    nswp=nswp+1
                  enddo
                  if (nswp.gt.ncon) write(0,*) 'shit',ncon,nswp
                  nswp=ntry(nswp)
                  mswp=nswp/iabs(nswp)
                  nswp=iabs(nswp)
C angle metropolis
                  i2=iabs(path(nswp-1))
                  j2=path(nswp-1)/i2
                  i3=iabs(path(nswp))
                  j3=path(nswp)/i3
                  if(j2.eq.j3) then
                    prob=rhat(1,i2,n)*rhat(1,n,i3)+
     &                   rhat(2,i2,n)*rhat(2,n,i3)+
     &                   rhat(3,i2,n)*rhat(3,n,i3)
                  else
                    prob=1.d0
                  endif
                  if(nswp.ne.2) then
                    i1=iabs(path(nswp-2))
                    j1=path(nswp-2)/i1
                    if(j1.eq.mswp) then
                      prob=prob*(rhat(1,i1,i2)*rhat(1,i2,n )+
     &                           rhat(2,i1,i2)*rhat(2,i2,n )+
     &                           rhat(3,i1,i2)*rhat(3,i2,n ))
                    endif
                    if(j1.eq.j3) then
                      prob=prob/(rhat(1,i1,i2)*rhat(1,i2,i3)+
     &                           rhat(2,i1,i2)*rhat(2,i2,i3)+
     &                           rhat(3,i1,i2)*rhat(3,i2,i3))
                    endif
                  endif
                  if(nswp.ne.npath) then
                    i4=iabs(path(nswp+1))
                    j4=path(nswp+1)/i4
                    if(mswp.eq.j4) then
                      prob=prob*(rhat(1,n ,i3)*rhat(1,i3,i4)+
     &                           rhat(2,n ,i3)*rhat(2,i3,i4)+
     &                           rhat(3,n ,i3)*rhat(3,i3,i4))
                    endif
                    if(j2.eq.j4) then
                      prob=prob/(rhat(1,i2,i3)*rhat(1,i3,i4)+
     &                           rhat(2,i2,i3)*rhat(2,i3,i4)+
     &                           rhat(3,i2,i3)*rhat(3,i3,i4))
                    endif
                  endif
C
                  prob=dabs(prob)**bet
                  if(prob.lt.1.d0) then
                    try=ran3()
                    if(try.lt.prob) then
                      igo=1
                    else
                      igo=0
                    endif
                  else
                    igo=1
                  endif
C swap
                  if(igo.eq.1) then
                    do i=npath,nswp,-1
                      j=iabs(path(i))
                      k=path(i)/j
                      ipath(j)=k*(i+1)
                      path(i+1)=k*j
                    enddo
                    path(nswp)=mswp*n
                    ipath(n)=mswp*nswp
                    npath=npath+1
                  endif
                endif
              endif
C remove from path
            elseif(jgo.eq.1) then
              prob=dexp(-wmax*bet)/wtot
              if (prob.lt.1.d0) then
                try=ran3()
                if(try.lt.prob) then
                  igo=1
                else
                  igo=0
                endif
              else
                igo=1
              endif
              if (igo.eq.1) then
                nswp=iabs(ipath(n))
                mswp=ipath(n)/nswp
C angle metropolis
                i2=iabs(path(nswp-1))
                j2=path(nswp-1)/i2
                i3=iabs(path(nswp+1))
                j3=path(nswp+1)/i3
                if(j2.eq.j3) then
                  prob=rhat(1,i2,n)*rhat(1,n,i3)+
     &                 rhat(2,i2,n)*rhat(2,n,i3)+
     &                 rhat(3,i2,n)*rhat(3,n,i3)
                else
                  prob=1.d0
                endif
                if(nswp.ne.2) then
                  i1=iabs(path(nswp-2))
                  j1=path(nswp-2)/i1
                  if(j1.eq.mswp) then
                    prob=prob*(rhat(1,i1,i2)*rhat(1,i2,n )+
     &                         rhat(2,i1,i2)*rhat(2,i2,n )+
     &                         rhat(3,i1,i2)*rhat(3,i2,n ))
                  endif
                  if(j1.eq.j3) then
                    prob=prob/(rhat(1,i1,i2)*rhat(1,i2,i3)+
     &                         rhat(2,i1,i2)*rhat(2,i2,i3)+
     &                         rhat(3,i1,i2)*rhat(3,i2,i3))
                  endif
                endif
                if(nswp.ne.npath-1) then
                  i4=iabs(path(nswp+2))
                  j4=path(nswp+2)/i4
                  if(mswp.eq.j4) then
                    prob=prob*(rhat(1,n ,i3)*rhat(1,i3,i4)+
     &                         rhat(2,n ,i3)*rhat(2,i3,i4)+
     &                         rhat(3,n ,i3)*rhat(3,i3,i4))
                  endif
                  if(j2.eq.j4) then
                    prob=prob/(rhat(1,i2,i3)*rhat(1,i3,i4)+
     &                         rhat(2,i2,i3)*rhat(2,i3,i4)+
     &                         rhat(3,i2,i3)*rhat(3,i3,i4))
                  endif
                endif
C
                prob=dabs(prob)**(-bet)
                if(prob.lt.1.d0) then
                  try=ran3()
                  if(try.lt.prob) then
                    igo=1
                  else
                    igo=0
                  endif
                else
                  igo=1
                endif
Cswap
                if(igo.eq.1) then
                  ipath(n)=0
                  do i=nswp+1,npath
                    j=iabs(path(i))
                    k=path(i)/j
                    path(i-1)=k*j
                    ipath(j)=k*(i-1)
                  enddo
                  npath=npath-1
                endif
              endif
            endif
          enddo
C wgt
          i1=iabs(path(1))
          j1=path(1)/i1
          wgt=0.d0
          do i=2,npath
            i2=iabs(path(i))
            j2=path(i)/i2
            if(j1.eq.j2) then
              wgt=wgt+adj(i1,i2)
            endif
            i1=i2
            j1=j2
          enddo
          i1=iabs(path(1))
          j1=path(1)/i1
          i2=iabs(path(2))
          j2=path(2)/i2
          wgt2=0.d0
          do i=3,npath
            i3=iabs(path(i))
            j3=path(i)/i3
            if(j1.eq.j3) then
              wgt2=wgt2-dlog(dabs(rhat(1,i1,i2)*rhat(1,i2,i3)+
     &                            rhat(2,i1,i2)*rhat(2,i2,i3)+
     &                            rhat(3,i1,i2)*rhat(3,i2,i3)))
            endif
            i1=i2
            j1=j2
            i2=i3
            j2=j3
          enddo
          write(20,888) wgt+wgt2,wgt,wgt2,npath
          do i=1,npath
            write(20,887) path(i)
          enddo
          write(20,*) ''
888       format(3(e18.12,1X),i4,$)
887       format(1X,i4,$)
        enddo
        close(20)
C
C  Put random number seed
        open(20,FILE='ran3.dat',STATUS='unknown')
        write(20,*) inext,inextp
        do i=1,55
          write(20,*) ma(i)
        enddo
        close(20)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function ran3()
        integer mbig,mseed,mz,mbig1,inext,inextp,mj
        double precision fac
        parameter(mbig=1000000000,mseed=161803398,mz=0)
        parameter(mbig1=mbig-1,fac=1.d0/mbig)
        integer ma(55)
C
        common /random/ ma,inext,inextp
C$OMP threadprivate(/random/)
C
        inext=inext+1
        if (inext.eq.56) inext=1
        inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if (mj.le.mz)mj=mj+mbig1
        ma(inext)=mj
        ran3=mj*fac
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

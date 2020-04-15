CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MD_CHX
        integer npar
        parameter(npar=454)
C
        integer path(npar),ipath(npar)
C        integer minp(npar),minn
        double precision x(3,npar)
        double precision minl
        integer ntry(npar)
        double precision wtry(npar)
        double precision adj(npar,npar)
        double precision anow,ran3,try,prob,wgt,dwgt,T
        integer i,j,k,isrc,isnk,i1,j1,n,npath,nswp,ngo,mgo
        integer ncon,i2,igo,jgo,nmax
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
C
        call dijkstra(adj,path,ipath,npath,isrc,isnk,minl)
        write(0,*) minl,npath
        write(0,*) path(1)
        i1=path(1)
        do i=2,npath
          i2=path(i)
          write(0,*) i2,adj(i1,i2)
          i1=i2
        enddo

C
        open(20,FILE=fxyz,STATUS='old')
        do i=1,npar
          read(20,*)  x(1,i),x(2,i),x(3,i)
        enddo
        close(20)
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
            i1=path(1)
            jgo=1
            wmax=-1.d12
            do i=2,npath
              i2=path(i)
              if(i2.ne.n) then
                if(adj(n,i1).gt.-0.1d0.and.adj(n,i2).gt.-0.1d0) then
                  ncon=ncon+1
                  ntry(ncon)=i
                  wtry(ncon)=adj(i1,i2)-adj(n,i1)-adj(n,i2)
                  if(wtry(ncon).gt.wmax) then
                    wmax=wtry(ncon)
                  endif
                endif
                i1=i2
              else
                i2=path(i+1)
                if(adj(i1,i2).lt.-0.1d0) jgo=0
              endif
            enddo
            wtot=0.d0
            do i=1,ncon
              wtot=wtot+dexp((wtry(i)-wmax)/T)
              wtry(i)=wtot
            enddo
C add to path
            if(ipath(n).eq.0) then
C swap
              if(ncon.ne.0) then
                prob=dexp(wmax/T)*wtot
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
                  do i=npath,nswp,-1
                    j=path(i)
                    ipath(j)=i+1
                    path(i+1)=j
                  enddo
                  path(nswp)=n
                  ipath(n)=nswp
                  npath=npath+1
                endif
              endif
C remove from path
            elseif(jgo.eq.1) then
              prob=dexp(-wmax/T)/wtot
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
Cswap
              if (igo.eq.1) then
                nswp=ipath(n)
                ipath(n)=0
                do i=nswp+1,npath
                  j=path(i)
                  path(i-1)=j
                  ipath(j)=i-1
                enddo
                npath=npath-1
              endif
            endif
          enddo
C wgt
          i1=path(1)
          wgt=0.d0
          do i=2,npath
            i2=path(i)
            wgt=wgt+adj(i1,i2)
            i1=i2
          enddo
          write(20,888) wgt,npath
          do i=1,npath
            write(20,887) path(i)
          enddo
          write(20,*) ''
888       format(e18.12,1X,i3,$)
887       format(1X,i3,$)
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

      subroutine dijkstra(adj,path,ipath,npath,isrc,isnk,mdis)
        integer npar
        parameter(npar=454)
        double precision adj(npar,npar)
        integer path(npar),ipath(npar)
        integer isrc,isnk,npath
        double precision dist(npar)
        double precision mdis,alt
        integer prev(npar),q(npar)
        integer i,mpar,j,n,m
C
        mpar=npar
        do i=1,npar
          q(i)=i
          dist(i)=9.d99
          prev(i)=0
        enddo
        dist(isrc)=0.d0

        do while(mpar.ne.0)
          n=1
          m=q(n)
          mdis=dist(m)
          do i=2,mpar
            j=q(i)
            if(dist(j).lt.mdis) then
              mdis=dist(j)
              n=i
              m=j
            endif
          enddo
C
          if(m.eq.isnk) then
            mpar=0
          else
            q(n)=q(mpar)
            mpar=mpar-1
C
            do i=1,npar
              if(i.ne.m) then
                if(adj(m,i).gt.0.d0) then
                  alt=mdis+adj(m,i)
                  if (alt.lt.dist(i)) then
                    dist(i)=alt
                    prev(i)=m
                  endif
                endif
              endif
            enddo
          endif
        enddo
C
        npath=0
        m=isnk
        do while (m.ne.0)
          npath=npath+1
          path(npath)=m
          m=prev(m)
        enddo
C
        do i=1,npath/2
          j=npath+1-i
          n=path(i)
          path(i)=path(j)
          path(j)=n
        enddo
C
        do i=1,npar
          ipath(i)=0
        enddo
        do i=1,npath
          j=path(i)
          ipath(j)=i
        enddo
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

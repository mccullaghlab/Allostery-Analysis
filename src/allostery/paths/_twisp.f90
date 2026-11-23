!=======================================================================
! Monte Carlo sampling of allosteric paths on a weighted graph.
!
! This is a consolidation of the original fixed-size programs
!     paths_TWISP.f, paths_TWISP_3D.f, paths_TWISP_broken.f and
!     paths_TWISP_broken_3D.f
! into a single library routine callable from Python (via ctypes).
! The control flow and the order of random-number calls are identical
! to the original programs, so for the same RNG state the same paths
! are produced (see tests/test_legacy_regression.py).
!
! Differences from the original programs:
!   * no hard-coded system size (npar) and no file I/O;
!   * the RNG state (Numerical Recipes ran3) is passed explicitly, so
!     several independent samplers may run concurrently;
!   * sampling is done in chunks: the current path is passed in/out
!     and the routine stops early when the output buffer is full.
!
! Conventions (all arrays are passed from C-ordered numpy arrays):
!   cost(j,i) = C[i,j]  edge cost -ln(A[i,j]) of the normalized adjacency A,
!                       -1 for absent edges (computed by the Python caller,
!                       so that every backend uses bit-identical costs)
!   xyz(k,i)  = X[i,k]  node positions (used only if angular /= 0)
!   path      = 1-based node indices.  For broken paths the sign marks
!               the segment: +  source segment,  -  sink segment.
!=======================================================================
module allostery_twisp
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   implicit none
   private
   public :: twisp_sample, twisp_dijkstra

   integer(c_int), parameter :: MBIG = 1000000000_c_int
   integer(c_int), parameter :: MZ = 0_c_int
   integer(c_int), parameter :: MBIG1 = MBIG - 1_c_int
   real(c_double), parameter :: FAC = 1.0d0/MBIG

contains

   !--------------------------------------------------------------------
   ! Numerical Recipes ran3 (as modified in the original code).
   ! state(1:55) = ma, state(56) = inext, state(57) = inextp
   !--------------------------------------------------------------------
   function ran3(state) result(r)
      integer(c_int), intent(inout) :: state(57)
      real(c_double) :: r
      integer(c_int) :: mj

      state(56) = state(56) + 1
      if (state(56) == 56) state(56) = 1
      state(57) = state(57) + 1
      if (state(57) == 56) state(57) = 1
      mj = state(state(56)) - state(state(57))
      if (mj <= MZ) mj = mj + MBIG1
      state(state(56)) = mj
      r = mj*FAC
   end function ran3


   !--------------------------------------------------------------------
   ! O(n^2) Dijkstra shortest path from isrc to isnk (1-based) on the
   ! cost graph; identical to the subroutine of paths_TWISP*.f.
   ! On return path(1:npath) holds the path and mdis its length;
   ! mdis >= 9e99 means that the sink is unreachable.
   !--------------------------------------------------------------------
   subroutine twisp_dijkstra(n, cost, isrc, isnk, path, npath, mdis) &
      bind(C, name="allostery_twisp_dijkstra")
      integer(c_int), value, intent(in) :: n, isrc, isnk
      real(c_double), intent(in) :: cost(n, n)
      integer(c_int), intent(out) :: path(n)
      integer(c_int), intent(out) :: npath
      real(c_double), intent(out) :: mdis

      real(c_double), allocatable :: dist(:)
      integer(c_int), allocatable :: prev(:), q(:)
      real(c_double) :: alt
      integer :: i, mpar, j, k, m

      allocate (dist(n), prev(n), q(n))

      mpar = n
      do i = 1, n
         q(i) = i
         dist(i) = 9.d99
         prev(i) = 0
      end do
      dist(isrc) = 0.d0

      do while (mpar /= 0)
         k = 1
         m = q(k)
         mdis = dist(m)
         do i = 2, mpar
            j = q(i)
            if (dist(j) < mdis) then
               mdis = dist(j)
               k = i
               m = j
            end if
         end do

         if (m == isnk) then
            mpar = 0
         else
            q(k) = q(mpar)
            mpar = mpar - 1
            do i = 1, n
               if (i /= m) then
                  if (cost(m, i) > 0.d0) then
                     alt = mdis + cost(m, i)
                     if (alt < dist(i)) then
                        dist(i) = alt
                        prev(i) = m
                     end if
                  end if
               end if
            end do
         end if
      end do

      npath = 0
      m = isnk
      do while (m /= 0)
         npath = npath + 1
         path(npath) = m
         m = prev(m)
      end do
      do i = 1, npath/2
         j = npath + 1 - i
         k = path(i)
         path(i) = path(j)
         path(j) = k
      end do
      deallocate (dist, prev, q)
   end subroutine twisp_dijkstra

   !--------------------------------------------------------------------
   ! cosine of the angle between bonds (a->b) and (b->c)
   !--------------------------------------------------------------------
   pure function bond_cos(rhat, n, a, b, c) result(v)
      integer(c_int), intent(in) :: n, a, b, c
      real(c_double), intent(in) :: rhat(3, n, n)
      real(c_double) :: v
      v = rhat(1, a, b)*rhat(1, b, c) + &
          rhat(2, a, b)*rhat(2, b, c) + &
          rhat(3, a, b)*rhat(3, b, c)
   end function bond_cos

   !--------------------------------------------------------------------
   ! Sample paths.
   !
   ! n            number of nodes
   ! cost(n,n)    edge costs (see the module header)
   ! xyz(3,n)     node coordinates (referenced only if angular /= 0)
   ! angular      0: path weight = sum of edge costs ("TWISP")
   !              1: add angular term -sum log|cos(theta)| ("3D")
   ! (broken paths are selected by the signs of the initial path:
   !  e.g. path = (isrc, -isnk) starts a broken-path simulation)
   ! isrc, isnk   1-based source / sink
   ! temperature  MC temperature
   ! nsweeps      maximum number of paths to generate in this call
   ! path, npath  current path (in/out)
   ! state(57)    ran3 state (in/out)
   ! cap          size of nodes_out
   ! nodes_out    concatenated output paths
   ! lengths_out  number of nodes of each output path
   ! weights_out  (total, edge, angular) weight of each output path
   ! ndone        number of paths generated
   ! nnodes       number of entries written to nodes_out
   ! status       0 ok, 1 inconsistent insertion index encountered
   !--------------------------------------------------------------------
   subroutine twisp_sample(n, cost, xyz, angular, isrc, isnk, &
                           temperature, nsweeps, path, npath, state, &
                           cap, nodes_out, lengths_out, weights_out, &
                           ndone, nnodes, status) &
      bind(C, name="allostery_twisp_sample")
      integer(c_int), value, intent(in) :: n, angular
      integer(c_int), value, intent(in) :: isrc, isnk, nsweeps, cap
      real(c_double), value, intent(in) :: temperature
      real(c_double), intent(in) :: cost(n, n)
      real(c_double), intent(in) :: xyz(3, *)
      integer(c_int), intent(inout) :: path(n)
      integer(c_int), intent(inout) :: npath
      integer(c_int), intent(inout) :: state(57)
      integer(c_int), intent(out) :: nodes_out(cap)
      integer(c_int), intent(out) :: lengths_out(nsweeps)
      real(c_double), intent(out) :: weights_out(3, nsweeps)
      integer(c_int), intent(out) :: ndone, nnodes, status

      real(c_double), allocatable :: rhat(:, :, :), wtry(:)
      integer(c_int), allocatable :: ipath(:), ntry(:)
      real(c_double) :: T, bet, wtot, wmax, prob, try, wgt, wgt2, x2, dx
      integer(c_int) :: i, j, k, nn, i1, i2, i3, i4, j1, j2, j3, j4
      integer(c_int) :: ngo, mgo, ncon, nswp, mswp, igo, jgo
      logical :: is3d

      status = 0
      ndone = 0
      nnodes = 0
      T = temperature
      bet = 1.d0/T
      is3d = (angular /= 0)

      allocate (ipath(n), ntry(2*n + 2), wtry(2*n + 2))

      if (is3d) then
         allocate (rhat(3, n, n))
         do i = 1, n
            do j = 1, n
               if (i /= j) then
                  x2 = 0.d0
                  do k = 1, 3
                     dx = xyz(k, j) - xyz(k, i)
                     rhat(k, j, i) = dx
                     x2 = x2 + dx*dx
                  end do
                  dx = x2**0.5d0
                  rhat(1, j, i) = rhat(1, j, i)/dx
                  rhat(2, j, i) = rhat(2, j, i)/dx
                  rhat(3, j, i) = rhat(3, j, i)/dx
               else
                  rhat(:, j, i) = 0.d0
               end if
            end do
         end do
      else
         allocate (rhat(1, 1, 1))
      end if

      ! rebuild the position index of the current path
      ipath = 0
      do i = 1, npath
         j = iabs(path(i))
         ipath(j) = sign(i, path(i))
      end do

      do ngo = 1, nsweeps
         if (nnodes + n > cap) exit
         do mgo = 1, n
            nn = int(ran3(state)*dble(n)) + 1
            do while (nn == isrc .or. nn == isnk)
               nn = int(ran3(state)*dble(n)) + 1
            end do
            ! enumerate the positions where nn may be inserted
            ncon = 0
            i1 = iabs(path(1))
            j1 = path(1)/i1
            jgo = 1
            wmax = -1.d12
            do i = 2, npath
               i2 = iabs(path(i))
               j2 = path(i)/i2
               if (j1 == j2) then
                  if (i2 /= nn) then
                     if (cost(nn, i1) > -0.1d0 .and. cost(nn, i2) > -0.1d0) then
                        ncon = ncon + 1
                        ntry(ncon) = j1*i
                        wtry(ncon) = cost(i1, i2) - cost(nn, i1) - cost(nn, i2)
                        if (wtry(ncon) > wmax) wmax = wtry(ncon)
                     end if
                     i1 = i2
                     j1 = j2
                  else
                     i2 = iabs(path(i + 1))
                     j2 = path(i + 1)/i2
                     if (j1 == j2 .and. cost(i1, i2) < -0.1d0) jgo = 0
                  end if
               else
                  if (i2 /= nn) then
                     if (cost(nn, i1) > -0.1d0) then
                        ncon = ncon + 1
                        ntry(ncon) = j1*i
                        wtry(ncon) = -cost(nn, i1)
                        if (wtry(ncon) > wmax) wmax = wtry(ncon)
                     end if
                     if (cost(nn, i2) > -0.1d0) then
                        ncon = ncon + 1
                        ntry(ncon) = j2*i
                        wtry(ncon) = -cost(nn, i2)
                        if (wtry(ncon) > wmax) wmax = wtry(ncon)
                     end if
                     i1 = i2
                     j1 = j2
                  end if
               end if
            end do
            wtot = 0.d0
            do i = 1, ncon
               wtot = wtot + dexp(scaled(wtry(i) - wmax))
               wtry(i) = wtot
            end do

            if (ipath(nn) == 0) then
               !------------------------------ insertion move
               if (ncon /= 0) then
                  prob = dexp(scaled(wmax))*wtot
                  igo = accept(prob, state)
                  if (igo == 1) then
                     try = ran3(state)*wtot
                     nswp = 1
                     do while (try > wtry(nswp))
                        nswp = nswp + 1
                        if (nswp > ncon) exit
                     end do
                     if (nswp > ncon) then
                        status = 1
                        nswp = ncon
                     end if
                     nswp = ntry(nswp)
                     mswp = nswp/iabs(nswp)
                     nswp = iabs(nswp)
                     if (is3d) then
                        i2 = iabs(path(nswp - 1))
                        j2 = path(nswp - 1)/i2
                        i3 = iabs(path(nswp))
                        j3 = path(nswp)/i3
                        if (j2 == j3) then
                           prob = bond_cos(rhat, n, i2, nn, i3)
                        else
                           prob = 1.d0
                        end if
                        if (nswp /= 2) then
                           i1 = iabs(path(nswp - 2))
                           j1 = path(nswp - 2)/i1
                           if (j1 == mswp) prob = prob*bond_cos(rhat, n, i1, i2, nn)
                           if (j1 == j3) prob = prob/bond_cos(rhat, n, i1, i2, i3)
                        end if
                        if (nswp /= npath) then
                           i4 = iabs(path(nswp + 1))
                           j4 = path(nswp + 1)/i4
                           if (mswp == j4) prob = prob*bond_cos(rhat, n, nn, i3, i4)
                           if (j2 == j4) prob = prob/bond_cos(rhat, n, i2, i3, i4)
                        end if
                        prob = dabs(prob)**bet
                        igo = accept(prob, state)
                     end if
                     if (igo == 1) then
                        do i = npath, nswp, -1
                           j = iabs(path(i))
                           k = path(i)/j
                           ipath(j) = k*(i + 1)
                           path(i + 1) = k*j
                        end do
                        path(nswp) = mswp*nn
                        ipath(nn) = mswp*nswp
                        npath = npath + 1
                     end if
                  end if
               end if
            elseif (jgo == 1) then
               !------------------------------ removal move
               prob = dexp(scaled(-wmax))/wtot
               igo = accept(prob, state)
               if (igo == 1) then
                  nswp = iabs(ipath(nn))
                  mswp = ipath(nn)/nswp
                  if (is3d) then
                     i2 = iabs(path(nswp - 1))
                     j2 = path(nswp - 1)/i2
                     i3 = iabs(path(nswp + 1))
                     j3 = path(nswp + 1)/i3
                     if (j2 == j3) then
                        prob = bond_cos(rhat, n, i2, nn, i3)
                     else
                        prob = 1.d0
                     end if
                     if (nswp /= 2) then
                        i1 = iabs(path(nswp - 2))
                        j1 = path(nswp - 2)/i1
                        if (j1 == mswp) prob = prob*bond_cos(rhat, n, i1, i2, nn)
                        if (j1 == j3) prob = prob/bond_cos(rhat, n, i1, i2, i3)
                     end if
                     if (nswp /= npath - 1) then
                        i4 = iabs(path(nswp + 2))
                        j4 = path(nswp + 2)/i4
                        if (mswp == j4) prob = prob*bond_cos(rhat, n, nn, i3, i4)
                        if (j2 == j4) prob = prob/bond_cos(rhat, n, i2, i3, i4)
                     end if
                     prob = dabs(prob)**(-bet)
                     igo = accept(prob, state)
                  end if
                  if (igo == 1) then
                     ipath(nn) = 0
                     do i = nswp + 1, npath
                        j = iabs(path(i))
                        k = path(i)/j
                        path(i - 1) = k*j
                        ipath(j) = k*(i - 1)
                     end do
                     npath = npath - 1
                  end if
               end if
            end if
         end do

         !------------------------------ path weight
         i1 = iabs(path(1))
         j1 = path(1)/i1
         wgt = 0.d0
         do i = 2, npath
            i2 = iabs(path(i))
            j2 = path(i)/i2
            if (j1 == j2) wgt = wgt + cost(i1, i2)
            i1 = i2
            j1 = j2
         end do
         wgt2 = 0.d0
         if (is3d .and. npath >= 3) then
            i1 = iabs(path(1))
            j1 = path(1)/i1
            i2 = iabs(path(2))
            j2 = path(2)/i2
            do i = 3, npath
               i3 = iabs(path(i))
               j3 = path(i)/i3
               if (j1 == j3) then
                  wgt2 = wgt2 - dlog(dabs(bond_cos(rhat, n, i1, i2, i3)))
               end if
               i1 = i2
               j1 = j2
               i2 = i3
               j2 = j3
            end do
         end if

         ndone = ndone + 1
         weights_out(1, ndone) = wgt + wgt2
         weights_out(2, ndone) = wgt
         weights_out(3, ndone) = wgt2
         lengths_out(ndone) = npath
         do i = 1, npath
            nodes_out(nnodes + i) = path(i)
         end do
         nnodes = nnodes + npath
      end do

      deallocate (ipath, ntry, wtry, rhat)

   contains

      ! The original 1D programs divide by T, the 3D programs multiply by
      ! 1/T; both are reproduced so that results are bit-for-bit identical.
      pure function scaled(x) result(y)
         real(c_double), intent(in) :: x
         real(c_double) :: y
         if (is3d) then
            y = x*bet
         else
            y = x/T
         end if
      end function scaled
   end subroutine twisp_sample

   !--------------------------------------------------------------------
   ! Metropolis acceptance exactly as in the original code: a random
   ! number is drawn only if prob < 1 (NaN and Inf are accepted).
   !--------------------------------------------------------------------
   function accept(prob, state) result(igo)
      real(c_double), intent(in) :: prob
      integer(c_int), intent(inout) :: state(57)
      integer(c_int) :: igo
      if (prob < 1.d0) then
         if (ran3(state) < prob) then
            igo = 1
         else
            igo = 0
         end if
      else
         igo = 1
      end if
   end function accept

end module allostery_twisp

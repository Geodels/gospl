! #include "petsc/finclude/petsc.h"
! #include "petsc/finclude/petscmat.h"
! #include "petscversion.h"

! #undef  CHKERRQ
! #define CHKERRQ(n) if ((n) .ne. 0) return;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                 INTERNAL MODULE                  !!
!!           Fortran Module: meshparams             !!
!!                                                  !!
!!  - Define mesh variables                         !!
!!  - Set the main functions for priority queues    !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module meshparams

  implicit none

  integer, dimension(:), allocatable :: FVnNb
  integer, dimension(:,:), allocatable :: FVnID
  double precision, dimension(:), allocatable :: FVarea
  double precision, dimension(:,:), allocatable :: FVeLgt
  double precision, dimension(:,:), allocatable :: FVvDist
  double precision, dimension(:, :), allocatable :: lcoords

  ! Queue node definition: index and elevation
  type node
    integer :: id
    double precision :: Z
  end type

  ! Simple queue node definition: index
  type snode
    integer :: id
  end type

  ! Pit queue node definition: index pit1 and pit2
  type pnode
    integer :: p1
    integer :: p2
  end type

  ! Watershed node definition: index of connected watersheds and lowest elevation
  type wnode
    integer :: id
    integer :: w1
    integer :: w2
    double precision :: Z
  end type

  ! Definition of plain queue (no priority)
  type queue
    type(node), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: pop
    procedure :: push
  end type
  type (queue) :: plainqueue

  ! Definition of simple plain queue (no priority)
  type squeue
    type(snode), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: spop
    procedure :: spush
  end type
  type (squeue) :: simplequeue

  ! Definition of pit plain queue (no priority)
  type ptqueue
    type(pnode), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: ppop
    procedure :: ppush
  end type
  type (ptqueue) :: pitqueue

  type wgraph
    type(wnode), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: wpop
    procedure :: wpush
  end type
  type (wgraph) :: graph

  ! Definition of priority queue functions
  ! The priority is based on the mesh elevation
  type pqueue
    type(node), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: PQpop
    procedure :: PQpush
    procedure :: shiftdown
  end type pqueue
  type(pqueue) :: priorityqueue

  contains

    ! Move the new element down the priority queue
    subroutine shiftdown(this, a)
      class (pqueue)  :: this
      integer :: a, parent, child
      associate (x => this%buf)
      parent = a
      do while(parent*2 <= this%n)
        child = parent*2
        if (child + 1 <= this%n) then
          if (x(child+1)%Z < x(child)%Z ) then
            child = child +1
          end if
        end if
        if (x(parent)%Z > x(child)%Z) then
          x([child, parent]) = x([parent, child])
          parent = child
        else
          exit
        end if
      end do
      end associate
    end subroutine shiftdown
    ! Pop the top element in the priority queue
    function PQpop(this) result (res)
      class(pqueue) :: this
      type(node) :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
      call this%shiftdown(1)
    end function PQpop
    ! Add a new element to the priority queue
    subroutine PQpush(this, Z, id)
      class(pqueue), intent(inout) :: this
      double precision :: Z
      integer  :: id
      type(node)  :: x
      type(node), allocatable  :: tmp(:)
      integer :: ii
      x%Z = Z
      x%id = id
      this%n = this%n +1
      if (.not.allocated(this%buf)) allocate(this%buf(1))
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if
      this%buf(this%n) = x
      ii = this%n
      do
        ii = ii / 2
        if (ii==0) exit
        call this%shiftdown(ii)
      end do
    end subroutine PQpush

    ! Pops first values in the watershed graph
    function wpop(this) result (res)
      class(wgraph) :: this
      type(wnode)   :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
    end function wpop
    ! Pushes new values in the  watershed graph
    subroutine wpush(this, w1, w2, Z, id)
      class(wgraph), intent(inout) :: this
      double precision :: Z
      integer  :: w1, w2, id, k
      type(wnode)  :: x
      type(wnode), allocatable  :: tmp(:)
      logical :: add
      if (.not.allocated(this%buf)) allocate(this%buf(1))
      add = .True.
      lp: do k = 1, this%n
        x = this%buf(k)
        if(w1 == x%w1 .and. w2 == x%w2)then
          if(Z < x%Z)then
            x%Z = Z
            x%id = id
            this%buf(k) = x
          endif
          add = .False.
          exit lp
        endif
      enddo lp
      if(add .or. this%n == 0)then
        x%Z = Z
        x%id = id
        x%w1 = w1
        x%w2 = w2
        this%n = this%n+1
        if (size(this%buf)<this%n) then
          allocate(tmp(2*size(this%buf)))
          tmp(1:this%n-1) = this%buf
          call move_alloc(tmp, this%buf)
        end if
        this%buf(this%n) = x
      endif
    end subroutine wpush

    ! Pops first values in pit plain queue
    function ppop(this) result (res)
      class(ptqueue) :: this
      type(pnode)   :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
    end function ppop
    ! Pushes new values in pit plain queue
    subroutine ppush(this, p1, p2)
      class(ptqueue), intent(inout) :: this
      integer  :: p1
      integer  :: p2
      type(pnode)  :: x
      type(pnode), allocatable  :: tmp(:)
      x%p1 = p1
      x%p2 = p2
      this%n = this%n +1
      if (.not.allocated(this%buf)) allocate(this%buf(1))
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if
      this%buf(this%n) = x
    end subroutine ppush

    ! Pops first values in a simple plain queue
    function spop(this) result (res)
      class(squeue) :: this
      type(snode)   :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
    end function spop
    ! Pushes new values in a simple plain queue
    subroutine spush(this, id)
      class(squeue), intent(inout) :: this
      integer  :: id
      type(snode)  :: x
      type(snode), allocatable  :: tmp(:)
      x%id = id
      this%n = this%n +1
      if (.not.allocated(this%buf)) allocate(this%buf(1))
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if
      this%buf(this%n) = x
    end subroutine spush

    ! Pops first values in a plain queue
    function pop(this) result (res)
      class(queue) :: this
      type(node)   :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
    end function pop
    ! Pushes new values in a plain queue
    subroutine push(this, Z, id)
      class(queue), intent(inout) :: this
      double precision :: Z
      integer  :: id
      type(node)  :: x
      type(node), allocatable  :: tmp(:)
      x%Z = Z
      x%id = id
      this%n = this%n +1
      if (.not.allocated(this%buf)) allocate(this%buf(1))
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if
      this%buf(this%n) = x
    end subroutine push

end module meshparams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!               INTERNAL FUNCTIONS                 !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euclid( p1, p2, norm)
!*****************************************************************************
! Computes the Euclidean vector norm between 2 points along the 3 dimensions
! based on fortran norm2 function.

  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(out) :: norm
  double precision :: vec(3)

  vec = p2 - p1
  norm = norm2(vec)

  return

end subroutine euclid

recursive subroutine quicksort(array, first, last, indices)
!*****************************************************************************
! quicksort routine from http://www.cgd.ucar.edu/pubsoft/TestQuicksort.html
! Reference:
! Nyhoff & Leestma, Fortran 90 for Engineers & Scientists
! (New Jersey: Prentice Hall, 1997), pp 575-577.

  double precision, dimension(:), intent(inout) :: array
  integer, intent(in)  :: first, last
  integer, dimension(:), intent(inout) :: indices

  interface
       subroutine split(array, low, high, mid, indices)
          double precision, dimension(:), intent(inout) :: array
          integer, intent(in) :: low, high
          integer, intent(out) :: mid
          integer, dimension(:), intent(inout) :: indices
       end subroutine split
  end interface

  integer :: mid

  if(first < last)then
    call split(array, first, last, mid, indices)
    call quicksort(array, first, mid-1, indices)
    call quicksort(array, mid+1, last,  indices)
  endif

end subroutine quicksort

subroutine split(array, low, high, mid, indices)
!*****************************************************************************
! used by quicksort from http://www.cgd.ucar.edu/pubsoft/TestQuicksort.html
! Reference:
! Nyhoff & Leestma, Fortran 90 for Engineers & Scientists
! (New Jersey: Prentice Hall, 1997), pp 575-577.

  double precision, dimension(:), intent(inout) :: array
  integer, intent(in) :: low, high
  integer, intent(out) :: mid
  integer, dimension(:), intent(inout) :: indices

  integer :: left, right
  double precision ::  pivot, swap
  integer :: ipivot, iswap

  left = low
  right = high
  pivot = array(low)
  ipivot = indices(low)

  do
    if( left >= right ) exit
    do
      if( left >= right .or. array(right) < pivot ) exit
      right = right - 1
    enddo
    do
      if(left > 8) exit
      if(array(left) > pivot) exit
      left = left + 1
    enddo

    if( left < right )then
      swap  = array(left)
      array(left)  = array(right)
      array(right) = swap
      iswap = indices(left)
      indices(left)  = indices(right)
      indices(right) = iswap
    endif
  enddo

  array(low) = array(right)
  array(right) = pivot
  mid = right
  indices(low) = indices(right)
  indices(right) = ipivot

end subroutine split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!         HILLSLOPE PROCESSES FUNCTIONS            !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setmaxnb(nb, maxnb)
!*****************************************************************************
! Get the maximum number of neighbours for each mesh vertice.
! This function is called once during the initialisation phase.

    use meshparams
    implicit none

    integer :: nb

    integer, intent(out) :: maxnb

    integer :: k

    maxnb  = 0
    do k = 1, nb
      if(FVarea(k)>0)then
        maxnb = max(maxnb,FVnNb(k))
      endif
    enddo

    return

end subroutine setmaxnb

subroutine sethillslopecoeff(nb, Kd, dcoeff)
!*****************************************************************************
! Define hillslope coefficients based on a finite volume spatial
! discretisation as proposed in Tucker et al. (2001).

    use meshparams
    implicit none

    integer :: nb

    double precision, intent(in) :: Kd(nb)
    double precision, intent(out) :: dcoeff(nb,9)

    integer :: k, p, n
    double precision :: s1, c, ck, cn, v

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        ck = Kd(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            cn = Kd(n)
            c = 0.5*(ck+cn)/FVarea(k)
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            s1 = s1 + v
            dcoeff(k,p+1) = -v
          endif
        enddo
        dcoeff(k,1) = 1.0 + s1
      endif
    enddo

    return

end subroutine sethillslopecoeff

subroutine fctcoeff(h, Kd, dcoeff, nb)
!*****************************************************************************
! Define the nonlinear diffusion coefficients based on a finite volume spatial
! discretisation as proposed in Tucker et al. (2001).

    use meshparams
    implicit none

    integer :: nb
    double precision, intent(in) :: h(nb)
    double precision, intent(in) :: Kd(nb)
    double precision, intent(out) :: dcoeff(nb)

    integer :: k, p, n
    double precision :: s1, c, ck, cn, v

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        ck = Kd(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            cn = Kd(n)
            c = 0.5*(ck+cn)/FVarea(k)
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            dcoeff(k) = dcoeff(k) - v*h(n)
            s1 = s1 + v
          endif
        enddo
        dcoeff(k) = dcoeff(k) + s1*h(k)
      endif
    enddo

    return

end subroutine fctcoeff

subroutine jacobiancoeff(h, Kd, Kp, dcoeff, nb)
!*****************************************************************************
! Define the nonlinear jacobian coefficients based on a finite volume spatial
! discretisation as proposed in Tucker et al. (2001).

    use meshparams
    implicit none

    integer :: nb
    double precision, intent(in) :: h(nb)
    double precision, intent(in) :: Kd(nb)
    double precision, intent(in) :: Kp(nb)
    double precision, intent(out) :: dcoeff(nb,9)

    integer :: k, p, n
    double precision :: s1, c, ck, cn, cpk, cpn, v

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        ck = Kd(k)
        cpk = Kp(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            cn = Kd(n)
            cpn = Kp(n)
            c = 0.5*(ck+cn+cpk*(h(k)-h(n)))/FVarea(k)
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            dcoeff(k,1) = dcoeff(k,1) + v
            c = 0.5*(ck+cn+cpn*(h(n)-h(k)))/FVarea(k)
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            dcoeff(k,p+1) = -v
          endif
        enddo
      endif
    enddo

    return

end subroutine jacobiancoeff

subroutine distocean(nrcv, sid, flux, rcv, wght, area, depth, dep, nb, nbi)
!*****************************************************************************
! Distribute marine sediment downstream in open water.

  use meshparams
  implicit none

  integer :: nb
  integer :: nbi
  integer :: nrcv
  integer,intent(in) :: sid(nbi)
  double precision,intent(in) :: flux(nb)
  integer,intent(in) :: rcv(nb,nrcv)
  double precision,intent(in) :: wght(nb,nrcv)
  double precision,intent(in) :: area(nb)
  double precision,intent(in) :: depth(nb)
  double precision,intent(out) :: dep(nb)

  integer :: k, i, p, n
  double precision :: flx(nb), vol(nb), nflx

  dep = 0.
  vol = area*depth
  flx = flux

  do k = 1, nbi
    i = sid(k)+1
    if(flx(i)>0.)then
      nflx = 0.
      if(flx(i)>vol(i))then
        dep(i) = depth(i)
        nflx = flx(i)-vol(i)
        vol(i) = 0.
        do p = 1, nrcv
          if(wght(i,p)>0.)then
            n = rcv(i,p)+1
            flx(n) = flx(n)+wght(i,p)*nflx
          endif
        enddo
        flx(i) = 0.
      else
        dep(i) = dep(i)+flx(i)/area(i)
        vol(i) = vol(i)-flx(i)
        flx(i) = 0.
      endif
    endif
  enddo

  return

end subroutine distocean

subroutine scale_volume(ids, vscale, scale, nb, n)
!*****************************************************************************
! Scale depression volume based on deposited one

  use meshparams
  implicit none

  integer :: nb
  integer :: n

  integer, intent(in) :: ids(nb)
  double precision, intent(in) :: vscale(n)
  double precision, intent(out) :: scale(nb)

  integer :: pitid, k
  scale = 0.

  do k = 1, nb
    pitid = ids(k)
    if(pitid>-1)then
      scale(k) = vscale(pitid+1)
    endif
  enddo

  return

end subroutine scale_volume


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!            FLOW DIRECTION FUNCTION               !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mfdreceivers( nRcv, exp, inIDs, elev, sl, rcv, dist, wgt, nb)
!*****************************************************************************
! Compute receiver characteristics based on multiple flow direction
! algorithm.

  use meshparams
  implicit none

  interface
    recursive subroutine quicksort(array, first, last, indices)
      double precision, dimension(:), intent(inout) :: array
      integer, intent(in)  :: first, last
      integer, dimension(:), intent(inout) :: indices
    end subroutine quicksort
  end interface

  integer :: nb

  integer, intent(in) :: nRcv
  double precision, intent(in) :: exp
  integer, intent(in) :: inIDs(nb)
  double precision,intent(in) :: sl
  double precision, intent(in) :: elev(nb)

  integer, intent(out) :: rcv(nb,nRcv)
  double precision, intent(out) :: dist(nb,nRcv)
  double precision, intent(out) :: wgt(nb,nRcv)

  integer :: k, n, p, kk
  double precision :: slp(8),dst(8),val,slope(8)
  integer :: id(8)

  rcv = -1
  dist = 0.
  wgt = 0.

  do k = 1, nb
    if(inIDs(k)>0)then
      if(elev(k)<=sl)then
        rcv(k,1:nRcv) = k-1
      else
        slp = 0.
        id = 0
        val = 0.
        kk = 0
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          if(n>0 .and. FVeLgt(k,p)>0.)then
            val = (elev(k) - elev(n))**exp/FVeLgt(k,p)
            if(val>0.)then
              kk = kk + 1
              slp(kk) = val
              id(kk) = n-1
              dst(kk) = FVeLgt(k,p)
            endif
          endif
        enddo

        if(kk == 0)then
          rcv(k,1:nRcv) = k-1
        elseif(kk <= nRcv)then
          val = 0.
          rcv(k,1:nRcv) = k-1
          do p = 1, kk
            rcv(k,p) = id(p)
            dist(k,p) = dst(p)
            val = val + slp(p)
          enddo
          do p = 1, nRcv
            wgt(k,p) = slp(p) / val
          enddo
        else
          rcv(k,1:nRcv) = k-1
          call quicksort(slp,1,kk,id)
          n = 0
          val = 0.
          slope = 0.
          do p = kk,kk-nRcv+1,-1
            n = n + 1
            slope(n) = slp(p)
            rcv(k,n) = id(p)
            dist(k,n) = dst(p)
            val = val + slp(p)
          enddo
          do p = 1, nRcv
            wgt(k,p) = slope(p)/val
          enddo
        endif
      endif
    endif
  enddo

  return

end subroutine mfdreceivers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!            STRATIGRAPHIC FUNCTIONS               !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine strataonesed(n, stratnb, ids, weights, strath, stratz, phis, &
                        nstrath, nstratz, nphis, nb)
!*****************************************************************************
! Record stratigraphic layers through time.

  implicit none

  integer :: nb
  integer, intent(in) :: n
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)

  double precision, intent(in) :: stratz(n,stratnb)
  double precision, intent(in) :: strath(n,stratnb)
  double precision, intent(in) :: phis(n,stratnb)

  double precision, intent(out) :: nstratz(nb,stratnb)
  double precision, intent(out) :: nstrath(nb,stratnb)
  double precision, intent(out) :: nphis(nb,stratnb)

  integer :: k, p, kk
  double precision :: tmp1, tmp2, tmp3, sum_weight

  do k = 1, nb
    sum_weight = weights(k,1) + weights(k,2) + weights(k,3)
    do kk = 1, stratnb
      tmp1 = 0.0
      tmp2 = 0.0
      tmp3 = 0.0
      do p = 1, 3
        tmp1 = tmp1 + weights(k,p)*stratz(ids(k,p)+1,kk)
        tmp2 = tmp2 + weights(k,p)*strath(ids(k,p)+1,kk)
        tmp3 = tmp3 + weights(k,p)*phis(ids(k,p)+1,kk)
      enddo
      nstratz(k,kk) = tmp1/sum_weight
      nstrath(k,kk) = tmp2/sum_weight
      nphis(k,kk) = tmp3/sum_weight
    enddo
  enddo

  return

end subroutine strataonesed

subroutine stratathreesed(n, stratnb, ids, weights, strath, stratz, stratf, &
                       stratw, phis, phif, phiw, nstrath, nstratz, nstratf, &
                       nstratw, nphis, nphif, nphiw, nb)
!*****************************************************************************
! Record stratigraphic layers through time with multiple lithologies.

  implicit none

  integer :: nb
  integer, intent(in) :: n
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)

  double precision, intent(in) :: stratz(n,stratnb)
  double precision, intent(in) :: strath(n,stratnb)
  double precision, intent(in) :: stratf(n,stratnb)
  double precision, intent(in) :: stratw(n,stratnb)
  double precision, intent(in) :: phis(n,stratnb)
  double precision, intent(in) :: phif(n,stratnb)
  double precision, intent(in) :: phiw(n,stratnb)

  double precision, intent(out) :: nstratz(nb,stratnb)
  double precision, intent(out) :: nstrath(nb,stratnb)
  double precision, intent(out) :: nstratf(nb,stratnb)
  double precision, intent(out) :: nstratw(nb,stratnb)
  double precision, intent(out) :: nphis(nb,stratnb)
  double precision, intent(out) :: nphif(nb,stratnb)
  double precision, intent(out) :: nphiw(nb,stratnb)

  integer :: k, p, kk
  double precision :: tmp1, tmp2, tmp3, tmp4, tmp5
  double precision :: tmp6, tmp7, sum_weight, tot

  do k = 1, nb
    sum_weight = weights(k,1) + weights(k,2) + weights(k,3)
    do kk = 1, stratnb
      tmp1 = 0.0
      tmp2 = 0.0
      tmp3 = 0.0
      tmp4 = 0.0
      tmp5 = 0.0
      tmp6 = 0.0
      tmp7 = 0.0
      do p = 1, 3
        tmp1 = tmp1 + weights(k,p)*stratz(ids(k,p)+1,kk)
        tmp2 = tmp2 + weights(k,p)*strath(ids(k,p)+1,kk)
        tmp3 = tmp3 + weights(k,p)*stratf(ids(k,p)+1,kk)
        tmp4 = tmp4 + weights(k,p)*phis(ids(k,p)+1,kk)
        tmp5 = tmp5 + weights(k,p)*phif(ids(k,p)+1,kk)
        tmp6 = tmp6 + weights(k,p)*stratw(ids(k,p)+1,kk)
        tmp7 = tmp7 + weights(k,p)*phiw(ids(k,p)+1,kk)
      enddo
      nstratz(k,kk) = tmp1/sum_weight
      nstrath(k,kk) = tmp2/sum_weight
      nstratf(k,kk) = tmp3/sum_weight
      nphis(k,kk) = tmp4/sum_weight
      nphif(k,kk) = tmp5/sum_weight
      nstratw(k,kk) = tmp6/sum_weight
      nphiw(k,kk) = tmp7/sum_weight
      if(nstratf(k,kk)<0.) nstratf(k,kk) = 0.0
      if(nstratw(k,kk)<0.) nstratw(k,kk) = 0.0
      tot = nstratw(k,kk)+nstratf(k,kk)
      if(tot>1.0)then
        nstratf(k,kk) = nstratf(k,kk)/tot
        nstratw(k,kk) = nstratw(k,kk)/tot
      endif
    enddo
  enddo

  return

end subroutine stratathreesed

subroutine stratafullsed(n, stratnb, ids, weights, strath, stratz, stratf, &
                           stratw, stratc, phis, phif, phiw, phic, nstrath, &
                           nstratz, nstratf, nstratw, nstratc, nphis, nphif, &
                           nphiw, nphic, nb)
!*****************************************************************************
! Record stratigraphic layers through time with carbonate module turned on.

  implicit none

  integer :: nb
  integer, intent(in) :: n
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)

  double precision, intent(in) :: stratz(n,stratnb)
  double precision, intent(in) :: strath(n,stratnb)
  double precision, intent(in) :: stratf(n,stratnb)
  double precision, intent(in) :: stratw(n,stratnb)
  double precision, intent(in) :: stratc(n,stratnb)
  double precision, intent(in) :: phis(n,stratnb)
  double precision, intent(in) :: phif(n,stratnb)
  double precision, intent(in) :: phiw(n,stratnb)
  double precision, intent(in) :: phic(n,stratnb)

  double precision, intent(out) :: nstratz(nb,stratnb)
  double precision, intent(out) :: nstrath(nb,stratnb)
  double precision, intent(out) :: nstratf(nb,stratnb)
  double precision, intent(out) :: nstratw(nb,stratnb)
  double precision, intent(out) :: nstratc(nb,stratnb)
  double precision, intent(out) :: nphis(nb,stratnb)
  double precision, intent(out) :: nphif(nb,stratnb)
  double precision, intent(out) :: nphiw(nb,stratnb)
  double precision, intent(out) :: nphic(nb,stratnb)

  integer :: k, p, kk
  double precision :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
  double precision :: tmp7, tmp8, tmp9, sum_weight, tot

  do k = 1, nb
    sum_weight = weights(k,1) + weights(k,2) + weights(k,3)
    do kk = 1, stratnb
      tmp1 = 0.0
      tmp2 = 0.0
      tmp3 = 0.0
      tmp4 = 0.0
      tmp5 = 0.0
      tmp6 = 0.0
      tmp7 = 0.0
      tmp8 = 0.0
      tmp9 = 0.0
      do p = 1, 3
        tmp1 = tmp1 + weights(k,p)*stratz(ids(k,p)+1,kk)
        tmp2 = tmp2 + weights(k,p)*strath(ids(k,p)+1,kk)
        tmp3 = tmp3 + weights(k,p)*stratf(ids(k,p)+1,kk)
        tmp8 = tmp8 + weights(k,p)*stratw(ids(k,p)+1,kk)
        tmp4 = tmp4 + weights(k,p)*stratc(ids(k,p)+1,kk)
        tmp5 = tmp5 + weights(k,p)*phis(ids(k,p)+1,kk)
        tmp6 = tmp6 + weights(k,p)*phif(ids(k,p)+1,kk)
        tmp7 = tmp7 + weights(k,p)*phic(ids(k,p)+1,kk)
        tmp9 = tmp9 + weights(k,p)*phiw(ids(k,p)+1,kk)
      enddo
      nstratz(k,kk) = tmp1/sum_weight
      nstrath(k,kk) = tmp2/sum_weight
      nstratf(k,kk) = tmp3/sum_weight
      nstratc(k,kk) = tmp4/sum_weight
      nstratw(k,kk) = tmp8/sum_weight
      if(nstratf(k,kk)<0.) nstratf(k,kk) = 0.0
      if(nstratc(k,kk)<0.) nstratc(k,kk) = 0.0
      if(nstratw(k,kk)<0.) nstratw(k,kk) = 0.0
      tot = nstratw(k,kk)+nstratf(k,kk)+nstratc(k,kk)
      if(tot>1.)then
        nstratf(k,kk) = nstratf(k,kk)/tot
        nstratw(k,kk) = nstratw(k,kk)/tot
        nstratc(k,kk) = nstratc(k,kk)/tot
      endif
      nphis(k,kk) = tmp5/sum_weight
      nphif(k,kk) = tmp6/sum_weight
      nphic(k,kk) = tmp7/sum_weight
      nphiw(k,kk) = tmp9/sum_weight
    enddo
  enddo

  return

end subroutine stratafullsed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!             PIT FILLING FUNCTIONS                !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_dir(spill, pitids, h, ptdir, m, n)
!*****************************************************************************
! Perform pit filling on each depression.

  use meshparams
  implicit none

  integer :: m
  integer :: n
  integer, intent(in) :: spill(m)
  integer, intent(in) :: pitids(n)
  double precision, intent(in) :: h(n)
  integer, intent(out) :: ptdir(n)

  integer :: i, k, c
  type (node)  :: ptID
  logical :: flag(n)

  double precision :: one

  one = 1.
  ptdir = -1

  ! Push depression spill over nodes to priority queue
  flag = .False.
  do c = 1, m
    i = spill(c)+1
    flag(i) = .True.
    ptdir(i) = 1
    call priorityqueue%PQpush(one, i)
  enddo

  ! Find closest neighbours to available priority queue node id
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    do k = 1, FVnNb(i)
      c = FVnID(i,k)+1
      if(c>0)then
        if(.not.flag(c) .and. pitids(c) > -1)then
          if(h(c)>=h(i))then
            flag(c) = .True.
            if(ptdir(c)<ptdir(i)+1)then
              ptdir(c) = ptdir(i) + 1
              call priorityqueue%PQpush(ptID%Z+one, c)
            endif
          endif
        endif
      endif
    enddo
  enddo

  return

end subroutine fill_dir

subroutine fill_rcvs(pitids, h, ptdir, rcvs, n)
!*****************************************************************************
! Perform pit filling on each depression.

  use meshparams
  implicit none

  integer :: n
  integer, intent(in) :: pitids(n)
  double precision, intent(in) :: h(n)
  integer, intent(in) :: ptdir(n)

  integer, intent(out) :: rcvs(n)

  integer :: i, k, c, v

  rcvs = -1

  ! Find receivers id for each depression
  do i = 1, n
    if(pitids(i)>-1)then
      rcvs(i)  = i-1
      v = ptdir(i)
      do k = 1, FVnNb(i)
        c = FVnID(i,k)+1
        if(c>0)then
          if(ptdir(c) < v .and. ptdir(c)>-1)then
            if(h(c)<=h(i))then
              rcvs(i)  = c-1
              v = ptdir(c)
            endif
          endif
        endif
      enddo
    endif
  enddo

  return

end subroutine fill_rcvs

subroutine nghb_dir(pitids, h, ptdir, nptdir, n)
!*****************************************************************************
! Perform pit filling on depression crossing different processors.

  use meshparams
  implicit none

  integer :: n
  integer, intent(in) :: pitids(n)
  integer,intent(in) :: ptdir(n)
  double precision, intent(in) :: h(n)

  integer, intent(out) :: nptdir(n)

  integer :: i, k, c
  type (node)  :: ptID
  logical :: flag(n)

  nptdir = -1

  ! Push existing directions to priority queue
  flag = .False.
  do c = 1, n
    if(ptdir(c)>-1)then
      flag(c) = .True.
      nptdir(c) = ptdir(c)
      call priorityqueue%PQpush(dble(ptdir(c)), c)
    endif
  enddo

  ! Find closest neighbours to available priority queue node id
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    do k = 1, FVnNb(i)
      c = FVnID(i,k)+1
      if(c>0)then
        if(.not.flag(c) .and. pitids(c) > -1)then
          if(h(c)>=h(i))then
            flag(c) = .True.
            nptdir(c) = nptdir(i)+1
            call priorityqueue%PQpush(dble(nptdir(c)), c)
          endif
        endif
      endif
    enddo
  enddo

  return

end subroutine nghb_dir

subroutine getpitvol(hlvl,elev,pit,id,vol,m,n)
!*****************************************************************************
! Compute the volume of the depressions at different elevations.

  use meshparams
  implicit none

  integer :: m
  integer :: n
  double precision,intent(in) :: hlvl(n,4)
  double precision,intent(in) :: elev(m)
  integer,intent(in) :: pit(m)
  integer,intent(in) :: id(m)
  double precision,intent(out) :: vol(n,4)


  integer :: i, p, k

  vol = 0.

  do i = 1, m
    p = pit(i)
    if(p>-1 .and. id(i)>0)then
      do k = 1, 4
        if(hlvl(p+1,k)>elev(i))then
          vol(p+1,k) = vol(p+1,k)+(hlvl(p+1,k)-elev(i))*FVarea(i)
        endif
      enddo
    endif
  enddo

  return

end subroutine getpitvol

subroutine pits_cons(pitids, pitnb, npitids, m, n)
!*****************************************************************************
! Define consecutive indices for each pits.

  use meshparams
  implicit none

  integer :: m
  integer :: n
  integer,intent(in) :: pitids(m)
  integer,intent(in) :: pitnb(n)
  integer,intent(out) :: npitids(m)

  integer :: i, k

  npitids = pitids
  do i = 1, m
    if(pitids(i)>-1)then
      lp: do k = 1, n
        if(pitids(i)==pitnb(k))then
          npitids(i) = k
          exit lp
        endif
      enddo lp
    endif
  enddo

  return

end subroutine pits_cons

subroutine edge_tile(lvl,border,elev,ledge,nb)
!*****************************************************************************
! Define edges of tile based on provided level, elevation and borders.

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: lvl
  integer,intent(in) :: border(nb, 2)
  double precision,intent(in) :: elev(nb)
  integer,intent(out) :: ledge(nb)

  integer :: c, p, nc
  double precision :: maxe
  ledge = -1

  ! Local edges
  do c = 1, nb
    if(elev(c) < lvl)then
      maxe = lvl - 1
      lp: do p = 1, FVnNb(c)
        nc = FVnID(c,p)+1
        if(nc > 0)then
          maxe = max(elev(nc),maxe)
          if(maxe >= lvl)then
            exit lp
          endif
        endif
      enddo lp
      if(maxe >= lvl .and. border(c,2) > -1)then
        ledge(c) = border(c,2)+1
      elseif(maxe >= lvl)then
        ledge(c) = 0
      else
        ledge(c) = -2
      endif
    elseif(border(c, 2) > -1)then
      ledge(c) = border(c,2)+1
    endif
  enddo

  return

end subroutine edge_tile

subroutine fill_tile(edge, elev, inids, fillz, labels, graphnb, m, nb)
!*****************************************************************************
! Perform pit filling using a priority queue approach following Barnes (2015).
! This algorithm is ran on each mesh belonging to a single processor.

  use meshparams
  implicit none

  integer :: m
  integer :: nb
  integer,intent(in) :: edge(m,2)
  double precision,intent(in) :: elev(nb)
  integer,intent(in) :: inids(nb)
  double precision,intent(out) :: fillz(nb)
  integer,intent(out) :: labels(nb)
  integer, intent(out) :: graphnb

  integer :: i, k, c, lb1, lb2, nblab
  type (node)  :: ptID

  labels = -1
  fillz = elev

  ! Push edge nodes
  do i = 1, m
    c = edge(i,1) +  1
    call priorityqueue%PQpush(fillz(c), c)
  enddo

  nblab = 1
  do while(priorityqueue%n>0 .or. plainqueue%n>0)

    if(plainqueue%n>0)then
      ptID = plainqueue%pop()
    else
      ptID = priorityqueue%PQpop()
    endif

    i = ptID%id
    if(labels(i)<=0)then
      do k = 1, FVnNb(i)
        c = FVnID(i,k)+1
        if(c>0 .and. inids(c)>0)then
          if(labels(c)>0 .and. fillz(c) <= fillz(i))then
            labels(i) = labels(c)
          endif
        endif
      enddo
      if(labels(i) <= 0)then
        labels(i) = nblab
        nblab = nblab + 1
      endif
    endif

    do k = 1, FVnNb(i)
      c = FVnID(i,k)+1
      if(c>0 .and. inids(c)>0)then
        if(labels(c)>0)then
          if(labels(c) .ne. labels(i))then
            lb1 = labels(c)
            lb2 = labels(i)
            if(labels(c)>labels(i))then
              lb1 = labels(i)
              lb2 = labels(c)
            endif
            if(fillz(c)>fillz(i))then
              call graph%wpush(lb1, lb2, fillz(c), c)
            else
              call graph%wpush(lb1, lb2, fillz(i), i)
            endif
          endif
        else
          labels(c) = labels(i)
          if(fillz(c) <= fillz(i))then
            fillz(c) = fillz(i)
            call plainqueue%push(fillz(c), c)
          else
            call priorityqueue%PQpush(fillz(c), c)
          endif
        endif
      endif
    enddo
  enddo

  ! Push border nodes
  do i = 1, m
    if(edge(i,2) == 0 .or. edge(i,2) == 2 )then
      c = edge(i,1) +  1
      call graph%wpush(labels(c), 0, elev(c), c)
    endif
  enddo

  graphnb = graph%n

  return

end subroutine fill_tile

subroutine fill_edges(nb, cgraph, maxnghbs, nelev, spillrank, spillnodes, spillid, m)
!*****************************************************************************
! This function returns filled graph based on priority flood algorithm.

  use meshparams
  implicit none

  integer :: m
  integer, intent(in) :: nb, maxnghbs
  double precision,intent(in) :: cgraph(m,5)

  integer,intent(out) :: spillrank(nb)
  integer,intent(out) :: spillnodes(nb)
  integer,intent(out) :: spillid(nb)
  double precision,intent(out) :: nelev(nb)

  integer :: k, c, nc, n1, n2, p
  logical :: inFlag(nb)
  integer :: ngbNb(nb)
  type (node)  :: ptID

  integer :: rank(nb,maxnghbs)
  integer :: ranknode(nb)
  integer :: tmp(nb)
  integer :: ngbhArr(nb,maxnghbs)
  integer :: spillnode(nb,maxnghbs)
  double precision :: spill(nb,maxnghbs)
  double precision :: spillz(nb)

  ! Initialise graph as a mesh
  inFlag = .False.
  ngbNb = 0
  ngbhArr = -1
  spillnodes = -1
  spillrank = -1
  nelev = 1.e8
  nelev(1) = -1.e8
  spillnode = -1
  spillz = 1.e8
  ranknode = -1
  spillid = -1
  rank = -1
  tmp = -1
  do k = 1, m
    n1 = int(cgraph(k,1))+1
    n2 = int(cgraph(k,2))+1
    ranknode(n1) = int(cgraph(k,5))
    ngbNb(n1) = ngbNb(n1)+1
    ngbNb(n2) = ngbNb(n2)+1
    ngbhArr(n1,ngbNb(n1)) = n2
    ngbhArr(n2,ngbNb(n2)) = n1
    spill(n1,ngbNb(n1)) = cgraph(k,3)
    spill(n2,ngbNb(n2)) = cgraph(k,3)
    spillnode(n1,ngbNb(n1)) = int(cgraph(k,4))
    spillnode(n2,ngbNb(n2)) = int(cgraph(k,4))
    rank(n1,ngbNb(n1)) = int(cgraph(k,5))
    rank(n2,ngbNb(n2)) = int(cgraph(k,5))
  enddo

  ! Perform pit filling using priority flood algorithm
  inFlag = .False.
  call priorityqueue%PQpush(nelev(1), 1)
  p = 0
  do while(priorityqueue%n > 0)
    ptID = priorityqueue%PQpop()
    c = ptID%id
    if(.not.inFlag(c))then
      p = p+1
      tmp(p) = c
      nelev(c) = ptID%Z
      inFlag(c) = .True.
      do k = 1, ngbNb(c)
        nc = ngbhArr(c,k)
        if(nc>0)then
          if(.not.inFlag(nc))then
            call priorityqueue%PQpush(max(spill(c,k),nelev(c)), nc)
            if(spillid(nc)>=0 .and. max(spill(c,k),nelev(c))<spillz(nc))then
              if(ranknode(c)==rank(c,k))then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              endif
            elseif(spillid(nc)<0)then
              if(c==1)then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              elseif(ranknode(c)==rank(c,k))then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo

  return

end subroutine fill_edges

subroutine fill_depressions(lvl, dem, fillp, wsh, ggraph, elev, m, nb)
!*****************************************************************************
! Find the fill elevation for each depressions per processors, returning the
! global mesh solution.

  use meshparams
  implicit none

  integer :: m, nb
  integer, intent(in) :: wsh(m)
  double precision, intent(in) :: lvl
  double precision, intent(in) :: dem(m)
  double precision, intent(in) :: fillp(m)
  double precision, intent(in) :: ggraph(nb)

  double precision, intent(out) :: elev(m)

  integer :: k, n

  elev = dem

  do k = 1, m
    if(dem(k)>lvl)then
      n = wsh(k)+1
      if(dem(k) < fillp(k) .and. fillp(k) > ggraph(n))then
        elev(k) = fillp(k)
      elseif(dem(k) <= ggraph(n))then
        elev(k) = ggraph(n)
      endif
    endif
  enddo

  return

end subroutine fill_depressions

subroutine combine_edges(elev, labels, ins, outs, newgraph, graphnb, m, n)
!*****************************************************************************
! Combine unstructured grids along each edges based on watershed numbers and elevations

  use meshparams
  implicit none

  integer :: m, n
  double precision, intent(in) :: elev(m)
  integer, intent(in) :: labels(m)
  integer, intent(in) :: ins(n)
  integer, intent(in) :: outs(m)

  integer,intent(out) :: graphnb
  double precision,intent(out) :: newgraph(n*8,4)

  integer :: i, c, p, nc, lb1, lb2

  type(wnode)  :: wID
  double precision :: eo

  ! Local edges
  do i = 1, n
    c = ins(i)+1
    do p = 1, FVnNb(c)
      nc = FVnID(c,p)+1
      if(nc > 0)then
        if(outs(nc) > 0)then
          if(labels(c) .ne. labels(nc))then
            eo = max(elev(c),elev(nc))
            if(labels(c)<labels(nc))then
              lb1 = labels(c)
              lb2 = labels(nc)
            else
              lb2 = labels(nc)
              lb1 = labels(c)
            endif
            if(elev(c)>elev(nc))then
              call graph%wpush(lb1, lb2, eo, c)
            else
              call graph%wpush(lb1, lb2, eo, nc)
            endif
          endif
        endif
      endif
    enddo
  enddo

  i = 1
  graphnb = graph%n
  do while(graph%n >0)
    wID = graph%wpop()
    newgraph(i,1) = wID%w1
    newgraph(i,2) = wID%w2
    newgraph(i,3) = wID%Z
    newgraph(i,4) = wID%id-1
    i = i+1
  enddo

  return

end subroutine combine_edges

subroutine graph_nodes(graphnb, newwgraph)
!*****************************************************************************
! This function returns local pit/depression graph.

  use meshparams
  implicit none

  integer,intent(in) :: graphnb
  double precision,intent(out) :: newwgraph(graphnb,4)
  type(wnode)  :: wID
  integer :: p

  p = 1
  do while(graph%n >0)
    wID = graph%wpop()
    newwgraph(p,1) = wID%w1
    newwgraph(p,2) = wID%w2
    newwgraph(p,3) = wID%Z
    newwgraph(p,4) = wID%id-1
    p = p+1
  enddo

  return

end subroutine graph_nodes

subroutine label_pits(lvl, fill, label, nb)
!*****************************************************************************
! This function finds local depression ids based on neighbors elevation.

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: lvl
  double precision,intent(in) :: fill(nb)

  integer, intent(out) :: label(nb)

  integer :: i, k, c, p, lbl
  logical :: noflow, flag(nb)
  type (snode)  :: ptID

  ! noflows = 0
  label = -1
  flag = .False.
  lbl = 0

  ! Resolve flats
  do i = 1, nb
    if(fill(i)>lvl .and. .not. flag(i))then
      flag(i) = .True.
      noflow = .True.
      lp: do k = 1, FVnNb(i)
        c = FVnID(i,k)+1
        if(c>0)then
          if(fill(i)>fill(c))then
            noflow = .False.
            exit lp
          endif
        endif
      enddo lp
      if(noflow)then
        lbl = lbl+1
        label(i) = lbl
        call simplequeue%spush(i)
        do while(simplequeue%n>0)

          ptID = simplequeue%spop()
          p = ptID%id
          do k = 1, FVnNb(p)
            c = FVnID(p,k)+1
            if(c>0)then
              if(fill(c)==fill(p) .and. .not. flag(c))then
                call simplequeue%spush(c)
                label(c) = lbl
                flag(c) = .True.
              elseif(fill(c)>fill(p))then
                flag(c) = .True.
              endif
            endif
          enddo
        enddo
      endif
    else
      flag(i) = .True.
    endif
  enddo

  return

end subroutine label_pits

subroutine spill_pts(mpirk, pitsnb, elev, pitids, border, spill, lspill, rank, m)
!*****************************************************************************
! Find for each depression its spillover point id.

  use meshparams
  implicit none

  integer :: m
  integer,intent(in) :: mpirk
  integer,intent(in) :: pitsnb
  double precision,intent(in) :: elev(m)
  integer,intent(in) :: pitids(m)
  integer,intent(in) :: border(m)

  ! Spill point ID
  integer, intent(out) :: spill(pitsnb)
  integer, intent(out) :: rank(pitsnb)
  integer,intent(out) :: lspill(m)

  integer :: i, k, c, pid

  double precision :: hmax

  spill = -1
  rank = -1
  lspill = -1

  do i = 1, m
    if(pitids(i)>-1)then
      pid = pitids(i)
      if(spill(pid+1) == -1)then
        hmax = elev(i)
        lp: do k = 1, FVnNb(i)
          c = FVnID(i,k)+1
          if(c>0)then
            if(hmax>elev(c))then
              spill(pid+1) = i-1
              rank(pid+1) = mpirk
              lspill(i) = 1
              exit lp
            elseif(hmax==elev(c))then
              if(border(c)>0)then
                spill(pid+1) = c-1
                rank(pid+1) = mpirk
                lspill(c) = 1
                exit lp
              elseif(pitids(c) == -1)then
                spill(pid+1) = c-1
                rank(pid+1) = mpirk
                lspill(c) = 1
                exit lp
              endif
            endif
          endif
        enddo lp
      endif
    endif
  enddo

  return

end subroutine spill_pts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!          MESH DECLARATION FUNCTIONS              !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine definetin( coords, cells_nodes, cells_edges, edges_nodes, &
                      area, circumcenter, ngbID, edgemax, n, nb, m)
!*****************************************************************************
! Compute for a specific triangulation the characteristics of each node and
! associated voronoi for finite volume discretizations

  use meshparams
  implicit none

  integer :: m, n, nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: cells_edges(n,3)
  integer, intent(in) :: edges_nodes(m, 2)

  double precision, intent(in) :: coords(nb,3)
  double precision, intent(in) :: area(nb)
  double precision, intent(in) :: circumcenter(3,n)

  integer, intent(out) :: ngbID(nb, 8)
  double precision, intent(out) :: edgemax

  integer :: i, n1, n2, k, l, p, eid, cid
  integer :: e, id, nbngbs, ngbIDv
  integer :: nid(2), nc(3), edge(nb, 8)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 8)

  double precision :: coords0(3), coordsID(3)
  double precision :: midpoint(3), dist


  logical :: inside

  ! Define mesh parameters
  if(allocated(FVarea)) deallocate(FVarea)
  if(allocated(FVnID)) deallocate(FVnID)
  if(allocated(FVnNb)) deallocate(FVnNb)
  if(allocated(FVeLgt)) deallocate(FVeLgt)
  if(allocated(FVvDist)) deallocate(FVvDist)
  if(allocated(lcoords)) deallocate(lcoords)

  allocate(FVarea(nb))
  allocate(FVnNb(nb))
  allocate(FVnID(nb,8))
  allocate(FVeLgt(nb,8))
  allocate(FVvDist(nb,8))
  allocate(lcoords(nb,3))

  cell_ids = -1
  edge = -1
  FVnNb = 0
  ngbID = -1
  FVeLgt = 0.
  FVvDist = 0.
  edgemax = 0.

  lcoords = coords
  FVarea = area

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells_nodes(i,1:3)+1
    do p = 1, 3
      inside = .False.
      lp: do k = 1, 8
        if( cell_ids(nc(p),k) == i-1 )then
          exit lp
        elseif( cell_ids(nc(p),k) == -1 )then
          inside = .True.
          exit lp
        endif
      enddo lp
      if( inside )then
        cell_ids(nc(p),k)  = i-1
      endif
    enddo
  enddo

  ! Find all edges connected to a given vertice
  do i = 1, m
    n1 = edges_nodes(i,1)+1
    n2 = edges_nodes(i,2)+1
    inside = .False.
    lp0: do k = 1, 8
      if(edge(n1,k) == i-1)then
        exit lp0
      elseif(edge(n1,k) == -1)then
        inside = .True.
        exit lp0
      endif
    enddo lp0
    if( inside )then
      edge(n1,k)  = i-1
      FVnNb(n1) = FVnNb(n1) + 1
    endif
    inside = .False.
    lp1: do k = 1, 8
      if(edge(n2,k) == i-1)then
        exit lp1
      elseif(edge(n2,k) == -1)then
        inside = .True.
        exit lp1
      endif
    enddo lp1
    if( inside )then
      edge(n2,k)  = i-1
      FVnNb(n2) = FVnNb(n2) + 1
    endif
  enddo

  do k = 1, nb
    ! Get triangulation edge lengths
    coords0 = coords(k,1:3)
    l = 0
    nbngbs = FVnNb(k)

    do eid = 1, nbngbs
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        ngbID(k,l) = nid(2)
        coordsID = coords(nid(2)+1,1:3)
      else
        l = l + 1
        ngbID(k,l) = nid(1)
        coordsID = coords(nid(1)+1,1:3)
      endif
      call euclid( coords0, coordsID, FVeLgt(k,l) )
      edgemax = max(edgemax, FVeLgt(k,l))
    enddo

    ! Get voronoi edge lengths
    lp2: do cid = 1, 8
      if( cell_ids(k,cid) == -1 ) exit lp2
      edgeNb(1:3) = cells_edges( cell_ids(k,cid)+1,1:3 )
      do e = 1, 3
        edges(e,1:2) = edges_nodes(edgeNb(e)+1,1:2)
        if( k-1 == edges(e,1) .or. k-1 == edges(e,2))then
          midpoint(1:3) = 0.5 * (coords(edges(e,1)+1,1:3)+coords(edges(e,2)+1,1:3))
          id = -1
          if( edges(e,1) == k-1 )then
            lp3: do i = 1, nbngbs
              ngbIDv = ngbID(k,i)
              if(ngbIDv == edges(e,2))then
                id = i
                exit lp3
              endif
            enddo lp3
          else
            lp4: do i = 1, nbngbs
              ngbIDv = ngbID(k,i)
              if(ngbIDv == edges(e,1))then
                id = i
                exit lp4
              endif
            enddo lp4
          endif
          call euclid(midpoint(1:3), circumcenter(1:3, cell_ids(k,cid)+1),  dist)
          FVvDist(k,id) = FVvDist(k,id) + dist
        endif
      enddo
    enddo lp2
  enddo


  FVnID = ngbID

end subroutine definetin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!            PREPROCESSING FUNCTIONS               !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine definegtin( nb, cells_nodes, edges_nodes, ngbNb, ngbID, n, m)
!*****************************************************************************
! Compute for the global mesh, the characteristics of each node and
! associated voronoi for finite volume discretizations

  use meshparams
  implicit none

  integer :: m, n
  integer, intent(in) :: nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: edges_nodes(m, 2)

  integer, intent(out) :: ngbID(nb, 12)
  integer, intent(out) :: ngbNb(nb)

  integer :: i, n1, n2, k, l, eid, p
  integer :: nid(2), nc(3), edge(nb, 12)
  integer :: cell_ids(nb, 12)

  logical :: inside

  cell_ids(:,:) = -1
  edge(:,:) = -1
  ngbNb(:) = 0
  ngbID(:,:) = -1

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells_nodes(i,1:3)+1
    do p = 1, 3
      inside = .False.
      lp: do k = 1, 12
        if( cell_ids(nc(p),k) == i-1 )then
          exit lp
        elseif( cell_ids(nc(p),k) == -1 )then
          inside = .True.
          exit lp
        endif
      enddo lp
      if( inside )then
        cell_ids(nc(p),k)  = i-1
      endif
    enddo
  enddo

  ! Find all edges connected to a given vertice
  do i = 1, m
    n1 = edges_nodes(i,1)+1
    n2 = edges_nodes(i,2)+1
    inside = .False.
    lp0: do k = 1, 12
      if(edge(n1,k) == i-1)then
        exit lp0
      elseif(edge(n1,k) == -1)then
        inside = .True.
        exit lp0
      endif
    enddo lp0
    if( inside )then
      edge(n1,k)  = i-1
      ngbNb(n1) = ngbNb(n1) + 1
    endif
    inside = .False.
    lp1: do k = 1, 12
      if(edge(n2,k) == i-1)then
        exit lp1
      elseif(edge(n2,k) == -1)then
        inside = .True.
        exit lp1
      endif
    enddo lp1
    if( inside )then
      edge(n2,k)  = i-1
      ngbNb(n2) = ngbNb(n2) + 1
    endif
  enddo

  do k = 1, nb
    ! Get triangulation edge lengths
    l = 0
    do eid = 1, ngbNb(k)
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        ngbID(k,l) = nid(2)
      else
        l = l + 1
        ngbID(k,l) = nid(1)
      endif
    enddo
  enddo

end subroutine definegtin

subroutine gfill(sl, elev, ngbid, hmax, fillz, nb)
!*****************************************************************************
! Perform pit filling using a priority queue approach following Barnes (2015).
! This function is done on a single processors.

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: sl
  double precision,intent(in) :: elev(nb)
  integer, intent(in) :: ngbid(nb, 12)
  double precision,intent(in) :: hmax
  double precision,intent(out) :: fillz(nb)

  logical :: flag(nb)

  integer :: i, k, c

  type (node)  :: ptID
  double precision :: h, limitz(nb)

  fillz = elev
  limitz = elev

  ! Push marine edges nodes to priority queue
  flag = .False.
  do i = 1, nb
    if(fillz(i)<sl)then
      flag(i) = .True.
      lp: do k = 1, 12
        c = ngbid(i,k)+1
        if(c>0)then
          if(fillz(c)>=sl)then
            call priorityqueue%PQpush(fillz(i), i)
            exit lp
          endif
        endif
      enddo lp
    endif
  enddo

  ! Perform pit filling using priority total queue
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    do k = 1, 12
      c = ngbid(i,k)+1
      if(c>0)then
        if(.not.flag(c))then
          flag(c) = .True.
          h = nearest(fillz(i), 1.0)
          ! Not a depression
          if(fillz(c)>h)then
            call priorityqueue%PQpush(fillz(c), c)
          ! Find a depression
          else
            fillz(c) = h
            limitz(c) = h
            if(fillz(c)-elev(c)>hmax) limitz(c) = elev(c)+hmax
            call priorityqueue%PQpush(fillz(c), c)
          endif
        endif
      endif
    enddo
  enddo

  fillz = limitz

  return

end subroutine gfill

subroutine filllabel(sl, elev, ngbid, fillz, labels, nb)
!*****************************************************************************
! Perform pit filling and watershed labeling using a variant of the priority
! queue approach following Barnes (2015).

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: sl
  double precision,intent(in) :: elev(nb)
  integer, intent(in) :: ngbid(nb, 8)
  double precision,intent(out) :: fillz(nb)
  ! Pit number and overspilling point ID
  integer,intent(out) :: labels(nb)
  logical :: flag(nb)

  integer :: i, k, c, label

  type (node)  :: ptID
  double precision :: h

  fillz = elev
  labels = -1

  ! Push marine edges nodes to priority queue
  flag = .False.
  do i = 1, nb
    if(fillz(i)<sl)then
      flag(i) = .True.
      labels(i) = 0
      lp: do k = 1, 6
        c = ngbid(i,k)+1
        if(c>0)then
          if(fillz(c)>=sl)then
            call priorityqueue%PQpush(fillz(i), i)
            exit lp
          endif
        endif
      enddo lp
    endif
  enddo

  ! Perform pit filling using priority total queue
  label = 0
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    if(labels(i)==0)then
      label = label + 1
      labels(i) = label
    endif
    do k = 1, 6
      c = ngbid(i,k)+1
      if(c>0)then
        if(.not.flag(c))then
          flag(c) = .True.
          h = nearest(fillz(i), 1.0)
          ! Not a depression
          if(fillz(c)>h)then
            call priorityqueue%PQpush(fillz(c), c)
          ! Find a depression
          else
            ! This is a new one update information
            fillz(c) = h
            call priorityqueue%PQpush(fillz(c), c)
          endif
          labels(c) = labels(i)
        endif
      endif
    enddo
  enddo

  return

end subroutine filllabel

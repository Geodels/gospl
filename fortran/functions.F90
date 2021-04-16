! #include "petsc/finclude/petsc.h"
! #include "petsc/finclude/petscmat.h"
! #include "petscversion.h"

#undef  CHKERRQ
#define CHKERRQ(n) if ((n) .ne. 0) return;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                 INTERNAL MODULE                  !!
!!           Fortran Module: meshparams             !!
!!                                                  !!
!!  - Define global variables for mesh parameters   !!
!!  - Set the main functions for priority queues    !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module meshparams

  implicit none

  integer :: nGlobal
  integer :: nLocal

  integer, dimension(:), allocatable :: FVnNb
  integer, dimension(:,:), allocatable :: gnID
  integer, dimension(:,:), allocatable :: FVnID

  double precision, dimension(:), allocatable :: FVarea
  double precision, dimension(:,:), allocatable :: FVeLgt
  double precision, dimension(:,:), allocatable :: FVvDist

  ! Queue node definition: index and elevation
  type node
    integer :: id
    double precision :: Z
  end type

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

    ! Move the new element down the stack
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

    ! Pop the top element in the stack
    function PQpop(this) result (res)
      class(pqueue) :: this
      type(node) :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
      call this%shiftdown(1)
    end function PQpop

    ! Add a new element to the stack
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

subroutine marinecoeff(nb, Ks, dcoeff)
!*****************************************************************************
! Define marine diffusion coefficients based on a finite volume spatial
! discretisation as proposed in Tucker et al. (2001).

    use meshparams
    implicit none

    integer :: nb

    double precision, intent(in) :: Ks(nb)
    double precision, intent(out) :: dcoeff(nb,9)

    integer :: k, p, n
    double precision :: s1, c, ck, cn, v

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        ck = Ks(k)
        if(ck>0.)then
          do p = 1, FVnNb(k)
            n = FVnID(k,p)+1
            cn = Ks(n)
            if(cn>0. .and. FVvDist(k,p)>0.)then
              c = 0.5*(ck+cn)/FVarea(k)
              v = c*FVvDist(k,p)/FVeLgt(k,p)
              s1 = s1 + v
              dcoeff(k,p+1) = -v
            endif
          enddo
          dcoeff(k,1) = 1.0 + s1
        else
          dcoeff(k,1) = 1.0
        endif
      endif
    enddo

    return

end subroutine marinecoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!            FLOW DIRECTION FUNCTION               !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mfdreceivers( nRcv, inIDs, elev, sl, rcv, dist, wgt, nb)
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
            val = (elev(k) - elev(n))/FVeLgt(k,p)
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


subroutine stratasimple(nb, stratnb, ids, weights, strath, stratz, phis, &
                        nstrath, nstratz, nphis)
!*****************************************************************************
! Record stratigraphic layers through time.

  implicit none

  integer, intent(in) :: nb
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)
  double precision, intent(in) :: stratz(nb,stratnb)
  double precision, intent(in) :: strath(nb,stratnb)
  double precision, intent(in) :: phis(nb,stratnb)

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

end subroutine stratasimple

subroutine stratabuild(nb, stratnb, ids, weights, strath, stratz, stratf, &
                       stratw, phis, phif, phiw, nstrath, nstratz, nstratf, &
                       nstratw, nphis, nphif, nphiw)
!*****************************************************************************
! Record stratigraphic layers through time with multiple lithologies.

  implicit none

  integer, intent(in) :: nb
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)
  double precision, intent(in) :: stratz(nb,stratnb)
  double precision, intent(in) :: strath(nb,stratnb)
  double precision, intent(in) :: stratf(nb,stratnb)
  double precision, intent(in) :: stratw(nb,stratnb)
  double precision, intent(in) :: phis(nb,stratnb)
  double precision, intent(in) :: phif(nb,stratnb)
  double precision, intent(in) :: phiw(nb,stratnb)

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

end subroutine stratabuild

subroutine stratabuildcarb(nb, stratnb, ids, weights, strath, stratz, stratf, &
                           stratw, stratc, phis, phif, phiw, phic, nstrath, &
                           nstratz, nstratf, nstratw, nstratc, nphis, nphif, &
                           nphiw, nphic)
!*****************************************************************************
! Record stratigraphic layers through time with carbonate module turned on.

  implicit none

  integer, intent(in) :: nb
  integer, intent(in) :: stratnb

  integer, intent(in) :: ids(nb,3)
  double precision,intent(in) :: weights(nb,3)
  double precision, intent(in) :: stratz(nb,stratnb)
  double precision, intent(in) :: strath(nb,stratnb)
  double precision, intent(in) :: stratf(nb,stratnb)
  double precision, intent(in) :: stratw(nb,stratnb)
  double precision, intent(in) :: stratc(nb,stratnb)
  double precision, intent(in) :: phis(nb,stratnb)
  double precision, intent(in) :: phif(nb,stratnb)
  double precision, intent(in) :: phiw(nb,stratnb)
  double precision, intent(in) :: phic(nb,stratnb)

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

end subroutine stratabuildcarb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!             PIT FILLING FUNCTIONS                !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fillpit(sl, elev, hmax, fillz, pits, nb)
!*****************************************************************************
! Perform pit filling using a priority queue approach following Barnes (2015).

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: sl
  double precision,intent(in) :: elev(nb)
  double precision,intent(in) :: hmax
  double precision,intent(out) :: fillz(nb)
  ! Pit number and overspilling point ID
  integer,intent(out) :: pits(nb,2)
  logical :: flag(nb)

  integer :: i, k, c, pitNb

  type (node)  :: ptID
  double precision :: h, limitz(nb)

  fillz = elev
  limitz = elev
  pits = -1

  ! Push marine edges nodes to priority queue
  flag = .False.
  do i = 1, nb
    if(fillz(i)<sl)then
      flag(i) = .True.
      lp: do k = 1, 6
        c = gnID(i,k)
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
  pitNb = 0
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    do k = 1, 6
      c = gnID(i,k)
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
            if(pits(i,1)==-1)then
              fillz(c) = h
              limitz(c) = h
              if(fillz(c)-elev(c)>hmax) limitz(c) = elev(c)+hmax
              pitNb = pitNb+1
              pits(i,1) = pitNb
              pits(c,1) = pitNb
              pits(i,2) = i-1
              pits(c,2) = i-1
            ! This is an existing one: add nodes to the depression
            else
              fillz(c) = h
              limitz(c) = h
              if(fillz(c)-elev(c)>hmax) limitz(c) = elev(c)+hmax
              pits(c,1) = pits(i,1)
              pits(c,2) = pits(i,2)
            endif
            call priorityqueue%PQpush(fillz(c), c)
          endif
        endif
      endif
    enddo
  enddo

  fillz = limitz

  return

end subroutine fillpit

subroutine filllabel(sl, elev, fillz, labels, nb)
!*****************************************************************************
! Perform pit filling and watershed labeling using a variant of the priority
! queue approach following Barnes (2015).

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: sl
  double precision,intent(in) :: elev(nb)
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
        c = gnID(i,k)
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
      c = gnID(i,k)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!          MESH DECLARATION FUNCTIONS              !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ngbglob( nb, ngbIDs )
!*****************************************************************************
! Set global mesh neighbours indices.

  use meshparams
  implicit none

  integer, intent(in) :: nb
  integer, intent(in) :: ngbIDs(nb, 8)

  nGlobal = nb
  if(allocated(gnID)) deallocate(gnID)
  allocate(gnID(nGlobal,8))
  gnID = ngbIDs+1

  return

end subroutine ngbglob

subroutine definetin( nloc, coords, lgIDs, cells_nodes, cells_edges, edges_nodes, &
                      area, circumcenter, ngbID, edgemax, n, nb, m)
!*****************************************************************************
! Compute for a specific triangulation the characteristics of each node and
! associated voronoi for finite volume discretizations

  use meshparams
  implicit none

  integer :: m, n, nb
  integer, intent(in) :: nloc
  integer, intent(in) :: lgIDs(nb)
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: cells_edges(n,3)
  integer, intent(in) :: edges_nodes(m, 2)

  double precision, intent(in) :: coords(nb,3)
  double precision, intent(in) :: area(nb)
  double precision, intent(in) :: circumcenter(3,n)

  integer, intent(out) :: ngbID(nloc, 8)
  double precision, intent(out) :: edgemax

  integer :: i, n1, n2, k, l, p, eid, cid
  integer :: e, id, lid, lid2, nbngbs, ngbIDv
  integer :: nid(2), nc(3), edge(nb, 8)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 8)

  double precision :: coords0(3), coordsID(3)
  double precision :: midpoint(3), dist

  integer, dimension(:), allocatable :: gFVnNb
  integer, dimension(:,:), allocatable :: gngbID
  ! integer, dimension(:,:), allocatable :: FVnID
  double precision, dimension(:,:), allocatable :: gFVeLgt
  double precision, dimension(:,:), allocatable :: gFVvDist

  logical :: inside

  ! Define fortran global parameters
  nLocal = nloc

  if(nLocal < nb)then
    if(allocated(gFVnNb)) deallocate(gFVnNb)
    if(allocated(gngbID)) deallocate(gngbID)
    if(allocated(gFVeLgt)) deallocate(gFVeLgt)
    if(allocated(gFVvDist)) deallocate(gFVvDist)

    allocate(gFVnNb(nb))
    allocate(gngbID(nb,8))
    allocate(gFVeLgt(nb,8))
    allocate(gFVvDist(nb,8))

    gFVnNb = 0
    gngbID = -1
    gFVeLgt = 0.
    gFVvDist = 0.
  endif

  if(allocated(FVarea)) deallocate(FVarea)
  if(allocated(FVnID)) deallocate(FVnID)
  if(allocated(FVnNb)) deallocate(FVnNb)
  if(allocated(FVeLgt)) deallocate(FVeLgt)
  if(allocated(FVvDist)) deallocate(FVvDist)

  allocate(FVarea(nLocal))
  allocate(FVnNb(nLocal))
  allocate(FVnID(nLocal,8))
  allocate(FVeLgt(nLocal,8))
  allocate(FVvDist(nLocal,8))

  cell_ids = -1
  edge = -1
  FVnNb = 0
  ngbID = -1
  FVeLgt = 0.
  FVvDist = 0.
  edgemax = 0.

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
      if(nLocal < nb)then
        gFVnNb(n1) = gFVnNb(n1) + 1
      else
        FVnNb(n1) = FVnNb(n1) + 1
      endif
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
      if(nLocal < nb)then
        gFVnNb(n2) = gFVnNb(n2) + 1
      else
        FVnNb(n2) = FVnNb(n2) + 1
      endif
    endif
  enddo

  do k = 1, nb

    ! Get triangulation edge lengths
    coords0 = coords(k,1:3)
    l = 0

    if(nLocal < nb)then
      nbngbs = gFVnNb(k)
    else
      nbngbs = FVnNb(k)
    endif

    do eid = 1, nbngbs
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        if(nLocal < nb)then
          gngbID(k,l) = nid(2)
        else
          ngbID(k,l) = nid(2)
        endif
        coordsID = coords(nid(2)+1,1:3)
      else
        l = l + 1
        if(nLocal < nb)then
          gngbID(k,l) = nid(1)
        else
          ngbID(k,l) = nid(1)
        endif
        coordsID = coords(nid(1)+1,1:3)
      endif
      if(nLocal < nb)then
        call euclid( coords0, coordsID, gFVeLgt(k,l) )
        edgemax = max(edgemax, gFVeLgt(k,l))
      else
        call euclid( coords0, coordsID, FVeLgt(k,l) )
        edgemax = max(edgemax, FVeLgt(k,l))
      endif
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
              if(nLocal < nb)then
                ngbIDv = gngbID(k,i)
              else
                ngbIDv = ngbID(k,i)
              endif
              if(ngbIDv == edges(e,2))then
                id = i
                exit lp3
              endif
            enddo lp3
          else
            lp4: do i = 1, nbngbs
              if(nLocal < nb)then
                ngbIDv = gngbID(k,i)
              else
                ngbIDv = ngbID(k,i)
              endif
              if(ngbIDv == edges(e,1))then
                id = i
                exit lp4
              endif
            enddo lp4
          endif
          call euclid(midpoint(1:3), circumcenter(1:3, cell_ids(k,cid)+1),  dist)
          if(nLocal < nb)then
            gFVvDist(k,id) = gFVvDist(k,id) + dist
          else
            FVvDist(k,id) = FVvDist(k,id) + dist
          endif
        endif
      enddo
    enddo lp2

    ! In case of parallel simulation only allocate local arrays
    if(nLocal < nb)then
      ! Is the vertex in the partition
      if(lgIDs(k)>-1)then
        lid = lgIDs(k) + 1
        FVnNb(lid) = gFVnNb(k)
        FVarea(lid) = area(k)
        p = 1
        do i = 1, gFVnNb(k)
          ngbIDv = gngbID(k, i)
          lid2 = lgIDs(ngbIDv + 1) + 1
          if(lid2 > 0)then
            ngbID(lid, p) =  lid2 - 1
            FVeLgt(lid, p) = gFVeLgt(k, i)
            FVvDist(lid, p) = gFVvDist(k, i)
            FVnID(lid, p) = ngbID(lid, p)
            p = p+1
          else
            FVnNb(lid) = FVnNb(lid) - 1
          endif
        enddo
      endif
    endif
  enddo


  if(nLocal < nb)then
    deallocate(gFVnNb)
    deallocate(gngbID)
    deallocate(gFVeLgt)
    deallocate(gFVvDist)
  else
    FVnID = ngbID
    FVarea = area
  endif

end subroutine definetin

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

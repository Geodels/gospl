#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscmat.h"

#include "petscversion.h"

#undef  CHKERRQ
#define CHKERRQ(n) if ((n) .ne. 0) return;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INTERNAL FUNCTIONS                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module meshparams

  implicit none

  integer :: nGlobal
  integer :: nLocal
  real(kind=8) :: eps

  integer, dimension(:,:), allocatable :: gnID
  integer, dimension(:,:), allocatable :: FVnID
  integer, dimension(:), allocatable :: FVnNb

  real( kind=8 ), dimension(:), allocatable :: FVarea
  real( kind=8 ), dimension(:,:), allocatable :: FVeLgt
  real( kind=8 ), dimension(:,:), allocatable :: FVvDist

  ! Queue node definition: index and elevation
  type node
    integer :: id
    real(kind=8) :: Z
  end type

  ! Definition of priority queue (1 priority: elevation)
  type pqueue
    type(node), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: PQpop
    procedure :: PQpush
    procedure :: shiftdown
  end type

  type(pqueue) :: priorityqueue

  contains

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
    function PQpop(this) result (res)
      class(pqueue) :: this
      type(node)   :: res
      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1
      call this%shiftdown(1)
    end function PQpop
    subroutine PQpush(this, Z, id)
      class(pqueue), intent(inout) :: this
      real(kind=8) :: Z
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

subroutine euclid( p1, p2, norm)
!*****************************************************************************
! Computes the Euclidean vector norm between 2 points
! along the 3 dimensions

  implicit none

  real( kind=8 ), intent(in) :: p1(3)
  real( kind=8 ), intent(in) :: p2(3)
  real( kind=8 ), intent(out) :: norm
  real( kind=8 ) :: vec(3)

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

  real( kind=8 ), dimension(:), intent(inout) :: array
  integer, intent(in)  :: first, last
  integer, dimension(:), intent(inout) :: indices

  interface
       subroutine split(array, low, high, mid, indices)
          real( kind=8 ), dimension(:), intent(inout) :: array
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
! used by quicksort  from http://www.cgd.ucar.edu/pubsoft/TestQuicksort.html
! Reference:
! Nyhoff & Leestma, Fortran 90 for Engineers & Scientists
! (New Jersey: Prentice Hall, 1997), pp 575-577.

  real( kind=8 ), dimension(:), intent(inout) :: array
  integer, intent(in) :: low, high
  integer, intent(out) :: mid
  integer, dimension(:), intent(inout) :: indices

  integer :: left, right
  real( kind=8 ) ::  pivot, swap
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HILLSLOPE PROCESSES FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setHillslopeCoeff(nb, Kd, dcoeff, maxnb)
!*****************************************************************************
! Define hillslope coefficients

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd
    integer, intent(out) :: maxnb
    real( kind=8 ), intent(out) :: dcoeff(nb,8)

    integer :: k, p
    real( kind=8 ) :: s1, c, v

    dcoeff = 0.
    maxnb  = 0
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        c = Kd/FVarea(k)
        maxnb = max(maxnb,FVnNb(k))
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            s1 = s1 + v
            dcoeff(k,p+1) = -v
          endif
        enddo
        dcoeff(k,1) = 1.0 + s1
      endif
    enddo

    return

end subroutine setHillslopeCoeff

subroutine setDiffusionCoeff(Kd, limit, elev, dh, dcoeff, nb)
!*****************************************************************************
! Define freshly deposited sediments diffusion implicit matrix coefficients

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd
    real( kind=8 ), intent(in) :: limit
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: dh(nb)

    real( kind=8 ), intent(out) :: dcoeff(nb,13)

    integer :: k, p, n
    real( kind=8 ) :: s1, c, v, limiter

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        c = Kd/FVarea(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            n = FVnID(k,p)+1
            limiter = 0.
            if(elev(n)>elev(k))then
              limiter = dh(n)/(dh(n)+limit)
            elseif(elev(n)<elev(k))then
              limiter = dh(k)/(dh(k)+limit)
            endif
            s1 = s1 + v*limiter
            dcoeff(k,p+1) = -v*limiter
          endif
        enddo
        dcoeff(k,1) = 1.0 + s1
      endif
    enddo

    return

end subroutine setDiffusionCoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FLOW DIRECTION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MFDreceivers( nRcv, inIDs, elev, rcv, slope, dist, wgt, nb)
!*****************************************************************************
! Compute receiver characteristics based on multiple flow direction algorithm

  use meshparams
  implicit none

  interface
    recursive subroutine quicksort(array, first, last, indices)
      real( kind=8 ), dimension(:), intent(inout) :: array
      integer, intent(in)  :: first, last
      integer, dimension(:), intent(inout) :: indices
    end subroutine quicksort
  end interface

  integer :: nb

  integer, intent(in) :: nRcv
  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: elev(nb)

  integer, intent(out) :: rcv(nb,nRcv)
  real( kind=8 ), intent(out) :: slope(nb,nRcv)
  real( kind=8 ), intent(out) :: dist(nb,nRcv)
  real( kind=8 ), intent(out) :: wgt(nb,nRcv)

  integer :: k, n, p, kk
  real( kind=8 ) :: slp(8),dst(8),val
  integer :: id(8)

  rcv = -1
  slope = 0.
  dist = 0.
  wgt = 0.

  do k = 1, nb
    if(inIDs(k)>0)then
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
          slope(k,p) = slp(p)
          dist(k,p) = dst(p)
          val = val + slp(p)
        enddo
        do p = 1, nRcv
          wgt(k,p) = slp(p)/val
        enddo
      else
        rcv(k,1:nRcv) = k-1
        call quicksort(slp,1,kk,id)
        n = 0
        val = 0.
        do p = kk,kk-nRcv+1,-1
          n = n + 1
          rcv(k,n) = id(p)
          slope(k,n) = slp(p)
          dist(k,n) = dst(p)
          val = val + slp(p)
        enddo
        do p = 1, nRcv
          wgt(k,p) = slope(k,p)/val
        enddo
      endif
    endif
  enddo

  return

end subroutine MFDreceivers

subroutine fillPIT(sl, elev, fillz, nb)
  use meshparams
  implicit none

  integer :: nb
  real(kind=8),intent(in) :: sl
  real(kind=8),intent(in) :: elev(nb)
  real(kind=8),intent(out) :: fillz(nb)
  logical :: flag(nb)

  integer :: i, k, c

  type (node)  :: ptID

  fillz = elev

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

  ! Perform pit filling using priority flood algorithm variant from Barnes 2014
  ! Here we use only one priority total queue as the plain queue doesn't ensure
  ! a consistent pit+epsilon filling in our case... not sure why?
  do while(priorityqueue%n>0)
    ptID = priorityqueue%PQpop()
    i = ptID%id
    do k = 1, 6
      c = gnID(i,k)
      if(c>0)then
        if(.not.flag(c))then
          flag(c) = .True.
          fillz(c) = max(fillz(c),fillz(i)+eps)
          call priorityqueue%PQpush(fillz(c), c)
        endif
      endif
    enddo
  enddo

  return

end subroutine fillPIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MESH DECLARATION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ngbGlob( nb, ngbIDs )

  use meshparams
  implicit none

  integer, intent(in) :: nb
  integer, intent(in) :: ngbIDs(nb, 6)

  ! integer :: k, p

  nGlobal = nb
  if(allocated(gnID)) deallocate(gnID)
  allocate(gnID(nGlobal,6))
  gnID = ngbIDs+1
  eps = 0.0000001

  return

end subroutine ngbGlob

subroutine defineTIN( coords, cells_nodes, cells_edges, edges_nodes, area, &
                      circumcenter, ngbID, n, nb, m)
!*****************************************************************************
! Compute for a specific triangulation the characteristics of each node and
! associated voronoi for finite volume discretizations

  use meshparams
  implicit none

  integer :: m, n, nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: cells_edges(n,3)
  integer, intent(in) :: edges_nodes(m, 2)

  real( kind=8 ), intent(in) :: coords(nb,3)
  real( kind=8 ), intent(in) :: area(nb)
  real( kind=8 ), intent(in) :: circumcenter(3,n)

  integer, intent(out) :: ngbID(nb, 8)

  integer :: i, n1, n2, k, l, p, eid, cid, e, id
  integer :: nid(2), nc(3), edge(nb, 8)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 8)

  real( kind=8 ) :: coords0(3), coordsID(3)
  real( kind=8 ) :: midpoint(3), dist

  logical :: inside


  ! Define fortran global parameters
  nLocal = nb
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
  FVarea = area
  cell_ids = -1
  edge = -1
  FVnNb = 0
  ngbID = -1
  FVeLgt = 0.
  FVvDist = 0.

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
    do eid = 1, FVnNb(k)
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
            lp3: do i = 1, FVnNb(k)
              if(ngbID(k,i) == edges(e,2))then
                id = i
                exit lp3
              endif
            enddo lp3
          else
            lp4: do i = 1, FVnNb(k)
              if(ngbID(k,i) == edges(e,1))then
                id = i
                exit lp4
              endif
            enddo lp4
          endif
          call euclid( midpoint(1:3), circumcenter(1:3,cell_ids(k,cid)+1),  dist)
          FVvDist(k,id) = FVvDist(k,id) + dist
        endif
      enddo
    enddo lp2

  enddo

  FVnID = ngbID

end subroutine defineTIN

subroutine defineGTIN( nb, cells_nodes, edges_nodes, ngbNb, ngbID, n, m)
!*****************************************************************************
! Compute for global triangulation the characteristics of each node

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

end subroutine defineGTIN

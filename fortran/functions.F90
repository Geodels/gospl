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

  integer, dimension(:), allocatable :: FVnNb           ! Number of vertex neighbors
  integer, dimension(:), allocatable :: FVgnNb           ! Number of global vertex neighbors
  integer, dimension(:,:), allocatable :: FVnID          ! Index of vertex neighbors
  integer, dimension(:,:), allocatable :: FVgnID          ! Index of global vertex neighbors
  integer, dimension(:,:), allocatable :: FVnIDfNb         ! Index of vertex neighbors connected face nb  
  double precision, dimension(:), allocatable :: FVarea    ! Voronoi area
  double precision, dimension(:,:), allocatable :: FVeLgt   ! Length btw vertex
  double precision, dimension(:,:), allocatable :: FVvDist   ! Voronoi edge length
  double precision, dimension(:, :), allocatable :: lcoords     ! Vertex coordinates
  double precision, dimension(:, :, :), allocatable :: faceVec  ! Voronoi edge normal vector
  double precision, dimension(:, :, :), allocatable :: midFace  ! Voronoi edge mid face coordinates
  double precision, dimension(:, :), allocatable :: faceVel  ! Voronoi edge mid face velocity

  integer, dimension(:), allocatable :: stencilNb        
  integer, dimension(:,:), allocatable :: stencilNgb      

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

  ! Definition of watershed queue (no priority)
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

subroutine spherical2lonlat(xyz, long, lat) 

  double precision, intent(in) :: xyz(3)
  double precision, intent(out) :: long
  double precision, intent(out) :: lat

  double precision :: r
  double precision :: pi = 3.1415926535897932_8

  r = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))
  
  lat = asin(xyz(3)/r)*(180.0_8/pi)
  long = atan2(xyz(2),xyz(1))*(180.0_8/pi)

  lat = lat * pi / 180.0_8
  long = long * pi / 180.0_8

  return

end subroutine spherical2lonlat

subroutine great_circle_distance(long1, lat1, long2, lat2, radius, d) 
!*****************************************************************************
! Computes the great circle distance on a spherical body, using the 
! Vincenty algorithm.

  implicit none

  double precision, intent(in) :: long1    !! longitude of first site [rad]
  double precision, intent(in) :: lat1     !! latitude of the first site [rad]
  double precision, intent(in) :: long2    !! longitude of the second site [rad]
  double precision, intent(in) :: lat2     !! latitude of the second site [rad]
  double precision, intent(in) :: radius   !! radius sphere in meters
  double precision, intent(out) :: d       !! great circle distance from 1 to 2 [km]

  double precision :: c1,s1,c2,s2,dlon,clon,slon

  c1   = cos(lat1)
  s1   = sin(lat1)
  c2   = cos(lat2)
  s2   = sin(lat2)
  dlon = long1-long2
  clon = cos(dlon)
  slon = sin(dlon)
  
  d = radius*atan2(sqrt((c2*slon)**2+(c1*s2-s1*c2*clon)**2), &
      (s1*s2+c1*c2*clon))

  return

end subroutine great_circle_distance

subroutine distance(p1, p2, norm)
!*****************************************************************************
! Computes the distance between 2 points along the 2 or 3 dimensions.
  
  use meshparams
  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(out) :: norm

  double precision :: vec(3), lon1, lon2, lat1, lat2, radius

  if(p2(3) == 0. .and. p1(3) == 0)then
    vec = p2 - p1
    norm = norm2(vec)
  else
    radius = norm2(p1)
    call spherical2lonlat(p1, lon1, lat1) 
    call spherical2lonlat(p2, lon2, lat2)
    call great_circle_distance(lon1, lat1, lon2, lat2, radius, norm)
  endif 

  return

end subroutine distance

subroutine projecttangentplane(p, v, pv)
!*****************************************************************************
! Project a vector v at point P onto the tangent plane of a sphere at P.
! Parameters:
!   p (array-like): The point on the sphere (x, y, z).
!   v (array-like): The vector to be projected (vx, vy, vz).
! Returns:
!   pv array: The projected vector on the tangent plane.

  implicit none

  double precision, intent(in) :: p(3)
  double precision, intent(in) :: v(3)
  double precision, intent(out) :: pv(3)

  double precision :: n(3), v_dot_n

  ! Calculate the normal vector n at point p (normalized P)
  n = p/norm2(p)
    
  ! Calculate the dot product of v and n
  v_dot_n = dot_product(v, n)
    
  ! Calculate the projection of v onto the tangent plane
  pv = v - v_dot_n * n
    
  return 

end subroutine projecttangentplane
  
subroutine slerpmidpt(p1, p2, mpt)
!*****************************************************************************
! Computes the mid point of an arc based on spherical linear interpolation 
! between 2 points along the 3 dimensions.
! see https://en.wikipedia.org/wiki/Slerp#Quaternion_slerp and 
! https://brsr.github.io/2021/05/01/vector-spherical-geometry.html

  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(out) :: mpt(3)

  double precision :: v1(3), v2(3), pt(3), a

  ! Set on the unit sphere
  v1 = p1/norm2(p1)
  v2 = p2/norm2(p2)

  ! Angle subtended by the arc
  a = acos(dot_product(v1, v2))

  pt = sin(0.5*a)*v1/sin(a) + sin(0.5*a)*v2/sin(a)
  mpt = norm2(p1)*pt/norm2(pt)

  return

end subroutine slerpmidpt

subroutine parallel_transport(pt1, pt2, v1, v2)
!*****************************************************************************
! Performs parallel transport along the unique minimizing geodesic connecting p1 and p2, 
! if it exists, on the sphere.    
! pt1: A vector (or column matrix) representing a point on the manifold.
! pt2: A vector (or column matrix) representing a point on the manifold.
! v1: A vector (or column matrix) tangent to p1. 
! v2: A vector tangent to p2.
! Follows: https://rdrr.io/cran/GeodRegr/src/R/others.R

  implicit none

  double precision, intent(in) :: pt1(3)
  double precision, intent(in) :: pt2(3)
  double precision, intent(in) :: v1(3)
  double precision, intent(out) :: v2(3)

  double precision :: p1(3), p2(3), w(3), invar(3), pv1(3), theta
  double precision :: e1(3), e2(3), a, t, magp1, magp2, tang(3)

  call projecttangentplane(pt1, v1, pv1)

  ! Set on the unit sphere
  magp1 = norm2(pt1)
  magp2 = norm2(pt2)
  p1 = pt1/magp1
  p2 = pt2/magp2

  ! Log map
  a = max(min(dot_product(p1, p2), 1.0), -1.0)
  theta = acos(a)
  tang = p2 - a * p1
  t = norm2(tang)
  if(t == 0)then
    w(1:3) = 0
  else
    w = theta * (tang / t)
  endif

  t = norm2(w)
  if(t == 0)then
    v2 = pv1
  else
    e1 = p1
    e2 = w / t
    a = dot_product(pv1, e2)
    invar = pv1 - a * e2
    v2 = a * (cos(t) * e2 - sin(t) * e1) + invar
  endif

  return

end subroutine parallel_transport

subroutine slerpvec(p1, p2, t, vec)
!*****************************************************************************
! Computes the face vector direction based on spherical linear interpolation 
! between 2 points along the 3 dimensions.
! see https://en.wikipedia.org/wiki/Slerp#Quaternion_slerp and 
! https://brsr.github.io/2021/05/01/vector-spherical-geometry.html

  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(in) :: t
  double precision, intent(out) :: vec(3)

  double precision :: v1(3), v2(3), newpt(3), a

  ! Set on the unit sphere
  v1 = p1/norm2(p1)
  v2 = p2/norm2(p2)

  ! Angle subtended by the arc
  a = acos(dot_product(v1, v2))

  newpt = sin((1.-t)*a)*v1/sin(a) + sin(t*a)*v2/sin(a)
  vec = (v2-newpt)/norm2(v2-newpt)

  return

end subroutine slerpvec

subroutine striarea( v1, v2, v3, area )
!*****************************************************************************
! Spherical triangle areas computation on the unit sphere.

  implicit none

  double precision, intent(in) :: v1(3)
  double precision, intent(in) :: v2(3)
  double precision, intent(in) :: v3(3)
  double precision, intent(out) :: area


  double precision :: lat1, lon1, lat2, lon2, lat3, lon3
  double precision :: a, b, c
  double precision :: alpha, beta, gamma, E
  double precision :: pi = 3.1415926535897932_8
  double precision :: r

  r = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))

  call spherical2lonlat(v1, lon1, lat1) 
  call spherical2lonlat(v2, lon2, lat2) 
  call spherical2lonlat(v3, lon3, lat3) 
  
  ! Calculate the sides of the spherical triangle (arc lengths)
  a = acos(sin(lat2)*sin(lat3) + cos(lat2)*cos(lat3)*cos(lon3-lon2))
  b = acos(sin(lat1)*sin(lat3) + cos(lat1)*cos(lat3)*cos(lon3-lon1))
  c = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))

  ! Calculate the angles of the spherical triangle
  alpha = acos((cos(a) - cos(b)*cos(c)) / (sin(b)*sin(c)))
  beta  = acos((cos(b) - cos(a)*cos(c)) / (sin(a)*sin(c)))
  gamma = acos((cos(c) - cos(a)*cos(b)) / (sin(a)*sin(b)))

  ! Calculate the spherical excess
  E = alpha + beta + gamma - pi

  ! Calculate the area of the spherical triangle
  area = E * r**2

end subroutine striarea

subroutine gtriarea(p1, p2, p3, area)
!*****************************************************************************
! Heron's Formula for the area of the triangle

  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(in) :: p3(3)
  double precision, intent(out) :: area

  double precision :: x1, x2, x3, y1, y2, y3
  double precision :: l1, l2, l3, p

  x1 = p1(1)
  x2 = p2(1)
  x3 = p3(1)
  y1 = p1(2)
  y2 = p2(2)
  y3 = p3(2)
  
  if(x1 == x2 .and. y1 == y2)then
    area = 0.0
  elseif(x1 == x3 .and. y1 == y3)then
    area = 0.0
  elseif(x2 == x3 .and. y2 == y3)then
    area = 0.0
  else
    l1 = sqrt((x1 - x2)**2 + (y1 - y2)**2)
    l2 = sqrt((x2 - x3)**2 + (y2 - y3)**2)
    l3 = sqrt((x3 - x1)**2 + (y3 - y1)**2)
    p = (l1 + l2 + l3)/2
    area = sqrt(p * (p - l1) * (p - l2) * (p - l3))
  endif

  return

end subroutine gtriarea

subroutine gtriareasphere(p1, p2, p3, area)
!*****************************************************************************
! Area of a triangle on the sphere

  implicit none

  double precision, intent(in) :: p1(3)
  double precision, intent(in) :: p2(3)
  double precision, intent(in) :: p3(3)
  double precision, intent(out) :: area

  double precision :: x1, x2, x3, y1, y2, y3, z1, z2, z3
  double precision :: v1x, v1y, v1z, v2x, v2y, v2z
  double precision :: cross_x, cross_y, cross_z

  x1 = p1(1)
  x2 = p2(1)
  x3 = p3(1)
  y1 = p1(2)
  y2 = p2(2)
  y3 = p3(2)
  z1 = p1(3)
  z2 = p2(3)
  z3 = p3(3)
  
  ! Compute the vectors v1 = P2 - P1 and v2 = P3 - P1
  v1x = x2 - x1
  v1y = y2 - y1
  v1z = z2 - z1

  v2x = x3 - x1
  v2y = y3 - y1
  v2z = z3 - z1

  ! Compute the cross product v1 x v2
  cross_x = v1y * v2z - v1z * v2y
  cross_y = v1z * v2x - v1x * v2z
  cross_z = v1x * v2y - v1y * v2x

  ! Compute the area of the triangle
  area = 0.5 * sqrt(cross_x**2 + cross_y**2 + cross_z**2)

  return

end subroutine gtriareasphere

subroutine euclid(p1, p2, norm)
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
      if(left > 12) exit
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
!!   HILLSLOPE AND ADVECTION PROCESSES FUNCTIONS    !!
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
    double precision, intent(out) :: dcoeff(nb,13)

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

subroutine getfacevelocity(nb, vel)
!*****************************************************************************
! Update the face velocity of each vertex-based voronoi and compute the dot product

  use meshparams
  implicit none

  integer :: nb
  double precision, intent(in) :: vel(nb, 3)

  integer :: k, p, n
  double precision :: v(3), v1(3), v2(3), radius

  if(.not. allocated(faceVel)) allocate(faceVel(nb,12))
  faceVel(:,:) = 0.0
  radius = 0.
  if(lcoords(1,3) .ne. 0) radius = norm2(lcoords(1,1:3))

  if(radius == 0.0)then
    do k = 1, nb
      if(FVarea(k)>0)then
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            ! Velocity at the face is taken to be the linear interpolation for each vertex (in a vertex-centered discretisation the dual of the delaunay triangulation i.e. the voronoi mesh has its edges on the middle of the nodes edges)
            v = 0.5 * (vel(k, 1:3) + vel(n, 1:3))
            faceVel(k,p) = dot_product(v(1:3), faceVec(k,p,1:3))
          endif
        enddo
      endif
    enddo
  else
    do k = 1, nb
      if(FVarea(k)>0)then
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            ! Velocity at the face is based on parallel transport on the sphere from the vertex center to the middle of the voronoi edge.
            call parallel_transport(lcoords(k, 1:3), midFace(k,p,1:3), vel(k, 1:3), v1(1:3))
            n = FVnID(k,p)+1
            call parallel_transport(lcoords(n, 1:3), midFace(k,p,1:3), vel(n, 1:3), v2(1:3))
            v = 0.5 * (v1(1:3) + v2(1:3))
            faceVel(k,p) = dot_product(v(1:3), faceVec(k,p,1:3))
          endif
        enddo
      endif
    enddo

  endif

  return

end subroutine getfacevelocity

subroutine advecupwind(nb, dt, lcoeff)
!*****************************************************************************
! Define advection coefficients based on a finite volume spatial
! discretisation using first order upwind scheme solved implicitly.
! See https://link.springer.com/article/10.1007/s13344-016-0039-1

    use meshparams
    implicit none

    integer :: nb

    double precision, intent(in) :: dt
    double precision, intent(out) :: lcoeff(nb,13)

    integer :: k, p
    double precision :: dotv

    lcoeff = 0.

    ! Upwind Scheme for Solving Advection Equations
    do k = 1, nb
      if(FVarea(k)>0)then
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            dotv = faceVel(k,p)
            ! Matrix coefficients
            lcoeff(k,1) = lcoeff(k,1) - 0.5 * dt * (dotv - abs(dotv)) * FVvDist(k,p) / FVarea(k)
            lcoeff(k,p+1) = 0.5 * dt * (dotv - abs(dotv)) * FVvDist(k,p) / FVarea(k)
          endif
        enddo
      endif
      lcoeff(k,1) = 1.0 + lcoeff(k,1)
    enddo

    return

end subroutine advecupwind

subroutine getrange(nb, data, dmin, dmax)
!*****************************************************************************
! For a considered variable get its range based on its neighbors values.

    use meshparams
    implicit none

    integer :: nb
    double precision, intent(in) :: data(nb)

    double precision, intent(out) :: dmin(nb)
    double precision, intent(out) :: dmax(nb)

    integer :: k, p, n

    dmin(1:nb) = data(1:nb)
    dmax(1:nb) = data(1:nb)
    do k = 1, nb
      if(FVarea(k)>0)then
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            dmin(k) = min(dmin(k),data(n))
            dmax(k) = max(dmax(k),data(n))
          endif
        enddo
      endif
    enddo

    return

end subroutine getrange

subroutine adveciioe(nb, dt, nbout, lcoeff, rcoeff)
!*****************************************************************************
! Define advection coefficients based on a finite volume spatial
! discretisation using the Inflow Implicit and Outflow Explicit scheme (I2EO).
! The approach is based on https://www.math.sk/mikula/mo-FVCA6.pdf

  use meshparams
  implicit none

  integer :: nb

  double precision, intent(in) :: dt

  integer, intent(out) :: nbout(nb)
  double precision, intent(out) :: lcoeff(nb,13)
  double precision, intent(out) :: rcoeff(nb,13)

  integer :: k, p
  double precision :: dotv, fluxin, fluxout

  lcoeff = 0.
  rcoeff = 0.
  nbout = 0
  
  ! The inflow implicit and outflow explicit approach is preferred for 2 main reasons: its application to bigger time steps and its non-dependence on the variables to advect in the matrix construction.

  ! Inflow-Implicit/Outflow-Explicit Scheme for Solving Advection Equations
  do k = 1, nb
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        if(FVvDist(k,p)>0.)then
          ! Defining the inward flux (inverse of the dot product btw face normal and velocity) 
          dotv = -faceVel(k,p)
          
          ! Number of nonzero outflows
          if(dotv < 0.) nbout(k) = nbout(k) + 1

          ! Similarly we consider that the advected variable at the face is defined by linear interpolation from each connected vertex.
          fluxin = 0.5 * dt * max(dotv, 0.) * FVvDist(k,p) / FVarea(k)
          fluxout = 0.5 * dt * min(dotv, 0.) * FVvDist(k,p) / FVarea(k)

          ! Left-hand side matrix coefficients (implicit solution)
          lcoeff(k,1) = lcoeff(k,1) + fluxin
          lcoeff(k,p+1) = -fluxin

          ! Right-hand side matrix coefficients (explicit solution)
          rcoeff(k,1) = rcoeff(k,1) - fluxout
          rcoeff(k,p+1) = fluxout
          
        endif
      enddo
    endif
    lcoeff(k,1) = 1.0 + lcoeff(k,1)
    rcoeff(k,1) = 1.0 + rcoeff(k,1)
  enddo

  return

end subroutine adveciioe

subroutine adveciioe2(nb, dt, nbout, var, vmin, vmax, lcoeff, rcoeff)
!*****************************************************************************
! Define advection coefficients based on a finite volume spatial
! discretisation using the Inflow Implicit and Outflow Explicit scheme (I2EO).
! The approach is based on https://www.sciencedirect.com/science/article/pii/S0168927414001032

  use meshparams
  implicit none

  integer :: nb

  double precision, intent(in) :: dt
  integer, intent(in) :: nbout(nb)
  double precision, intent(in) :: var(nb)
  double precision, intent(in) :: vmin(nb)
  double precision, intent(in) :: vmax(nb)
  double precision, intent(out) :: lcoeff(nb,13)
  double precision, intent(out) :: rcoeff(nb,13)

  integer :: k, p, n, q
  double precision :: thetain(nb,12), thetaout(nb,12)
  double precision :: aout, fluxin, fluxout
  double precision :: dotv, val

  lcoeff = 0.
  rcoeff = 0.

  ! Weighting parameter for flux-corrected transport for the stabilized IIOE scheme
  thetaout = 0.5
  do k = 1, nb
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        dotv = -faceVel(k,p)
        aout = min(dotv, 0.) * FVvDist(k,p) 
        if(aout*(var(n) - var(k))>0)then
          val = FVarea(k) * (vmax(k) - var(k))
          val = val / (aout * dt * nbout(k) * (var(n) - var(k)))
          thetaout(k,p) = min(0.5, val)
        elseif(aout*(var(n) - var(k))<0)then
          val = FVarea(k) * (vmin(k) - var(k))
          val = val / (dt * aout * nbout(k) * (var(n) - var(k)))
          thetaout(k,p) = min(0.5, val)
        else
          thetaout(k,p) = 0.5
        endif
      enddo
    endif
  enddo

  thetain = 0.5
  do k = 1, nb
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        q = FVnIDfNb(k,p)
        thetain(k,p) = 1.0 - thetaout(n,q)
      enddo
    endif
  enddo

  ! Inflow-Implicit/Outflow-Explicit Scheme for Solving Advection Equations
  do k = 1, nb
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        if(FVvDist(k,p)>0.)then
          ! Defining the inward flux (inverse of the dot product btw face normal and velocity) 
          dotv = -faceVel(k,p)

          ! Similarly we consider that the advected variable at the face is defined by linear interpolation from each connected vertex.
          fluxin = dt * thetain(k,p) * max(dotv, 0.) * FVvDist(k,p) / FVarea(k)
          fluxout = dt * thetaout(k,p) * min(dotv, 0.) * FVvDist(k,p) / FVarea(k)

          ! Left-hand side matrix coefficients (implicit solution)
          lcoeff(k,1) = lcoeff(k,1) + fluxin
          lcoeff(k,p+1) = -fluxin

          ! Right-hand side matrix coefficients (explicit solution)
          rcoeff(k,1) = rcoeff(k,1) - fluxout
          rcoeff(k,p+1) = fluxout
          
        endif
      enddo
    endif
    lcoeff(k,1) = 1.0 + lcoeff(k,1)
    rcoeff(k,1) = 1.0 + rcoeff(k,1)
  enddo

  return

end subroutine adveciioe2

subroutine fitedges(h, nh, nb)
!*****************************************************************************
! Set edges elevations based on local neighbours.

  use meshparams
  implicit none

  integer :: nb
  double precision, intent(in) :: h(nb)
  double precision, intent(out) :: nh(nb)

  integer :: f, k, p, n, m, t, l, id(nb), id2(nb)
  double precision :: s, val

  nh = 0.0

  id = -1
  m = 1
  do k = 1, nb
    if(h(k) <= -1.e7)then
      s = 0.0
      val = 0.0
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        if(h(n)>-1.e7)then
          s = s + FVeLgt(k,p)
          val = val + FVeLgt(k,p) * h(n)
        endif
      enddo
      if(s > 0.0)then
        nh(k) = val/s
      else
        id(m) = k
        m = m + 1
        nh(k) = -1.e7
      endif
    else
      nh(k) = h(k)
    endif
  enddo

  f = 1
  if(m == 1) f = 0
  
  do while(f == 1)
    l = 1
    id2 = -1
    do t = 1, m
      k = id(t)
      if(k>0)then
        s = 0.0
        val = 0.0
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          if(nh(n)>-1.e7)then
            s = s + FVeLgt(k,p)
            val = val + FVeLgt(k,p) * nh(n)
          endif
        enddo
        if(s > 0.0)then
          nh(k) = val/s
        else
          id2(l) = k
          l = l + 1
        endif
      endif
    enddo

    if(l == 1)then
      f = 0
    else
      m = l
      l = 1
      id = id2
    endif
  enddo

  return
  
end subroutine fitedges

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
    double precision :: c, ck, cn, v

    dcoeff = 0.
    do k = 1, nb
      if(FVarea(k)>0)then
        ck = Kd(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            n = FVnID(k,p)+1
            cn = Kd(n)
            c = 0.5*(ck+cn)/FVarea(k)
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            dcoeff(k) = dcoeff(k) - v*h(n)
            dcoeff(k) = dcoeff(k) + v*h(k)
          endif
        enddo
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
    double precision, intent(out) :: dcoeff(nb,13)

    integer :: k, p, n
    double precision :: c, ck, cn, cpk, cpn, v

    dcoeff = 0.
    do k = 1, nb
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

subroutine donorslist(nrcv, inIDs, rcvs, donors, nb)
!*****************************************************************************
! Compute the list of donors based on mesh connectivity.

  use meshparams
  implicit none

  integer :: nb 

  integer, intent(in) :: nrcv
  integer, intent(in) :: inIDs(nb)
  integer, intent(in) :: rcvs(nb,nrcv)
  integer, intent(out) :: donors(nb,12)

  integer :: k, i, p
  integer :: nbdonors(nb)

  donors = -1
  nbdonors = 0

  do k = 1, nb
    if(inIDs(k)>0) p = 1
    do p = 1, nrcv
      i = rcvs(k,p) + 1
      if(i .ne. k .and. i > 0)then
        nbdonors(i) = nbdonors(i) + 1
        if(nbdonors(i) <= 12)then
          donors(i,nbdonors(i)) = k - 1
        endif
      endif
    enddo 
  enddo

  return

end subroutine donorslist

subroutine donorsmax(dat, donors, valmax, nb)
!*****************************************************************************
! Compute the donors maximum based on a specific variable.

  use meshparams
  implicit none

  integer :: nb 

  double precision, intent(in) :: dat(nb)
  integer, intent(in) :: donors(nb,12)
  double precision, intent(out) :: valmax(nb)

  integer :: k, i, p

  valmax = -1.e8 

  do k = 1, nb
    do p = 1, 12
      i = donors(k,p) + 1
      if(i > 0)then
        valmax(k) = max(dat(i),valmax(k))
      endif
    enddo
  enddo

  return

end subroutine donorsmax

subroutine mfdreceivers(nRcv, exp, elev, sl, rcv, dist, wgt, nb)
!*****************************************************************************
! Compute receiver characteristics based on multiple flow direction algorithm.
! The exponent is referred to as the flow‐partition exponent following: Quin et al., 2007
! An adaptive approach to selecting a flow‐partition exponent for a multiple‐flow‐direction algorithm
! The larger the value of the exponent, the more similar MFD is to SFD.
! The smaller the value of the exponent, the more spread is the distribution of water across downstream nodes.

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
  double precision,intent(in) :: sl
  double precision, intent(in) :: elev(nb)

  integer, intent(out) :: rcv(nb,nRcv)
  double precision, intent(out) :: dist(nb,nRcv)
  double precision, intent(out) :: wgt(nb,nRcv)

  integer :: k, n, p, kk
  double precision :: slp(12),dst(12),val,slope(12)
  double precision :: e, fexp
  integer :: id(12)

  rcv = -1
  dist = 0.
  wgt = 0.

  do k = 1, nb
    if(elev(k)<=sl)then
      rcv(k,1:nRcv) = k-1
    else
      ! Determination of flow-partition exponent using the maximum downslope gradient (Qin et al. 2007)
      e = 0.
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        if(n>0 .and. FVeLgt(k,p)>0.)then
          val = (elev(k) - elev(n))/FVeLgt(k,p) 
          e = max(val,e) 
        endif
      enddo
      if(e>0)then
        fexp = 5.0 * min(e,1.0) + exp
      else
        fexp = 1.0
      endif
      ! fexp = exp
      slp = 0.
      id = 0
      val = 0.
      kk = 0
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        if(n>0 .and. FVeLgt(k,p)>0.)then
          val = (elev(k) - elev(n))**fexp/FVeLgt(k,p)
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
  enddo

  return

end subroutine mfdreceivers

subroutine mfdrcvrs(nRcv, exp, elev, sl, rcv, dist, wgt, nb)
!*****************************************************************************
! Compute receiver characteristics based on multiple flow direction algorithm for marine environment.
! The exponent is referred to as the flow‐partition exponent following: Quin et al., 2007
! An adaptive approach to selecting a flow‐partition exponent for a multiple‐flow‐direction algorithm
! The larger the value of the exponent, the more similar MFD is to SFD.
! The smaller the value of the exponent, the more spread is the distribution of water across downstream nodes.

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
  double precision,intent(in) :: sl
  double precision, intent(in) :: elev(nb)

  integer, intent(out) :: rcv(nb,12)
  double precision, intent(out) :: dist(nb,12)
  double precision, intent(out) :: wgt(nb,12)

  integer :: k, n, p, kk, ngbs
  double precision :: fexp,slp(12),dst(12),val,slope(12)
  integer :: id(12)

  rcv = -1
  dist = 0.
  wgt = 0.

  do k = 1, nb
    if(elev(k)>sl)then
      ngbs = nRcv
      fexp = exp
    else
      ngbs = 12
      fexp = 0.01
    endif
    slp = 0.
    id = 0
    val = 0.
    kk = 0
    do p = 1, FVnNb(k)
      n = FVnID(k,p)+1
      if(n>0 .and. FVeLgt(k,p)>0.)then
        val = (elev(k) - elev(n))**fexp/FVeLgt(k,p)
        if(val>0.)then
          kk = kk + 1
          slp(kk) = val
          id(kk) = n-1
          dst(kk) = FVeLgt(k,p)
        endif
      endif
    enddo

    if(kk == 0)then
      rcv(k,1:12) = k-1
    elseif(kk <= ngbs)then
      val = 0.
      rcv(k,1:ngbs) = k-1
      do p = 1, kk
        rcv(k,p) = id(p)
        dist(k,p) = dst(p)
        val = val + slp(p)
      enddo
      do p = 1, ngbs
        wgt(k,p) = slp(p) / val
      enddo
    else
      rcv(k,1:ngbs) = k-1
      call quicksort(slp,1,kk,id)
      n = 0
      val = 0.
      slope = 0.
      do p = kk,kk-ngbs+1,-1
        n = n + 1
        slope(n) = slp(p)
        rcv(k,n) = id(p)
        dist(k,n) = dst(p)
        val = val + slp(p)
      enddo
      do p = 1, ngbs
        wgt(k,p) = slope(p)/val
      enddo
    endif
  enddo

  return

end subroutine mfdrcvrs

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
  double precision,intent(out) :: newgraph(n*12,4)

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

subroutine sort_ids(df1, df2, id2, m)
!*****************************************************************************
! Sort pit ids.

  use meshparams
  implicit none

  integer :: m
  integer,intent(in) :: df1(m)
  integer,intent(in) :: df2(m)

  integer, intent(out) :: id2(m)
  integer :: k

  id2 = -1

  do k = 1, m
    if(k == 1)then
      id2(k) = df2(1)
    else
      if(df2(k) == df2(k-1))then
        id2(k) = df1(k-1)
      else
        id2(k) = df2(k)
      endif
    endif
  enddo

  return

end subroutine sort_ids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!          MESH DECLARATION FUNCTIONS              !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine globalngbhs(nt, cells, n)

  use meshparams
  implicit none

  integer :: n
  integer, intent(in) :: nt
  integer, intent(in) :: cells(n, 3)

  integer :: i, k, nb, n1, n2, n3
  integer :: nc11, nc12, nc21, nc22, nc31, nc32
  integer :: nc(3)

  logical :: in11, in12, in21, in22, in31, in32

  ! Define mesh parameters
  if(allocated(FVgnID)) deallocate(FVgnID)
  if(allocated(FVgnNb)) deallocate(FVgnNb)
  allocate(FVgnID(nt,12))
  allocate(FVgnNb(nt))

  FVgnID = -1
  FVgnNb = 0

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells(i,1:3)+1

    n1 = nc(1)
    nc11 = nc(2)
    nc12 = nc(3)
    n2 = nc(2)
    nc21 = nc(1)
    nc22 = nc(3)
    n3 = nc(3)
    nc31 = nc(1)
    nc32 = nc(2)

    in11 = .True.
    in12 = .True.
    in21 = .True.
    in22 = .True.
    in31 = .True.
    in32 = .True.
    do k = 1, 12
      if(FVgnID(n1,k)==nc11) in11 = .False.
      if(FVgnID(n1,k)==nc12) in12 = .False.
      if(FVgnID(n2,k)==nc21) in21 = .False.
      if(FVgnID(n2,k)==nc22) in22 = .False.
      if(FVgnID(n3,k)==nc31) in31 = .False.
      if(FVgnID(n3,k)==nc32) in32 = .False.
    enddo

    if(in11)then
      nb = FVgnNb(n1) + 1
      FVgnNb(n1) = nb
      FVgnID(n1, nb) = nc11
    endif

    if(in12)then
      nb = FVgnNb(n1) + 1
      FVgnNb(n1) = nb
      FVgnID(n1, nb) = nc12
    endif

    if(in21)then
      nb = FVgnNb(n2) + 1
      FVgnNb(n2) = nb
      FVgnID(n2, nb) = nc21
    endif

    if(in22)then
      nb = FVgnNb(n2) + 1
      FVgnNb(n2) = nb
      FVgnID(n2, nb) = nc22
    endif

    if(in31)then
      nb = FVgnNb(n3) + 1
      FVgnNb(n3) = nb
      FVgnID(n3, nb) = nc31
    endif

    if(in32)then
      nb = FVgnNb(n3) + 1
      FVgnNb(n3) = nb
      FVgnID(n3, nb) = nc32
    endif
  enddo

  return

end subroutine globalngbhs 

subroutine epsfill(sl, elev, fillz, nb)
!*****************************************************************************
! Perform pit filling using a priority queue approach following Barnes (2015).
! This function is done on a single processors.

  use meshparams
  implicit none

  integer :: nb
  double precision,intent(in) :: sl
  double precision,intent(in) :: elev(nb)
  double precision,intent(out) :: fillz(nb)

  logical :: flag(nb)

  integer :: i, k, c

  type (node)  :: ptID
  double precision :: h

  fillz = elev

  ! Push marine edges nodes to priority queue
  flag = .False.
  do i = 1, nb
    if(fillz(i)<sl)then
      flag(i) = .True.
      lp: do k = 1, FVgnNb(i)
        c = FVgnID(i,k)
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
    do k = 1, FVgnNb(i)
      c = FVgnID(i,k)
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
            call priorityqueue%PQpush(fillz(c), c)
          endif
        endif
      endif
    enddo
  enddo

  return

end subroutine epsfill

subroutine updatearea(narea, nb)
!*****************************************************************************
! Update voronoi areas.

  use meshparams
  implicit none

  integer :: nb
  double precision, intent(in) :: narea(nb)

  FVarea = narea

  return

end subroutine updatearea 

subroutine definetin( coords, cells_nodes, cells_edges, edges_nodes, &
                      circumcenter, ngbID, narea, n, nb, m)
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
  double precision, intent(in) :: circumcenter(3,n)

  integer, intent(out) :: ngbID(nb, 12)
  double precision, intent(out) :: narea(nb)

  integer :: i, n1, n2, k, l, p, eid, cid
  integer :: e, id, nbngbs, ngbIDv
  integer :: nid(2), nc(3), edge(nb, 12)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 12)

  double precision :: coords0(3), coordsID(3), cc(3), triarea
  double precision :: midpoint(3), dist, vnorm, radius

  logical :: inside

  ! Define mesh parameters
  if(allocated(FVarea)) deallocate(FVarea)
  if(allocated(FVnID)) deallocate(FVnID)
  if(allocated(FVnNb)) deallocate(FVnNb)
  if(allocated(FVeLgt)) deallocate(FVeLgt)
  if(allocated(FVvDist)) deallocate(FVvDist)
  if(allocated(lcoords)) deallocate(lcoords)
  if(allocated(faceVec)) deallocate(faceVec)
  if(allocated(FVnIDfNb)) deallocate(FVnIDfNb)
  if(allocated(midFace)) deallocate(midFace)

  allocate(FVarea(nb))
  allocate(FVnNb(nb))
  allocate(FVnID(nb,12))
  allocate(FVeLgt(nb,12))
  allocate(FVvDist(nb,12))
  allocate(lcoords(nb,3))
  allocate(faceVec(nb,12,3))
  allocate(FVnIDfNb(nb,12))

  cell_ids = -1
  edge = -1
  FVnNb = 0
  ngbID = -1
  FVeLgt = 0.
  FVvDist = 0.
  FVarea = 0.
  FVnIDfNb = -1

  lcoords = coords
  radius = 0.
  if(lcoords(1,3) .ne. 0.0) radius = norm2(lcoords(1,1:3))
  if(radius > 0.)then
    allocate(midFace(nb,12,3))
    midFace = -1.0
  endif

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
      FVnNb(n1) = FVnNb(n1) + 1
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
      call distance(coords0, coordsID, FVeLgt(k,l))
    enddo

    ! Get voronoi edge lengths
    lp2: do cid = 1, 12
      if( cell_ids(k,cid) == -1 ) exit lp2
      edgeNb(1:3) = cells_edges( cell_ids(k,cid)+1,1:3 )
      do e = 1, 3
        edges(e,1:2) = edges_nodes(edgeNb(e)+1,1:2)
        if( k-1 == edges(e,1) .or. k-1 == edges(e,2))then
          ! Midpoint on an arc of the sphere or on flat mesh
          if(radius > 0)then 
            call slerpmidpt(coords(edges(e,1)+1,1:3), coords(edges(e,2)+1,1:3), &
                    midpoint(1:3))
          else
            midpoint(1:3) = 0.5 * (coords(edges(e,1)+1,1:3)+coords(edges(e,2)+1,1:3))
          endif
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
          cc = circumcenter(1:3, cell_ids(k,cid)+1)
          call distance(midpoint(1:3), cc,  dist)
          FVvDist(k,id) = FVvDist(k,id) + dist
          ! Get the tangent vector either from euclidian or spherical formula
          if(radius == 0.0)then
            call euclid(coords0(1:3), midpoint(1:3),  vnorm)
            faceVec(k,id,1:3) = (midpoint(1:3) - coords0(1:3)) / vnorm
            call gtriarea(coords0(1:3), midpoint(1:3), cc, triarea)
            FVarea(k) = FVarea(k) + triarea
          else
            call slerpvec(coords0(1:3), midpoint(1:3), 0.9_8, faceVec(k,id,1:3))
            ! call gtriareasphere(coords0(1:3), midpoint(1:3), cc, triarea)
            call striarea(coords0(1:3), midpoint(1:3), cc, triarea)
            FVarea(k) = FVarea(k) + triarea
            midFace(k,id,1:3) = midpoint(1:3) / norm2(midpoint(1:3))
          endif
        endif
      enddo
    enddo lp2
  enddo

  FVnID = ngbID
  narea = FVarea

  do k = 1, nb
    if(isnan(FVarea(k))) FVarea(k) = 0.0
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        inner: do l = 1, FVnNb(n)
          if(k == FVnID(n,l)+1)then
            FVnIDfNb(k,p) = l
            exit inner
          endif
        enddo inner
      enddo
    endif
  enddo

  return

end subroutine definetin

subroutine stencil(nb,ngbid,maxnb)
!*****************************************************************************
! Compute the neighbors of the neighborhood used for the flexural equation
! solution

  use meshparams
  implicit none
  
  integer, intent(in) :: nb
  integer, intent(out) :: ngbid(nb, 41)
  integer, intent(out) :: maxnb

  integer :: k, p, n, q, m, l, i, res(41)
  integer :: stenNb(nb), stenNgb(nb, 65)
  logical :: add

  if(allocated(stencilNb)) deallocate(stencilNb)
  allocate(stencilNb(nb))
  if(allocated(stencilNgb)) deallocate(stencilNgb)
  allocate(stencilNgb(nb,41))
  stenNb = 0
  stenNgb = -1

  stencilNb = 0
  stencilNgb = -1

  ! Find all cells surrounding a given vertice
  do k = 1, nb
    stenNgb(k,1:FVnNb(k)) = FVnID(k,1:FVnNb(k))
    stenNb(k) = FVnNb(k)
    if(FVarea(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        do l = 1, FVnNb(n)
          m = FVnID(n,l)+1
          add = .True.
          if(add)then
            stenNb(k) = stenNb(k)+1
            stenNgb(k,stenNb(k)) = m-1
          endif
        enddo
      enddo
    endif
  enddo

  ngbid = -1
  maxnb = 0
  do q = 1, nb
    res(1) = q-1
    res(2:FVnNb(q)+1) = FVnID(q,1:FVnNb(q))
    res(FVnNb(q)+2:41) = -1
    k = FVnNb(q)+1
    do i=1,65
        if (any( res == stenNgb(q,i) )) cycle
        k = k + 1
        if(k>41)then
          write(*,*)'There is more than 40 neighbors on the second-order stencil,'
          write(*,*)' you will need to increase the allocation array accordingly!'
          exit
        endif
        res(k) = stenNgb(q,i)
    enddo
    ngbid(q,1:k) = res(1:k)
    maxnb = max(k,maxnb)
    stencilNb(q) = k
  enddo
  
  stencilNgb = ngbid

  return

end subroutine stencil

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
  integer,intent(in) :: ngbid(nb, 12)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                  !!
!!           FLEXURE PROCESSES FUNCTIONS            !!
!!                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine four1(data1,nn,isign)
!*****************************************************************************
! This function is from Fastscape.

  integer :: isign,nn
  double precision :: data1(2*nn)
  integer i,istep,j,m,mmax,n
  double precision :: tempi,tempr
  double precision :: theta,wi,wpi,wpr,wr,wtemp
  
  n=2*nn
  j=1
  do i=1,n,2
    if(j>i)then
      tempr=data1(j)
      tempi=data1(j+1)
      data1(j)=data1(i)
      data1(j+1)=data1(i+1)
      data1(i)=tempr
      data1(i+1)=tempi
    endif
    m=n/2
    do while((m>=2).and.(j>m)) 
      j=j-m
      m=m/2
    enddo
    j=j+m
  enddo
  
  mmax=2
  do while(n>mmax) 
    istep=2*mmax
    theta=6.28318530717959d0/(isign*mmax)
    wpr=-2.d0*dsin(0.5d0*theta)**2
    wpi=dsin(theta)
    wr=1.d0
    wi=0.d0
    do m=1,mmax,2
      do i=m,n,istep
        j=i+mmax
        tempr=sngl(wr)*data1(j)-sngl(wi)*data1(j+1)
        tempi=sngl(wr)*data1(j+1)+sngl(wi)*data1(j)
        data1(j)=data1(i)-tempr
        data1(j+1)=data1(i+1)-tempi
        data1(i)=data1(i)+tempr
        data1(i+1)=data1(i+1)+tempi
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    enddo
    mmax=istep
  enddo

end subroutine four1

subroutine realft(data1,n,isign)
!*****************************************************************************
! This function extracts the real part of the FFT, the routine is from Fastscape.

  double precision :: wr,wi,wpi,wpr,wtemp,theta
  double precision :: data1(2*n)
  integer :: n,isign,i1,i2,i3,i4,i
  double precision :: c2,h2r,h2i,wis,h1r,h1i
  
  theta=6.28318530717959d0/2.0d0/dfloat(n)
  c1=0.5
  if(isign==1)then
    c2=-0.5
    call four1(data1,n,+1)
  else
    c2=0.5
    theta=-theta
  endif
  
  wpr=-2.0d0*dsin(0.5d0*theta)**2
  wpi=dsin(theta)
  wr=1.0d0+wpr
  wi=wpi
  n2p3=2*n+3
  
  do i=2,n/2+1
    i1=2*i-1
    i2=i1+1
    i3=n2p3-i2
    i4=i3+1
    wrs=sngl(wr)
    wis=sngl(wi)
    h1r=c1*(data1(i1)+data1(i3))
    h1i=c1*(data1(i2)-data1(i4))
    h2r=-c2*(data1(i2)+data1(i4))
    h2i=c2*(data1(i1)-data1(i3))
    data1(i1)=h1r+wrs*h2r-wis*h2i
    data1(i2)=h1i+wrs*h2i+wis*h2r
    data1(i3)=h1r-wrs*h2r+wis*h2i
    data1(i4)=-h1i+wrs*h2i+wis*h2r
    wtemp=wr
    wr=wr*wpr-wi*wpi+wr
    wi=wi*wpr+wtemp*wpi+wi
  enddo
  
  if(isign==1)then
    h1r=data1(1)
    data1(1)=h1r+data1(2)
    data1(2)=h1r-data1(2)
  else
    h1r=data1(1)
    data1(1)=c1*(h1r+data1(2))
    data1(2)=c1*(h1r-data1(2))
    call four1(data1,n,-1)
  endif

end subroutine realft

subroutine sinft(y,n)
!*****************************************************************************
! This function computes FFT of weights, the routine is from Fastscape.

  integer :: n
  double precision, intent(inout) :: y(n)
  
  integer :: j
  double precision :: sum,y1,y2
  double precision :: theta,wi,wpi,wpr,wr,wtemp
  
  theta=3.141592653589793d0/dble(n)
  wr=1.0d0
  wi=0.0d0
  wpr=-2.0d0*sin(0.5d0*theta)**2
  wpi=sin(theta)
  y(1)=0.0
  
  do j = 1,n/2
    wtemp=wr
    wr=wr*wpr-wi*wpi+wr
    wi=wi*wpr+wtemp*wpi+wi
    y1=wi*(y(j+1)+y(n-j+1))
    y2=0.5*(y(j+1)-y(n-j+1))
    y(j+1)=y1+y2
    y(n-j+1)=y1-y2
  enddo
  
  call realft(y,n/2,+1)
  
  sum=0.0
  y(1)=0.5*y(1)
  y(2)=0.0
  
  do j= 1, n-1, 2
    sum=sum+y(j)
    y(j)=y(j+1)
    y(j+1)=sum
  enddo
  
end subroutine sinft

subroutine addw(w,nxflex,nyflex,iflexmin,iflexmax,jflexmin,jflexmax,cbc)
!*****************************************************************************
! This function add flexural weights, the routine is from Fastscape

  implicit none

  integer :: nxflex,nyflex,i,j,iflexmin,iflexmax,jflexmin,jflexmax
  double precision w(nxflex,nyflex)
  character cbc*4

  if (cbc(1:1).eq.'0') w(:,jflexmin)=w(:,jflexmin+1)
  if (cbc(2:2).eq.'0') w(iflexmax,:)=w(iflexmax-1,:)
  if (cbc(3:3).eq.'0') w(:,jflexmax)=w(:,jflexmax-1)
  if (cbc(4:4).eq.'0') w(iflexmin,:)=w(iflexmin+1,:)

  do j=jflexmin,jflexmax
    do i=iflexmin,iflexmax
      if (cbc(1:1).eq.'0') w(i,jflexmin-(j-jflexmin+1))=w(i,jflexmin-(j-jflexmin+1))+w(i,j)
      if (cbc(3:3).eq.'0') w(i,jflexmax+(jflexmax-j+1))=w(i,jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),j)=w(iflexmin-(i-iflexmin+1),j)+w(i,j)
      if (cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),j)=w(iflexmax+(iflexmax-i+1),j)+w(i,j)
      if (cbc(1:1).eq.'0'.and.cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))= &
      w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))+w(i,j)
      if (cbc(2:2).eq.'0'.and.cbc(3:3).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))= &
      w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(3:3).eq.'0'.and.cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))= &
      w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(4:4).eq.'0'.and.cbc(1:1).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))= &
      w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))+w(i,j)
    enddo
  enddo

  return

end subroutine addw

subroutine flexure(dh,nx,ny,xl,yl,young,nu,rhos,rhoa,eet,ibc,newh)
!*****************************************************************************
! Routine to compute the flexural response of erosion from Fastscape
! in input:
! dh(nx,ny), the difference in topography after and before erosion, in m,
! rhos(nx,ny), surface rock density, in kg/m^3
! rhoa, asthenospheric density, in kg/m^3
! eet, effective elastic thickness, in m
! nx,ny, resolution of the input topography
! xl,yl, horizontal dimensions of the input topography, in m

! Here fixed values are assumed for:
! Young modulus, 1.d11 Pa
! Poisson ratio, 0.25
! g, 9.81 m/s^2

! the flexural, biharmonic equation is solved by FFT method (see Nunn and Aires, 1988)
! on a 1024x1024 mesh

  implicit none

  integer, intent(in) :: nx,ny,ibc
  double precision, intent(in), dimension(nx,ny) :: dh
  double precision, intent(in) :: xl,yl,rhoa,eet,rhos,young,nu
  double precision, intent(out), dimension(nx,ny) :: newh

  integer nxflex,nyflex,i,j,ii,jj
  double precision, dimension(:,:), allocatable :: w
  double precision hx,hy,dflex,d,xk,pihx,pihy,g,fi,fj,tij,dx,dy,r,s,h1,h2,h3,h4
  double precision ddxf,ddyf,xloc,yloc,dw,xflexloc,yflexloc
  integer iflexmin,iflexmax,jflexmin,jflexmax
  character cbc*4

  double precision, dimension(:,:,:), allocatable :: hw
  integer, dimension(:,:), allocatable :: iiw,jjw

  write (cbc,'(i4)') ibc
  if (cbc(1:1).eq.'') cbc(1:1)='0'
  if (cbc(2:2).eq.'') cbc(2:2)='0'
  if (cbc(3:3).eq.'') cbc(3:3)='0'
  if (cbc(4:4).eq.'') cbc(4:4)='0'

  ! allocate memory

  nxflex=1
  do while (nxflex.lt.nx)
    nxflex=nxflex*2
  enddo

  nyflex=1
  do while (nyflex.lt.ny)
    nyflex=nyflex*2
  enddo

  allocate (hw(4,nx,ny), iiw(nx,ny), jjw(nx,ny))
  allocate (w(nxflex,nyflex))

  ! compute relevant geometrical, flexural and spectral parameters

  iflexmin=nxflex/2-nxflex/8
  iflexmax=nxflex/2+nxflex/8
  jflexmin=nyflex/2-nyflex/8
  jflexmax=nyflex/2+nyflex/8

  dx=xl/(nx-1)
  dy=yl/(ny-1)
  ddxf=xl/(iflexmax-iflexmin)
  ddyf=yl/(jflexmax-jflexmin)
  hx=ddxf*(nxflex-1)
  hy=ddyf*(nyflex-1)
  dflex=young/12.d0/(1.d0-nu**2)
  d=dflex*eet**3
  g=9.81d0
  xk=rhoa*g
  pihx=3.141592654d0/hx
  pihy=3.141592654d0/hx

  ! compute weigths corresponding to the increase in topography by interpolation
  ! from the nx,ny grid to the nflex, nflex grid, using a bilinear interpolation scheme

  w=0.d0
  jj=jflexmin
  yflexloc=0.d0
  do j=1,ny
    yloc=(j-1)*dy
    if (yloc.gt.yflexloc+ddyf) jj=jj+1
    yflexloc=(jj-jflexmin)*ddyf
    ii=iflexmin
    xflexloc=0.d0
    do i=1,nx
      xloc=(i-1)*dx
      if (xloc.gt.xflexloc+ddxf) ii=ii+1
      xflexloc=(ii-iflexmin)*ddxf
      r=(xloc-xflexloc)/ddxf*2.d0-1.d0
      s=(yloc-yflexloc)/ddyf*2.d0-1.d0
      h1=(1.d0-r)*(1.d0-s)/4.d0
      h2=(1.d0+r)*(1.d0-s)/4.d0
      h3=(1.d0-r)*(1.d0+s)/4.d0
      h4=(1.d0+r)*(1.d0+s)/4.d0
      iiw(i,j)=ii
      jjw(i,j)=jj
      hw(1,i,j)=h1
      hw(2,i,j)=h2
      hw(3,i,j)=h3
      hw(4,i,j)=h4
      dw=-dh(i,j)*rhos*dx*dy*g
      w(ii,jj) = w(ii,jj) + dw*h1
      w(ii+1,jj) = w(ii+1,jj) + dw*h2
      w(ii,jj+1) = w(ii,jj+1) + dw*h3
      w(ii+1,jj+1) = w(ii+1,jj+1) + dw*h4
    enddo
  enddo

  call addw(w,nxflex,nyflex,iflexmin,iflexmax,jflexmin,jflexmax,cbc)

  ! compute FFT of weights

  do j=1,nyflex
    call sinft(w(:,j),nxflex)
  enddo

  w=transpose(w)

  do i=1,nxflex
    call sinft(w(:,i),nyflex)
  enddo

  ! apply filter to FFT of weights to simulated flexure (see Nunn and Aires, 1988)

  w=w*4./hx/hy

  do j=1,nyflex
    fj=(j*pihx)**2
    do i=1,nxflex
      fi=(i*pihx)**2
      tij=d/xk*(fi**2+2.d0*fi*fj+fj**2)+1.d0
      w(j,i)=w(j,i)/xk/tij
    enddo
  enddo

  ! compute inverse FFT of filtered weights to obtain deflection
  do i=1,nxflex
    call sinft (w(:,i),nyflex)
  enddo

  w=transpose(w)

  do j=1,nyflex
    call sinft (w(:,j),nxflex)
  enddo

  ! add  deflection by interpolation from the nflex,nflex grid to the nx,ny grid
  ! by bilinear interpolation
  do j=1,ny
    do i=1,nx
      ii=iiw(i,j)
      jj=jjw(i,j)
      h1=hw(1,i,j)
      h2=hw(2,i,j)
      h3=hw(3,i,j)
      h4=hw(4,i,j)
      newh(i,j)=w(ii,jj)*h1+w(ii+1,jj)*h2+w(ii,jj+1)*h3+w(ii+1,jj+1)*h4
    enddo
  enddo

  ! deallocate memory
  deallocate (w,iiw,jjw,hw)

end subroutine flexure

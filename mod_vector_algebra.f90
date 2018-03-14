! © 2017 DANG Truong
module mod_vector_algebra
  implicit none

! define here the vector type and overload operators
  type :: vector_t
      real(8) :: x = 0.0d0
      real(8) :: y = 0.0d0
      real(8) :: z = 0.0d0
  end type vector_t

  interface operator(+)
      module procedure vector_sum
  end interface

  interface operator(-)
      module procedure vector_subs
  end interface

  interface operator(*)
    module procedure vector_dot_product
  end interface

  interface operator(.x.)
    module procedure vector_cross_product
  end interface

  interface operator(.abs.)
    module procedure vector_norm
  end interface
  
  interface operator(.times.)
    module procedure vector_times_scalar
  end interface

contains
! implement the functions

  function vector_sum(v1, v2) result(v3)
    type(vector_t), intent(in) :: v1, v2
    type(vector_t)             :: v3

    v3%x = v1%x + v2%x
    v3%y = v1%y + v2%y
    v3%z = v1%z + v2%z
!
  end function vector_sum

  function vector_subs(v1, v2) result(v3)
    type(vector_t), intent(in) :: v1, v2
    type(vector_t)             :: v3

    v3%x = v1%x - v2%x
    v3%y = v1%y - v2%y
    v3%z = v1%z - v2%z
!
  end function vector_subs

  function vector_dot_product(v1, v2) result(d)
    type(vector_t), intent(in) :: v1, v2
    real(8) :: d

    d = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
!
  end function vector_dot_product

  function vector_cross_product(v1, v2) result(v3)
    type(vector_t), intent(in) :: v1, v2
    type(vector_t)             :: v3

    v3%x = v1%y * v2%z - v1%z * v2%y
    v3%y = v1%z * v2%x - v1%x * v2%z
    v3%z = v1%x * v2%y - v1%y * v2%x
!
  end function vector_cross_product

  function vector_norm(v) result(l)
    type(vector_t), intent(in) :: v
    real(8) :: l

    l = sqrt(v%x**2 + v%y**2 + v%z**2 )
!
  end function vector_norm
  
  function vector_times_scalar(v1, v2) result(v3)
    type(vector_t), intent(in) :: v1
    real(8),        intent(in) :: v2
    type(vector_t)             :: v3

    v3%x = v1%x * v2
    v3%y = v1%y * v2
    v3%z = v1%z * v2
!
  end function vector_times_scalar

end module mod_vector_algebra

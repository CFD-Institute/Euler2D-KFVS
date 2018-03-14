module mod_point_util
    use mod_read_gmsh, only: point
    implicit none
    
    interface operator(.eq.)
        module procedure equal
    end interface
    
    contains
!----------------------------------------------------------------------     
        function equal(p1, p2) result(res)
        implicit none 
        logical :: res
        type(point), intent(in) :: p1, p2
        
        if (abs(p1%x - p2%x) <= 1.0d-9 .and. abs(p1%y - p2%y) <= 1.0d-9 .and. abs(p1%z -p2%z) <= 1.0d-9) then 
            res = .true.
        else
            res = .false.
        endif
        
        end function equal
!---------------------------------------------------------------------- 
        subroutine swap_point(p1,p2)
        implicit none 
        type(point), intent(inout) :: p1, p2
        type(point)                :: temp
        
        temp = p1
        p1   = p2
        p2   = temp
        end subroutine swap_point
end module mod_point_util
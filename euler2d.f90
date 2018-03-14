!  euler2d.f90 

!****************************************************************************
!
!  PROGRAM: euler2d
!
!  PURPOSE:  Illustrate the using of gmsh-to-vtk-and-tecplot-Fortran2003.
!
!****************************************************************************

    program euler2d

    implicit none

    call pre_processing ! refactoring code of gmsh-to-vtk-and-tecplot
    
    call solver

    end program euler2d


subroutine solver
    implicit none
    
    call solver_kfvs 
end subroutine solver
    
subroutine solver_kfvs
    use mod_solver_kfvs
    use mod_struct_to_array
    use mod_read_gmsh,     only: fname
    implicit none 
    integer :: n
    
    call struct_to_array 
    
    call donnee_initiale
    
    if (fname(1:3) == 'c31') then 
        call allocate_vardummy_multi_elem_airfoil
        call conditions_aux_limites_multi_elem_airfoil
        call assign_lr_cell_multi_elem_airfoil 
        call calcul_derived_quantities_multi_elem_airfoil 
    else
        call allocate_vardummy
        call conditions_aux_limites
        call assign_lr_cell 
        call calcul_derived_quantities 
    endif
    
    call calcul_conservative_vector 
    
    call timestep 
    
    !--- temps de simulation et parametres de sorties
    tmax=0.3828823925d-2
    nmax=floor(tmax/dt)
    
    ! debug
    nmax = 5000000
    
    !--- evolution
    do n=1,nmax
        
        !-- calcul des flux
        if (fname(1:3) == 'c31') then 
            call calcul_flux_multi_elem_airfoil 
        else
            call calcul_flux
        endif 
        
        !-- iteration en temps
        call calcul_rhs
        call euler_time_iteration
        
        !-- mise a jour
        vect_u=vect_unew
        
        !-- calcul de rho,ux,uy,t
        call calcul_rho_ux_uy_t
        
        !-- mise a jour cl
        if (fname(1:3) == 'c31') then 
            call conditions_aux_limites_multi_elem_airfoil 
        else
            call conditions_aux_limites
        endif 
        
        !-- mise à jour des quantités dérivées
        if (fname(1:3) == 'c31') then 
            call calcul_derived_quantities_multi_elem_airfoil 
        else
            call calcul_derived_quantities
        endif 
        
        !-- sauvegarde resultats (format vtk, lisible par Paraview)
        if (mod(n,1000) == 0) then
            write(*,*) 'Writing solution file at iteration ', n, '...'
            call write_solution_vtk(n)
            call write_pressure_coefficient(n)
        endif 
    enddo 
    
end subroutine solver_kfvs
    
subroutine write_solution_vtk(iter)
      use mod_solver_kfvs
      implicit none

      integer, intent(in) :: iter
      integer(4)          :: i
      character(300)      :: foutput
      character(7)        :: citer
      
      write(citer,'(I7.7)') iter
      foutput = trim(fname)//'_'//trim(citer)//'.vtk'

      open(unit = 21, file = foutput, status = 'replace')

      write(21,'(a)') '# vtk DataFile Version 2.0'
      write(21,'(a)') 'VTK format for unstructured mesh'
      write(21,'(a)') 'ASCII'
      write(21,'(a)') 'DATASET POLYDATA'
      write(21,1) list_cell%nbnode

      do i = 1, list_cell%nbnode
        write(21,*) list_cell%coord_nodes(i)%p%x, list_cell%coord_nodes(i)%p%y, list_cell%coord_nodes(i)%p%z
      end do

      write(21,2) list_cell%nbelm, 5*list_cell%nbelm

      do i = 1, list_cell%nbelm
        write(21,*) 4, list_cell%id_nodes(i)%pn%id_node(6:9)-1
      end do

      write(21,3) list_cell%nbelm
      write(21,'(a)') 'SCALARS CELL_IDENT integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, list_cell%nbelm
        write(21,*) list_cell%id_nodes(i)%pn%id_node(1)
      end do

      write(21,'(a)') 'SCALARS Density float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) rho(i)
      enddo 
      
      write(21,'(a)') 'SCALARS Temperature float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) t(i)
      enddo 
      
      write(21,'(a)') 'SCALARS Pressure float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) p(i)
      enddo
      
      write(21,'(a)') 'SCALARS Velocity_u float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) ux(i)
      enddo
      
      write(21,'(a)') 'SCALARS Velocity_v float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) uy(i)
      enddo
      
      write(21,'(a)') 'SCALARS Velocity_magnitude float'
      write(21,'(a)') 'LOOKUP_TABLE default '
      do i = 1, list_cell%nbelm
          write(21,*) sqrt(ux(i)**2 + uy(i)**2)
      enddo

      close(unit = 21)

1     format('POINTS',i9,' float')
2     format('POLYGONS ',2i9)
3     format('CELL_DATA',i9)

end subroutine    

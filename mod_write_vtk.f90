module mod_write_vtk
    use mod_read_gmsh
    use mod_cell_2D
    implicit none
    
    type :: obj_msh_vis
        class(obj_cell_2D), pointer :: smart_pointer
    contains
        procedure :: write_vtk     => write_mesh_vtk
        procedure :: write_tecplot => write_mesh_tecplot
    end type

    contains
!----------------------------------------------------------------------    
      subroutine write_mesh_vtk(this)
      implicit none

      class(obj_msh_vis), intent(in) :: this
      integer(4)                 :: i
      character(300)             :: foutput
      
      foutput = trim(fname)//'.vtk'

      open(unit = 21, file = foutput, status = 'replace')

      write(21,'(a)') '# vtk DataFile Version 2.0'
      write(21,'(a)') 'VTK format for unstructured mesh'
      write(21,'(a)') 'ASCII'
      write(21,'(a)') 'DATASET POLYDATA'
      write(21,1) this%smart_pointer%nbnode

      do i = 1, this%smart_pointer%nbnode
        write(21,*) this%smart_pointer%coord_nodes(i)%p%x, this%smart_pointer%coord_nodes(i)%p%y, this%smart_pointer%coord_nodes(i)%p%z
      end do

      write(21,2) this%smart_pointer%nbelm, 5*this%smart_pointer%nbelm

      do i = 1, this%smart_pointer%nbelm
        write(21,*) 4,this%smart_pointer%id_nodes(i)%pn%id_node(6:9)-1
      end do

      write(21,3) this%smart_pointer%nbelm
      write(21,'(a)') 'SCALARS CELL_IDENT integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, this%smart_pointer%nbelm
        write(21,*) this%smart_pointer%id_nodes(i)%pn%id_node(1)
      end do

      write(21,'(a)') 'SCALARS NEIGHBOR1 integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, this%smart_pointer%nbelm
          if (associated(this%smart_pointer%cell(i)%p%neighbor1)) then
              write(21,*) this%smart_pointer%cell(i)%p%neighbor1%ident
          else
              write(21,*) 0 
          endif 
      end do

      write(21,'(a)') 'SCALARS NEIGHBOR2 integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, this%smart_pointer%nbelm
          if (associated(this%smart_pointer%cell(i)%p%neighbor2)) then
              write(21,*) this%smart_pointer%cell(i)%p%neighbor2%ident
          else
              write(21,*) 0 
          endif 
      end do

      write(21,'(a)') 'SCALARS NEIGHBOR3 integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, this%smart_pointer%nbelm
          if (associated(this%smart_pointer%cell(i)%p%neighbor3)) then
              write(21,*) this%smart_pointer%cell(i)%p%neighbor3%ident
          else
              write(21,*) 0 
          endif 
      end do

      write(21,'(a)') 'SCALARS NEIGHBOR4 integer 1'
      write(21,'(a)') 'LOOKUP_TABLE default '

      do i = 1, this%smart_pointer%nbelm
          if (associated(this%smart_pointer%cell(i)%p%neighbor4)) then
              write(21,*) this%smart_pointer%cell(i)%p%neighbor4%ident
          else
              write(21,*) 0 
          endif 
      end do

      close(unit = 21)

1     format('POINTS',i9,' float')
2     format('POLYGONS ',2i9)
3     format('CELL_DATA',i9)

      end subroutine
 !----------------------------------------------------------------------     
      subroutine write_mesh_tecplot(this)
      implicit none

      class(obj_msh_vis), intent(in) :: this
      integer(4)                 :: i
      character(300)             :: foutput
      
      foutput = trim(fname)//'.dat'

      open(unit = 21, file = foutput, status = 'replace')

      !write(21,'(a)') 'VARIABLES=X,Y,CELL_IDENT,NEIGHBOR1,NEIGHBOR2,NEIGHBOR3,NEIGHBOR4' ! coming soon
      write(21,'(a)') 'VARIABLES=X,Y'
      write(21,'(a)') 'ZONE T="UNSTRUCTURED-COUNTOUR"'
      write(21,'(a)') 'ZONETYPE=FEPOLYGON'
      write(21,*) 'NODES=',this%smart_pointer%nbnode
      write(21,*) 'ELEMENTS=',this%smart_pointer%nbelm
      write(21,*) 'FACES=',this%smart_pointer%nbelm*4
      write(21,*) 'NumConnectedBoundaryFaces=0'
      write(21,*) 'TotalNumBoundaryConnections=0'

      do i = 1, this%smart_pointer%nbnode
        write(21,*) this%smart_pointer%coord_nodes(i)%p%x
      end do

      do i = 1, this%smart_pointer%nbnode
        write(21,*) this%smart_pointer%coord_nodes(i)%p%y
      end do 
      
      ! Node indexes
      do i = 1, this%smart_pointer%nbelm
        write(21,*) this%smart_pointer%id_nodes(i)%pn%id_node(6), this%smart_pointer%id_nodes(i)%pn%id_node(7)
        write(21,*) this%smart_pointer%id_nodes(i)%pn%id_node(7), this%smart_pointer%id_nodes(i)%pn%id_node(8)
	    write(21,*) this%smart_pointer%id_nodes(i)%pn%id_node(8), this%smart_pointer%id_nodes(i)%pn%id_node(9)
	    write(21,*) this%smart_pointer%id_nodes(i)%pn%id_node(9), this%smart_pointer%id_nodes(i)%pn%id_node(6)
      enddo
      
      do i = 1, this%smart_pointer%nbelm
        write(21,*) i,i,i,i
      enddo  
  
      do i = 1, this%smart_pointer%nbelm
        write(21,*) 0,0,0,0
      enddo

      
      close(unit = 21)

      end subroutine
!----------------------------------------------------------------------    
end module

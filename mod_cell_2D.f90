module mod_cell_2D
    use mod_read_gmsh
    use mod_point_util
    implicit none

    type face
        type(point)            :: p1,p2
        integer(4)             :: bc_typ = 0
        integer(4)             :: idface = 0
        real(8)                :: area = 0.0d0
        type(point)            :: centroid
    end type face

    type cell_2D
        integer(4)  :: ident = 0
        type(point) :: vertex(4)

        type(face),dimension(4) :: faces

        type(cell_2D),pointer :: neighbor1 => null()
        type(cell_2D),pointer :: neighbor2 => null()
        type(cell_2D),pointer :: neighbor3 => null()
        type(cell_2D),pointer :: neighbor4 => null()
        
        real(8)     :: vol = 0.0d0

    end type cell_2D
     
    type list_cell_2D
        type(cell_2D),pointer :: p
    end type

    type, extends(obj_gmsh_reader) :: obj_cell_2D
        type(list_cell_2D), dimension(:), allocatable :: cell
    contains
        procedure :: build_cells      => construct_cells
        procedure :: calcul_vol       => calcul_vol_cells
        procedure :: calcul_face_info => calcul_area_cent_faces
        procedure :: assign_face_id   => assign_id_face_v2
        procedure :: find_neighbor    => detect_neighbor
    end type obj_cell_2D
    
    integer(4)      :: nbfaces
    
    private :: sort_vertex

    contains
!----------------------------------------------------------------------    
       subroutine construct_cells(this)
       implicit none

       class(obj_cell_2D), intent(inout) :: this
       integer(4) :: i,j
       integer(4) :: idnode
       
       integer(4) :: idnode1, idnode2
       
       type(cell_2D), pointer        :: pcell
       type(node_ident_msh), pointer :: pmsh

       allocate(this%cell(1:this%nbelm))

       do i = 1, this%nbelm
          allocate(this%cell(i)%p)
       end do

       ! Assign vertexes of a cell

       do i = 1, this%nbelm
         this%cell(i)%p%ident = i
         this%cell(i)%p%vertex(1)%ident = this%id_nodes(i)%pn%id_node(6)
         idnode = this%cell(i)%p%vertex(1)%ident
         do j = 1, this%nbnode
            if (idnode == this%coord_nodes(j)%p%ident) then
                this%cell(i)%p%vertex(1)%x = this%coord_nodes(j)%p%x
                this%cell(i)%p%vertex(1)%y = this%coord_nodes(j)%p%y
            end if
         end do

         this%cell(i)%p%vertex(2)%ident = this%id_nodes(i)%pn%id_node(7)
         idnode = this%cell(i)%p%vertex(2)%ident
         do j = 1, this%nbnode
            if (idnode == this%coord_nodes(j)%p%ident) then
                this%cell(i)%p%vertex(2)%x = this%coord_nodes(j)%p%x
                this%cell(i)%p%vertex(2)%y = this%coord_nodes(j)%p%y
            end if
         end do

         this%cell(i)%p%vertex(3)%ident = this%id_nodes(i)%pn%id_node(8)
         idnode = this%cell(i)%p%vertex(3)%ident
         do j = 1, this%nbnode
            if (idnode == this%coord_nodes(j)%p%ident) then
                this%cell(i)%p%vertex(3)%x = this%coord_nodes(j)%p%x
                this%cell(i)%p%vertex(3)%y = this%coord_nodes(j)%p%y
            end if
         end do

         this%cell(i)%p%vertex(4)%ident = this%id_nodes(i)%pn%id_node(9)
         idnode = this%cell(i)%p%vertex(4)%ident
         do j = 1, this%nbnode
            if (idnode == this%coord_nodes(j)%p%ident) then
                this%cell(i)%p%vertex(4)%x = this%coord_nodes(j)%p%x
                this%cell(i)%p%vertex(4)%y = this%coord_nodes(j)%p%y
            end if
         end do

       end do
       
       ! Sort vertexes for each cell
       ! sort_vertex failed
       !do i = 1, nbelm 
       !    call sort_vertex(cell(i)%p%vertex)
       !enddo 

       ! Assign faces for each cell

       do i = 1, this%nbelm
         this%cell(i)%p%faces(1)%p1%x = this%cell(i)%p%vertex(1)%x
         this%cell(i)%p%faces(1)%p1%y = this%cell(i)%p%vertex(1)%y

         this%cell(i)%p%faces(1)%p2%x = this%cell(i)%p%vertex(2)%x
         this%cell(i)%p%faces(1)%p2%y = this%cell(i)%p%vertex(2)%y

         this%cell(i)%p%faces(2)%p1%x = this%cell(i)%p%vertex(2)%x
         this%cell(i)%p%faces(2)%p1%y = this%cell(i)%p%vertex(2)%y

         this%cell(i)%p%faces(2)%p2%x = this%cell(i)%p%vertex(3)%x
         this%cell(i)%p%faces(2)%p2%y = this%cell(i)%p%vertex(3)%y

         this%cell(i)%p%faces(3)%p1%x = this%cell(i)%p%vertex(3)%x
         this%cell(i)%p%faces(3)%p1%y = this%cell(i)%p%vertex(3)%y

         this%cell(i)%p%faces(3)%p2%x = this%cell(i)%p%vertex(4)%x
         this%cell(i)%p%faces(3)%p2%y = this%cell(i)%p%vertex(4)%y

         this%cell(i)%p%faces(4)%p1%x = this%cell(i)%p%vertex(4)%x
         this%cell(i)%p%faces(4)%p1%y = this%cell(i)%p%vertex(4)%y

         this%cell(i)%p%faces(4)%p2%x = this%cell(i)%p%vertex(1)%x
         this%cell(i)%p%faces(4)%p2%y = this%cell(i)%p%vertex(1)%y

       end do
       
       ! Assign boundary condition for each face of a cell 
       
       do i = 1, this%nbelm
           pcell => this%cell(i)%p 
           idnode1 = pcell%vertex(1)%ident 
           idnode2 = pcell%vertex(2)%ident 
           
           loop_msh1: do j = 1, this%nbel_msh
               pmsh => this%id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(1)%bc_typ = pmsh%tag1
                       cycle loop_msh1
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(1)%bc_typ = pmsh%tag1
                       cycle loop_msh1
                   endif 
               endif
           enddo loop_msh1
           
           idnode1 = pcell%vertex(2)%ident 
           idnode2 = pcell%vertex(3)%ident 
           
           loop_msh2: do j = 1, this%nbel_msh
               pmsh => this%id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(2)%bc_typ = pmsh%tag1
                       cycle loop_msh2
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(2)%bc_typ = pmsh%tag1
                       cycle loop_msh2
                   endif 
               endif
           enddo loop_msh2
           
           idnode1 = pcell%vertex(3)%ident 
           idnode2 = pcell%vertex(4)%ident 
           
           loop_msh3: do j = 1, this%nbel_msh
               pmsh => this%id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(3)%bc_typ = pmsh%tag1
                       cycle loop_msh3
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(3)%bc_typ = pmsh%tag1
                       cycle loop_msh3
                   endif 
               endif
           enddo loop_msh3
           
           idnode1 = pcell%vertex(4)%ident 
           idnode2 = pcell%vertex(1)%ident 
           
           loop_msh4: do j = 1, this%nbel_msh
               pmsh => this%id_nodes_msh(j)%pn 
               if (pmsh%elem_typ == 27) then ! I write a solver for NACA 0012 HIOCFD firstly
                   if (idnode1 == pmsh%id_node(6) .and. idnode2 == pmsh%id_node(7)) then
                       pcell%faces(4)%bc_typ = pmsh%tag1
                       cycle loop_msh4
                   else if (idnode1 == pmsh%id_node(7) .and. idnode2 == pmsh%id_node(6)) then
                       pcell%faces(4)%bc_typ = pmsh%tag1
                       cycle loop_msh4
                   endif 
               endif
           enddo loop_msh4
           
       enddo
       
       end subroutine
!----------------------------------------------------------------------  
       subroutine calcul_vol_cells(this)
       implicit none 
       class(obj_cell_2D), intent(inout) :: this
       integer                           :: i
       type(cell_2D), pointer            :: pc
       real(8)                           :: vol, x1, x2, x3, x4, y1, y2, y3, y4
       
       do i = 1, this%nbelm
           pc => this%cell(i)%p
           
           x1     = pc%vertex(1)%x
           x2     = pc%vertex(2)%x
           x3     = pc%vertex(3)%x
           x4     = pc%vertex(4)%x
           
           y1     = pc%vertex(1)%y
           y2     = pc%vertex(2)%y
           y3     = pc%vertex(3)%y
           y4     = pc%vertex(4)%y
           
           vol    = 0.5d0 * ( (x1 - x3)*(y2 - y4) + (x4 - x2)*(y1 - y3) )
           this%cell(i)%p%vol = vol
       enddo 
       end subroutine calcul_vol_cells
!----------------------------------------------------------------------
       subroutine calcul_area_cent_faces(this)
       use mod_vector_algebra
       implicit none 
       class(obj_cell_2D), intent(inout) :: this
       integer                           :: i, j
       type(cell_2D), pointer            :: pcel
       type(face),    pointer            :: pfac
       type(vector_t)                    :: vec
       
       do i = 1, this%nbelm
           pcel => this%cell(i)%p
           do j = 1, 4
               pfac => pcel%faces(j)
               vec                                = vector_t( pfac%p1%x - pfac%p2%x, pfac%p1%y - pfac%p2%y, 0.0d0 )
               this%cell(i)%p%faces(j)%area       = .abs.(vec)
               this%cell(i)%p%faces(j)%centroid%x =  0.5d0 * (pfac%p1%x + pfac%p2%x)
               this%cell(i)%p%faces(j)%centroid%y =  0.5d0 * (pfac%p1%y + pfac%p2%y)
           enddo 
       enddo 
       end subroutine calcul_area_cent_faces
!----------------------------------------------------------------------  
       !subroutine assign_id_face
       !use mod_point_util
       !implicit none
       !integer                :: i, j, k
       !type(cell_2D), pointer :: pc, pc2
       !type(face),    pointer :: pfac, pfac2
       !
       !nbfaces = 0
       !do i = 1, nbelm
       !    pc => cell(i)%p
       !        pfac => pc%faces(1)
       !        if (pfac%idface == 0) then ! a face that has not been indexed yet
       !        if (associated(pc%neighbor1)) then
       !           pc2 => pc%neighbor1
       !           if (associated(pc%neighbor1%neighbor3, pc)) then 
       !             if (pfac%bc_typ == 0) then ! normal face
       !                 pfac2 => pc2%faces(3)
       !                 if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
       !                     if (pfac2%idface == 0) then
       !                         nbfaces = nbfaces + 1
       !                         pfac%idface = nbfaces
       !                         pfac2%idface = nbfaces
       !                     endif
       !                 endif
       !             endif
       !           endif 
       !         else if (.not. associated(pc%neighbor1)) then 
       !            if (pfac%bc_typ /= 0) then ! boundary face
       !                 nbfaces = nbfaces + 1
       !                 pfac%idface = nbfaces 
       !            endif
       !         endif
       !        endif 
       !        
       !        pfac => pc%faces(2)
       !        if (pfac%idface == 0) then ! a face that has not been indexed yet
       !        if (associated(pc%neighbor2)) then
       !           pc2 => pc%neighbor2
       !           if (associated(pc%neighbor2%neighbor4, pc)) then 
       !             if (pfac%bc_typ == 0) then ! normal face
       !                 pfac2 => pc2%faces(4)
       !                 if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
       !                     if (pfac2%idface == 0) then
       !                         nbfaces = nbfaces + 1
       !                         pfac%idface = nbfaces
       !                         pfac2%idface = nbfaces
       !                     endif
       !                 endif
       !             endif
       !           endif 
       !         else if (.not. associated(pc%neighbor2)) then 
       !            if (pfac%bc_typ /= 0) then ! boundary face
       !                 nbfaces = nbfaces + 1
       !                 pfac%idface = nbfaces 
       !            endif
       !         endif
       !        endif
       !        
       !        pfac => pc%faces(3)
       !        if (pfac%idface == 0) then ! a face that has not been indexed yet
       !        if (associated(pc%neighbor3)) then
       !           pc2 => pc%neighbor3
       !           if (associated(pc%neighbor3%neighbor1, pc)) then 
       !             if (pfac%bc_typ == 0) then ! normal face
       !                 pfac2 => pc2%faces(1)
       !                 if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
       !                     if (pfac2%idface == 0) then
       !                         nbfaces = nbfaces + 1
       !                         pfac%idface = nbfaces
       !                         pfac2%idface = nbfaces
       !                     endif
       !                 endif
       !             endif
       !           endif 
       !         else if (.not. associated(pc%neighbor3)) then 
       !            if (pfac%bc_typ /= 0) then ! boundary face
       !                 nbfaces = nbfaces + 1
       !                 pfac%idface = nbfaces 
       !            endif
       !         endif
       !        endif
       !        
       !        pfac => pc%faces(4)
       !        if (pfac%idface == 0) then ! a face that has not been indexed yet
       !        if (associated(pc%neighbor4)) then
       !           pc2 => pc%neighbor4
       !           if (associated(pc%neighbor4%neighbor2, pc)) then 
       !             if (pfac%bc_typ == 0) then ! normal face
       !                 pfac2 => pc2%faces(2)
       !                 if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
       !                     if (pfac2%idface == 0) then
       !                         nbfaces = nbfaces + 1
       !                         pfac%idface = nbfaces
       !                         pfac2%idface = nbfaces
       !                     endif
       !                 endif
       !             endif
       !           endif 
       !         else if (.not. associated(pc%neighbor4)) then 
       !            if (pfac%bc_typ /= 0) then ! boundary face
       !                 nbfaces = nbfaces + 1
       !                 pfac%idface = nbfaces 
       !            endif
       !         endif
       !         endif
       !enddo 
       !end subroutine assign_id_face
!----------------------------------------------------------------------  
       subroutine sort_vertex(vertices)
       !Reference: http://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90                  
       implicit none 
       
       type(point) :: vertices(:)
       type(point) :: vertices_temp(4)
       
       integer     :: i, location
       
       do i = 1, size(vertices) - 1
           location = findminimum(i,size(vertices))
           call swap_point(vertices(i),vertices(location))
       enddo 
       
       vertices_temp = vertices
       
       if (vertices_temp(1)%y < vertices_temp(2)%y) then
           vertices(1) = vertices_temp(1)
           vertices(4) = vertices_temp(2)
       else if (vertices_temp(1)%y >= vertices_temp(2)%y) then 
           vertices(1) = vertices_temp(2)
           vertices(4) = vertices_temp(1)
       endif
       
       if (vertices_temp(3)%y < vertices_temp(4)%y) then
           vertices(2) = vertices_temp(3)
           vertices(3) = vertices_temp(4)
       else if (vertices_temp(3)%y >= vertices_temp(4)%y) then 
           vertices(2) = vertices_temp(4)
           vertices(3) = vertices_temp(3)
       endif 
       
       contains
          integer function findminimum(start,endp)
          implicit none
          integer, intent(in) :: start, endp
          
          integer             :: location, i
          type(point)         :: minimum
          
          minimum = vertices(start)
          location = start
          do i = start + 1, endp
              if (vertices(i)%x < minimum%x) then
                  minimum = vertices(i)
                  location = i
              endif 
          enddo 
          findminimum = location
          
          end function findminimum
       end subroutine sort_vertex
!----------------------------------------------------------------------
       subroutine assign_id_face_v2(this)
       use mod_point_util
       implicit none
       class(obj_cell_2D), intent(inout) :: this
       integer                           :: i, j, k
       type(cell_2D), pointer            :: pc, pc2
       type(face),    pointer            :: pfac, pfac2
       
       nbfaces = 0
       do i = 1, this%nbelm
           pc => this%cell(i)%p
           do j = 1, 4
               pfac => pc%faces(j)
               if (pfac%idface == 0) then ! a face that has not been indexed yet
                   pc2 => pc%neighbor1
                   if (associated(pc%neighbor1)) then
                       do k = 1, 4
                           pfac2 => pc2%faces(k)
                           if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                this%cell(i)%p%faces(j)%idface           = nbfaces
                                this%cell(i)%p%neighbor1%faces(k)%idface = nbfaces
                            endif
                           endif
                       enddo
                   else if (.not. associated(pc%neighbor1)) then 
                       if (pfac%bc_typ /= 0) then ! boundary face
                           if (pfac%idface > 0) cycle 
                        nbfaces = nbfaces + 1
                        this%cell(i)%p%faces(j)%idface = nbfaces 
                       endif
                   endif 
                   pc2 => pc%neighbor2
                   if (associated(pc%neighbor2)) then
                       do k = 1, 4
                           pfac2 => pc2%faces(k)
                           if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                this%cell(i)%p%faces(j)%idface           = nbfaces
                                this%cell(i)%p%neighbor2%faces(k)%idface = nbfaces
                            endif
                           endif
                       enddo 
                   else if (.not. associated(pc%neighbor2)) then 
                       if (pfac%bc_typ /= 0) then ! boundary face
                           if (pfac%idface > 0) cycle
                        nbfaces = nbfaces + 1
                        this%cell(i)%p%faces(j)%idface = nbfaces 
                       endif
                   endif 
                   pc2 => pc%neighbor3
                   if (associated(pc%neighbor3)) then
                       do k = 1, 4
                           pfac2 => pc2%faces(k)
                           if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                this%cell(i)%p%faces(j)%idface           = nbfaces
                                this%cell(i)%p%neighbor3%faces(k)%idface = nbfaces
                            endif
                           endif
                       enddo
                   else if (.not. associated(pc%neighbor3)) then 
                       if (pfac%bc_typ /= 0) then ! boundary face
                           if (pfac%idface > 0) cycle
                        nbfaces = nbfaces + 1
                        this%cell(i)%p%faces(j)%idface = nbfaces 
                       endif
                   endif 
                   pc2 => pc%neighbor4
                   if (associated(pc%neighbor4)) then
                       do k = 1, 4
                           pfac2 => pc2%faces(k)
                           if ((pfac%centroid .eq. pfac2%centroid) .and. abs(pfac%area - pfac2%area) <= 1.0d-9) then
                            if (pfac2%idface == 0) then
                                nbfaces = nbfaces + 1
                                this%cell(i)%p%faces(j)%idface           = nbfaces
                                this%cell(i)%p%neighbor4%faces(k)%idface = nbfaces
                            endif
                           endif
                       enddo 
                   else if (.not. associated(pc%neighbor4)) then 
                       if (pfac%bc_typ /= 0) then ! boundary face
                           if (pfac%idface > 0) cycle
                        nbfaces = nbfaces + 1
                        this%cell(i)%p%faces(j)%idface = nbfaces 
                       endif
                   endif 
               endif 
           enddo 
       enddo 
       
       !begin debug
       !do i = 1, nbelm
       !    pc => cell(i)%p
       !    do j = 1, 4
       !        pfac => pc%faces(j)
       !        if (pfac%idface == 4409) then
       !            print*,'check cell', i
       !        endif 
       !    enddo
       !enddo 
       !end debug 
       end subroutine assign_id_face_v2
!----------------------------------------------------------------------  
    subroutine detect_neighbor(this)
      implicit none
      class(obj_cell_2D), intent(inout) :: this
      integer(4)                        :: i,j,k,t
      integer(4)                        :: idnode1,idnode2
      integer(4)                        :: cnt
      type(cell_2D),pointer             :: pc => null()
      type(cell_2D),pointer             :: pr => null()

!      do i = 1, nbelm
!        do k = 1, 4
!            write(*,*) cell(i)%p%vertex(k)%ident
!        end do
!        write(*,*)
!      end do

      do i = 1, this%nbelm
          idnode1 = this%cell(i)%p%vertex(1)%ident
          idnode2 = this%cell(i)%p%vertex(2)%ident
          pc => this%cell(i)%p

          do j = 1, this%nbelm
            pr => this%cell(j)%p
            cnt = 0
            if (associated(pc,pr)) cycle
            do t = 1, 4
              if (idnode1 == pr%vertex(t)%ident) cnt = cnt + 1
              if (idnode2 == pr%vertex(t)%ident) cnt = cnt + 1
            end do
            if (cnt == 2) then
              this%cell(i)%p%neighbor1 => pr
            end if
          end do

       end do

       do i = 1, this%nbelm

          idnode1 = this%cell(i)%p%vertex(2)%ident
          idnode2 = this%cell(i)%p%vertex(3)%ident
          pc => this%cell(i)%p

          do j = 1, this%nbelm
            pr => this%cell(j)%p
            cnt = 0
            if (associated(pc,pr)) cycle
            do t = 1, 4
              if (idnode1 == pr%vertex(t)%ident) cnt = cnt + 1
              if (idnode2 == pr%vertex(t)%ident) cnt = cnt + 1
            end do
            if (cnt == 2) then
              this%cell(i)%p%neighbor2 => pr
            end if
           end do

      end do

      do i = 1, this%nbelm

          idnode1 = this%cell(i)%p%vertex(3)%ident
          idnode2 = this%cell(i)%p%vertex(4)%ident
          pc => this%cell(i)%p

          do j = 1, this%nbelm
            pr => this%cell(j)%p
            cnt = 0
            if (associated(pc,pr)) cycle
            do t = 1, 4
              if (idnode1 == pr%vertex(t)%ident) cnt = cnt + 1
              if (idnode2 == pr%vertex(t)%ident) cnt = cnt + 1
            end do
            if (cnt == 2) then
              this%cell(i)%p%neighbor3 => pr
            end if
           end do

       end do

       do i = 1, this%nbelm


           idnode1 = this%cell(i)%p%vertex(4)%ident
           idnode2 = this%cell(i)%p%vertex(1)%ident

           pc => this%cell(i)%p
           do j = 1, this%nbelm
             pr => this%cell(j)%p
             cnt = 0
             if (associated(pc,pr)) cycle
             do t = 1, 4
               if (idnode1 == pr%vertex(t)%ident) cnt = cnt + 1
               if (idnode2 == pr%vertex(t)%ident) cnt = cnt + 1
             end do
             if (cnt == 2) then
               this%cell(i)%p%neighbor4 => pr
             end if
            end do

      end do
      end subroutine

end module

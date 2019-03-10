module mod_read_gmsh

    type point
        integer(4) :: ident
        real(8)    :: x = 0.0d0
        real(8)    :: y = 0.0d0
        real(8)    :: z = 0.0d0
    end type

    type list_point
        type(point),pointer :: p => null()
    end type

    type node_ident
        integer(4),dimension(:),allocatable:: id_node
    end type

    type,extends(node_ident) :: node_ident_msh
        integer(4) :: ident
        integer(4) :: elem_typ
        integer(4) :: nb_tags
        integer(4) :: tag1
        integer(4) :: tag2
    end type

    type ptr_node_ident
        type(node_ident),pointer :: pn => null()
    end type

    type ptr_node_ident_msh
        type(node_ident_msh),pointer :: pn => null()
    end type

    type obj_gmsh_reader

        integer(4) :: nbnode
        integer(4) :: nbel_msh
        integer(4) :: nbelm
        type(ptr_node_ident), dimension(:), allocatable     :: id_nodes
        type(ptr_node_ident_msh), dimension(:), allocatable :: id_nodes_msh
        type(list_point), dimension(:), allocatable         :: coord_nodes

    contains
        procedure :: read_msh => read_mesh
        procedure :: assign_id_nodes => construct_id_nodes
    end type obj_gmsh_reader

    character(300), save :: fname
    real(8),        save :: rho_init, ux_init, uy_init, t_init

    contains
!----------------------------------------------------------------------
       subroutine read_mesh(this)
       implicit none

       namelist /mesh_file/ fname
       namelist /init_valu/ rho_init, ux_init, uy_init, t_init

       class(obj_gmsh_reader), intent(inout) :: this
       character(300) :: finput
       character(30) :: cline
       integer(4) :: i
       integer(4) :: elem_typ

       open(unit = 100, file = 'input.dat')
       read(100, nml = mesh_file)
       read(100, nml = init_valu)
       close(100)

       finput = trim(fname)//'.msh'

       write(6,'(a)') 'Input mesh file: ',finput

       open(unit = 20, file = finput, status = 'old')
       read(20,*)
       read(20,*)
       read(20,*)
       read(20,*)
       read(20,*) this%nbnode

       allocate(this%coord_nodes(1:this%nbnode))
       do i = 1, this%nbnode
         allocate(this%coord_nodes(i)%p)
       end do

       do i = 1, this%nbnode
        read(20,*) this%coord_nodes(i)%p%ident, this%coord_nodes(i)%p%x, this%coord_nodes(i)%p%y, this%coord_nodes(i)%p%z
       end do



       read(20,*)
       read(20,*)

       read(20,*) this%nbel_msh

       allocate(this%id_nodes_msh(1:this%nbel_msh))

       do i = 1, this%nbel_msh
        allocate(this%id_nodes_msh(i)%pn)
       end do

       this%nbelm = 0

       do i = 1, this%nbel_msh
        read(20,*) this%id_nodes_msh(i)%pn%ident, this%id_nodes_msh(i)%pn%elem_typ, this%id_nodes_msh(i)%pn%nb_tags, &
            & this%id_nodes_msh(i)%pn%tag1, this%id_nodes_msh(i)%pn%tag2
        elem_typ = this%id_nodes_msh(i)%pn%elem_typ
        select case (elem_typ)
            case (1) ! 2-node line.
                allocate(this%id_nodes_msh(i)%pn%id_node(1:7))
            case (3) ! 4-node quadrangle.
                allocate(this%id_nodes_msh(i)%pn%id_node(1:9))
                this%nbelm = this%nbelm + 1
            case (15) ! 1-node point.
                allocate(this%id_nodes_msh(i)%pn%id_node(1:6))
            case (37) ! 5-node edge quadrangle.
                allocate(this%id_nodes_msh(i)%pn%id_node(1:30))
                this%nbelm = this%nbelm + 1
            case (27) ! boundary 5-node edge.
                allocate(this%id_nodes_msh(i)%pn%id_node(1:10))
            case default
                print*,'comming soon !'
                stop
        end select
       end do

       close(unit = 20)

       open(unit = 20, file = finput, status = 'old')
       do
           read(20,*) cline
           cline = trim(cline)
           if (cline == '$Elements') exit
       enddo

       read(20,*)
       do i = 1, this%nbel_msh
           read(20,*) this%id_nodes_msh(i)%pn%id_node(:)
       enddo

       close(unit = 20)

       end subroutine
!----------------------------------------------------------------------
       subroutine construct_id_nodes(this)
       implicit none
       class(obj_gmsh_reader), intent(inout) :: this
       integer(4) :: i,j

       allocate(this%id_nodes(1:this%nbelm))
       do i = 1, this%nbelm
           allocate(this%id_nodes(i)%pn)

           if (this%id_nodes_msh(i)%pn%elem_typ == 3) then
                allocate(this%id_nodes(i)%pn%id_node(1:9))
           endif

           if (this%id_nodes_msh(i)%pn%elem_typ == 37) then
                allocate(this%id_nodes(i)%pn%id_node(1:30))
           endif

       enddo

       j = 0
       do i = 1, this%nbel_msh
           if (this%id_nodes_msh(i)%pn%elem_typ == 3) then
               j = j + 1
               this%id_nodes(j)%pn%id_node(:) = this%id_nodes_msh(i)%pn%id_node(:)
           endif

           if (this%id_nodes_msh(i)%pn%elem_typ == 37) then
               j = j + 1
               this%id_nodes(j)%pn%id_node(:) = this%id_nodes_msh(i)%pn%id_node(:)
           endif
       enddo

       end subroutine
!----------------------------------------------------------------------
end module

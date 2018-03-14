module mod_fvm_face_2D
    use mod_cell_2D
    use mod_vector_algebra
    implicit none 
    
    type, extends(face) :: fvm_face_2D
        type(vector_t)         :: normal_vector
        type(vector_t)         :: tangent_vector
        type(cell_2D), pointer :: left_cell => null()
        type(cell_2D), pointer :: right_cell => null()
        real(8)                :: len_nor 
        real(8)                :: len_tan
    end type
    
    type :: pt_fvm_face_2D 
        type(fvm_face_2D), pointer :: f => null()
    end type 
    
    type :: obj_fvm_face_2D
        type(pt_fvm_face_2D), dimension(:), allocatable :: face_2D
        class(obj_cell_2D), pointer                     :: smart_pointer
    contains
        procedure :: alloc => allocate_face_2D
        procedure :: get   => assign_face_2D
    end type
    
    contains 
!----------------------------------------------------------------------
        subroutine allocate_face_2D(this)
        implicit none 
        
        class(obj_fvm_face_2D), intent(inout) :: this
        integer :: i
        
        if (.not. allocated(this%face_2D)) then
            allocate(this%face_2D(1:nbfaces))
        endif 
        
        do i = 1, nbfaces
            allocate(this%face_2D(i)%f)
        enddo 
        end subroutine allocate_face_2D
!----------------------------------------------------------------------
        subroutine assign_face_2D(this)
        implicit none 
        
        class(obj_fvm_face_2D), intent(inout) :: this
        integer                    :: icel, ifac, jfac, jfac_nbor
        type(cell_2D),     pointer :: pcel, pcel_nbor
        type(fvm_face_2D), pointer :: pfac
        type(face),        pointer :: pfac2, pfac2_nbor
        real(8)                    :: norm_vec
        
        do ifac = 1, nbfaces
            
            !if (ifac == 4409) then
            !    print*,'4409'
            !endif 
            
            pfac => this%face_2D(ifac)%f
            do icel = 1, this%smart_pointer%nbelm
                pcel => this%smart_pointer%cell(icel)%p
                do jfac = 1, 4
                    pfac2 => pcel%faces(jfac)
                    
                    !if (pfac2%idface == 4409) then
                    !    print*,'4409 2' 
                    !endif 
                    
                    if (pfac2%idface == ifac) then
                        
                        if (this%face_2D(ifac)%f%idface > 0) cycle ! just to index face_2D one time
                        
                        ! Copy data from face to fvm_face_2D
                        this%face_2D(ifac)%f%p1       = this%smart_pointer%cell(icel)%p%faces(jfac)%p1
                        this%face_2D(ifac)%f%p2       = this%smart_pointer%cell(icel)%p%faces(jfac)%p2
                        this%face_2D(ifac)%f%bc_typ   = this%smart_pointer%cell(icel)%p%faces(jfac)%bc_typ
                        this%face_2D(ifac)%f%idface   = this%smart_pointer%cell(icel)%p%faces(jfac)%idface
                        this%face_2D(ifac)%f%area     = this%smart_pointer%cell(icel)%p%faces(jfac)%area
                        this%face_2D(ifac)%f%centroid = this%smart_pointer%cell(icel)%p%faces(jfac)%centroid
                        
                        ! Assign left cell and right cell for each face_2D
                        pcel_nbor => pcel%neighbor1
                        if (associated(pcel_nbor)) then
                            do jfac_nbor = 1, 4
                                pfac2_nbor => pcel_nbor%faces(jfac_nbor)
                                if (pfac2_nbor%idface == ifac) then
                                    if (pcel%ident <= pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p%neighbor1
                                    else if (pcel%ident > pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p%neighbor1
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p
                                    endif 
                                endif
                            enddo
                        else if (.not. associated(pcel_nbor)) then 
                            if (.not. associated(this%face_2D(ifac)%f%left_cell)) then
                                this%face_2D(ifac)%f%left_cell => this%smart_pointer%cell(icel)%p
                            endif 
                        endif 
                        
                        pcel_nbor => pcel%neighbor2
                        if (associated(pcel_nbor)) then
                            do jfac_nbor = 1, 4
                                pfac2_nbor => pcel_nbor%faces(jfac_nbor)
                                if (pfac2_nbor%idface == ifac) then
                                    if (pcel%ident <= pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p%neighbor2
                                    else if (pcel%ident > pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p%neighbor2
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p
                                    endif 
                                endif
                            enddo 
                        else if (.not. associated(pcel_nbor)) then 
                            if (.not. associated(this%face_2D(ifac)%f%left_cell)) then
                                this%face_2D(ifac)%f%left_cell => this%smart_pointer%cell(icel)%p
                            endif
                        endif 
                        
                        pcel_nbor => pcel%neighbor3
                        if (associated(pcel_nbor)) then
                            do jfac_nbor = 1, 4
                                pfac2_nbor => pcel_nbor%faces(jfac_nbor)
                                if (pfac2_nbor%idface == ifac) then
                                    if (pcel%ident <= pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p%neighbor3
                                    else if (pcel%ident > pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p%neighbor3
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p
                                    endif 
                                endif
                            enddo
                        else if (.not. associated(pcel_nbor)) then 
                            if (.not. associated(this%face_2D(ifac)%f%left_cell)) then
                                this%face_2D(ifac)%f%left_cell => this%smart_pointer%cell(icel)%p
                            endif
                        endif 
                        
                        pcel_nbor => pcel%neighbor4
                        if (associated(pcel_nbor)) then
                            do jfac_nbor = 1, 4
                                pfac2_nbor => pcel_nbor%faces(jfac_nbor)
                                if (pfac2_nbor%idface == ifac) then
                                    if (pcel%ident <= pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p%neighbor4
                                    else if (pcel%ident > pcel_nbor%ident) then
                                        this%face_2D(ifac)%f%left_cell  => this%smart_pointer%cell(icel)%p%neighbor4
                                        this%face_2D(ifac)%f%right_cell => this%smart_pointer%cell(icel)%p
                                    endif 
                                endif
                            enddo
                        else if (.not. associated(pcel_nbor)) then 
                            if (.not. associated(this%face_2D(ifac)%f%left_cell)) then
                                this%face_2D(ifac)%f%left_cell => this%smart_pointer%cell(icel)%p
                            endif
                        endif 
                        
                        ! Assign tangent and normal vectors
                        this%face_2D(ifac)%f%tangent_vector = vector_t(this%face_2D(ifac)%f%p2%x - this%face_2D(ifac)%f%p1%x, &
                            & this%face_2D(ifac)%f%p2%y - this%face_2D(ifac)%f%p1%y, 0.0d0)
                        norm_vec                       = .abs. (this%face_2D(ifac)%f%tangent_vector)
                        this%face_2D(ifac)%f%len_tan        = norm_vec
                        this%face_2D(ifac)%f%tangent_vector = this%face_2D(ifac)%f%tangent_vector .times. (1.0d0 / norm_vec)
                        
                        if (associated(this%face_2D(ifac)%f%right_cell, this%smart_pointer%cell(icel)%p)) then
                            this%face_2D(ifac)%f%normal_vector  = vector_t(this%face_2D(ifac)%f%p1%y - this%face_2D(ifac)%f%p2%y, &
                                & this%face_2D(ifac)%f%p2%x - this%face_2D(ifac)%f%p1%x, 0.0d0)
                        else if (associated(this%face_2D(ifac)%f%left_cell, this%smart_pointer%cell(icel)%p)) then 
                            this%face_2D(ifac)%f%normal_vector =  vector_t(-this%face_2D(ifac)%f%p1%y + this%face_2D(ifac)%f%p2%y, &
                                & -this%face_2D(ifac)%f%p2%x + this%face_2D(ifac)%f%p1%x, 0.0d0)
                        endif
                        norm_vec                       = .abs. (this%face_2D(ifac)%f%normal_vector)
                        this%face_2D(ifac)%f%len_nor        = norm_vec
                        this%face_2D(ifac)%f%normal_vector  = this%face_2D(ifac)%f%normal_vector .times. (1.0d0 / norm_vec) 
                    endif 
                enddo
            enddo
        enddo 
        end subroutine assign_face_2D
end module mod_fvm_face_2D
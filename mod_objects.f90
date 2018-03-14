module mod_objects
    use mod_read_gmsh
    use mod_write_vtk
    use mod_cell_2D
    use mod_fvm_face_2D
    implicit none
    
    type(obj_cell_2D), target :: list_cell
    type(obj_msh_vis)         :: visualizer
    type(obj_fvm_face_2D)     :: faces_fvm
    
end module
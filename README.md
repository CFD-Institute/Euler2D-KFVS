# euler2D-kfvs-Fortran2003

Euler2D-kfvs-Fortran2003 is a compressible Euler equations solver using kinetic flux vector splitting (KFVS) scheme. It is the successor of my Fortran 95 code euler-2D-kfvs

Some features:

- Cell centered finite volume method.

- Unstructured quadrilateral meshes.

- Kinetic flux vector splitting (KFVS) scheme.

- Written in Modern Fortran (Fortran 2003). Object-oriented programming has been implemented, mainly in the mesh reader. The solver is kept in Fortran 95 now, to retain its execution speed, which is the major advantage of the Fortran programming language. For further development, maybe I will upgrade the solver to Fortran 2003, in order to integrate new numerical schemes.

Euler2D-kfvs-Fortran2003 is inspired by the structured finite volume solver of [Prof. Luc Mieussens](https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK3/SIMULATION_NUMERIQUE_ECOULEMENTS_FLUIDES/simulations.html)

[More details](https://github.com/truongd8593/euler2D-kfvs-Fortran2003/wiki)

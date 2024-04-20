#gfortran -c mod_solid.f90 -o mod_solid.o
#gfortran -c mod_io.f90 -o mod_io.o
#gfortran -c mod_solid.o fem2d.f90 -o fem2d.o
#gfortran -c -fPIC mod_io.o fem2d.o soft_particles.f90 -o soft_particles.o
#gfortran -shared -o libsoftparticles.so soft_particles.o

# Compile individual modules
gfortran -c mod_solid.f90 -o mod_solid.o
gfortran -c mod_krod.f90 -o mod_krod.o
gfortran -c mod_io.f90 -o mod_io.o

# Compile fem2d.f90 and mod_solid.f90 into object files
gfortran -c mod_cilia.f90 -o mod_cilia.o

# Compile soft_particles.f90 and link with object files to create shared library
#gfortran -c -fPIC soft_cilia.f90 -o soft_cilia.o
#gfortran -shared -o libsoftcilia.so soft_cilia.o mod_cilia.o mod_solid.o mod_io.o

#gfortran ibmc.f90 -L. -lsoftcilia -o ibmc
gfortran check_mod.f90 mod_io.o mod_cilia.o mod_solid.o mod_krod.o -o check_mod

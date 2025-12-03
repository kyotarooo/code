### OPTIMIZATION ###
OPT	= -O3
#-- Enbale OpenMP for gcc compiler
#OPT	+= -fopenmp
#-- Enable OpenMP for Intel compiler
#OPT	+= -qopenmp

CFLAGS	= -MMD -MP
GFLAGS	= -Wall -g

### DEFINITIONS ###
#-- If the simulation is in a finite volume, and dislocations will path through the free surface
#DFLAGS	+= -D_FREE_SURFACE_
#DFLAGS	+= -D_NOTERMINATION_
#-- If the simulation is not finite volume
DFLAGS	+= -D_NOFEM_
#-- If _PARDISO_ is defined, PARDISO is enabled
#-- If not, MUMPS solver will be selected
#DFLAGS	+= -D_PARDISO_
#DFLAGS += -D_OUTPUT_VTU_

### LIBRARY ###
LFLAGS	= -lm
#-- Library links for PARDISO solver
#LFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#-- Library links for MUMPS solver
#LFLAGS	= -lmumps_common -ldmumps -Wl,-no_compact_unwind -lm
#INCLUDEDIR	= -I/opt/homebrew/opt/brewsci-mumps/include/
#LIBDIR	= -L/opt/homebrew/opt/brewsci-mumps/lib/

### COMPILER ###
CC	= gcc
#-- Intel compiler 
#CC	= icc
#-- If you use MUMPS solver, mpirun is required
#CC	= mpicc

PROGRAM	= ebpdd
OBJS	= dd_annihilation.o\
	dd_constraint.o\
	dd_direct_interaction.o\
	dd_dislocation_core.o\
	dd_dislocation_stress.o\
	dd_element_handle.o\
	dd_grid_stress.o\
	dd_junction.o\
	dd_node_handle.o\
	dd_periodic_bc.o\
	dd_pk_force.o\
	dd_plot.o\
	dd_read.o\
  dd_terminate.o\
	dd_rearrange.o\
	dd_restart.o\
	dd_shape.o\
	dd_solve.o\
	dd_time_step.o\
	dd_write.o\
	dd_fem_correction.o\
	dd_fem_elastic_interaction.o\
	dd_fem_main.o\
	dd_fem_mechanical_behavior.o\
	dd_fem_rk_integral.o\
	dd_fem_trapezoidal_integral.o\
	dd_fem_velocity.o\
	dd_fem_write.o\
	fem_b_matrix.o\
	fem_d_matrix.o\
	fem_read.o\
	fem_shape.o\
	fem_solve.o\
	fem_stiffness.o\
	fem_stress.o\
	fem_surface.o\
	fem_util.o\
	fem_write.o
.c.o	:
	$(CC) $(OPT) $(CFLAGS) $(GFLAGS) $(DFLAGS) -c $< $(INCLUDEDIR)
all	: $(PROGRAM)
$(PROGRAM)	: $(OBJS)
	$(CC) $(OPT) $(GFLAGS) $(DFLAGS) -o $@ $^ $(INCLUDEDIR) $(LIBDIR) $(LFLAGS)
clean	:
	rm -f $(OBJS) $(PROGRAM) ./*~ ./*.d
-include *.d

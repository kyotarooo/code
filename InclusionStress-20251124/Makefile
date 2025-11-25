OPT	= -O3
CFLAGS	= -MMD -MP
GFLAGS	= -Wall -g
LFLAGS	= -lm
CC	= gcc
PROGRAM	= incstrgen
OBJS	= incstr_axial_cylinder.o\
	incstr_bessel.o\
	incstr_cuboid.o\
	incstr_cylinder.o\
	incstr_dilatation_center.o\
	incstr_elastic_dipole.o\
	incstr_elliptic.o\
	incstr_grid.o\
	incstr_legendre.o\
	incstr_main.o\
	incstr_math.o\
	incstr_periodic.o\
	incstr_read.o\
	incstr_sphere.o\
	incstr_stress.o\
	incstr_truncated_sphere.o\
	incstr_write.o
.c.o	:
	$(CC) $(OPT) $(DFLAGS) $(CFLAGS) $(GFLAGS) -c $<
all	: $(PROGRAM)
$(PROGRAM)	: $(OBJS)
	$(CC) $(OPT) $(GFLAGS) -o $@ $(OBJS) $(LFLAGS)
clean	:
	rm -f $(OBJS) $(PROGRAM) ./*~ ./*.d
-include *.d

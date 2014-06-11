#cpl = ifort
#lflags_precise = -fp-model precise -fp-model source -fimf-arch-consistency=true
#lflags_debug = -debug -save-temps -fpic -heap-arrays -O2 -g -traceback -check all -fp-stack-check -vec-report5 $(lflags_precise)
#lflags_fast = -fast
#lflags_profiling = $(lflags_fast) -profile-functions -profile-loops=all -profile-loops-report=2

cpl	= gfortran
lflags_precise = -g -Wall -fcheck=all -ffpe-trap=zero,invalid,overflow,underflow -pedantic-errors
lflags_debug = $(lflags_precise)
lflags_fast = -O3
lflags_profiling = $(lflags_fast) -pg

OBJS = chemistry.o configure.o disk.o heating_cooling.o load_Neufeld_cooling_CO.o data_struct.o grid.o  \
		load_Neufeld_cooling_H2.o load_Neufeld_cooling_H2O.o main.o opkda1.o opkda2.o \
		opkdmain.o sub_global_variables.o sub_trivials.o load_Bethell_Xray.o load_Visser_CO_selfshielding.o \
		interpolation.o statistic_equilibrium.o quick_sort.o \
		montecarlo.o load_Draine_dusts.o ray_tracing.o voigt.o data_dump.o \
        hitran.o lamda.o cdms.o ray_propagating.o vertical_structure.o

LIBS = -lm -lcfitsio

exe_name   = rac
exe_alt    = rac_alt
debug_name = rac_debug
prof_name  = rac_prof

all: lflags = $(lflags_fast)
all: $(exe_name)

alt: lflags = $(lflags_fast)
alt: $(exe_alt)

precise: lflags = $(lflags_fast) $(lflags_precise)
precise: $(exe_name)

debug: lflags = $(lflags_debug)
debug: $(debug_name)

prof: lflags = $(lflags_profiling)
prof: $(prof_name)

cflags = $(lflags) -c

$(exe_name):  $(OBJS)
	$(cpl) $(lflags) -o $(exe_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(debug_name):  $(OBJS)
	$(cpl) $(lflags) -o $(debug_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(prof_name):  $(OBJS)
	$(cpl) $(lflags) -o $(prof_name) $(OBJS) $(LIBS)
	echo -e '\a'

$(exe_alt):  $(OBJS)
	$(cpl) $(lflags) -o $(exe_alt) $(OBJS) $(LIBS)
	echo -e '\a'

main.o: main.f90 configure.o disk.o sub_trivials.o montecarlo.o ray_tracing.o hitran.o lamda.o cdms.o ray_propagating.o
	$(cpl) $(cflags) main.f90

disk.o: disk.f90 chemistry.o heating_cooling.o load_Visser_CO_selfshielding.o data_struct.o grid.o montecarlo.o data_dump.o load_Draine_dusts.o ray_propagating.o vertical_structure.o
	$(cpl) $(cflags) disk.f90

vertical_structure.o: vertical_structure.f90 data_struct.o sub_trivials.o grid.o
	$(cpl) $(cflags) vertical_structure.f90

data_dump.o: data_dump.f90 data_struct.o sub_trivials.o grid.o
	$(cpl) $(cflags) data_dump.f90

hitran.o: hitran.f90 data_struct.o sub_trivials.o sub_global_variables.o quick_sort.o
	$(cpl) $(cflags) hitran.f90

lamda.o: lamda.f90 data_struct.o sub_trivials.o sub_global_variables.o quick_sort.o
	$(cpl) $(cflags) lamda.f90

cdms.o: cdms.f90 data_struct.o sub_trivials.o sub_global_variables.o quick_sort.o
	$(cpl) $(cflags) cdms.f90

grid.o: grid.f90 data_struct.o sub_trivials.o interpolation.o
	$(cpl) $(cflags) grid.f90

chemistry.o: chemistry.f90 sub_trivials.o sub_global_variables.o data_struct.o quick_sort.o
	$(cpl) $(cflags) chemistry.f90

heating_cooling.o: heating_cooling.f90 sub_global_variables.o load_Neufeld_cooling_H2.o chemistry.o \
		load_Neufeld_cooling_H2O.o load_Neufeld_cooling_CO.o load_Bethell_Xray.o data_struct.o statistic_equilibrium.o lamda.o
	$(cpl) $(cflags) heating_cooling.f90

data_struct.o: data_struct.f90
	$(cpl) $(cflags) data_struct.f90

configure.o: configure.f90 disk.o chemistry.o heating_cooling.o sub_trivials.o grid.o ray_tracing.o
	$(cpl) $(cflags) configure.f90

load_Neufeld_cooling_H2.o: load_Neufeld_cooling_H2.f90 sub_trivials.o
	$(cpl) $(cflags) load_Neufeld_cooling_H2.f90

load_Neufeld_cooling_H2O.o: load_Neufeld_cooling_H2O.f90 sub_trivials.o
	$(cpl) $(cflags) load_Neufeld_cooling_H2O.f90

load_Neufeld_cooling_CO.o: load_Neufeld_cooling_CO.f90 sub_trivials.o
	$(cpl) $(cflags) load_Neufeld_cooling_CO.f90

load_Bethell_Xray.o: load_Bethell_Xray.f90 sub_global_variables.o
	$(cpl) $(cflags) load_Bethell_Xray.f90

load_Visser_CO_selfshielding.o: load_Visser_CO_selfshielding.f90 sub_trivials.o
	$(cpl) $(cflags) load_Visser_CO_selfshielding.f90

sub_trivials.o: sub_trivials.f90 sub_global_variables.o
	$(cpl) $(cflags) sub_trivials.f90

sub_global_variables.o: sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

interpolation.o: interpolation.f90
	$(cpl) $(cflags) interpolation.f90

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o
	$(cpl) $(cflags) statistic_equilibrium.f90

montecarlo.o: montecarlo.f90 sub_trivials.o sub_global_variables.o data_struct.o grid.o quick_sort.o load_Bethell_Xray.o voigt.o ray_propagating.o
	$(cpl) $(cflags) montecarlo.f90

ray_tracing.o: ray_tracing.f90 sub_trivials.o sub_global_variables.o data_struct.o grid.o quick_sort.o lamda.o hitran.o cdms.o ray_propagating.o
	$(cpl) $(cflags) ray_tracing.f90

ray_propagating.o: ray_propagating.f90 sub_trivials.o data_struct.o grid.o
	$(cpl) $(cflags) ray_propagating.f90

voigt.o: voigt.f90
	$(cpl) $(cflags) voigt.f90

load_Draine_dusts.o: load_Draine_dusts.f90 sub_trivials.o sub_global_variables.o
	$(cpl) $(cflags) load_Draine_dusts.f90

quick_sort.o: quick_sort.f90
	$(cpl) $(cflags) quick_sort.f90

opkdmain.o: opkdmain.f opkda1.o opkda2.o
	$(cpl) $(cflags) opkdmain.f

opkda1.o: opkda1.f
	$(cpl) $(cflags) opkda1.f

opkda2.o: opkda2.f
	$(cpl) $(cflags) opkda2.f

analytical_H2O_OH.o: analytical_H2O_OH.f90
	$(cpl) $(cflags) analytical_H2O_OH.f90

sync:
	rsync -acvz *.f90 fdu@shire.astro.lsa.umich.edu:~/now/src/

syncall:
	rsync -acvz *.f90 fdu@shire.astro.lsa.umich.edu:~/now/src/
	rsync -acvz makefile fdu@shire.astro.lsa.umich.edu:~/now/src/
clean:
	 rm ./*.mod ./*.o

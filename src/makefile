#cpl	= mpif90
#cpl = ifort
#cflags = -debug -save-temps -fpic -Wl,-no_pie -heap-arrays -O2 -g -traceback -check all -fp-stack-check -c
#lflags = -debug -save-temps -fpic -heap-arrays -O2 -g -traceback -check all -fp-stack-check
#cflags = -debug -c
#lflags = -debug
#cflags = -fast -c
#lflags = -fast
cpl	= gfortran
lflags = -g -W -Wall -fbounds-check #-pedantic-errors -ffpe-trap=zero,invalid,overflow,underflow
cflags = $(lflags) -c

OBJS = chemistry.o configure.o disk.o heating_cooling.o load_Neufeld_cooling_CO.o data_struct.o grid.o dust.o \
		load_Neufeld_cooling_H2.o load_Neufeld_cooling_H2O.o main.o opkda1.o opkda2.o \
		opkdmain.o sub_global_variables.o sub_trivials.o load_Bethell_Xray.o load_Visser_CO_selfshielding.o \
		interpolation.o statistic_equilibrium.o quick_sort.o \
		montecarlo.o

OBJSANA = analytical_H2O_OH.o opkdmain.o opkda1.o opkda2.o sub_global_variables.o

exe_name = a.out
exe_alt = a1.out
exe_ana = ana.out

all: $(exe_name)

alt: $(exe_alt)

ana: $(exe_ana)

$(exe_name):  \
		$(OBJS)
	$(cpl) $(lflags) $(OBJS)

$(exe_alt):  \
		$(OBJS)
	$(cpl) $(lflags) -o $(exe_alt) $(OBJS)

$(exe_ana): \
        $(OBJSANA)
	$(cpl) $(lflags) -o $(exe_ana) $(OBJSANA)

main.o: main.f90 configure.o disk.o sub_trivials.o dust.o interpolation.o statistic_equilibrium.o montecarlo.o
	$(cpl) $(cflags) main.f90

disk.o: disk.f90 chemistry.o heating_cooling.o load_Visser_CO_selfshielding.o data_struct.o grid.o
	$(cpl) $(cflags) disk.f90

dust.o: dust.f90 sub_trivials.o
	$(cpl) $(cflags) dust.f90

grid.o: grid.f90 data_struct.o sub_trivials.o interpolation.o
	$(cpl) $(cflags) grid.f90

chemistry.o: chemistry.f90 sub_trivials.o sub_global_variables.o data_struct.o quick_sort.o
	$(cpl) $(cflags) chemistry.f90

heating_cooling.o: heating_cooling.f90 sub_global_variables.o load_Neufeld_cooling_H2.o \
		load_Neufeld_cooling_H2O.o load_Neufeld_cooling_CO.o load_Bethell_Xray.o data_struct.o statistic_equilibrium.o
	$(cpl) $(cflags) heating_cooling.f90

data_struct.o: data_struct.f90
	$(cpl) $(cflags) data_struct.f90

configure.o: configure.f90 disk.o chemistry.o heating_cooling.o sub_trivials.o grid.o
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

sub_trivials.o: sub_trivials.f90
	$(cpl) $(cflags) sub_trivials.f90

sub_global_variables.o: sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

interpolation.o: interpolation.f90
	$(cpl) $(cflags) interpolation.f90

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o
	$(cpl) $(cflags) statistic_equilibrium.f90

montecarlo.o: montecarlo.f90
	$(cpl) $(cflags) montecarlo.f90

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
	#rsync -acvz *.dat fdu@shire.astro.lsa.umich.edu:~/now/
	#
	#rsync -acvz *.f90 ~/SkyDrive/backup/work/protoplanetary_disk/
	#rsync -acvz *.py  ~/SkyDrive/backup/work/protoplanetary_disk/
	#rsync -acvz *.dat ~/SkyDrive/backup/work/protoplanetary_disk/
	#
	#rsync -acvz *.f90 ~/Dropbox/work/protoplanetary_disk/
	#rsync -acvz *.py  ~/Dropbox/work/protoplanetary_disk/
	#
	#rsync -acvz *.f90 fdu@sftp.itd.umich.edu:~/now/
	#rsync -acvz *.py  fdu@sftp.itd.umich.edu:~/now/
	#rsync -acvz *.dat fdu@sftp.itd.umich.edu:~/now/

clean:
	 rm $(OBJS)

PROG =  art_mpi.exe

SRCS =  ART_Advance.f ART_Auxiliary.f ART_Density2.f \
        ART_FFT.f ART_Gravity.f \
        ART_Init.f ART_MPI_IO2.f ART_MPI_Comm1.f \
        ART_MPI_Main.f ART_Modify.f \
        ART_Particles.f ART_Tools.f ART_MPI_LG1.f

OBJS =  ART_Advance.o ART_Auxiliary.o ART_Density2.o \
        ART_FFT.o ART_Gravity.o \
        ART_Init.o ART_MPI_IO2.o ART_MPI_Comm1.o \
        ART_MPI_Main.o ART_Modify.o \
        ART_Particles.o ART_Tools.o ART_MPI_LG1.o

##CC = cc
##CFLAGS = -O
##      SP3
#FC = mpxlf90_r
#MODE =  -q64  -O3 -C -g -qsigtrap -qsmp  -qarch=pwr3  -qfixed=72
MODE =  -q64 -O3 -qstrict  -qsmp -qmaxmem=-1 -qarch=pwr4  -qfixed=72
#MODE =    -O3  -qsmp  -qarch=pwr3  -qfixed=72 -bmaxdata:0x78000000
#FC = mpif90
#FC = mpif90
#MODE = +DD64 +DSitanium2  +U77 +O3 +Ofltacc=relaxed +Oinitcheck  +Oopenmp +Odataprefetch
#LDFLAGS = +DD64 +DSitanium2 +U77 +O3 +Ofltacc=relaxed +Oinitcheck +Oopenmp +FPD +Odataprefetch -Wl,+pd1M -lmtmpi

#LDFLAGS = -O3 -tpp2 -ipo -IPF_fma -ftz -IPF_fltacc -openmp -stack_temps -Vaxlib -w -g -WB -fpp -lmpi
#FFLAGS  = -O3 -tpp2 -ipo -IPF_fma -ftz -IPF_fltacc -openmp -stack_temps -Vaxlib -w -g -WB -fpp
#FC = efc

LDFLAGS = -O2 -openmp -traceback   -mcmodel=medium  -i-dynamic -axT   -IPF_fma -ftz    -lifport -w -convert big_endian  -fpp -lmpi
FFLAGS  = -O2 -openmp -traceback   -mcmodel=medium  -i-dynamic -axT   -IPF_fma -ftz   -w -convert big_endian  -fpp
#DFLAGS = -O2 -C -g -traceback -ftz    -w -convert big_endian   -lmpi
#FFLAGS  = -O2 -C -g -traceback -ftz    -w -convert big_endian -lmpi
FC = gfortan


#MODE = +DD64 +DSitanium2  +U77 +O3 +Ofltacc=relaxed +Oinitcheck
#LDFLAGS = +DD64 +DSitanium2 +U77 +O3 +Ofltacc=relaxed +Oinitcheck +FPD -Wl,+pd1M

##    i n t e l    s c o r e
## -compiler intel is for score mode (/opt/score/bin/mpif90)
##MODE = -w  -compiler intel -O2 -w90 -w95
##MODE =  -C -w -compiler intel -O2 -w90 -w95
##MODE =  -C  -compiler intel -O2 -w90 -w95


##     i n t e l
## setenv  /opt/share/non-score/intel/mpich-1.2.4/bin
## set path  /opt/share/non-score/intel/compiler60/ia32/bin
##MODE = -w -O2 -w90     #   intel compiler


##     p g i
#FC = pgf90
## set path = ( /opt/share/non-score/pgi/mpich-1.2.4/bin $path )
## setenv LD_LIBRARY_PATH /opt/share/pgi/compiler60/ia32/lib/
##MODE = -Mbounds -mp -Bstatic -byteswapio -fast -lmpich -L/opt/share/pgi/linux86/lib
##MODE = -Mbounds -Knoieee -Bstatic -byteswapio -fast -lmpich -L/opt/share/pgi/linux86/lib
#MODE = -Knoieee -Bstatic -byteswapio -fast -lmpich -L/opt/share/pgi/linux86/lib
##MODE =  -byteswapio -fast -lmpich -L/opt/share/pgi/linux86/lib

#FFLAGS = $(MODE)

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

MPI_toART2:
	$(FC) $(FFLAGS)  MPI_toART2.f90 -o MPI_toART2.exe
clean:
	rm -f $(PROG) $(OBJS)
	rm -f *.il

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

ART_Advance.o: a_tree.h a_control.h
ART_Auxiliary.o:  a_tree.h a_control.h
ART_Density1.o: a_tree.h a_control.h
ART_Density2.o: a_tree.h a_control.h
ART_FFT.o: a_tree.h
ART_Gravity.o: a_tree.h a_control.h
ART_Init.o: a_tree.h a_control.h
ART_MPI_IO.o: a_control.h a_constant.h a_mpi.h a_setup.h
ART_MPI_IO1.o: a_control.h a_constant.h a_mpi.h a_setup.h
ART_MPI_Comm.o: a_control.h a_constant.h a_mpi.h a_setup.h
ART_MPI_Comm1.o: a_control.h a_constant.h a_mpi.h a_setup.h
ART_MPI_Main.o: a_mpi.h a_setup.h a_tree.h a_control.h
ART_Modify.o: a_tree.h a_control.h
ART_Particles.o: a_tree.h a_control.h
ART_Tools.o: a_tree.h
ART_MPI_LG.o: a_lg.h a_mpi.h a_setup.h a_tree.h a_control.h
ART_MPI_LG1.o: a_lg.h a_mpi.h a_setup.h a_tree.h a_control.h

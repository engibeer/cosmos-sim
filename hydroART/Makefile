this = .
PROG = hartCSF.x
O = .

SRCS =	a_dnew.F \
	ART_Main.F ART_Advance.F ART_Auxiliary.F ART_Control.F \
	ART_Conversions.F ART_Cooling.F  ART_Density1.F \
	ART_EoS.F ART_FFT.F ART_FU.F ART_FluxH.F \
	 ART_Gravity.F ART_IO.F ART_Init.F ART_Magic.F \
	ART_Modify.F ART_ON.F \
	ART_ParticlesF.F ART_Riemann.F ART_SF.F ART_CF.F \
	ART_SplitUpdate.F ART_Step.F ART_Sweep.F ART_TimeStep.F \
	ART_Tools.F Start_Galaxy3.F \
        ART_RadPre.F \
	hpcwall.c nrutil.c 

OBJS =	a_dnew.o \
	$(O)/ART_Main.o $(O)/ART_Advance.o $(O)/ART_Auxiliary.o \
        $(O)/ART_Control.o $(O)/ART_Conversions.o $(O)/ART_Cooling.o \
        $(O)/ART_Density1.o $(O)/ART_EoS.o $(O)/ART_FFT.o $(O)/ART_FU.o \
        $(O)/ART_FluxH.o   \
	$(O)/ART_Gravity.o  $(O)/ART_IO.o \
	$(O)/ART_Init.o $(O)/ART_Magic.o $(O)/ART_Modify.o \
	$(O)/ART_ON.o $(O)/ART_ParticlesF.o $(O)/ART_Riemann.o \
	$(O)/ART_SF.o $(O)/ART_CF.o $(O)/ART_SplitUpdate.o \
	$(O)/ART_Step.o $(O)/ART_Sweep.o $(O)/ART_TimeStep.o \
	$(O)/ART_Tools.o $(O)/Start_Galaxy3.o   \
        $(O)/ART_RadPre.o \
	$(O)/hpcwall.o   $(O)/nrutil.o 

LIBS =	


# GNU Fortran compiler
CC = gcc
CFLAGS = -O2 -w 
FC = gfortran
FFLAGS = -O2 -w -openmp 
LDFLAGS = -O2 -w -openmp 
LIB = 
CLIB = 

depall = a_def.h

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(this)/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

.SUFFIXES: $(SUFFIXES)

$(O)/ART_Advance.o: ART_Advance.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Advance.F $(LIB)
$(O)/ART_Auxiliary.o: ART_Auxiliary.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Auxiliary.F $(LIB)
$(O)/ART_Control.o: ART_Control.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Control.F $(LIB)
$(O)/ART_Conversions.o: a_numbers.h a_control.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Conversions.F $(LIB)
$(O)/ART_Cooling.o: ART_Cooling.F a_tree.h a_control.h a_clcool.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Cooling.F $(LIB)
$(O)/ART_Density1.o: ART_Density1.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Density1.F $(LIB)
$(O)/ART_EoS.o: ART_EoS.F a_tree.h a_control.h  a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_EoS.F $(LIB)
$(O)/ART_FFT.o: ART_FFT.F a_tree.h  a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS2) -c ART_FFT.F $(LIB)
$(O)/ART_FU.o: ART_FU.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_FU.F $(LIB)
$(O)/ART_FluxH.o: ART_FluxH.F a_tree.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_FluxH.F $(LIB)
$(O)/ART_Gravity.o: ART_Gravity.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Gravity.F $(LIB)
$(O)/ART_IO.o: ART_IO.F a_control.h a_tree.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_IO.F $(LIB)
$(O)/ART_Init.o: ART_Init.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Init.F $(LIB)
$(O)/ART_Magic.o: ART_Magic.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Magic.F $(LIB)
$(O)/ART_Main.o: ART_Main.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Main.F $(LIB)
$(O)/ART_Modify.o: ART_Modify.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Modify.F $(LIB)
$(O)/ART_ON.o: ART_ON.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_ON.F $(LIB)
$(O)/ART_ParticlesF.o: ART_ParticlesF.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_ParticlesF.F $(LIB)
$(O)/ART_Riemann.o: ART_Riemann.F a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Riemann.F $(LIB)
$(O)/ART_SF.o: ART_SF.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_SF.F $(LIB)
$(O)/ART_CF.o: ART_CF.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_CF.F $(LIB)
$(O)/ART_SplitUpdate.o: ART_SplitUpdate.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_SplitUpdate.F $(LIB)
$(O)/ART_Step.o: ART_Step.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Step.F $(LIB)
$(O)/ART_Sweep.o: ART_Sweep.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Sweep.F $(LIB)
$(O)/ART_TimeStep.o: ART_TimeStep.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_TimeStep.F $(LIB)
$(O)/ART_Tools.o: ART_Tools.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Tools.F $(LIB)
$(O)/ART_RadPre.o: ART_RadPre.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_RadPre.F $(LIB)
$(O)/Start_Galaxy3.o: Start_Galaxy3.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c Start_Galaxy3.F $(LIB)
$(O)/density_alloc.mod: a_dnew.F
	$(FC) $(FFLAGS) -c a_dnew.F $(LIB)

# C routines called from Fortran
$(O)/hpcwall.o: hpcwall.c 
	$(CC) $(CFLAGS) -c hpcwall.c $(CLIB)
$(O)/nrutil.o: nrutil.c nrutil.h
	$(CC) $(CFLAGS) -c nrutil.c $(CLIB) 

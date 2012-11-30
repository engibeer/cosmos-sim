this = .
PROG = hartCSF.x
O = build

SRCS =	a_dnew.F \
	ART_Main.F ART_Advance.F ART_Auxiliary.F ART_Control.F \
	ART_Conversions.F ART_Cooling.F  ART_Density1.F \
	ART_EoS.F ART_FFT.F ART_FU.F ART_FluxH.F \
	ART_Gravity.F ART_IO.F ART_Init.F ART_Magic.F \
	ART_Modify.F ART_ON.F \
	ART_ParticlesF.F ART_Riemann.F ART_SF.F ART_CF.F \
	ART_SplitUpdate.F ART_Step.F ART_Sweep.F ART_TimeStep.F \
	ART_Tools.F Start_Galaxy3.F ART_RadPre.F \
	hpcwall.c nrutil.c 

OBJS =	a_dnew.o \
	$(O)/ART_Main.o $(O)/ART_Advance.o $(O)/ART_Auxiliary.o \
    $(O)/ART_Control.o $(O)/ART_Conversions.o $(O)/ART_Cooling.o \
    $(O)/ART_Density1.o $(O)/ART_EoS.o $(O)/ART_FFT.o $(O)/ART_FU.o \
    $(O)/ART_FluxH.o $(O)/ART_Gravity.o  $(O)/ART_IO.o \
	$(O)/ART_Init.o $(O)/ART_Magic.o $(O)/ART_Modify.o \
	$(O)/ART_ON.o $(O)/ART_ParticlesF.o $(O)/ART_Riemann.o \
	$(O)/ART_SF.o $(O)/ART_CF.o $(O)/ART_SplitUpdate.o \
	$(O)/ART_Step.o $(O)/ART_Sweep.o $(O)/ART_TimeStep.o \
	$(O)/ART_Tools.o $(O)/Start_Galaxy3.o $(O)/ART_RadPre.o \
	$(O)/hpcwall.o $(O)/nrutil.o 

# GNU Fortran compiler
CC = gcc
CFLAGS = -O2 -w 
FC = gfortran
FFLAGS = -O2 -w -openmp 
LDFLAGS = -O2 -w -openmp 

depall = a_def.h

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

$(O)/ART_Advance.o: ART_Advance.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Advance.F
$(O)/ART_Auxiliary.o: ART_Auxiliary.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Auxiliary.F
$(O)/ART_Control.o: ART_Control.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Control.F
$(O)/ART_Conversions.o: a_numbers.h a_control.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Conversions.F
$(O)/ART_Cooling.o: ART_Cooling.F a_tree.h a_control.h a_clcool.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Cooling.F
$(O)/ART_Density1.o: ART_Density1.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Density1.F
$(O)/ART_EoS.o: ART_EoS.F a_tree.h a_control.h  a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_EoS.F
$(O)/ART_FFT.o: ART_FFT.F a_tree.h  a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS2) -c ART_FFT.F
$(O)/ART_FU.o: ART_FU.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_FU.F
$(O)/ART_FluxH.o: ART_FluxH.F a_tree.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_FluxH.F
$(O)/ART_Gravity.o: ART_Gravity.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Gravity.F
$(O)/ART_IO.o: ART_IO.F a_control.h a_tree.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_IO.F
$(O)/ART_Init.o: ART_Init.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Init.F
$(O)/ART_Magic.o: ART_Magic.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Magic.F
$(O)/ART_Main.o: ART_Main.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Main.F
$(O)/ART_Modify.o: ART_Modify.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Modify.F
$(O)/ART_ON.o: ART_ON.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_ON.F
$(O)/ART_ParticlesF.o: ART_ParticlesF.F a_tree.h a_control.h a_dnew.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_ParticlesF.F
$(O)/ART_Riemann.o: ART_Riemann.F a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Riemann.F
$(O)/ART_SF.o: ART_SF.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_SF.F
$(O)/ART_CF.o: ART_CF.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_CF.F
$(O)/ART_SplitUpdate.o: ART_SplitUpdate.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_SplitUpdate.F
$(O)/ART_Step.o: ART_Step.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Step.F
$(O)/ART_Sweep.o: ART_Sweep.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Sweep.F
$(O)/ART_TimeStep.o: ART_TimeStep.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_TimeStep.F
$(O)/ART_Tools.o: ART_Tools.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_Tools.F
$(O)/ART_RadPre.o: ART_RadPre.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c ART_RadPre.F
$(O)/Start_Galaxy3.o: Start_Galaxy3.F a_tree.h a_control.h a_numbers.h a_setup.h a_constant.h $(depall)
	$(FC) $(FFLAGS) -c Start_Galaxy3.F
$(O)/density_alloc.mod: a_dnew.F
	$(FC) $(FFLAGS) -c a_dnew.F

# C routines called from Fortran
$(O)/hpcwall.o: hpcwall.c 
	$(CC) $(CFLAGS) -c hpcwall.c
$(O)/nrutil.o: nrutil.c nrutil.h
	$(CC) $(CFLAGS) -c nrutil.c

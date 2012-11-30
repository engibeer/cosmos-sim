!!
#define PRESSUREFLOOR
!!
!! include N-body particles
!!
#define NBODY
!!
!! include gasdynamics
!!
#define HYDRO
!!
!! include self-gravity
!!
#define GRAVITY
!!
!! compile for cosmological runs
!!
#define COSMOLOGY
!!
!! allow mesh refinements
!!
#define REFINE
!!
!! include Lapidus viscosity
!!
#define LAPIDUS
!!
!! smoother density gradients
!!
#define DENSGRADSMOOTH
!!
!! advect extra species (defined by nchem > 0 in a_setup.h)
!!
#define ADVECT_SPECIES
!!
!! include gas cooling
!!
#define COOLING
!!
!! allow for metallicity-dependent cooling
!!
#define METALCOOLING
!!
!! use Cloudy cooling/heating table
!!
#define CLOUDY_COOLING
!!
!! allow starformation
!!
#define STARFORM
!!
!! allow for feedback from SNII
!!
#define FEEDBACK
!!
!! allow metal enrichment by SNII
!!
#define ENRICH
!!
!! allow for feedback from SNIa
!!
#define FEEDBACK_SNIa
!!
!! allow metal enrichment by SNIa
!!
#define ENRICH_SNIa
!!
!! allow stellar mass loss
!!
#define STELLARMASSLOSS
!!
!! allow AGN feedback
!!
!!#define AGNFEEDBACK
!!
!! compile for AIX OS
!!
#define OS_AIX
!!
!! include work progress output
!!
#define WORKOUT

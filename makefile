# Import variables/options/rules from PETSc.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

PETSC_MAKE_STOP_ON_ERROR=

MYFLAGS = -I.

###############################################################################
# Assign additional compiler/preprocessor flags
###############################################################################

# These flags are supplemental to the PETSc flags
CFLAGS   =
FFLAGS   = ${USER_FFLAGS} ${INC_NETCDF} -fallow-argument-mismatch
LDFLAGS  = -L${LIB_NETCDF} -lnetcdf -lnetcdff
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES = mksurfdat_petsc

mksurfdat_petsc_obj = mksurfdat_petsc.o \
  shr_kind_mod.o \
  mkvarctl.o \
  mkvarpar.o \
  mkdomainMod.o \
  nanMod.o \
  shr_kind_mod.o \
  mkutilsMod.o \
  mkncdio.o \
  shr_sys_mod.o \
  shr_log_mod.o \
  mkgridmapMod.o \
  mkpftMod.o \
  shr_string_mod.o \
  shr_timer_mod.o \
  fileutils.o \
  shr_file_mod.o \
  mkglcmecMod.o \
  mksoilMod.o \
  mkurbanparMod.o \
  mkurbanparCommonMod.o \
  mkindexmapMod.o \
  mklanwatMod.o \
  mkdiagnosticsMod.o \
  mkchecksMod.o \
  mkgdpMod.o \
  mkpeatMod.o \
  mkagfirepkmonthMod.o


mksurfdat_petsc : $(mksurfdat_petsc_obj)
	${FLINKER} -o mksurfdat_petsc $(mksurfdat_petsc_obj) ${PETSC_LIB} ${LIBS} ${FFLAGS} ${LDFLAGS}

mkindexmapMod.o : \
  shr_kind_mod.o \
  mkncdio.o

mkurbanparMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkindexmapMod.o \
  mkurbanparCommonMod.o

mkurbanparCommonMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o

mksoilMod.o : \
  shr_kind_mod.o \
  mkdomainMod.o

mkglcmecMod.o : \
  shr_kind_mod.o \
  mkdomainMod.o

mkutilsMod.o : \
  shr_kind_mod.o \
  mkncdio.o

mkncdio.o : \
  shr_sys_mod.o

mkvarctl.o :

mkvarpar.o : \
  shr_const_mod.o

mkpftMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkvarctl.o \
  mkvarpar.o \
  mkdomainMod.o \
  shr_string_mod.o \
  fileutils.o

nanMod.o :

mkpeatMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkagfirepkmonthMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkdomainMod.o : \
  shr_kind_mod.o \
  mkvarpar.o \
  mkutilsMod.o \
  nanMod.o \
  mkgridmapMod.o

mkgridmapMod.o : \
  mkvarctl.o \
  shr_kind_mod.o

mksurfdat_petsc.o : \
  shr_kind_mod.o \
  mkdomainMod.o \
  mkpftMod.o \
  mkglcmecMod.o \
  mksoilMod.o \
  mkvarctl.o \
  mkurbanparMod.o \
  mkurbanparCommonMod.o \
  mklanwatMod.o \
  mkgdpMod.o \
  mkpeatMod.o \
  mkagfirepkmonthMod.o

mklanwatMod.o	: \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o \
  mkdiagnosticsMod.o \
  mkchecksMod.o

mkchecksMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o

mkgdpMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkdiagnosticsMod.o : \
  shr_kind_mod.o

shr_const_mod.o : \
  shr_kind_mod.o

shr_kind_mod.o :

shr_log_mod.o :

shr_sys_mod.o : \
  shr_kind_mod.o \
  shr_log_mod.o

shr_string_mod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  shr_timer_mod.o \
  shr_log_mod.o

shr_timer_mod.o : \
  shr_kind_mod.o \
  shr_log_mod.o

fileutils.o : \
  shr_sys_mod.o \
  shr_file_mod.o

shr_file_mod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  shr_log_mod.o

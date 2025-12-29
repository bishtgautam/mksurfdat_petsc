# Sanity checks for required environment variables
$(if $(strip $(PETSC_DIR)),,$(error PETSC_DIR is not set. Set PETSC_DIR=/path/to/petsc_main))
$(if $(strip $(PETSC_ARCH)),,$(error PETSC_ARCH is not set. Set PETSC_ARCH=your-arch))
$(if $(strip $(SCORPIO_INSTALL)),,$(error SCORPIO_INSTALL is not set))

# Import variables/options/rules from PETSc.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

PETSC_MAKE_STOP_ON_ERROR=

MYFLAGS = -I.

###############################################################################
# Assign additional compiler/preprocessor flags
###############################################################################

# These flags are supplemental to the PETSc flags
MYFLAGS = -I${SCORPIO_INSTALL}/include 
SCORPIO_LIB=-L${SCORPIO_INSTALL}/lib -lpiof -lgptl -lpioc

CFLAGS   =
FFLAGS   = ${USER_FFLAGS} ${INC_NETCDF} -fallow-argument-mismatch
LDFLAGS  = -L${LIB_NETCDF_C} -lnetcdf -L${LIB_NETCDF_FORTRAN} -lnetcdff ${SCORPIO_LIB}
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES = mksurfdat_petsc

mksurfdat_petsc_obj = mksurfdat_petsc.o \
  shr_kind_mod.o \
  mkvarctl.o \
  mkvarpar.o \
  mkdomainMod.o \
  mkdomainPIOMod.o \
  nanMod.o \
  shr_kind_mod.o \
  mkutilsMod.o \
  mkncdio.o \
  shr_sys_mod.o \
  shr_log_mod.o \
  mkgridmapMod.o \
  mkgridmapPIOMod.o \
  piofileutils.o \
  mkpftMod.o \
  mkpftPIOMod.o \
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
  mkagfirepkmonthMod.o \
  mkVICparamsMod.o \
  mkCH4inversionMod.o \
  mkvocefMod.o \
  mkSedMod.o \
  mksoilphosphorusMod.o \
  mkfileMod.o \
  mkharvestMod.o \
  mklaiMod.o \
  mkdataPIOMod.o \
  spmdMod.o

mksurfdat_petsc : $(mksurfdat_petsc_obj)
	${FLINKER} -o mksurfdat_petsc $(mksurfdat_petsc_obj) ${PETSC_LIB} ${LIBS} ${FFLAGS} ${LDFLAGS}

spmdMod.o :

mkdataPIOMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkchecksMod.o

mkfileMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  fileutils.o \
  mkvarpar.o \
  mkvarctl.o \
  mkurbanparMod.o \
  mkglcmecMod.o \
  mkpftMod.o \
  mksoilMod.o \
  mkSedMod.o \
  mkharvestMod.o \
  mkncdio.o \
  mkdomainMod.o \
  spmdMod.o

mklaiMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o \
  mkvarctl.o

mkharvestMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkSedMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mksoilphosphorusMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkvocefMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkVICparamsMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkCH4inversionMod : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkindexmapMod.o : \
  shr_kind_mod.o \
  mkncdio.o

mkurbanparMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkindexmapMod.o \
  mkurbanparCommonMod.o \
  mkdataPIOMod.o

mkurbanparCommonMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o

mksoilMod.o : \
  shr_kind_mod.o \
  mkdataPIOMod.o \
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

mkpftPIOMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkvarctl.o \
  mkvarpar.o \
  mkdomainPIOMod.o \
  mkgridmapPIOMod.o \
  shr_string_mod.o \
  piofileutils.o \
  fileutils.o

nanMod.o :

mkpeatMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkagfirepkmonthMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkdomainPIOMod.o : \
  shr_kind_mod.o \
  mkvarpar.o \
  spmdMod.o \
  nanMod.o

mkdomainMod.o : \
  shr_kind_mod.o \
  mkvarpar.o \
  mkutilsMod.o \
  nanMod.o \
  mkgridmapMod.o

mkgridmapMod.o : \
  mkvarctl.o \
  shr_kind_mod.o \
  mkgridmapPIOMod.o

mkgridmapPIOMod.o : \
  mkvarctl.o \
  shr_kind_mod.o \
  piofileutils.o \
  spmdMod.o

mksurfdat_petsc.o : \
  shr_kind_mod.o \
  mkdomainMod.o \
  mkdomainPIOMod.o \
  mkpftMod.o \
  mkpftPIOMod.o \
  mkglcmecMod.o \
  mksoilMod.o \
  mkvarctl.o \
  mkurbanparMod.o \
  mkurbanparCommonMod.o \
  mklanwatMod.o \
  mkgdpMod.o \
  mkpeatMod.o \
  mkagfirepkmonthMod.o \
  mkVICparamsMod.o \
  mkCH4inversionMod.o \
  mkvocefMod.o \
  mkSedMod.o \
  mksoilphosphorusMod.o \
  mkfileMod.o \
  mkharvestMod.o \
  mklaiMod.o \
  mkdataPIOMod.o \
  spmdMod.o

mkdataPIOMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainPIOMod.o \
  mkgridmapMod.o \
  mkvarpar.o \
  mkchecksMod.o \
  piofileutils.o

mklanwatMod.o	: \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o \
  mkdiagnosticsMod.o \
  mkchecksMod.o \
  piofileutils.o

mkchecksMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o

mkgdpMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdataPIOMod.o \
  mkdomainMod.o

mkdiagnosticsMod.o : \
  shr_kind_mod.o

piofileutils.o : \
  mkdomainPIOMod.o \
  shr_kind_mod.o \
  spmdMod.o

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

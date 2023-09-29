# Import variables/options/rules from PETSc.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

PETSC_MAKE_STOP_ON_ERROR=

MYFLAGS = -I.

###############################################################################
# Assign additional compiler/preprocessor flags
###############################################################################

# These flags are supplemental to the PETSc flags

MYFLAGS=-DCPRAPPLECLANG -DCPRGNU -DHAVE_MPI -DINCLUDE_CMAKE_FCI -DMPICH_SKIP_MPICXX -DMPI_NO_CPPBIND -DOMPI_SKIP_MPICXX -DTIMING -DUSE_PNETCDF_VARN -DUSE_PNETCDF_VARN_ON_READ -D_NONETCDF -D_PNETCDF -I/Users/bish218/projects/e3sm/scorpio/src/flib -I/Users/bish218/projects/e3sm/scorpio/build/src/flib -I/Users/bish218/projects/petsc/petsc_99e66fd434959f7d2/clang-openmpi/include -I/Users/bish218/projects/e3sm/scorpio/src/gptl -I/Users/bish218/projects/e3sm/scorpio/build/src/gptl -I/Users/bish218/projects/e3sm/scorpio/src/clib -I/Users/bish218/projects/e3sm/scorpio/build/src/clib -g -O0 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX13.3.sdk

SCORPIO_LIB=/Users/bish218/projects/e3sm/scorpio/build/src/flib/libpiof.a /Users/bish218/projects/e3sm/scorpio/build/src/gptl/libgptl.a /Users/bish218/projects/e3sm/scorpio/build/src/clib/libpioc.a /Users/bish218/projects/petsc/petsc_99e66fd434959f7d2/clang-openmpi/lib/libpnetcdf.dylib -lmpi_cxx -lc++

CFLAGS   =
FFLAGS   = ${USER_FFLAGS} ${INC_NETCDF} -fallow-argument-mismatch
LDFLAGS  = -L${LIB_NETCDF} -lnetcdf -lnetcdff ${SCORPIO_LIB}
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
  mkagfirepkmonthMod.o \
  mkVICparamsMod.o \
  mkCH4inversionMod.o \
  mkvocefMod.o \
  mkSedMod.o \
  mksoilphosphorusMod.o \
  mkfileMod.o \
  mkharvestMod.o \
  mklaiMod.o \
  spmdMod.o

mksurfdat_petsc : $(mksurfdat_petsc_obj)
	${FLINKER} -o mksurfdat_petsc $(mksurfdat_petsc_obj) ${PETSC_LIB} ${LIBS} ${FFLAGS} ${LDFLAGS}

spmdMod.o :

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
  mkdomainMod.o

mksoilphosphorusMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkvocefMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkVICparamsMod.o : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

mkCH4inversionMod : \
  shr_kind_mod.o \
  shr_sys_mod.o \
  mkdomainMod.o

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

mkdomainPIOMod.o : \
  shr_kind_mod.o \
  mkvarpar.o \
  nanMod.o

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
  mkdomainPIOMod.o \
  mkpftMod.o \
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
  spmdMod.o

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

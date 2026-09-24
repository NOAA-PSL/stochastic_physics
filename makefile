.SUFFIXES:
.SUFFIXES: .F90 .o
.SUFFIXES: .F .o


FC = mpif90
FFLAGS_STOCH = $(FFLAGS_PROMOTION) $(FFLAGS) 
FCINCLUDES_STOCH = $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90 
RM = rm
# CORE is a make variable passed via the command line (e.g. "make CORE=atmosphere")
# in the main call to make to build MPAS.  It specifies which of the various MPAS
# cores to build.  The MPAS macro passed via the -DMPAS flag specifies to the
# pre-processor that stochastic_physics is being built for MPAS (as opposed to
# another dycore such as FV3).
ifeq ($(CORE),atmosphere)
COREDEF = -DMPAS
endif

all: stoch_physics_lib

dummy:
	echo "****** compiling stochastic physics ******"

OBJS = \
	kinddef.o \
	mpi_wrapper.o \
	stochy_internal_state_mod.o \
	stochy_namelist_def.o \
	spectral_transforms.o \
	mersenne_twister.o \
	stochy_patterngenerator.o \
	compns_stochy.o \
	stochy_nml_rec.o \
	stochy_data_mod.o \
	get_stochy_pattern.o \
	stochastic_physics_m.o \
	mpas_stochastic_physics.o	

stoch_physics_lib: $(OBJS)
	ar -ru libstochphys.a $(OBJS)

phys_interface: $(OBJS)

# DEPENDENCIES:
stochy_internal_state_mod.o: kinddef.o
stochy_namelist_def.o: kinddef.o
spectral_transforms.o: kinddef.o mpi_wrapper.o stochy_internal_state_mod.o stochy_namelist_def.o
mersenne_twister.o: kinddef.o
stochy_patterngenerator.o: kinddef.o spectral_transforms.o mersenne_twister.o mpi_wrapper.o
compns_stochy.o: stochy_namelist_def.o
stochy_nml_rec.o: kinddef.o stochy_namelist_def.o
stochy_data_mod.o: kinddef.o stochy_nml_rec.o spectral_transforms.o stochy_namelist_def.o \
	mpi_wrapper.o stochy_patterngenerator.o stochy_internal_state_mod.o \
	mersenne_twister.o compns_stochy.o
get_stochy_pattern.o: kinddef.o spectral_transforms.o stochy_namelist_def.o \
	stochy_data_mod.o stochy_patterngenerator.o stochy_internal_state_mod.o \
	mpi_wrapper.o mersenne_twister.o
stochastic_physics_m.o: kinddef.o stochy_data_mod.o stochy_namelist_def.o \
	spectral_transforms.o mpi_wrapper.o get_stochy_pattern.o
mpas_stochastic_physics.o: kinddef.o stochastic_physics_m.o

clean:
	$(RM) *.o *.mod libstochphys.a

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS_STOCH) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES_STOCH) 

.F.o:
	echo $(FFLAGS_STOCH)
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS_STOCH) -c $*.F $(CPPINCLUDES) $(FCINCLUDES_STOCH) 

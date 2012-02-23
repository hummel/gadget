OPT   +=  -DPERIODIC
OPT   +=  -DUNEQUALSOFTENINGS
OPT   +=  -DPEANOHILBERT
OPT   +=  -DWALLCLOCK
OPT   +=  -DPMGRID=128
#OPT   +=  -DPLACEHIGHRESREGION=3
#OPT   +=  -DENLARGEREGION=1.2
#OPT   +=  -DASMTH=1.25
#OPT   +=  -DRCUT=4.5
OPT   +=  -DDOUBLEPRECISION
#OPT   +=  -DDOUBLEPRECISION_FFTW
#OPT   +=  -DSYNCHRONIZATION
OPT   +=  -DFLEXSTEPS
#OPT   +=  -DPSEUDOSYMMETRIC
OPT   +=  -DNOSTOP_WHEN_BELOW_MINTIMESTEP
#OPT   +=  -DNOPMSTEPADJUSTMENT
#OPT   +=  -DHAVE_HDF5
#OPT   +=  -DOUTPUTPOTENTIAL
#OPT   +=  -DOUTPUTACCELERATION
#OPT   +=  -DOUTPUTCHANGEOFENTROPY
#OPT   +=  -DOUTPUTTIMESTEP
#OPT   +=  -DOUTPUTPRESSURE
#OPT   +=  -DOUTPUTCOLUMN
#OPT   +=  -DNOGRAVITY
#OPT   +=  -DNOTREERND
#OPT   +=  -DNOTYPEPREFIX_FFTW
#OPT   +=  -DLONG_X=3
#OPT   +=  -DLONG_Y=2
#OPT   +=  -DLONG_Z=1
#OPT   +=  -DTWODIMS
#OPT   +=  -DSPH_BND_PARTICLES
#OPT   +=  -DNOVISCOSITYLIMITER
#OPT   +=  -DCOMPUTE_POTENTIAL_ENERGY
#OPT   +=  -DLONGIDS
#OPT   +=  -DISOTHERM_EQS
OPT   +=  -DADAPTIVE_GRAVSOFT_FORGAS
#OPT   +=  -DSELECTIVE_NO_GRAVITY=2+4+8+16
#OPT   +=  -DFORCETEST=0.1
#OPT   +=  -DDEBUG_CALC_PHOTO
OPT   +=  -DDEBUG_RATE_EQ
#OPT   +=  -DMAKEGLASS=262144
#OPT   +=  -DPOLYTROPE
#OPT   +=  -DTURBULENCE
#OPT   += -DJUMP
OPT   += -DCHEMCOOL
OPT   += -DCHEMISTRYNETWORK=1
#OPT   += -DRAYTRACE
#OPT   += -DCO_SHIELDING
OPT   += -DRAYTRACE_TG
#OPT   += -DMETALS_TG
OPT += -DSINKVAL

EXEC = g2

COBJS =	accel.o \
		accrete.o \
		allocate.o \
		allvars.o \
                alpha_calc.o \
		begrun.o \
		chemcool.o \
		density.o  \
		domain.o \
		driftfac.o  \
		endrun.o \
		forcetree.o \
                ghost.o \
		global.o  \
		gravtree.o \
		gravtree_forcetest.o \
                heat_ion_rates.o \
		hydra.o \
		init.o \
		iohr.o    \
		longrange.o \
                lum_calc.o  \
		main.o \
                mdot_calc.o \
		ngb.o  \
		peano.o \
		pm_nonperiodic.o \
		pm_periodic.o \
		potential.o  \
		predict.o \
		ray.o \
		raytrace.o \
		read_ichr.o \
		restart.o \
		rsk_turbdriving_field.o \
		rsk_turbdriving_NGP.o \
		run.o \
                sink.o \
		system.o \
		timestep.o \
                compute_gamma.o

FOBJS =	calc_photo.o \
		cheminmo.o \
		const_rates.o \
		cool_func.o \
		cool_util.o \
		coolinmo.o    \
		dvode.o \
		evolve_abundances.o \
		jac.o \
		photoinit_lowZ.o \
		rate_eq.o \
		rate_eq_GMC.o \
		rate_eq_highn.o \
		rate_eq_mol.o \
		rate_eq_primordial.o \
		rate_eq_simple.o \
		spline.o \
		validate_rates.o

INCL =	allvars.h \
		chemcool_consts.h \
		proto.h \
		tags.h \
		turbulence.h \
		makefile

FINCL =	cool.h \
		chemcool_consts.h \
		fs_data.h \
		mol_data.h \
		non_eq.h \
		sd_metal.h \
		shield_data.h \
		makefile

CC        = /home/r900-1/pawlik/sw/mpich2/bin/mpicc
FC        = /home/r900-1/pawlik/sw/mpich2/bin/mpif90 -nofor-main #-Mnomain
OPTIMIZE  = -Wall -g -O3
OPTIONS = $(OPTIMIZE) $(OPT)
GSL_INCL  = -I/home/r900-1/pawlik/sw/gsl/include  
GSL_LIBS  = -L/home/r900-1/pawlik/sw/gsl/lib -lgsl -lgslcblas -lm -lfrtbegin -lg2c
FFTW_INCL = -I/home/r900-1/pawlik/sw/fftw2mpich2/include
FFTW_LIBS = -L/home/r900-1/pawlik/sw/fftw2mpich2/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw 
MPICHINCL = -I/home/r900-1/pawlik/sw/mpich2/include 
MPICHLIB  = -L/home/r900-1/pawlik/sw/mpich2/lib -lmpich
HDF5INCL  = -I/home/r900-1/pawlik/sw/hdf5/include
HDF5LIB   = -L/home/r900-1/pawlik/sw/hdf5/lib -lhdf5
#IMFLIB    = -L/opt/hpc/intel/fc_i386/10.1.015/lib -limf
#IMFINCL   = -I/opt/hpc/intel/fc_i386/10.1.015/include


CFLAGS = $(OPTIONS)  $(GSL_INCL) $(FFTW_INCL) $(MPICHINCL)

FFLAGS = $(OPTIONS)

LIBS = $(MPICHLIB) -g -lm $(GSL_LIBS) $(FFTW_LIBS) 

$(EXEC): $(FOBJS) $(COBJS)
	$(FC) $(FOBJS) $(COBJS) $(LIBS) -o $(EXEC)

$(COBJS): $(INCL)

$(FOBJS): $(FINCL)

clean:
	rm -f $(COBJS) $(FOBJS) $(EXEC)

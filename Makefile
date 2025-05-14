FC = ifort 
OMPFC = mpif90
MPIFC = mpiifort
LIBPATH = -L/data1/share/software/intel/oneapi/mkl/2022.0.2/lib/intel64 
SCALAPACK = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
LIBS =  -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64  -lmkl_intel_thread  -lmkl_core -lpthread -liomp5
#COPTS =   -qopenmp -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
COPTS =  -fopenmp
#all: shale01 shale02 shale03 shale04 shale05 shale06 shalecg02 shalecg03 shalecg04 shalecg05 shalecg06 shalecg07 shalecg08
all: shale01 shale02 shale03 shale04 shale05 shale06 shalecg02 shalecg03 shalecg04 shalecg05 shalecg06 shalecg07 shalecg08
#ifeq ($(FC),gfortran)
#	COPTS = -fopenmp
#	LIBS = -llapack -lblas
#endif
shalecg08: shale_cg_v08.f90 pnm.f90
	$(FC) $(COPTS) -o shalecg08 shale_cg_v08.f90 pnm.f90 $(LIBPATH) $(LIBS)

shalecg07: shale_cg_v07.f90 pnm.f90
	$(MPIFC) $(COPTS) -o shalecg07 shale_cg_v07.f90 pnm.f90 $(LIBPATH) $(LIBS) $(SCALAPACK)

shalecg06: shale_cg_v06.f90 pnm.f90
	$(MPIFC) $(COPTS) -o shalecg06 shale_cg_v06.f90 pnm.f90 $(LIBPATH) $(LIBS)

shale06: shale_v06.f90 pnm.f90
	$(MPIFC) $(COPTS) -o shale06 shale_v06.f90 pnm.f90 $(LIBPATH) $(SCALAPACK) $(LIBS)

shalecg05: shale_cg_v05.f90 pnm.f90
	$(FC) $(COPTS) -o shalecg05 shale_cg_v05.f90 pnm.f90 $(LIBPATH) $(LIBS)

shale05: shale_v05.f90 pnm.f90
	$(OMPFC) $(COPTS) -o shale05 shale_v05.f90 pnm.f90 $(LIBPATH) $(SCALAPACK) $(LIBS)

shale04: shale_v04.f90 pnm.f90
	$(FC) $(COPTS) -o shale04 shale_v04.f90 pnm.f90 $(LIBPATH) $(LIBS)

shalecg04: shale_cg_v04.f90 pnm.f90
	$(FC) $(COPTS) -o shalecg04 shale_cg_v04.f90 pnm.f90 $(LIBPATH) $(LIBS)

shale03: shale_v03.f90 pnm.f90
	$(FC) $(COPTS) -o shale03 shale_v03.f90 pnm.f90 $(LIBPATH) $(LIBS)

shalecg03: shale_cg_v03.f90 pnm.f90
	$(FC) $(COPTS) -o shalecg03 shale_cg_v03.f90 pnm.f90 $(LIBPATH) $(LIBS)

shalecg02: shale_cg_v02.f90 pnm.f90
	$(FC) $(COPTS) -o shalecg02 shale_cg_v02.f90 pnm.f90 $(LIBPATH) $(LIBS)

shale02: shale_v02.f90 pnm.f90
	$(FC) $(COPTS) -o shale02 shale_v02.f90 pnm.f90 $(LIBPATH) $(LIBS)

shale01: shale_v01.f90 pnm.f90
	$(FC) $(COPTS) -o shale01 shale_v01.f90 pnm.f90 $(LIBPATH) $(LIBS)
clean:
	rm -f *.o shale01 shale02 shale03 shale04 shale05 shalecg02 shalecg03 shalecg04 shalecg05 shalecg06 shalecg07 shalecg08


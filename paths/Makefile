
prefix = ~/bin
mpif90 = /opt/local/lib/openmpi/bin/mpif90
f90 = gfortran
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
lapacklib = /usr/local/lib

all: paths_TWISP.f paths_TWISP_3D.f paths_TWISP_broken.f paths_TWISP_broken_3D.f paths_analyze_3D_multi.f paths_analyze_3D_multi.f paths_analyze_broken_3D_multi.f paths_analyze_broken_multi.f paths_analyze_multi.f paths_deg_3D_multi.f paths_deg_broken_3D_multi.f paths_deg_broken_multi.f paths_deg_multi.f paths_distill_3D_multi.f paths_distill_broken_3D_multi.f paths_distill_broken_multi.f paths_distill_multi.f paths_partition.f paths_partition_3D.f paths_partition_broken.f paths_partition_broken_3D.f 
	$(f90) paths_TWISP.f -o paths_TWISP.e  $(flags) #-L$(lapacklib) -llapack -lblas
	$(f90) paths_TWISP_3D.f -o paths_TWISP_3D.e  $(flags) 
	$(f90) paths_TWISP_broken.f -o paths_TWISP_broken.e  $(flags) 
	$(f90) paths_TWISP_broken_3D.f -o paths_TWISP_broken_3D.e  $(flags) 
	$(f90) paths_analyze_3D_multi.f -o paths_analyze_3D_multi.e $(flags)
	$(f90) paths_analyze_broken_3D_multi.f -o paths_analyze_broken_3D_multi.e $(flags)
	$(f90) paths_analyze_broken_multi.f -o paths_analyze_broken_multi.e $(flags)
	$(f90) paths_analyze_multi.f -o paths_analyze_multi.e $(flags)
	$(f90) paths_deg_3D_multi.f -o paths_deg_3D_multi.e $(flags)
	$(f90) paths_deg_broken_3D_multi.f -o paths_deg_broken_3D_multi.e $(flags)
	$(f90) paths_deg_broken_multi.f -o paths_deg_broken_multi.e $(flags)
	$(f90) paths_deg_multi.f -o paths_deg_multi.e $(flags)
	$(f90) paths_distill_3D_multi.f -o paths_distill_3D_multi.e $(flags)
	$(f90) paths_distill_broken_3D_multi.f -o paths_distill_broken_3D_multi.e $(flags)
	$(f90) paths_distill_broken_multi.f -o paths_distill_broken_multi.e $(flags)
	$(f90) paths_distill_multi.f -o paths_distill_multi.e $(flags)
	$(f90) paths_partition.f -o paths_partition.e $(flags)
	$(f90) paths_partition_3D.f -o paths_partition_3D.e $(flags) 
	$(f90) paths_partition_broken.f -o paths_partition_broken.e $(flags)
	$(f90) paths_partition_broken_3D.f -o paths_partition_broken_3D.e $(flags)
	


# Makefile for Elmer
# ----------------------------------------
# Use external Geometry to create mesh
# Calculate Depth and Height for Paraview

EXECUTABLES = src/DistanceSolverRD 


NumProcs=80
NumProcx=1
NumProcy=1
NumProcz=1

InputSif=WriteOutGeo.sif



.SUFFIXES: .f90

all: clean ini grid submit

grid:
	#ElmerGrid 2 2 Mesh -partition 1 1 1 -autoclean
	#ElmerGrid 2 2 Mesh -metis 80 4  -autoclean

submit: ini

	mpirun -n 80 ElmerSolver_mpi

compile:
	elmerf90 src/WriteOutput2File.f90 -o src/WriteOutput2File
	elmerf90 src/DistanceSolveRD.f90 -o src/DistanceSolveRD
	elmerf90 src/BedrockBump.f90 -o src/BedrockBump
	elmerf90 src/GroundedMaskSolver.f90 -o src/GroundedMaskSolver
	elmerf90 src/USF_BMB.f90 -o src/USF_BMB; elmerf90 -o src/USF_Contact src/USF_Contact.f90 src/USF_Sliding.f90

clean:
	rm -fr results/*

ini:
	echo WriteOutGeo.sif > ELMERSOLVER_STARTINFO

.f90:
	elmerf90  -o WriteOutGeo.sif 1 1 1 80 $<
.c:
	gcc  -o WriteOutGeo.sif 1 1 1 80 $< -lm

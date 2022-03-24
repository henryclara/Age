# Makefile for Elmer
# ----------------------------------------
# Use external Geometry to create mesh
# Calculate Depth and Height for Paraview

EXECUTABLES = src/DistanceSolverRD 


NumProcs=1
NumProcx=1
NumProcy=1

InputSif=IceRiseRemesh.sif



.SUFFIXES: .f90

all: clean ini grid submit

grid:
	gmsh Mesh.geo -1 -2
	ElmerGrid 14 2 Mesh.msh  -autoclean

submit: ini

	mpirun -n 1 ElmerSolver_mpi

compile:
	elmerf90 src/AgeSolverRD.f90 -o src/AgeSolverRD  
	elmerf90 src/DistanceSolveRD.f90 -o src/DistanceSolveRD
	elmerf90 src/BedrockBump.f90 -o src/BedrockBump
	elmerf90 src/GroundedMaskSolver.f90 -o src/GroundedMaskSolver
	elmerf90 src/SeaLevel.f90 -o src/SeaLevel
	elmerf90 src/USF_BMB.f90 -o src/USF_BMB;  elmerf90 -o src/USF_Contact src/USF_Contact.f90 src/USF_Sliding.f90

clean:
	rm -fr results/*

ini:
	echo IceRiseRemesh.sif > ELMERSOLVER_STARTINFO

.f90:
	elmerf90  -o IceRiseRemesh.sif $<
.c:
	gcc  -o IceRiseRemesh.sif $< -lm

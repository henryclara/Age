#!/bin/bash
#SBATCH -o /work/bm1164/m300792/AgeExperiments//IceRiseAgeTestNewNewTunedIntNew/WriteOutGeo/SLURM_job.%j.%N.out
#SBATCH -e /work/bm1164/m300792/AgeExperiments//IceRiseAgeTestNewNewTunedIntNew/WriteOutGeo/SLURM_job.%j.%N.err
#SBATCH -D /work/bm1164/m300792/AgeExperiments//IceRiseAgeTestNewNewTunedIntNew/WriteOutGeo
#SBATCH -J WriteOutGeo
#SBATCH --get-user-env
#SBATCH --account=bm1164
#SBATCH --ntasks=80
#SBATCH --time=00:15:00
#SBATCH --partition=compute,compute2
#=================================================================================================================
set -e
echo Here comes the Nodelist:
echo $SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo $SLURM_JOB_PARTITION
cd $SLURM_SUBMIT_DIR

source ModulesPlusPathsMistralGCC71.sh

cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so src/MyFreeSurfaceSolver.so
echo $ELMER_HOME
echo $ELMER_SOLVER_HOME

make compile
make ini
srun -l --export=ALL --cpu_bind=cores --distribution=block:cyclic -n 80 ElmerSolver_mpi
mkdir -p GeoOut
cd src
#./ProcessSaveLine.sh
#./ProcessWriteOutput.sh

#cd ..

#cp GeoOut/* ../Remesh/DEM/
#cp GeoOut/* ../AgeForward/DEM/
#cp Mesh/mesh* ../Remesh/Mesh/

#cd ../Remesh
#sbatch Submit.sh

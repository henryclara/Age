#!/bin/bash
#SBATCH -o /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Remesh/SLURM_job.%j.%N.out
#SBATCH -e /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Remesh/SLURM_job.%j.%N.err
#SBATCH -D /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Remesh
#SBATCH -J RemeshInit
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --time=00:15:00
#SBATCH --partition=compute,compute2
#SBATCH --account=bm1164
#=================================================================================================================
set -e
source ModulesPlusPathsMistral.sh
echo Here comes the Nodelist:
echo $SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo $SLURM_JOB_PARTITION
cd $SLURM_SUBMIT_DIR

make compile
make ini
make grid
make submit

cp -r square_iso_N6/mesh*  /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Init/Mesh/
#cp -r square_iso_N6/mesh*  /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Forward/Mesh/

cd /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersAgeSolver/Init
sbatch Submit.sh 0


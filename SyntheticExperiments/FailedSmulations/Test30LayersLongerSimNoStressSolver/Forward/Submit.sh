#!/bin/bash
#SBATCH -o /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersLongerSimNoStressSolver/Forward/SLURM_job.%j.%N.out
#SBATCH -e /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersLongerSimNoStressSolver/Forward/SLURM_job.%j.%N.err
#SBATCH -D /work/bm1164/m300832/Age/SyntheticExperiments//Test30LayersLongerSimNoStressSolver/Forward
#SBATCH -J IceRiseForward
#SBATCH --get-user-env
#SBATCH --account=bm1164
#SBATCH --ntasks=80
#SBATCH --time=08:00:00
#SBATCH --partition=compute,compute2
#=================================================================================================================
set -e
echo Here comes the Nodelist:
echo $SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo $SLURM_JOB_PARTITION
cd $SLURM_SUBMIT_DIR

source ModulesPlusPathsMistral.sh

cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so src/MyFreeSurfaceSolver.so
YearCounter=$1
echo YearCounter is: $YearCounter
if [ "${YearCounter}" -lt "2" ]; then
	YearCounterFormatted=$(printf %06d $YearCounter)
	YearCounter=$(($YearCounter+50))
	YearCounterFormattedNew=$(printf %06d $YearCounter)
	cp Forward.sif.bak Forward.sif
	sed -i "s/START/${YearCounterFormatted}/g" Forward.sif
	sed -i "s/END/${YearCounterFormattedNew}/g" Forward.sif
	echo $YearCounter
	make compile
	make ini
	make grid
srun -l --export=ALL --cpu_bind=cores --distribution=block:cyclic -n 80 ElmerSolver_mpi
	mv Mesh Output${YearCounterFormattedNew}
	mkdir -p Mesh
	cp Output${YearCounterFormattedNew}/mesh* Mesh/
	cp Output${YearCounterFormattedNew}/Forward${YearCounterFormattedNew}*.result* 		 Mesh
	sbatch Submit.sh $YearCounter
fi

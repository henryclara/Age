module load gcc/7.1.0
module load openmpi/1.8.4-gcc71

export ELMER_HOME="/home/mpim/m300792/ElmerIce/elmericeScatteredGCC91/Elmer_devel_12-11-20"
export ELMER_SOLVER_HOME="$ELMER_HOME/bin"
export PATH=/home/mpim/m300792/ElmerIce/elmericeScatteredGCC91/Elmer_devel_12-11-20/bin:$PATH
#export \
#PATH=/home/mpim/m300792/ElmerIce/NETCDF_INST/include:/home/mpim/m300792/ElmerIce/NETCDF_INST/lib/libnetcdf.so:/home/mpim/m300792/ElmerIce/NETCDF_INST/lib/libnetcdff.so:/home/mpim/m300792/ElmerIce/NETCDF_INST/lib/:$PATH
export LD_LIBRARY_PATH=/home/mpim/m300792/ElmerIce/elmericeScatteredGCC91/Elmer_devel_12-11-20/share/elmersolver/lib/:$LD_LIBRARY_PATH

export \
LD_LIBRARY_PATH=/home/mpim/m300792/ElmerIce/elmericeScatteredGCC91/ThirdPartySoftware/netcdf-4.3.3/Install/include:/home/mpim/m300792/ElmerIce/elmericeScatteredGCC91/ThirdPartySoftware/netcdf-4.3.3/Install/lib:$LD_LIBRARY_PATH


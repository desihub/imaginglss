#PBS -S /bin/bash
#PBS -l mppwidth=256,walltime=00:15:00
#PBS -N Example
#PBS -o Example.out
#PBS -e Example.err
#PBS -q debug
#PBS -A m779

#
# An example "submit" script, to work on Edison at NERSC,
# to read a catalog of Decals objects and print some
# basic statics about them, with 256 mpi ranks.
#
# the script has been optimized for faster, consistent python startup
#

# Set the ImagingLSS configuration
export DECALS_PY_CONFIG=/project/projectdirs/m779/imaginglss/dr1/dr1.conf.py

cd $PBS_O_WORKDIR
#
# Load the modules we need.
# use the faster $SCRATCH version of python
module use /project/projectdirs/mpccc/usg/modulefiles/edison
module load python mpi4py
# but these packages have been pre-packaged to a payload to be
# broadcast by python-mpi -bcast
export PYTHONPAYLOAD=/project/projectdirs/m779/usg-python-2.7.9.tar.gz

# loacation of the broadcast
export PYTHON_MPI_HOME=/dev/shm

# patch up python locations
export PYTHONPATH=`echo $PYTHONPATH|sed -e "s;/scratch2/usg-python/;${PYTHON_MPI_HOME};g"`
export PYTHONHOME=${PYTHON_MPI_HOME}/2.7.9


# This is important to avoid python stating source files in $HOME
# Also means all ~/.local packages are unavailable
export PYTHONUSERBASE=$SCRATCH

# No OMP threading
export OMP_NUM_THREADS=1

aprun -n 256 -d 1 /project/projectdirs/m779/python-mpi/python-mpi \
    -bcast $PYTHONPAYLOAD \
    ./select_elgs.py
#

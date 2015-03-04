#PBS -S /bin/bash
#PBS -l mppwidth=24,walltime=00:10:00
#PBS -N Example
#PBS -o Example.out
#PBS -e Example.err
#PBS -q debug
#PBS -A m779
#
# An example "submit" script, to work on Carver at NERSC,
# to read a catalog of Decals objects and print some
# basic statics about them.
#
cd $PBS_O_WORKDIR
#
# Load the modules we need.
source /project/projectdirs/desi/software/modules/desi_environment.sh
#
# Set the path to the imaging data.
export DECALS_IMAGING=/project/projectdirs/cosmo/work/decam/cats/edr3
export DECALS_IMAGING=/project/projectdirs/cosmo/work/decam/release/edr/
export DECALS_CACHE=$GSCRATCH/desicache
#
# and set up the Python path...this is an example.
export PYTHONPATH=${PYTHONPATH}:${PBS_O_WORKDIR}
#
python << EOF
import numpy
from model.datarelease import DataRelease
dr    = DataRelease(version='EDR')
ramin = dr.observed_bricks[ 'RA1'].min()
ramax = dr.observed_bricks[ 'RA2'].max()
decmin= dr.observed_bricks['DEC1'].min()
decmax= dr.observed_bricks['DEC2'].max()
print ramin,ramax,decmin,decmax
EOF
#

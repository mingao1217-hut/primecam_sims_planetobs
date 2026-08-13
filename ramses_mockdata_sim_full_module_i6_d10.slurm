#!/bin/bash

################################################
###     SLURM script for primecam_sims      ####
###  10-det, full 3-array module, I6-        ####
###  translated focal plane, Ben's 5deg/     ####
###  +4min schedule, scan_rate_az=0.5 deg/s   ####
###					                        ####
###     Tested on Uni Koeln Ramses Cluster  ####
###             v1.0 , 08.13.2026           ####
################################################

### Stage 7 of the Ankur/Ben full-module validation:
### generate the actual TOAST boresight trajectory with a very small
### detector count (10), using:
###   - MODULE_ARRAYS=3 + --fp-source fp_f280_dettable_I6.h5, i.e. the
###     focalplane already hardware-translated to tube I6 via Ben's
###     validated method (dxi=-92.5', deta=-53.4' additive xi/eta shift),
###     NOT the old --tube-xi-off-deg/--tube-eta-off-deg runtime transform
###     -- do not pass those flags here, or the offset would apply twice.
###   - Ben's current wider Jupiter schedule
###     (schedule_Jupiter_5_0deg_plus4min.txt, the 11-min extension that
###     gave max dmin=2.050', frac<=1FWHM=59.9% for the single w2 array
###     at 0.2 deg/s).
###   - --scan-rate-az-deg-per-s 0.5, per Ben's request.
###   - atmosphere/noise left off (--sim-atm/--sim-noise both default off).
###
### Only 10 detectors are actually simulated here -- this run exists only
### to produce the boresight trajectory (times + boresight_azel). The
### dmin coverage check against all 10,308 translated truth detectors is
### done separately by diagnose_detector_proximity.py against the truth
### focalplane file directly, independent of how many detectors were
### simulated.
###
### --run-tag full_module_i6_d10 keeps this in its own output directory,
### separate from every single-array (w2-only) run so far.

### section 1 - SLURM params

#SBATCH --partition=smp
#SBATCH --account=ag-riechers # For Ramses/ CCAT

#SBATCH --ntasks-per-node 2
#SBATCH --cpus-per-task 1

#SBATCH --nodes 1

#SBATCH --job-name=mgao_jup_fullmod_i6_d10

#SBATCH --mem=50G
### Maximum requested time (days-hrs:min:sec)
#SBATCH --time 0-01:00:00

#SBATCH --output ./logs_sims/%A_%a_%x.out
#SBATCH --error ./logs_sims/%A_%a_%x.err
### To handle dir ownership issues
umask 0002
mkdir -p logs_sims
chmod g+rws logs_sims

### section 2 - Runtime env set-up
echo ""
echo "*************"
echo "Loading Modules and Env..."
module purge
### load all modules needed for Job
module load lib/mpi4py/3.1.5-gompi-2023b
module load lang/Miniconda3/23.9.0-0
### initiate conda
conda deactivate
conda activate toast3
echo "*************"

### section 3 - Job Logging
echo ""
echo "*************"
echo "Running Job..."
echo "Starting at `date`"
echo "Hostname $HOSTNAME"
echo "Job Name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Node List: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Slurm Ntasks: $SLURM_NTASKS"
echo "Number of Tasks per Node: $SLURM_NTASKS_PER_NODE"
echo "Number of CPUs per Task: $SLURM_CPUS_PER_TASK"
echo "Cores per Node: $SLURM_CPUS_ON_NODE"
echo "Total Number of Nodes: $SLURM_JOB_NUM_NODES"
echo "Current working directory is `pwd`"
echo "Python path: $(which python)"
echo "Using MPI lib: $(which mpirun)"
echo "mpi4py import version: $(python -c 'import mpi4py; print(mpi4py.__version__)')"
echo "Using GCC lib: $(which gcc)"
echo ""

### section 4 - Set Job parameters
### Set number of detectors (kept very small -- only the boresight
### trajectory is needed from this run, not per-detector timestreams)
ndets=10
### Full 3-array module
MODULE_ARRAYS=3
### Already-translated (tube I6) full-module focal plane
FP_SOURCE=fp_f280_dettable_I6.h5

### section 5 - Job Run
echo "***** EXEC SCRIPT *****"
echo `date '+%F %H:%M:%S'`
echo "***********************"
echo ""

schedule_file="schedule_Jupiter_5_0deg_plus4min.txt"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo ""
echo "Running schedule file: $schedule_file"
echo "Number of dets: $ndets"
echo "Module arrays: $MODULE_ARRAYS"
echo "FP source: $FP_SOURCE"
echo "Scan rate: 0.5 deg/s"
echo ""

### Making detector files
python -m scripts.fp_scripts.fp_trim $ndets --module-arrays $MODULE_ARRAYS --fp-source $FP_SOURCE
if [ $? -ne 0 ]; then
    echo "Error occurred in FP file generation script. Exiting."
    exit 1
fi

### Simulating timestream data: full module, I6-translated focal plane
### baked in via --fp-source (no --tube-xi-off-deg/--tube-eta-off-deg here),
### 0.5 deg/s scan rate, atmosphere/noise off (defaults).
mpirun -np $SLURM_NTASKS python3 sim_data_primecam_mpi.py \
    -s $schedule_file -d $ndets --module-arrays $MODULE_ARRAYS \
    --fp-source $FP_SOURCE --run-tag full_module_i6_d10 \
    --scan-rate-az-deg-per-s 0.5

echo ""
echo "******** DONE *********"
echo `date '+%F %H:%M:%S'`
echo "***********************"

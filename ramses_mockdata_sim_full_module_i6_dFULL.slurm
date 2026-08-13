#!/bin/bash

################################################
###     SLURM script for primecam_sims      ####
###  FULL (10,308-det), full 3-array module,  ####
###  I6-translated focal plane, Ben's 5deg/   ####
###  +4min schedule, scan_rate_az=0.5 deg/s   ####
###					                        ####
###     Tested on Uni Koeln Ramses Cluster  ####
###             v1.0 , 08.13.2026           ####
################################################
###
### Final stage of the Ankur/Ben full-module validation: the actual full
### 10,308-detector (all of w1+w2+w3) timestream run, following the clean
### 10-det (stage 7/8) and 2000-det results against this same schedule/
### geometry -- both showed strong coverage (98.6% of all 10,308 truth
### detectors within 1 FWHM; per-wafer w1=96.7%, w2=100.0%, w3=99.2%).
### This run exists to feed fit_all_detectors_beam.py for a real
### recovered-vs-truth position check across the full module, not just
### the geometric dmin proxy.
###
### Same geometry/schedule/tube handling as the stage-7/9 runs:
###   - MODULE_ARRAYS=3 + --fp-source fp_f280_dettable_I6.h5 (Ben's
###     validated I6 hardware translation baked into the focalplane file)
###     -- NOT --tube-xi-off-deg/--tube-eta-off-deg, which would double
###     apply the offset.
###   - schedule_Jupiter_5_0deg_plus4min.txt (Ben's 5deg schedule
###     geometry, +4min duration extension from earlier in the session).
###   - --scan-rate-az-deg-per-s 0.5.
###   - atmosphere/noise left off (defaults) -- clean pointing/fitting
###     check, matching every prior fit_all_detectors_beam.py validation
###     run this session.
###
### Resources bumped to Ankur's own full-module template default
### (ntasks-per-node=24, cpus-per-task=8) rather than the 12-task setup
### used for the 2000-det run: at ~10,308 dets x 322k samples each, each
### rank's share of the per-observation HDF5 write needs to stay
### comfortably under the 2GiB parallel-write limit that caused failures
### earlier this session at lower rank counts; 24 ranks keeps each rank's
### share to roughly ~1.1 GB.
###
### --run-tag full_module_i6_dFULL keeps this in its own output directory,
### separate from the stage-7 (10-det) and stage-9 (2000-det) runs.
###
### NOT YET SUBMITTED -- prepared pending review.

### section 1 - SLURM params

#SBATCH --partition=smp
#SBATCH --account=ag-riechers # For Ramses/ CCAT

#SBATCH --ntasks-per-node 24
#SBATCH --cpus-per-task 8

#SBATCH --nodes 1

#SBATCH --job-name=mgao_jup_fullmod_i6_dFULL

#SBATCH --mem=200G
### Maximum requested time (days-hrs:min:sec)
#SBATCH --time 0-01:30:00

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
### FULL selects every detector in the module (no trimming)
ndets=FULL
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
    --fp-source $FP_SOURCE --run-tag full_module_i6_dFULL \
    --scan-rate-az-deg-per-s 0.5

echo ""
echo "******** DONE *********"
echo `date '+%F %H:%M:%S'`
echo "***********************"

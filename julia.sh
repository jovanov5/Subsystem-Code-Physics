#!/bin/bash -l
echo =========================================================   
echo Job submitted  date = Tue  1 Oct 11:05:58 BST 2024      
date_start=`date +%s`
echo $SLURM_JOB_NUM_NODES nodes \( $SLURM_CPUS_ON_NODE processes per node \)        
echo $SLURM_JOB_NUM_NODES hosts used: $SLURM_JOB_NODELIST      
# Set this otherwise a different transport gets selected on some nodes and things break in strange ways
export OMPI_MCA_pml=^cm
echo Job output begins                                           
echo -----------------                                           
echo   
#hostname

# Need to set the max locked memory very high otherwise IB can't allocate enough and fails with "UCX  ERROR Failed to allocate memory pool chunk: Input/output error"
ulimit -l unlimited

export OMP_NUM_THEADS=1
 /usr/local/shared/slurm/bin/srun -u -n 1 --mpi=pmix_v4 --mem-per-cpu=6144 nice -n 10 /usr/local/shared/julia/julia-1.10.0/bin/julia /mnt/users/jovanovic/GitHub/Subsystem-Code-Physics/Remote-Hydra/Toric_Code/toric_code-Zassenhausen.jl /mnt/users/jovanovic/GitHub/Subsystem-Code-Physics/Remote-Hydra/Toric_Code/data_z/args/argfile_exp:99.json
  echo ---------------                                           
  echo Job output ends                                           

  date_end=`date +%s`
  seconds=$((date_end-date_start))
  minutes=$((seconds/60))
  seconds=$((seconds-60*minutes))
  hours=$((minutes/60))
  minutes=$((minutes-60*hours))
  echo =========================================================   
  echo PBS job: finished   date = `date`   
  echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
  echo =========================================================

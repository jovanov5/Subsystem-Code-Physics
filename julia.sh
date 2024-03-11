#!/bin/bash -l
echo =========================================================   
echo Job submitted  date = Wed  6 Mar 16:31:39 GMT 2024      
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
 /usr/local/shared/slurm/bin/srun -u -n 1 --mpi=pmix_v4 --mem-per-cpu=4096 nice -n 10 /usr/local/shared/julia/julia-1.10.0/bin/julia /mnt/users/jovanovic/GitHub/Subsystem-Code-Physics/Remote-Hydra/Toric_Code/toric_code-TEE-FermiBoseCond.jl 15 3 0.5 11 0.5 11 100000.0 10 11 400
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

#!/bin/bash
#
# Project/Account (use your own)
#SBATCH -A stabile
#
# Number of MPI tasks
##SBATCH -n 2
#
# Number of tasks per node
##SBATCH --tasks-per-node=1
#
# Runtime of this jobs is less then 12 hours.
##SBATCH --time=12:00:00
#
# Name
#SBATCH -J "allgather"
#
# Partition
##SBATCH --partition=mpi
#
##SBATCH --output=bandwidth_%a.out
##SBATCH --nodelist=nodo[06-07]
#SBATCH --distribution=cyclic

echo $SLURM_JOB_NODELIST

export NODELIST=outputfiles/nodelist.$$
exe="mesh"
dir="mesh_output"
#parameters
size_matrix_row="131072"
size_matrix_col="131072"
procs_row="8"
procs_col="8"
procs=$(($procs_row*$procs_col))

mkdir -p ${dir}
srun -l bash -c 'hostname' | sort | awk '{print $2}' > $NODELIST

cat $NODELIST > outputfiles/myhostfile.$$ 
echo "-----------------------------------------------"

for i in ${procs}    
do
    mpirun -np $i --map-by node --mca btl '^openib' --mca pml ucx \
    --oversubscribe --mca coll_tuned_use_dynamic_rules 1 \
    --hostfile outputfiles/myhostfile.$$ --mca mpi_warn_on_fork 0 \
        ./${exe} $size_matrix_row $size_matrix_col $procs_row $procs_col  > ${dir}/mesh_p${i}.dat 
done

# End of submit file
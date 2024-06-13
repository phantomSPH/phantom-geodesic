#!/bin/bash
#
#SBATCH --job-name=test_emb_arr
#SBATCH --partition=sstar
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00
#SBATCH --mem-per-cpu=100
#SBATCH --output=geodesic.in.out
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=megha.sharma@monash.edu

ulimit -s unlimited
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "HOSTNAME = $HOSTNAME"
echo "HOSTTYPE = $HOSTTYPE"
echo Time is `date`
echo Directory is `pwd`

echo "starting geodesic run..."

start=100.000  # Modify the start value as per your requirement
end=105.000    # Modify the end value as per your requirement
grid_points=100  # Modify the number of grid points as per your requirement

log_base=$(awk -v start="$start" -v end="$end" -v grid_points="$grid_points" 'BEGIN { print exp(log(end/start) / (grid_points-1)) }')

rp_newton_values=()

for ((i=0; i<grid_points; i++)); do
    value=$(awk -v start="$start" -v i="$i" -v log_base="$log_base" 'BEGIN { print start * (log_base^i) }')
    rp_newton_values+=($value)
done

rp_newton_values=(100.049)

echo "rp_newton,i" >> rp_i_values.txt
i=0
for rp_newton in "${rp_newton_values[@]}"; do
    srun ./new_run_geo.sh "$rp_newton" $i
    echo "$rp_newton, $i" >> rp_i_values.txt
    ((i++))
done




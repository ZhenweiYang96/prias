#!/bin/bash
#SBATCH -p normal
#SBATCH -n 250
#SBATCH -t 24:00:00
module load pre2019
module load R
for seed in {2001..2010}
do
	for testId in {751..1000}
	do
	(
	  export MKL_NUM_THREADS=1
	  Rscript controller_schedule_start_year1.R $seed $testId > "/home/at38680/log/pat_${seed}_${testId}.out" 2>&1
	) &
	done
	wait
done

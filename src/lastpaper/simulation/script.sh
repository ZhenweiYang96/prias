#!/bin/bash
batch_size=5
max_batch_number=250/batch_size

for ((seed=2001; seed<=2010; seed++))
do
	for ((batch_number=1; batch_number<=max_batch_number; batch_number++))
	do
		for ((pat_counter=1; pat_counter<=batch_size; pat_counter++))
		do
		(
		  testId=$(( 750 + pat_counter + batch_size * (batch_number-1) ))
		  echo $testId
		  sleep 5
		  # Rscript controller_schedule_start_year1.R $seed $testId > "/home/atomer/log/pat_${seed}_${testId}.out" 2>&1
		) &
		done
		wait
	done
done

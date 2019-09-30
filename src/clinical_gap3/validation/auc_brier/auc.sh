#!/bin/bash
for i in {1..10}
do
	Rscript auc_calculator.R Hopkins $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {11..20}
do
	Rscript auc_calculator.R Hopkins $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {21..30}
do
	Rscript auc_calculator.R Hopkins $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {1..10}
do
	Rscript auc_calculator.R MSKCC $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {11..20}
do
	Rscript auc_calculator.R MSKCC $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {21..30}
do
	Rscript auc_calculator.R MSKCC $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {1..10}
do
	Rscript auc_calculator.R PRIAS $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {11..20}
do
	Rscript auc_calculator.R PRIAS $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait

for i in {21..30}
do
	Rscript auc_calculator.R PRIAS $i > "/home/atomer/log/hopkins${i}.out" 2>&1 &
done

wait





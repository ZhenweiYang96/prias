#!/bin/bash
for batch in {0..10}
do
  for i in {1..5}
  do
	  Rscript auc_calculator.R $((batch*5 + i)) > "/home/atomer/log/PRIAS_$((batch*5 + i)).out" 2>&1 &
  done
  wait
done





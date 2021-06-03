#! /bin/bash


for f in ./*job.sh; do
	echo $f
	$f; #sbatch

done
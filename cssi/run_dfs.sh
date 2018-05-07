#!/bin/sh
clear

n_cores="$1"	# number of processors
max="$2"			# number of random instances to be used


echo "from math import *" > ./dfs_results.py;
echo "from numpy import median" >> ./dfs_results.py;
echo "from numpy import mean" >> ./dfs_results.py;
echo "import sys" >> ./dfs_results.py;
echo "def printf(format, *args):" >> ./dfs_results.py;
echo "    sys.stdout.write(format % args)" >> ./dfs_results.py;
echo "" >> ./dfs_results.py;
echo "" >> ./dfs_results.py;

echo "print \"\"" >> ./dfs_results.py;
echo "clock_cycles_dfs = []" >> ./dfs_results.py;
echo "running_time_dfs_dbls = []" >> ./dfs_results.py;
echo "running_time_dfs_adds = []" >> ./dfs_results.py;
echo "running_time_dfs_comp = []" >> ./dfs_results.py;
echo "running_time_dfs_eval = []" >> ./dfs_results.py;

for i in `seq 1 $max`
do
	sh config.sh;	# generating a new random instances
	sh build.sh;	# compiling the new parameters

	# dfs
	(./bin/dfs "$n_cores" A) &> "./dfs_outputs/$i.dfs.out";
	magma "./dfs_outputs/$i.dfs.out";
done

echo "print \"Here, M. R-T means Measured running time, and \"" >> ./dfs_results.py;
echo "print \"      M. CC means Measured Clock Cycles\"" >> ./dfs_results.py;
echo "print \"\"" >> ./dfs_results.py;
echo "print \"\"" >> ./dfs_results.py;

#dfs
echo "print \"Median (DFS):\"" >> ./dfs_results.py;
echo "printf(\"M. R-T DOUBLINGS: %3.2lf\\n\", log(median(running_time_dfs_dbls), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T ADDITIONS: %3.2lf\\n\", log(median(running_time_dfs_adds), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T 2-ISOGENY COMPUTATIONS: %3.2lf\\n\", log(median(running_time_dfs_comp), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T 2-ISOGENY EVALUATIONS: %3.2lf\\n\", log(median(running_time_dfs_eval), 2))" >> ./dfs_results.py;
echo "printf(\"M. CC : %3.2lf\\n\", log(median(clock_cycles_dfs), 2))" >> ./dfs_results.py;
echo "print \"\"" >> ./dfs_results.py;

echo "print \"Mean (DFS):\"" >> ./dfs_results.py;
echo "printf(\"M. R-T DOUBLINGS: %3.2lf\\n\", log(mean(running_time_dfs_dbls), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T ADDITIONS: %3.2lf\\n\", log(mean(running_time_dfs_adds), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T 2-ISOGENY COMPUTATIONS: %3.2lf\\n\", log(mean(running_time_dfs_comp), 2))" >> ./dfs_results.py;
echo "printf(\"M. R-T 2-ISOGENY EVALUATIONS: %3.2lf\\n\", log(mean(running_time_dfs_eval), 2))" >> ./dfs_results.py;
echo "printf(\"M. CC : %3.2lf\\n\", log(mean(clock_cycles_dfs), 2))" >> ./dfs_results.py;
echo "print \"\"" >> ./dfs_results.py;
echo "print \"\"" >> ./dfs_results.py;

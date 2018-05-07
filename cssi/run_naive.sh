#!/bin/sh
clear

n_cores="$1"	# number of processors
max="$2"			# number of random instances to be used


echo "from math import *" > ./naive_results.py;
echo "from numpy import median" >> ./naive_results.py;
echo "from numpy import mean" >> ./naive_results.py;
echo "import sys" >> ./naive_results.py;
echo "def printf(format, *args):" >> ./naive_results.py;
echo "    sys.stdout.write(format % args)" >> ./naive_results.py;
echo "" >> ./naive_results.py;
echo "" >> ./naive_results.py;

echo "print \"\"" >> ./naive_results.py;
echo "clock_cycles_naive = []" >> ./naive_results.py;
echo "running_time_naive = []" >> ./naive_results.py;

for i in `seq 1 $max`
do
	sh config.sh;	# generating a new random instances
	sh build.sh;	# compiling the new parameters

	# NAIVE
	(./bin/naive "$n_cores" A) &> "./naive_outputs/$i.naive.out";
	magma "./naive_outputs/$i.naive.out";
done

echo "print \"Here, M. R-T means Measured running time, and \"" >> ./naive_results.py;
echo "print \"      M. CC means Measured Clock Cycles\"" >> ./naive_results.py;
echo "print \"\"" >> ./naive_results.py;
echo "print \"\"" >> ./naive_results.py;

#NAIVE
echo "print \"Median (Naive):\"" >> ./naive_results.py;
echo "printf(\"M. R-T : %3.2lf\\n\", log(median(running_time_naive), 2))" >> ./naive_results.py;
echo "printf(\"M. CC : %3.2lf\\n\", log(median(clock_cycles_naive), 2))" >> ./naive_results.py;
echo "print \"\"" >> ./naive_results.py;

echo "print \"Mean (Naive):\"" >> ./naive_results.py;
echo "printf(\"M. R-T : %3.2lf\\n\", log(mean(running_time_naive), 2))" >> ./naive_results.py;
echo "printf(\"M. CC : %3.2lf\\n\", log(mean(clock_cycles_naive), 2))" >> ./naive_results.py;
echo "print \"\"" >> ./naive_results.py;
echo "print \"\"" >> ./naive_results.py;


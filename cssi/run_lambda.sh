#!/bin/sh
clear

n_cores="$1"	# number of processors
var_w="$2"		# w : number of elements in memory
beta="$3"		# \beta
max="$4"			# number of random instances to be used
Fs="$5"			# Number of functions to be used


echo "from math import *" > ./lambda_results.py;
echo "from numpy import median" >> ./lambda_results.py;
echo "from numpy import mean" >> ./lambda_results.py;
echo "import sys" >> ./lambda_results.py;
echo "def printf(format, *args):" >> ./lambda_results.py;
echo "    sys.stdout.write(format % args)" >> ./lambda_results.py;
echo "" >> ./lambda_results.py;
echo "" >> ./lambda_results.py;

echo "print \"\"" >> ./lambda_results.py;
echo "w = 2**$var_w" >> ./lambda_results.py;
echo "print \"w : 2^$var_w\\n\"" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;

echo "running_time_dist_pnts = []" >> ./lambda_results.py;
echo "running_time_collision = []" >> ./lambda_results.py;

echo "clock_cycles_lambda = []" >> ./lambda_results.py;
echo "running_time_lambda = []" >> ./lambda_results.py;
echo "diffunctions_lambda = []" >> ./lambda_results.py;

echo "c = []" >> ./lambda_results.py;
echo "t = []" >> ./lambda_results.py;

for i in `seq 1 $max`
do
	sh config.sh;	# generating a new random instances
	sh build.sh;	# compiling the new parameters
		
	# LAMBDA 
	(./bin/lambda "$n_cores" A "$var_w" "$beta" "$Fs") &> "./lambda_outputs/$i.lambda.out";
	magma "./lambda_outputs/$i.lambda.out";
done

echo "print \"Here, M. R-T means Measured running time, and \"" >> ./lambda_results.py;
echo "print \"      M. CC means Measured Clock Cycles\"" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;

#LAMBDA
echo "print \"Median (VW):\"" >> ./lambda_results.py;
echo "printf(\"M. R-T : %3.2lf\\n\", log(median(running_time_lambda), 2))" >> ./lambda_results.py;
echo "printf(\"M. R-T for all dist. points reached : %3.2lf\\n\", log(median(running_time_dist_pnts), 2))" >> ./lambda_results.py;
echo "printf(\"M. R-T for all collisions reached : %3.2lf\\n\", log(median(running_time_collision), 2))" >> ./lambda_results.py;
echo "printf(\"Number of functions : %3.2lf\\n\", median(diffunctions_lambda))" >> ./lambda_results.py;
echo "printf(\"M. CC : %3.2lf\\n\", log(median(clock_cycles_lambda), 2))" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;

echo "print \"Mean (VW):\"" >> ./lambda_results.py;
echo "printf(\"M. R-T : %3.2lf\\n\", log(mean(running_time_lambda), 2))" >> ./lambda_results.py;
echo "printf(\"M. R-T for all dist. points reached : %3.2lf\\n\", log(mean(running_time_dist_pnts), 2))" >> ./lambda_results.py;
echo "printf(\"M. R-T for all collisions reached : %3.2lf\\n\", log(mean(running_time_collision), 2))" >> ./lambda_results.py;
echo "printf(\"Number of functions : %3.2lf\\n\", mean(diffunctions_lambda))" >> ./lambda_results.py;
echo "printf(\"M. CC : %3.2lf\\n\",  log(mean(clock_cycles_lambda), 2))" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;

echo "print \"Different collisions reached:\""  >> ./lambda_results.py;
echo "seq_c_median = []" >> ./lambda_results.py;
echo "seq_c_mean = []" >> ./lambda_results.py;
echo "for c_i in c:" >> ./lambda_results.py;
echo "    seq_c_median.append(median(c_i))" >> ./lambda_results.py;
echo "    seq_c_mean.append(mean(c_i))" >> ./lambda_results.py;
echo "" >> ./lambda_results.py;
echo "printf(\"Medians: %1.2lfw\\n\", median(seq_c_median) / float(w))" >> ./lambda_results.py;
echo "printf(\"Means: %1.2lfw\\n\", mean(seq_c_mean) / float(w))" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;

echo "print \"Total collisions reached:\""  >> ./lambda_results.py;
echo "seq_t_median = []" >> ./lambda_results.py;
echo "seq_t_mean = []" >> ./lambda_results.py;
echo "for t_i in t:" >> ./lambda_results.py;
echo "    seq_t_median.append(median(t_i))" >> ./lambda_results.py;
echo "    seq_t_mean.append(mean(t_i))" >> ./lambda_results.py;
echo "" >> ./lambda_results.py;
echo "printf(\"Medians: %1.2lfw\\n\", median(seq_t_median) / float(w))" >> ./lambda_results.py;
echo "printf(\"Means: %1.2lfw\\n\", mean(seq_t_mean) / float(w))" >> ./lambda_results.py;
echo "print \"\"" >> ./lambda_results.py;


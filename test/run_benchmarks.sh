#!/bin/bash

efforts="1 3 6 9"
benchmarks=$(find ISPD06 -name "*.aux" | sort)
names=""

for benchmark in $benchmarks
do
	name=$(basename "${benchmark}" .aux)
	names="${names} ${name}"
done

echo "Running placement benchmarks"
echo "Benchmarks: ${names}"
echo "Efforts: ${efforts}"
echo

echo

mkdir -p logs

for effort in $efforts
do
	resfile="results_${effort}.csv"
	echo "Benchmark Effort WL" > "${resfile}"
	for benchmark in $benchmarks
	do
		name=$(basename "${benchmark}" .aux)
		logfile="logs/log_${name}_${effort}.txt"

		# Run placement
		echo "  Running benchmark ${name} effort ${effort}"
		coloquinte ${benchmark} --effort $effort 2>&1 > "${logfile}" || echo "      Failure"

		# Write the results (WL) in a csv
		echo -n "${name} ${effort} " >> "${resfile}"
		grep "Detailed placement done" "${logfile}" | sed "s/.*(WL \(.*\)).*/\1/" >> "${resfile}"
	done
done

#! /bin/sh
#$ -N CH_final_spins
#$ -cwd
#$ -l h_rt=160:00:00
#$ -l h_vmem=16G
#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
g++ -std=c++11 -g -O2 -mcmodel=large ../CH_MC_heatbath.cpp -o outp.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files

>average_spin_final.txt

for i in $(seq 1 1 10);
do
	>energy.txt
	>energy2.txt
	>spin_total.txt
	>spin_total_even.txt
	>spin_total_odd.txt
	>spin2_total.txt
	>spin2_total_even.txt
	>spin2_total_odd.txt
	>spins.txt
	>spins2.txt
	>spin_corr.txt

	#temperature step size
	dT=0.1

	#Run simulation

	#initial high temperature run
	./outp.o none 5 1
	#cooling system
	for j in $(seq 49 -1 1); do
	kT=$(echo "scale=4; $j*$dT"|bc)
	./outp.o spins_after.txt $kT 1
	done

	#Copy results from current simulation into new subdirectory
	mkdir result$i
	cp energy.txt result$i
	cp energy2.txt result$i
	cp spin_total.txt result$i
	cp spin_total_even.txt result$i
	cp spin_total_odd.txt result$i
	cp spin2_total.txt result$i
	cp spin2_total_even.txt result$i
	cp spin2_total_odd.txt result$i
	cp spins.txt result$i
	cp spins2.txt result$i
	cp spin_corr.txt result$i
	v=$(tail -n 1 spin_total.txt)
	echo $(cut -c 4- <<< $v) >> average_spin_final.txt


done
#! /bin/sh
#$ -N HL_cooling
#$ -cwd
#$ -l h_rt=100:00:00
#$ -l h_vmem=16G
#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
g++ -std=c++11 -g -O2 -mcmodel=medium ../HL_MC_heatbath.cpp -o outp.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files

#temperature step size
dT=0.1
u0=(0.5 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 50 100 200)

#initial high temperature run

for i in ${u0[*]}; do

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

	./outp.o none 5 1 $i

	#cooling system
	for j in $(seq 49 -1 1); do
	kT=$(echo "scale=4; $j*$dT"|bc)
	./outp.o spins_after.txt $kT 1 $i
	done
	
	mkdir result$i
	cp energy.txt result$i
	cp energy2.txt result$i
	cp spin_total.txt result$i
	#cp spin_total_even.txt result$i
	#cp spin_total_odd.txt result$i
	cp spin2_total.txt result$i
	#cp spin2_total_even.txt result$i
	#cp spin2_total_odd.txt result$i
	cp spins.txt result$i
	cp spins2.txt result$i
	#cp spin_corr.txt result$i



done



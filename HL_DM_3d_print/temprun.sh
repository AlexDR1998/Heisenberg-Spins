#! /bin/sh
#$ -N CH_cooling
#$ -cwd
#$ -l h_rt=04:00:00
#$ -l h_vmem=8G
#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
g++ -std=c++11 -g -O2 -mcmodel=medium ../HL_MC_heatbath.cpp -o out.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files
>debug.txt
>energy.txt
>energy2.txt
>spin_total.txt
>spin_total_even.txt
>spin_total_even.txt
>spin2_total.txt
>spin2_total_even.txt
>spin2_total_odd.txt
>spins.txt
>spins2.txt
>spin_corr.txt

#temperature step size
dT=0.1

#initial high temperature run
./out.o none 5 1 10



#cooling system
for j in $(seq 49 -1 1); do
kT=$(echo "scale=4; $j*$dT"|bc)
./out.o spins_after.txt $kT 1 10
done







#! /bin/sh
#$ -N CH_cooling
#$ -cwd
#$ -l h_rt=04:00:00
#$ -l h_vmem=8G
#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
g++ -std=c++11 -g -O3 -mcmodel=large ../CH_MC_heatbath.cpp -o outp.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files
>energy.txt
>energy2.txt
>spin_total.txt
>spin2_total.txt
>spins.txt
>spins2.txt
>spin_corr.txt

#temperature step size
dT=0.02

#initial high temperature run
./outp.o none 10 1

#cooling system
for j in $(seq 500 -1 1); do
kT=$(echo "scale=4; $j*$dT"|bc)
./outp.o spins_after.txt $kT 1
done

#! /bin/sh
#$ -N CH_cooling
#$ -cwd
#$ -l h_rt=100:00:00
#$ -l h_vmem=16G
#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files
mkdir result$SGE_TASK_ID
cd result$SGE_TASK_ID

g++ -std=c++11 -g -O2 -mcmodel=medium ../../HL_MC_heatbath.cpp -o outp.o
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

#Set parameters
dT=0.05
#u0=$SGE_TASK_ID

#initial high temperature run
./outp.o none 5 1 $SGE_TASK_ID

#cooling system
for j in $(seq 99 -1 1); do
kT=$(echo "scale=4; $j*$dT"|bc)
./outp.o spins_after.txt $kT 1 $SGE_TASK_ID
done


#Store results

#cp energy.txt result$SGE_TASK_ID
#cp energy2.txt result$SGE_TASK_ID
#cp spin_total.txt result$SGE_TASK_ID
#cp spin_total_even.txt result$SGE_TASK_ID
#cp spin_total_odd.txt result$SGE_TASK_ID
#cp spin2_total.txt result$SGE_TASK_ID
#cp spin2_total_even.txt result$SGE_TASK_ID
#cp spin2_total_odd.txt result$SGE_TASK_ID
#cp spins.txt result$SGE_TASK_ID
#cp spins2.txt result$SGE_TASK_ID

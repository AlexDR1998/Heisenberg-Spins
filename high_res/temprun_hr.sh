#! /bin/sh
#$ -N Cooling_high_res
#$ -cwd
#$ -l h_rt=100:00:00
#$ -l h_vmem=16G


mkdir result$SGE_TASK_ID
#cp outp result$1
cp system.txt result$SGE_TASK_ID
cd result$SGE_TASK_ID


#compile c++ script first
#g++ -std=c++11 -g -O3 -mcmodel=medium source.cpp -o outp.o
g++ -std=c++11 -g -O2 -mcmodel=medium ../HL_MC_heatbath.cpp -o outp.o
#g++ -std=c++11 -g -O3 IL_MC_heatbath_split.cpp -o outp.o
#empties the output files
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
dT=0.2

#initial high temperature run
./outp.o none 5 1 $SGE_TASK_ID

#cooling system at low temperature resolution
for j in $(seq 24 -1 1); do
kT=$(echo "scale=4; $j*$dT"|bc)
./outp.o spins_after.txt $kT 1 $SGE_TASK_ID
done



peak=$(python3 peak.py $dT)
t_start=$(echo "$peak+$dT"|bc)
t_end=$(echo "$peak-$dT"|bc)
echo $t_start
echo $peak
echo $t_end

#Empty output files

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

#Re-run at higher temperature resolution, only around the peak
./outp.o none $t_start 1 $SGE_TASK_ID
dT=0.01
for j in $(seq 39 -1 1); do
kT=$(echo "scale=4; $j*$dT+$t_end"|bc)
./outp.o spins_after.txt $kT 1 $SGE_TASK_ID
done

peak=$(python3 peak.py $dT)
t_start=$(echo "$peak+2*$dT"|bc)
t_end=$(echo "$peak-2*$dT"|bc)
echo $t_start
echo $peak
echo $t_end


#Empty output files

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

#Re-run at higher temperature resolution, only around the peak
./outp.o none $t_start 1 $SGE_TASK_ID
dT=0.001
for j in $(seq 39 -1 1); do
kT=$(echo "scale=4; $j*$dT+$t_end"|bc)
./outp.o spins_after.txt $kT 1 $SGE_TASK_ID
done


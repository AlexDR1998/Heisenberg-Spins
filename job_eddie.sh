#! /bin/sh
#$ -N Heisenberg_MC
#$ -cwd
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G

. /etc/profile.d/modules.sh

module load intel/2016

#compile c++ script first
g++ -std=c++11 -g -O3 -mcmodel=medium CH_MC_heatbath.cpp -o outp.o

#high temperature run
./outp.o none 10.0 1


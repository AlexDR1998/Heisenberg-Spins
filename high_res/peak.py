import numpy as np
import scipy.signal as sp
import matplotlib.pyplot as plt 
import sys

#short script that just finds the peak temperature in heat capacity plot
size=10

def peak_find():
	dT = sys.argv[1]
	mean_energies = np.loadtxt("energy.txt").transpose()
	d = -(np.gradient(mean_energies[1])/size**3)/float(dT)
	peaks,_ = sp.find_peaks(d,height=1,distance=10)
	#print(len(peaks))
	return mean_energies[0][peaks]

def main():
	print(peak_find()[0])

main()
import numpy as np 
import sys

#short script that just finds the peak temperature in heat capacity plot
size=20

def peak_find():
	dT = sys.argv[1]
	means = np.loadtxt("energy.txt").transpose()
	means_squared = np.loadtxt("energy2.txt").transpose()
	#d = -(np.gradient(means[1])/size**3)/float(dT)
	d = (means_squared[1]-(means[1]*means[1]))/(means[0]*means[0]*size**3)
	peak = np.argmax(d)
	return mean_energies[0][peak]

def main():
	print(peak_find())

main()
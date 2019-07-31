import numpy as np 
import sys

#short script that just finds the peak temperature in heat capacity plot
size=20

def peak_find():
	dT = sys.argv[1]
	mean_energies = np.loadtxt("energy.txt").transpose()
	d = -(np.gradient(mean_energies[1])/size**3)/float(dT)
	peak = np.argmax(d)
	return mean_energies[0][peak]

def main():
	print(peak_find())

main()
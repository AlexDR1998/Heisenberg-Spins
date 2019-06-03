#Simple python script to plot results of IL_MC_heatbath_split from txt files
import numpy as np
import matplotlib.pyplot as plt 
import sys

def main():
	n=1
	file_available=True
	while file_available:
		try:
			file = sys.argv[n]
			n=n+1
		except Exception as e:
			file_available=False
			break

		data = np.loadtxt(file)
		plt.plot(data)
	plt.show()

main()
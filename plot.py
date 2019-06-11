#Simple python script to plot results of IL_MC_heatbath_split from txt files
import numpy as np
import matplotlib.pyplot as plt 
import sys

def multiplot():
	#plots multiple simple (1D) data sets together
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


def single_plot():
	file = sys.argv[1]
	data = np.loadtxt(file).transpose()
	#data = data.transpose()
	plt.scatter(data[0],data[1])
	plt.show()


def heat_capacity():
	#Input files should be: energy.txt energy2.txt
	f1 = sys.argv[1]
	f2 = sys.argv[2]
	means = np.loadtxt(f1).transpose()
	squared_means = np.loadtxt(f2).transpose()
	heat_capacity = (squared_means[1]-(means[1]*means[1]))/(means[0]*means[0])
	plt.scatter(means[0],heat_capacity)
	plt.show()

def difference_plot():
	#Plots difference of 2 files
	f1 = sys.argv[1]
	f2 = sys.argv[2]
	d1 = np.loadtxt(f1).transpose()
	d2 = np.loadtxt(f2).transpose()
	plt.plot(d1-d2)
	plt.show()

def main():
	#multiplot()
	single_plot()
	#heat_capacity()
	#difference_plot()
main()

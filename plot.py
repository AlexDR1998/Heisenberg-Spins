#Simple python script to plot results of IL_MC_heatbath_split from txt files
import numpy as np
#import scipy
import scipy.signal as sp
#from scipy.signal import find_peaks
import matplotlib.pyplot as plt 
import sys

size=24

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


def image():
	file = sys.argv[1]
	data = np.loadtxt(file).transpose()
	plt.matshow(data)

def single_plot():
	file = sys.argv[1]
	data = np.loadtxt(file).transpose()
	#data = data.transpose()
	plt.scatter(data[0],data[1])
	#plt.scatter(data[0],data[1])
	plt.show()



def heat_capacity(folder,n=1,spins=True,norm=size**3):
	#Input files should be: energy.txt energy2.txt
	#f1 = sys.argv[n]
	#f2 = sys.argv[m]
	means = np.loadtxt(str(folder)+"/result"+str(n)+"/energy.txt").transpose()
	#means = np.loadtxt("high_res/HL_flux/result"+str(n)+"/energy.txt").transpose()
	squared_means = np.loadtxt(str(folder)+"/result"+str(n)+"/energy2.txt").transpose()
	#heat_capacity = (squared_means[1]-(means[1]*means[1]))/((means[0]*means[0]*norm))
	#plt.plot(means[0],heat_capacity,label="Variance")
	heat_capacity = -40*(np.gradient(means[1])/norm)
	if n=="CH":
		plt.scatter(means[0],heat_capacity,label="Classical")	
	else:
		plt.plot(means[0],heat_capacity,label="u0 = "+str(n),linewidth=2)
	if (n=="10") and spins:
	#	print("hello")
		plt.scatter([means[0][140],means[0][160],means[0][180]],[heat_capacity[140],heat_capacity[160],heat_capacity[180]],label="Spin structures",color="orange",linewidth=2,marker="x",s=100)
		#plt.scatter([means[0][110],means[0][127],means[0][160]],[heat_capacity[110],heat_capacity[127],heat_capacity[160]],label="Spin structures",color="orange")
	#peaks,_ = sp.find_peaks(d,height=1,distance=10)
	#plt.plot(means[0][peaks],d[peaks],"x",label="Peak")
	
	#plt.show()


def heat_capacity2(n=1,m=2,norm=size**3):
	#Input files should be: energy.txt energy2.txt
	f1 = sys.argv[n]
	f2 = sys.argv[m]
	means = np.loadtxt(f1).transpose()
	squared_means = np.loadtxt(f2).transpose()
	heat_capacity = (squared_means[1]-(means[1]*means[1]))/((means[0]*means[0]*norm))
	#plt.plot(means[0],heat_capacity,label="Variance")
	d = -40*(np.gradient(means[1])/norm)
	plt.plot(means[0],d,label="Numerical derivative")
	plt.plot(means[0],heat_capacity,label="Variance")
	#peaks,_ = sp.find_peaks(d,height=1,distance=10)
	#plt.plot(means[0][peaks],d[peaks],"x",label="Peak")
	plt.xlabel("Temperature")
	plt.ylabel("Heat Capacity")
	plt.title("Comparison of methods (J=1,D=0.75,u0=20)")
	plt.legend()
	plt.show()


def difference_plot():
	#Plots difference of 2 files
	f1 = sys.argv[1]
	f2 = sys.argv[2]
	d1 = np.loadtxt(f1).transpose()
	d2 = np.loadtxt(f2).transpose()
	plt.plot(d1-d2)
	#plt.show()

def plot_from_csv():
	f1 = sys.argv[3]
	d = np.genfromtxt(f1,delimiter=',').transpose()
	plt.plot(d[0],d[1],label="miyatake")

def plot_spin():
	f1 = sys.argv[1]
	d1 = np.loadtxt(f1).transpose()
	plt.plot(d1[0],np.sqrt(d1[1]**2+d1[2]**2+d1[3]**2)/size**3)
	
	#plt.plot(d1[0],d1[1]/size**3)
	#plt.plot(d1[0],d1[2]/size**3)
	#plt.plot(d1[0],d1[3]/size**3)


def final_spins():
	#plots theta and phi components of final spin vectors to check
	#they are uncorrelated
	f1 = sys.argv[1]
	d1 = np.loadtxt(f1).transpose()
	theta = np.arccos(d1[2]/size**3)
	phi = np.arctan2(d1[0],d1[1])
	plt.scatter(theta,phi)

def hist():
	f1 = sys.argv[1]
	d1 = np.loadtxt(f1)
	plt.hist(d1,bins=50)


def peak_find(n,folder):
	#Loads file from "HL_fluxuations/result#n/energy.txt"
	means = np.loadtxt(str(folder)+"/result"+str(n)+"/energy.txt").transpose()
	squared_means = np.loadtxt(str(folder)+"/result"+str(n)+"/energy2.txt").transpose()
	#heat_capacity = (squared_means[1]-(means[1]*means[1]))/((means[0]*means[0]*size**3))
	heat_capacity = -40*(np.gradient(means[1])/size**3)
	#peaks,_ = sp.find_peaks(heat_capacity,threshold=0,distance=500)
	peaks = np.argmax(heat_capacity)
	#print(len(peaks))
	return means[0][peaks]
	#Finds location of peak in data array. Returns location of peak



def u0_tc(folder):

	#u0_str = np.array(["1","2","3","4","5","6","7","8","9","10","12","14","16","18","20","30","40","50","60","70","80","90","100","200","300","400","500"])
	u0_str = np.array(["1","10","30","50","100"])
	u0 = np.array([1,10,30,50,100])
	
	#u0 = np.array([1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100,200,300,400,500])

	ch_line = peak_find("CH")
	plt.hlines(ch_line,np.min(u0),np.max(u0))
	v = np.vectorize(lambda x:peak_find(x,folder))
	#plt.semilogx(u0,(v(u0_str)))
	plt.scatter(u0,(v(u0_str)))
	plt.xscale('log')

	plt.title("Heisenberg Landau Tc vs u0 (J=1,D=0)")
	plt.ylabel("Transition Temperature")
	plt.xlabel("u0")
	plt.show()


def multi_heat_capacity(folder,spins):

	#u0_str = np.array(["1","2","3","4","5","6","7","8","9","10","12","14","16","18","20","30","40","50","60","70","80","90","100","200","300","400","500"])
	u0_str = np.array(["1","10","500","CH"])
	
	for u in u0_str:
		heat_capacity(folder,u,spins)
	plt.title("Heisenberg-Landau Heat capacity (J=1,D=0)")
	plt.ylabel("Heat capacity")
	plt.xlabel("Temperature")
	plt.legend()
	plt.show()


def multi_heat_capacity2():
	folder1 = "HL_DM_flux5"
	folder2 = "HL_flux4"
	norm=size**3
	fig = plt.figure()
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(223)
	ax3 = fig.add_subplot(222)
	ax4 = fig.add_subplot(224)
	u0_str = np.array(["1","2","10","CH"])
	u0 = np.array([1,2,3,4,5,6,7,8,9,10,12,16,18,20,30,40,50,100])
	
	for n in u0_str:
		means = np.loadtxt(str(folder1)+"/result"+str(n)+"/energy.txt").transpose()
		squared_means = np.loadtxt(str(folder1)+"/result"+str(n)+"/energy2.txt").transpose()
		#heat_capacity = (squared_means[1]-(means[1]*means[1]))/((means[0]*means[0]*norm))
		heat_capacity = -40*(np.gradient(means[1])/norm)
		if n=="CH":
			ax1.scatter(means[0],heat_capacity,label="Classical")	
		else:
			ax1.plot(means[0],heat_capacity,label="u0 = "+str(n),linewidth=2)
		if n=="10":
			ax1.scatter([means[0][140],means[0][159],means[0][180]],[heat_capacity[140],heat_capacity[158],heat_capacity[180]],label="Spin structures",color="orange",linewidth=2,marker="x",s=100)
			#plt.scatter([means[0][110],means[0][127],means[0][160]],[heat_capacity[110],heat_capacity[127],heat_capacity[160]],label="Spin structures",color="orange")
		

		means = np.loadtxt(str(folder2)+"/result"+str(n)+"/energy.txt").transpose()
		squared_means = np.loadtxt(str(folder2)+"/result"+str(n)+"/energy2.txt").transpose()
		#heat_capacity = (squared_means[1]-(means[1]*means[1]))/((means[0]*means[0]*norm))
		heat_capacity = -40*(np.gradient(means[1])/norm)
		if n=="CH":
			ax3.scatter(means[0],heat_capacity,label="Classical")	
		else:
			ax3.plot(means[0],heat_capacity,label="u0 = "+str(n),linewidth=2)
		
	#u0 = np.array([1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100,200,300,400,500])

	ch_line = peak_find("CH",folder1)
	ax2.hlines(ch_line,np.min(u0)/10,np.max(u0)*10,label="Classical limit")
	v = np.vectorize(lambda x:peak_find(x,folder1))
	ax2.scatter(u0,(v(u0)))
	ax2.set_xscale('log')


	ch_line = peak_find("CH",folder2)
	ax4.hlines(ch_line,np.min(u0)/10,np.max(u0)*10,label="Classical limit")
	v = np.vectorize(lambda x:peak_find(x,folder2))
	ax4.scatter(u0,(v(u0)))
	ax4.set_xscale('log')
	#plt.title("Heisenberg Landau Tc vs u0 (J=1,D=0)")
	#plt.ylabel("Transition Temperature")
	#plt.xlabel("u0")
	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()

	ax1.set_title("J=1,D=0.75",size=30)
	ax1.set_ylabel("Heat Capacity",size=20)
	ax1.set_xlabel("Temperature",size=20)
	
	ax3.set_title("J=1,D=0",size=30)
	#ax3.set_ylabel("Heat Capacity",size=20)
	ax3.set_xlabel("Temperature",size=20)

	ax2.set_ylabel("Transition Temperature",size=20)
	ax2.set_xlabel("u0",size=20)

	#ax4.set_ylabel("Transition Temperature")
	ax4.set_xlabel("u0",size=20)
	plt.show()
	
	#plt.show()
def potential():
	#u0=1
	for u0 in [0.5,1,10,100]:
		f = np.vectorize(lambda x:u0*(1-2*x**2+x**4))
		xs = np.linspace(0,2,1000)
		ys = f(xs)
		plt.plot(xs,ys,label="u0="+str(u0))

	plt.legend(loc=1)

	plt.xlabel("m")
	plt.ylabel("u(m)")
	plt.show()


def main():
	#print(scipy.__version__)
	#image()
	#heat_capacity2()
	#hist()
	#single_plot()

	#heat_capacity2()
	#single_plot()
	#plt.show()


	"""
	u0_str = np.array(["1","2","3","4","5","6","7","8","9","10","12","14","16","18","20","30","40","50","60","70","80","90","100","200","300","400","500"])
	#u0_str = np.array(["1","2","3","10","50","100","1000"])
	
	u0 = np.array([1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100,200,300,400,500])
	v = np.vectorize(peak_find)
	#plt.semilogx(u0,(v(u0_str)))
	plt.scatter(u0,(v(u0_str)))
	plt.xscale('log')
	#heat_capacity(1,2)
	#for u in u0_str:

	#	heat_capacity(u)
	#plt.title("Heisenberg-Landau Heat capacity (J=0,D=1)")
	#plt.ylabel("Heat capacity")
	#plt.xlabel("Temperature")
	plt.title("Heisenberg Landau Tc vs u0 (J=0,D=1)")
	plt.ylabel("Transition Temperature")
	plt.xlabel("u0")
	#plt.legend()
	plt.show()
	"""
	#u0_tc()
	#multi_heat_capacity("HL_DM_flux5",True)
	#multi_heat_capacity2()
	#plot_spin()
	#plt.show()
	heat_capacity2(1,2)
	#hist()
	plt.show()
	#potential()
	#plt.show()
	#heat_capacity(3,4)
	#heat_capacity(5,6)
	#heat_capacity(7,8)
	#final_spins()
	#plot_spin()
	#plot_from_csv()
	#difference_plot()
	
main()

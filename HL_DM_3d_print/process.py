import numpy as np
import sys



size=12
no_of_temps=49

def process_spin_data(filename):
	#load spins
	data = np.loadtxt(filename)
	print(data.shape)
	#remove temperatuers written to same file
	mask = np.ones(len(data),dtype=bool)
	indices = [x*size*size*size*3+x for x in range(no_of_temps)]
	mask[indices] = False
	#reshape
	data = np.reshape(data[mask],(no_of_temps,size,size,size,3))
	return data

def main():
	c = 180/np.pi
	t=45
	spins = process_spin_data("spins.txt")[t]
	mags = np.linalg.norm(spins,axis=3)
	phi = c*np.arctan(spins[:,:,:,2],spins[:,:,:,1]).reshape(1,spins.shape[0]*spins.shape[1]*spins.shape[2])
	theta = c*np.arccos(spins[:,:,:,0]/mags).reshape(1,spins.shape[0]*spins.shape[1]*spins.shape[2])
	mags = mags.reshape(1,spins.shape[0]*spins.shape[1]*spins.shape[2])
	print(spins.shape)
	print(mags.shape)
	print(phi.shape)
	print(theta.shape)

	np.savetxt('mags_out.txt',mags,delimiter=',',fmt="%5.5f")
	np.savetxt('phi_out.txt',phi,delimiter=',',fmt="%0.1i")
	np.savetxt('theta_out.txt',theta,delimiter=',',fmt="%0.1i")


main()
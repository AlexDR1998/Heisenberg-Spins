import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import sys




size = 24
no_of_temps = 98

#size = 30
#no_of_temps = 9

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

def fft_2d_sum(data,ax=0):
	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	#log might scale the image better
	#q_sum = np.log(np.sum(q,axis=ax))
	q_sum = np.sum(q,axis=ax)
	plt.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#cb = fig.colorbar(im)
	plt.xlabel("qx")
	plt.ylabel("qy")
	plt.show()

def fft_2d_sum_sub(data,ax=0):
	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	#log might scale the image better
	#q_sum = np.log(np.sum(q,axis=ax))
	q_sum = np.sum(q,axis=ax)

	#find peak locations
	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)

	#see how data at peak varies along summation axis
	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	if ax==0:
		xs1 = q[:,peak1_location[0],peak1_location[1]]
		xs2 = q[:,peak2_location[0],peak2_location[1]]
		print(np.sum(xs1))
	elif ax==1:
		xs1 = q[peak1_location[0],:,peak1_location[1]]
		xs2 = q[peak2_location[0],:,peak2_location[1]]
		print(np.sum(xs1))
	else:
		xs1 = q[peak1_location[0],peak1_location[1],:]
		xs2 = q[peak2_location[0],peak2_location[1],:]
		print(np.sum(xs1))
	
	print(peak1_location-np.array([int(size/2),int(size/2)]))
	print(peak2_location-np.array([int(size/2),int(size/2)]))

	ax1.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#cb = fig.colorbar(im)
	ax1.set_xlabel("qx")
	ax1.set_ylabel("qy")
	#plt.colorbar()
	#ax2.semilogy(range(-10,10),xs1)
	#ax2.semilogy(range(-10,10),xs2)
	ax2.plot(range(-int(size/2),int(size/2)),xs1)
	ax2.plot(range(-int(size/2),int(size/2)),xs2)

	ax2.set_xlabel("qz")
	ax1.set_title("J=1,D=0.75,u0=30")
	plt.show()



def fft_2d_sum_animate(data,ax=0):
	#projects vectors from a 2d slice onto total spin vector at given temperature
	def update(i):
		#screendata = data[i]
		sx = data[i,:,:,:,0]
		sy = data[i,:,:,:,1]
		sz = data[i,:,:,:,2]
		qx = (np.fft.fftshift(np.fft.fftn(sx)))
		qy = (np.fft.fftshift(np.fft.fftn(sy)))
		qz = (np.fft.fftshift(np.fft.fftn(sz)))
		q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
		screendata = np.log(np.sum(q,axis=0))
		#screendata = data[i,:,:,c,0]*axis[i,0]+data[i,:,:,c,1]*axis[i,1]+data[i,:,:,c,2]*axis[i,2]
		matrix.set_array(screendata)
	
	#screendata = data[0,:,:,c,0]*axis[0,0]+data[0,:,:,c,1]*axis[0,1]+data[0,:,:,c,2]*axis[0,2]
	sx = data[-1,:,:,:,0]
	sy = data[-1,:,:,:,1]
	sz = data[-1,:,:,:,2]
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	screendata = np.log(np.sum(q,axis=0))

	#data_flat = screendata[:,:,c,0]*axis[0]+screendata[:,:,c,1]*axis[1]+screendata[:,:,c,2]*axis[2]
	fig, ax = plt.subplots()            

	#print(data_flat.shape)
	matrix = ax.matshow(screendata,cmap="nipy_spectral")

	#plt.colorbar(matrix)
	ani = animation.FuncAnimation(fig,update,frames=no_of_temps,interval=200)
	plt.show()


def fft_2d_slice(data,c=10):

	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))


	plt.matshow((np.abs(qx[:,:,c]))**2)
	plt.colorbar()
	plt.show()

	plt.matshow((np.abs(qy[:,:,c]))**2)
	plt.colorbar()
	plt.show()

	plt.matshow((np.abs(qz[:,:,c]))**2)
	plt.colorbar()
	plt.show()



def vec_plot(data,thresh=10):
	#Doesn't look very good for large grid size
	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	qx = np.abs(np.fft.fftshift(np.fft.fftn(sx)))
	qy = np.abs(np.fft.fftshift(np.fft.fftn(sy)))
	qz = np.abs(np.fft.fftshift(np.fft.fftn(sz)))
	X,Y,Z = np.meshgrid(np.arange(size),np.arange(size),np.arange(size))

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	#for x in range(10):

	ax.quiver(X,Y,Z,np.where(qx>thresh,np.log(qx),0),np.where(qy>thresh,np.log(qy),0),np.where(qz>thresh,np.log(qz),0),linewidths=1)
	plt.show()


def voxel_plot(data,thresh=20):

	X,Y,Z = np.meshgrid(np.arange(size),np.arange(size),np.arange(size))
	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	qx = np.abs(np.fft.fftshift(np.fft.fftn(sx)))
	qy = np.abs(np.fft.fftshift(np.fft.fftn(sy)))
	qz = np.abs(np.fft.fftshift(np.fft.fftn(sz)))

	#qx = np.abs((np.fft.fftn(sx)))
	#qy = np.abs((np.fft.fftn(sy)))
	#qz = np.abs((np.fft.fftn(sz)))

	ft = (qx>thresh)|(qy>thresh)|(qz>thresh)
	colours = np.zeros(ft.shape +(3,))
	colours[...,0]=np.log(qx)/np.max(np.log(qx))
	colours[...,1]=np.log(qy)/np.max(np.log(qy))
	colours[...,2]=np.log(qz)/np.max(np.log(qz))
	#colours[...,3]=(colours[...,0]+colours[...,1]+colours[...,2])/3
	print(colours.shape)
	#print(qx>thresh)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.voxels(ft,facecolors=colours,edgecolor='k')
	plt.xlabel("qx")
	plt.ylabel("qy")
	plt.show()

def main():
	f1 = sys.argv[1]
	spins = process_spin_data(f1)
	
	#fft_2d_slice(spins[39])
	#fft_2d_slice(spins[48])
	t = int(sys.argv[2])
	fft_2d_sum(spins[t],0)
	#fft_2d_sum(spins[t],1)
	#fft_2d_sum(spins[t],2)
	#fft_2d_sum_animate(spins)
	#a = np.array([[1,2,3],[4,5,3],[9,0,2]])
	#b = np.ravel(a)
	#print(b)
	#c = np.unravel_index(np.argmax(np.ravel(a)),a.shape)
	#print(a[c])
	
	

main()
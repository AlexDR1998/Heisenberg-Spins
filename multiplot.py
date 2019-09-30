import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import sys

no_of_temps = 191
size = 24
c = 15
ax = 0


def flat_vec1(data,c=0):
	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data[c,:,:,0]
	V = data[c,:,:,1]
	col = data[c,:,:,2]
	plt.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#plt.imshow(col,origin="lower",interpolation="none")
	cbar = plt.colorbar()
	cbar.set_label("z",rotation=0)
	
	plt.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	#plt.title("Spin structure (J=1,D=0.75,u0=10)")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.show()

def fft_2d_sum(data,ax=0):
	sx = data[:,:,:,0]
	sy = data[:,:,:,1]
	sz = data[:,:,:,2]

	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2

	q_sum = np.sum(q,axis=ax)
	plt.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	plt.xlabel("qx")
	plt.ylabel("qy")
	plt.show()

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


def plot3():
	filename = sys.argv[1]
	data1 = process_spin_data(filename)[160]
	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)

	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data1[c,:,:,0]
	V = data1[c,:,:,1]
	col = data1[c,:,:,2]
	sx = data1[:,:,:,0]
	sy = data1[:,:,:,1]
	sz = data1[:,:,:,2]
	ax1.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax1.imshow(col,origin="lower",interpolation="none")
	ax1.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	ax1.set_title("Real space spin structure")
	#ax1.set_xlabel("x")
	ax1.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)
	
	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	




	ax2.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#ax2.set_xlabel("qx")
	ax2.set_title("2D projected diffraction")
	ax2.set_ylabel("qy")
	#plt.show()
	ax3.set_title("Z axis diffraction location")
	ax3.plot(range(-int(size/2),int(size/2)),xs1)
	ax3.plot(range(-int(size/2),int(size/2)),xs2)

	plt.show()


def plot6():

	filename = sys.argv[1]
	data1 = process_spin_data(filename)[120]
	data2 = process_spin_data(filename)[159]
	fig = plt.figure()
	ax1 = fig.add_subplot(231)
	ax2 = fig.add_subplot(232)
	ax3 = fig.add_subplot(233)
	ax4 = fig.add_subplot(234)
	ax5 = fig.add_subplot(235)
	ax6 = fig.add_subplot(236)

	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data1[c,:,:,0]
	V = data1[c,:,:,1]
	col = data1[c,:,:,2]
	sx = data1[:,:,:,0]
	sy = data1[:,:,:,1]
	sz = data1[:,:,:,2]
	ax1.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax1.imshow(col,origin="lower",interpolation="none")
	ax1.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	ax1.set_title("Real space spin structure")
	#ax1.set_xlabel("x")
	ax1.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)
	
	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	




	ax2.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#ax2.set_xlabel("qx")
	ax2.set_title("2D projected diffraction")
	ax2.set_ylabel("qy")
	#plt.show()
	ax3.set_title("Z axis diffraction intensity")
	ax3.plot(range(-int(size/2),int(size/2)),xs1)
	#ax3.plot(range(-int(size/2),int(size/2)),xs2)




	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data2[c,:,:,0]
	V = data2[c,:,:,1]
	col = data2[c,:,:,2]
	sx = data2[:,:,:,0]
	sy = data2[:,:,:,1]
	sz = data2[:,:,:,2]
	ax4.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax3.imshow(col,origin="lower",interpolation="none")
	ax4.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	#ax4.set_title("At phase transition")
	ax4.set_xlabel("x")
	ax4.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)

	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	
	ax5.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	ax5.set_xlabel("qx")
	ax5.set_ylabel("qy")
	ax6.plot(range(-int(size/2),int(size/2)),xs1)
	ax6.set_xlabel("qz")
	#ax6.plot(range(-int(size/2),int(size/2)),xs2)

	plt.show()
def plot9():

	filename = sys.argv[1]
	data1 = process_spin_data(filename)[50]#[140]
	data2 = process_spin_data(filename)[115]#[159]
	data3 = process_spin_data(filename)[140]#[180]
	fig = plt.figure()
	ax1 = fig.add_subplot(331)
	ax2 = fig.add_subplot(332)
	ax3 = fig.add_subplot(333)
	ax4 = fig.add_subplot(334)
	ax5 = fig.add_subplot(335)
	ax6 = fig.add_subplot(336)
	ax7 = fig.add_subplot(337)
	ax8 = fig.add_subplot(338)
	ax9 = fig.add_subplot(339)

	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data1[c,:,:,0]
	V = data1[c,:,:,1]
	col = data1[c,:,:,2]
	sx = data1[:,:,:,0]
	sy = data1[:,:,:,1]
	sz = data1[:,:,:,2]
	ax1.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax1.imshow(col,origin="lower",interpolation="none")
	ax1.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	ax1.set_title("Real space spin structure")
	#ax1.set_xlabel("x")
	ax1.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)
	
	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	




	ax2.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#ax2.set_xlabel("qx")
	ax2.set_title("2D projected diffraction")
	ax2.set_ylabel("qy")
	#plt.show()
	ax3.set_title("Z axis diffraction intensity")
	ax3.plot(range(-int(size/2),int(size/2)),xs1)
	ax3.plot(range(-int(size/2),int(size/2)),xs2)




	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data2[c,:,:,0]
	V = data2[c,:,:,1]
	col = data2[c,:,:,2]
	sx = data2[:,:,:,0]
	sy = data2[:,:,:,1]
	sz = data2[:,:,:,2]
	ax4.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax3.imshow(col,origin="lower",interpolation="none")
	ax4.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	#ax4.set_title("At phase transition")
	#ax4.set_xlabel("x")
	ax4.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)

	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	
	ax5.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	#ax5.set_xlabel("qx")
	ax5.set_ylabel("qy")
	ax6.plot(range(-int(size/2),int(size/2)),xs1)
	ax6.plot(range(-int(size/2),int(size/2)),xs2)



	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data3[c,:,:,0]
	V = data3[c,:,:,1]
	col = data3[c,:,:,2]
	sx = data3[:,:,:,0]
	sy = data3[:,:,:,1]
	sz = data3[:,:,:,2]
	ax7.contourf(X,Y,col,levels=np.linspace(np.min(col),np.max(col),100))
	#ax5.imshow(col,origin="lower",interpolation="none")
	ax7.quiver(X,Y,U,V,linewidths=1,pivot="middle")
	#ax7.set_title("Below phase transition")
	ax7.set_xlabel("x")
	ax7.set_ylabel("y")
	#fft stuff
	qx = (np.fft.fftshift(np.fft.fftn(sx)))
	qy = (np.fft.fftshift(np.fft.fftn(sy)))
	qz = (np.fft.fftshift(np.fft.fftn(sz)))
	q = np.abs(qx)**2+np.abs(qy)**2+np.abs(qz)**2
	q_sum = np.sum(q,axis=ax)

	peak1_location = np.unravel_index(np.argmax(np.ravel(q_sum)),q_sum.shape)
	q_copy = np.copy(q_sum)
	q_copy[peak1_location] = 0
	peak2_location = np.unravel_index(np.argmax(np.ravel(q_copy)),q_copy.shape)
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
	
	ax8.imshow(q_sum,cmap="magma",extent=[-int(size/2),int(size/2),-int(size/2),int(size/2)],interpolation='none')
	ax8.set_xlabel("qx")
	ax8.set_ylabel("qy")
	ax9.plot(range(-int(size/2),int(size/2)),xs1)
	ax9.plot(range(-int(size/2),int(size/2)),xs2)
	ax9.set_xlabel("qz")



	plt.show()



def main():
	plot9()
main()
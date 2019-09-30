import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import sys



"""
Code to visualise spin structures from HL_MC or CH_MC simulation runs
"""

size = 12
no_of_temps = 49


#size = 30
#no_of_temps = 9
def vec_plot(data):
	
	X,Y,Z = np.meshgrid(np.arange(size),np.arange(size),np.arange(size))
	U = data[:,:,:,0]
	V = data[:,:,:,1]
	W = data[:,:,:,2]

	#c = np.sqrt(U**2+V**2+W**2)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	#for x in range(10):

	ax.quiver(X,Y,Z,U,V,W,linewidths=1)
	plt.show()

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

def flat_vec2(data,c=0):
	X,Y,Z = np.meshgrid(np.arange(size),np.arange(size),np.arange(2))
	U = data[:,:,:2,0]
	V = data[:,:,:2,1]
	W = data[:,:,:2,2]

	#c = np.sqrt(U**2+V**2+W**2)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	#ax.autoscale(enable=False,axis="z")
	#for x in range(10):
	ax.quiver(X,Y,Z,U,V,W,pivot="middle")
	plt.show()

def flat_vec3(data,c=0):
	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data[:,c,:,0]
	V = data[:,c,:,1]
	col = data[:,c,:,2]
	
	plt.quiver(X,Y,U,V,col,linewidths=0.1,pivot="middle",headwidth=5,headlength=5)
	plt.title("Spin structures with DM interactions")
	plt.xlabel("x")
	plt.ylabel("y")
	cbar = plt.colorbar()
	cbar.set_label("z axis",rotation=0)
	plt.show()

def flat_plot(data,axis,c=0):
	#projects onto axis

	
	data_flat = data[:,:,c,0]*axis[0]+data[:,:,c,1]*axis[1]+data[:,:,c,2]*axis[2]
	plt.matshow(data_flat)
	plt.show()

def stream(data,c=0):
	X,Y = np.meshgrid(np.arange(size),np.arange(size))
	U = data[:,:,c,0]
	V = data[:,:,c,1]
	#col = np.sqrt(U**2+V**2)
	col = data[:,:,c,2]
	plt.streamplot(X,Y,U,V,color=col)
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

def flat_animate1(data,axis,c=0):
	#projects vectors from a 2d slice onto total spin vector at given temperature
	def update(i):
		#screendata = data[i]
		screendata = data[i,:,:,c,0]*axis[i,0]+data[i,:,:,c,1]*axis[i,1]+data[i,:,:,c,2]*axis[i,2]
		matrix.set_array(screendata)
	
	screendata = data[0,:,:,c,0]*axis[0,0]+data[0,:,:,c,1]*axis[0,1]+data[0,:,:,c,2]*axis[0,2]
	#data_flat = screendata[:,:,c,0]*axis[0]+screendata[:,:,c,1]*axis[1]+screendata[:,:,c,2]*axis[2]
	fig, ax = plt.subplots()            

	#print(data_flat.shape)
	matrix = ax.matshow(screendata)

	#plt.colorbar(matrix)
	ani = animation.FuncAnimation(fig,update,frames=no_of_temps,interval=200)
	plt.show()

def animated_slices(data,axis):
	def update(i):
		#screendata = data[i]
		screendata = data[:,:,i,0]*axis[0]+data[:,:,i,1]*axis[1]+data[:,:,i,2]*axis[2]
		matrix.set_array(screendata)
	
	screendata = data[:,:,0,0]*axis[0]+data[:,:,0,1]*axis[1]+data[:,:,0,2]*axis[2]
	#data_flat = screendata[:,:,c,0]*axis[0]+screendata[:,:,c,1]*axis[1]+screendata[:,:,c,2]*axis[2]
	fig, ax = plt.subplots()            

	#print(data_flat.shape)
	matrix = ax.matshow(screendata)

	#plt.colorbar(matrix)
	ani = animation.FuncAnimation(fig,update,frames=size,interval=200)
	plt.show()

def flat_animate2(data,axis,c=0):
	#projects vectors from a 2d slice onto specific axis
	def update(i):
		#screendata = data[i]
		screendata = data[i,:,:,c,0]*axis[0]+data[i,:,:,c,1]*axis[1]+data[i,:,:,c,2]*axis[2]
		matrix.set_array(screendata)
	
	screendata = data[0,:,:,c,0]*axis[0]+data[0,:,:,c,1]*axis[1]+data[0,:,:,c,2]*axis[2]
	#data_flat = screendata[:,:,c,0]*axis[0]+screendata[:,:,c,1]*axis[1]+screendata[:,:,c,2]*axis[2]
	fig, ax = plt.subplots()            

	#print(data_flat.shape)
	matrix = ax.matshow(screendata)

	#plt.colorbar(matrix)
	ani = animation.FuncAnimation(fig,update,frames=no_of_temps,interval=200)
	plt.show()

def complex_array_to_rgb(X, theme='dark', rmax=None):
	'''
    Maps array of complex numbers to colours. Taken from stack overflow:
	https://stackoverflow.com/questions/15207255/is-there-any-way-to-use-bivariate-colormaps-in-matplotlib
    '''
	absmax = rmax or np.abs(X).max()
	Y = np.zeros(X.shape + (3,), dtype='float')
	Y[..., 0] = np.angle(X)/(2 * np.pi) % 1
	if theme == 'light':
		Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
		Y[..., 2] = 1
	elif theme == 'dark':
		Y[..., 1] = 1
		Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
	Y = matplotlib.colors.hsv_to_rgb(Y)
	return Y


def angular_animate1(data,axis,c=0):
	#Like flat_animate1 but maps theta and phi angles of spins from total spin to colours
	def update(i):
		#screendata = data[i]
		#screendata = data[i,:,:,c,0]*axis[i,0]+data[i,:,:,c,1]*axis[i,1]+data[i,:,:,c,2]*axis[i,2]
		screendata = complex_array_to_rgb(data[i,:,:,c,0]*1j + data[i,:,:,c,1])
		matrix.set_array(screendata)
	
	#screendata = data[0,:,:,c,0]*axis[0,0]+data[0,:,:,c,1]*axis[0,1]+data[0,:,:,c,2]*axis[0,2]
	screendata = complex_array_to_rgb(data[0,:,:,c,0]*1j + data[0,:,:,c,1])
	#data_flat = screendata[:,:,c,0]*axis[0]+screendata[:,:,c,1]*axis[1]+screendata[:,:,c,2]*axis[2]
	fig, ax = plt.subplots()            

	#print(data_flat.shape)
	matrix = ax.imshow(screendata)

	#plt.colorbar(matrix)
	ani = animation.FuncAnimation(fig,update,frames=no_of_temps,interval=200)
	plt.show()




def vector_demo():
	fmx = np.array([0,0])
	fmy = np.array([1,1])

	fmdmx = np.array([1,0])
	fmdmy = np.array([0,1])

	
	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	X,Y = np.meshgrid(np.arange(2),np.arange(1))
	ax1.quiver(X,Y,fmx,fmy,scale=0.005,scale_units='dots',pivot='mid',width=0.01)
	ax2.quiver(X,Y,fmdmx,fmdmy,scale=0.005,scale_units='dots',pivot='mid',width=0.01)
	ax1.set_title("J minimising",size=30)
	ax2.set_title("D minimising",size=30)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	ax2.axes.get_xaxis().set_visible(False)
	ax2.axes.get_yaxis().set_visible(False)
	plt.show()

def main():
	#filename = sys.argv[1]
	#data = process_spin_data(filename)
	#print(data.shape)
	#t = int(sys.argv[2])
	#axis=[1,0,0]
	#fn = sys.argv[2]
	#axis = np.loadtxt(fn)
	#t=88
	#flat_plot(data[40],[1,0,0])
	#axis = (axis/size**3)[:,1:]
	#animated_slices(data[t],[0,1,0])
	#flat_animate2(data,[1,0,0],0)
	
	#axis = (axis/size**3)[-1,1:]
	#flat_animate2(data,axis,0)
	#angular_animate1(data,axis,0)
	

	#stream(data[t])
	#for t in [87,88,89,90,91,92,93,94,95,96,97,98]:
	#for c in range(12):
	#	flat_vec1(data[t],c)
	vector_demo()



	#flat_vec3(data[t],0)
	#flat_animate2(data,axis)
	#stream(data[0])
	#vec_plot(data[0])
	#print(axis)
	#vec_plot(data[0])
	#for c in range(size):
	#flat_vec(data[40])
	#for c in range(10):
	#	flat_plot(data,axis,c)
	#flat_plot(data,[0,1,0])
	#flat_plot(data,[0,0,1])
main()
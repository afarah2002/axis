import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def read_data(file):
	'''
	Arg: txt file as str
	Function: Takes in file, returns list of data lists 
	File format: time is final string, makes n-sized
				 list of n-elements before 
	'''
	xs = []
	ys = []
	zs = []

	with open(file, 'r') as file:
		full_data = []
		for line in file:
			instance = []
			for element in line.split():
				instance.append(float(element))
			full_data.append(instance)
		return np.array(full_data).flatten()

def plotter(densities, energies, SRList):

	# print(self.stoppingRanges)
	# print(densities)

	xs = densities
	ys = energies
	zs = SRList

	# xs = [0.001, 0.011, 0.020999999999999998, 0.030999999999999996, 0.040999999999999995, 0.05099999999999999, 0.06099999999999999, 0.071, 0.08099999999999999, 0.09099999999999998]
	# ys = [0.1, 4.6 ,9.1]
	# zs = [[0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

	print(len(xs),len(ys), len(zs[0]))

	xlabel = r'$\rho$' + "(g/cm" r'$^3$' + ")"
	ylabel = "Alpha energy (MeV)"
	zlabel = "Stopping Range (mm)"
	# plt.xlabel(xlabel)
	# plt.ylabel(ylabel)
	# plt.title(ylabel + " vs " + xlabel)
	# plt.grid()
	# plt.scatter(densities, self.stoppingRanges)
	# plt.show()

	fig = plt.figure()
	# first subplot: a 3D scatter plot of positions
	ax = fig.add_subplot(111, projection='3d')

	axes = plt.gca()
	axes.set_xlim([min(xs),max(xs)])
	axes.set_ylim([min(ys),max(ys)])
	axes.set_zlim([0,20])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_zlabel(zlabel)

	# print(xs)
	# print(ys)
	# print(zs)

	for i in np.arange(0, len(ys)): # energies is on the y
		ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])



	plt.title(zlabel + " vs " + xlabel + "vs" + ylabel)
	plt.grid()
	plt.legend(("Po210","Po209", "Po208", "Pu238", "Cm244", "Am241"),loc=6)
	# plt.draw() 
	plt.show()

if __name__ == '__main__':
	energyRange = [5.3,4.9,5.1,5.6,5.8,5.5]
	# energyRange = [5.3,4.9]
	densityRange = list(np.arange(0.01000001,1., .001)) # g/cm3
	# energy = 5.3
	SRList = []
	for energy in energyRange:
		file = "data_01/"+str(energy)+".txt"
		stopping_ranges = read_data(file)
		SRList.append(stopping_ranges)
		print stopping_ranges
	plotter(densityRange, energyRange, SRList)

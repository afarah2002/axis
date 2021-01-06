import matplotlib as mpl
from matplotlib import rc
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.linear_model import LinearRegression
plt.rcParams["font.family"] = "Times New Roman"
mpl.rcParams.update({'font.size': 14})
rc('text', usetex=True)

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


prop_list = ["cesium", "bismuth","mercury","xenon","iodine"]
for prop in prop_list:
	data_file = prop + "_data/mean_ionization.txt"
	ionization = read_data(data_file)[0] 
	print prop, ionization
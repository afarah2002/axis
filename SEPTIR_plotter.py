#!/usr/bin/python2
#-----------------imports and setup-----------------#
import itertools
import matplotlib as mpl
from matplotlib import rc
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
# plt.rcParams["font.family"] = "Times New Roman"
# mpl.rcParams.update({'font.size': 14})
# rc('text', usetex=True)

#-----------------classes-----------------#
class Constants(object):
	#-----------------CONSTANTS USED IN SCRIPT-----------------#
	#-----------------assumes 100% ionization in chamber-----------------#
	pct_chamber_ionized = 1
	#-----------------num of alphas used in G4 sim-----------------#
	num_alphas = 50
	#-----------------chmaber volume from prelim CAD model, 1P unit-----------------#
	Chamber_volume = pct_chamber_ionized*(np.pi*2**2)*4 # cm3
	#-----------------mission duration, 20 years in seconds-----------------#
	duration = np.arange(0,10*3.1536e7,86400) # s
	# #-----------------travel time, years to seconds-----------------#
	# travel_time = 5*3.1536e7
	#-----------------initial radioisotope mass (g), proportional to propellant mass flow rate-----------------#
	RI_m0 = .1
	#-----------------spacecraft initial mass (g)-----------------#
	S_m0 = 750. 
	#-----------------zeta, propellant mass fraction-----------------#
	zeta = .3
	#-----------------avogadro's number (atoms/mol)-----------------#
	N_A = 6.02e23 
	#-----------------first/mean ionization efficiencies-----------------#
	first_or_mean = "mean"
	#-----------------plot bounds-----------------#
	xbounds = [0,max(duration)/3.1536e7] # years
	ybounds = [0,3] # mg/s
	ybounds_Bool = False 
	#-----------------Isp (s)-----------------#
	Isp_BIT1 = 1600
	Isp_BIT3 = 2100
	Isp_BIT7 = 3300
	Isp_Tigris = 3048
	#-----------------earth gravitational acc (m/s)-----------------#
	g0 = 9.8 
	#-----------------plotted data-----------------#
	plotted_data_list = ["mdot", "Fthrust"]
	plotted_data = "Fthrust"
	#-----------------saved file directory-----------------#
	saved_file_directory = "/home/nasa01/Documents/howe_internship/axis/AXIS_Report_II_Plots/"
	save_bool = False

class Radioisotope(object):

	def __init__(self, name, atomic_number, molar_mass, half_life, specific_activity, alpha_energy):
		#-----------------properties of a radioisotope-----------------#
		self.name = name
		self.Z = atomic_number
		self.M = molar_mass
		self.HL = half_life*3.1536e7 # ** years to seconds
		self.SA = specific_activity
		self.energy = alpha_energy

class Propellant(object):

	def __init__(self, name, atomic_number, molar_mass):
		self.name = name
		self.Z = atomic_number
		self.M = molar_mass

class Calculator(object):

	def __init__(self):
		pass

	def calc(self, prop_molar_mass, max_ionizations_per_alpha, spec_rad, half_life, Isp):
		#-----------------common throughout calcs, declared here-----------------#
		coefficient = prop_molar_mass*max_ionizations_per_alpha*spec_rad*Constants().RI_m0/Constants().N_A

		#-----------------propellant mass flow rate (g/s)-----------------#
		mdot = np.multiply(coefficient, np.exp(np.multiply(Constants().duration, np.float(-0.693/half_life))))

		#-----------------total propellant mass that can be ionized and ejected (g)-----------------#
		mTotal = -coefficient*half_life*(np.exp(-0.693*max(Constants().duration)/half_life)-1)/0.693

		#-----------------time-dependent thrust (N)-----------------#
		Fthrust = mdot*Isp*Constants().g0/1000 # convert mdot from g/s to kg/s

		#-----------------delta v calc from propellant mass-----------------#
		del_v = Isp*Constants().g0*np.log(Constants().S_m0/(Constants().S_m0 - Constants().S_m0*Constants().zeta))
		
		#-----------------time when prop runs out (years)-----------------#
		t_empty = -half_life/0.693*np.log(1 - (Constants().S_m0*Constants().zeta*0.693/(half_life*coefficient)))*1.15741e-5

		#-----------------travel time (s)----------------#
		travel_time = max(Constants().duration) - t_empty

		#-----------------extra RI mass (g)-----------------#
		RI_m_extra = Constants().RI_m0/(np.exp(-0.693*travel_time/half_life)) - Constants().RI_m0

		results_dict = dict(zip(Constants().plotted_data_list,
								[mdot, Fthrust]))


		# print "Time to empty ", t_empty, " days"
		# print "del v ", del_v, " m/s"
		# print "Extra RI mass", RI_m_extra, " g"

		return results_dict[Constants().plotted_data]

	def calc2(self, prop_molar_mass, 
			   max_ionizations_per_alpha, 
			   spec_rad, 
			   half_life, 
			   S_m0, 
			   RI_mT, #RI mass at launch
			   Isp, 
			   zeta, 
			   duration, 
			   plotted_data,
			   transit_time):


		#-----------------RI_m0, RI mass at burn-----------------#
		RI_m0 = RI_mT*np.exp(-0.693*transit_time/half_life) 

		#-----------------common throughout calcs, declared here-----------------#
		coefficient = prop_molar_mass*max_ionizations_per_alpha*spec_rad*RI_m0/Constants().N_A

		#-----------------propellant mass flow rate (g/s)-----------------#
		mdot = np.multiply(coefficient, np.exp(np.multiply(duration, np.float(-0.693/half_life))))

		#-----------------total propellant mass that can be ionized and ejected (g)-----------------#
		# mTotal = -coefficient*half_life*(np.exp(-0.693*max(duration)/half_life)-1)/0.693

		#-----------------time-dependent thrust (N)-----------------#
		Fthrust = mdot*Isp*Constants().g0/1000 # convert mdot from g/s to kg/s

		#-----------------delta v calc from propellant mass-----------------#
		# del_v = Isp*Constants().g0*np.log(S_m0/(S_m0 - S_m0*zeta))
		
		#-----------------time when prop runs out (years)-----------------#
		# t_empty = -half_life/0.693*np.log(1 - (S_m0*zeta*0.693/(half_life*coefficient)))*1.15741e-5

		#-----------------travel time (s)----------------#
		# travel_time = max(duration) - t_empty

		#-----------------extra RI mass (g)-----------------#
		# RI_m_extra = RI_m0/(np.exp(-0.693*travel_time/half_life)) - RI_m0

		results_dict = dict(zip(Constants().plotted_data_list,
								[mdot, Fthrust]))

		return results_dict[plotted_data]

class MaterialAssignment():

	def __init__(self):
		pass

	def assign(self):
		#-----------------MATERIAL ASSIGNMENT-----------------#
		#-----------------radioisotopes-----------------#
		# *** Half life input is years, output is seconds
		Po209 = Radioisotope("Po209",84,208.9824,109,0.63e12,4.9)
		Po208 = Radioisotope("Po208",84,207.9812,2.9,21.8e12,5.1)
		Am241 = Radioisotope("Am241",95,241.05683,432.2,0.126e12,5.5)
		Pu238 = Radioisotope("Pu238",94,238.04956,97.8,0.643e12,5.6)
		Cm244 = Radioisotope("Cm244",96,244.06275,18.1,3.03e12,5.8)
		#-----------------propellants-----------------#
		cesium = Propellant("cesium",55,132.91)
		bismuth = Propellant("bismuth",83,208.98)
		xenon = Propellant("xenon",54,131.29)
		mercury = Propellant("mercury",80,200.59)
		iodine = Propellant("iodine",53,126.90)
		#-----------------------------------------------------#
		radioisotope_list = [Po209,Po208,Am241,Pu238,Cm244]
		propellants_list = [cesium,bismuth,xenon,mercury,iodine]

		return radioisotope_list, propellants_list

class XYPlotter(object):
	def __init__(self, plot_data_list):
		#-----------------initialize plot-----------------#
		self.fig, self.subplot = plt.subplots()
		#-----------------unpack plot data list-----------------#
		propellant = plot_data_list[0] # propellant instance
		data = np.array(plot_data_list[1:]).reshape(len(plot_data_list)-1,2) # 
		#-----------------duration, plotted in years-----------------#
		x = np.divide(Constants().duration,3.1536e7) # years
		#-----------------data for all RIs plotted with this propellant-----------------#
		RIs = np.take(data,0,axis=1)
		ydata = np.take(data,1,axis=1)

		#-----------------plot each RI data set-----------------#
		for y in ydata:
			y *= 1000 # g/s to mg/s, N to mN
			self.subplot.plot(x,y)

		#-----------------plot parameters-----------------#
		xlabel = "Time (years)"
		if Constants().plotted_data == "mdot":
			ylabel = propellant.name.capitalize() + " " r'$\dot{m}_{P}$'+ " (mg/s)"
		if Constants().plotted_data == "Fthrust":
			ylabel = propellant.name.capitalize() + " " r'${F}_{T}$'+ " (mN)"

		self.subplot.set_xlabel(xlabel)
		self.subplot.set_xlim(Constants().xbounds)
		self.subplot.set_ylabel(ylabel)
		if Constants().ybounds_Bool:
			self.subplot.set_ylim(Constants().ybounds)
		self.subplot.set_title(ylabel + " vs " + xlabel + " (" + Constants().first_or_mean + ")")
		self.subplot.grid()

		self.subplot.legend((RIs[0].name + ", " + str(RIs[0].energy) + " MeV",
							 RIs[1].name + ", " + str(RIs[1].energy) + " MeV",
							 RIs[2].name + ", " + str(RIs[2].energy) + " MeV",
							 RIs[3].name + ", " + str(RIs[3].energy) + " MeV",
							 RIs[4].name + ", " + str(RIs[4].energy) + " MeV"),loc=1)

		#-----------------save figure-----------------#
		if Constants().save_bool:
			file_loc = Constants().saved_file_directory + \
					   propellant.name + "_" + \
					   Constants().plotted_data + "_" + \
					   Constants().first_or_mean
			self.fig.savefig(file_loc,bbox_inches='tight')
			# print file_loc

class XYPlotter2(object):
	def __init__(self, plot_data_list, plotted_data):
		#-----------------initialize plot-----------------#
		self.fig, self.subplot = plt.subplots()
		#-----------------unpack plot data list-----------------#
		propellant = plot_data_list[0] # propellant instance
		data = np.array(plot_data_list[1:]).reshape(len(plot_data_list)-1,2) # 
		#-----------------duration, plotted in years-----------------#
		x = np.divide(Constants().duration,3.1536e7) # years
		#-----------------data for all RIs plotted with this propellant-----------------#
		RIs = np.take(data,0,axis=1)
		ydata = np.take(data,1,axis=1)

		#-----------------plot each RI data set-----------------#
		for y in ydata:
			y *= 1000 # g/s to mg/s, N to mN
			self.subplot.plot(x,y)

		#-----------------plot parameters-----------------#
		xlabel = "Time (years)"
		if Constants().plotted_data == "mdot":
			ylabel = propellant.name.capitalize() + " " r'$\dot{m}_{P}$'+ " (mg/s)"
		if Constants().plotted_data == "Fthrust":
			ylabel = propellant.name.capitalize() + " " r'${F}_{T}$'+ " (mN)"

		self.subplot.set_xlabel(xlabel)
		self.subplot.set_xlim(Constants().xbounds)
		self.subplot.set_ylabel(ylabel)
		if Constants().ybounds_Bool:
			self.subplot.set_ylim(Constants().ybounds)
		self.subplot.set_title(ylabel + " vs " + xlabel + " (" + Constants().first_or_mean + ")")
		self.subplot.grid()

		self.subplot.legend((RIs[0].name + ", " + str(RIs[0].energy) + " MeV",
							 RIs[1].name + ", " + str(RIs[1].energy) + " MeV",
							 RIs[2].name + ", " + str(RIs[2].energy) + " MeV",
							 RIs[3].name + ", " + str(RIs[3].energy) + " MeV",
							 RIs[4].name + ", " + str(RIs[4].energy) + " MeV"),loc=1)

		#-----------------save figure-----------------#
		if Constants().save_bool:
			file_loc = Constants().saved_file_directory + \
					   propellant.name + "_" + \
					   Constants().plotted_data + "_" + \
					   Constants().first_or_mean
			self.fig.savefig(file_loc,bbox_inches='tight')
			# print file_loc

class DataScanner(object):
	#-----------------reads data from files-----------------#
	def __init__(self):
		pass
	def read_data(self, file):
		'''
		Arg: txt file as str
		Function: Takes in file, returns list of data lists 
		File format: time is final string, makes n-sized
					 list of n-elements before 
		'''
		with open(file, 'r') as file:
			full_data = []
			for line in file:
				instance = []
				for element in line.split():
					instance.append(float(element))
				full_data.append(instance)
			return np.array(full_data).flatten()

class DataProcessor(object):
	#-----------------collects data and runs calcs funcs on it-----------------#
	def __init__(self):
		pass

	def process1(self, combination):
		#-----------------RI and prop from combo-----------------#
		radioisotope = combination[0]
		propellant = combination[1]
		#-----------------read file with IE data-----------------#
		if Constants().first_or_mean == "first":
			tag = "_first"
		else:
			tag = ""

		totalIonizNumList = []
		total_Ioniz_file = "data/" + propellant.name + "2/" + \
								str(radioisotope.energy) + \
								tag + "_ionization_num.txt"

		ioniz_nums = DataScanner().read_data(total_Ioniz_file)
		alpha_ioniz_nums = ioniz_nums[::2].flatten() # all even indexed counts
		elec_ioniz_nums = ioniz_nums[1::2].flatten() # all odd indexed counts
		total_ioniz_nums = np.add(alpha_ioniz_nums, elec_ioniz_nums)

		totalIonizNumList.append(total_ioniz_nums)
		#---------------------------------------------------------#
		#-----------------calculate ionizations/alpha-----------------#
		IE = max(total_ioniz_nums)/Constants().num_alphas
		# print "IE: ", IE
		#-----------------calc stuff w/ IE -----------------#
		calc_data = Calculator().calc(propellant.M,
								 IE,
								 radioisotope.SA,
								 radioisotope.HL,
								 Constants().Isp_Tigris)

		# print propellant.name, radioisotope.name, "\n"

		# return IE
		return propellant, [radioisotope, calc_data]

	def process2(self, combination, S_m0, RI_mT, Isp, zeta, duration, plotted_data, transit_time):
		#-----------------RI and prop from combo-----------------#
		radioisotope = combination[0]
		propellant = combination[1]
		#-----------------read file with IE data-----------------#
		if Constants().first_or_mean == "first":
			tag = "_first"
		else:
			tag = ""

		totalIonizNumList = []
		total_Ioniz_file = "data/" + propellant.name + "2/" + \
								str(radioisotope.energy) + \
								tag + "_ionization_num.txt"

		ioniz_nums = DataScanner().read_data(total_Ioniz_file)
		alpha_ioniz_nums = ioniz_nums[::2].flatten() # all even indexed counts
		elec_ioniz_nums = ioniz_nums[1::2].flatten() # all odd indexed counts
		total_ioniz_nums = np.add(alpha_ioniz_nums, elec_ioniz_nums)

		totalIonizNumList.append(total_ioniz_nums)
		#---------------------------------------------------------#
		#-----------------calculate ionizations/alpha-----------------#
		IE = max(total_ioniz_nums)/Constants().num_alphas
		# print "IE: ", IE
		#-----------------calc stuff w/ IE -----------------#
		calc_data = Calculator().calc2(propellant.M,
								 IE,
								 radioisotope.SA,
								 radioisotope.HL,
								 S_m0, 
			   					 RI_mT, 
			    				 Isp, 
			  					 zeta, 
			 					 duration,
			 					 plotted_data,
			 					 transit_time)

		# print propellant.name, radioisotope.name, "\n"

		# return IE
		return calc_data

	def process3(self, combination):
		#-----------------jsut gets the ionization efficiency IE-----------------#

		#-----------------RI and prop from combo-----------------#
		radioisotope = combination[0]
		propellant = combination[1]
		#-----------------read file with IE data-----------------#
		if Constants().first_or_mean == "first":
			tag = "_first"
		else:
			tag = ""

		totalIonizNumList = []
		total_Ioniz_file = "data/" + propellant.name + "2/" + \
								str(radioisotope.energy) + \
								tag + "_ionization_num.txt"

		ioniz_nums = DataScanner().read_data(total_Ioniz_file)
		alpha_ioniz_nums = ioniz_nums[::2].flatten() # all even indexed counts
		elec_ioniz_nums = ioniz_nums[1::2].flatten() # all odd indexed counts
		total_ioniz_nums = np.add(alpha_ioniz_nums, elec_ioniz_nums)

		totalIonizNumList.append(total_ioniz_nums)
		#---------------------------------------------------------#
		#-----------------calculate ionizations/alpha-----------------#
		IE = max(total_ioniz_nums)/Constants().num_alphas
		# print "IE: ", IE
		#-----------------calc stuff w/ IE -----------------#
		# calc_data = Calculator().calc2(propellant.M,
		# 						 IE,
		# 						 radioisotope.SA,
		# 						 radioisotope.HL,
		# 						 S_m0, 
		# 	   					 RI_mT, 
		# 	    				 Isp, 
		# 	  					 zeta, 
		# 	 					 duration,
		# 	 					 plotted_data,
		# 	 					 transit_time)

		# print propellant.name, radioisotope.name, "\n"

		# return IE
		return IE


def main():

	#-----------------generate materials-----------------#
	radioisotope_list, propellants_list = MaterialAssignment().assign()

	#-----------------generate array of RI/prop combos-----------------#
	prop_RI_combos = np.array(list(itertools.product(radioisotope_list,propellants_list)), dtype=object)

	#-----------------apply calculations for each combo-----------------#
	prop_RI_combos_calc = np.array(np.apply_along_axis(DataProcessor().process1, 1, prop_RI_combos))

	#-----------------sort combos alphabetically for easy reading-----------------#
	prop_RI_combos_calc_sorted = np.array(sorted(prop_RI_combos_calc, key=lambda x: x[0])).reshape(5,5,2)

	#-----------------generate plot for each propellant-----------------#
	for prop in prop_RI_combos_calc_sorted:
		propellant = prop[0][0]
		prop_plot_data = [propellant] + list(np.concatenate(prop)[1::2])
		XYPlotter(prop_plot_data)

	#-----------------show plots-----------------#
	# plt.show()

if __name__ == '__main__':
	main()
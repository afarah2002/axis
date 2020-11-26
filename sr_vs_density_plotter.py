import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
plt.rcParams["font.family"] = "Times New Roman"


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

def plotter(densities, RI_densities, energies, SRList):

	# print(self.stoppingRanges)
	# print(densities)

	xs = densities
	ys = energies
	zs = SRList

	# xs = [0.001, 0.011, 0.020999999999999998, 0.030999999999999996, 0.040999999999999995, 0.05099999999999999, 0.06099999999999999, 0.071, 0.08099999999999999, 0.09099999999999998]
	# ys = [0.1, 4.6 ,9.1]
	# zs = [[0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

	print(len(xs),len(ys), len(zs[0]))

	xlabel = propellant + " " + r'$\rho$' + " (g/cm" r'$^3$' + ")"
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
	# ax = fig.add_subplot(111, projection='3d')
	ax = fig.add_subplot(111)

	# RI_masses = np.multiply(RI_densities, np.multiply(4*np.pi/3, np.power(np.subtract(2,zs[i]),3)))
	Chamber_volume = (np.pi*2**2)*4# cm3
	Propellant_density = xs

	decay_num = 1e5

	for i in np.arange(0, len(ys)): # energies is on the y
		# ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])
		alpha_energy = ys[i]
		ax.plot(Propellant_density, zs[i])
		# plot the RI masses vs Propellant density
		# RI_masses = np.multiply(RI_densities[i], np.multiply(4*np.pi/3, np.power(np.subtract(2,np.divide(zs[i],10)),3))))
		# RI_radius = np.subtract(2,np.divide(zs[i],10)) # cm
		# RI_volume = np.multiply(4*np.pi/3, np.power(RI_radius,3)) # cm3
		# RI_mass   = np.multiply(RI_densities[i], RI_volume) # g
		# RI_surface_area = np.multiply(4*np.pi, np.power(RI_radius, 2)) # cm2
		# Propellant_volume = np.subtract(Chamber_volume, RI_volume) #cm3
		# Propellant_mass   = np.multiply(Propellant_density,Propellant_volume) # g
		# Propellant_moles = Propellant_mass/prop_molar_mass # mol
		# Propellant_total_iEnergy = Propellant_moles*iEnergy_T
		# total_alpha = Propellant_total_iEnergy/alpha_energy # total number of alphas need for 100% ionization of the propellant
		# alpha_mass = total_alpha*6.643e-24 # change in RI mass (g) due to 100% ionization (loss of total_alpha)
		
		# ion_num_per_SR = np.divide(total_alpha,zs[i])

		# RI_mass = decay_num*alpha_mass
		# RI_power = np.divide(np.multiply(RI_mass, [141,.4855,18,.56, 2.65,.11][i]),1000) # kW
		# total_mass = np.add(Propellant_mass,RI_mass)
		# RI_power_per_total_mass = np.divide(RI_power,total_mass) # kW/cm2
		# RI_power_per_Propellant_mass = np.divide(RI_power,Propellant_mass,RI_mass) # kW/cm2
		# RI_power_by_Propellant_volume = np.multiply(RI_power,Propellant_volume)
		# RI_power_by_Propellant_mass = np.multiply(RI_power,Propellant_mass)

		# ax2.plot(Propellant_density,alpha_mass)
		# ax2.plot(xs, RI_radius)


	ax.set_xlim([min(xs),max(xs)])
	# axes.set_zlim([0,20])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(zlabel)
	ax.set_title(zlabel + " vs " + xlabel)
	ax.grid()


	# ax2.set_xlabel("Propellant_density (g/cm3)")	
	# ax2.set_ylabel("alpha_mass (g)"	)

	# ax2.grid()
	ax.legend(("Po210, energy = " + str(ys[0]),"Po209, energy = " + str(ys[1]), "Po208, energy = " + str(ys[2]), "Pu238, energy = " + str(ys[3]), "Cm244, energy = " + str(ys[4]), "Am241, energy = " + str(ys[5])),loc=5)

	# plt.draw() 
	plt.show()

def ioniz_plotter(densities, energies, totalIonizNumList):

	xs = densities
	ys = energies
	zs = totalIonizNumList


	# plt.xlabel(xlabel)
	# plt.ylabel(ylabel)
	# plt.title(ylabel + " vs " + xlabel)
	# plt.grid()
	# plt.scatter(densities, self.stoppingRanges)
	# plt.show()

	fig = plt.figure()
	# first subplot: a 3D scatter plot of positions
	# ax = fig.add_subplot(111, projection='3d')
	ax = fig.add_subplot(111)


	# RI_masses = np.multiply(RI_densities, np.multiply(4*np.pi/3, np.power(np.subtract(2,zs[i]),3)))
	pct_chamber_ionized = 1
	Chamber_volume = pct_chamber_ionized*(np.pi*2**2)*4# cm3
	Propellant_densities = xs
	specific_radioactivities = [0.63e12,21.8e12,0.643e12,3.03e12,0.126e12] # Bcq/g Po210, Po209, Po208, Pu238, Cm244, Am241
	# time_elapsed = 1e-9
	num_alphas = 50
	RI_molar_masses = [208.9824,207.9812,238.04956,244.06275,241.05683] # g/mol

	num_full_ionizations = 1e12 # number of full-chamber ionizations the RI should be capable of producing

	for i in np.arange(0, len(ys)): # energies is on the y
		# ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])
		alpha_energy = ys[i]
		ionizations_per_alpha = np.divide(zs[i],num_alphas)
		# ax.plot(Propellant_densities, ionizations_per_alpha)
		spec_rad = specific_radioactivities[i] # Bcq/g, decays/gs
		prop_particles_in_chamber = np.multiply(Propellant_densities,Chamber_volume*6.02e23/prop_molar_mass)
		num_alphas_full_chamber = np.divide(prop_particles_in_chamber,ionizations_per_alpha)
		RI_molar_mass = RI_molar_masses[i]
		RI_mass_full_chamber = np.multiply(num_alphas_full_chamber,RI_molar_mass/(6.02e23)) # mass in g
		time_to_ionize = np.divide(num_alphas_full_chamber,np.multiply(RI_mass_full_chamber,spec_rad*num_full_ionizations)) # time in s

		Propellant_mass = np.multiply(Propellant_densities,Chamber_volume)
		Propellant_mass_per_RI_mass = np.divide(Propellant_mass,RI_mass_full_chamber)
		# ax.plot(Propellant_densities,np.multiply(RI_mass_full_chamber,1e6)) # mass in ug
		ax.plot(Propellant_densities, np.multiply(time_to_ionize,1e3)) # time in ms
		# ax.plot(Propellant_densities,Propellant_mass_per_RI_mass)
		# ax.plot(RI_mass_full_chamber,Propellant_mass)
		# ax.plot(Propellant_densities, num_alphas_full_chamber)


	xlabel = propellant + " " + r'$\rho$' + " (g/cm" r'$^3$' + ")"
	# xlabel = "Num of propellant particles in chamber"
	ylabel = "Alpha energy (MeV)"
	# zlabel = "Ionizations per alpha at first ionization energy"
	# zlabel = "RI mass to fully ionize the chamber once (" + r'$\mu$' + "g)"
	zlabel = "Time to fully ionize chamber (ms)"
	# zlabel = "Propellant mass ionized per gram of RI"
	# zlabel = "Number of alphas needed for full ionization"

	# ax.set_ylim([0,4])
	# axes.set_zlim([0,20])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(zlabel)
	ax.set_title(zlabel + " vs " + xlabel)
	ax.grid()
	ax.legend(("Po209, energy = " + str(ys[0]), 
			   "Po208, energy = " + str(ys[1]), 
			   "Pu238, energy = " + str(ys[2]), 
			   "Cm244, energy = " + str(ys[3]), 
			   "Am241, energy = " + str(ys[4])),
				loc=4)

	# plt.draw() 
	mass_plt_file_name = propellant + "_plots/RI_mass_vs_density.pdf"
	first_IE_plt_file_name = propellant + "/first_IE_vs_density.pdf"
	mean_IE_plt_file_name = propellant + "/mean_IE_vs_density.pdf"

	# fig.savefig(mass_plt_file_name, bbox_inches='tight') # save as pdf
	plt.show()





if __name__ == '__main__':
	energyRange = [4.9,5.1,5.6,5.8,5.5]
	# ioniz_energyRange = [4.9]# ,5.1,5.6,5.8,5.5]

	RI_densities = [9.196, 9.196, 9.196, 19.8, 13.51, 12.] # g/cm3    Po210, Po209, Po208, Pu238, Cm244, Am241
	# energyRange = [5.3,4.9]
	# densityRange =  list(np.arange(0.006,.05, .0001))
	# densityRange_2 = list(np.arange(0.003,.05, .0001)) #use for SR, mean
	densityRange_2 = list(np.arange(0.003,.015, .0001)) # use for first ioniz


	# energy = 5.3
	SRList = []
	totalIonizNumList = []



	propellant = raw_input("Which propellant do you want to use? ")

	for energy in energyRange:
		# ALL OF THE PREVIOUS DATA FROM THE FIRST BATCH OF FOLDERS
		# IS BAD!!! IT WAS RECORDING ELECTRONS, NOT ALPHAS
		# MAYBE IT COULD BE USEFUL THO......
		# SR_file = propellant + "2/"+str(energy)+"_SR.txt"
		totalIoniz_file = propellant + "2/"+str(energy)+"_first_ionization_num.txt"

		# using total ionization energies:
		# Source: # https://en.wikipedia.org/wiki/Molar_ionization_energies_of_the_elements ##
		if propellant == "cesium":
			# iEnergy_T = 3.7511e+19 # total ionization energy per mol of propellant MeV
			prop_molar_mass = 132.90545 # g/mol
		if propellant == "bismuth":
			# iEnergy_T = 1.43985e+20 # total ionization energy per mol of propellant MeV
			prop_molar_mass = 208.9804 # g/mol
		if propellant == "mercury":
			# iEnergy_T = 3.817994e+19 # total ionization energy per mol of propellant MeV
			prop_molar_mass = 200.59 # g/mol
		if propellant == "xenon":
			# iEnergy_T = 3.942262e+19 # total ionization energy per mol of propellant MeV
			prop_molar_mass = 131.293 # g/mol
		if propellant == "iodine":
			# iEnergy_T = 3.766314e+19 # total ionization energy per mol of propellant MeV
			prop_molar_mass = 126.90447 # g/mol

		# stopping_ranges = read_data(SR_file)
		# SRList.append(stopping_ranges)

		ioniz_nums = read_data(totalIoniz_file)
		alpha_ioniz_nums = ioniz_nums[::2].flatten()
		elec_ioniz_nums = ioniz_nums[1::2].flatten()
		# total_ioniz_nums = np.sum(alpha_ioniz_nums, elec_ioniz_nums)
		print ioniz_nums
		totalIonizNumList.append(alpha_ioniz_nums)
		# print stopping_ranges
	# plotter(densityRange, RI_densities, energyRange, SRList)
	ioniz_plotter(densityRange_2, energyRange, totalIonizNumList)
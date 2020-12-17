import matplotlib as mpl
from matplotlib import rc
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.linear_model import LinearRegression
plt.rcParams["font.family"] = "Times New Roman"
mpl.rcParams.update({'font.size': 24})
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

def plotter(densities, RI_densities, energies, SRList):

	# `self.stoppingRanges)
	# print(densities)

	xs = densities
	ys = energies
	zs = SRList

	# xs = [0.001, 0.011, 0.020999999999999998, 0.030999999999999996, 0.040999999999999995, 0.05099999999999999, 0.06099999999999999, 0.071, 0.08099999999999999, 0.09099999999999998]
	# ys = [0.1, 4.6 ,9.1]
	# zs = [[0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

	# print(len(xs),len(ys), len(zs[0]))

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
	ax1 = fig.add_subplot(111)

	# RI_masses = np.multiply(RI_densities, np.multiply(4*np.pi/3, np.power(np.subtract(2,zs[i]),3)))
	Chamber_volume = (np.pi*2**2)*4# cm3
	Propellant_density = xs

	decay_num = 1e5

	for i in np.arange(0, len(ys)): # energies is on the y
		# ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])
		alpha_energy = ys[i]
		ax1.plot(Propellant_density, zs[i])
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


	ax1.set_xlim([min(xs),max(xs)])
	# axes.set_zlim([0,20])
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(zlabel)
	ax1.set_title(zlabel + " vs " + xlabel)
	ax1.grid()


	# ax2.set_xlabel("Propellant_density (g/cm3)")	
	# ax2.set_ylabel("alpha_mass (g)"	)

	# ax2.grid()
	ax1.legend(("Po209, energy = " + str(ys[0]) + " MeV", 
			   "Po208, energy = " + str(ys[1]) + " MeV", 
			   "Am241, energy = " + str(ys[2]) + " MeV",
			   "Pu238, energy = " + str(ys[3]) + " MeV", 
			   "Cm244, energy = " + str(ys[4]) + " MeV"),
				loc=1)
	# plt.draw() 
	# fig.savefig(propellant + "_plots/" + propellant + "_SR_vs_density.pdf", bbox_inches='tight') # save as pdf

	plt.show()


# bracket dimensions based on Cm244 (max) and Po209 (min)
y1s = []
y2s = []


def ioniz_plotter(densities, energies, totalIonizNumList, prop_molar_mass, ax_plt):

	xs = densities
	ys = energies
	zs = totalIonizNumList


	# plt.xlabel(xlabel)
	# plt.ylabel(ylabel)
	# plt.title(ylabel + " vs " + xlabel)
	# plt.grid()
	# plt.scatter(densities, self.stoppingRanges)
	# plt.show()

	# first subplot: a 3D scatter plot of positions
	# ax = fig.add_subplot(111, projection='3d')


	fig = plt.figure()
	ax = fig.add_subplot(111)


	# RI_masses = np.multiply(RI_densities, np.multiply(4*np.pi/3, np.power(np.subtract(2,zs[i]),3)))
	pct_chamber_ionized = 1
	num_alphas = 50
	Chamber_volume = pct_chamber_ionized*(np.pi*2**2)*4# cm3
	Propellant_densities = xs
	# specific_radioactivities = [0.63e12,21.8e12,0.643e12,3.03e12,0.126e12] # Bcq/(g*s) Po209, Po208, Pu238, Cm244, Am241
	specific_radioactivities = [0.63e12,21.8e12,0.126e12,0.643e12,3.03e12] # Bcq/(g*s) Po209, Po208, Am241, Pu238, Cm244
	# time_elapsed = 1e-9
	# RI_molar_masses = [208.9824,207.9812,238.04956,244.06275,241.05683] # g/mol Po209, Po208, Pu238, Cm244, Am241
	RI_molar_masses = [208.9824, 207.9812, 241.05683, 238.04956, 244.06275] # g/mol Po209, Po208, Am241, Pu238, Cm244

	half_lives = np.multiply([109, 2.9, 432.2, 97.8, 18.1], 3.1536e7)	# s Po209, Po208,  Am241, Pu238, Cm244
	num_full_ionizations = 1e6 				# number of full-chamber ionizations the RI should be capable of producing
	time_to_ionize = 1 					# time to fully ionize chamber, s
	duration = np.arange(0,20*3.1536e7,86400) 	#
	RI_m0 = 1									# inital mass, g

	RI_mass_vs_density_slopes = []
	densityRange_test = np.arange(0,max(xs),0.0001)

	for i in np.arange(0, len(ys)): # energies is on the y
		# ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])
		alpha_energy = ys[i]
		spec_rad = specific_radioactivities[i] # Bcq/g, decays/gs
		half_life = half_lives[i]


		ionizations_per_alpha = np.divide(zs[i],num_alphas)
		print(propellant, alpha_energy, max(ionizations_per_alpha))

		# if alpha_energy == 5.8: # Cm244, max
		# 	y1s.append(max(ionizations_per_alpha)*1.e-6)
		# if alpha_energy == 4.9: # Po209, min
		# 	y2s.append(min(ionizations_per_alpha)*1.e-6)

		# ax_plt.plot(Propellant_densities, np.divide(ionizations_per_alpha,1e6))
		# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

		prop_particles_in_chamber = np.multiply(Propellant_densities,Chamber_volume*6.02e23/prop_molar_mass)
		num_alphas_full_chamber = np.divide(prop_particles_in_chamber,ionizations_per_alpha)
		RI_molar_mass = RI_molar_masses[i]
		RI_mass_full_chamber = np.multiply(num_alphas_full_chamber,RI_molar_mass/(6.02e23)) # mass in g
		# RI_mass_full_chamber_SPECACTIVITY = np.divide(num_alphas_full_chamber,spec_rad) # uses specific activity, I don't think I can use it like this lol
		# get slope of RI mass vs density line
		lm = LinearRegression()
		lm.fit(np.array(Propellant_densities[30:]).reshape(-1,1),np.array(RI_mass_full_chamber[30:])) # omits the litte jag at the beginning [30:120]
		RI_mass_vs_density_slopes.append(lm.coef_[0]) 
		# ####################################
		# Prop density vs t throughout mission duration
		Prop_mass_mission = np.multiply(np.exp(np.multiply(duration,np.float(-0.693/half_life))), np.float(0.693*RI_m0*time_to_ionize/(half_life*np.float(lm.coef_[0]))))

		Prop_mass_mission_test = np.multiply(np.exp(np.multiply(duration,np.float(-0.693/half_life))), 
											 np.float(0.693*RI_m0*prop_molar_mass*max(ionizations_per_alpha)*time_to_ionize/(RI_molar_mass*half_life)))
		# time_to_ionize = np.divide(num_alphas_full_chamber,np.multiply(RI_mass_full_chamber,spec_rad*num_full_ionizations)) # time in s

		# Prop_mass_mission = np.multiply(Prop_density_mission,Chamber_volume)
		RI_mass_full_chamber_test = np.multiply(densityRange_test,lm.coef_[0])
		# Prop_mass_full_chamber_test = np.multiply(densityRange_test,Chamber_volume)

		# ax_plt.plot(Prop_mass_full_chamber_test, np.multiply(RI_mass_full_chamber_test,1e6))
		# Propellant_mass_per_RI_mass = np.divide(Propellant_mass,RI_mass_full_chamber)
		ax_plt.plot(densityRange_test,np.multiply(RI_mass_full_chamber_test,1e6)) # mass in ug
		# ax.plot(Propellant_densities, np.multiply(time_to_ionize,1e3)) # time in ms
		# ax.plot(Propellant_densities,Propellant_mass_per_RI_mass)
		# ax.plot(RI_mass_full_chamber,Propellant_mass)
		ax.plot(np.divide(duration,3.1536e7), np.multiply(Prop_mass_mission_test,1e6/Chamber_volume)) # prop mass in ug

		if alpha_energy == 5.5: # Am241, max
			y1s.append(max(RI_mass_full_chamber_test)*1e6)
		if alpha_energy == 5.1: # Po208, min
			y2s.append(max(RI_mass_full_chamber_test)*1e6)


	# xlabel = "Num of propellant particles in chamber"
	# zlabel = "Millions of ionizations per alpha at first ionization energy"
	# zlabel = "Time to fully ionize chamber (ms)"
	# zlabel = "Propellant mass ionized per gram of RI"
	ylabel = "Alpha energy (MeV)"
	# zlabel =  propellant + " " + r'$\rho$' + " (g/cm" r'$^3$' + ")"
	# zlabel = propellant + " mass (" + r'$\mu$' + "g)"

	xlabel = "t (years)"
	zlabel = propellant.capitalize() + " m" r'$_{P}$'+ " (" + r'$\mu$' + "g)"
	ax.set_xlabel(xlabel)
	ax.set_ylabel(zlabel)
	ax.set_title(zlabel + " vs " + xlabel + " (first-ionized)")
	ax.legend(("Po209, energy = " + str(ys[0]) + " MeV", 
			   "Po208, energy = " + str(ys[1]) + " MeV", 
			   "Am241, energy = " + str(ys[2]) + " MeV",
			   "Pu238, energy = " + str(ys[3]) + " MeV", 
			   "Cm244, energy = " + str(ys[4]) + " MeV"),
				loc=1)

	# axes.set_zlim([0,20])
	xlabel = r'$\rho_{P}$' + " (g/cm" r'$^3$' + ")"
	# xlabel = propellant + " (propellant) mass (g)"
	# zlabel = "Number of alphas needed for full ionization"
	zlabel = r'$\Delta$' + "m" r'$_{RI}$' + " (" + r'$\mu$' + "g)"
	ax_plt.set_xlabel(xlabel)
	ax_plt.set_ylabel(zlabel)
	ax_plt.set_title(zlabel + " vs " + xlabel + " (first-ionized)")
	ax_plt.legend(("Po209, energy = " + str(ys[0]) + " MeV", 
			   "Po208, energy = " + str(ys[1]) + " MeV", 
			   "Am241, energy = " + str(ys[2]) + " MeV",
			   "Pu238, energy = " + str(ys[3]) + " MeV", 
			   "Cm244, energy = " + str(ys[4]) + " MeV"),
				loc=2)
	ax.grid()
	ax.text(12,.8*np.mean(ax.get_ylim()), r'$M_{0}$'+ ' = ' + str(RI_m0) + ' g' + "\n" + r'$\Delta$' + r'$t_{ionize}$' + ' = ' + str(time_to_ionize) + ' s',fontsize=20)

	# plt.draw() 
	RI_mass_plt_file_name = propellant + "_plots/" + propellant + "_RI_mass_vs_density.pdf"
	first_IE_plt_file_name = propellant + "/" + propellant + "_first_IE_vs_density.pdf"
	mean_IE_plt_file_name = propellant + "/" + propellant + "_mean_IE_vs_density.pdf"
	# prop_mass_plt_file_name = propellant + "_plots/" + propellant + "_prop_mass_vs_mission_time__first.pdf"
	prop_mass_plt_file_name = "/home/nasa01/Documents/howe_internship/axis/AXIS_Report_I_Plots/" + propellant + "_prop_mass_vs_mission_time__first.pdf"
	# fig.savefig(prop_mass_plt_file_name, bbox_inches='tight') # save as pdf
	print(RI_mass_vs_density_slopes)
	# plt.show()

def alpha_en_vs_IE():
	energyRange = [4.9,5.1,5.5,5.6,5.8]
	prop_names = ["Cesium", "Bismuth", "Xenon", "Mercury", "Iodine"]
	IE_list = [[1.26e6, 1.31e6, 1.41e6, 1.44e6, 1.49e6], 
			   [6.73e5, 7.00e5, 7.55e5, 7.69e5, 7.96e5],
			   [4.04e5, 4.20e5, 4.53e5, 4.62e5, 4.78e5],
			   [4.69e5, 4.88e5, 5.27e5, 5.36e5, 5.56e5],
			   [4.69e5, 4.88e5, 5.26e5, 5.36e5, 5.55e5]]
	
	fig2 = plt.figure()
	ax3 = fig2.add_subplot(111)



	for i in range(len(prop_names)):
		ax3.scatter(energyRange,np.multiply(IE_list[i],1e-6))


	xlabel = r'$E_{\alpha}$'
	ylabel = "I"
	ax3.set_xlabel(xlabel)
	ax3.set_ylabel(ylabel + " (millions of ionizations)")
	ax3.set_title(ylabel + " vs " + xlabel + " (first-ionized)")
	ax3.grid()
	ax3.legend(prop_names, loc=4)

def main(propellant, ax_plt, densityRange):
	energyRange = [4.9,5.1,5.5,5.6,5.8]
	# ioniz_energyRange = [4.9]# ,5.1,5.6,5.8,5.5]
	# energyRange = [5.3,4.9]



	# energy = 5.3
	SRList = []
	totalIonizNumList = []



	# propellant = raw_input("Which propellant do you want to use? ")

	for energy in energyRange:
		# ALL OF THE PREVIOUS DATA FROM THE FIRST BATCH OF FOLDERS
		# IS BAD!!! IT WAS RECORDING ELECTRONS, NOT ALPHAS
		# MAYBE IT COULD BE USEFUL THO......
		SR_file = propellant + "2/"+str(energy)+"_SR.txt"
		total_first_Ioniz_file = propellant + "2/"+str(energy)+"_first_ionization_num.txt"
		total_mean_Ioniz_file = propellant + "2/"+str(energy)+"_ionization_num.txt"

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

		ioniz_nums = read_data(total_first_Ioniz_file)
		# ioniz_nums = read_data(total_mean_Ioniz_file)
		alpha_ioniz_nums = ioniz_nums[::2].flatten() # all even indexed counts
		elec_ioniz_nums = ioniz_nums[1::2].flatten() # all odd indexed counts
		total_ioniz_nums = np.add(alpha_ioniz_nums, elec_ioniz_nums)
		# print ioniz_nums
		stopping_ranges = read_data(SR_file)
		SRList.append(stopping_ranges)
		totalIonizNumList.append(total_ioniz_nums)
		# print stopping_ranges

	# plotter(densityRange, RI_densities, energyRange, SRList)
	ioniz_plotter(densityRange, energyRange, totalIonizNumList, prop_molar_mass, ax_plt)

def CurlyBrace(color, ll_corner=(0, 0), width=1, height=1): # from https://stackoverflow.com/questions/52767053/making-annotations-on-axis-of-heatmaps 
	Path = mpath.Path
	verts = np.array([(0, 0), (.5, 0), (.5, .2), (.5, .3), (.5, .5), (1, .5), (.5, .5), (.5, .7), (.5, .8), (.5, 1), (0, 1)]) 
	verts[:, 0] *= width
	verts[:, 1] *= height
	verts[:, 0] += ll_corner[0]
	verts[:, 1] += ll_corner[1]

	cb_patch = mpatches.PathPatch(
		Path(verts,
			 [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3]),
		fc="none", ec=color, clip_on=False, transform=ax_plt.transData)
	return cb_patch




if __name__ == '__main__':

	# alpha_en_vs_IE()
	propellants_list = ["cesium", "bismuth", "xenon", "mercury", "iodine"]
	NUM_COLORS = 5 # num of radioisotopes
	cm = plt.get_cmap('flag')
	# cm = plt.get_cmap('Paired')
	fig = plt.figure()
	ax_plt = fig.add_subplot(111)
	ax_plt.grid()
	ax_plt.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

  	# densityRange =  list(np.arange(0.006,.05, .0001))
	# densityRange_2 = list(np.arange(0.003,.05, .0001)) #use for SR, mean
	densityRange_2 = list(np.arange(0.003,.015, .0001)) # use for first ioniz

	for propellant in propellants_list:
		main(propellant, ax_plt, densityRange_2)

	# widths = [3,1.5,3,1.5,1]
	# widths = [1,1.5,1,1.5,1]
	widths = [.5,.5,.5,.5,.5]
	# colors = ["red", "black", "red", "black", "black"]
	colors = ["black", "black", "black", "black", "red"]
	# print(y1s, "\n", y2s)
	duration = np.arange(0,20*3.1536e7,86400)
	x = np.divide(duration,3.1536e7)
	ax_plt.set_ylim([0,max(y1s)+0.5])
	for y, h, w, c, prop in zip(y1s, np.subtract(y2s,y1s), widths, colors, propellants_list): # <--------------------mess with this
		cb = CurlyBrace("black", [densityRange_2[-1], y+h*0.025], w/1000., h*.95)
		ax_plt.add_patch(cb)
		ax_plt.set_xlim([0,densityRange_2[-1]+0.0025+w/1000.])
		ax_plt.text(densityRange_2[-1]+0.0001+w/1000., y+h/2, prop, va='center', color="black")

	plt.show()
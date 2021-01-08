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

def plotter(densities, energies, SRList):

	xs = densities
	ys = energies
	zs = SRList

	xlabel = propellant.capitalize() + " " + r'$\rho$' + " (g/cm" r'$^3$' + ")"
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


	ax1.set_xlim([min(xs),max(xs)])
	ax1.set_ylim([0,40])
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(zlabel)
	ax1.set_title(zlabel + " vs " + xlabel)
	ax1.grid()

	ax1.legend(("Po209, energy = " + str(ys[0]) + " MeV", 
			   "Po208, energy = " + str(ys[1]) + " MeV", 
			   "Am241, energy = " + str(ys[2]) + " MeV",
			   "Pu238, energy = " + str(ys[3]) + " MeV", 
			   "Cm244, energy = " + str(ys[4]) + " MeV"),
				loc=1)
	# plt.draw() 
	# fig.savefig("/home/nasa01/Documents/howe_internship/axis/AXIS_Report_I_Plots/" + propellant + "_SR_vs_density.pdf", bbox_inches='tight') # save as pdf

	# plt.show()


# bracket dimensions based on Cm244 (max) and Po209 (min)
y1s = []
y2s = []

def ioniz_plotter(propellant, densities, energies, totalIonizNumList, prop_molar_mass, dRI_vs_dens, mdot_vs_time):
	###################################################################################################################################################	
	################################################################### GIVEN DATA ####################################################################
	###################################################################################################################################################
	# ASSUMES 100% IONIZATION IN CHAMBER
	pct_chamber_ionized = 1
	# NUM OF ALPHAS USED IN G4 SIM
	num_alphas = 50
	# CHAMBER VOLUME BASED ON PRELIM CAD MODEL, 1P UNIT
	Chamber_volume = pct_chamber_ionized*(np.pi*2**2)*4 # cm3
	# PROPELLANT DENSITIES, DETERMINED TO HAVE NO EFFECT ON IONIZATION EFFICIENCY
	Propellant_densities = densities # g/cm3
	# RI SPECIFIC ACTIVITIES, FROM FIRESTORM PHASE I REPORT
	specific_radioactivities = [0.63e12,21.8e12,0.126e12,0.643e12,3.03e12] # Bcq/(g*s) 	Po209, Po208, Am241, Pu238, Cm244
	# RI MOLAR MASSES, FROM PERIODIC TABLE
	RI_molar_masses = [208.9824, 207.9812, 241.05683, 238.04956, 244.06275] # g/mol 	Po209, Po208, Am241, Pu238, Cm244
	# RI HALF LIVES, FROM FIRESTORM PHASE I REPORT
	half_lives = np.multiply([109, 2.9, 432.2, 97.8, 18.1], 3.1536e7)	# s 	 		Po209, Po208, Am241, Pu238, Cm244
	# MISSION DURATION, SET TO 20 YEARS
	duration = np.arange(0,20*3.1536e7,86400) # s
	# INITIAL MASS OF RADIOISOTOPE, PROPORTIONAL TO MASS FLOW RATE
	RI_m0 = 1 # g
	# AVOGADRO'S NUMBER
	N_A = 6.02e23 # atoms/mol

	RI_mass_vs_density_slopes = []
	# MORE CONTINUOUS DENSITY RANGE
	densityRange_test = np.arange(0,max(Propellant_densities),0.0001)
	###################################################################################################################################################
	###################################################################################################################################################

	###################################################################################################################################################
	################################################# SET UP DEL_RI_MASS VS PROP MASS (DENSITY) PLOT ##################################################
	###################################################################################################################################################
	xlabel = r'$\rho_{P}$' + " (g/cm" r'$^3$' + ")"
	ylabel = r'$\Delta$' + "m" r'$_{RI}$' + " (" + r'$\mu$' + "g)"
	dRI_vs_dens.set_xlabel(xlabel)
	dRI_vs_dens.set_ylabel(ylabel)
	dRI_vs_dens.set_title(ylabel + " vs " + xlabel + " (first-ionized)")
	dRI_vs_dens.legend(("Po209, energy = " + str(energies[0]) + " MeV", 
			   	   "Po208, energy = " + str(energies[1]) + " MeV", 
			   	   "Am241, energy = " + str(energies[2]) + " MeV",
			   	   "Pu238, energy = " + str(energies[3]) + " MeV", 
			   	   "Cm244, energy = " + str(energies[4]) + " MeV"),
					loc=2)

	###################################################################################################################################################
	###################################################################################################################################################



	for RI in range(len(energies)): 

		# GIVEN DATA FOR A PARTICULAR RI
		alpha_energy = energies[RI]
		spec_rad = specific_radioactivities[RI] # Bcq/g, decays/gs
		half_life = half_lives[RI]
		RI_molar_mass = RI_molar_masses[RI]

		##############################################################################################################################################
		################################### CALCULATES AND PLOTS TIME-DEPENDENT PROPELLANT MASS FLOW RATE ############################################
		##############################################################################################################################################
		# PER ALPHA IONIZATION EFFICIENCY, FROM G4 SIMS
		ionizations_per_alpha = np.divide(totalIonizNumList[RI],num_alphas)

		# IONIZATIONS PER ALPHA EVENTUALLY MAX OUT, LEVEL OFF
		max_ionizations_per_alpha = max(ionizations_per_alpha)

		# TIME-DEPENDENT PROPELLANT MASS FLOW RATE USING SPECIFIC ACTIVITY
		Prop_mdot_mission_A = np.multiply(prop_molar_mass*max_ionizations_per_alpha*spec_rad*RI_m0/N_A,
										  np.exp(np.multiply(duration,
										  					 np.float(-0.693/half_life))))

		# PLOT TIME-DEPENDENT PROPELLANT MASS FLOW RATE USING SPECIFIC ACTIVITY
		mdot_vs_time.plot(np.divide(duration,3.1536e7), np.multiply(Prop_mdot_mission_A,1e3)) # prop mdot mg/s
		##############################################################################################################################################
		##############################################################################################################################################

		##############################################################################################################################################	
		############################################## CALCULATES AND PLOTS DEL RI MASS VS PROP DENSITY ##############################################
		##############################################################################################################################################	
		prop_particles_in_chamber = np.multiply(Propellant_densities,Chamber_volume*N_A/prop_molar_mass)
		num_alphas_full_chamber = np.divide(prop_particles_in_chamber,ionizations_per_alpha)
		RI_mass_full_chamber = np.multiply(num_alphas_full_chamber,RI_molar_mass/(N_A)) # mass in g
		lm = LinearRegression()
		lm.fit(np.array(Propellant_densities[30:]).reshape(-1,1),np.array(RI_mass_full_chamber[30:])) # omits the litte jag at the beginning [30:120]
		RI_mass_vs_density_slopes.append(lm.coef_[0]) 
		RI_mass_full_chamber_test = np.multiply(densityRange_test,lm.coef_[0])

		# GIVES UPPER AND LOWER BOUND FOR BRACKETS SET
		dRI_vs_dens.plot(densityRange_test,np.multiply(RI_mass_full_chamber_test,1e6)) # mass in ug
		if alpha_energy == 5.8: # Cm244, max
			y1s.append(max(RI_mass_full_chamber_test)*1e6)
		if alpha_energy == 5.1: # Po208, min
			y2s.append(max(RI_mass_full_chamber_test)*1e6)
		##############################################################################################################################################		
		##############################################################################################################################################		

	##################################################################################################################################################
	############################################################ SET UP MDOT VS TIME PLOT ############################################################ 
	##################################################################################################################################################
	xlabel = "t (years)"
	ylabel = propellant.capitalize() + " " r'$\dot{m}_{P}$'+ " (mg/s)"
	# ylabel = propellant.capitalize() + " m" r'$_{P}$'+ " (kg)"
	# ylabel = propellant.capitalize() + " m" r'$_{P}$'+ " (" + r'$\mu$' + "g)"
	mdot_vs_time.set_xlim([0,20])	
	mdot_vs_time.set_ylim([0,8])
	mdot_vs_time.set_xlabel(xlabel)
	mdot_vs_time.set_ylabel(ylabel)
	mdot_vs_time.set_title("       " + ylabel + " vs " + xlabel + " (first-ionized)", y=1.04)
	mdot_vs_time.grid()
	mdot_vs_time.text(12,.8*np.mean(mdot_vs_time.get_ylim()), r'$M_{0}$'+ ' = ' + str(RI_m0) + ' g',fontsize=20)
	mdot_vs_time.legend(("Po209, energy = " + str(energies[0]) + " MeV", 
			   	   "Po208, energy = " + str(energies[1]) + " MeV", 
			   	   "Am241, energy = " + str(energies[2]) + " MeV",
			   	   "Pu238, energy = " + str(energies[3]) + " MeV", 
			   	   "Cm244, energy = " + str(energies[4]) + " MeV"),
					loc=1)
	##################################################################################################################################################
	##################################################################################################################################################

	RI_mass_plt_file_name = propellant + "_plots/" + propellant + "_RI_mass_vs_density.pdf"

	# prop_mass_plt_file_name = propellant + "_plots/" + propellant + "_prop_mass_vs_mission_time__first.pdf"
	prop_mdot_plt_file_name = "/home/nasa01/Documents/howe_internship/axis/AXIS_Report_II_Plots/" + propellant + "_prop_mdot_vs_mission_time__first.pdf"
	# fig.savefig(prop_mdot_plt_file_name, bbox_inches='tight') # save as pdf
	print(RI_mass_vs_density_slopes)
	# plt.show()



def run(propellant, dRI_vs_dens, densityRange, mdot_vs_time, prop_molar_mass):
	energyRange = [4.9,5.1,5.5,5.6,5.8] # MeV; Po209, Po208, Am241, Pu238, Cm244

	SRList = []
	totalIonizNumList = []


	for energy in energyRange:

		SR_file = propellant + "2/"+str(energy)+"_SR.txt"
		total_first_Ioniz_file = propellant + "2/"+str(energy)+"_first_ionization_num.txt"
		total_mean_Ioniz_file = propellant + "2/"+str(energy)+"_ionization_num.txt"

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

	# plotter(densityRange, energyRange, SRList)
	ioniz_plotter(propellant, densityRange, energyRange, totalIonizNumList, prop_molar_mass, dRI_vs_dens, mdot_vs_time)

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
		fc="none", ec=color, clip_on=False)
	return cb_patch




def main():

	# alpha_en_vs_IE()
	propellants_dict = {"cesium":132.90545, 
						"bismuth":208.9804, 
						"xenon":131.293, 
						"mercury":200.59, 
						"iodine":126.90447}

	NUM_COLORS = 5 # num of radioisotopes
	cm = plt.get_cmap('flag')
	# cm = plt.get_cmap('Paired')
	fig, dRI_vs_dens = plt.subplots()
	# dRI_vs_dens = fig.add_subplot(111)
	dRI_vs_dens.grid()
	dRI_vs_dens.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])



	# densityRange_2 = list(np.arange(0.003,.05, .0001)) #use for SR, mean
	densityRange_2 = list(np.arange(0.003,.015, .0001)) # use for first ioniz

	for propellant, prop_molar_mass in propellants_dict.items():
		fig2 = plt.figure()
		mdot_vs_time = fig2.add_subplot(111)
		run(propellant, dRI_vs_dens, densityRange_2,mdot_vs_time,prop_molar_mass)
		# plt.show()

	# widths = [3,1.5,3,1.5,1]
	# widths = [1,1.5,1,1.5,1]
	widths = [.5,.5,.5,.5,.5]
	# colors = ["red", "black", "red", "black", "black"]
	colors = ["black", "black", "black", "black", "red"]
	# print(y1s, "\n", y2s)
	duration = np.arange(0,20*3.1536e7,86400)
	x = np.divide(duration,3.1536e7)
	dRI_vs_dens.set_ylim([0,max(y1s)+0.5])
	for y, h, w, c, prop in zip(y1s, np.subtract(y2s,y1s), widths, colors, list(propellants_dict.keys())): # <--------------------mess with this
		cb = CurlyBrace("black", [densityRange_2[-1], y+h*0.025], w/1000., h*.95)
		dRI_vs_dens.add_patch(cb)
		dRI_vs_dens.set_xlim([0,densityRange_2[-1]+0.0025+w/1000.])
		dRI_vs_dens.text(densityRange_2[-1]+0.0001+w/1000., y+h/2, prop, va='center', color="black")

	plt.show()

if __name__ == '__main__':
	main()
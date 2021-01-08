#----------GEANT4 IMPORTS----------#
from Geant4 import * 
from Geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex
import g4py.ParticleGun, g4py.MedicalBeam

# import g4py.GeneralParticleSource



# particleTable = G4ParticleTable.GetParticleTable()

# print particleTable.FindParticle("anttriton")

# --------PYTHON IMPORTS ----------#
from collections import Counter
import itertools
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import pandas as pd
import random
from scipy.signal import find_peaks
import scipy.stats as ss
import seaborn as sns  # for nicer graphics
import time
plt.rcParams["font.family"] = "Times New Roman"
# -------- FILE IMPORTS ------- #
# from arrow_generator import Arrow3D
global stoppingRangesList
stoppingRangesList = []
global elecDataList
elecDataList = []
global elecDataList_final
elecDataList_final = []
global elecPositionsList
elecPositionsList = []
global stoppingPowerList
stoppingPowerList = []
global ionizationDataList
ionizationDataList = []
global alphaIonizNum
alphaIonizNum = []
global elecIonizNum
elecIonizNum = []

#----------code starts here!----------#
class SetGlobalData(object):

	def __init__(self):
		pass

	def set_ionization_energies(self, propellant):
		global ioniz_energies_list
		if propellant == "cesium":
			ioniz_energies_list = [3.894e-6,23.16e-6,35.24e-6] # MeV
		if propellant == "bismuth":
			ioniz_energies_list = [7.286e-6,16.69e-6,25.56e-6] # MeV
		if propellant == "mercury":
			ioniz_energies_list = [10.44e-6,18.76e-6,34.2e-6,45.29e-6,55.97e-6,88.3e-6] # MeV
		if propellant == "xenon":
			ioniz_energies_list = [12.13e-6,21.21e-6,32.12e-6] # MeV
		if propellant == "iodine":
			ioniz_energies_list = [10.45e-6,19.13e-6,32.96e-6] # MeV
		if propellant =="lithium":
			ioniz_energies_list = [5.3915e-6,75.64e-6,122.45e-6] # MeV

class GetGlobalData(object):

	def __init__(self):
		pass

	def get_stopping_range(self):
		return stoppingRange_final

	def get_electron_num(self):
		return electron_number_final # returns num of electron within the stopping range

	def get_total_ioniz_num(self):
		# global alpha_ioniz_num
		# global elec_ioniz_num
		alpha_ioniz_num = len(alphaIonizNum) 
		elec_ioniz_num = len(elecIonizNum)
		print alpha_ioniz_num, elec_ioniz_num
		# time.sleep(2)
		alphaIonizNum[:] = []
		elecIonizNum[:] = []
		

		return alpha_ioniz_num, elec_ioniz_num

class DataAnalysis(object):

	'''
	Analyzes the stopping ranges of alphas of different 
		energies in different densities of Xe gas
	Produces a 3D plot OR 2D plots (range vs density) for
		chosen radioisotopes
	'''

	def __init__(self):

		stoppingRangesList[:] = []

	def data_collection(self, rng):

		stoppingRangesList.append(rng)
		# energies.append(energy)

	def elec_data_collection(self, elec_data):

		elecDataList.append(elec_data)

	def compute_stopping_range(self):

		# print stoppingRangesList

		# plt.hist(stoppingRangesList,100)
		# # ax2.scatter(energies, stoppingRangesList)
		# plt.show()
		amounts, vals = np.histogram(stoppingRangesList,100)
		index = list(amounts).index(max(amounts))
		peak_val = vals[index]
		print("\n")
		print("\n")
		print("\n")
		print(len(stoppingRangesList))
		print "Peak stopping range: ", peak_val, "mm"
		print("\n")
		print("\n")
		print("\n")
		global stoppingRange_final
		stoppingRange_final = peak_val

		# print(peak_val)

	def elec_gen_analysis(self):

		c = Counter()
		energies = (np.array(elecDataList).T[:2]).T
		# print len(energies) 
		# time.sleep(2)
		# dup_trackIDs = [item for item, count in Counter(elecDataList.T[0]).items() if count > 1]
		for trackID, eDep in energies:
			c.update({trackID:eDep})
		energies = list(c.items())
		genDeps = []
		# print("\n",max(np.array(elecDataList).T[0]))
		# print(set(np.array(elecDataList).T[0]))
		# print len(set(np.array(elecDataList).T[0]))

		trackID_filtered = dict(enumerate(pd.Index(np.array(elecDataList).T[0]).duplicated(keep='first')))
		trackID_filtered = {key:val for key, val in trackID_filtered.items() if val == True}
		trackID_filtered_indices_bad = list(trackID_filtered.keys())
		genDeps = np.delete(np.array(elecDataList).T[1],trackID_filtered_indices_bad)

		# for trackID in set(np.array(elecDataList).T[0]):
		# 	print "HELLO"
		# 	genDeps.append(np.array(elecDataList).T[1][list(np.array(elecDataList).T[0]).index(trackID)])

		print len(genDeps), len(energies)

		# elecDataList_final.append(energies) 
		# elecDataList_final.append(genDeps)

		# print(np.array(genDeps).shape())

		plt.hist(np.array(genDeps),100)
		plt.show()

	def electron_displayer(self):

		parentIds, trackIds, kineticEnergys, totalEnergyDeposits, positions = np.array(elecPositionsList).T
		
		# for p in positions

		positions_all = list(itertools.chain(*positions))
		print "num of e- trackIds = ", len(set(trackIds)), "num of unique tracks = ", max(trackIds) - min(trackIds)

		number_of_electrons = len(set(trackIds))

		trackID_filtered = dict(enumerate(pd.Index(trackIds).duplicated(keep='first')))
		trackID_filtered = {key:val for key, val in trackID_filtered.items() if val == True}
		trackID_filtered_indices_bad = list(trackID_filtered.keys())

		positions_filtered = np.delete(positions,trackID_filtered_indices_bad)
		print "num of filtered e- = ", len(positions_filtered)
		
		positions_filtered_all = list(itertools.chain(*positions_filtered))
		xs, ys, zs = positions_filtered_all[0::3], positions_filtered_all[1::3], positions_filtered_all[2::3]

		positions_filtered_combined = np.array([xs,ys,zs]).T
		for i in range(0,len(positions_filtered_combined)):
			pos = positions_filtered_combined[i]
			if all(i <= stoppingRange_final for i in pos) == False:
				np.delete(positions_filtered_combined, i)
		print "num of filtered e- within stopping range = ", len(positions_filtered_combined)
		global electron_number_final
		electron_number_final = len(positions_filtered_combined)

		# fig = plt.figure()
		# # first subplot: a 3D scatter plot of positions
		# ax = fig.add_subplot(111, projection='3d')
		# axmin = -10
		# axmax = 10
		# axes = plt.gca()
		# axes.set_xlim([axmin,axmax])
		# axes.set_ylim([axmin,axmax])
		# axes.set_zlim([axmin,axmax])
		# ax.set_xlabel('X (mm)')
		# ax.set_ylabel('Y (mm)')
		# ax.set_zlabel('Z (mm)')
		# ax.scatter(xs, ys, zs)
		# plt.title("Electron positions") 
		# plt.show()
		# pass

	def calc_ionization(self, procName, particleName, deltaEnergy, ionizingEnergy, energyLevel):

		'''
		procname = process (ionIoni, eIoni)
		particleName = alpha, e-
		deltaEnergy = change in energy due to ionization
		ionizingEnergy = change in energy due to ionization
		energyLevel = which ionization level is used (1st, 2nd, 3rd, "mean")
		'''
		if energyLevel != "mean":
			energyLevel -= 1 # python is indexed 0
			unit_ionization_energy = ioniz_energies_list[energyLevel]
		else:
			unit_ionization_energy = np.mean(ioniz_energies_list)

		ionization_num = int(ionizingEnergy/unit_ionization_energy)
		if procName == "ionIoni":
			for i in range(0,ionization_num):
				alphaIonizNum.append("alpha")
		if procName == "eIoni":
			for i in range(0,ionization_num):
				elecIonizNum.append("e-")

	def prop_stopping_power(self):

		stoppingPowerData = np.array(stoppingPowerList).T
		depths = stoppingPowerData[0]
		kineticEnergys = stoppingPowerData[1]
		totalEnergyDeposits = stoppingPowerData[2]
		totalEnergyDeposit_fracs = stoppingPowerData[3]
		totalEnergyDeposits_per_depth = np.divide(totalEnergyDeposits,depths)
		totalEnergyDeposit_fracs_per_depth = np.divide(totalEnergyDeposit_fracs, depths)

		zipped_x = depths
		zipped_y = totalEnergyDeposits_per_depth

		zipped_data = zip(zipped_x, zipped_y)
		sorted_zipped_data = sorted(zipped_data)
		sorted_y = [element for _, element in sorted_zipped_data]
		sorted_x = sorted(zipped_x)
		# print totalEnergyDeposit_fracs, "\n", sorted_depths

		plt.plot(sorted_x, sorted_y)
		plt.grid()
		plt.xlabel("Depths (mm)") 
		plt.ylabel("Energy deposited per depth (MeV/mm)")
		plt.title("Stopping power of propellant")
		plt.show()

DA = DataAnalysis()

'''																									     											   '''
   #                         								   ____________   _  ____________ 								                           #
   ########################################################## / ___/ __/ _ | / |/ /_  __/ / / ##########################################################
   ##########################################################/ (_ / _// __ |/    / / / /_  _/ ##########################################################
   ##########################################################\___/___/_/ |_/_/|_/ /_/   /_/   ##########################################################
   #																																				   #
'''																									    											   '''

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self, particle, energy, energyUnit, emitterCount, momentumArray):
		G4VUserPrimaryGeneratorAction.__init__(self)
		# self.ionGun = G4IonGun(1)
		self.pg = G4ParticleGun(1)
		self.pg.SetParticleByName(self.particle)
		g4py.ParticleGun.Construct()
		self.energy = energy
		self.energyUnit = energyUnit
		self.emitterCount = emitterCount
		self.momentumArray = momentumArray
		self.particle = particle


	def GeneratePrimaries(self, event):
		# Particle param
		#################################################
		# locationArray = [0, 0, 0]

		energyUnit = self.energyUnit 
		dimensionUnit = cm
		energy = self.energy
		print("\n\n\n\n\n\n")
		time.sleep(1)
		# emitter locations on outer/inner surfaces of radioisotope sphere
		# outer, 
		r = 10
		phi = np.linspace(0, np.pi, 20)
		theta = np.linspace(0, 2 * np.pi, 40)
		xs = []
		ys = []
		zs = []
		for p in phi:
			for t in theta:
				x = r*np.sin(t)*np.cos(p)
				y = r*np.sin(t)*np.sin(p)
				z = r*np.cos(t)

				xs.append(x)
				ys.append(y)
				zs.append(z)

		## Random generation from a point
		locationArray = [0,0,0]


		for i in range(0, self.emitterCount): # creates random momentum vectors originating from [0, 0, 0]
			mx = random.uniform(-1,1)
			my = random.uniform(-1,1)
			mz = random.uniform(-1,1)
			momentumArray = [mx, my, mz]
			gApplyUICommand("/gun/particle ion")
			gApplyUICommand("/gun/ion 92 233 0")
			gApplyUICommand("/gun/energy 90 MeV")
			gApplyUICommand("/gun/position 0. 0. 0. m")
			gApplyUICommand("/gun/direction " + str(mx) + " " 
											  + str(my) + " " 
											  + str(mz) + " ")
			# gRunManager.BeamOn(1)
			time.sleep(0.01)
			self.pg.GeneratePrimartVertex(event)

#-------------------------------------------------------------------

class MyRunAction(G4UserRunAction):
	"My Run Action"

	def EndOfRunAction(self, run):
		print "*** End of Run"
		print "- Run sammary : (id= %d, #events= %d)" \
		% (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
		stoppingRangesList[:] = []
		stoppingPowerList[:] = []
		elecPositionsList[:] = []
		pass

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		# DA.compute_stopping_range()
		# DA.elec_gen_analysis()
		# DA.prop_stopping_power()
		# if len(elecPositionsList) != 0:
		# DA.electron_displayer()
		# else:
		# 	global electron_number_final
		# 	electron_number_final = 0
		# total_alpha_ionizations = 0
		# total_elec_ionizations = 0
		# print "alphas = ", alpha_ioniz_num, "e- =" , elec_ioniz_num
		print "*** End of Event"
		pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
	"My Stepping Action"

	def UserSteppingAction(self, step):
		# stoppingRangesList[:] = []
		# energies[:] = []
		preStepPoint = step.GetPreStepPoint()
		postStepPoint = step.GetPostStepPoint()
		procName = postStepPoint.GetProcessDefinedStep().GetProcessName()
		track = step.GetTrack()
		parentId = track.GetParentID()
		trackId = track.GetTrackID()
		particleName = track.GetDefinition().GetParticleName() 
		touchable = track.GetTouchable()
		kineticEnergy = track.GetKineticEnergy()


		totalEnergyDeposit = step.GetTotalEnergyDeposit()
		nonIonizingEnergyDeposit = step.GetNonIonizingEnergyDeposit()
		deltaEnergy = step.GetDeltaEnergy() 
		ionizingEnergy = totalEnergyDeposit - nonIonizingEnergyDeposit
		# print parentId, trackId, "dE = ", deltaEnergy, "  ionizing energy = ", ionizingEnergy, " ", particleName


		# --- calculate number of ionizations --- #
		energyLevel = "mean"
		# energyLevel = 1
		DA.calc_ionization(procName, particleName, deltaEnergy, ionizingEnergy, energyLevel)


		# stoppedThreshold = .1 # 10 eV, less than the 1st ionization energy of Xe
		# print totalEnergyDeposit

		# if totalEnergyDeposit > 0.999*5.3:
			# if kineticEnergy < stoppedThreshold:
			# print particleName, " ", totalEnergyDeposit, " ", 25+postStepPoint.GetPosition().x
		# if kineticEnergy == 0.0 :
		# 	print parentId, trackId, stoppingRange, kineticEnergy, particleName

		# --- alphas, stopping range --- #
		# if parentId == 0 and kineticEnergy == 0.0:
		# 	position = [float(postStepPoint.GetPosition().x), float(postStepPoint.GetPosition().y), float(postStepPoint.GetPosition().z)]
		# 	stoppingRange  = np.sqrt(np.sum(np.power(position,2)))
		# 	# print parentId, trackId, stoppingRange, kineticEnergy, particleName
		# 	# DA.data_collection(stoppingRange)

		# --- alphas, ionization (process) data -- #
		# if parentId == 0:
		# 	print particleName, deltaEnergy, totalEnergyDeposit, procName
		# if procName == "ionIoni" or procName == "eIoni":
		# 	print particleName, parentId, trackId, deltaEnergy, totalEnergyDeposit, procName
		# if procName == "ionIoni":
		# 	alphaIonizNum.append([particleName, parentId, trackId, deltaEnergy, totalEnergyDeposit, procName])
		# 	print particleName, parentId, trackId, deltaEnergy, totalEnergyDeposit, procName
		# if procName == "eIoni":
		# 	elecIonizNum.append([particleName, parentId, trackId, deltaEnergy, totalEnergyDeposit, procName])
		# 	print particleName, parentId, trackId, deltaEnergy, totalEnergyDeposit, procName


		# # --- alphas, stopping power --- #
		# if particleName == "alpha":
		# 	if kineticEnergy == 0:
		# 		totalEnergyDeposit_frac = 0
		# 	else:
		# 		totalEnergyDeposit_frac = float(totalEnergyDeposit/kineticEnergy)
		# 	print parentId, trackId, 1250 + preStepPoint.GetPosition().x, kineticEnergy, totalEnergyDeposit, totalEnergyDeposit_frac
		# 	stoppingPowerList.append([1250 + preStepPoint.GetPosition().x, kineticEnergy, totalEnergyDeposit, totalEnergyDeposit_frac])

		# # --- electrons, position plotting --- #
		# if particleName == "e-":
		# 	position = [float(postStepPoint.GetPosition().x), float(postStepPoint.GetPosition().y), float(postStepPoint.GetPosition().z)]
		# 	print parentId, trackId, kineticEnergy, totalEnergyDeposit, position, procName
		# 	elecPositionsList.append([parentId, trackId, kineticEnergy, totalEnergyDeposit, position])

			# DA.elec_data_collection([trackId, totalEnergyDeposit, 1250 + preStepPoint.GetPosition().x])

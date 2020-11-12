#----------GEANT4 IMPORTS----------#
from Geant4 import * 
from Geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex



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


#----------code starts here!----------#
class GetGlobalData(object):

	def __init__(self):
		
		pass

	def get_stopping_range(self):
		return stoppingRange_final

	def get_electron_num(self):
		return electron_number_final # returns num of electron within the stopping range

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
		time.sleep(2)
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

		parentIds, trackIds, KEs, edeps, positions = np.array(elecPositionsList).T
		
		# for p in positions

		positions_all = list(itertools.chain(*positions))
		print "num of e- trackIds = ", len(set(trackIds)), "num of unique tracks = ", max(trackIds) - min(trackIds)

		number_of_electrons = len(set(trackIds))





		# c = Counter()
		# positions_unfiltered = np.concatenate(trackIds, positions)
		# for trackID, position in positions_unfiltered:
		# 	c.update({trackID:position})
		# positions = list(c.items())

		trackID_filtered = dict(enumerate(pd.Index(trackIds).duplicated(keep='first')))
		trackID_filtered = {key:val for key, val in trackID_filtered.items() if val == True}
		trackID_filtered_indices_bad = list(trackID_filtered.keys())

		positions_filtered = np.delete(positions,trackID_filtered_indices_bad)
		print "num of filtered e- = ", len(positions_filtered)
		
		positions_filtered_all = list(itertools.chain(*positions_filtered))
		xs, ys, zs = positions_all[0::3], positions_all[1::3], positions_all[2::3]

		fig = plt.figure()
		# first subplot: a 3D scatter plot of positions
		ax = fig.add_subplot(111, projection='3d')
		axmin = -10
		axmax = 10
		axes = plt.gca()
		axes.set_xlim([axmin,axmax])
		axes.set_ylim([axmin,axmax])
		axes.set_zlim([axmin,axmax])
		ax.set_xlabel('X (mm)')
		ax.set_ylabel('Y (mm)')
		ax.set_zlabel('Z (mm)')
		ax.scatter(xs, ys, zs)
		plt.title("Electron positions") 
		plt.show()
		pass

	def prop_stopping_power(self):

		stoppingPowerData = np.array(stoppingPowerList).T
		depths = stoppingPowerData[0]
		KEs = stoppingPowerData[1]
		edeps = stoppingPowerData[2]
		edep_fracs = stoppingPowerData[3]
		edeps_per_depth = np.divide(edeps,depths)
		edep_fracs_per_depth = np.divide(edep_fracs, depths)

		zipped_x = depths
		zipped_y = edeps_per_depth

		zipped_data = zip(zipped_x, zipped_y)
		sorted_zipped_data = sorted(zipped_data)
		sorted_y = [element for _, element in sorted_zipped_data]
		sorted_x = sorted(zipped_x)
		# print sorted_edep_fracs, "\n", sorted_depths

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
		self.particleGun_outer = G4ParticleGun(1)
		# self.particleGun_inner = G4ParticleGun(1)
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
		self.particleGun_outer.SetParticleByName(self.particle) # define particle
		self.particleGun_outer.SetParticleEnergy(energy*energyUnit) # define particle energy 
		# self.particleGun_inner.SetParticleByName(self.particle) # define particle
		# self.particleGun_inner.SetParticleEnergy(energy*energyUnit) # define particle energy 

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
			self.particleGun_outer.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*dimensionUnit) # define first particle generator location
			self.particleGun_outer.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0], momentumArray[1], momentumArray[2])*dimensionUnit) # define first particle generator momentum
			self.particleGun_outer.GeneratePrimaryVertex(event)

		# generate emitter locations
		# emitterNumber = self.particleCount
		# emitterNumber = int(np.sqrt(self.emitterCount))
		# emitterLocations = np.indices((emitterNumber, emitterNumber)).T.flatten().reshape(emitterNumber**2,2)

		# print(emitterLocations)
		# time.sleep(2)

		# for loc in emitterLocations:
		# 	locationArray = [-1250,loc[0]-int(emitterNumber/2), loc[1]-int(emitterNumber/2)]
		# 	locationArray = [-1250,0,0]
		# 	momentumArray = self.momentumArray
		# 	# locationArray = [loc[0], loc[1], loc[2]] # spherical
		# 	# momentumArray_outer = np.subtract(locationArray,0)
		# 	self.particleGun_outer.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*mm) # define first particle generator location
		# 	self.particleGun_outer.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0], momentumArray[1], momentumArray[2])*dimensionUnit) # define first particle generator momentum

		# 	self.particleGun_outer.GeneratePrimaryVertex(event)
		


		# emitterLocations = np.array([xs,ys,zs]).T
		# for loc in emitterLocations:
		# 	# locationArray = [-25,loc[0]-int(emitterNumber/2), loc[1]-int(emitterNumber/2)]
		# 	# momentumArray = self.momentumArray
		# 	locationArray = [loc[0], loc[1], loc[2]] # spherical
		# 	momentumArray_outer = np.subtract(locationArray,0)
		# 	self.particleGun_outer.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*mm) # define first particle generator location
		# 	self.particleGun_outer.SetParticleMomentumDirection(G4ThreeVector(momentumArray_outer[0], momentumArray_outer[1], momentumArray_outer[2])*dimensionUnit) # define first particle generator momentum

		# 	self.particleGun_outer.GeneratePrimaryVertex(event)

		# for loc in emitterLocations:

		# 	locationArray = [loc[0], loc[1], loc[2]] # spherical
		# 	momentumArray_inner = np.subtract(0,locationArray)
		# 	self.particleGun_inner.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*mm) # define first particle generator location
		# 	self.particleGun_inner.SetParticleMomentumDirection(G4ThreeVector(momentumArray_inner[0], momentumArray_inner[1], momentumArray_inner[2])*dimensionUnit) # define first particle generator momentum

		# 	self.particleGun_inner.GeneratePrimaryVertex(event)

		#################################################

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
		DA.compute_stopping_range()
		# DA.elec_gen_analysis()
		# DA.prop_stopping_power()
		DA.electron_displayer()
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
		track = step.GetTrack()
		parentId = track.GetParentID()
		trackId = track.GetTrackID()
		particleName = track.GetDefinition().GetParticleName() 
		touchable = track.GetTouchable()
		KE = track.GetKineticEnergy()

		edep = step.GetTotalEnergyDeposit()



		# stoppedThreshold = .1 # 10 eV, less than the 1st ionization energy of Xe
		# print edep

		# if edep > 0.999*5.3:
			# if KE < stoppedThreshold:
			# print particleName, " ", edep, " ", 25+postStepPoint.GetPosition().x
		# if KE == 0.0 :
		# 	print parentId, trackId, stoppingRange, KE, particleName

		# --- alphas, stopping range --- #
		if parentId == 0 and KE == 0.0:
			position = [float(postStepPoint.GetPosition().x), float(postStepPoint.GetPosition().y), float(postStepPoint.GetPosition().z)]
			stoppingRange  = np.sqrt(np.sum(np.power(position,2)))
			print parentId, trackId, stoppingRange, KE, particleName
			DA.data_collection(stoppingRange)

		# # --- alphas, stopping power --- #
		# if particleName == "alpha":
		# 	if KE == 0:
		# 		edep_frac = 0
		# 	else:
		# 		edep_frac = float(edep/KE)
		# 	print parentId, trackId, 1250 + preStepPoint.GetPosition().x, KE, edep, edep_frac
		# 	stoppingPowerList.append([1250 + preStepPoint.GetPosition().x, KE, edep, edep_frac])

		# # --- electrons, position plotting --- #
		if particleName == "e-":
			position = [float(postStepPoint.GetPosition().x), float(postStepPoint.GetPosition().y), float(postStepPoint.GetPosition().z)]
			print parentId, trackId, KE, edep, position
			elecPositionsList.append([parentId, trackId, KE, edep, position])

			# DA.elec_data_collection([trackId, edep, 1250 + preStepPoint.GetPosition().x])

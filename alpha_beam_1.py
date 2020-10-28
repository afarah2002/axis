#----------GEANT4 IMPORTS----------#
from Geant4 import * 
from Geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex



# particleTable = G4ParticleTable.GetParticleTable()

# print particleTable.FindParticle("anttriton")

# --------PYTHON IMPORTS ----------#
import itertools
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import random
from scipy.signal import find_peaks
import scipy.stats as ss
import seaborn as sns  # for nicer graphics
import time
plt.rcParams["font.family"] = "Times New Roman"
# -------- FILE IMPORTS ------- #
# from arrow_generator import Arrow3D

#----------code starts here!----------#
class GetGlobalData(object):

	def __init__(self):
		pass

	def get_stopping_range(self):
		return stoppingRange_final

class DataAnalysis(object):

	'''
	Analyzes the stopping ranges of alphas of different 
		energies in different densities of Xe gas
	Produces a 3D plot OR 2D plots (range vs density) for
		chosen radioisotopes
	'''

	def __init__(self):

		self.stoppingRangesList = []

	def data_collection(self, range):

		self.stoppingRangesList.append(range)

	def compute_stopping_range(self):

		amounts, vals = np.histogram(self.stoppingRangesList,100)
		index = list(amounts).index(max(amounts))
		peak_val = vals[index]
		# print(len(self.stoppingRangesList))
		global stoppingRange_final
		stoppingRange_final = peak_val

		# print(peak_val)


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
		self.particleGun_inner = G4ParticleGun(1)
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
		self.particleGun_inner.SetParticleByName(self.particle) # define particle
		self.particleGun_inner.SetParticleEnergy(energy*energyUnit) # define particle energy 

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


		# generate emitter locations
		# emitterNumber = self.particleCount
		emitterNumber = int(np.sqrt(self.emitterCount))
		emitterLocations = np.indices((emitterNumber, emitterNumber)).T.flatten().reshape(emitterNumber**2,2)

		# print(locations)

		for loc in emitterLocations:
			locationArray = [-25,loc[0]-int(emitterNumber/2), loc[1]-int(emitterNumber/2)]
			momentumArray = self.momentumArray
			# locationArray = [loc[0], loc[1], loc[2]] # spherical
			# momentumArray_outer = np.subtract(locationArray,0)
			self.particleGun_outer.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*mm) # define first particle generator location
			self.particleGun_outer.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0], momentumArray[1], momentumArray[2])*dimensionUnit) # define first particle generator momentum

			self.particleGun_outer.GeneratePrimaryVertex(event)
		
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
		# print "*** End of Run"
		# print "- Run sammary : (id= %d, #events= %d)" \
		# % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
		DA.compute_stopping_range()
		pass

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
	"My Stepping Action"

	def UserSteppingAction(self, step):
		preStepPoint = step.GetPreStepPoint()
		postStepPoint = step.GetPostStepPoint()
		track = step.GetTrack()
		parentId = track.GetParentID()
		particleName = track.GetDefinition().GetParticleName() 
		touchable = track.GetTouchable()
		KE = track.GetKineticEnergy()

		stoppedThreshold = 0.1 # MeV

		if parentId == 0:
			# print particleName, " ", KE, " ", 25+postStepPoint.GetPosition().x
			if KE < stoppedThreshold:
				stoppingRange = 25 + postStepPoint.GetPosition().x
				DA.data_collection(stoppingRange)


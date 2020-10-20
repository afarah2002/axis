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

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self, particle, energy, energyUnit, particleCount, momentumArray):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		self.energy = energy
		self.energyUnit = energyUnit
		self.particleCount = particleCount
		self.momentumArray = momentumArray
		self.particle = particle


	def GeneratePrimaries(self, event):
		# Particle param
		#################################################
		# locationArray = [0, 0, 0]

		energyUnit = self.energyUnit 
		dimensionUnit = cm
		energy = self.energy
		self.particleGun.SetParticleByName(self.particle) # define particle
		self.particleGun.SetParticleEnergy(energy*energyUnit) # define particle energy 

		# generate emitter locations
		emitterNumber = 10
		emitterLocations = np.indices((emitterNumber, emitterNumber)).T.flatten().reshape(emitterNumber**2,2)

		for loc in emitterLocations:
			locationArray = [-25,loc[0]-int(emitterNumber/2), loc[1]-int(emitterNumber/2)]
			self.particleGun.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*mm) # define first particle generator location
			self.particleGun.SetParticleMomentumDirection(G4ThreeVector(self.momentumArray[0], self.momentumArray[1], self.momentumArray[2])*dimensionUnit) # define first particle generator momentum
			self.particleGun.GeneratePrimaryVertex(event)
		#################################################

#-------------------------------------------------------------------

class MyRunAction(G4UserRunAction):
	"My Run Action"

	def EndOfRunAction(self, run):
		# print "*** End of Run"
		# print "- Run sammary : (id= %d, #events= %d)" \
		# % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
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

		if parentId == 0:
			print particleName, " ", KE, " ", postStepPoint.GetPosition().x
		# kinetic energy in MeV - PRE
		# initialKE = preStepPoint.GetKineticEnergy() 
		# kinetic energy in MeV - POST
		# finalKE = postStepPoint.GetKineticEnergy()
		# if particleName == 'e+':
		# 	p_test = [step.GetDeltaPosition().x,step.GetDeltaPosition().y,step.GetDeltaPosition().z]
		# 	p = [postStepPoint.GetPosition().x, postStepPoint.GetPosition().y, postStepPoint.GetPosition().z] # (mm)
		# 	# p and p_test are the SAME 
		# 	t = step.GetDeltaTime()
		# 	# t = track.GetGlobalTime() # (ns)
		# 	# t and t_test are the SAME
		# 	m = [postStepPoint.GetMomentum().x, postStepPoint.GetMomentum().y, postStepPoint.GetMomentum().z]
		# 	# m = [step.GetDeltaMomentum().x, step.GetDeltaMomentum().y, step.GetDeltaMomentum().z] # equal to the postStepPoint momentum
		# 	mm = np.sqrt((m[0])**2 + (m[1])**2 + (m[2])**2)
		# 	# momenta - PRE
		# 	initialMomentum = [preStepPoint.GetMomentum().x, preStepPoint.GetMomentum().y, preStepPoint.GetMomentum().z]
		# 	# momenta - POST
		# 	# print KE, "\n", p, "\n", initialMomentum, "\n", finalMomentum, "\n\n" 
		# 	# energy = step.GetTotalEnergyDeposit()
		# 	DA.dataCollection(p, m, t) # calls data collection and analysis on final positions and momenta
		# 	# return initialMomentum, finalMomentum 

### Raw G4 imports

from Geant4 import *
from Geant4 import FTFP_BERT, FTFP_BERT_HP, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP, QGSP_BIC_HP
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex

### G4PY imports

import g4py.ExN03geom
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

### python imports
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
plt.rc('font',family='Times New Roman')
import numpy as np
import pymouse
import time
### python file imports

from visualizer import Visualizer as VIS

### global variables
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
global ionIonizNum
ionIonizNum = []
global elecIonizNum
elecIonizNum = []

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
		# global ion_ioniz_num
		# global elec_ioniz_num
		ion_ioniz_num = len(ionIonizNum) 
		elec_ioniz_num = len(elecIonizNum)
		print ion_ioniz_num, elec_ioniz_num
		# time.sleep(2)
		ionIonizNum[:] = []
		elecIonizNum[:] = []
		

		return ion_ioniz_num, elec_ioniz_num

class DataAnalysis(object):

	'''
	Analyzes the stopping ranges of ions of different 
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
		particleName = ion, e-
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
				ionIonizNum.append("ion")
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
SGD = SetGlobalData()
GGD = GetGlobalData()

'''																									     											   '''
   #                         								   ____________   _  ____________ 								                           #
   ########################################################## / ___/ __/ _ | / |/ /_  __/ / / ##########################################################
   ##########################################################/ (_ / _// __ |/    / / / /_  _/ ##########################################################
   ##########################################################\___/___/_/ |_/_/|_/ /_/   /_/   ##########################################################
   #																																				   #
'''																																					   '''

def setup():
	# ------------------------------------------------------------------
	# setup for materials
	# ------------------------------------------------------------------
	# simple materials for Qgeom
	g4py.NISTmaterials.Construct()

	# ------------------------------------------------------------------
	# setup for geometry
	# ------------------------------------------------------------------
	#g4py.Qgeom.Construct()
	g4py.ezgeom.Construct()  # initialize

	# ------------------------------------------------------------------
	# setup for physics list
	# ------------------------------------------------------------------
	# g4py.EMSTDpl.Construct()
	physicsList = QGSP_BERT()
	physicsList = FTFP_BERT()
	gRunManager.SetUserInitialization(physicsList)

	# ------------------------------------------------------------------
	# setup for primary generator action
	# ------------------------------------------------------------------
	g4py.ParticleGun.Construct()

	# ------------------------------------------------------------------
	# Resize world
	# ------------------------------------------------------------------
	# prop_gas = G4Material("prop_cesium", 55., 132.90545*g/mole, .002*g/cm3)

def set_prop(prop_gas):
	g4py.ezgeom.SetWorldMaterial(prop_gas)
	g4py.ezgeom.ResizeWorld(5000000.*m, 5000000.*m, 5000000.*m)


class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
  "My Primary Generator Action"

  def __init__(self):
	G4VUserPrimaryGeneratorAction.__init__(self)
	self.particleGun= G4ParticleGun(1)

  def GeneratePrimaries(self, event):
	#dx= random.gauss(0., 0.1)
	dx=0.
	self.particleGun.SetParticleMomentumDirection(G4ThreeVector(dx, 0., 1.))
	self.particleGun.GeneratePrimaryVertex(event)

# ------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
  "My Run Action"

  def BeginOfRunAction(self, run):
	print "*** #event to be processed (BRA)=",
	run.GetNumberOfEventToBeProcessed()

  def EndOfRunAction(self, run):
	print "*** run end run(ERA)=", run.GetRunID()

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
  "My Event Action"

  #def BeginOfEventAction(self, event):
	#print "*** current event (BEA)=", event.GetEventID()
  #  pass

  def EndOfEventAction(self, event):
   # print "*** current event (EEA)=", event.GetEventID()

	for pv_idx in range(event.GetNumberOfPrimaryVertex()):
		pvtx = event.GetPrimaryVertex(pv_idx)
		for pp_idx in range(pvtx.GetNumberOfParticle()):
			ppart = pvtx.GetPrimary(pp_idx)
			pmass = ppart.GetMass()
			print pmass
			# time.sleep(1)
# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
  "My Stepping Action"

  def UserSteppingAction(self, step):
	#print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
	# track= step.GetTrack()
	# touchable= track.GetTouchable()
	# pv= touchable.GetVolume()

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

	# print "process", procName, ionizingEnergy, deltaEnergy

	energyLevel = "mean"
	# energyLevel = 1
	DA.calc_ionization(procName, particleName, deltaEnergy, ionizingEnergy, energyLevel)

	#print pv.GetCopyNo()
	#print touchable.GetReplicaNumber(0)

# ------------------------------------------------------------------

# ==================================================================
# main
# ==================================================================

# qDC= gtest01.QDetectorConstruction()
# gRunManager.SetUserInitialization(qDC)

# qPL= gtest01.QPhysicsList()
# gRunManager.SetUserInitialization(qPL)


# set user actions...
#qPGA= gtest01.QPrimaryGeneratorAction()
# myPGA= MyPrimaryGeneratorAction()
# gRunManager.SetUserAction(myPGA)

#gRunManager.SetUserAction(myRA)

# myEA= MyEventAction()
# gRunManager.SetUserAction(myEA)

# mySA= MySteppingAction()
# gRunManager.SetUserAction(mySA)

# set particle gun
#ApplyUICommand("/control/execute gun.mac")
#pg= qPGA.GetParticleGun()
# pg= myPGA.particleGun
# pg.SetParticleByName("GenericIon")
# pg.SetParticleEnergy(200.*MeV)
# pg.SetParticlePosition(G4ThreeVector(0.,0.,0)*cm)

# magnetic field
# fieldMgr= gTransportationManager.GetFieldManager()
# myField= G4UniformMagField(G4ThreeVector(0.,10.*tesla,0.))
# #myField= MyField()
# fieldMgr.SetDetectorField(myField)
# fieldMgr.CreateChordFinder(myField)


# beamOn
def main():

	prop_list = ["cesium", "bismuth","mercury","xenon","iodine"]
	density_list = list(np.arange(0.02,.015, .0001))
	prop_dict = {"cesium" : 	(55., 132.90545), 
				 "bismuth" : 	(83., 208.9804), 
				 "mercury" : 	(80., 200.59), 
				 "xenon" : 		(54., 131.293), 
				 "iodine" : 	(53., 126.90447)}


	setup()
	gRunManager.Initialize()
	for prop in ["cesium", "bismuth","mercury","xenon","iodine"]:
		SGD.set_ionization_energies(prop)
		
		ion_ioniz_filename = prop + "_data/mean_ionization_Rb98.txt" # <<<<<<<------- mean or 1st IE
		# os.system("rm " + ion_ioniz_filename)
		os.system("touch " + ion_ioniz_filename)
		# time.sleep(5)
		# open(ion_ioniz_filename).close()
		ion_ioniz_file = open(ion_ioniz_filename, "a")

		density = 0.00000001

		prop_gas = G4Material(prop + "_" + str(density), 
					  prop_dict[prop][0], 
					  prop_dict[prop][1]*g/mole, 
					  density*g/cm3)
		set_prop(prop_gas)

		VIS().visualizer(45, 45, "viewer")
		# gRunManager.BeamOn(1)

		myEA= MyEventAction()
		gRunManager.SetUserAction(myEA)

		mySA= MySteppingAction()
		gRunManager.SetUserAction(mySA)

		gControlExecute("gun.mac")

		ion_ioniz, elec_ioniz = GGD.get_total_ioniz_num()
		total_ioniz = ion_ioniz + elec_ioniz
		ion_ioniz_file.writelines(str(total_ioniz)+"\n")
		print ion_ioniz
			# time.sleep(2)

if __name__ == '__main__':
	main()
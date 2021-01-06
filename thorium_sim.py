#----------GEANT4 imports----------#
from Geant4 import *
import g4py.ExN03geom
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
from Geant4 import FTFP_BERT, FTFP_BERT_HP, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP, QGSP_BIC_HP
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex

#----------PYTHON imports----------#
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
plt.rc('font',family='Times New Roman')
import numpy as np
import pymouse
import random
import time
#-----------FILE imports------------#
from geom_constructor import GeomConstructor 
from visualizer import Visualizer
from read_gdml_test import MyDetectorConstruction
from alpha_beam_2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, DataAnalysis, GetGlobalData, SetGlobalData

GC = GeomConstructor()
VIS = Visualizer()
DA = DataAnalysis()
GGD = GetGlobalData()
SGD = SetGlobalData()

######### CODE STARTS HERE ##########
# physicsList = QGSP_BERT()
# physicsList = FTFP_BERT_HP()

# ==================================================================
# intialize
# ==================================================================
def Configure():
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
	# gRunManager.SetUserInitialization(physicsList)

	# ------------------------------------------------------------------
	# setup for primary generator action
	# ------------------------------------------------------------------
	g4py.ParticleGun.Construct()

	# ------------------------------------------------------------------
	# Resize world
	# ------------------------------------------------------------------
	g4py.ezgeom.ResizeWorld(500.*mm, 500.*mm, 500.*mm)

	gApplyUICommand("/gun/number 1")
	gApplyUICommand("/gun/particle e-")
	gApplyUICommand("/gun/energy 0.2 MeV")
	gApplyUICommand("/gun/direction 0. 0.  1.")
	gApplyUICommand("/gun/position 0. -5. -30. cm")


class SpaceConstructor(object):
	def __init__(self):
		g4py.NISTmaterials.Construct()
		g4py.ezgeom.Construct()
		exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(exN03PL)
		xe = gNistManager.FindOrBuildMaterial("G4_Xe") # USE THIS MATERIAL MANAGER
		au = G4Material.GetMaterial("G4_Au", 1)

		# density = 5.56e-2
		radioisotope = G4Material("U235", 92,  235.0439299*g/mole, 19.1*g/cm3)
		material1 = gNistManager.FindOrBuildMaterial("G4_C")

		g4py.ezgeom.ResizeWorld(500.*mm, 500.*mm, 500.*mm)



# -----------CLASS ASSIGNMENTS----------- #

SC = SpaceConstructor()


class Main(object):

	def __init__(self, viewer_name):
		self.stoppingRanges = []
		self.name = viewer_name

	def run(self, densities, particleEnergy, particle, particleCount, viz_theta, viz_phi, prop):
		self.stoppingRanges = []
		for d in densities:
			DA
			print(d)
			if prop == "cesium":
				prop_gas = G4Material("prop_cesium", 55., 132.90545*g/mole, d*g/cm3)
			if prop == "bismuth":
				prop_gas = G4Material("prop_bismuth", 83., 208.9804*g/mole, d*g/cm3)
			if prop == "mercury":
				prop_gas = G4Material("prop_mercury", 80., 200.59*g/mole, d*g/cm3)
			if prop == "iodine":
				prop_gas = G4Material("prop_iodine", 53., 126.90447*g/mole, d*g/cm3)
			if prop == "xenon":
				prop_gas = G4Material("prop_xenon", 54., 131.293*g/mole, d*g/cm3)
			if prop == "lithium":
				prop_gas = G4Material("prop_lithium", 3., 6.941*g/mole, d*g/cm3)

			print prop, "\n"

			# vacuum = gNistManager.FindOrBuildMaterial("vacuum") 
			# g4py.ezgeom.SetWorldMaterial(vacuum)
			SGD.set_ionization_energies(prop)
			g4py.ezgeom.SetWorldMaterial(prop_gas)
			# alphaEnergy = 5.49 # MeV

			# gApplyUICommand("/run/initialize")
			gRunManager.Initialize()

			# PGA_1 = MyPrimaryGeneratorAction(particle, particleEnergy, MeV, particleCount, [1,0,0])
			# gRunManager.SetUserAction(PGA_1)

			for i in range(0, particleCount): # creates random momentum vectors originating from [0, 0, 0]
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

				gRunManager.BeamOn(1)
			myEA = MyEventAction()
			gRunManager.SetUserAction(myEA)
			mySA = MySteppingAction()
			gRunManager.SetUserAction(mySA)
			VIS.visualizer(viz_theta, viz_phi, self.name)
			myRA = MyRunAction()
			gRunManager.SetUserAction(myRA)

			gRunManager.Initialize()

			# stoppingRange = GGD.get_stopping_range()
			# SR_file.writelines(str(stoppingRange)+"\n")

			# elec_num = GGD.get_electron_num()
			# elec_file.writelines(str(elec_num)+"\n")

			alpha_ioniz, elec_ioniz = GGD.get_total_ioniz_num()
			# print "Total num of ionizations = ", alpha_ioniz + elec_ioniz
			print alpha_ioniz, elec_ioniz
			# total_ioniz_file.writelines(str(alpha_ioniz)+ " " + str(elec_ioniz)+"\n")
			
			# self.stoppingRanges.append(stoppingRange)

		return self.stoppingRanges

		# print("")
		# print("Density = ", density)
		# print("SR = ", stoppingRange)
		# print(" ")

	def plotter(self, densities, energies, SRList):

		# print(self.stoppingRanges)
		# print(densities)

		xs = densities
		ys = energies
		zs = SRList

		print(len(xs),len(ys), len(zs[0]))

		xlabel = r'$\rho$' + "(g/cm" r'$^3$' + ")"
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
		ax = fig.add_subplot(111, projection='3d')

		axes = plt.gca()
		axes.set_xlim([min(xs),max(xs)])
		axes.set_ylim([min(ys),max(ys)])
		axes.set_zlim([0,50])
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_zlabel(zlabel)

		print(xs)
		print(ys)
		print(zs)

		for i in np.arange(0, len(ys)): # energies is on the y
			ax.plot(xs, list(ys[i]*np.ones(len(xs))), zs[i])



		plt.title(zlabel + " vs " + xlabel + "vs" + ylabel)
		plt.grid()
		plt.show()


if __name__ == '__main__':

	# setup
	Configure()

	# run
	gRunManager.Initialize()

	propellants = ["cesium", "bismuth","mercury","xenon","iodine"]
	# propellants = ["lithium"]
	energyRange = [4.9,5.1,5.6,5.8,5.5] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241

	radioactivities = [0.63e12,21.8e12,0.643e12,3.03e12,0.126e12] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241
	# densityRange = list(np.arange(0.00000001,.005,0.0001)) # g/cm3

	# densityRange =  list(np.arange(0.006,.05, .0001))
	densityRange_2 = list(np.arange(0.003,.015, .0001))
	# densityRange_2 = list(np.arange(0.003,.03, .0001))
	# densityRange_2 = [0.01]

	stoppingRangesList = []

	viewer_name = "alpha_beam"
	particle = "GenericIon"
	# particle = "alpha"
	# particle = "e+"
	# particle = "proton"




	# G4ParticleGun(1).SetParticleByName("GenericIon")

	# gApplyUICommand("/gun/particle ion")
	# gApplyUICommand("/gun/ion 92 233 0")
	# gApplyUICommand("/gun/energy 90 MeV")
	# gApplyUICommand("/gun/position 0. 0. 0. m")

	time.sleep(2)

	MN = Main(viewer_name)
	GGD

	SC
	for prop in propellants:

		## Holds stopping range data for each RI
		# SR_filename = prop + "2/"+str(e)+"_SR.txt"
		# os.system("touch " + SR_filename)
		# open(SR_filename).close()
		# SR_file = open(SR_filename, "a")

		## Holds number of electrons generated within the SR radius for each RI
		# elec_num_filename = prop + "2/"+str(e)+"_elec_num.txt"
		# os.system("touch " + elec_num_filename)
		# open(elec_num_filename).close()
		# elec_file = open(elec_num_filename, "a")

		# Holds number of alpha and e- induced ionizations
		e = 90 # MeV
		theta = 45
		phi = 45
		# total_ionization_num_filename = prop + "2/"+str(e)+"_ionization_num.txt"
		# os.system("touch " + total_ionization_num_filename)
		# open(total_ionization_num_filename).close()
		# total_ioniz_file = open(total_ionization_num_filename, "a")

		particleCount = 50
		MN
		stoppingRanges = MN.run(densityRange_2, e, particle, particleCount, theta, phi, prop)
		stoppingRangesList.append(stoppingRanges)
	# time.sleep(.005)


	# MN.plotter(densityRange, energyRange, stoppingRangesList)
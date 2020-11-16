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
import time
#-----------FILE imports------------#
from geom_constructor import GeomConstructor 
from visualizer import Visualizer
from read_gdml_test import MyDetectorConstruction
from alpha_beam_1 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, DataAnalysis, GetGlobalData, SetGlobalData

GC = GeomConstructor()
VIS = Visualizer()
DA = DataAnalysis()
GGD = GetGlobalData()
SGD = SetGlobalData()

######### CODE STARTS HERE ##########
physicsList = QGSP_BERT()
# physicsList = FTFP_BERT_HP()

# mouse stuff
mouse = pymouse.PyMouse()
mouse_x = [0, mouse.position()[0]]
mouse_y = [0, mouse.position()[1]]

theta = 45
phi = 45

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
		# GC.ConstructBox("carbon_plate", material1, [0,0,0], mm,[1,20,20])
		g4py.ezgeom.ResizeWorld(500.*mm, 500.*mm, 500.*mm)
		# GC.ConstructBox("uranium", radioisotope, [-49.5/2,0,0], mm, [1,49.5,49.5])


# -----------CLASS ASSIGNMENTS----------- #

SC = SpaceConstructor()

# myDC = MyDetectorConstruction()
# gRunManager.SetUserInitialization(myDC)

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
				prop_gas = G4Material("gaseousProp", 54., 131.293*g/mole, d*g/cm3)

			print prop, "\n"

			# vacuum = gNistManager.FindOrBuildMaterial("vacuum") 
			# g4py.ezgeom.SetWorldMaterial(vacuum)
			SGD.set_ionization_energies(prop)
			g4py.ezgeom.SetWorldMaterial(prop_gas)
			# alphaEnergy = 5.49 # MeV
			PGA_1 = MyPrimaryGeneratorAction(particle, particleEnergy, MeV, particleCount, [1,0,0])
			gRunManager.SetUserAction(PGA_1)
			myEA = MyEventAction()
			gRunManager.SetUserAction(myEA)
			mySA = MySteppingAction()
			gRunManager.SetUserAction(mySA)
			VIS.visualizer(viz_theta, viz_phi, self.name)
			myRA = MyRunAction()
			gRunManager.SetUserAction(myRA)
			gRunManager.Initialize()
			gRunManager.BeamOn(1)

			# stoppingRange = GGD.get_stopping_range()
			# SR_file.writelines(str(stoppingRange)+"\n")

			# elec_num = GGD.get_electron_num()
			# elec_file.writelines(str(elec_num)+"\n")

			alpha_ioniz, elec_ioniz = GGD.get_total_ioniz_num()
			# print "Total num of ionizations = ", alpha_ioniz + elec_ioniz
			print alpha_ioniz, elec_ioniz
			total_ioniz_file.writelines(str(alpha_ioniz)+ " " + str(elec_ioniz)+"\n")
			
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

		# xs = [0.001, 0.011, 0.020999999999999998, 0.030999999999999996, 0.040999999999999995, 0.05099999999999999, 0.06099999999999999, 0.071, 0.08099999999999999, 0.09099999999999998]
		# ys = [0.1, 4.6 ,9.1]
		# zs = [[0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.5828279651355435, 0.19427598837851448, 0.14570699128388587, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.04856899709462862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

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

	# energyRange = np.arange(0.1,10,.1) # MeV
	# propellants = ["cesium"]
	# energyRange = [5.5]
	propellants = ["cesium", "bismuth","mercury","xenon","iodine"]
	energyRange = [4.9,5.1,5.6,5.8,5.5] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241

	# energyRange = [5.3,4.9,5.1,5.6,5.8,5.5] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241
	# radioactivities = [166e12,0.63e12,21.8e12,0.643e12,3.03e12,0.126e12] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241

	# energyRange = [5.1,5.6,5.8,5.5] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241
	radioactivities = [0.63e12,21.8e12,0.643e12,3.03e12,0.126e12] # Bcq Po210, Po209, Po208, Pu238, Cm244, Am241
	# densityRange = list(np.arange(0.00000001,.005,0.0001)) # g/cm3

	# densityRange =  list(np.arange(0.006,.05, .0001))
	densityRange_2 = list(np.arange(0.003,.015, .0001))
	# densityRange_2 = list(np.arange(0.003,.03, .0001))
	# densityRange_2 = [0.01]

	stoppingRangesList = []

	gRunManager.SetUserInitialization(physicsList)
	# gApplyUICommand("/process/em/fluo true")
	# gApplyUICommand("/process/em/auger true")
	# gApplyUICommand("/process/em/augerCascade truegApplyUICommand")
	# gApplyUICommand("/process/em/pixe true")
	viewer_name = "alpha_beam"
	# particle = "GenericIon"
	particle = "alpha"
	# particle = "e+"
	# particle = "proton"
	MN = Main(viewer_name)
	GGD
	# gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
	# gApplyUICommand("/vis/viewer/create OGLSX " + viewer_name)

	SC
		# theta += 10
		# phi += 10

		#### 	MOUSE CONTROL OF GEANT4 SPACE 	####
		# mouse_x[0] = mouse_x[-1]
		# mouse_y[0] = mouse_y[-1]
		# mouse_x[-1] = mouse.position()[0]
		# mouse_y[-1] = mouse.position()[1]

		# dphi = 0.5*(mouse_x[-1] - mouse_x[0])
		# dtheta = 0.5*(mouse_y[-1] - mouse_y[0])

		# theta += dtheta
		# phi += dphi
		############################################


	# for e in energyRange:

	# energy = energyRange[0]


	# energy = 5.49

	for prop in propellants:

		for i in range(0,len(energyRange)):

			e = energyRange[i]

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
			total_ionization_num_filename = prop + "2/"+str(e)+"_first_ionization_num.txt"
			os.system("touch " + total_ionization_num_filename)
			open(total_ionization_num_filename).close()
			total_ioniz_file = open(total_ionization_num_filename, "a")

			# e = 5.49
			# rad = radioactivities[i]
			# time_elapsed = 1e-9
			# particleCount = int(rad*time_elapsed)
			particleCount = 50
			MN
			stoppingRanges = MN.run(densityRange_2, e, particle, particleCount, theta, phi, prop)
			stoppingRangesList.append(stoppingRanges)
	# time.sleep(.005)


	# MN.plotter(densityRange, energyRange, stoppingRangesList)
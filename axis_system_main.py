#----------GEANT4 imports----------#
from Geant4 import *
import g4py.ExN03geom
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
from Geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QGSP_BIC_AllHP, QGSP_BIC_HP
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
from alpha_beam_1 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, DataAnalysis, GetGlobalData

GC = GeomConstructor()
VIS = Visualizer()
DA = DataAnalysis()
GGD = GetGlobalData()

# exN03geom= g4py.ExN03geom.ExN03DetectorConstruction()
# gRunManager.SetUserInitialization(exN03geom)

######### CODE STARTS HERE ##########

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
		g4py.ezgeom.ResizeWorld(50.*mm, 50.*mm, 50.*mm)
		# GC.ConstructBox("uranium", radioisotope, [-49.5/2,0,0], mm, [1,49.5,49.5])


# -----------CLASS ASSIGNMENTS----------- #

SC = SpaceConstructor()

# myDC = MyDetectorConstruction()
# gRunManager.SetUserInitialization(myDC)

class Main(object):

	def __init__(self, viewer_name):
		self.stoppingRanges = []
		self.name = viewer_name

	def run(self, densities, particleEnergy, particle, viz_theta, viz_phi):
		self.stoppingRanges = []
		for d in densities:

			xe_highPressure = G4Material("gaseousXenon", 54., 131.293*g/mole, d*g/cm3)
			g4py.ezgeom.SetWorldMaterial(xe_highPressure)
			# alphaEnergy = 5.49 # MeV
			PGA_1 = MyPrimaryGeneratorAction(particle, particleEnergy, MeV, 100, [1,0,0])
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

			stoppingRange = GGD.get_stopping_range()
			self.stoppingRanges.append(stoppingRange)
			file.writelines(str(stoppingRange)+"\n")

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
		# plt.draw() 
		plt.show()


if __name__ == '__main__':

	# energyRange = np.arange(0.1,10,.1) # MeV
	energyRange = [5.3,4.9,5.1,5.6,5.8,5.5]
	densityRange = list(np.arange(0.000001,.1,0.001)) # g/cm3
	stoppingRangesList = []

	physicsList = QGSP_BIC_HP()
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
	DA
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

	energy = energyRange[5]

	open("data/"+str(energy)+".txt").close()
	file = open("data/"+str(energy)+".txt", "a")

	MN
	stoppingRanges = MN.run(densityRange, energy, particle, theta, phi)
	stoppingRangesList.append(stoppingRanges)
	# time.sleep(.005)


	MN.plotter(densityRange, energyRange, stoppingRangesList)
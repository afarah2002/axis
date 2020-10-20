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
from alpha_beam_1 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction

GC = GeomConstructor()
VIS = Visualizer()

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

		density = .1
		xe_highPressure = G4Material("gaseousXenon", 54., 131.293*g/mole, density*g/cm3)
		radioisotope = G4Material("U235", 92,  235.0439299*g/mole, 19.1*g/cm3)
		g4py.ezgeom.SetWorldMaterial(xe_highPressure)

		g4py.ezgeom.ResizeWorld(50.*mm, 50.*mm, 50.*mm)
		# GC.ConstructBox("uranium", radioisotope, [-49.5/2,0,0], mm, [1,49.5,49.5])


# -----------CLASS ASSIGNMENTS----------- #
SC = SpaceConstructor()


# myDC = MyDetectorConstruction()
# gRunManager.SetUserInitialization(myDC)

class Main(object):

	def __init__(self, viewer_name):
		self.name = viewer_name

	def run(self, particle, viz_theta, viz_phi):

		PGA_1 = MyPrimaryGeneratorAction(particle, 5.3, MeV, 500, [1,0,0])
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





if __name__ == '__main__':


	physicsList = QGSP_BIC_HP()
	gRunManager.SetUserInitialization(physicsList)
	gApplyUICommand("/process/em/fluo true")
	gApplyUICommand("/process/em/auger true")
	gApplyUICommand("/process/em/augerCascade truegApplyUICommand")
	gApplyUICommand("/process/em/pixe true")
	viewer_name = "alpha_beam"
	# particle = "GenericIon"
	particle = "alpha"
	# particle = "e+"
	MN = Main(viewer_name)
	# gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
	# gApplyUICommand("/vis/viewer/create OGLSX " + viewer_name)

	while True:
		SC
		# theta += 10
		# phi += 10
		mouse_x[0] = mouse_x[-1]
		mouse_y[0] = mouse_y[-1]
		mouse_x[-1] = mouse.position()[0]
		mouse_y[-1] = mouse.position()[1]

		dphi = 0.5*(mouse_x[-1] - mouse_x[0])
		dtheta = 0.5*(mouse_y[-1] - mouse_y[0])


		theta += dtheta
		phi += dphi

		MN
		MN.run(particle, theta, phi)
		# time.sleep(.005)
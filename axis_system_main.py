#----------GEANT4 imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
#----------PYTHON imports----------#
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
plt.rc('font',family='Times New Roman')
import numpy as np
import time
#-----------FILE imports------------#
from geom_constructor import GeomConstructor 
from visualizer import Visualizer
from read_gdml_test import MyDetectorConstruction


######### CODE STARTS HERE ##########
class SpaceConstructor(object):
	def __init__(self):
		g4py.NISTmaterials.Construct()
		g4py.ezgeom.Construct()
		exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(exN03PL)
		# xenon = G4Material.GetMaterial("G4_XE", 1)
		# g4py.ezgeom.SetWorldMaterial(xenon)
		g4py.ezgeom.ResizeWorld(20.*cm, 20.*cm, 20.*cm)
		# gold = G4Material.GetMaterial("G4_Au", 1)
		# GC.ConstructSphere("sphere1", 
		# 					gold, 
		# 					[10,10,10], 
		# 					cm, 
		# 					0, 
		# 					10, 
		# 					0, 
		# 					360, 
		# 					0, 
		# 					360)

# -----------CLASS ASSIGNMENTS----------- #
SC = SpaceConstructor()
GC = GeomConstructor()
VIS = Visualizer()


myDC = MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

class Main(object):

	def __init__(self, viz_theta, viz_phi, viewer_name):
		self.viz_theta = viz_theta
		self.viz_phi = viz_phi
		self.name = viewer_name
		VIS.visualizer(self.viz_theta, self.viz_phi, self.name)
		gRunManager.Initialize()


if __name__ == '__main__':
	SC
	theta = 45
	phi = 45
	viewer_name = "main"
	gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
	gApplyUICommand("/vis/viewer/create OGLSX " + viewer_name)
	while True:
		theta += 1
		phi += 1
		Main(theta, phi, viewer_name)
		time.sleep(.005)
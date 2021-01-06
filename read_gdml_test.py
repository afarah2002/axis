#!/usr/bin/python
# ==================================================================
# An example for reading a GDML file
#
# ==================================================================
from Geant4 import *
import g4py.ezgeom
import g4py.ExN01pl, g4py.ParticleGun
from geom_constructor import GeomConstructor 
import glob

GC = GeomConstructor()


# ==================================================================
# user actions in python
# ==================================================================
class MyDetectorConstruction(G4VUserDetectorConstruction):
	"My Detector Construction"

	def __init__(self):
		G4VUserDetectorConstruction.__init__(self)
		self.world= None
		self.gdml_parser= G4GDMLParser()

	# -----------------------------------------------------------------
	def __del__(self):
		pass
		
	# -----------------------------------------------------------------
	def Construct(self):
		# for gdml in glob.glob("	*m.gdml"):
		gdml = "../gdmls/fullmodelV10102inletbulkheadassembly10101propellantplenum2010103electricalisolator1_Aluminum.gdml"
			# print gdml
		self.gdml_parser.Read(gdml)
		# gdml = "fullmodelV10202ionizerassembly1020201hollowdischargeinjectioncathode1_Aluminum.gdml"
		# self.gdml_parser.Read("bulkhead.gdml")
		self.world = self.gdml_parser.GetWorldVolume()
		# xenon = G4Material.GetMaterial("G4_Xe")
		# au = G4Material.GetMaterial("G4_Au", 1)
		# GC.ConstructBox("xenon",
		# 								au,
		# 								[10,0,0],
		# 								cm,
		# 								[1,5,5])
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

		return self.world


# ==================================================================
# main
# ==================================================================
# set geometry
myDC= MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)
g4py.ezgeom.ResizeWorld(10.*m, 10.*m, 10.*m)

# minimal physics list
g4py.ExN01pl.Construct()

# set primary generator action
g4py.ParticleGun.Construct()

# initialize
gRunManager.Initialize()

gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
gApplyUICommand("/vis/open OGLSX")
t = 45.
p = 45.
while True:
	t += .5
	p += .5
	# visualization
	gApplyUICommand("/vis/scene/create")
	gApplyUICommand("/vis/scene/add/volume")
	gApplyUICommand("/vis/sceneHandler/attach")
	gApplyUICommand("/vis/viewer/set/viewpointThetaPhi " + str(t) + " " + str(p))

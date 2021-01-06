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

import numpy as np
import time 

### python file imports

from visualizer import Visualizer as VIS

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
	prop_gas = G4Material("prop_cesium", 55., 132.90545*g/mole, .002*g/cm3)
	g4py.ezgeom.SetWorldMaterial(prop_gas)
	g4py.ezgeom.ResizeWorld(500.*mm, 500.*mm, 500.*mm)


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
			time.sleep(1)
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

	#print pv.GetCopyNo()
	#print touchable.GetReplicaNumber(0)

# ------------------------------------------------------------------
class MyField(G4MagneticField):
  "My Magnetic Field"

  def GetFieldValue(self, pos, time):
	bfield= G4ThreeVector()
	bfield.x= 0.
	bfield.y= 5.*tesla
	bfield.z= 0.
	return bfield

# ==================================================================
# main
# ==================================================================

# qDC= gtest01.QDetectorConstruction()
# gRunManager.SetUserInitialization(qDC)

# qPL= gtest01.QPhysicsList()
# gRunManager.SetUserInitialization(qPL)

setup()

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

gRunManager.Initialize()

# beamOn

for i in np.arange(0,10000,1):
	VIS().visualizer(45, 45, "viewer")
	# gRunManager.BeamOn(1)

	myEA= MyEventAction()
	gRunManager.SetUserAction(myEA)

	mySA= MySteppingAction()
	gRunManager.SetUserAction(mySA)

	gControlExecute("gun.mac")
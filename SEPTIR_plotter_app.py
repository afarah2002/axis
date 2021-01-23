#-----------------imports-----------------#
import numpy as np
from tkinter import *

#-----------------app, main class-----------------#
from SEPTIR_plotter import MaterialAssignment


class App(Frame):

	def init(self):
		#-----------------title-----------------#
		title = Label(self, text="FIRESTORM: SEPTIR Performance Plotter")
		title.pack()

		#-----------------create RI and propellant dropdown-----------------#

		#-----------------set RI and propellant attributes-----------------#
		radioisotopesList, propellantsList = MaterialAssignment().assign()
		self.radioisotope = radioisotopesList[0]
		self.propellant = propellantsList[0]

		#-----------------get list of RI and prop names for labelling-----------------#
		radioisotopeNames = [RI.name for RI in radioisotopesList]
		propellantNames = [prop.name for prop in propellantsList]

		#-----------------name-wise dictionaries of RIs and Props-----------------#
		self.radioisotopeDict = dict(zip(radioisotopeNames, radioisotopesList))
		self.propellantDict = dict(zip(propellantNames, propellantsList))

		#-----------------drop downs-----------------#
		self.radioisotopeLabel = StringVar(self)
		self.propellantLabel = StringVar(self)
		self.radioisotopeLabel.set(radioisotopeNames[0]) 
		self.propellantLabel.set(propellantNames[0])

		RI_dropdown = OptionMenu(self, self.radioisotopeLabel, *radioisotopeNames)
		RI_dropdown.config(width=90, font=('Helvetica', 12))

		prop_dropdown = OptionMenu(self, self.propellantLabel, *propellantNames)
		prop_dropdown.config(width=90, font=('Helvetica', 12))

		RI_dropdown.pack()
		prop_dropdown.pack()

		#-----------------Button: Get RI-----------------#
		print_RI_Button = Button(self, text="Print RI", command=self.get_radioisotope)
		print_RI_Button.pack()


	def get_radioisotope(self):
		self.radioisotope = self.radioisotopeDict.get(self.radioisotopeLabel.get())

	def get_propellant(self):
		self.propellant = self.propellantDict.get(self.propellantLabel.get())



	def __init__(self, master=None):
		Frame.__init__(self, master)
		self.init()
		self.grid()

#-----------------run app-----------------#
app = App()
app.mainloop()
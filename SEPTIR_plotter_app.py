#-----------------imports-----------------#
import numpy as np
from tkinter import *

#-----------------app, main class-----------------#
from SEPTIR_plotter import MaterialAssignment


class App(Frame):

	def init(self):
		#title and header    row=0, 1
		title = Label(self, text="FIRESTORM: SEPTIR Performance Plotter")
		title.pack()
		# title.grid(row=0, column=1, columnspan=5)
		header = Label(self, text="Measurement of Mass Attenuation Coefficient")
		# header.grid(row=1, column=1, columnspan=5)

		radioisotopesList, propellantsList = MaterialAssignment().assign()
		radioisotope = StringVar(self)
		radioisotope.set(radioisotopesList[0])


		RI_dropdown = OptionMenu(self, radioisotope, *radioisotopesList)
		RI_dropdown.config(width=90, font=('Helvetica', 12))

		RI_dropdown.pack()

	def __init__(self, master=None):
		Frame.__init__(self, master)
		self.init()
		self.grid()

#-----------------run app-----------------#
app = App()
app.mainloop()
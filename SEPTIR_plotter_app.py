#-----------------imports-----------------#
import numpy as np
from tkinter import *
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.backends.backend_tkagg
from matplotlib.figure import Figure

#-----------------app, main class-----------------#
from SEPTIR_plotter import MaterialAssignment, DataProcessor, Constants


class App(Frame):

	def init(self):
		#-----------------title-----------------#
		title = Label(self, text="FIRESTORM: SEPTIR Performance Plotter")
		title.grid(row=1, column=1, sticky=W)

		#-----------------set up plot-----------------#
		self.f = Figure(figsize=(5,5), dpi= 100)
		self.a = self.f.add_subplot(111)
		self.a.grid()
		#-----------------add plot-----------------#
		self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.f, self)
		self.canvas.draw()
		self.canvas.get_tk_widget().grid(row=3, column=6)
		#-----------------add MPL toolbar-----------------#
		# toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2TkAgg(canvas, self)
		# toolbar.update()
		# canvas._tkcanvas.grid(row=4, column=5)	
		#-----------------create RI and propellant dropdown-----------------#


		#-----------------get list of RI and prop names for labelling-----------------#
		radioisotopesList, propellantsList = MaterialAssignment().assign()


		#-----------------name-wise dictionaries of RIs and Props-----------------#




		#-----------------drop downs-----------------#
		#-----------------RI drop down-----------------#

		#-----------------set RI and propellant attributes-----------------#
		radioisotopeNames = [RI.name for RI in radioisotopesList]
		self.radioisotopeDict = dict(zip(radioisotopeNames, radioisotopesList))

		self.radioisotopeLabel = StringVar(self)
		self.radioisotopeLabel.set(radioisotopeNames[0]) 
		RI_dropdown = OptionMenu(self, self.radioisotopeLabel, *radioisotopeNames)
		RI_dropdown.config(width=10, font=('Helvetica', 12))
		RI_dropdown.grid(row=2, column=1, sticky=W)

		#-----------------Propellant drop down-----------------#
		propellantNames = [prop.name for prop in propellantsList]
		self.propellantDict = dict(zip(propellantNames, propellantsList))

		self.propellantLabel = StringVar(self)
		self.propellantLabel.set(propellantNames[0])
		prop_dropdown = OptionMenu(self, self.propellantLabel, *propellantNames)
		prop_dropdown.config(width=10, font=('Helvetica', 12))
		prop_dropdown.grid(row=3, column=1, sticky=W)

		#-----------------YData drop down-----------------#
		#-----------------mdot or Fthrust-----------------#
		self.ydataLabel = StringVar(self)
		self.ydataLabel.set("")
		ydata_dropdown = OptionMenu(self, self.ydataLabel, *["mdot", "Fthrust"])
		ydata_dropdown.config(width=10, font=('Helvetica', 12))
		ydata_dropdown.grid(row=4, column=1, sticky=W)



		#-----------------Button: Plots mdot or thrust----------------#
		mdot_thrust_button = Button(self, text="Plot mdot/thrust", command=self.mdot_thrust)
		mdot_thrust_button.grid(row=9, column=6)

		#-----------------Button: mRI vs dv-----------------#
		plot_button = Button(self, text="Plot mRI vs dV", command=self.mri_vs_dv)
		plot_button.grid(row=10, column=6)
		
		#-----------------entry boxes-----------------#
		#-----------------spacecraft mass, g-----------------#
		spacecraft_mass_entry_Label = Label(self, text="Spacecraft mass (g)")
		spacecraft_mass_entry_Label.grid(row=2, column=2, sticky=W)
		self.spacecraft_mass_entry = Entry(self)
		self.spacecraft_mass_entry.grid(row=2, column=3, sticky=W)
		#-----------------RI initial mass, g-----------------#
		RI_launch_mass_entry_Label = Label(self, text="RI launch mass (g)")
		RI_launch_mass_entry_Label.grid(row=3, column=2, sticky=W)
		# self.RI_launch_mass_entry = Entry(self)
		# self.RI_launch_mass_entry.grid(row=3, column=3, sticky=W)
		#-----------------Isp, s-----------------#
		Isp_entry_Label = Label(self, text="Isp (s)")
		Isp_entry_Label.grid(row=4, column=2, sticky=W)
		self.Isp_entry = Entry(self)
		self.Isp_entry.grid(row=4, column=3, sticky=W)
		#-----------------delta v (m/s)-----------------#
		dv_entry_Label = Label(self, text="dv (m/s)")
		dv_entry_Label.grid(row=5, column=2, sticky=W)
		self.dv_entry = Entry(self)
		self.dv_entry.grid(row=5, column=3, sticky=W)
		#-----------------burn duration, years-----------------#
		burn_duration_entry_Label = Label(self, text="Burn duration (years)")
		burn_duration_entry_Label.grid(row=6, column=2, sticky=W)
		self.burn_duration_entry = Entry(self)
		self.burn_duration_entry.grid(row=6, column=3, sticky=W)
		#-----------------transit time, years-----------------#
		transit_time_entry_Label = Label(self, text="Transit time (years)")
		transit_time_entry_Label.grid(row=7, column=2, sticky=W)
		self.transit_time_entry = Entry(self)
		self.transit_time_entry.grid(row=7, column=3, sticky=W)
		#---------------------------------------------#


	def mdot_thrust(self):
		#-----------------processes a specific RI-propellant combination-----------------#
		self.update()
		self.a.plot(np.divide(self.burn_duration,3.1536e7), np.multiply(self.plotdata,1000))
		#-----------------plot formatting-----------------#
		xlabel = "Time (years)"
		if self.ydata == "mdot":
			ylabel = self.propellant.name.capitalize() + " " r'$\dot{m}_{P}$'+ " (mg/s)"
		if self.ydata == "Fthrust":
			ylabel = self.propellant.name.capitalize() + " " r'${F}_{T}$'+ " (mN)"

		self.a.set_xlabel(xlabel)
		self.a.set_ylabel(ylabel)
		self.a.set_title(ylabel + " vs " + xlabel)
		self.a.grid()
		self.canvas.draw()



	def mri_vs_dv(self):
		#-----------------plots mRI vs dV -----------------#
		self.update()
		dvList = np.arange(0,50000,500)
		mRI_launch_List = -np.multiply(np.exp(0.693*self.transit_time/self.radioisotope.HL)*Constants().N_A/(self.propellant.M*self.IE*self.radioisotope.SA*self.radioisotope.HL)* \
					self.S_m0/(np.exp(-0.693*(max(self.burn_duration)/self.radioisotope.HL))-1), \
					np.subtract(1, np.exp(np.divide(dvList,(-self.Isp*Constants().g0)))))

		self.a.plot(dvList, mRI_launch_List)
		xlabel = r'$\Delta$'+ 'v (m/s)'
		ylabel = self.radioisotope.name + " launch mass (g)"
		self.a.set_xlabel(xlabel)
		self.a.set_ylabel(ylabel)
		self.a.set_title(ylabel + " vs " + xlabel)
		self.a.grid()
		self.canvas.draw()

	def update(self):
		#-----------------clear plot-----------------#
		self.a.cla()

		#-----------------updates GUI with newly entered values-----------------#
		self.S_m0 = float(self.spacecraft_mass_entry.get())
		# self.RI_mT = float(self.RI_launch_mass_entry.get())
		self.Isp = float(self.Isp_entry.get())
		# self.zeta = float(self.zeta_entry.get())
		self.dv = float(self.dv_entry.get())
		self.zeta = float(1-np.exp(-self.dv/(self.Isp*Constants().g0)))
		self.burn_duration = np.arange(0,float(self.burn_duration_entry.get())*3.1536e7,86400) # **years to seconds
		self.transit_time = float(self.transit_time_entry.get())*3.1536e7 # **years to seconds		
		#-----------------updates GUI with new selections from dropdowns-----------------#
		self.radioisotope = self.radioisotopeDict.get(self.radioisotopeLabel.get())
		self.propellant = self.propellantDict.get(self.propellantLabel.get())
		self.ydata = self.ydataLabel.get()
		#-----------------get ionization efficiency of RI/prop combo-----------------#
		self.IE = DataProcessor().process3([self.radioisotope, self.propellant])
		#-----------------update in-app text-----------------#
		#-----------------update m_RI_T based on selected dv-----------------#
		self.RI_mT = -(np.exp(0.693*self.transit_time/self.radioisotope.HL)*Constants().N_A/(self.propellant.M*self.IE*self.radioisotope.SA*self.radioisotope.HL)* \
					self.S_m0/(np.exp(-0.693*(max(self.burn_duration)/self.radioisotope.HL))-1)*(1 - np.exp(self.dv/(-self.Isp*Constants().g0))))
		#-----------------update m_RI_T text-----------------#
		self.RI_mT_Label = Label(self, text=str(self.RI_mT))
		self.RI_mT_Label.grid(row=3, column=3, sticky=W)

		#-----------------run calc process with new data-----------------#
		self.plotdata = DataProcessor().process2([self.radioisotope, self.propellant],
									self.S_m0,
									self.RI_mT,
									self.Isp,	
									self.zeta,
									self.burn_duration,
									self.ydata,
									self.transit_time)




	def __init__(self, master=None):
		Frame.__init__(self, master)
		self.init()
		self.grid()

app = App()
app.mainloop()
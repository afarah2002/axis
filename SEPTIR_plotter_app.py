#-----------------imports-----------------#
import numpy as np
from tkinter import *
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.backends.backend_tkagg
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from scipy.optimize import curve_fit

#-----------------app, main class-----------------#
from SEPTIR_plotter import MaterialAssignment, DataProcessor, Constants, DataScanner


class App(Frame):

	def init(self):
		#-----------------title-----------------#
		self.font = ("Times New Roman", 20)
		title = Label(self, text="FIRESTORM: SEPTIR Performance Plotter",font=self.font)
		title.grid(row=1, column=1, sticky=W)

		#-----------------set up plot-----------------#
		self.f = Figure(figsize=(5,5), dpi= 100)
		self.a = self.f.add_subplot(111)
		# self.f.tight_layout(pad=1.08)
		self.a.grid()
		#-----------------add plot-----------------#
		self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.f, self)
		self.canvas.draw()
		self.canvas.get_tk_widget().grid(row=2, column=1)

		#-----------------create RI and propellant dropdown-----------------#
		#-----------------get list of RI and prop names for labelling-----------------#
		radioisotopesList, propellantsList = MaterialAssignment().assign()
		#-----------------drop downs-----------------#
		#-----------------RI drop down-----------------#

		#-----------------set RI and propellant attributes-----------------#
		radioisotopeNames = [RI.name for RI in radioisotopesList]
		self.radioisotopeDict = dict(zip(radioisotopeNames, radioisotopesList))

		self.radioisotopeLabel = StringVar(self)
		self.radioisotopeLabel.set(radioisotopeNames[0]) 
		RI_dropdown = OptionMenu(self, self.radioisotopeLabel, *radioisotopeNames)
		RI_dropdown.config(width=10, font=self.font)
		RI_dropdown.grid(row=4, column=1, sticky=W)

		#-----------------Propellant drop down-----------------#
		propellantNames = [prop.name for prop in propellantsList]
		self.propellantDict = dict(zip(propellantNames, propellantsList))

		self.propellantLabel = StringVar(self)
		self.propellantLabel.set(propellantNames[0])
		prop_dropdown = OptionMenu(self, self.propellantLabel, *propellantNames)
		prop_dropdown.config(width=10, font=self.font)
		prop_dropdown.grid(row=5, column=1, sticky=W)

		#-----------------YData drop down-----------------#
		#-----------------mdot or Fthrust-----------------#
		self.ydataLabel = StringVar(self)
		self.ydataLabel.set("")
		ydata_dropdown = OptionMenu(self, self.ydataLabel, *["mdot", "Fthrust"])
		ydata_dropdown.config(width=10, font=self.font)
		ydata_dropdown.grid(row=6, column=1, sticky=W)



		#-----------------Button: Plots mdot or thrust----------------#
		mdot_thrust_button = Button(self, text="Plot mdot/thrust", command=self.mdot_thrust, font=self.font)
		mdot_thrust_button.grid(row=8, column=1)

		#-----------------Button: mRI vs dv-----------------#
		plot_button = Button(self, text="Plot mRI vs dV", command=self.mri_vs_dv, font=self.font)
		plot_button.grid(row=9, column=1)

		#-----------------Button: Plots SR vs rho w/ fit-----------------#
		plot_button = Button(self, text="Plot SR vs density", command=self.SR_fitter, font=self.font)
		plot_button.grid(row=10, column=1)
		#-------------------------------------------------------------#

		#-----------------entry boxes-----------------#
		#-----------------spacecraft mass, g-----------------#
		spacecraft_mass_entry_Label = Label(self, text="Spacecraft mass (g)", font=self.font)
		spacecraft_mass_entry_Label.grid(row=4, column=2, sticky=W)
		self.spacecraft_mass_entry = Entry(self)
		self.spacecraft_mass_entry.grid(row=4, column=3, sticky=W)
		#-----------------RI initial mass text box, g-----------------#
		RI_launch_mass_text_Label = Label(self, text="RI launch mass (g)", font=self.font)
		RI_launch_mass_text_Label.grid(row=9, column=2, sticky=W)
		#-----------------propellant mass text box, g-----------------#
		propellant_mass_text_Label = Label(self, text="Total propellant mass (g)", font=self.font)
		propellant_mass_text_Label.grid(row=10, column=2, sticky=W)
		#-----------------Isp, s-----------------#
		Isp_entry_Label = Label(self, text="Isp (s)", font=self.font)
		Isp_entry_Label.grid(row=5, column=2, sticky=W)
		self.Isp_entry = Entry(self)
		self.Isp_entry.grid(row=5, column=3, sticky=W)
		#-----------------delta v (m/s)-----------------#
		dv_entry_Label = Label(self, text="dv (m/s)", font=self.font)
		dv_entry_Label.grid(row=6, column=2, sticky=W)
		self.dv_entry = Entry(self)
		self.dv_entry.grid(row=6, column=3, sticky=W)
		#-----------------burn duration, years-----------------#
		burn_duration_entry_Label = Label(self, text="Burn duration (years)", font=self.font)
		burn_duration_entry_Label.grid(row=7, column=2, sticky=W)
		self.burn_duration_entry = Entry(self)
		self.burn_duration_entry.grid(row=7, column=3, sticky=W)
		#-----------------transit time, years-----------------#
		transit_time_entry_Label = Label(self, text="Transit time (years)", font=self.font)
		transit_time_entry_Label.grid(row=8, column=2, sticky=W)
		self.transit_time_entry = Entry(self)
		self.transit_time_entry.grid(row=8, column=3, sticky=W)
		#---------------------------------------------#


	def mdot_thrust(self):
		#-----------------processes a specific RI-propellant combination-----------------#
		self.update()
		self.a.plot(np.divide(self.burn_duration,3.1536e7), np.multiply(self.plotdata,1000))
		#-----------------plot formatting-----------------#
		xlabel = "Burn time (years)"
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

	#-----------------FITTING CURVES-----------------#
	def rational(self, x, p, q):
		return np.polyval(p, x) / np.polyval(q + [1.0], x)

	def rational3_3(self, x, p0, p1, p2, q1, q2):
		return self.rational(x, [p0, p1, p2], [q1, q2])

	def SR_fitter(self):
		#-----------------plots a fit of the SR vs prop rho curve-----------------#
		#-----------------plot raw data-----------------#
		self.update()
		rho_list = np.arange(0.006,.05, .0001)
		data_file = "data/" + self.propellant.name + "2/"+str(self.radioisotope.energy)+"_SR.txt"
		SR_list = DataScanner().read_data(data_file)
		self.a.plot(rho_list, SR_list)
		#-----------------plot fit curve-----------------#
		# popt, pcov = curve_fit(self.rational3_3, rho_list, SR_list, p0=(0.2, 0.3, 0.5, -1.0, 2.0), maxfev=5000)
		# self.a.plot(rho_list, rational3_3(rho_list, *popt), 'r-', label='fit')

		xlabel = r'$\rho$'+ ' (g/cm' + r'$^{3}$' + ')'
		ylabel = r'$\alpha^{+}$' + ' Stopping Range (mm)'
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
		self.Isp = float(self.Isp_entry.get())
		self.dv = float(self.dv_entry.get())
		self.zeta = float(1-np.exp(-self.dv/(self.Isp*Constants().g0)))
		self.burn_duration = np.arange(0,float(self.burn_duration_entry.get())*3.1536e7,86400.) # **years to seconds
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
					self.S_m0/(np.exp(-0.693*(float(max(self.burn_duration))/self.radioisotope.HL))-1)*(1 - np.exp(self.dv/(-self.Isp*Constants().g0))))
		self.RI_mT_Label = Label(self, text=str(round(self.RI_mT, 3)), font=self.font)
		self.RI_mT_Label.grid(row=9, column=3, sticky=W)
		#-----------------update propellant mass text based on dv-----------------#
		self.m_prop = self.S_m0*(1 - np.exp(-self.dv/(self.Isp*Constants().g0)))
		self.propellant_mass_Label = Label(self, text=str(round(self.m_prop,3)), font=self.font)
		self.propellant_mass_Label.grid(row=10, column=3, sticky=W)
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
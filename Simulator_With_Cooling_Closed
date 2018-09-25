#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 08:43:27 2018

@author: test
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 16:50:05 2018

@author: test
"""
'''
This code is a gekko simulation create to model the dynamic temperature of the 
Aquila HALE UAV during a flight (data provided below in excel file). This is a 
collaboration of Christian Abbott - christian.t.abbott@gmail.com and Nathaniel
Gates -  nsgates@gmail.com.

SimWithCooling Version 1.0
Last updated: 9/7/2018

MAJOR ASSUMPTIONS IN MODEL
- The aircraft in question is similar to the Aquila HALE UAV
- The battery network behaves comparable to a single large cylinderic battery
- The ambient air and air let into the battery pod is completely turbulent and 
    the Nu number is constant accross the battery pod
- The battery resistance is a funtion of battery temperature and follows an
    Arrhenius-based expression and is accurately approximated with the average
    battery temperature even though the battery experiences a temperature 
    gradient.
- The battery temperature only changes axially, not radially
- The density of the cooling air is constant axially throughout the annulus. And
    is calculated assuming the ambient air conditions (Will provide inflated cooling)
'''

#import pip
#pip.main(['install','gekko'])
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from gekko import GEKKO

#Data Import
excel_file = 'Facebook_data2.xlsx'
data = pd.read_excel(excel_file)
v_data = data['v'] #m/s
rho_data = data['rho'] #kg/m3
t_air_data = data['t_air'] #K
re_data = data['re']
time_data = data['time'] #seconds
motor_data = data['p_n'] #W 
h_data = data['h'] #height (m)
charge_data = data['e_batt'] #MJ
mu_data = data['mu']
motor_efficiency_data = data['eta_motor']
motor_voltage_data = data['motor_voltage'] #V
motor_current_data = data['motor_current'] #A


#%% Compiling Data

####################### CHANGE ME #############################################

#Selecting Data to use in Gekko model from original 
n = len(data)-1  #Number of data points to calculate.
shift = 1  #Shift data that is selected by this many (max data point = 8635)
Optimal_temp = 298 #K (Assumed optimal battery temperature)
space = 4 #only include every (space)th data point
data_pieces = int(n/space)-2

###############################################################################

#Storage for loop below
v = np.zeros(data_pieces) #velocity// storage
rho = np.zeros(data_pieces) #air density storage
t_air = np.zeros(data_pieces) #ambient temperature storage
re2 = np.zeros(data_pieces) #reynolds number storage
time = np.zeros(data_pieces) #time storage 
motor = np.zeros(data_pieces) #motor power storage
motor[0] = motor_data[0+shift] #motor power initial value
T_SP = np.ones(data_pieces)*Optimal_temp #Creating an array of optimal temperatures (K)
mu = np.zeros(data_pieces) #viscosity storage
motor_efficiency = np.zeros(data_pieces)
motor_voltage = np.zeros(data_pieces)
motor_current = np.zeros(data_pieces)


#Loop to select data from large data set
for d in range(data_pieces):
    v[d] = float(v_data[d*space+shift])
    rho[d] = float(rho_data[d*space+shift])
    t_air[d] = float(t_air_data[d*space+shift])
    re2[d] = float(re_data[d*space+shift])
    time[d] = float(time_data[d*space+shift]-time_data[0+shift]) #makes the first set of data start at t=0
    motor[d] = float(motor_data[d*space+shift])
    mu[d] = float(mu_data[d*space+shift])
    motor_efficiency[d] = float(motor_efficiency_data[d*space+shift])
    motor_voltage[d] = float(motor_voltage_data[d*space+shift])
    motor_current[d] = float(motor_current_data[d*space+shift])
    


#%% Gekko model

m = GEKKO(server='http://byu.apmonitor.com')
m.time = time 

#Material property constants and assumptions about the batter/motor pod shape and size (in gekko objects)
Optimal_Temp = m.Const(value = Optimal_temp,name='Optimal_Temp')
MW = m.Const(value=0.0289644,name='MW') #kg/mol
Cp_air = m.Const(value=975.4,name='Cp_air') #J/kgK (found using HYSYS calculations)
kf_air = m.Const(value=.01855,name='kf_air') #W/mK (found using HYSYS calculations)
mass_batt = m.Const(value=142/4,name='mass_batt') #kg
Cp_batt = m.Const(value=900,name='Cp_batt') #J/kgK
kf_batt = m.Const(value=6,name='kf_batt') #W/mK
radius_batt = m.Const(value=.5/2,name='radius_batt') #m (Estimation based on photographs)
length_total = m.Const(value=2,name='length_total') #m (Estimation based on photographs)
length_motor = m.Const(value=.5,name='length_motor') #m (Guess)
V_open_circuit = m.Const(value=44,name='V_open_circuit') #Volts
mu2 = m.Param(value=mu,name='mu2')

#Intermediates
length_batt = m.Intermediate(length_total - length_motor,name='length_batt') #m
Pr = m.Intermediate(Cp_air*mu2/kf_air,name='Pr') #Prandlt number assumed to be the same for both internal and external air

#Parameters necessary for the model made into gekko model objects
rho_air - m.Param(value=rho,name='rho_air') #density of air (changes due to temperature and pressure variations during flight)
motor_power = m.Param(value=motor,name='motor_power') #W (Not used)
t_air2 = m.Param(value=t_air,name='t_air2') #K
T_SP2 = m.Param(value=T_SP,name='T_SP2') #K
re_external = m.Param(value=re2,name='re_external')
re_internal = m.Param(value=re2,name='re_internal')
eff_motor = m.Param(value=motor_efficiency,name='eff_motor') #efficiency of the motor
I_load = m.Param(value=motor_current,name='I_load') #Amps used by motor
V_load = m.Param(value=motor_voltage,name='V_load') #Volts dropped by motor useage 

############ CHANGE ME ###############
# How many slices to break up the battery into (length-wise)
discretize = 11 #Number of discretizations spatially
thick_insulation = m.Const(value=.0207,name='thick_insulation') #m
#.008037601 optimized thickness with constant internal battery resistance
#.02389585 optimized thickness with variable internal battery resistance
######################################

# More constants and assumptions about the model made into gekko objects 
thick_air = m.Const(value=.02,name='thick_air') #m      #Thickness of an inside tube of air used for insulation and cooling (tube can be opened to the ambient air)
thick_motor_insulation = m.Const(value=0,name='thick_motor_insulation') #m
thick_carbon = m.Const(value=.0004*2,name='thick_carbon')        #m
num_slices = m.Const(value=discretize,name='num_slices') #Inputting the number of slices into Gekko model
kf_insulation = m.Const(value=.03,name='kf_insulation') #W/mK
kf_carbon = m.Const(value=100.5,name='kf_carbon') #W/mK
Cp_motor = m.Const(value=900,name='Cp_motor') #J/kgK
Cp_insulation = m.Const(value=1331,name='Cp_insulation') #J/kgK Calculated at 298K
Cp_carbon = m.Const(value=1130,name='Cp_carbon') #J/kgK
dens_insulation = m.Const(value=.205,name='dens_insulation') #kg/m3
dens_carbon = m.Const(value=1.800,name='dens_carbon') #kg/m3
slice_length = m.Intermediate(length_batt/num_slices,name='slice_length') #Length of slice
Nu_external = m.Intermediate( .0296*re_external**(4/5)*Pr**(1/3) ,name='Nu_external') #Nusselt number correlation for completely turbulent air flow at length = length of the whole battery pod
Nu_internal = m.Intermediate( .0296*re_internal**(4/5)*Pr**(1/3) ,name='Nu_internal') #Nusselt number correlation for completely turbulent air flow at length = length of the whole battery pod
D_hydraulic = m.Const(value=thick_air,name='D_hydraulic') #hydraulic diameter used to calculate convective heat transfer coefficient inside an annulus
h_external = m.Intermediate( Nu_external*kf_air/length_batt,name='h_external') #convective heat transfer coefficient for external air
h_internal = m.Intermediate( Nu_internal*kf_air/length_batt,name='h_internal') #convective heat transfer coefficient for air passing through annulus

#Geometries necessary to calculate the heat transfer through discretized segments
D1 = m.Intermediate( 2*radius_batt,name='D1') #Diameter of the battery/motor
D2 = m.Intermediate( 2*thick_air+D1,name='D2') #Outer diameter for the air
D3 = m.Intermediate( 2*thick_insulation+D2,name='D3') #Outer diameter of the insulation
D4 = m.Intermediate( 2*thick_carbon+D3,name='D4') #Outer diameter of the carbon fiber
motor_face = m.Intermediate( np.pi*D1*length_motor+np.pi*radius_batt**2 ,name='A1_inner') #Surface area of motor casing touching air
A1_face = m.Intermediate( np.pi*(D1/2)**2 ,name='A1_face') #area of the battery face
A2_face = m.Intermediate( np.pi*(D2/2)**2 - A1_face,name='A2_face') #area of the air annulus
A3_face = m.Intermediate( np.pi*(D3/2)**2 - A2_face,name='A3_face') #area of the insulation annulus
A4_face = m.Intermediate( np.pi*(D4/2)**2 - A3_face,name='A4_face') #area of the carbon fiber annulus
A1_shell = m.Intermediate( np.pi*D1*slice_length,name='A1_shell') #area of battery shell
A2_shell = m.Intermediate( np.pi*D2*slice_length,name='A2_shell') #area of air shell
A3_shell = m.Intermediate( np.pi*D3*slice_length,name='A3_shell') #area of insulation shell
A4_shell = m.Intermediate( np.pi*D4*slice_length,name='A4_shell') #area of carbon fiber shell (very outside)
#Important intermediate that describes the resistance of heat loss from battery to the shell
U_slice = m.Intermediate( 1/(1/(h_internal*A1_shell)+1/(h_internal*A2_shell)+m.log(D3/D2)/(2*np.pi*kf_insulation*length_batt/num_slices)+m.log(D4/D3)/(2*np.pi*kf_carbon*length_batt/num_slices)+1/(h_external*A4_shell)))

vol_annulus = m.Intermediate((A2_face-A1_face)*slice_length)
mass_air_annulus = m.Intermediate( rho_air*vol_annulus)
############################# Equation Section ################################

"""
Calculating an internal battery resistance that varies with temperature is important for this
model because it provides the optimizer a push away from lower temperatures (where energy is greatly lost)
All constants used to calculate the battery internal resistance (R_ref, T_ref, and beta)
were taken from the following literature:
Motapon, Souleman Njoya, et al. "A Generic electrothermal li-ion battery model for 
rapid evaluation of cell temperature temporal evolution." IEEE Transactions on Industrial
 Electronics 64.2 (2017): 998-1008.
"""
R_ref = m.Const(value=0.0126,name='R_ref')
T_ref = m.Const(value=299.15,name='T_ref')
beta = m.Const(value=2.836*10**3,name='beta')
air_valve = m.Const(value=,lb=0,ub=1,name='air_valve') #IMPORTANT: This is the actuator that lets in air to cool the battery (CV -> simulation, MV -> controller)

R_batt = m.SV(name='R_batt')
Heater_Energy = m.SV(name='Heater_Energy') #calculating total heat over the time period devoted to heating the battery (assuming all heat from heaters goes to battery)
T_motor = m.SV(value=273+10,name='T_motor') # Setting the motor temperature as a state variable (not to be controlled at optimal temperature)
T = [m.SV(value=273+20,name='Tbatt_'+str(i)) for i in range(discretize)] #Create a battery temperature variable for (discretize) number of sections 
Ta = [m.SV(value=273+20,name='Tair_'+str(i)) for i in range(discretize)] #Create an air temperature variable for (discretize) number of sections 
Ta[0] = m.MV(lb=215,ub=300)
#Q_heater = [m.MV(value=0,name='Q_'+str(i)) for i in range(discretize)] #Creating a heater variable for each slice
Q_heater = m.Const(value=0,name='Q_heater')
T_avg = m.SV(value=273+20,name='T_avg') #Average battery temperature (needed to use a 1 heater temperature management system) -- also good assumption

m.Equation( T_avg == sum(T)/discretize)
m.Equation( R_batt == R_ref*m.exp(beta*(1/T_avg-1/T_ref)))
m.Equation( Heater_Energy.dt() == Q_heater) 
for j in range(discretize):
    if j == 0:
        # Motor (Section 0)
        m.Equation( (mass_batt/4*Cp_motor)*T_motor.dt() ==  motor_power*(1-eff_motor) #Heat generation (From Motor) 
                                                                 -kf_batt*A2_inner*(T_motor-T[j])  #Conduction to next slice
                                                                 -h*A1_inner*(T_motor-t_air2)) #Convection from motor 
#                                                                 -U_motor*(T[j]-t_air2))       #Convection from slice through carbon 
        if open_air = True:
            Ta[j] = 215
        
        # Calculating temperature for slice touching motor
        m.Equation( (mass_batt*Cp_batt)/num_slices*T[j].dt() == kf_batt*A2_inner*(T_motor-T[j]) #Conduction from previous slice
                                                            + I_load**2*R_batt/num_slices #Heat generation (Joule heating)
                                                            -kf_batt*A2_inner*(T[j]-T[j+1]) #Conduction to next slice
                                                            -U_slice*(T[j]-t_air2) #Conduction through both insulation and carbon with convection on outside
                                                            +Q_heater/num_slices) #Heat from external coil
        
        m.Equation( (mass_air_annulus*Cp_air)/num_slices*Ta[j].dt() == )                              

                               
    elif j > 0 and j < discretize-1:
        # Battery slice (Section 1 --> Section n-1)
        m.Equation( (mass_batt*Cp_batt)/num_slices*T[j].dt() == kf_batt*A2_inner*(T[j-1]-T[j]) #Conduction from previous slice
                                                                + I_load**2*R_batt/num_slices #Heat generation (Joule heating)
                                                                -kf_batt*A2_inner*(T[j]-T[j+1]) #Conduction to next slice
                                                                -U_slice*(T[j]-t_air2) #Conduction through both insulation and carbon with convection on outside
                                                                +Q_heater/num_slices) #Heat from external coil 
    else:
        # Battery (Last section)
        m.Equation( (mass_batt*Cp_batt)/num_slices*T[j].dt() == kf_batt*A2_inner*(T[j-1]-T[j]) #Conduction from previous slice
                                                                + I_load**2*R_batt/num_slices #Heat generation (Joule heating)
                                                                -U_slice*(T[j]-t_air2) #Conduction through both insulation and carbon with convection on outside
                                                                +Q_heater/num_slices) #Heat from external coil
        
############################ Global Options ###################################


m.options.imode = 4
m.options.nodes = 3
m.options.solver = 3
#m.options.EV_Type = 1

m.solve() #Solve simulation

#Load results
sol = m.load_results()
solData = pd.DataFrame.from_dict(sol)

#%% Plot results
plt.figure(figsize=(7,7))
plt.title('Battery Temperature Profile '+str(discretize)+' slices '+str(data_pieces)+' points')
plt.ylabel('Temperature (C)')
plt.xlabel('Time (hr)')
plt.plot(solData.time/3600,solData['t_motor']-273.15,label='T_motor')
for i in range(discretize):
    plt.plot(solData.time/3600, solData['t_'+str(i)]-273.15,label='T_'+str(i))
if discretize < 20:
    plt.legend(["T_"+str(q) for q in range(discretize)])
        
    #Save results
    #m.CSV_WRITE = 1
#    plt.savefig('/Users/test/Desktop/School/Winter 2018/Dynamic_Control/FACEBOOK_Project/Grid_Refinement/'+str(discretize)+' slices '+str(data_pieces)+' points')
#    plt.savefig('/Users/test/Desktop/Work/FACEBOOK/Added_Resistance'+str(thicknesses[step])+' m insulation'+' with resistance.png')

"""
Created on Sat Nov 2018
Example of the NTempObject version
We consider the 2 layers SRO=Strontium Oxide and STO=Strontium Titanium Oxide
with respective parameters for length, k(T), C(T), rho (and coupling G)
We look at two cases: 
    1) only electron temperature
    2) electron and lattice temperature
Dependencies: Numpy, Matplotlib, B-Spline Package, (SI- units package)
@author: Lukas Alber & Valentino Scalera
"""
from NTempObjectGIT import * 
import numericalunits as u
u.reset_units('SI')


#%% Define a Source
s                           = source()
s.optical_penetration_depth = [45*u.nm,50*u.nm]  #optical penetration depth for each layer 
s.FWHM                      = 1*u.ps 
s.fluence                   = 15*u.mJ/u.cm**2 
s.t0                        = 5*u.ps            #Peake of Gaussian
#%%1 Temperature model: Set up simulation 
#Define numer of Temperatures taken into consideration 1 => only electron temperature
#Pass on the source object s as an input parameter
sim = simulation(1,s) 
#add layers (Length,conductivity,heatCapacity,density,coupling)
sim.addLayer(40e-9,[6],[lambda Te: 0.112*Te],6500) #SRO Layer
sim.addLayer(80e-9,[12],[lambda Te: 0.025*Te],5100)#STO Layer
sim.final_time = 50*u.ps                           #final time of the simulation
#To get the raw output
[phi_E,x,t] = sim.run() 
#%%2 Temperature Model: 
#Two temperatures are considered, electron and lattice
sim = simulation(2,s)
#add parameters for both layers and both systems
sim.addLayer(40e-9,[6,1],[lambda Te: 0.112*Te,450],6500,5e17) #SRO Layer
sim.addLayer(80e-9,[12,1],[lambda Te: 0.025*Te,730],5100,5e17)#STO Layer
sim.final_time = 50*u.ps
#To get the raw output
[phi_E,phi_L,x,t] = sim.run() 
#%%Visualize result
#Create a visual object where the simulation gets passed on 
v = visual(sim)
#output of v.source is the full matrix of the source(x,t)
so = v.source()
#A contourplot of the electron system
v.contour('Electron')
v.average()
#This will only work in the two temperature case.
v.contour('Lattice') 

#%%Animation of dynamics
#animate the result input parameter is the speed of the animation
v.animation(1)





# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 2018
This is our code we use for the 2Temperatue solver. 
1)Multi layers and 2 temperatures can be solved. 
2)The timestep for the euler loop is automatically calculated.
3)We use an approch of B-splines and explicit euler to compute the solution
4)In the source class one can find two different approaches on how the energy
deposit is being calculated.
    a) Transfer matrix method (Reflection, incident angle, refractive index, polarization)
    b)  Lambert Beer law (exponential decay in space according to optical penetration depth)
@author: Lukas Alber & Valentino Scalera
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from bspline import Bspline
from bspline.splinelab import aptknt
import time
from matplotlib.animation import FuncAnimation as movie

from tqdm import tqdm


#==============================================================================
class temperature(object): 
        
    def __init__(self): 
        self.plt_points     = 30                    #number of points in x grid
        self.length         = np.array([0,0])       #length of x space,starting from 0
        self.Left_BC_Type   = 1                     #Boundary conditions Default is Neumann
        self.Right_BC_Type  = 1                     #1=> Neumann; 0=> Dirichlet
        self.init           = lambda x: 300+0*x     # initial temperature of probe
        self.conductivity   =  [1]                  #This gets deleted after initialisation
        self.heatCapacity   =  [1]                  #those values are just here to make space
        self.rho            = [1]                   #Actual values are given, when 'addLayer(length, conductivity,heatCapacity,rho)' is executed
        self.setup          = False                 #first time setup to not double calculated
    
    def getProperties(self):                        #to depict the properties of the object
        for i in (self.__dict__): 
            print(i,' : ',self.__dict__[i])
            
    def __repr__(self): 
        return('Temperature')
    #for every layer, a function to calculate the derivative of k(T)        
    def diff_conductivity(self,phi,num_of_material): 
        eps =1e-9
        dc = (self.conductivity[num_of_material](phi+eps)-self.conductivity[num_of_material](phi))/eps
        return(dc)
    #Creating the key matrices for B-splines. Those are A0,A1,A2
    #A0 => Zero derivative; A1 => 1st order derivative.... 
    #We create the matices for every layer, with respective length ect
    #then we put them together to Abig => Boundary and interface conditions are applied here.         
    def Msetup(self):
        #Deleting the ifrst element of the default initialization
        #After creating the element with 'addLayer' we dont need them!
        if not self.setup:
            self.length         = self.length[1:]
            self.conductivity   = self.conductivity[1:]
            self.heatCapacity   = self.heatCapacity[1:]
            self.rho            = self.rho[1:]
            self.setup          = True
        #Length and numper of grid points for each respective layer    
        length          = self.length
        plt_points      = self.plt_points
        num_of_points   = 12  #Number of points per layer used in the spline
        order           = 5   #order of the spline
        x               = np.array(np.zeros([np.size(length)-1,num_of_points]))
        x_plt           = np.array(np.zeros([np.size(length)-1,plt_points]))
        knot_vector     = np.array(np.zeros([np.size(length)-1,num_of_points+order+1]))
        basis           = np.array(np.zeros(np.size(length)-1))
        A0h = []; A1h = []; A2h = []; Ch = [];
        LayerMat        = np.array([np.zeros((num_of_points,num_of_points))])
        
        #Create all the big matices A0,A1,A2 & C. C is used to map on a fine mesh in x- space.
        #For every layer we set up splines between the boundaries
        for i in range(0,np.size(length)-1):
            x[i,:]   = np.linspace(length[i], length[i+1] , num_of_points)
            x_plt[i,:]       = np.linspace(length[i], length[i+1] , plt_points)
            knot_vector[i,:] = aptknt(x[i,:], order) #prepare for Spline matrix
            basis    = Bspline(knot_vector[i,:],order)
            A0hinter = basis.collmat(x[i,:], deriv_order = 0); A0hinter[-1,-1] = 1
            A1hinter = basis.collmat(x[i,:], deriv_order = 1); A1hinter[-1] = -np.flip(A1hinter[0],0)
            A2hinter = basis.collmat(x[i,:], deriv_order = 2); A2hinter[-1,-1] = 1
            Chinter  = basis.collmat(x_plt[i,:], deriv_order = 0); Chinter[-1,-1] = 1
            LayerMat = np.append(LayerMat,np.array([np.dot(A2hinter,np.linalg.inv(A0hinter))]),axis = 0)
            A0h      =  np.append(A0h,A0hinter)
            A1h      =  np.append(A1h,A1hinter)
            A2h      =  np.append(A2h,A2hinter)
            Ch       = np.append(Ch,Chinter)
        #Reshape the long string of appended Matrix, such that
        #rows: x-points; colums: i´th basis spline
        LayerMat = LayerMat[1:,:,:]
        A0h = np.reshape(A0h, (-1,num_of_points)) 
        A1h = np.reshape(A1h, (-1,num_of_points))
        A2h = np.reshape(A2h, (-1,num_of_points))
        Ch  = np.reshape(Ch,(-1,num_of_points)) 
        #Ch => More points in x, but same number of basis splines
        #Clearing the interface points, to not double count 
        N               = num_of_points 
        plp             = plt_points
        interfaces      = np.shape(x)[0]-1
        sizeA           = np.shape(x)[0]*N-interfaces
        sizeCb          = np.shape(x)[0]*plp-interfaces
        Abig            = np.zeros([sizeA,sizeA]) 
        A1b             = np.zeros([sizeA,sizeA])
        A2b             = np.zeros([sizeA,sizeA])
        Cb              = np.zeros([sizeCb,sizeA])
        #Clearing the double counts from the space grid 
        xflat = x.flatten()
        x_plt_flat = x_plt.flatten()
        #index of double counts        
        doublec = np.array([np.arange(1,len(length)-1)])*N
        doublec_plt = np.array([np.arange(1,len(length)-1)])*plp
        xflat = np.delete(xflat,doublec)
        x_plt_flat = np.delete(x_plt_flat,doublec_plt)
        #Filling the big matrices.
        startA = 0; endA = N-1
        startC = 0; endC = plp-1 
        for i in range(0,interfaces+1):
            Abig[startA:endA,startA:endA+1]  = A0h[startA+i:endA+i,:]   
            A1b[startA:endA+1,startA:endA+1] = A1h[startA+i:endA+i+1,:]
            A2b[startA:endA+1,startA:endA+1] = A2h[startA+i:endA+i+1,:]
            Cb[startC:endC+1,startA:endA+1]  = Ch[startC+i:endC+i+1,:] 
            startA += N-1; endA += N-1
            startC += plp-1; endC += plp-1
            
        #Create A00 with no interface condition to correctly compute phi in loop
        #The copy needs to be done befor interface conditions are applied in Abig
        A00 = Abig.copy() 
        A00[-1,-1] = 1;
        #Here we make init, conductivity & capacity all functions, in case they are
        # given as integeres or floats. Also thorw warinings if not every layer has a
        # conducitvity or capacity ============================================
        #Making init a function, in case it is given as a scalar
        if np.size(self.init) == 1 and isinstance(self.init,(int,float)):
            dummy   = self.init                                      
            self.init    = lambda x: dummy + 0*x        
        if len(length) > 2: #multilayer case
            if len(length)-1 !=(  len(self.heatCapacity) & len(self.conductivity) ): 
                print('--------------------------------------------------------')
                print('The number of different layers must match the number of number of' \
                      'inputs for Conductivity, heatCapacity, rho.')
                print('--------------------------------------------------------')
            if np.size(self.conductivity) is not interfaces+1:
                print('--------------------------------------------------------')
                print('Not every Layer has been given a conductivity function' \
                      'Adjust your input of the conductivity functions with respect to the layers.')
                print('--------------------------------------------------------')
            if np.size(self.heatCapacity) is not interfaces+1:
                print('--------------------------------------------------------')
                print('Not every Layer has been given a heatCapacity function value.'\
                      'Adjust your input of the heatCapacity functions with respect to the layers.')   
                print('--------------------------------------------------------')
        #Make Functions in case heat capacity/conductivity are given as variables
        if (all(self.conductivity) or all(self.heatCapacity) or all(self.init)) == False:
            print('No heatCapacity, conductivity or initial function given.')
            print('--------------------------------------------------------')
            #make the conductivity always a function
        if len(length) >2 or np.size(self.conductivity)>=2:
            for j in list(range (0,np.size(self.conductivity))):
                if isinstance(self.conductivity[j],(int,float,list)) :  
                    dummy3 = self.conductivity[j]
                    self.conductivity[j] = (lambda b: lambda a: b+0*a)(dummy3) 
            #make the conductivity always a function 
            for j in list(range (0,np.size(self.heatCapacity))):  
                if isinstance(self.heatCapacity[j],(int, float,list)) : 
                    dummy4 = self.heatCapacity[j]
                    self.heatCapacity[j] = (lambda b: lambda a: b+0*a)(dummy4)  
        else : 
            if isinstance(self.conductivity[0],(int,float)):
                dummy1 = self.conductivity
                self.conductivity = [lambda phi: dummy1 + 0*phi]
               
            if isinstance(self.heatCapacity[0],(int,float)):
                dummy2 = self.heatCapacity
                self.heatCapacity = lambda phi: dummy2 + 0*phi
                self.heatCapacity = [self.heatCapacity]          
        #End of function creation for init(x), conductivity[l](phi), heatCapacity[l](phi)
        # with respect to every layer 'l' =====================================
        def interconditions(phi,interfaces):
            N = num_of_points
            end_i = N-1
            intercondiL = np.zeros((interfaces,N)) 
            intercondiR = np.zeros((interfaces,N)) 
            for i in range(interfaces): 
                intercondiL[i] = self.conductivity[i](phi[end_i])*A1h[end_i+i]
                intercondiR[i] = self.conductivity[i+1](phi[end_i])*A1h[end_i+i+1]
                end_i += N-1
            return(intercondiL,intercondiR) 
        #Initial Electron temperature   
        initphi = self.init(xflat)       
        initphi_large = self.init(x_plt_flat)
        intercon = interconditions(initphi,interfaces)            
        #filling up Abig wiht the interface condition in the middle of the grid
        start_i = 0; end_i = N-1  
        for i in range(0,interfaces): 
            Abig[end_i,start_i:end_i]    =  intercon[0][i][:-1]#Lhs interface flow
            Abig[end_i,end_i+1:end_i+N]  = -intercon[1][i][1:]#Rhs interface flow
            Abig[end_i,end_i]            = intercon[0][i][-1] -intercon[1][i][0]
            start_i += N-1; end_i += N-1
        Abig[-1,-1] = 1 #to correct Cox algorithm        
        #Now Matrix Abig is completed and interface condition is applied.
        #Treating 2 types of boundary conditions: 0=> Dirichlet; 1=> Neumann, 
        # where 0´th and -1´th row need to be first order derivatives for flux.
        neumannBL = A1b[0].copy(); 
        neumannBR = A1b[-1].copy(); 
        if self.Left_BC_Type  == 1: Abig[0]  = -neumannBL
        if self.Right_BC_Type == 1: Abig[-1] = neumannBR    
        #Clear for BC! (first and last row need to be cleared to correctly apply BC)         
        A1b[0]     = 0; A2b[0]     = 0;
        A1b[-1]    = 0; A2b[-1]    = 0;
        #Get inital c coefficients for splines using init (=phi_init)
        c = np.dot(np.linalg.inv(A00),self.init(xflat))
        #Passed on properties to the simulation class
        return(c,A00,Abig,A1b,A2b,Cb,length,N,plp,xflat,x_plt_flat,initphi_large,interfaces,LayerMat,A1h)
        
    def addLayer(self,L,conductivity,heatCapacity,rho):
        self.length = np.append(self.length,self.length[-1]+L)
        self.conductivity.append(conductivity)
        self.heatCapacity.append(heatCapacity)
        self.rho = np.append(self.rho,rho)
        
#==============================================================================        
class simulation(object): 
    
    def __init__(self,num_of_temp,source): 
        self.temp_data          = temperature() #import the temperatuer object
        self.num_of_temp        = num_of_temp   #1 if only electron temp. 2 if electron and lattice temp.
        self.start_time         = 0             #starting time (can be negative)
        self.final_time         = 10            #time when simulation stops
        self.time_step          = []            #can either be given or is automatically calculated in stability
        self.left_BC            = 0             #function or constant what the boundary condition 
        self.right_BC           = 0             #on the left or right side of the problem is.
        self.temp_data_Lat      = []            #Default case is without lattice temperature
        if num_of_temp == 2:                    #if Lattice temp is considered
            self.temp_data_Lat  = temperature() #in case also a lattice module is given
            self.coupling       = []            #Coupling between the two systems
            self.left_BC_L      = 0
            self.right_BC_L     = 0
        self.source             = source        #object source can be passed on
    #to depict the properties of the object   
    def getProperties(self): 
        for i in (self.__dict__): 
            print(i,' : ',self.__dict__[i])
            
    def __repr__(self): 
        return('Simulation')
        
    def addLayer(self,L,conductivity,heatCapacity,rho,coupling=0):
        #check all input arguments and make them to lists, for the multi layer case
        #make list when given as int or float
        typecheck = np.array([])
        if type(conductivity) is not (list or type(typecheck)):
            conductivity = [conductivity]
        if type(heatCapacity) is not (list or type(typecheck)):
            heatCapacity = [heatCapacity]
        #do typecheck for the lattice system
        if self.num_of_temp == 2:
            if (np.size(conductivity) or np.size(heatCapacity))<2: 
                print('Lattice parameters are missing.\n Add parameters for Lattice system.')
                return(128)
            self.temp_data_Lat.addLayer(L,conductivity[1],heatCapacity[1],rho)
            self.coupling = np.append(self.coupling,coupling)   
        self.temp_data.addLayer(L,conductivity[0],heatCapacity[0],rho)
    #a function which gives back an array where the intereface condition is returned
    #for the left and right side of the interface. Gets called in the loop.     
    def interconditions(self,phi,interfaces,conductivity,N,A1h):
        end_i = N-1
        intercondiL = np.zeros((interfaces,N)) 
        intercondiR = np.zeros((interfaces,N)) 
        for i in range(interfaces): 
            intercondiL[i] = conductivity[i](phi[end_i])*A1h[end_i+i] 
            intercondiR[i] = conductivity[i+1](phi[end_i])*A1h[end_i+i+1]
            end_i += N-1        
        return(intercondiL,intercondiR) 
    #create mulit layer source
    def init_G_source(self,xflat,t,opt_pen,N,func):
        """
        First an empty array 'sourceM' is created.
        Then we iterate over the different layers and call
        func --> Gaussian. 
        This will create a 2D (time, space) Gaussian source grid
        with different lam[i].
        For each layer, the problem is receted, i.e. we have new
        Amplitude, new scope of the x-grid, new lambda. Only sigma stays the same.
        """
        x0 = self.temp_data.length#interfaces

        lam = 1/opt_pen 
        #create space for solution of source in matrix form
        xmg, tmg = np.meshgrid(xflat,t)
        sourceM  = np.zeros(np.shape(xmg))
        #Convert the input, fluence & FWHM given in 'source' class to Amplitude and sigma
        sigma2 = self.source.FWHM**2/(2*np.log(2))
        A = self.source.fluence/np.sqrt(2*np.pi*sigma2)
        #loop over all layers and change lambda, Amplitude and scope of the x-grid
        startL = 0; endL = N-1
        for i in range(2,len(opt_pen)+2):
            sourceM[:,startL:endL] = func(xmg[:,startL:endL],tmg[:,startL:endL],lam[i-2],A,sigma2,x0[i-2])            
            #at the end of each loop: the intensity of the end of previous layer
            #will be the starting intensity of the next layer
            A = A*np.exp(-(x0[i-1]-x0[i-2])*lam[i-2])
            startL = endL; endL = i*N-i+1
        return(sourceM)

    # This is the main Explicit Euler loop where the solution to T(x,t) is calculated.
    def run(self):
        idealtimestep = self.stability()
        if not self.time_step: 
            self.time_step  = idealtimestep
            print('-----------------------------------------------------------')            
            print(f' No specific time constant has been indicated. \n '\
                  f'The stability region has been calculated and an appropriate timestep has been chosen.\n '\
                  f'Timestep = {idealtimestep:.2e} s')
            print('-----------------------------------------------------------') 
        if (self.time_step-idealtimestep)/idealtimestep > 0.1: 
            print('-----------------------------------------------------------')            
            print(f'The manually chosen time step of {self.time_step:.2e} is eventually too big and could cause instabilities in the simulation.\n '\
                  f'We suggest a timestep of {idealtimestep:.2e} s')
            print('-----------------------------------------------------------')
        if(self.time_step-idealtimestep)/idealtimestep < -0.2: 
            print('-----------------------------------------------------------')  
            print(f'The maunually chosen time step of {self.time_step:.2e} is very small and will eventually cause a long simulation time.\n'\
                  f'We suggest a timestep of {idealtimestep:.2e} s')
            print('-----------------------------------------------------------')            
        #loading simulation relevant properties from the structural tmeperature object
        [c,A00,Abig,A1b,A2b,Cb,length,N,plp,xflat,x_plt_flat,initphi_large,interfaces,LayerMat,A1h]  = self.temp_data.Msetup()
        t = np.arange(self.start_time,self.final_time,self.time_step)
        #only if the injection would make the time grid smaller, to not move into instable regime
        if self.source.FWHM == True:
            if (6*self.source.FWHM/200 < idealtimestep):
                #inject 200 extra points around pulse to fully capture the shape of the pulse
                tinj    = np.linspace(self.source.t0 - 3*self.source.FWHM,self.source.t0 + 3*self.source.FWHM,200)
                smaller = np.where(t<self.source.t0 - 3*self.source.FWHM)[0]
                bigger  = np.where(t>self.source.t0 + 3*self.source.FWHM)[0]
            #new time grid with higher resolution
                t = np.concatenate((t[smaller],tinj,t[bigger]),axis=0)
        tstep = np.ones(len(t))
        tstep[:-1] = np.diff(t); tstep[-1] = np.diff(t)[-1]
        #If a more refined grid is choosen around t0. We inject a fine time grid around t0, to correctly capture the pulse shape
        if self.source.adjusted_grid is not False: 
            if self.source.dt0 == False:
                print('-----------------------------------------------------------')            
                print('The option for an adjusted grid is True, but no interval for a more refined grid has been given.'/
                      'Indicate dt0 (around which values the time grid should have higher resolution) in the source object')
                print('-----------------------------------------------------------')

            if 2*self.source.dt0/self.source.extra_points < idealtimestep:
                print('-----------------------------------------------------------')            
                print('A refined Grid around t0 has been applied')
                print('-----------------------------------------------------------')
                tinj    = np.linspace(self.source.t0-self.source.dt0,self.source.t0+self.source.dt0,self.source.extra_points)
                smaller = np.where(t<self.source.t0 - self.source.dt0)[0]
                bigger  = np.where(t>self.source.t0 + self.source.dt0)[0]
                #new time grid with higher resolution
                t = np.concatenate((t[smaller],tinj,t[bigger]),axis=0)
                tstep = np.ones(len(t))
                tstep[:-1] = np.diff(t); tstep[-1] = np.diff(t)[-1]
            else: 
                print('-----------------------------------------------------------')            
                print('No refined time grid is applied. The timestep is alerady very small.' \
                      'You can use the simulation class with the property self.time_step and '\
                      'assign it to a smaller value as the current time step.')
                print('-----------------------------------------------------------')
                
        #If the two temperature system is considered the matrix setup for the lattice is executed
        if self.temp_data_Lat:
            self.temp_data_Lat.Msetup()
            
        #The source gets loaded and made into a matrix=======
        """
        First we load from the source class the module for the 
        Gaussian source type. Then we do a check if for multiple
        layers, also multiple optical penetration depths are given. 
        Then the method 'init_G_source(xflat,t,opt_pen,N,func)' 
        gets called, where func is the Gaussian function, defined
        in the 'source' class. 
        """
        if self.source.sourcetype == 'Gaussian':
            print('-----------------------------------------------------------')
            print('Gaussian source type is applied')
            print('-----------------------------------------------------------')
            #typechecking
            typecheck = np.array([])
            #in case of a number (int or float)
            if type(self.source.optical_penetration_depth) is not (list or type(typecheck)): 
                self.source.optical_penetration_depth = np.array([self.source.optical_penetration_depth])
                #if given as a list
            if type(self.source.optical_penetration_depth) is list: 
                self.source.optical_penetration_depth = np.asarray(self.source.optical_penetration_depth)
            #To maintain form after multiple calls
            if np.shape(self.source.optical_penetration_depth) == (1,len(length)-1):
                self.source.optical_penetration_depth = self.source.optical_penetration_depth[0]
            #dimensioncheck    
            if (len(self.source.optical_penetration_depth) is not len(length)-1):
                self.source.optical_penetration_depth = self.source.optical_penetration_depth*np.ones(len(length)-1)
                print('-----------------------------------------------------------')
                print(f'Not every layer has a unique optical penetration depth. \n'\
                      f'=> optical penetration depth will be set to the value of the first layer = '\
                      f' {self.source.optical_penetration_depth[0]:.2e} \n  for all other layers.')
                print('-----------------------------------------------------------')
            #Source with different optical penetration depth defined on the xflat gird   
            sourceM = self.init_G_source(xflat,t,self.source.optical_penetration_depth,N,self.source.Gaussian)
 
        if self.source.sourcetype == 'THz': 
            """
            here I want to consider the simplistic model of a THz pulse. 
            Data of the time evolution and amplitude has to be loaded into the source class
            exponential decay in space is considered
            """
            if self.source.loadData is False: 
                print('-----------------------------------------------------------')
                print('Import an array, containing the data of the THz pulse. arr[0,:] = time; arr[1,:] = amplitude')
                print('-----------------------------------------------------------')
            else:  
                print('-----------------------------------------------------------')
                print('THz time profile of source is taken into account.')
                print('-----------------------------------------------------------')
                lam = 1/self.source.optical_penetration_depth
                [ttime,amplitude] = self.source.loadData
                sourceM = self.source.THz(t,xflat,ttime,amplitude,lam)
        if self.source.sourcetype == 'TMM': 
            """
            This will implement a transfer matrix approach to local absorption
            instead as using the Lambert Beer´s law considered in the Gaussian
            source type
            """
            print('-----------------------------------------------------------')
            print('Transfer matrix absorption profile of source is taken into account.'\
                  f'Length of every layer has to be given in units of meter.')
            print('-----------------------------------------------------------')
            #!!I am assuming the input for the disstances is given in nm
            #!!Then I am multiplying with 1e9, since the absorption()-function
            #!!In the source module only works if length is in units of nm!!
            sourceM = self.source.createTMM(xflat,t,self.temp_data.length*1e9)
# =============================================================================
#             xmg, tmg = np.meshgrid(xflat,t)
#             x0      = self.temp_data.length*1e9 #All the interfaces
#             #Adding infinite layer at the beginning and at the end of the distance vector
#             d_vec   = np.diff(x0); d_vec = np.insert(d_vec,(0,len(d_vec)),np.inf) 
#             #Calculating the 1D absorption profile according to TMM
#             [absorption,grid] = self.source.absorption(self.source.theta_in,self.source.lambda_vac,self.source.n_vec,d_vec,self.source.polarization,np.shape(xflat)[0])
#             #Convert the input, fluence & FWHM given in 'source' class to Amplitude and sigma
#             sigma2  = self.source.FWHM**2/(2*np.log(2))
#             A       = self.source.fluence/np.sqrt(2*np.pi*sigma2)
#             #Evaluate the Gaussian time profile
#             sourceM = A*np.exp(-(tmg-self.source.t0)**2/(2*sigma2))
#             #Multiplying it with the absorption profile to obtain 2D (space-time source map)
#             sourceM *= absorption/1e-9
# =============================================================================

        #Making the boundary conditions a function of t, in case they are given as scalars
        if isinstance(self.left_BC,(int,float)): 
            dummy = self.left_BC
            self.left_BC = lambda t: dummy + 0*t
        if isinstance(self.right_BC,(int,float)): 
            dummy1 = self.right_BC
            self.right_BC = lambda t: dummy1 + 0*t
        #Makint the boundary conditions a matrix for the electron case      
        BC_E      = np.zeros((len(c),len(t))) 
        BC_E[0]   = self.left_BC(t)
        BC_E[-1]  = self.right_BC(t)
        
        if self.temp_data_Lat:
            if isinstance(self.left_BC_L,(int,float)): 
                dummy2 = self.left_BC_L
                self.left_BC_L = lambda t: dummy2 + 0*t 
            if isinstance(self.right_BC_L,(int,float)): 
                dummy3 = self.right_BC_L
                self.right_BC_L = lambda t: dummy3 + 0*t 
        #Makint the boundary conditions a matrix for the lattice case
            BC_L      = np.zeros((len(c),len(t))) 
            BC_L[0]   = self.left_BC_L(t)
            BC_L[-1]  = self.right_BC_L(t)
            
            if np.size(self.coupling)<np.size(length)-1:
                self.coupling = self.coupling*np.ones(np.size(self.temp_data.length)-1)
                print('-----------------------------------------------------------')
                print(f'Not every layer has a unique coupling constant \'G \'. \n => G will be set to the value of the first layer = {self.coupling[0]:.2e}\n  for all other layers.')
                print('-----------------------------------------------------------')
        
        if self.temp_data_Lat: #The two temperature model is considered
            #Setup arrays for electron temperature
            c_E = np.copy(c)
            phi_E = np.zeros((len(t),len(x_plt_flat))); phi_E[0] = initphi_large
            Flow_1E = np.zeros(len(c))
            Flow_2E = np.zeros(len(c)) 
            dphi_E = np.zeros(len(c))
            intphi_E = np.zeros(len(c))
            #Setup arrays for lattice temperature
            c_L = np.copy(c)
            phi_L = np.zeros((len(t),len(x_plt_flat))); phi_L[0] = initphi_large
            Flow_1L = np.zeros(len(c))
            Flow_2L = np.zeros(len(c)) 
            dphi_L = np.zeros(len(c))
            intphi_L = np.zeros(len(c))
            #General setup for E.E. loop
            condi       = np.array([np.arange(1,len(length)-1)])*(N-1) #Index to apply interface condition
            cnfill      = np.array([np.arange(1,len(length)-1)])*(plp-1)#correct interface condition with real value for phi
            A00[0]      = 1; A00[-1] = 1 #Avoide devide through 0 in dphi_L! Clar for BC before intphi calc.
            Abig_E      = np.copy(Abig)#Since Abig can change due to interconditions we split it here
            Abig_L      = np.copy(Abig)# The interface conditions are applied on every time step 
            start_EL    = time.time()
            for i in tqdm(range(1,len(t))):
                phi0_E = np.dot(A00,c_E); phi1_E = np.dot(A1b,c_E); phi2_E = np.dot(A2b,c_E)
                phi0_L = np.dot(A00,c_L); phi1_L = np.dot(A1b,c_L); phi2_L = np.dot(A2b,c_L)
                intercon_E = self.interconditions(phi_E[i-1],interfaces,self.temp_data.conductivity,N,A1h)    #Applying interface conditions
                intercon_L = self.interconditions(phi_L[i-1],interfaces,self.temp_data_Lat.conductivity,N,A1h)
                startf = 0;endf = N-1
                #construct all picewise flows and piecewise dphi
                for j in range(0,interfaces+1): 
                    #electron
                    Flow_1E[startf:endf] = self.temp_data.diff_conductivity(phi0_E[startf:endf],j)
                    Flow_2E[startf:endf] = self.temp_data.conductivity[j](phi0_E[startf:endf])
                    Flow_1E[startf:endf] *=phi1_E[startf:endf]**2
                    Flow_2E[startf:endf] *= phi2_E[startf:endf] 
                    #lattice
                    Flow_1L[startf:endf] = self.temp_data_Lat.diff_conductivity(phi0_L[startf:endf],j)
                    Flow_2L[startf:endf] = self.temp_data_Lat.conductivity[j](phi0_L[startf:endf])
                    Flow_1L[startf:endf] *=phi1_L[startf:endf]**2
                    Flow_2L[startf:endf] *= phi2_L[startf:endf]  
                    #calculate delta phi for electron and lattice
                    dphi_E[startf:endf] = 1/(self.temp_data.heatCapacity[j](phi0_E)[startf:endf]*self.temp_data.rho[j])*(Flow_1E[startf:endf]+Flow_2E[startf:endf]+sourceM[i,startf:endf] + self.coupling[j]*(phi0_L[startf:endf]-phi0_E[startf:endf]))  
                    dphi_L[startf:endf] = 1/(self.temp_data_Lat.heatCapacity[j](phi0_L)[startf:endf]*self.temp_data_Lat.rho[j])*(Flow_1L[startf:endf]+Flow_2L[startf:endf] + self.coupling[j]*(phi0_E[startf:endf]-phi0_L[startf:endf])) 
                    startf += N-1; endf +=N-1 
                #filling up Abig wiht the interface condition in the middle of the grid
                start_i = 0; end_i = N-1  
                for k in range(0,interfaces): #Apply interface conditions for all layers in every time step 
                    #for the electron system 
                    Abig_E[end_i,start_i:end_i]    =  intercon_E[0][k][:-1]#Lhs interface flow
                    Abig_E[end_i,end_i+1:end_i+N]  = -intercon_E[1][k][1:]#Rhs interface flow
                    Abig_E[end_i,end_i]            = intercon_E[0][k][-1] -intercon_E[1][k][0]
                    #for the lattice system
                    Abig_L[end_i,start_i:end_i]    =  intercon_L[0][k][:-1]#Lhs interface flow
                    Abig_L[end_i,end_i+1:end_i+N]  = -intercon_L[1][k][1:]#Rhs interface flow
                    Abig_L[end_i,end_i]            = intercon_L[0][k][-1] -intercon_L[1][k][0]
                    start_i += N-1; end_i += N-1
                #computing the flux for every time step for the boundaries
                Flux_E = BC_E[:,i]
                Flux_E[0] /= self.temp_data.conductivity[0](c_E[0])   + 1e-12
                Flux_E[1] /= self.temp_data.conductivity[-1](c_E[-1]) + 1e-12
                Flux_L = BC_L[:,i]
                Flux_L[0] /= self.temp_data_Lat.conductivity[0](c_L[0])   + 1e-12
                Flux_L[1] /= self.temp_data_Lat.conductivity[-1](c_L[-1]) + 1e-12
                #Clear for boundary conditions at the edgeds of the grid
                dphi_E[0] = 0; dphi_E[-1] = 0; dphi_L[0] = 0; dphi_L[-1] = 0
                phi0_E[0] = 0; phi0_E[-1]  = 0; phi0_L[0] = 0; phi0_L[-1] = 0;   
                #intermediate phi with low resolution in space according to explicit euler
                intphi_E        = phi0_E + tstep[i] * dphi_E + Flux_E
                intphi_E[condi] = 0
                intphi_L        = phi0_L + tstep[i] * dphi_L + Flux_L
                intphi_L[condi] = 0                               #Interface condition: Q_1 -Q_2 = 0
                #electron: use c to map on high resolution x-grid
                #since in Abig k(T(t)) is inserted we have to solve the system for every step
                c_E                 = np.linalg.solve(Abig_E,intphi_E)  # c(t) for every timestep
                phi_E[i]            = np.dot(Cb,c_E)                    # map spline coefficients to fine Cb grid
                phi_E[i,cnfill]     = c_E[condi]                        #correct the values for phi at interface
                #lattice
                c_L                 = np.linalg.solve(Abig_L,intphi_L)   
                phi_L[i]            = np.dot(Cb,c_L)
                phi_L[i,cnfill]     = c_L[condi]             
            end_EL = time.time()
            print('-----------------------------------------------------------')  
            print('Heat diffusion in a coupled electron-lattice system has been simulated')
            print('Eleapsed time in E.E.- loop:', end_EL-start_EL)
            print('-----------------------------------------------------------')  
            return(phi_E,phi_L,x_plt_flat,t)
        else: #this is the single temperature case. (Only electron temperature)
            #prepare space to store phi solution on fine plt grid. And Flow_1,2 vectors
            phi = np.zeros((len(t),len(x_plt_flat))); phi[0] = initphi_large
            Flow_1 = np.zeros(len(c))
            Flow_2 = np.zeros(len(c)) 
            dphi = np.zeros(len(c))
            intphi = np.zeros(len(c))
            condi = np.array([np.arange(1,len(length)-1)])*(N-1)    #Index to apply interface condition
            cnfill = np.array([np.arange(1,len(length)-1)])*(plp-1) #correct interface condition with real value for phi
            A00[0] = 1; A00[-1] = 1 #Avoid 1/0 division in dphi calculation. See E.E. loop
            
            startE = time.time()
            for i in range(1,len(t)):
                phi0 = np.dot(A00,c); phi1 = np.dot(A1b,c); phi2 = np.dot(A2b,c)  
                intercon_E = self.interconditions(phi[i-1],interfaces,self.temp_data.conductivity,N,A1h) #get interface conditions for every time step
                startf = 0;endf = N-1
                #construct all picewise flows and piecewise dphi
                for j in range(0,interfaces+1): 
                    Flow_1[startf:endf] = self.temp_data.diff_conductivity(phi0[startf:endf],j)
                    Flow_2[startf:endf] = self.temp_data.conductivity[j](phi0[startf:endf])
                    Flow_1[startf:endf] *=phi1[startf:endf]**2
                    Flow_2[startf:endf] *= phi2[startf:endf] 
                    dphi[startf:endf] = 1/(self.temp_data.heatCapacity[j](phi0)[startf:endf]*self.temp_data.rho[j])*(Flow_1[startf:endf]+Flow_2[startf:endf]+sourceM[i,startf:endf])  
                    startf += N-1; endf +=N-1 
                #filling up Abig wiht the interface condition in the middle of the grid
                start_i = 0; end_i = N-1  
                for k in range(0,interfaces): #Apply interface conditions for all layers in every time step 
                    Abig[end_i,start_i:end_i]    =  intercon_E[0][k][:-1]#Lhs interface flow
                    Abig[end_i,end_i+1:end_i+N]  = -intercon_E[1][k][1:]#Rhs interface flow
                    Abig[end_i,end_i]            = intercon_E[0][k][-1] -intercon_E[1][k][0]
                start_i += N-1; end_i += N-1
                #computing the flux for every time step for the boundaries
                Flux_E = BC_E[:,i]
                Flux_E[0] /= self.temp_data.conductivity[0](c[0])   + 1e-12
                Flux_E[1] /= self.temp_data.conductivity[-1](c[-1]) + 1e-12
                #Make space for BC    
                dphi[0] = 0; dphi[-1] = 0
                phi0[0] = 0; phi0[-1]  = 0
                
                intphi = phi0 + tstep[i] * dphi + Flux_E
                intphi[condi] = 0                               #Interface condition: Q_1 -Q_2 = 0
                # c(t) for every timestep
                #since in Abig k(T(t)) is inserted we have to solve the system for every step
                c               = np.linalg.solve(Abig,intphi)  #this system has to be solved in every time step
                phi[i]          = np.dot(Cb,c)                  # map spline coefficients to fine Cb grid
                phi[i,cnfill]   = c[condi]                      #correct the values for phi at interface
            endE = time.time()
            print('-----------------------------------------------------------')  
            print('Electron temperature heat diffusion has been simulated.')
            print('Eleapsed time in E.E.- loop:', endE-startE)
            print('-----------------------------------------------------------')  
        return(phi,x_plt_flat,t)
        
    def stability(self):
        """
        If only the electron temperature system is under consideration, we only 
        compute the eigenvalues of lambda_i = k/(C*rho)*A00^-1*A2b. This is 
        we consider the minimum Eigenvalue for each layer to represent the time konstant.
        The time constant for E.E. is then given by -2/min(lambda_i), which is 
        the criterion for stability for E.E. loops, to obtain convergence.
        """
        [c,A00,Abig,A1b,A2b,Cb,length,N,plp,xflat,x_plt_flat,initphi_large,interfaces,LayerMat,A1h]  = self.temp_data.Msetup()
        A00[0,0] = 1; A00[-1,-1] = 1
        rho_E = self.temp_data.rho
        conductivity_E = self.temp_data.conductivity
        conductivity_E = np.asarray(conductivity_E)
        typecheck      = np.array([1])[0]
        for i in range(0,len(conductivity_E)):
            #In case conductivity is a function k(T) we compute a worst case scenario
            #this is because we can only compare integers. 
            if not isinstance(conductivity_E[i],(int,float,type(typecheck))): 
                testT = np.linspace(270,2000,50)
                conductivity_E[i] = max(conductivity_E[i](testT))
        
        heatCapacity_E = self.temp_data.heatCapacity
        heatCapacity_E = np.asarray(heatCapacity_E)
        for i in range(0,len(heatCapacity_E)): 
            #In case heatCapacity is a function C(T) we compute a worst case scenario
            #and take an integer value to compare
            if not isinstance(heatCapacity_E[i],(int,float,type(typecheck))): 
                testT = np.linspace(270,2000,50)
                heatCapacity_E[i] = min(heatCapacity_E[i](testT))
        Eval = np.zeros(interfaces+1) #for each layer there will be an eigenvalue    
        koeff1   = conductivity_E/(heatCapacity_E*rho_E)
        for i in range(0,interfaces+1):
            Lambda   = koeff1[i]*LayerMat[i]
            Eval[i]  = min(np.real(np.linalg.eig(Lambda)[0]))
        tkonst_E = -2/Eval
        
        if self.num_of_temp == 2:
            """
            In the multy temperature case, we also consider the lattice dynamics,
            with respective k_L/(C_L*rho_L) dynamics. In addition, we also have to
            consider, the coupling between those two layers. I.e. G_mat. 
            with coefficients koeff_2 = G/(heatCapacity*rho)
            Therefor we compute eigenvalues of the combined system:
            lambda_i = eval(Lambda + G_mat) for each layer.
            The time constant is again -2/min(lambda_i)
            """
            [c,A00_L,Abig,A1b,A2b_L,Cb,length,N,plp,xflat,x_plt_flat,initphi_large,interfaces,LayerMat,A1h]  = self.temp_data_Lat.Msetup()
            A00_L[0,0] = 1; A00_L[-1,-1] = 1  
            rho_L = self.temp_data_Lat.rho
            G = self.coupling
            G = np.asarray(G)
             
            conductivity_L = self.temp_data_Lat.conductivity
            conductivity_L = np.asarray(conductivity_L)
            #In case conductivity is a function k(T) we compute a worst case scenario
            for i in range(0,len(conductivity_L)): 
                if not isinstance(conductivity_L[i],(int ,float,type(typecheck))):
                    testT = np.linspace(270,2000,50)
                    conductivity_L[i] = max(conductivity_L[i](testT))
                    
            heatCapacity_L = self.temp_data_Lat.heatCapacity
            heatCapacity_L = np.asarray(heatCapacity_L)
            #In case heatCapacity is a function C(T) we compute a worst case scenario
            for i in range(0,len(heatCapacity_L)): 
                if not isinstance(heatCapacity_L[i],(int,float,type(typecheck))): 
                    testT = np.linspace(270,2000,50)
                    heatCapacity_L[i] = min(heatCapacity_L[i](testT))
            #M: for every layer we loade the respective matrix from the temperature class                    
            M                       = np.shape(LayerMat)[1] 
            Lambda                  = np.zeros((2*M,2*M))
            Eval                    = np.zeros(interfaces+1)
            G_mat                   = np.zeros((2*M,2*M))          
            koeff1_E                = conductivity_E/(heatCapacity_E*rho_E)
            koeff1_L                = conductivity_L/(heatCapacity_L*rho_L)
            koeff2_E                = G/(heatCapacity_E*rho_E)            
            koeff2_L                = G/(heatCapacity_L*rho_L)
            
            for i in range(0,interfaces+1):
                Lambda[0:M,0:M]     = koeff1_E[i]*LayerMat[i]
                Lambda[M:,M:]       = koeff1_L[i]*LayerMat[i]
                G_mat[0:M,0:M]      = -koeff2_E[i]*np.eye(M)
                G_mat[0:M,M:]       = koeff2_E[i]*np.eye(M)
                G_mat[M:,0:M]       = koeff2_L[i]*np.eye(M)
                G_mat[M:,M:]        = -koeff2_L[i]*np.eye(M)
                Eval[i]             = min(np.real(np.linalg.eig(Lambda+G_mat)[0]))
            tkonst = -2/Eval
            return(min(tkonst))            
        else: 
            #if there is only electron temperature, only those dynamics will be
            #considered, when time step for the E.E. loop is calculated.
            return(min(tkonst_E))
        

class source(object):

    def __init__(self): 
        self.sourcetype     = 'Gaussian'
        self.fluence        = []
        self.optical_penetration_depth = []
        self.t0             = 0
        self.FWHM           = False
        self.loadData       = False
        self.multipulse     = False
        self.frequency      = False
        self.num_of_pulses  = False
        self.adjusted_grid  = False
        self.dt0            = False
        self.extra_points   = 200
        self.theta_in       = 0# 0 is perpendicular to surface/ pi/2 is grazing
        self.n_vec          = False
        self.lambda_vac     = False
        self.polarization   = 'p' 
        
               
    def getProperties(self): # to depict the properties of the object
        for i in (self.__dict__): 
            print(i,' : ',self.__dict__[i])
            
    def __repr__(self): 
        return('Source')
            
    def Gaussian(self,xmg,tmg,lam,A,sigma2,x0):
        if not (self.fluence or self.optical_penetration_depth or self.FWHM):
            print('------------------------------------------------------------')
            print('Create a pulse with defining pulse properties. \n ' +\
                  '.fluence, .optical_penetration_depth, .FWHM')
            print('------------------------------------------------------------')
        #Create a source with respect to each lam of every layer
        Gauss = A*np.exp(-(tmg-self.t0)**2/(2*sigma2))
        Gauss *= lam*np.exp(-lam*(xmg-x0))
        
        if (self.multipulse is not False and self.num_of_pulses is False):
            time_range  = tmg[-1,-1]-self.t0
            pulses      = int(round(time_range * self.frequency)) 
            #Add up Gaussian pulses with different t0, according to the frequency given
            #from t0 onwards, until the end of the time grid
            Gauss = np.zeros(np.shape(tmg))
            for i in range(0,pulses): 
                t00 = self.t0 + i/self.frequency
                Gauss += self.fluence/self.FWHM*np.sqrt(np.log(2)/np.pi)*np.exp(-(tmg-t00)**2*np.log(2)/(self.FWHM**2))*lam*np.exp(-lam*xmg) 
                 
        if (self.multipulse is not False and self.num_of_pulses is not False):
            #Creating a certain number of pulses according to self.num_of_pulses
            time_range = tmg[-1,-1]-self.t0
            pulses = self.num_of_pulses
            #If num_of_pulses is bigger too big to fit in the timerange [t0,t_end] throw warning
            if (pulses > int(round(time_range * self.frequency))):
                pulses      = int(round(time_range * self.frequency)) 
                print('Number of pulses is too big to fit in the timerange under consideration. \n'\
                      'Adjust t_end or consider a smaller number of pulses.')
            Gauss = np.zeros(np.shape(tmg))
            for i in range(0,pulses):
                t00 = self.t0 +i/self.frequency
                Gauss += self.fluence/self.FWHM*np.sqrt(np.log(2)/np.pi)*np.exp(-(tmg-t00)**2*np.log(2)/(self.FWHM**2))*lam*np.exp(-lam*xmg) 
        
        return(Gauss)
    #mytime is the timegrid of the simulation
    #time, amplitude are the timegrids of the inputdata collected from the lab    
    def THz(self,mytime,myspace,ttime,amplitude,lam):
        #Mapping the obtained data to the simulation time grid via interpolation
        ampl1D = np.interp(mytime,ttime,amplitude**2)
        #Compute the right amplitude. using the area under the curve
        integr = np.trapz(ampl1D,mytime,np.diff(mytime))
        #scling factore to get the amplitude right
        scaling = self.fluence/integr
        xmg,tmg = np.meshgrid(myspace,mytime)
        ampl2D    = np.interp(tmg,ttime,amplitude**2)
        ampl2D *= scaling
        #exponential decay in space
        lam = 1/self.optical_penetration_depth
        ampl2D *= lam*np.exp(-lam*xmg)
        return(ampl2D)
        
        
        
        
        
    def fresnel(self,theta_in,n_in,n_out,pol): 
        n_in = complex(n_in); n_out = complex(n_out)
        theta_out = np.arcsin(n_in*np.sin(theta_in)/n_out)
        if pol == 's': 
            rs = (n_in*np.cos(theta_in) - n_out*np.cos(theta_out))/\
                 (n_in*np.cos(theta_in) + n_out*np.cos(theta_out))
            ts = 2*n_in*np.cos(theta_in)/(n_in*np.cos(theta_in)+n_out*np.cos(theta_out))
            return(theta_out,rs,ts)
        if pol == 'p':
            rp = (n_out*np.cos(theta_in)-n_in*np.cos(theta_out))/\
                 (n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
            tp = 2* n_in*np.cos(theta_in)/(n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
            return(theta_out,rp,tp) 
    
    
    def TM(self,theta_in,lambda0,n_vec,d_vec,pol): 
        #create complex arrays for variables
        theta   = np.zeros(len(n_vec), dtype = complex); theta[0] = theta_in
        phi     = np.zeros(len(n_vec),dtype = complex)
        rn      = np.zeros(len(n_vec)-1, dtype = complex) 
        tn      = np.zeros_like(rn,dtype = complex)
        M_n     = np.zeros((len(n_vec),2,2), dtype = complex)
        M       = np.eye(2,dtype = complex)
        for i in range(len(n_vec)-1): # to obtian all angels/rn/tn for each layer
            [theta[i+1],rn[i],tn[i]] = self.fresnel(theta[i],n_vec[i],n_vec[i+1],pol)
        #M = M0*M1*M2*M4*....
        for k in range(1,len(n_vec)-1):#loop over all interfaces except 1st
            phi[k]  = 2*np.pi*n_vec[k]*np.cos(theta[k])*d_vec[k]/lambda0
            Tn      = np.array([[np.exp(-1j*phi[k]),0],[0,np.exp(1j*phi[k])]],dtype = complex)/tn[k]
            Pn      = np.array([[1,rn[k]],[rn[k],1]],dtype = complex)
            M_n[k]     = np.dot(Tn,Pn)
            M = np.dot(M,M_n[k])
        #compute for the first interface: 
        trans0 = np.array([[1,rn[0]],[rn[0],1]],dtype= complex)/tn[0]
        M = np.dot(trans0,M)
        #Complex transmission/reflection amplitude
        t = 1/M[0,0]
        r = M[1,0]/M[0,0]
        #Fraction of power transmitted
        if pol == 's': #s-polarized
            T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(theta[-1]))/\
                np.real(n_vec[0]*np.cos(theta[0]))
        elif pol == 'p': #p-polarized
            T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(np.conj(theta[-1])))/\
                np.real(n_vec[0]*np.cos(np.conj(theta[0])))
        #Fraction of power reflected
        R = np.abs(r)**2
        A = 1.-T-R
        return(M,M_n,t,r,T,R,A,theta)
    
    def layerAmpl(self,theta_in,lambda0,n_vec,d_vec,pol): 
        """
        After r & t have been calculated and all the respective matrices M_n
        for each layer are known, we can go 'backwards', i.e. from the last to the
        first layer, and compute all the amplituedes for the forward v_n and 
        backward w_n traveling wave. -> [v_n,w_n].T = M_n @ [v_{n+1},w_{n+1}].T
        """
        [M,M_n,t,r,T,R,A,theta] = self.TM(theta_in,lambda0,n_vec,d_vec,pol)
        vw_list = np.zeros((len(n_vec),2),dtype = complex)
        vw =np.array([[t],[0]])
        vw_list[-1,:] = vw.T
        for i in range(len(n_vec)-2,0,-1):
            vw = np.dot(M_n[i],vw)
            vw_list[i,:] = vw.T
        return(vw_list,theta)
        
    def absorption(self,theta_in,lambda0,n_vec,d_vec,pol,points):
        #reload the forward and backward wave coefficients for every layer
        [vw_n,theta]  = self.layerAmpl(theta_in,lambda0,n_vec,d_vec,pol)
        total_len = np.sum(d_vec[1:-1])
        pointcount = 0
        #a is an array where the normalized absorbtion for the entire grid is stored 
        a = []
        for i in range(1,len(n_vec)-1):
            kz      = 2*np.pi*n_vec[i]*np.cos(theta[i])/lambda0
            #How many points, out of the total 'points' does each layer get, with respect to 
            #the length of the layer
            if i == 1: 
                points_per_layer = int(np.round(points/(len(d_vec)-2)))+1
            else:
                points_per_layer = int(np.round(points/(len(d_vec)-2)))
            #for every layer, the space grid is reset. I.e. every layer starts at 0
            layer   = np.linspace(0,d_vec[i],points_per_layer)
            v = vw_n[i,0]; w = vw_n[i,1];#complex wave amplitudes for every layer
            Ef = v * np.exp(1j * kz * layer)#forward traveling wave 
            Eb = w * np.exp(-1j * kz *layer)#backward traveling wave
            if pol == 'p':#p-polarized
                a_layer = (n_vec[i]*np.conj(np.cos(theta[i]))*(kz*np.abs(Ef-Eb)**2-np.conj(kz)*np.abs(Ef+Eb)**2)).imag /\
                        (n_vec[0]*np.conj(np.cos(theta[0]))).real
            if pol == 's': 
                a_layer = (np.abs(Ef+Eb)**2 *np.imag(kz*n_vec[i]*np.cos(theta[i])))/\
                        (np.real(n_vec[0]*np.cos(theta[0])))
            #for every layer calculate the absorbtion grid and append it to the total 
            a   = np.append(a,a_layer)
            #to get the right amount of points considered in the grid, since we round.
            pointcount += points_per_layer 
        a = np.real(a)
        grid = np.linspace(0,total_len,pointcount)
        return(a,grid)
        
        
    def createTMM(self,xflat,t,x0): 
        #The distances and hence x0 have to be given in units of nm!
        xmg, tmg = np.meshgrid(xflat,t)
        #Adding infinite layer at the beginning and at the end of the distance vector
        d_vec   = np.diff(x0); d_vec = np.insert(d_vec,(0,len(d_vec)),np.inf) 
        #Calculating the 1D absorption profile according to TMM
        [absorption,grid] = self.absorption(self.theta_in,self.lambda_vac,self.n_vec,d_vec,self.polarization,np.shape(xflat)[0])
        #Convert the input, fluence & FWHM given in 'source' class to Amplitude and sigma
        sigma2  = self.FWHM**2/(2*np.log(2))
        A       = self.fluence/np.sqrt(2*np.pi*sigma2)
        #Evaluate the Gaussian time profile
        sourceM = A*np.exp(-(tmg-self.t0)**2/(2*sigma2))
        #Multiplying it with the absorption profile to obtain 2D (space-time source map)
        sourceM *= absorption/1e-9
        return(sourceM)
    
        
        

        
        
class visual(object): 
    
    def __init__(self,*args):
        self.data = False
        typecheck = np.array([])
        if isinstance(args[0],(float,int,list,type(typecheck)))== True:
            #If arrays are passed on
            if len(args) == 4:
                self.T_E = args[0]
                self.T_L = args[1]
                self.x = args[2]
                self.t = args[3]
            if len(args) == 3:
                self.T_E = args[1]
                self.x = args[2]
                self.t = args[3]
            if len(args) < 3: 
                print('Not enough input arguments are given. \n Pass on \'Electron temperature\', \'x-grid\', \'time grid\'. \n They are output of simulation.run(). ')
            if len(args)>5: 
                print('Too many input arguments are given. \n Pass on \'Electron temperature\', \'x-grid\', \'time grid\'. \n They are output of simulation.run(). ')
        #If the simulation object is passed on          
        else:
            self.sim = args[0]
            print('------------------------------------------------------------')
            print('The simulation object of the'+str(self.sim.num_of_temp)+' temerature system has been passed on to the visual class.')
            print('------------------------------------------------------------')    
            if (self.sim.num_of_temp == 2 and self.data == False):  #two temperature case  
                [self.T_E, self.T_L, self.x, self.t]    = self.sim.run()
                self.so                                 = self.sim.source
                self.data                               = True
            if (self.sim.num_of_temp == 1 and self.data == False):  #1 temperature case  
                [self.T_E, self.x, self.t]              = self.sim.run()
                self.so                                 = self.sim.source 
                self.data                               = True
    def __repr__(self): 
        return('Visual')
        
    def source(self):
        opt_pen = self.so.optical_penetration_depth
        if self.so.sourcetype == "Gaussian":
        #Source with different optical penetration depths
        #defined on the x_plt and time grid. 
        #Call the initalisation 'init_G_source' = method in 'simulation' class
        #The 'init_G_source' needs the Gauss function of the 'source' class in the source class
            init_Gauss_is = self.sim.init_G_source
            N             = self.sim.temp_data.plt_points
            s             = init_Gauss_is(self.x,self.t,opt_pen,N,self.so.Gaussian)
        if self.so.sourcetype == "THz": 
            ttime       = self.so.loadData[0] 
            amplitude   = self.so.loadData[1]
            s           = self.so.THz(self.t,self.x,ttime,amplitude,1/opt_pen)
        if self.so.sourcetype == "TMM": 
            s           = self.so.createTMM(self.x,self.t,self.sim.temp_data.length*1e9)
            
        xx,tt = np.meshgrid(self.x,self.t)
        fig  = plt.figure()
        ax   = fig.gca(projection='3d')
        if (np.mean(self.x) <1e-7 and np.mean(self.t)<1e-10):
            surf = ax.plot_surface(xx/1e-9,tt/1e-12,s,cmap = 'jet')
            plt.xlabel('Depth from surface into material in nm')
            plt.ylabel('time  in ps')
        else: 
            surf = ax.plot_surface(xx,tt,s,cmap = 'jet')
            plt.xlabel('Depth from surface into material in m')
            plt.ylabel('time in s')   
        fig.colorbar(surf,shrink=0.7, aspect=5)
        plt.title(r'S(x,t)')
        plt.show()
        #A contour plot
        plt.figure()
        plt.title("Contour plot of input source")
        plt.xlabel("Depth into material from surface in m")
        plt.ylabel("Time profile in s")
        plt.contourf(xx,tt,s,cmap = 'jet')
        plt.colorbar()
        plt.show()
        return(s)
        
    def contour(self,*args):
        if 'Lattice' in args:
            T = self.T_L
            name = args
            if self.sim.num_of_temp != 2:
                print('The lattice temperature can not be considered, since it has not been simulated yet. \n Initialize lattice parameters and simulate them first.')
        else: 
            T = self.T_E
            if args:
                name = args
        plt.figure()
        plt.title(r'$T(x,t)$ in K')
        plt.xlabel('x-Space in m')
        plt.ylabel('time in s')
        xx,tt = np.meshgrid(self.x,self.t)
        if (np.mean(self.x) <1e-7 and np.mean(self.t)<1e-10):
            xx,tt = np.meshgrid(self.x/1e-9,self.t/1e-12) 
            plt.xlabel('x-Space in nm')
            plt.ylabel('time in ps')
        plt.contourf(xx,tt,T,50,cmap = 'plasma')
        plt.colorbar( orientation='vertical', shrink=0.8)
        if args:
            plt.title(r'The temperature of '+str(name)+' $T(x,t)$ in K')    
        plt.show()
    
    def average(self):
        #load Data for weighted means
        points          = self.sim.temp_data.plt_points
        tot_length      = self.sim.temp_data.length[-1]
        len_of_layer    = np.diff(self.sim.temp_data.length)
        if self.sim.num_of_temp == 2:#electron and lattice temperature under consideration
            #Take weighted averages with respect to the length of each layer
            avT_E = len_of_layer[0]/tot_length*np.mean(self.T_E[:,0:points-1],1)
            avT_L = len_of_layer[0]/tot_length*np.mean(self.T_L[:,0:points-1],1)
            for i in range(1,len(len_of_layer)):
                avT_E += len_of_layer[i]/tot_length*np.mean(self.T_E[:,i*points:(i+1)*points],1)
                avT_L += len_of_layer[i]/tot_length*np.mean(self.T_L[:,i*points:(i+1)*points],1)
            plt.figure()
            plt.xlabel(r'Time in s')
            plt.ylabel(r'Temperature in K')
            if np.mean(self.t) < 1e-10: 
                self.t/=1e-12
                plt.xlabel(r'Time in ps')
                plt.ylabel(r'Temperature in K')   
            plt.plot(self.t,avT_L,c = 'k',label = 'Lattice')
            plt.plot(self.t,avT_E,c = 'r',label = 'Electron')
            plt.title('Temperature averaged in space vs time')
            plt.grid()
            plt.legend()
            plt.show()
            return(self.t,avT_E,avT_L)
        else:#only electron temperature under consideration
            #Take the weighted mean with respect to length of layer only for T_E
            avT_E = len_of_layer[0]/tot_length*np.mean(self.T_E[:,0:points-1],1)
            for i in range(1,len(len_of_layer)):
                avT_E += len_of_layer[i]/tot_length*np.mean(self.T_E[:,i*points:(i+1)*points],1)
                
            plt.figure()
            plt.xlabel(r'Time in s')
            plt.ylabel(r'Temperature in K')
            if np.mean(self.t) < 1e-10: 
                self.t/=1e-12
                plt.xlabel(r'Time in ps')
                plt.ylabel(r'Temperature in K')   
            plt.plot(self.t,avT_E,c = 'r',label = 'Electron')
            plt.title('Temperage averaged in space vs time')
            plt.grid()
            plt.legend()
            plt.show()
            #return(self.t,avT_E)
            
    #IN: vector to block, time grid, num_blocks,   
    def blocking(self,vector,tt,blocks):
        blockvec = np.zeros(blocks)
        elements = len(vector) 
        rest     = elements % blocks
        if rest != 0: #truncate vector if number of blocks dont fit in vector
            vector   = vector[0:-rest]
            tt       = tt[0:-rest]
            elements = len(vector)   
        meanA  = np.mean(vector)        
        bdata  = int((elements/blocks))#how many points per block
        sigBsq = 0; 
        for k in range(0,blocks):
            blockvec[k] = np.average(vector[k*bdata : (k+1)*bdata]) 
            sigBsq      = sigBsq + (blockvec[k]-meanA)**2    
        sigBsq *= 1/(blocks-1); 
        sigmaB = np.sqrt(sigBsq)
        error  = 1/np.sqrt(blocks)*sigmaB
        blocktt = tt[0:-1:bdata]
        return(blockvec,blocktt, error , sigmaB)
            
    def timegrid(self):
        """
        Plotting delta t_i vs t_i. This is the width of the timesteps
        with respect to the point of the timegrid where we are right now.:
        This is relevant, once adjusted_grid (in the source class )is True.
        """
        tstep = np.ones(len(self.t))
        tstep[:-1] = np.diff(self.t); tstep[-1] = np.diff(self.t)[-1]
        plt.figure()
        plt.plot(self.t,tstep,'-o')
        plt.xlabel('Time grid $t_i$'); plt.ylabel('$\Delta t_i$')
        plt.title('$\Delta t_i = t_{i+1}-t_i$ vs t_i')
            
        
            
    def animation(self,speed,sve = 0): 
        if len(self.t) % speed != 0:
            speed = speed-(len(self.t)%speed)
            print('Speed of animation',speed)
         
        fig, ax = plt.subplots()
        lineE, = ax.plot([], [],'r' , animated=True,label = 'Electron Temp. in [K]')
        #if there is only one temperature, then the lineL is just a placeholder for the update function
        if self.sim.num_of_temp == 2:
            lineL, = ax.plot([],[],'b',animated = True,label = 'Lattice Temp. in [K]')
            ax.set_xlim(0, self.x[-1]); ax.set_ylim(np.min((self.T_L,self.T_E))-(np.mean(self.T_E)-np.min((self.T_L,self.T_E)))/2, np.max((self.T_E,self.T_L))+(np.max(self.T_E)-np.mean(self.T_E))/2)
        else : 
            lineL, = ax.plot([],[],'r',animated = True)
            self.T_L = self.T_E
            ax.set_xlim(0, self.x[-1]); ax.set_ylim(np.min((self.T_L,self.T_E))-(np.mean(self.T_E)-np.min((self.T_L,self.T_E)))/2, np.max((self.T_E,self.T_L))+(np.max(self.T_E)-np.mean(self.T_E))/2)
        time_text = ax.text(0.02, 0.95, "", transform = ax.transAxes)
        plt.xlabel('Depth of Material'); plt.ylabel('Temperature')
        plt.title('Evolution of Temperature in space and time')
        plt.legend()
        
        def update(frame):
            lineE.set_data(self.x,self.T_E[speed*frame,:])
            lineL.set_data(self.x,self.T_L[speed*frame,:])
            time_text.set_text(f"time =  {self.t[speed*frame]:.2e}  s")
            if speed*frame >= len(self.t)-5: 
                ani.event_source.stop()
            return lineL,lineE, time_text
        
        ani = movie(fig, update, blit=True,save_count = 2000,interval = 15, repeat = True) 
        plt.grid()
        plt.show()
        if sve ==1: 
            ani.save("animation.mp4")
        return(ani)            
            
            
   
        
                        
        
        
        
            
    

     








        
        
        
        
        


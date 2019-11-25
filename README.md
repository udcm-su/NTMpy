# NTMpy - N-Temperature Model solver


#### To download                     `pip install NTMpy`

#### To call routines in a script:   `from NTMpy import NTMpy as ntm`

#### To update to lates version:     `pip install --upgrade NTMpy`

Further information on the solver package itselfe can be found here: [NTMpy](https://github.com/udcm-su/heat-diffusion-1D/tree/master/NTMpy).

Documentation and example sessions can be found in the [Wiki](https://github.com/udcm-su/heat-diffusion-1D/wiki).

------------------------------------------------------------------------------------------------------------------

This is a code providing a solution to the heat diffusion in a 1D structure in a 2-temperature model approximation.

We consider:
* Material comprising piecewise homogeneous layers
* Heating of electron system with an energy source with Gaussian shape (i.e. an ultrashort laser pulse)
* Three temperature model: electron, lattice and spin systems are considered and coupled
* Transfer Matrix Method to account for local absorbed energy in the multi layer material

The equation under consideration is: 

 <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/master/Pictures/Equation.PNG" width="750" height="120" />
 
 where *C* = specific heat, *k* = conductivity, *G* is the coupling constant between the two systems (Electron and Lattice)
  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiE.png" width="30" height="20" /> and <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiL.png" width="30" height="20" /> 
  are the respective temperatures of the electron and the lattice system with respect to space *x* and time *t*. The superscripts  *L* and *E* indicate wether a parameter belongs to the electron (E) or lattice (L) system and the sub index *i* denotes to which layer the parameter belongs.

 Our approach is to use a combination of **finite elements** (B-Splines) to approximate the derivation in space and **Explicit Euler** to approximate the evolution in time.
 To stay within the **stability** region of the Explicit Euler algorithm, a well suited time step is automatically computed.
 
  ### Example:
  In this case the Material under consideration is one layer of strontium ruthenate (SrRuO<sub>3</sub> or SRO) on top of a layer of strontium titanate (SrTiO<sub>3</sub> or STO).
    <p align="center"> 
   <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/SROSTO1.PNG" width="520" height="60" />     
   </p>
 The output of the solver is the temperature evolution of the electron and lattice system in space and time.
 
  Temperature evolution of electron- lattice system |  Gaussian laser pulse S(x,t) hitting sample
:-------------------------:|:-------------------------:
 <img src="https://media.giphy.com/media/RHE9DS2kPSdobin3hv/giphy.gif" width="320" height="300" />  |  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/Source.png" width="320" height="300" />
 
Here we consider an energy source *S(x,t)* = *lambda* * *exp(-lambda* * *(x-x0))* * *G(t-t0)* representing a laser pulse hitting a material from the left. 
 
The animation clearly shows the following effects: 
  1. The electron system heats up immediately.
  2. The heat in the electronic system diffuese along the x-axis.
  3. The heat transfers to the lattice system on longer time scales.

#### Documentation
Documentation and example sessions can be found in the [Wiki](https://github.com/udcm-su/heat-diffusion-1D/wiki)

#### With: 
This is a project from the Ultrafast Condensed Matter physics groupe in Stockholm. The main contributors are: 
* [Alber Lukas](https://github.com/luksen99) <img align="right" width="100" height="100" src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/SU.jpg">  <img align="right" width="100" height="100" src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/UDCM_logo.png">
* [Scalera Valentino](https://github.com/VaSca92) 
* [Vivek Unikandanunni](https://github.com/VivekUUnni)
* [UDCM Group of SU](http://udcm.fysik.su.se/)

You can directly contact us via mail: [Lukas Alber](mailto:lukas.alber@fysik.su.se)


#### Cite 
   <p align="center"> 
   <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3333493.svg" width="300" height="30" />   
   </p>



#### How to contribute : 

Since finding the thermophysical parameters needed to solve the equation under consideration (see above) is current subject to research, we want to encourage our users to contribute with their knowledge of parameters. That is, please **send us the data of the parameters in use, together with references from literature** ( [Lukas Alber](mailto:lukas.alber@fysik.su.se)). We are working on providing an open source data base. Also see [NTMpy](https://github.com/udcm-su/heat-diffusion-1D/edit/master/NTMpy/README.md) page.

Fork from the `Developer`- branch and pull request to merge back into the original `Developer`- branch. 
Working updates and improvements will then be merged into the `Master` branch, which will always contain the latest working version.


#### Dependencies:

[Numpy](http://www.numpy.org/)

[Matplotlib](https://matplotlib.org/)

[B-splines](https://github.com/johntfoster)

[Progressbar](https://pypi.org/project/tqdm/2.2.3/)

Note that by downloading the package via `pip install NTMpy` a routine, which automatically checks if all the required packages are existent on the local machine is implemented. If one of the dependent pip packages,listed here, is missing an automatic download is initiated.

  

       
  

 

# NTMpy - N-Temperature Model solver


#### To download                     `pip install NTMpy`

#### To call routines in a script:   `from NTMpy import NTMpy as ntm`

#### To update to lates version:     `pip install --upgrade NTMpy`

Information and citation refere to [NTMpy: An open source package for solving coupled parabolic differential equations in the framework of the three-temperature model](https://arxiv.org/abs/2002.04559) paper.

Further information on the solver package itself can be found here: [NTMpy](https://github.com/udcm-su/heat-diffusion-1D/tree/master/NTMpy).

Documentation and example sessions can be found in the [Wiki](https://github.com/udcm-su/heat-diffusion-1D/wiki).

------------------------------------------------------------------------------------------------------------------

This is a code providing a solution to the heat diffusion in a 1D structure in a 3-temperature model approximation.

We consider:
* Material comprising piecewise homogeneous layers
* Heating of electron system with an energy source with Gaussian or custom shape (i.e. an ultrashort laser pulse)
* Three temperature model: electron, lattice and spin systems are considered and coupled
* Transfer Matrix Method to account for local absorbed energy in the multi layer material

The equation under consideration is: 

 <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/master/Pictures/Equation.PNG" width="750" height="120" />
 
 where *C* = specific heat, *k* = conductivity, *G* is the coupling constant between the two systems (Electron and Lattice)
  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiE.png" width="30" height="20" /> and <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiL.png" width="30" height="20" /> 
  are the respective temperatures of the electron lattice abd spin system with respect to space *x* and time *t*. The superscripts  *L* and *E* indicate whether a parameter belongs to the electron (E) lattice (L) or spin (S) system and the sub index *i* denotes to which layer the parameter belongs.
  
  A sketch of the model: 
  
  <p align="center">
  <img src="https://github.com/udcm-su/NTMpy/blob/master/Pictures/ThreeTMscetch.png" width="420" height="320" />
  </p>
 
 Our approach is to use a combination of **finite elements** (B-Splines) to approximate the derivation in space and **Explicit Euler** to approximate the evolution in time.
 To stay within the **stability** region of the Explicit Euler algorithm, a well suited time step is automatically computed.
 
  ### Example:
Using the material parameters, reported in [Ultrafast Spin Dynamics in Ferromagnetic Nickel](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.76.4250) we consider a 3- temperature system in a Nickel layer. We are able to simulate the reported temperatuer map for all systems, in accordance with the paper.
 
  Temperature evolution of 3-coupled systems |  Gaussian laser pulse S(x,t) hitting sample
:-------------------------:|:-------------------------:
 <img src="https://media.giphy.com/media/RHE9DS2kPSdobin3hv/giphy.gif" width="320" height="300" />  |  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/Source.png" width="320" height="300" />
 
Here we consider an energy source *S(x,t)*, that is Gaussian in time and decays exponentially, representing a laser pulse hitting the material from the left. Note, that the implemented Transfer Matrix Method makes it easily able to simulate realistic laser sources, taking multiple reflection, different incident angles and polarization into account.
 
The animation clearly shows the following effects: 
  1. The electron system heats up immediately, as expected looking at the equation above.
  2. The heat in the electronic system diffuses along the x-axis.
  3. The heat transfers to the subsystems on longer time scales.

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

`@misc{alber2020ntmpy,
    title={NTMpy: An open source package for solving coupled parabolic differential equations in the framework of the three-temperature model},
    author={Lukas Alber and Valentino Scalera and Vivek Unikandanunni and Daniel Schick and Stefano Bonetti},
    year={2020},
    eprint={2002.04559},
    archivePrefix={arXiv},
    primaryClass={cond-mat.other}
}`



#### How to contribute : 

Since finding the thermophysical parameters needed to solve the equation under consideration (see above) is currently subject to research, we want to encourage our users to contribute with their knowledge of parameters. That is, please **send us the data of the parameters in use, together with references from literature** ( [Lukas Alber](mailto:lukas.alber@fysik.su.se)). We are working on providing an open source data base. Also see [NTMpy](https://github.com/udcm-su/heat-diffusion-1D/edit/master/NTMpy/README.md) page.

Fork from the `Developer`- branch and pull request to merge back into the original `Developer`- branch. 
Working updates and improvements will then be merged into the `Master` branch, which will always contain the latest working version.


#### Dependencies:

[Numpy](http://www.numpy.org/)

[Matplotlib](https://matplotlib.org/)

[B-splines](https://github.com/johntfoster)

[Progressbar](https://pypi.org/project/tqdm/2.2.3/)

Note that by downloading the package via `pip install NTMpy` a routine, which automatically checks if all the required packages are existent on the local machine is implemented. If one of the dependent pip packages, listed here, is missing an automatic download is initiated.

  

       
  

 

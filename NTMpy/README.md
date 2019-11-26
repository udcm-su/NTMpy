# NTMpy solver

### Introduction
The NTMpy software combines routines which allow the user to solve parabolic differential equations. After being subject to an external stimulus, e.g. a laser pulse, the diffusion dynamics of e.g. heat are modeled for up to three coupled systems in multiple piece wise homogeneous one dimensional layers. 
In order to provide a user friendly interface, the software is written in an object oriented way, allowing it to customize the input parameters in an easy to handle manner, such that material specific parameters as well as boundary- , initial conditions and the source can be adjusted for every simulation. To make data visualization and analysis more easy, ready made plotting routines are as well provided to the user. In this paper the mathematical background of the software, focusing on the implementation of the finite elements used to describe the solution in space, is presented. Further more we lay out the work flow, and the most important commands of the program. To show that this program really is suited to accurately model real life events, we demonstrate the physical relevance through a comparison of an experiment carried out in the lab and an independent simulation. 

### How to use
The solver is designed in the follwing way: 
**Source** --> **Simulation** --> **Visualization**

### Suggestions for further reading: 
* Non-equilibrium transient thermal grating relaxation in metal by K. Nelson et al. (Journal of Applied Physics 109)
* UDKM1DSIM - A simulation toolkit for 1D ultrafast dynamics in condensed matter by D. Schick et al. (Computer Physics Communications 185)
*    Multilayer optical calculations by Steven J. Byrnes (arXiv:1603.02720v3)
*    Electron and lattice dynamics following optical excitation of metals by j. Hohlfeld et al. (Chemical Physics 251)
*    Electron emission from metal surfaces exposed to ultrashort laser pulses by S. I. Anisimov et al. (Sov. Phys.-JETP Vol. 39 No. 2)
*    Explaining the paradoxical diversity of ultrafast laser- induced dematnetization (B. Koopmans et al. )


#### How to contribute : 

Since finding the thermophysical parameters needed to solve the equation under consideration is current subject to research, we want to encourage our users to contribute with their knowledge of parameters. That is, please **send us the data of the parameters in use, together with references from literature** ( [Lukas Alber](mailto:lukas.alber@fysik.su.se)). We are working on providing an open source data base. Also see [NTMpy](https://github.com/udcm-su/heat-diffusion-1D/edit/master/NTMpy/README.md) page.


#### Cite 
   <p align="center"> 
   <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3333493.svg" width="300" height="30" />   
   </p>

[Main page](https://github.com/udcm-su/heat-diffusion-1D)

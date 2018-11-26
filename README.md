# Heat-Diffusion-1D
This is a software providing a solution to the heat diffusion in a 1D structure. 

We consider:
* Multiple piecwise homogenious materials behind each other. 
* Electron system heating with a Gaussian like source
* Two temperature model: electron & lattice temperature

The equation under consideration is: 

 <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/DiffusionEq.png" width="700" height="100" />
 
 where *C* = heat capacity, *k* = conductivity, *G* is the coupling constant between the two systems (Electron and Lattice)
  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiE.png" width="30" height="20" /> and <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/phiL.png" width="30" height="20" /> 
  are the respective temperatures of the electron and the lattice system with respect to space *x* and time *t* .

 Our approach is to use a combination of **finite elements** (B-Splines) to approximate the derivation in space and **Explicit Euler** to approximate the evolution in time.
 To stay within the **stability** region of the Explicit Euler algorithm, a well suited time step is automatically computed.
 
  ### Example showing our output:
 Temperature evolution in time and space after a laser pulse hits the probe
 
  Temperature evolution of electron- lattice system |  Gaussian laser pulse S(x,t) hitting probe
:-------------------------:|:-------------------------:
 <img src="https://media.giphy.com/media/dIUAz7xfof5N8B8tUy/giphy.gif" width="320" height="300" />  |  <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/Source.png" width="320" height="300" />
 
 Here we consider a gaussian in time and exponentially decaying in space laser source *S(x,t)* = *lambda* * *exp(-lambda* * *(x-x0))* * *G(t-t0)* hitting a material from the left. 
 
 First the electron system emediately gets heated up. 
 Then we can see two effects: 
 
  1. Diffusion of the electron heating along the x-axis
  2. Heat gets transported to the lattice system
  
  
   | The material under consideration is  | 
   | ------------------------------------ | 
   | <img src="https://github.com/udcm-su/heat-diffusion-1D/blob/Developer/Pictures/SroSto.PNG" width="320" height="200" />     |   
       
  
  .
 
 

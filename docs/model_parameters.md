# Model parameters

The joint inversion program, `joint_inv`,  evaluates the posterior probability distributions of model parameters. The model parameters are:
 
* Number of layers (<img src="https://latex.codecogs.com/gif.latex?k" title="k" />)
* Depth of each layer bottom (<img src="https://latex.codecogs.com/gif.latex?z_i,&space;i=1&space;\cdots&space;k" title="z_i, i=1 \cdots k" />)
* S-wave velocity of each layer (<img src="https://latex.codecogs.com/gif.latex?\beta_i,&space;i=1&space;\cdots&space;k" title="\beta_i, i=1 \cdots k" />)    
* P-wave velocity of each layer (optional; <img src="https://latex.codecogs.com/gif.latex?\alpha_i,&space;i=1&space;\cdots&space;k" title="\alpha_i, i=1 \cdots k" />)
* Standard deviation of data noise (optional; <img src="https://latex.codecogs.com/gif.latex?\sigma_j,&space;j=1,&space;\cdots,&space;N_{data}" title="\sigma_j, j=1, \cdots, N_{data}" />)

Alternatively, one can solve for velocity anomalies centered on the given reference velocity, as is done by [Akuhara et al. (2020)](https://doi.org/10.1029/2020GL088280). This mode is activated by setting [solve_anomaly](parameter_list.md#solve_anomaly) to be 'T' in [the main parameter file](https://github.com/akuhara/SEIS_FILO/wiki/Main-Parameter-File). 


### Note 
* The ocean and bottom half-space layers are not counted for <img src="https://latex.codecogs.com/gif.latex?k" title="k" />.
* The properties of the ocean layer are fixed during the inversion (P-wave velocity: 1.5 km/s and density: 1.0 g/cm^3, by default). These default values can be changed through the parameters [vp_ocean](parameter_list.md#vp_ocean) and [rho_ocean](parameter_list.md#rho_ocean).
* The properties of the bottom half-sapce are fixed during the inversion (P-wave velocity: 8.1 km/s, S-wave velocity: 4.6 km/s and density: 3.3 g/cm^3, by default). These default values can be changed through the parameters [vp_bottom](parameter_list.md#vp_bottom), [vs_bottom](parameter_list.md#vs_bottom), and [rho_bottom](parameter_list.md#rho_bottom). If you want to solve for the velocities of this layer, set [vp_bottom](parameter_list.md#vp_bottom) and [vs_bottom](parameter_list.md#vs_bottom) to be negative (__NOT RECOMMENDED WHEN INVERTING DISPERSION CURVES__ because low velocity at the half space often results in the failure of dispersion computation.)
* The density of is calculated from the P-wave velocity using the empirical relationship by [Brocher (2005)](
https://doi.org/10.1785/0120050077): 
    * <img src="https://latex.codecogs.com/gif.latex?\rho=1.6612*V_p&space;-0.4721*V_p^2&space;&plus;&space;0.0671*V_p^3&space;-&space;0.0043*V_p^4&space;&plus;&space;0.000106*V_p^5" title="\rho=1.6612*V_p -0.4721*V_p^2 + 0.0671*V_p^3 - 0.0043*V_p^4 + 0.000106*V_p^5">

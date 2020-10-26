# Model parameters

The joint inversion program, `joint_inv`,  samples the posterior probability distributions of model parameters regarding 1-D velocity structure. The model parameters are:
 
* Number of layers: \\(k\\)
* Depth of each layer bottom: \\(z_i (i=1,\cdots,k)\\)
* S-wave velocity of each layer: \\(\beta_i (i=1,\cdots,k)\\)
* P-wave velocity of each layer (optional): \\(\alpha_i (i=1,\cdots,k)\\)
* Standard deviation of data noise (optional): \\(\sigma_i (i=1,\cdots,N_{data})\\)

Alternatively, one can solve for velocity anomalies centered on the given reference velocity, as is done by [Akuhara et al. (2020)](https://doi.org/10.1029/2020GL088280). This mode is activated by setting [solve_anomaly](parameter_list.md#solve_anomaly) to be 'T' in [the main parameter file](main_parameter_file.md). 


!!! NOTE
    * The ocean and bottom half-space layers are not counted for \\(k\\).
    * The properties of the ocean layer are fixed during the inversion (P-wave velocity: 1.5 km/s and density: 1.0 g/cm^3, by default). These default values can be changed through the parameters [vp_ocean](parameter_list.md#vp_ocean) and [rho_ocean](parameter_list.md#rho_ocean).
    * The properties of the bottom half-sapce are also fixed during the inversion (P-wave velocity: 8.1 km/s, S-wave velocity: 4.6 km/s and density: 3.3 g/cm^3, by default). These default values can be changed through the parameters [vp_bottom](parameter_list.md#vp_bottom), [vs_bottom](parameter_list.md#vs_bottom), and [rho_bottom](parameter_list.md#rho_bottom). If you want to solve for the velocities of this layer, set [vp_bottom](parameter_list.md#vp_bottom) and [vs_bottom](parameter_list.md#vs_bottom) to be negative (__NOT RECOMMENDED WHEN INVERTING DISPERSION CURVES__ because low velocity at the half space often results in the failure of dispersion computation.)
    * The density of is calculated from the P-wave velocity using the empirical relationship by [Brocher (2005)](https://doi.org/10.1785/0120050077): \\( \rho = 1\.6612V_p-0\.4721V_p^2+0\.0671V_p^3-0\.0043V_p^4+0\.000106V_p^5\\).

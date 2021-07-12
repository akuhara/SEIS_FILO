# Parameter list

SEIS_FILO's programs involve a number of tuning parameters. The complete list is given here. These parameters should appear in the  [main parameter file](main_parameter_file.md) or [data summary file](data_summary_file.md). 

## Surface wave 

### freq_or_period

Input type of dispersion measurements. This parameter must be either 'freq' or 'period' (freq: dispersion curve as a function of frequency; period: as a functon of period).

* Type: Character(6)

### xmin

The minimum frequency or period of a dispersion curve (Hz or s).

* Type: Double precision

### xmax

The maximum frequency or period of a dispersion curve (Hz or s).

* Type: Double precision

### dx

The frequency or period interval of a dispersion curve (Hz or s).

* Type: Double precision

### cmin

The minimum phase velocity (km/s). This value is used as the starting point of the root search for dispersion curves.

* Type: Double precision

### cmax

The maximum phase velocity (km/s).

* Type: Double precision

### dc

The phase velocity interval (km/s) used for the root search. This interval is used in the first-step rough search for zero crossings. The second-step binary search follows this rough search.

* Type: Double precision


### disp_phase

The phase-type of surface waves. L: Love-wave wave (not supported yet), R: Rayleigh wave.

* Type: Character (length=1)

### n_mode

The number of mode for surface waves. 0: fundamental mode, 1: 1st overtone, ...

* Type: Integer

!!! NOTE
    If this parameter is set to be less than 0, full search mode is activated. See [here](dispersion_curve.md) for more details. 

### n_disp

The number of dipsersuion curves used for inversion. Need to specify in the first line of the [data summary file](data_summary_file.md). 

* Type: Integer

### sig_c_min

The minimum standard deviation for data noise in phase velocity (km/s). Only used for inversion.

* Type: Double precision

### sig_c_max

The maximum standard deviation for data noise in phase velocity (km/s). Only used for inversion.

* Type: Double precision

### dev_sig_c

The standard deviation of random perturbation amount added to the standard deviation of phase velocity noise (km/s). Only used for inversion.

* Type: Double precision

### sig_u_min

The minimum standard deviation for data noise in group velocity (km/s). Only used for inversion.

* Type: Double precision

### sig_u_max

The maximum standard deviation for data noise in group velocity (km/s). Only used for inversion.

* Type: Double precision

### dev_sig_u

The standard deviation of random perturbation amount added to the standard deviation of group velocity noise (km/s). Only used for inversion.

* Type: Double precision

### disper_out

The output filename for _disper_fwd_.

* Type: Character (length=200)

## Receiver function

### n_smp

The number of elements in receiver functions.

* Type: Integer

### delta

The sampling interval of receiver functions (sec).

* Type: Double precision

### rayp

Ray parameter (s*km^-1).

* Type: Double precision

### a_gauss

The parameter for the Gaussian low-pass filter.

* Type: Double precision

### t_pre

Time length (sec) preceding the direct P arrival on receiver functions.

* Type: Double precision

### damp

Damping factor for water-level deconvolution. Default is 0.

* Type: Double precision

### rf_phase

The phase-type of receiver functions: P (P-wave); S (S-wave).

* Type: Characeter (length=1)

### deconv_flag

Whether the forward computation accompanies deconvolution (T) or not (F). For deconv_flag=T, receiver functions are computed by the deconvolution of the radial-component deconvolved by the vertical component. For deconv_flag=F, the output is the radial-component impulse response. 

* Type: Logical (T or F)

### correct_amp

Whether the observed receiver functions (input data) is corrected for the energy loss by the Gaussian low-pass filter.
See [here](receiver_function_amplitude.md) for more details.

* Type: Logical (T or F)

### recv_func_out

The output SAC filename for _recv_func_fwd_.

* Type: Character (length=200)

### n_rf

The number of receiver functions use for inverions. Need to specify in the first line of the [data summary file](data_summary_file.md).

* Type: Integer

### sig_rf_min

The minimum standard deviation for data noise for receiver function amplitude. Only used for inversion.

* Type: Double precision

### sig_rf_max

The maximum standard deviation for data noise for receiver function amplitude. Only used for inversion.

* Type: Double precision

### dev_sig_rf

The standard deviation of random perturbation amount added to the standard deviation of receiver function noise. Only used for inversion.

* Type: Double precision


## Forward computation 

### vmod_in

Input filename that contains a velocity model for forward computations.

* Type: Character (length=200)

### noise_added

Standard deviation of noise added to forward computation results. 

* Type: Double precision 

## Inversion

### n_iter

The number of total iterations for the inversion, which includes the burn-in period.

* Type: Integer

### n_burn 

The number of iterations in the burn-in period. The sampled model is not saved for this period to forget the initial state. 

* Type: Integer

### n_corr

The interation interval at whith the sampled model is saved. This is useful to avoid artificial correlation to the previous sample. 

* Type: Integer

### n_proc

The number of processes. __NOTE:__ this parameter should be passed as command line arugment: `mpirun -np [n_proc] bin/joint_int joint_inv.in`, for example. 

!!! Warning
    DONOT put this parameter in any input files.

* Type: Integer

### n_chain

The number of MCMC chains per process.

* Type: Integer

### n_cool

The number of non-tempered MCMC chains per process.

* Type: Integer

### temp_high

The maximum temperature allowed for parallel tempering scheme. 

* Type: Double precision

### is_ocean

Whether the velocity structure has the sea (T) or not(F).

* Type: Logical (T or F)

### ocean_thick

The thickness of the ocean layer (km)

* Type: Double precision 

### k_min

The minimum number of layers (excluding the seawater and bottom half-space).

* Type: Integer

### k_max

The maximum number of layer interfaces (excluding the seawater and bottom half-space).

* Type: Integer

### z_min

The minimum depth of layer interfaces (km below the sea surface). __NOTE:__ z_min must be deeper than the sea bottom for ocean setting.

* Type: Double precision

### z_max

The maximum depth of layer interfaces (km below the sea surface).

* Type: Double precision

### solve_vp

Whether Vp is solved for (T) or not (F).

* Type: Logical (T or F).

### vp_min

The minimum Vp (km/s)

* Type: Double precision

### vp_max

The maximum Vp (km/s)

* Type: Double precision

### vs_min

The minimum Vs (km/s).

* Type: Double precision

### vs_max

The maximum Vs (km/s).

* Type: Double precision

### vp_ocean

Vp of the ocean layer (km/s). The default value is 1.5 km/s.

* Type: Double precision

### rho_ocean

The density of the ocean layer (g/cm^3). The default value is 1.0 g/cm^3.

* Type: Double precision

### vp_bottom

Vp of the bottom half-space layer (km/s). The default value is 8.1 km/s. Set this parameter to be negative (and solve_vp to be .true.) if you want to solve for the Vp of this layer.

* Type: Double precision

### vs_bottom

Vs of the bottom half-space layer (km/s). The default value is 4.6 km/s. Set this parameter to be negative if you want to solve for the Vs of this layer.

* Type: Double precision

### rho_bottom

The density of the bottom half-space layer (km/s). The default value is 4.6 km/s. If a negative value is given to this parameter, the density of the bottom layer is determined empirically from its Vp.   

* Type: Double precision

### solve_anomaly

If this parameter is set to ".true." (or "T"), velocity anomalies are solved for instead of the absolute  values. These anomalies are defined relative to a reference velocity model given by users. 

* Type: Logical

### dvp_sig

The standard deviation of Gaussian prior probability for Vp anomaly.
This parameter is valid only if [solve_anomaly](parameter_list.md#solve_anomaly) is true. 

* Type: double precision 

### dvs_sig

The standard deviation of Gaussian prior probability for Vs anomaly.
This parameter is valid only if [solve_anomaly](parameter_list.md#solve_anomaly)
is true. 

* Type: double precision 

### ref_vmod_in

A filename that contains a reference velocity model. This parameter is valid only if [solve_anomaly](parameter_list.md#solve_anomaly) is true. 

* Type: Character (length=200)


### dev_z

The standard deviation of random perturbation added to layer interfaces (km)

* Type: Double precision

### dev_vp

The standard deviation of random perturbation added to P-wave velocity (km/s)

* Type: Double precision

### dev_vs

The standard deviation of random perturbation added to S-wave velocity (km/s)

* Type: Double precision

### n_bin_z

The number of depth bins to construct the posterior probability.

* Type: Integer

### n_bin_vp

The number of Vp bins to construct the posterior probability.

* Type: Integer

### n_bin_vs

The number of Vs bins to construct the posterior probability.

* Type: Integer

### disper_in

The filename of the [data summary file](data_summary_file.md) for dispersion curves.

* Type: Character (length=200)

### recv_func_in

The filename of the [data summary file](data_summary_file.md) for receiver functions.

* Type: Character (length=200)

### diagnostic_mode

If .true. (or T), it activates diagnostic mode (see [diagnostic_mode](diagnostic_mode.md)). Default is .false..

* Type: logical


## Common for all

### is_sphere

Given (or output) velocity model is for spherical Earth (T) or not (F). When it is true, the earth-flattening transformation is applied to the given velocity model before forward computation. Consider setting this parameter to be "T" for studies on Mantle structure (>~ 100 km depth).

* Type: logical

### r_earth

The earth's radius in kilo-meter used for the earth-flattening transformation. The default value is 6371.

* Type: double precision  

### i_seed1

The number to initialize pseudo-random numbers. The Xorshift algorithm is used for the random number generation. 

* Type: Integer

### i_seed2

The number to initialize pseudo-random numbers. The Xorshift algorithm is used for the random number generation. 

* Type: Integer

### i_seed3

The number to initialize pseudo-random numbers. The Xorshift algorithm is used for the random number generation. 

* Type: Integer

### i_seed4

The number to initialize pseudo-random numbers. The Xorshift algorithm is used for the random number generation. 

* Type: Integer

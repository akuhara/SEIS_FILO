# Main parameter file

The main parameter file contains user's choice on tuning parameters. This filename must be passed to SEIS_FILO's programs as the **first command-line argument**. See [parameter list](parameter_list.md) for descriptions of each parameter.  

__EXAMPLE__

```
recv_func_fwd [path to parameter file]
disper_fwd [path to parameter file]
mpirun -np 20 joint_inv [path to parameter file]
```

## Format
* Specify one parameter per line
* Parameter's name should be on the left side of "=" and the value on the right side
* Comment out by "#"
* Order insensitive
* Empty lines are just ignored
```
# Snippet of the parameter file
n_iter = 100000
n_chain = 20
# The following lines give the same results.
n_burn = 50000
n_ bu rn= 500 0 0
```

## Required parameters for _disper_fwd_

* [freq_or_period](parameter_list.md#freq_or_period): dispersion measurement type (freq: dispersion curve as a function of frequency; period: as a function of period)
* [xmin](parameter_list.md#xmin): minimum frequecny or period of dispersion curve (Hz or s)
* [xmax](parameter_list.md#xmax): maximum frequency or period of dispersion curve (Hz or s)
* [dx](parameter_list.md#dx): freqeuncy or period interval of dispersion curve (Hz or s)
* [cmin](parameter_list.md#cmin): minimum phase velocity for root search (km/s)
* [cmax](parameter_list.md#cmax): maximum phase velocity for root search (km/s)
* [dc](parameter_list.md#dc): inverval of phase velocity (km/s)
* [disper_phase](parameter_list.md#disper_phase): phase type (R: Rayleigh, L: Love)
* [n_mode](parameter_list.md#n_mode): mode number (n_mode >= 0)
* [vmod_in](parameter_list.md#vmod_in): filename for input velocity model
* [disper_out](parameter_list.md#disper_out): output filename


!!! Note
    * disper_phase = L (Love wave) is nut supported for current version.
    * When n_mode < 0, [full search mode](dispersion_curve.md) is acctivated

## Required parameters for _recv_func_fwd_

* [n_smp](parameter_list.md#n_smp): number of samples in receiver function data
* [delta](parameter_list.md#delta): time interval of receiver function (s)
* [rayp](parameter_list.md#rayp): ray parameter (s/km)
* [a_gauss](parameter_list.md#a_gauss): Gaussian low-pass filter parameter
* [t_pre](parameter_list.md#t_pre): time before direct P arrival in receiver function (s)
* [rf_phase](parameter_list.md#rf_phase): phase type (P: P-wave, S: S-wave)
* [deconv_flag](parameter_list.md#deconv_flag): whether deconvolved by Z component or not (T: with deconvolution, F: without deconvolution)
* [correct_amp](parameter_list.md#correct_amp): amplitude correction to account for energy loss by filtering (T: with correction, F: without correction)
* [vmod_in](parameter_list.md#vmod_in): filename for input velocity model
* [recv_func_out](parameter_list.md#recv_func_out): output filename

## Required parameters for _joint_inv_

* [n_iter](parameter_list.md#n_iter): number of iterations
* [n_burn](parameter_list.md#n_burn): number of iterations for burn-in period
* [n_corr](parameter_list.md#n_corr): number of iterations each time before saving a model
* [n_chain](parameter_list.md#n_chain): number of MCMC chains per process
* [n_cool](parameter_list.md#n_cool): number of non-tempered MCMC chains per process
* [temp_high](parameter_list.md#temp_high): maximum temperature
* [is_ocean](parameter_list.md#is_ocean): whether model includes ocean (T: with ocean layer, F: without ocean layer)
* [ocean_thick](parameter_list.md#ocean_thick): ocean layer thickness (only required when is_ocean=T)
* [k_min](parameter_list.md#k_min): minimum number of layer interfaces
* [k_max](parameter_list.md#k_max): maximum number of layer interfaces
* [z_min](parameter_list.md#k_min): minimum interface depth (km)
* [z_max](parameter_list.md#k_max): maximum interface depth (km)
* [solve_vp](parameter_list.md#solve_vp): whether solve for Vp (T: solve Vp, F: not)
* [vp_min](parameter_list.md#vp_min): minimum Vp (km/s) (only required when solve_vp=T)
* [vp_max](parameter_list.md#vp_max): maximum Vp (km/s) (only required when solve_vp=T)
* [vs_min](parameter_list.md#vs_min): minimum Vs (km/s).
* [vs_max](parameter_list.md#vs_max): maximum Vs (km/s).
* [dev_z](parameter_list.md#dev_z): standard deviation of proposed layer depth perturbation (km).
* [dev_vp](parameter_list.md#dev_vp): standard deviation of proposed Vp perturbation (only required when solve_vp=T) (km/s).
* [dev_vs](parameter_list.md#dev_vs): standard deviation of proposed Vs perturbation (km/s).
* [n_bin_z](parameter_list.md#n_bin_z): number of depth bin for display.
* [n_bin_vp](parameter_list.md#n_bin_vp): number of Vp bin for display (only required when solve_vp=T).
* [n_bin_vs](parameter_list.md#n_bin_vs): number of Vs bin for display.
* [i_seed1](parameter_list.md#i_seed1): random seed 1
* [i_seed2](parameter_list.md#i_seed2): random seed 2
* [i_seed3](parameter_list.md#i_seed3): random seed 3
* [i_seed4](parameter_list.md#i_seed4): random seed 4









 
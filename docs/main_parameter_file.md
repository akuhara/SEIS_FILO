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

* [fmin](parameter_list.md#fmin)
* [fmax](parameter_list.md#fmax)
* [df](parameter_list.md#df)
* [cmin](parameter_list.md#cmin)
* [cmax](parameter_list.md#cmax)
* [dc](parameter_list.md#dc)
* [disper_phase](parameter_list.md#disper_phase)
* [n_mode](parameter_list.md#n_mode)
* [vmod_in](parameter_list.md#vmod_in)
* [disper_out](parameter_list.md#disper_out)

## Required parameters for _recv_func_fwd_

* [n_smp](parameter_list.md#n_smp)
* [delta](parameter_list.md#delta)
* [rayp](parameter_list.md#rayp)
* [a_gauss](parameter_list.md#a_gauss)
* [t_pre](parameter_list.md#t_pre)
* [rf_phase](parameter_list.md#rf_phase)
* [deconv_flag](parameter_list.md#deconv_flag)
* [correct_amp](parameter_list.md#correct_amp)
* [vmod_in](parameter_list.md#vmod_in)
* [recv_func_out](parameter_list.md#recv_func_out)

## Required parameters for _joint_inv_

* [n_iter](parameter_list.md#n_iter)
* [n_burn](parameter_list.md#n_burn)
* [n_corr](parameter_list.md#n_corr)
* [n_chain](parameter_list.md#n_chain)
* [n_cool](parameter_list.md#n_cool)
* [temp_high](parameter_list.md#temp_high)
* [is_ocean](parameter_list.md#is_ocean)
* [ocean_thick](parameter_list.md#ocean_thick) (only required when is_ocean=T)
* [k_min](parameter_list.md#k_min)
* [k_max](parameter_list.md#k_max)
* [z_min](parameter_list.md#k_min)
* [z_max](parameter_list.md#k_max)
* [solve_vp](parameter_list.md#solve_vp)
* [vp_min](parameter_list.md#vp_min) (only required when solve_vp=T)
* [vp_max](parameter_list.md#vp_max) (only required when solve_vp=T)
* [vs_min](parameter_list.md#vs_min)
* [vs_max](parameter_list.md#vs_max)
* [dev_z](parameter_list.md#dev_z)
* [dev_vp](parameter_list.md#dev_vp) (only required when solve_vp=T)
* [dev_vs](parameter_list.md#dev_vs)
* [n_bin_z](parameter_list.md#n_bin_z)
* [n_bin_vp](parameter_list.md#n_bin_vp) (only required when solve_vp=T)
* [n_bin_vs](parameter_list.md#n_bin_vs)
* [i_seed1](parameter_list.md#i_seed1)
* [i_seed2](parameter_list.md#i_seed2)
* [i_seed3](parameter_list.md#i_seed3)
* [i_seed4](parameter_list.md#i_seed4)









 
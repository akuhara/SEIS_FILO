# SEIS_FILO 
SEISmological transdimensional inversion tools for Flat and Isotropic Layered medium in the Ocean

## Package contents

* __recv_func_fwd:__ Synthetic calculation of receiver functions

* __disper_fwd:__ Synthetic calculation of surface wave dispersion curves

* __joint_inv:__ Transdimensional joint inversion of receiver functions and surf
ace wave dispersion curves

## Before use

* [Overview](overview.md)

* [Requirements](requirements.md)

* [Install](install.md)

* [Temrs of use](terms_of_use.md)


## How to use

* Inputs
    * [Overview of input files](overview_of_input_files.md)
    * [Main parameter file (required for `disper_fwd`, `recv_func_fwd`, and `joint_inv`)](main_parameter_file.md)
    * [Data summary file (required for `joint_inv`)](data_summary_file.md)
    * [Data file (required for `joint_inv`)](data_file.md)
    * [Reference velocity file (option for `joint_inv`)](reference_velocity_file.md)
    * [Velocity model file (required for `disper_fwd` and `recv_func_fwd`)](velocity_model_file.md)
    * [Parameter list](parameter_list.md)	

* Outputs
    * [Number of layers](number_of_layers.md)
    * [Velocity profile](velocity_profile.md)
    * [Dispersion curve ensemble](dispersion_curve_ensemble.md)
    * [Receiver function ensemble](receiver_function_ensemble.md)
    * [Standard deviation of noise](standard_deviation_of_noise.md)

* [Plot utilities](plot_utilities.md)

* Forward calculation
    * [Dispersion curve](dispersion_curve.md)
    * [Receiver function](receiver_function.md)

* Inversion analysis
    * [Model parameters](model_parameters.md)
    * [Prior probability](prior_probability.md)
    * [MCMC step](mcmc_step.md)
    * [Parallel tempering](parallel_tempering.md)

* Miscellaneous
    * [Receiver function amplitude](receiver_function_amplitude.md)
    * [Verification](verification.md)
    * [Spherical Earth mode](spherical_earth_mode.md)
    * [Diagnostic mode](diagnostic_mode.md)
    * [Docker container](docker_container.md)
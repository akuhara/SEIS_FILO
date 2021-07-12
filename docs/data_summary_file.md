# Data summary file


The data summary file passes the summary of input data to the inversion program `joint_inv`. Different file formats are used for surface wave dispersion curves and receiver functions. The filename can be arbitrary, which are to be specified in [the parameter file](https://github.com/akuhara/SEIS_FILO/wiki/Parameter-File) via the parameters [disper_in](parameter_list.md#disper_in) and [recv_func_in](parameter_list.md#recv_func_in).

## Surface wave dispersion curves 

__Format__

* Comment out by # is acceptable

* Empty lines are ignored.

* Line 1: [n_disp](parameter_list.md#n_disp) (The number of data files)

* Line 2 and after: Summary of each data file (must be repeated n_disp times)
    * Line i: data file name
    
    * Line ii: [disper_phase](parameter_list.md#disper_phase), [n_mode](parameter_list.md#n_mode), [freq_or_period](parameter_list.md#freq_or_period)
    
    * Line iii: [nx](parameter_list.md#nx), [xmin](parameter_list.md#xmin), [dx](parameter_list.md#dx)

    * Line iv: [cmin](parameter_list.md#cmin), [cmax](parameter_list.md#cmax), [dc](parameter_list.md#dc)

    * Line v: [sig_c_min](parameter_list.md#sig_c_min), [sig_c_max](parameter_list.md#sig_c_max), [dev_sig_c](parameter_list.md#dev_sig_c)

    * Line vi: [sig_u_min](parameter_list.md#sig_u_min), [sig_u_max](parameter_list.md#sig_u_max), [dev_sig_u](parameter_list.md#dev_sig_u) 

* The example is [here](https://github.com/akuhara/SEIS_FILO/blob/main/sample/joint_inv/disper.in).

 

## Receiver functions

__Format__

* Comment out by # is acceptable

* Empty lines are ignored.

* Line 1: [n_rf](parameter_list.md#n_rf) (The number of data files)

* Line 2 and after: Summary of each data file (must be repeated n_rf times)
    * Line i: data file name

    * Line ii: [a_gauss](parameter_list.md#a_gauss), [rayp](parameter_list.md#rayp), [delta](parameter_list.md#delta)

    * Line iii: [rf_phase](parameter_list.md#rf_phase)

    * Line iv: t_start (start time of analysis window), t_end (end time of the analysis window)

    * Line v: [sig_rf_min](parameter_list.md#sig_rf_min), [sig_rf_max](parameter_list.md#sig_rf_max), [dev_sig_rf](parameter_list.md#dev_sig_rf) 

    * Line vi: [deconv_flag](parameter_list.md#deconv_flag), [correct_amp](parameter_list.md#correct_amp)

* The example is [here](https://github.com/akuhara/SEIS_FILO/blob/main/sample/joint_inv/recv_func.in).
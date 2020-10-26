# Reference velocity file

When solving for velocity anomalies ([solve_anomaly](parameter_list.md#solve_anomaly) = T), you must pass a reference velocity model through this file. The filename can be arbitrary but needs to be specified in a [main parameter file](main_parameter_file.md) using the parameter [solve_anomaly](parameter_list.md#solve_anomaly).

## Format
|1st column |2nd column |3rd column |
|:---|:---|:---|
|Depth (km) |Vp (km/s) |Vs (km/s) |

An example is [here](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/case_2_rec_func_with_ref_vmod/ref_vmod.in). 

!!! Note 
    * DEPTH INCREMENTS MUST BE CONSTANT.
    * The part shallower than the seafloor is ignored. The inversion fixes the Vp of the ocean layer at value of the parameter [vp_ocean](parameter_list.md#vp_ocean) (1.5 km/s by default). 


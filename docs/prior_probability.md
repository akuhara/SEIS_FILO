# Prior probability

User can adjust lower and upper bounds of uniform prior probabilities via the parameters listed below.

|Model parameter | Lower bound|Upper bound|
|:---|:---|:---|
|# of layer|[k_min](parameter_list.md#k_min)| [k_max](parameter_list.md#k_max)|
|Layer bottom depth|[z_min](parameter_list.md#z_min)| [z_max](parameter_list.md#z_max)|
|Vs|[vs_min](parameter_list.md#vs_min)| [vs_max](parameter_list.md#vs_max)|
|Vp|[vp_min](parameter_list.md#vp_min)| [vp_max](parameter_list.md#vp_max)|
|Standard deviation of receiver function noise|[sig_rf_min](parameter_list.md#sig_rf_min)| [sig_rf_max](parameter_list.md#sig_rf_max)|
|Standard deviation of phase velocity noise|[sig_c_min](parameter_list.md#sig_c_min)| [sig_c_max](parameter_list.md#sig_c_max)|
|Standard deviation of group velocity noise|[sig_u_min](parameter_list.md#sig_u_min)| [sig_u_max](parameter_list.md#sig_u_max)|

If [solve_anomaly](parameter_list.md#solve_anomaly) is true, Gaussian-type priors with zero-mean are assumed for the velocity anomalies. The standard deviation can be adjusted by the parameters below.

|Model parameter | Standard deviation (km/s) |
|:---|:---|
|Vp anomaly|[dvp_sig](parameter_list.md#dvp_sig) |
|Vs anomaly|[dvs_sig](parameter_list.md#dvs_sig) |

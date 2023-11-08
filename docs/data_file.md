# Data file

The data file contains the input data used for inversion. Different formats are used for surface wave dispersion curves and receiver functions. The filename and summary of the data should appear in the [data summary files](https://github.com/akuhara/SEIS_FILO/wiki/Data-Summary-File). 

## Surface wave dispersion curves

__Format__

|1st col.              |2nd col.                   |3rd col.                |4th col.                    |5th col.                    |6th col.                | 7th col.                  | 8th col.               |
|:---------------------|:--------------------------|:-----------------------|:---------------------------|:--------------------------|:--------------------------|:--------------------------|:-----------------------|
|Phase velocity (km/s) |Data availability (T or F) | Group velocity (km/s)  |Data availability (T or F)  | H/V ratio                 |Data availability (T or F)   |Rayleigh admittance (km/GPa) |Data availability (T or F) |
 
 
 
!!! NOTE
    * Ascending order is assumed in terms of frequency/period, where the first line must contain measurements at [xmin](parameter_list.md#xmin)).
    * An example is available at [GitHub repository](https://github.com/akuhara/SEIS_FILO/blob/main/sample/joint_inv/disper_obs.txt).

## Receiver functions

* Use SAC data file format. See [the official guidance by IRIS](https://ds.iris.edu/files/sac-manual/manual/file_format.html). 



# Data file

The data file contains the input data for inversion. Two different formats are used for surface wave dispersion curves and receiver functions. The file name and detailed information about the data should appear in the [data summary files](https://github.com/akuhara/SEIS_FILO/wiki/Data-Summary-File). 

## Surface wave dispersion curves

__Format__

|1st col.              |2nd col.                   |3rd col.                |4th col.                    |
|:---------------------|:--------------------------|:-----------------------|:---------------------------|
|Phase velocity (km/s) |Data availability (T or F) | Group velocity (km/s)  |Data availability (T or F)  | 
 
 
 
!!! NOTE
    * Ascending order is assumed in terms of frequency (the first line must contain measurements at [fmin](parameter_list.md#fmin)).
    * Example is [here](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/case_1_recv_func_and_rayleigh/rayleigh.0th).

## Receiver functions

* Use [SAC data file format](https://ds.iris.edu/files/sac-manual/manual/file_format.html). 



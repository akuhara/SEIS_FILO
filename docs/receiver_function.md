# Receiver functions

## Synthetic Seismograms
Radial and vertical components of synthetic seismograms are calculated using the Propagator Matrix method of Thomson and Haskell. In the case of offshore settings, we adopt a boundary condition that requires zero traction on the seafloor and allows discontinuous horizontal displacement. 

## Deconvolution
We have two options for deconvolution:

### Water-level deconvolution ([deconv_flag](parameter_list.md#deconv_flag) = T)
Radial component records are deconvolved by vertical component records using water-level deconvolution. The resulting receiver functions are compared with observed waveforms to measure misfits.

### No deconvolution ([deconv_flag](parameter_list.md#deconv_flag) = F)
No deconvolution is carried out. Synthetic radial component seismograms are directly compared with observed waveforms. You may want to choose this option if you have Green's function estimations instead of receiver functions (e.g., [Kumar et al. 2010](https://doi.org/10.1111/j.1365-246X.2009.04469.x); [Akuhara et al. 2019](https://doi.org/10.1029/2018JB016499)).

## Filtering
A Gaussian low-pass filter is applied to the deconvolved seismograms: 

![codecogseqn](https://user-images.githubusercontent.com/31914302/53931636-677f7d00-40d9-11e9-8325-de8240db9c02.gif),

where _a_ is a tuning parameter, which is to be specified as the parameter [a_gauss](parameter_list.md#a_gauss).

## S receiver function
When computing for S wave incidence, resultant seismogram is subject to time and sign reversals, following the convention for S receiver functions. 

## Amplitude normalization
We have two options for normalization:

### With amplitude correction ([correct_amp](parameter_list.md#correct_amp) = T)
The output is normalized such that the vertical component of the direct P arrival (or the horizontal component of the direct S arrival for S receiver functions) shows unit amplitude after the filtering.    

### Without amplitude correction ([correct_amp](parameter_list.md#correct_amp) = F)

Receiver functions are calculated without normalization. If this option is chosen, the amplitude varies in accordance with the filter parameter. See also [this page](https://github.com/akuhara/SEIS_FILO/wiki/Receiver-Function-Amplitude). 
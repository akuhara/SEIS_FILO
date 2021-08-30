# Receiver function forward computation


## How to run

`../../bin/recv_func_fwd recv_func_fwd.in`

## Inputs

### Main parameter file

* The filename must be specified as the first command line argument.
* The following can be used as a template.

```
#-----------------------------------------------------------------------
# Sample input to recv_func_fwd of SEIS_FILO package
#-----------------------------------------------------------------------

#--- Required parameters -----------------------------------------------

# Output filename
recv_func_out   = recv_func_fwd.out

# Input velocity model
vmod_in         = ../vmod/vmod.in

# Incident wave
rayp            = 0.05  # Ray parameter [s/km]
rf_phase        = P     # phase type P or S

# Low-pass filter
a_gauss         = 8.0   # Gaussian parameter
correct_amp     = T     # Amplitude correction (T or F)

# Deconvolution
deconv_flag     = T     # Deconvolved by vertical component (T) or not (F)
damp            = 0.00  # Water-level damping for deconvolution		 

# Timewindow
t_pre           = 3.0   # Time length before direct arrival [s]
n_smp           = 2048  # Number of elements in time series
delta           = 0.05  # Time interval [s]

#-----------------------------------------------------------------------

# Additive
noise_added = 0.0      # Standard deviation of noise (default=0.0)

# Random number seed (used for noise generation)
i_seed1     = 1001
i_seed2     = 22002
i_seed3     = 333003
i_seed4     = 4444004
```
### Velocity model file

* The filename must be specified in the parameter file.
* The following can be used as a template.
* For the ocean layer, set Vs to a negative value. This is allowed only for topmost layer.
* Thickness of the bottom layer is ignored.

```
#-----------------------------------------------------------------------
# Velocity model file:
# Sample input to disper_fwd & recv_func_fwd of SEIS_FILO package
#-----------------------------------------------------------------------
# Number of layers (including ocean & bottom layers)
         4
# Vp(km/s) Vs(km/s) Density(g/cm^3) Thickness(km)
       1.5     -3.0             1.0           2.0
       5.0      3.0             2.5           5.0
       7.0      4.0             3.0          10.0
       8.0      4.6             3.3           0.0
```


## Output 

### Receiver function waveform file
* SAC format file is used.

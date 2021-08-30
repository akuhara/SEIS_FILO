# Dispersion curve forward computation


## How to run

`../../bin/disper_fwd disper_fwd.in`

## Inputs

### Main parameter file
```
#-----------------------------------------------------------------------
# Sample input to disper_fwd of SEIS_FILO package
#-----------------------------------------------------------------------

#---- Required parameters ----------------------------------------------

# Output filename
disper_out     = disper_fwd.out

# Input velocity model
vmod_in        = ../vmod/vmod.in

# Output frequency or period range 
freq_or_period = freq    # Output observables (Freq or Period)
xmin           = 0.01    # Lower bound [Hz] / [s]
xmax           = 1.0     # Upper bound [Hz] / [s]
dx             = 0.05    # Interval    [Hz] / [s]

# Phase velocity range for root search
cmin           = 1.5     # Lower bound [km/s]
cmax           = 4.6     # Upper bound [km/s]
dc             = 0.02    # Interval    [km/s]

# Phase type
disper_phase   = R       # Currently only 'R' (Rayleigh wave) is supported
n_mode         = 0       # Mode number 
                         # (0: fundamental mode, 1: 1st higher mode, 
                         #  2: 2nd higher mode, ...)

#---- Optional parameters ----------------------------------------------

# Additive noise 
noise_added    = 0.0     # Standard deviation of noise (default=0.0)

# Random number seeds (used for noise generation)
i_seed1        = 1001
i_seed2        = 22002
i_seed3        = 333003
i_seed4        = 4444004


```

* [../vmod/vmod.in](https://github.com/akuhara/SEIS_FILO/blob/master/sample/vmod/vmod.in): input velocity model


## Output 
* disper.out: ASCII format file

```

test 

  test

```
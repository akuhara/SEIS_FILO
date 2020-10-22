# Velocity model file

This file describes the velocity model that is passed to the forward computation programs (i.e., `disper_fwd` and `recv_func_fwd`).

## Format
* Line 1: the number of layers (including seawater and half-space).
* Line 2 and after: Vp (km/s),  Vs (km/s), density (g/cm^3), thickness (km)

__NOTE__

* For the ocean layer, set Vs to be lower than 0 (this is allowed only for the topmost layer).

* The thickness of the bottom layer (half-space) is not used.

* Empty lines are ignored.

* Comment out by "#" works fine.

* Example is [here](https://github.com/akuhara/SEIS_FILO/blob/master/sample/vmod/vmod.in).


  
#!/bin/bash

model_file=../vmod/CPS_ocean.mod #CPS_input.mod
distance_file=distance_file

#sprep96 -M $model_file -d $distance_file -R -NMOD 3 # output: sdisp96.dat
sprep96 -M $model_file -DT 0.3 -NPTS 800 -R -NMOD 3 # output: sdisp96.dat
sdisp96                                             # output: sdisp96.ray
sregn96 -HR 2.0 -DER                                # output: sregn96.der
sdpder96 -R -TXT                                    # output: SRDER.TXT, SRDER.PLT (not used)

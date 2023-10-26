#!/bin/bash

model_file=../vmod/CPS_ocean.mod #CPS_input.mod
distance_file=distance_file

sprep96 -M $model_file -d $distance_file -R -NMOD 3
sdisp96
sregn96  -HR 2.0
sdpsrf96 -R -TXT
sdpegn96 -R -E -TXT

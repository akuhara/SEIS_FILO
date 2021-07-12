#!/bin/bash

model_file=../vmod/CPS_land.mod #CPS_input.mod
distance_file=distance_file

sprep96 -M $model_file -d $distance_file -R -NMOD 5
sdisp96
sregn96
sdpsrf96 -R -TXT
sdpegn96 -R -E -TXT


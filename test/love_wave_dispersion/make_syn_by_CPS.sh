#!/bin/bash

set -e

model_file=../vmod/CPS_land.mod #CPS_input.mod
distance_file=distance_file

sprep96 -M $model_file -d $distance_file -L -NMOD 3
sdisp96
slegn96
sdpsrf96 -L -TXT -YMIN 2.5 -YMAX 4.7
sdpegn96 -L -C -TXT

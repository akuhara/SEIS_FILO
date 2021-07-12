# Dispersion curve


Surface wave dispersion curves are calculated by the method of Saito and Kabasawa (1993), an improved version of the Thomson-Haskell matrix method. As is the case for other methods, it searches dispersion relation in the frequency (f)-phase velocity (c) domain by looking for zeros in a characteristic function. The SEIS_FILO programs involve several parameters that control the search process. 


## Root search process
The search begins with the evaluation of the characteristic function at (f, c) = ([xmin](parameter_list.md#xmin), [cmin](parameter_list.md#cmin)). This evaluation proceeds with an increment for c by [dc](parameter_list.md#dc) until we obtain [n_mode](parameter_list.md#n_mode)+1 zero crossings. Once the desired number of zero crossings occur, a binary search is performed for more tightly constraining the location of the zero-crossing, which is treated as phase velocity estimation at the corresponding frequency. A group velocity is then calculated by numerical differentiation centered on the estimated phase velocity. If c reaches [cmax](parameter_list.md#cmax) before finding the zero-crossing, the forward calculation terminates with an error message and will not proceed for higher frequencies. 

The procedure above are repeated over a given frequency range from 
[xmin](parameter_list.md#xmin) to [xmax](parameter_list.md#xmax) with an interval of [dx](parameter_list.md#dx); but the only difference is initial phase velocity. Instead of using [cmin](parameter_list.md#cmin) as the starting point, we can start from a point close to the prediction by a dispersion curve slope, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline{\frac{\partial&space;c}&space;{\partial&space;f}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline{\frac{\partial&space;c}&space;{\partial&space;f}}" title="dc/df" /></a>. In specific, the start value is given by <a href="https://www.codecogs.com/eqnedit.php?latex=\inline{c_{start}&space;=c_{previous}&space;&plus;a\frac{\partial&space;c}&space;{\partial&space;f}&space;df}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline{c_{start}&space;=c_{previous}&space;&plus;a\frac{\partial&space;c}&space;{\partial&space;f}&space;df}" title="\inline{c_{start} =c_{previous} +a\frac{\partial c} {\partial f} df}" /></a>, where a is an adjustable parameter which we fix at 3.5. Starting from this value, the zero-crossing first encountered is adopted as the target mode. 

## Full search mode
If one wish to obtain all dispersion curves within a certain region in the f-c domain given by ([xmin](parameter_list.md#xmin), [xmax](parameter_list.md#xmax)) × ([cmin](parameter_list.md#cmin), [cmax](parameter_list.md#cmax)), the full search mode may provide a solution. This mode is activated by setting [n_mode](parameter_list.md#n_mode)<0 and outputs characteristic functions at all grid points with intervals of [df](parameter_list.md#dc) and [dc](parameter_list.md#dc). 


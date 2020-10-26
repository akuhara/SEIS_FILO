# Parallel tempering

Parallel tempering is a technique offering more efficient MCMC sampling. Multiple MCMC samplings are performed in parallel, with each MCMC chain tempered by different temperatures. The temperatures control the performance of exploration in multidimensional space in the same manner as a popular simulated annealing method (_i.e._, The higher temperature is, the higher exploration ability). A vital element of this method is swapping temperatures between MCMC chains, which allows a non-tempered MCMC chain to make a long jump in multi-dimensional space. A review of the parallel tempering technique is given by [Sambridge (2014)](https://doi.org/10.1093/gji/ggt342) in the context of geophysics.

## Tuning parallel tempering

The performance of the parallel tempering depends on the number of MCMC chains and temperatures given to the chains. SEIS_FILO offers four adjustable parameters for this purpose: [n_proc](parameter_list.md#n_proc), [n_chain](parameter_list.md#n_chain), [n_cool](parameter_list.md#n_cool), and 
[temp_high](parameter_list.md#temp_high). 

!!! Note
    * The number of chains is set to n_proc × n_chains, where n_proc is the number of processes for parallel computing and n_chains is the number of MCMC chains per process.
    * Out of n_proc × n_chains chains, n_proc × n_cool chains are given unit temperature (i.e., non-tempered). 
    * Temperatures of the other chains are randomly distributed in logarithmic scale between 1 and temp_high, while the temperature of one chain is fixed at temp_high. 

 




# MCMC step

At each iteration, a trial model is proposed by perturbing the current model. One of the following 6 proposals is randomly selected. 

* Brith proposal: add a layer
* Death proposal: delete a layer
* Depth proposal: perturb a layer bottom depth
* Vs proposal: perturb Vs of a layer
* Vp proposal: perturb Vp of a layer
* Sigma proposal: perturb data noise standard deviation of a dataset

## Tuning parameters

The amount of perturbation can be adjusted via tuning parameters, with the exception of birth and death proposals. The higher the amount of these parameters, the more often a large perturbation amount is proposed. 

|Proposal type|Tuning parameter|
|:---|:---|
|Depth proposal |[dev_z](parameter_list.md#dev_z) |
|Vs prlposal|[dev_vs](parameter_list.md#dev_vs)|
|Vp prlposal|[dev_vp](parameter_list.md#dev_vp)| 
|Sigma prlposal|[dev_sig_rf](parameter_list.md#dev_sig_rf), [dev_sig_c](parameter_list.md#dev_sig_c), [dev_sig_u](parameter_list.md#dev_sig_u)|

## On birth proposal

In the birth proposal, the algorithm first selects the depth at which the new interface is created, and then assigns velocity to the layer. Note that both depth and velocity are randomly extracted from the corresponding prior probabilities.

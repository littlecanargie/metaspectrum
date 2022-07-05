# metaspectrum
This repository hosts R codes for meta-analysis spectrum methods under the structural equation modeling (SEM) ecosystem. SEM has been proven to be a useful tool for meta-analysis spectrum methods from its ability to specify complex, multivariate models with latent variables. Our code runs under the R environment and depends on [OpenMx](https://openmx.ssri.psu.edu/) to carryout the SEM specification and estimation procedures. 

## side_splitting.R
- Contains standalone function for side-splitting and symmetric side-splitting models for inconsistency in network meta-analysis.
- The architecture of the side-splitting model, which is a frequentist counterpart of the node-splitting model, can be found in [this 2010 SIM paper](https://doi.org/10.1002/sim.3767).
- The result of node(side)-splitting models may depend on an arbitrarily chosen reference treatment.  The symmetric side-splitting model was proposed to fix this asymmetry, which can be found in [this 2015 STATA journal paper](https://doi.org/10.1177/1536867X1501500403)

## evidence_splitting.R
- Contains standalone function for evidence-splitting model, an approach to inconsistency in network meta-analysis that takes evidence structure into account.
- The details of this model and why it provides more sensible results than previous models can be found in [this 2021 RSM paper](https://doi.org/10.1002/jrsm.1480).

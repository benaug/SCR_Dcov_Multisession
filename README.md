# SCR_Dcov_Multisession
Nimble MCMC samplers for Spatial Capture Recapture with density covariates and multiple closed sessions.

Uses Reversible Jump MCMC instead of data augmentation. Similar to other RJMCMC approaches except process model dimensions do not change,
only observation model. This allows us to let nimble assign samplers to them, tune them, etc. Pragmatic approach.

2 versions, M0 and Mb, each using data summed over occasions, so no occasion effects on detection possible without some modifications
to model file, data initalizer, and custom updates.
# SCR_Dcov_Multisession
Nimble MCMC samplers for Spatial Capture Recapture with density covariates and multiple closed sessions.

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

2 versions, M0 and Mb, each using data summed over occasions, so no occasion effects on detection possible without some modifications
to model file, data initalizer, and custom updates.
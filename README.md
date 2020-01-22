# Efficiency
Implementation of various improved efficiency statistics for skewed and periodic data

Python script implementing various efficiency estimators as defined in Lamontagne et al. (in Review), including: LBE, LBE', NSE, KGE, LNSE, PVSE. Also implemented are mixture estimators for LBE and LBE’: LBEm, LBE’m. R script implements LBEm and LBE'm

The Python script has been tested on large datasets.  The R script is still in development.

The Python script has functions MSE, NSE, LNSE, KGE_prime, and PVSE_prime that accept Nx1 vectors of simulations (S) and observations (O) and return estimators of the variance efficiencies.  The Python script also has LBE that accept Nx1 vectors S and O, and a logical variable indicating if the LBE estimator of E or E' should be returned.  Similarly, LBEm implements the mixture estimator of either E or E'.  LBEm also takes an Nx1 vector period that indicates the period (season, month, etc.) of either E or E'.  Note that Lamontagne et al. (2020) assumed monthly periods, but this is not a general rule for all hydrologic contexts.  The R script assumes monthly periods.

Repository is still under development and will be updated periodically.

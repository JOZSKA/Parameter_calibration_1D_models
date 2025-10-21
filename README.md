Code to calibrate model parameters using 1D EAT framework and eatpy library (Bruggeman et. al. 2024, https://gmd.copernicus.org/articles/17/5619/2024/). 
The code assumes that observation files have variables structured in gridded time x depth format, where depth index corresponding to the real depth in m.
The basic functionality is in calibrate_models.py and the different tasks can be run from run.py, run_seasons.py and run_obs_reductions.py.

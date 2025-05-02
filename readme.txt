The file sims_run_final.R will reproduce Table 2 in the paper where the user must change "setting = 1"  to "setting = 2" or "setting = 3" to get the specific setting results. Note that the code sets seeds and thus, should be fully reproducible. 

The main functions, all which are used in the simulation code are:
plugin_delta
plugin_delta_s
tmle_delta
tmle_delta_s
estimate_R

Several R packages are required; you must install these packages first before running the code. 

############################
********************************
************ dmd.py ************
********************************

#### 8 arguments in total are passed, beside the name of Python script ####

example commends: 
$python dmd.py 0.9 0.05 256 3 1 5 5 8
or
$python dmd.py 2.25 0.01 128 6 1 5 14 1

argument[0] -- Name of Python script
argument[1] -- ky
argument[2] -- delt                     # temporal resolution
argument[3] -- ntheta                   # spatial resolution
argument[4] -- nperiod                  # another spatial resolution parameter, but 3 is enough, 6/8/12 do not seem to be helpful
argument[5] -- DMD-order                # >1 means higher-order DMD
argument[6] -- # of SVD-truncation
argument[7] -- snapshot_start_ratio     # means the snapshot_start_point = len(t)//snapshot_start_ratio
argument[8] -- snapshot_end_ratio       # means the snapshot_end_point = len(t)*(100-snapshot_end_ratio)//100



***********************************
************ dmd.ipynb ************
***********************************

aky     = 0.9
delt    = 0.05   # temporal resolution
ntheta  = 256    # spatial resolution
nperiod = 3      # another sparial resolution parameter, but 3 is sufficient, 6/8/12 do not seem to be helpful

ns      = 1    # DMD-order
index   = 5    # number of SVD-truncation
snapshot_start_ratio = 5    # means the start_point = len(t)//snapshot_start_ratio
snapshot_end_ratio   = 8    # means the end_point = len(t)*(100-snapshot_end_ratio)//100



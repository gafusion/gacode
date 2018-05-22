#--------------------------------------------
# Environment variable setup for gacode
#--------------------------------------------
#!/bin/tcsh

setenv PATH $GACODE_ROOT/tgyro/bin:${PATH}
setenv PATH $GACODE_ROOT/gyro/bin:${PATH}
setenv PATH $GACODE_ROOT/cgyro/bin:${PATH}
setenv PATH $GACODE_ROOT/neo/bin:${PATH}
setenv PATH $GACODE_ROOT/vgen/bin:${PATH}
setenv PATH $GACODE_ROOT/tglf/bin:${PATH}
setenv PATH $GACODE_ROOT/glf23/bin:${PATH}
setenv PATH $GACODE_ROOT/le3/bin:${PATH}
setenv PATH $GACODE_ROOT/shared/bin:${PATH}
setenv PATH $GACODE_ROOT/profiles_gen/bin:${PATH}

if ( $?GACODE_ADD_ROOT ) then
   setenv PATH $GACODE_ADD_ROOT/freya/bin:${PATH}
   setenv PATH $GACODE_ADD_ROOT/prefreya/bin:${PATH}
endif

if ( $?PYTHONPATH ) then
   setenv PYTHONPATH $GACODE_ROOT/python:${PYTHONPATH}
else
   setenv PYTHONPATH $GACODE_ROOT/python
endif

if ( $?IDL_PATH ) then
   setenv IDL_PATH $GACODE_ROOT/gyro/vugyro:${IDL_PATH}
else
   setenv IDL_PATH $GACODE_ROOT/gyro/vugyro
endif

if ( $?NEURAL_ROOT ) then
    :
else if ( -d $GACODE_ROOT/../neural ) then
    setenv NEURAL_ROOT $GACODE_ROOT/../neural
else
    setenv NEURAL_ROOT
endif

if ( $?NEURAL_ROOT ) then
    setenv NN_LIB '-L$(NEURAL_ROOT) -I$(NEURAL_ROOT) -lbrainfusetf'
    setenv EPEDNN_MODEL $NEURAL_ROOT/eped1nn/models/EPED_mb_128_pow_norm_common_30x10.pb
    setenv TGLFNN_MODEL $NEURAL_ROOT/tglfnn/models/nn_SAT0_mb_1024_abs_reg_common_stair2x2x6.pb
    setenv NEONN_MODEL $NEURAL_ROOT/neonn/models/NEO_mb_64_common_30.pb
endif

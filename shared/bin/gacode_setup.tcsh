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

setenv EPEDNN_MODEL_DIR $GACODE_ROOT/../neural/eped1nn/models/EPED1_H_superH/
setenv TGLFNN_MODEL $GACODE_ROOT/../neural/tglfnn/models/nn_SAT0_mb_1024_abs_reg_common_stair2x2x6.pb

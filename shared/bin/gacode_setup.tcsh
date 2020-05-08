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
   setenv PYTHONPATH $GACODE_ROOT/f2py:$GACODE_ROOT/f2py/pygacode:${PYTHONPATH}
else
   setenv PYTHONPATH $GACODE_ROOT/f2py:$GACODE_ROOT/f2py/pygacode
endif

if ( $?IDL_PATH ) then
   setenv IDL_PATH $GACODE_ROOT/gyro/vugyro:${IDL_PATH}
else
   setenv IDL_PATH $GACODE_ROOT/gyro/vugyro
endif

setenv EPEDNN_MODEL_DIR $GACODE_ROOT/../neural/eped1nn/models/EPED1_H_superH/
setenv TGLFNN_MODEL_DIR $GACODE_ROOT/../neural/tglfnn/models/DIIID_ion_stiffness_60_rotation/

#!/task/imd/anaconda/bin/python
import os, sys
def make_clean():
    ret = os.system('cd $GACODE_ROOT; make clean')
    if ret:
        print('GACODE failed to make clean')
        sys.exit(125)
def make():
    ret = os.system('cd $GACODE_ROOT; make')
    if ret:
        print('GACODE failed to make')
        sys.exit(125)
class RegressionError():
    pass
def run_regressions():
    for code in ['gyro','neo','tglf','tgyro']:
        print('Starting regression testing of the %s code'%code)
        parallel = '-n 4'
        if code=='neo':
            parallel = ''
        logfile = '%s_regression.log'%(code)
        if not os.path.exists(logfile):
            ret = os.system('%s -r %s > %s'%(code,parallel,logfile))/256
            if ret:
                sys.exit(ret)
        num_tests = len(open(os.environ['GACODE_ROOT']+'/%s/tools/input/reg_list'%code,'r').read().strip().splitlines())
        if open(logfile,'r').read().count('PASS') < num_tests:
            print('Regression testing of the %s code failed'%code)
            raise RegressionError()
        print('Regression testing of the %s code PASSED %d tests'%(code,num_tests))
if '-clean' in sys.argv:
    make_clean()
if '-nomake' not in sys.argv:
    make()
try:
    run_regressions()
except RegressionError:
    print('You might try to `make clean` then run this script again')
    sys.exit(1)

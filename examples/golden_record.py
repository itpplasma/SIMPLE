
## Example of testing current (local) version of SIMPLE against
## main branch.

### ATTENTION: Some older variants of SIMPLE are not built via `make`,
###            which means the build step in build_and_run_legacy()
###            below needs to be adapted.

import os
from subprocess import run, Popen, PIPE, STDOUT

def get_test_data() -> str:
    cwd = os.getcwd()
    p = Popen([cwd+'/../test/test_data/get_test_data.sh'], cwd=cwd+"/../test/test_data" ,stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)
    return p.communicate(input='main')[0]
    
def build_and_run_simple():
    print(f"Build current version of SIMPLE")
    p = Popen("make", cwd="../")
    p.communicate()
    run(["cp","./simple.in","../build/simple.in"])
    run(["cp","../test/test_data/wout.nc","../build/wout.nc"])
    print(f"Run current version of SIMPLE")
    p = Popen("./simple.x", cwd="../build/")
    p.communicate()    
    run(["cp","../build/times_lost.dat","../test/tests/times_lost_new.dat"])

def build_and_run_legacy():
    print(f"Build legacy version of SIMPLE")
    cwd = os.getcwd()
    p = Popen("make", cwd=cwd+"/../test/test_data/simple_main/")
    p.communicate()
    run(["cp","./simple.in","../test/test_data/simple_main/build/simple.in"])
    run(["cp","../test/test_data/wout.nc","../test/test_data/simple_main/build/wout.nc"])
    print(f"Run legacy version of SIMPLE")
    p = Popen("./simple.x", cwd=cwd+"/../test/test_data/simple_main/build/")
    p.communicate()
    run(["cp","../test/test_data/simple_main/build/times_lost.dat","../test/tests/times_lost_old.dat"])

def cleanup():
    cwd = os.getcwd()
    os.remove(cwd+"/../test/tests/simple.in")
    os.remove(cwd+"/../test/tests/start.dat")
    os.remove(cwd+"/../test/tests/wout.nc")
    pass
    
if __name__ == "__main__":
    ## Checkout old version
    print(f"Getting test data: {get_test_data()}")

    ## Run current (local) version
    build_and_run_simple()

    ## Set up test data
    run(["ln","-sf","../test/test_data/wout.nc","../test/tests/wout.nc"])
    run(["ln","-sf","./simple.in","../test/tests/simple.in"])
    run(["ln","-sf","../build/start.dat","../test/tests/start.dat"])
    
    ## Run old version from git
    build_and_run_legacy()

    ## Test compares outputs
    cwd = os.getcwd()
    p = Popen(["python",cwd+"/../test/tests/test_against_legacy_behaviour.py"], cwd=cwd+"/../test/tests/")
    p.communicate()
    ## Clean artifacts
    cleanup()



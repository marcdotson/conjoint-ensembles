import numpy as np
from PYTHON import conjoint

def test_simulated_data_hbmnl(x):
    D = conjoint.get_data(x)
    result = conjoint.hbmnl(D)
    print(result)

def test_real_data_hbmnl(x):
    D = conjoint.get_data(x)
    result = conjoint.hbmnl(D, mu=(0,1), alpha=(0,10), lkj_param=5)
    print(result)

def test_simulated_data_ensemble(x):
    D = conjoint.get_data(x)
    result = conjoint.ensemble(D)
    print(result)

def test_real_data_ensemble(x):
    D = conjoint.get_data(x)
    result = conjoint.ensemble(D)
    print(result)

if __name__ == "__main__":
    # comment indicates pass

    ### PASS ###
    #test_simulated_data_hbmnl()
    #test_simulated_data_ensemble()
    #test_real_data_ensemble("04")

    ### FAIL ###
    test_real_data_hbmnl("01") # passes on 02 and 05


import numpy as np
from PYTHON import conjoint

def test_simulated_data_hbmnl():
    D = conjoint.get_data("06")
    result = conjoint.hbmnl(D)
    print(result)

def test_real_data_hbmnl():
    D = conjoint.get_data("03")
    result = conjoint.hbmnl(D)
    print(result)

if __name__ == "__main__":
    # comment indicates pass
    #test_simulated_data_hbmnl()
    test_real_data_hbmnl()

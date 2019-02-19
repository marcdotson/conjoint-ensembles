import conjoint

def test_hbmnl(x):
    D = conjoint.get_data(x)
    result = conjoint.hbmnl(D)
    print(result)
    
def test_ensemble(x):
    D = conjoint.get_data(x)
    result = conjoint.ensemble(D)
    print(result)

if __name__ == "__main__":

    x = "01"
    test_hbmnl(x)
    #test_ensemble(x)


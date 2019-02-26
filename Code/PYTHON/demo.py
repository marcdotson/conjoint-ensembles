import conjoint
import pandas as pd

dataset_list = [
        "01_PathologyNone/",
        "02_PathologyANA/",
        "03_PathologyScreening/",
        "04_PathologyMultiple/",
        "R1_PremiumChocolate/",
        "R2_GroundBeef/",
        "R3_ArtificialFlowers/",
        "R4_FloorCleaningServices/",
        "R5_InteriorPaint/"]

RESULTS = pd.DataFrame()
hbmnl_score = list()
ensemble_score = list()
hbmnl_time = list()
ensemble_time = list()
run_list = list()

for dataset in dataset_list:
    for run in range(10):
    
        path_to_data = "../../Data/{0}".format(dataset)
    
        if dataset in dataset_list[4:]:
            hbmnl_result,ensemble_result = conjoint.model_comparison(path_to_data, holdout=2, init_r=1)
    
        else:
            hbmnl_result,ensemble_result = conjoint.model_comparison(path_to_data, holdout=5)
    
        hbmnl_score.append(hbmnl_result['SCORE'])
        hbmnl_time.append(hbmnl_result['TIME'])
        ensemble_score.append(ensemble_result['SCORE'])
        ensemble_time.append(ensemble_result['TIME'])
        run_list.append(run)
    
RESULTS['dataset'] = dataset_list
RESULTS['hbmnl_score'] = hbmnl_score
RESULTS['ensemble_score'] = ensemble_score
RESULTS['hbmnl_time'] = hbmnl_time
RESULTS['ensemble_time'] = ensemble_time
RESULTS['run'] = run_list
RESULTS.to_csv("./results.csv")

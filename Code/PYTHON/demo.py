import conjoint
import pandas as pd

dataset = "03_PathologyScreening"
#dataset = "04_PathologyMultiple"
#dataset = "R1_PremiumChocolate"
#dataset = "R2_GroundBeef"
#dataset = "R3_ArtificialFlowers"
#dataset = "R4_FloorCleaningServices"
#dataset = "R5_InteriorPaint"

RESULTS = pd.DataFrame()
dataset_list = list()
hbmnl_score = list()
ensemble_score = list()
hbmnl_time = list()
ensemble_time = list()
run_list = list()

for run in range(10):
    print(dataset,run)

    path_to_data = "../../Data/{0}/".format(dataset)

    hbmnl_result,ensemble_result = conjoint.model_comparison(
                path_to_data,
                holdout=5, 
                iter=1000, 
                chains=4)
                #holdout=2, 
                #iter=500, 
                #chains=4,
                #control={'adapt_delta':.9, 'max_treedepth':3}, 
                #init_r=1)

    dataset_list.append(dataset)
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

RESULTS.to_csv("./{0}_results.csv".format(dataset))

print(RESULTS)

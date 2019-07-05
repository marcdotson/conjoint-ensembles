import utils
import conjoint



data_dict = utils.generate_factorial_design(10, 2, 2)
data_dict = utils.compute_beta_response(data_dict)
data_dict['Xtrain'] = data_dict['X'].copy()
data_dict['Ytrain'] = data_dict['Y'].copy()
data_dict['Xtest'] = data_dict['X'].copy()
data_dict['Ytest'] = data_dict['Y'].copy()

control_params = {'adapt_delta':.9, 'max_treedepth':10}
print(data_dict['X'].shape)

results, FIT = conjoint.hbmnl(
    data_dict,
    mu=[0,3.5],
    alpha=[0,1],
    return_fit=True,
    iter=1000,
    chains=4,
#    control=control_params,
    init_r=1)


print(results)
print(FIT)



### END ###

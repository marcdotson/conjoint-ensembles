from . import utils
import numpy as np
import pytest

def test_get_data():
    datasets = ['01_PathologyNone',
                '02_PathologyANA',
                '03_PathologyScreening',
                '04_PathologyMultiple']

    for dataset in datasets:
        # define filepath to data directory
        path_to_data = "../../Data/{0}/".format(dataset)
        # load correctly formatted data
        true_data = np.load(path_to_data + "data_dict.npz")
        # test against get_data function
        data_dict = utils.get_data(path_to_data)
        assert np.allclose(data_dict['X'], true_data['X'])
        assert np.allclose(data_dict['Y'], true_data['Y'])
        assert np.allclose(data_dict['Xtrain'], true_data['Xtrain'])
        assert np.allclose(data_dict['Ytrain'], true_data['Ytrain'])
        assert np.allclose(data_dict['Xtest'], true_data['Xtest'])
        assert np.allclose(data_dict['Ytest'], true_data['Ytest'])



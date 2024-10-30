import pyreadr
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.cross_decomposition import PLSRegression
from scipy.stats import pearsonr

# import data, access pandas df via key "None"
X_Govaere = pyreadr.read_r("datasets/X_Govaere.rds")[None]
X_Hoang = pyreadr.read_r("datasets/X_Hoang.rds")[None]
X_Pantano = pyreadr.read_r("datasets/X_Pantano.rds")[None]

# actual, clinical results
Y_Govaere = pyreadr.read_r("datasets/Y_Govaere.rds")[None]
Y_Hoang = pyreadr.read_r("datasets/Y_Hoang.rds")[None]
Y_Pantano = pyreadr.read_r("datasets/Y_Pantano.rds")[None]

# viper metric, predicted results
Y_viper_Govaere = pyreadr.read_r("datasets/viper_X_Govaere.rds")[None]
Y_viper_Hoang = pyreadr.read_r("datasets/viper_X_Hoang.rds")[None]
Y_viper_Pantano = pyreadr.read_r("datasets/viper_X_Pantano.rds")[None]

## tuning procedure
# training: Govaere, 10-fold cross validation to find the ideal (num of latent variables) hyperparameter

# StratifiedKFold will split datasets evenly amongst the classes
# shuffle=True will shuffle each class's samples before splitting into batches
# random_state=5 will guarantee reproducible output across multiple function calls

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=5)
num_lvs = range(2, 11) # test number of latent variables as hyperparameter
pearson_coeff = []
# need to split into Fibrosis and NAS scores

# X = X_Govaere
# Y = Y_Govaere
X = np.array([['Sample 1', 3, 50], ['Sample 2', 5, 60], ['Sample 3', 7, 70], ['Sample 4', 7, 2], ['Sample 5', 5, 40], ['Sample 6', 7, 50]])
Y = np.array([4, 5, 6, 4, 5, 6])
X = pd.DataFrame(X, columns=["Sample", "Gene1", "Gene2"])
Y = pd.DataFrame(Y, columns=["Sample", "Score"])

for latent_var in num_lvs:

    for i, (train_index, test_index) in enumerate(skf.split(X, Y)):
        print(f'Fold {i}:')
        print(f'Training index: {train_index}')
        print(f'Testing index: {test_index}')
        X_train = X.loc[train_index]
        X_test = X.loc[test_index]
        Y_train = Y.loc[train_index]
        Y_test = Y.loc[test_index]

        model = PLSRegression(n_components=latent_var, scale=False)
        corr, _ = pearsonr(X, Y)

# some pickle thing to download the model
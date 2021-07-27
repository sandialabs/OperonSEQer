# Pandas and numpy for data manipulation
import pandas as pd
import numpy as np
np.random.seed(42)
import pickle
import sys
import configargparse

 
# Matplotlib and seaborn for plotting
#import matplotlib.pyplot as plt
#%matplotlib inline
print('importing')
#import matplotlib
#matplotlib.rcParams['font.size'] = 16
#matplotlib.rcParams['figure.figsize'] = (9, 9)

print('starting scipy')
# Scipy helper functions
from scipy.stats import percentileofscore
from scipy import stats
import scipy.stats as st
print('done with scipy')
# Standard ML Models for comparison
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
import xgboost
from xgboost import XGBClassifier


# Splitting data into training/testing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import scale
print('still importing')

from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from matplotlib import pyplot
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.preprocessing import label_binarize
from functools import reduce 

# Metrics
from sklearn.metrics import mean_squared_error, mean_absolute_error, median_absolute_error

# Distributions
import scipy
import pandas as pd
import numpy as np

import pymc3 as pm

print('ok done importing')

if __name__ == '__main__':
    if len (sys.argv) != 9 :
        print("Usage: python voting.py -n operon_file -g gff_file -f prediction_files(6) -o1 operon_output -o2 no_operon_output\nNumber of arguements is " + str(len(sys.argv)))
        sys.exit (1)
    p = configargparse.ArgParser(description='given the six predictions from OperonSEQer, this script issues a voting tally')
    p.add('-n', required=True, help='operon file (full path)',dest='op')
    p.add('-g', required=True, help='gff with genes (full path)',dest='gff')
    p.add('-f', required=True, nargs='+', help='six input prediction files',dest='preds')
    p.add('-o', required=True, help='output name',dest='out')
    #print('we made it in')
    args=p.parse_args()


    df1=pd.read_csv(args.preds[0], delimiter = '\t')# Select only relevant variables
    df2=pd.read_csv(args.preds[1], delimiter = '\t')# Select only relevant variables
    df3=pd.read_csv(args.preds[2], delimiter = '\t')# Select only relevant variables
    df4=pd.read_csv(args.preds[3], delimiter = '\t')# Select only relevant variables
    df5=pd.read_csv(args.preds[4], delimiter = '\t')# Select only relevant variables
    df6=pd.read_csv(args.preds[5], delimiter = '\t')# Select only relevant variables


    df1b = df1[df1.duplicated(subset=['SysName1','SysName2'], keep=False)]
    df2b = df2[df2.duplicated(subset=['SysName1','SysName2'], keep=False)]
    df3b = df3[df3.duplicated(subset=['SysName1','SysName2'], keep=False)]
    df4b = df4[df4.duplicated(subset=['SysName1','SysName2'], keep=False)]
    df5b = df5[df5.duplicated(subset=['SysName1','SysName2'], keep=False)]
    df6b = df6[df6.duplicated(subset=['SysName1','SysName2'], keep=False)]

    ##multiple calls, only use full calls
    df1a = df1b.groupby(['SysName1', 'SysName2'])['XGBoost-BayesOpNoSc_preds'].mean().reset_index()
    df2a = df2b.groupby(['SysName1', 'SysName2'])['RF-BayesOp_preds'].mean().reset_index()
    df3a = df3b.groupby(['SysName1', 'SysName2'])['LogisticRegressionL2_preds'].mean().reset_index()
    df4a = df4b.groupby(['SysName1', 'SysName2'])['MLP-BayesOp_preds'].mean().reset_index()
    df5a = df5b.groupby(['SysName1', 'SysName2'])['GaussianNB-BayesOp_preds'].mean().reset_index()
    df6a = df6b.groupby(['SysName1', 'SysName2'])['SVM-rbf-BayesOp_preds'].mean().reset_index()


    #filter for values that agree (remove anything between 0 and 1)

    df1m = df1a[(df1a['XGBoost-BayesOpNoSc_preds'] == 1) | (df1a['XGBoost-BayesOpNoSc_preds'] == 0)]
    df2m = df2a[(df2a['RF-BayesOp_preds'] == 1) | (df2a['RF-BayesOp_preds'] == 0)]
    df3m = df3a[(df3a['LogisticRegressionL2_preds'] == 1) | (df3a['LogisticRegressionL2_preds'] == 0)]
    df4m = df4a[(df4a['MLP-BayesOp_preds'] == 1) | (df4a['MLP-BayesOp_preds'] == 0)]
    df5m = df5a[(df5a['GaussianNB-BayesOp_preds'] == 1) | (df5a['GaussianNB-BayesOp_preds'] == 0)]
    df6m = df6a[(df6a['SVM-rbf-BayesOp_preds'] == 1) | (df6a['SVM-rbf-BayesOp_preds'] == 0)]



    data_frames = [df1m, df2m, df3m, df4m, df5m, df6m]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['SysName1','SysName2'],
                                                how='outer'), data_frames)

    dfm=df_merged.drop_duplicates()
    dfm2=dfm.dropna()

    col_list=['XGBoost-BayesOpNoSc_preds',
           'RF-BayesOp_preds', 'LogisticRegressionL2_preds', 'MLP-BayesOp_preds',
           'GaussianNB-BayesOp_preds', 'SVM-rbf-BayesOp_preds']

    dfm2['sum']=dfm2[col_list].sum(axis=1)

    dfm2['one'] = np.where(dfm2['sum']>= 1, 1, 0)
    dfm2['two'] = np.where(dfm2['sum']>= 2, 1, 0)
    dfm2['three'] = np.where(dfm2['sum']>= 3, 1, 0)
    dfm2['four'] = np.where(dfm2['sum']>= 4, 1, 0)
    dfm2['five'] = np.where(dfm2['sum']>= 5, 1, 0)
    dfm2['six'] = np.where(dfm2['sum']>= 6, 1, 0)

    dfm2.to_csv('OperonSEQer_vote_result.csv', sep='\t', index=False)



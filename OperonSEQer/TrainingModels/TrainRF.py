# Pandas and numpy for data manipulation
import pandas as pd
import numpy as np
np.random.seed(42)
import pickle
print('START RF')
# Matplotlib and seaborn for plotting
#import matplotlib.pyplot as plt
#%matplotlib inline

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

#import matplotlib
#matplotlib.rcParams['font.size'] = 16
#matplotlib.rcParams['figure.figsize'] = (9, 9)

#import seaborn as sns

#from IPython.core.pylabtools import figsize

# Scipy helper functions
from scipy.stats import percentileofscore
from scipy import stats

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

# Splitting data into training/testing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import scale

from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score

# Metrics
from sklearn.metrics import mean_squared_error, mean_absolute_error, median_absolute_error

# Distributions
import scipy
import pandas as pd
import numpy as np

import pymc3 as pm
print('initial imports done')
##put together data frames
dfa2=pd.read_csv('BpseuK96243_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfb2=pd.read_csv('CdiffR20291_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfc2=pd.read_csv('Synec7002_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfd2=pd.read_csv('EcoliMG1655_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfe2=pd.read_csv('Selon7942_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dff2=pd.read_csv('Synec6803_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfg2=pd.read_csv('SaureUSA300_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables
dfh2=pd.read_csv('Bsubt3610_mergedInts_pred_st_20.txt', delimiter = '\t')# Select only relevant variables


dfa=dfa2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfb=dfb2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfc=dfc2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfd=dfd2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfe=dfe2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dff=dff2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfg=dfg2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
dfh=dfh2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]

#numeric_cols = [col for col in dfa if dfa[col].dtype.kind != 'O']
#numeric_cols

dfa[['KWp']]+=1e-300
dfa[['KWAIp']]+=1e-300
dfa[['KWBIp']]+=1e-300
dfa[['KWABp']]+=1e-300

dfb[['KWp']]+=1e-300
dfb[['KWAIp']]+=1e-300
dfb[['KWBIp']]+=1e-300
dfb[['KWABp']]+=1e-300

dfc[['KWp']]+=1e-300
dfc[['KWAIp']]+=1e-300
dfc[['KWBIp']]+=1e-300
dfc[['KWABp']]+=1e-300

dfd[['KWp']]+=1e-300
dfd[['KWAIp']]+=1e-300
dfd[['KWBIp']]+=1e-300
dfd[['KWABp']]+=1e-300

dfe[['KWp']]+=1e-300
dfe[['KWAIp']]+=1e-300
dfe[['KWBIp']]+=1e-300
dfe[['KWABp']]+=1e-300

dff[['KWp']]+=1e-300
dff[['KWAIp']]+=1e-300
dff[['KWBIp']]+=1e-300
dff[['KWABp']]+=1e-300

dfg[['KWp']]+=1e-300
dfg[['KWAIp']]+=1e-300
dfg[['KWBIp']]+=1e-300
dfg[['KWABp']]+=1e-300

dfh[['KWp']]+=1e-300
dfh[['KWAIp']]+=1e-300
dfh[['KWBIp']]+=1e-300
dfh[['KWABp']]+=1e-300


dfa[['log_KWp']] = np.log(dfa[['KWp']])
dfa[['log_KWAIp']] = np.log(dfa[['KWAIp']])
dfa[['log_KWBIp']] = np.log(dfa[['KWBIp']])
dfa[['log_KWABp']] = np.log(dfa[['KWABp']])

dfb[['log_KWp']] = np.log(dfb[['KWp']])
dfb[['log_KWAIp']] = np.log(dfb[['KWAIp']])
dfb[['log_KWBIp']] = np.log(dfb[['KWBIp']])
dfb[['log_KWABp']] = np.log(dfb[['KWABp']])

dfc[['log_KWp']] = np.log(dfc[['KWp']])
dfc[['log_KWAIp']] = np.log(dfc[['KWAIp']])
dfc[['log_KWBIp']] = np.log(dfc[['KWBIp']])
dfc[['log_KWABp']] = np.log(dfc[['KWABp']])

dfd[['log_KWp']] = np.log(dfd[['KWp']])
dfd[['log_KWAIp']] = np.log(dfd[['KWAIp']])
dfd[['log_KWBIp']] = np.log(dfd[['KWBIp']])
dfd[['log_KWABp']] = np.log(dfd[['KWABp']])

dfe[['log_KWp']] = np.log(dfe[['KWp']])
dfe[['log_KWAIp']] = np.log(dfe[['KWAIp']])
dfe[['log_KWBIp']] = np.log(dfe[['KWBIp']])
dfe[['log_KWABp']] = np.log(dfe[['KWABp']])

dff[['log_KWp']] = np.log(dff[['KWp']])
dff[['log_KWAIp']] = np.log(dff[['KWAIp']])
dff[['log_KWBIp']] = np.log(dff[['KWBIp']])
dff[['log_KWABp']] = np.log(dff[['KWABp']])

dfg[['log_KWp']] = np.log(dfg[['KWp']])
dfg[['log_KWAIp']] = np.log(dfg[['KWAIp']])
dfg[['log_KWBIp']] = np.log(dfg[['KWBIp']])
dfg[['log_KWABp']] = np.log(dfg[['KWABp']])

dfh[['log_KWp']] = np.log(dfh[['KWp']])
dfh[['log_KWAIp']] = np.log(dfh[['KWAIp']])
dfh[['log_KWBIp']] = np.log(dfh[['KWBIp']])
dfh[['log_KWABp']] = np.log(dfh[['KWABp']])

#dfftemp=dff.dropna()
#dff=dfftemp
#dff.loc['log_KWp'] = np.log(dff['KWp']+1e-1000)
#dff.loc['log_KWAIp'] = np.log(dff['KWAIp']+1e-1000)
#dff.loc['log_KWBIp'] = np.log(dff['KWBIp']+1e-1000)
#dff.loc['log_KWABp'] = np.log(dff['KWABp']+1e-1000)

dfa3=dfa[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfb3=dfb[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfc3=dfc[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfd3=dfd[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfe3=dfe[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dff3=dff[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfg3=dfg[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]
dfh3=dfh[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]

frames = [dfa3, dfb3, dfc3, dfd3, dfe3, dff3, dfg3, dfh3]
#frames = [dfa3, dfb3, dfe3]
df = pd.concat(frames)
# df = dfh3
###BELOW FOR ONE DATA FRAME
#df2=pd.read_csv('Allno7002_mergedLite_pred_st.txt', delimiter = '\t')# Select only relevant variables
#print(df2.columns)
#category_df = df2.select_dtypes('object')
#print(category_df)
#df=df2[['meanCov1','meanCov2','meanInt','Distance','Strand1','Strand2','bOp','Sep','MOGScore','GOScore','ExprSim','MC2_Rank','MC1_Rank','MI_Rank','Ratio_AtoInt_x_Rank','Ratio_BtoInt_x_Rank','Ratio_AtoB_x_Rank','proOpDBpred']]
#df=df2[['sdCov1','sdCov2','sdInt','meanCov1','meanCov2','meanInt','Distance','Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred','bOp']]
#df=df2[['sdCov1','sdCov2','sdInt','meanCov1','meanCov2','meanInt','Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred','bOp']]
#print(len(df.columns))
#print(df)
# Correlations of numericalq values
#df.corr()['pred'].sort_values()

df3=df.dropna(subset=['Length1'])
df=df3
print('files loaded')
##remove non-same strand
# Get names of indexes for which column Stock has value No
indexNames = df[ df['KWs'] == 0 ].index
# Delete these row indexes from dataFrame
df.drop(indexNames , inplace=True)
#df.drop(columns=['strandMatch'], inplace=True)
#df

##subset variables
def format_data(df):
    # Targets are final grade of student
    import sklearn
    labels = df['pred']
    df2=df.drop(columns=['pred'])
    #print(df)
    #print(df['SMRT'])
    #print(labels)
    #df.drop(columns=['sdCov1','sdCov2','sdInt','meanCov1','meanCov2','meanInt','Length1','Length2','LengthInt','SMRT','mergedPred','bOp'])
    
    #ADD#df2=pd.DataFrame(sklearn.preprocessing.scale(df[['Distance','KWStat','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp']]))
    #print(df2)  
    #df2.index = df.index
    #print(df.index.equals(df2.index) )  
    #df[['Distance','Ratio_AtoAend','Ratio_BtoBend','Length1','Length2']]=sklearn.preprocessing.scale(df[['Distance','Ratio_AtoAend','Ratio_BtoBend','Length1','Length2']])
    # Split into training/testing sets with 25% split
    #ADD#df2.columns=['Distance','KWStat','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp']
    #print(labels)
    #print(len(labels))
    #df2['SMRT']=labels
    #df2['strandMatch']=df['strandMatch']
    #print(df2['SMRT'])
    #print(len(df2['SMRT']))
    #print(df2)
    #print(labels)
    #ADD#df[['Distance','KWStat','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp']]=df2[['Distance','KWStat','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp']]
    X_train, X_test, y_train, y_test = train_test_split(df2, labels, 
                                                        test_size = 0.3,
                                                        random_state=42)
    
    return X_train, X_test, y_train, y_test

#print(df['SMRT'])
#single species
X_train, X_test, y_train, y_test = format_data(df)
#multi species
#X_train, y_train = format_data_2(df)
#X_test, y_test = format_data_2(testdf)

print(X_train.head())
#X_train, X_test, y_train, y_test = format_data(df)
##scale
scaler = sklearn.preprocessing.StandardScaler().fit(X_train)


X_train2=pd.DataFrame(scaler.transform(X_train))
X_train2.columns=X_train.columns
X_test2=pd.DataFrame(scaler.transform(X_test))
X_test2.columns=X_test.columns

#X_test
#
from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix

import pickle

import sklearn.utils.fixes
from numpy.ma import MaskedArray

sklearn.utils.fixes.MaskedArray = MaskedArray

import skopt
from bayes_opt import BayesianOptimization

# example of bayesian optimization (bayesian-optimization)
# https://medium.com/shortcutnlp/07-hyperparameter-tuning-a-final-method-to-improve-model-accuracy-b98ba860f2a6
from sklearn.model_selection import cross_validate

def estimator(n_estimators,max_depth, min_samples_split):
    # initialize model
    model = RandomForestClassifier(n_estimators=int(n_estimators), max_depth=int(max_depth), min_samples_split=int(min_samples_split),n_jobs=20)
    # set in cross-validation
    result = cross_validate(model, X_train, y_train, cv=10)
    # result is mean of test_score
    return np.mean(result['test_score'])

# def estimator(parameters):
#     print('estimator iteration start')
#     parameters=parameters[0]
#     n_estimators=int(parameters[0])
#     max_depth=int(parameters[1])
#     min_samples_split=int(parameters[2])
#     print(n_estimators)
#     print(max_depth)
#     print(min_samples_split)
#     model = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, n_jobs=10)
#     print('model set')
#     # set in cross-validation
#     #result = cross_val_score(model, X_train2, y_train, cv=10, n_jobs=-1)
#     result = cross_validate(model, X_train, y_train, cv=10)
#     # result is mean of test_score
#     print('estimator iteration end')
#     return np.mean(result['test_score'])
#     #return np.mean(result)

    
#define hyperparametes
#hparams = {"C": (0.1, 10), "gamma": (0.1, 10)}
#RF
hparams = {"n_estimators": (10, 1000),"max_depth": (1, 150),"min_samples_split": (2, 10)}
#hparams=[{'name':'n_estimators', 'type':'continuous', 'domain': (10, 1000)},
#    {'name': 'max_depth', 'type': 'continuous', 'domain': (1, 150)},
#    {'name': 'min_samples_split', 'type': 'continuous', 'domain': (2, 10)}]

#add optimizer
print('ready for optimization')
#from GPyOpt.methods import BayesianOptimization
# give model and hyperparameter to optmizer
svc_bayesopt = BayesianOptimization(estimator, hparams)
print('starting optimizer')

svc_bayesopt.maximize(init_points=10, n_iter=20, acq='ucb')
#batch_size=10
#num_cores=10
#svc_bayesopt = BayesianOptimization(f=estimator, domain=hparams, model_type='GP',acquisition_type='EI',batch_size=batch_size,num_cores=num_cores,exact_feval=True,evaluator_type = 'local_penalization')
#print('about to run optimizer')
#svc_bayesopt.run_optimization(max_iter=100)
#print('optimizer done')
#svc_bayesopt.plot_acquisition(filename='/home/rkrishn/CERES/RFacquisition')
#print("plot 1 done")
#svc_bayesopt.plot_convergence(filename='/home/rkrishn/CERES/RFlinConvergence')
#print("plot 2 done")
#svc_bayesopt.maximize(init_points=2, n_iter=5, acq='ucb')
#print(svc_bayesopt.max)
#y_bo = np.maximum.accumulate(-optimizer.Y).ravel()
#arg=np.argmax(svc_bayesopt.Y)

#finTarget=svc_bayesopt.max['target']
#finEst=int(svc_bayesopt.X[arg][0])
#finDep=int(svc_bayesopt.X[arg][1])
#finSplit=int(svc_bayesopt.X[arg][2])

finTarget=svc_bayesopt.max['target']
# finV=svc_bayesopt.X[arg][0]
#finGamma=svc_bayesopt.X[arg][1]
#finkernel=kernelDict[int(svc_bayesopt.X[arg][2])]
# finTarget=svc_bayesopt.max['target']
finEst=int(svc_bayesopt.max['params']['n_estimators'])
finDep=svc_bayesopt.max['params']['max_depth']
finSplit=int(svc_bayesopt.max['params']['min_samples_split'])


model=RandomForestClassifier(n_estimators=finEst, max_depth=finDep, min_samples_split=finSplit,n_jobs=10)
model.fit(X_train, y_train)
filename = 'RF.p'
pickle.dump(model, open(filename, 'wb'))
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
precision = precision_score(y_test, predictions)
recall = recall_score(y_test, predictions)
f1 = f1_score(y_test, predictions)
tn, fp, fn, tp = confusion_matrix(y_test, predictions).ravel()
specificity = tn / (tn+fp)
#tn / (tn + fp)
model_name = "RF"
model_list=["RF"]
results = pd.DataFrame(columns=['accuracy', 'precision','recall','specificity','f1'], index = model_list)
preds = {}
results.loc[model_name, :] = [accuracy, precision, recall, specificity, f1]
preds[model_name]=predictions

final=X_test2
final['RF_preds']=preds['RF']

results.to_csv('RF-result.csv', sep='\t', index=False)
final.to_csv('RF-predictions.csv', sep='\t', index=False)



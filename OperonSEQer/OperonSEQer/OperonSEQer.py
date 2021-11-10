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

import seaborn as sns

from IPython.core.pylabtools import figsize
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
from xgboost import XGBClassifier


# Metrics
from sklearn.metrics import mean_squared_error, mean_absolute_error, median_absolute_error

# Distributions
import scipy
import pandas as pd
import numpy as np

print('ok done importing')

def CI(a):
    '''
    confidence interval function
    '''
    return st.t.interval(0.95, len(a)-1, loc=np.mean(a), scale=st.sem(a))


if __name__ == '__main__':

	if len (sys.argv) < 3 :
		print("Usage: python OperonSEQer.py -i input_data_file [-p]\nNumber of arguements is " + str(len(sys.argv)))
		sys.exit (1)
	p = configargparse.ArgParser(description='run a model on data with or wihtout existing predictions ')
	p.add('-i', required=True, help='input file',dest='file')
	p.add('-p', required=False, help='predictions in file',dest='pred',action='store_true') 

	args=p.parse_args()


    dfa2=pd.read_csv(input, delimiter = '\t')# Select only relevant variables
    dfb=dfa2.dropna(subset=['Length1'])
    dfa2=dfb
    if pred==True:
        try:
            dfa=dfa2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']]
        except:
            print("Check the columns of the input file, make sure the following are in the file: ['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch','pred']")
            print("Use FormatData.py to get the correct file format")
    else:
        try:
            dfa=dfa2[['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch']]
        except:
            print("Check the columns of the input file, make sure the following are in the file: ['Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp','strandMatch']")
            print("Use FormatData.py to get the correct file format")
        
    #numeric_cols = [col for col in dfa if dfa[col].dtype.kind != 'O']
    #numeric_cols

    dfa[['KWp']]+=1e-300
    dfa[['KWAIp']]+=1e-300
    dfa[['KWBIp']]+=1e-300
    dfa[['KWABp']]+=1e-300


    dfa[['log_KWp']] = np.log(dfa[['KWp']])
    dfa[['log_KWAIp']] = np.log(dfa[['KWAIp']])
    dfa[['log_KWBIp']] = np.log(dfa[['KWBIp']])
    dfa[['log_KWABp']] = np.log(dfa[['KWABp']])

    dfa3=dfa[['Length1','Length2','LengthInt','KWs','log_KWp','KWAIs','log_KWAIp','KWBIs','log_KWBIp','KWABs','log_KWABp','strandMatch','pred']]

    df=dfa3
    ##remove non-same strand
    # Get names of indexes for which column Stock has value No
    indexNames = df[ df['KWs'] == 0 ].index
    # Delete these row indexes from dataFrame
    df.drop(indexNames , inplace=True)
    dfa2.drop(indexNames , inplace=True)
    #df.drop(columns=['strandMatch'], inplace=True)
    #df

    ##subset variables
    def format_data(df):
        # Targets are final grade of student
        import sklearn
        labels = df['pred']
        df2=df.drop(columns=['pred'])
        return df2, labels

    #single species
    if pred==True:
        print('Predictions provided')
        X_test, y_test = format_data(df)
        #print(y_test)
    else:   
        X_test=df.copy()
    #multi species
    #X_train, y_train = format_data_2(df)
    #X_test, y_test = format_data_2(testdf)

    #X_train, X_test, y_train, y_test = format_data(df)
    ##scale
    X_test2=X_test.copy()
    X_test2.columns=X_test.columns
    scaler = sklearn.preprocessing.MinMaxScaler(feature_range = (0, 1)).fit(X_test)
    X_test3=pd.DataFrame(scaler.transform(X_test))
    X_test3.columns=X_test.columns

    # example of bayesian optimization (bayesian-optimization)
    # https://medium.com/shortcutnlp/07-hyperparameter-tuning-a-final-method-to-improve-model-accuracy-b98ba860f2a6
    from sklearn.model_selection import cross_validate

    models=['OperonSEQer/XGB.p','OperonSEQer/MLP.p','OperonSEQer/SVM.p','OperonSEQer/RF.p','OperonSEQer/LR.p','OperonSEQer/GNB.p']
    #models=['OperonSEQer/MLP.p']
    #model2=pickle.load( open( args.model, "rb" ) )

    def load_pkl(fname):
        with open(fname, 'rb') as f:
            obj = pickle.load(f)
        return obj
        
    for model in models:
        model2=load_pkl(model)
        model_name = str(model).split('.p')[0].split('/')[1]
        if pred==True:
            print('Running '+str(model_name))
            print('testing predicitons by boostrap')
            # configure bootstrap
            n_iterations = 100
            n_size = int(len(y_test) * 0.10)
            stats=[]
    # 	    print('this is X_test2')
    # 	    print(X_test2)
            if (model_name=='XGB') or (model_name=='RF'):
                tosample=X_test2.copy()
            else:
                tosample=X_test3.copy()
            #print('hi')
            #print(len(tosample['Length1']))
            #print(tosample)
    # 	    print('this is y_test')
    # 	    print(len(y_test))
    # 	    print(type(y_test))
            tosample['pred']=y_test.tolist()
            #print('here is tosample')
            #print(tosample)
            for i in range(n_iterations):
                # prepare train and test sets
                #print('hello before')
                data = sklearn.utils.resample(tosample, n_samples=n_size)
                #print('here is data')
                #print(data)
                labels=data['pred']
                testset=data.drop(columns=['pred'])
                #print('here is test')
                #print(testset)
                # evaluate model
                #print(model2)
                predictions = model2.predict(testset)
                accuracy = accuracy_score(labels, predictions)
                precision = precision_score(labels, predictions)
                recall = recall_score(labels, predictions)
                f1 = f1_score(labels, predictions)
                tn, fp, fn, tp = confusion_matrix(labels, predictions).ravel()
                specificity = tn / (tn+fp)
                stats.append([accuracy, precision, recall, specificity, f1])
            statsDF=pd.DataFrame(stats, columns=['accuracy', 'precision', 'recall','specificity','f1'])
            cols=statsDF.columns
            #print(model_name)
            #print('statsDF')
            #print(statsDF)
            means=statsDF.mean(axis=0)
            CIs=statsDF[cols].apply(CI)
            if (model_name=='XGB') or (model_name == 'RF'):
                predictions = model2.predict(X_test2)
                lr_probs = model2.predict_proba(X_test2)
            else:
                predictions = model2.predict(X_test3)
                lr_probs = model2.predict_proba(X_test3) 
            # keep probabilities for the positive outcome only
            lr_probs = lr_probs[:, 1]
            ns_probs = [0 for _ in range(len(y_test))]
            # calculate scores
            ns_auc = roc_auc_score(y_test, ns_probs)
            lr_auc = roc_auc_score(y_test, lr_probs)
            # summarize scores
            #print('No Skill: ROC AUC=%.3f' % (ns_auc))
            #print('Model: ROC AUC=%.3f' % (lr_auc))
            # calculate roc curves
            ns_fpr, ns_tpr, _ = roc_curve(y_test, ns_probs)
            lr_fpr, lr_tpr, _ = roc_curve(y_test, lr_probs)
            # plot the roc curve for the model
            pyplot.clf()
            fig = pyplot.figure(figsize=(5,5))
            ax = fig.add_subplot()
            ax.plot(ns_fpr, ns_tpr, linestyle='--', label='No Skill')
            ax.plot(lr_fpr, lr_tpr, marker='.', label='Model')
            ax.title.set_text(str(model_name)+' ROC curve')
            # axis labels
            ax.set(xlabel='False Positive Rate', ylabel='True Positive Rate')
            # show the legend
            ax.legend(bbox_to_anchor=(1, 1), loc="lower right")
            ax.text(0.55, 0.5,'No Skill: ROC AUC=%.3f' % (ns_auc))
            ax.text(0.55, 0.4,'Model: ROC AUC=%.3f' % (lr_auc))
            # show the plot
            if path!=None:
                fig.savefig(str(path)+str(model_name)+' ROC curve.png')
            else:
                fig.savefig(str(model_name)+' ROC curve.png')
            # Compute ROC curve and ROC area for each class
        else:
            print('Running '+str(model_name)+' - no predictions provided')
            if (model_name=='XGB') or (model_name == 'RF'):
                predictions = model2.predict(X_test2)
            else:
                predictions = model2.predict(X_test3)
        #print(predictions)
        #tn / (tn + fp)
        #model_name = str(args.model).split('.p')[0].split('/')[-1]
        model_list=[model_name]
        #print(model_list)
        if pred==True:
            results = pd.DataFrame(columns=['accuracy', 'precision','recall','specificity','f1'], index = model_list)
            results.loc['means'] = means.tolist()
            results.loc['lower CI'] = [a_tuple[0] for a_tuple in CIs]
            results.loc['upper CI'] = [a_tuple[1] for a_tuple in CIs]
            results.index.name = ''
            results.dropna(inplace=True)
            results.reset_index(inplace=True)
        preds = {}
        preds[model_name]=predictions
        final=pd.DataFrame(columns=['A'])
        final['SysName1']=dfa2['SysName1']
        final['SysName2']=dfa2['SysName2']
        #print(len(final['SysName1']))
        #print(len(preds[model_name]))
        final[str(model_name)+'_preds']=preds[model_name]
        final.drop(columns=['A'], inplace=True)
        if path!=None:
            nameRes=str(path)+str(model_name)+'_result.csv'
            namePred=str(path)+str(model_name)+'_predictions.csv'
        else:
            nameRes=str(model_name)+'_result.csv'
            namePred=str(model_name)+'_predictions.csv'            
        if pred==True:
            results.to_csv(nameRes, sep='\t', index=False)
        final.to_csv(namePred, sep='\t', index=False)
    print('done with OperonSEQer')

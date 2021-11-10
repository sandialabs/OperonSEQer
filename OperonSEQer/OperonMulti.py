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
    if len (sys.argv) < 10 :
        print("Usage: python OperonMulti.py -f prediction_file -o operon_output -t threshold -g gff_file -d gff_delimiter\nNumber of arguements is " + str(len(sys.argv)))
        sys.exit (1)
    p = configargparse.ArgParser(description='given the six predictions from OperonSEQer and a threshold, this script strings together multigene operons')
    p.add('-f', required=True, help='six input predictions',dest='Opreds')
    p.add('-o', required=True, help='output name',dest='out')
    p.add('-t', required=False, help='threshold for call', dest='thr', type=int, choices=range(1,7))
    p.add('-g', requied=False, help='gff file for stringing together operons (no header)', dest='gff')
    p.add('-d', requied=False, help='gff file delimiter (default as tab) - tab, space or comma', dest='deli')
    #print('we made it in')
    ##this is the format of the gff file:
    #0     Chromosome  ena  gene      190      255  .  +  .  b0001
    #1     Chromosome  ena  gene      337     2799  .  +  .  b0002

    args=p.parse_args()
    
    preds=pd.read_csv(args.Opreds, sep='\t')
    
    if args.thr:
        if args.gff:
            print('creating thresholded operon file')
        else:
            print('ERROR - if you want to string together operons, you need to provide a gff file')
            sys.exit(1)
        threshdict={1:'one',2:'two',3:'three',4:'four',5:'five',6:'six'}
        thresh=threshdict[args.thr]
        operons = preds[["SysName1", "SysName2", thresh]]
        operons.columns=['Gene1','Gene2','pred']
        if args.deli:
            if args.deli=='tab':
                gff=pd.read_csv('/home/rkrishn/projects/CERES/Raga/Genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37_lite.txt', sep='\t', header=None)
            elif args.deli=='comma':
                gff=pd.read_csv('/home/rkrishn/projects/CERES/Raga/Genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37_lite.txt', sep=',', header=None)
            elif args.deli=='space':
                gff=pd.read_csv('/home/rkrishn/projects/CERES/Raga/Genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37_lite.txt', sep=' ', header=None)
        else:
            gff=pd.read_csv('/home/rkrishn/projects/CERES/Raga/Genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37_lite.txt', sep='\t', header=None)
        genelist=list(gff.iloc[:, 8])
        multioperons={}
        n=0
        a=0
        print(len(genelist))
        while a < len(genelist)-1:
            inoperon=1
            templist=[genelist[a]]
            while inoperon==1:
                subset=operons[operons['Gene1'].str.match(genelist[a])]
                if not subset.empty:
                    if len(subset['pred'])!=1:
                        print('oops two')
                    if subset.iloc[0]['pred']==0:
                        inoperon=0
                        a+=1
                        continue
                    elif subset.iloc[0]['pred']!=0:
                        inoperon=1
                        if subset.iloc[0]['Gene2']==genelist[a+1]:
                            templist.append(genelist[a+1])
                            a+=1
                else:
                    a+=1
                    inoperon=0
            if len(templist)>1:
                multioperons[n]=templist
            n+=1

        df = pd.DataFrame.from_dict(multioperons, orient='index') #Convert dict to df
        df.to_csv(str(args.out) + '_OperonSEQer_voteThreshold_operonList.csv',header=False, index=False) #Convert df to csv

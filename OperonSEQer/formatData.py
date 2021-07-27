import pandas as pd
import numpy as np
np.random.seed(42)

# Distributions
import scipy
import pandas as pd
import numpy as np

import pymc3 as pm
import pickle
import csv

    
##format cross-species testing
def format_data_2(df):
    label = df['pred']
    if 'KWp' in df:
        df.drop(columns=['KWp','KWAIp','KWBIp','KWABp'])
    return df, label


if __name__ == '__main__':
    if len (sys.argv) < 3 :
        print("Usage: python FormatData.py -i input file  -p [optional argument to save as pickle]\nNumber of arguements is " + str(len(sys.argv)))
        sys.exit (1)
    p = configargparse.ArgParser(description='this script takes in a coverage file and outputs a table of pertinent features for operon prediction')
    p.add('-i', required=True, help='input data file',dest='input')
    p.add('-p', required=False, help='output file name',dest='pickle')
    #print('we made it in')
    args=p.parse_args()
    
    df=pd.read_csv(input, delimiter = '\t')
    
    df['KWp']+=1e-300
    df['KWAIp']+=1e-300
    df['KWBIp']+=1e-300
    df['KWABp']+=1e-300
    
    df['log_KWp'] = np.log(df['KWp'])
    df['log_KWAIp'] = np.log(df['KWAIp'])
    df['log_KWBIp'] = np.log(df['KWBIp'])
    df['log_KWABp'] = np.log(df['KWABp'])
    
    df3=df.dropna(subset=['Length1'])
    df=df3
    indexNames = df[ df['KWs'] == 0 ].index
    # Delete these row indexes from dataFrame
    df.drop(indexNames , inplace=True)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    if pred==True:
        df['pred'] = df['pred'].fillna(0)
    if path!=None:
        filename=str(path)+'Data_Formatted.p'
        outfile=str(path)+'Data_Formatted.txt'
    else:
        filename='Data_Formatted.p'
        outfile='Data_Formatted.txt'
    if pickleit==True:
        pickle.dump(df, open(filename, "wb"))
    df.to_csv(outfile, sep='\t')

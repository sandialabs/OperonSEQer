import csv
import sys
import configargparse
import pandas as pd


if __name__ == '__main__':

	if len (sys.argv) < 5 :
		print("Usage: python mergepred2021.py -f file\nNumber of arguements is " + str(len(sys.argv)))
		sys.exit (1)

	p = configargparse.ArgParser(description='this script merges predicitons with data - make sure format is Gene1 Gene2 Prediction (True or False)')
	p.add('-f', required=True, help='input file with full path',dest='file')
	p.add('-p', required=True, help='predictions in the format Gene1\tGene1\tPred with headers',dest='pred')
	args=p.parse_args()
	print('we have started')
	
    with open(input,"r") as f:
        reader = csv.reader(f, delimiter="\t")
        header1=next(reader)
        data = list(reader)
    
    with open(pred,"r") as g:
        reader2 = csv.reader(g, delimiter="\t")
        header2=next(reader2)
        data2 = list(reader2)
            
    df = pd.DataFrame(data, columns = ['SysName1','SysName2','meanCov1','sdCov1','meanCov2','sdCov2','meanInt','sdInt','Distance','Strand1','Strand2','Length1','Length2','LengthInt','KWs','KWp','KWAIs','KWAIp','KWBIs','KWBIp','KWABs','KWABp'])
    df.loc[(df['Strand1'] == df['Strand2']), 'strandMatch'] = 1
    df['strandMatch'] = df['strandMatch'].fillna(0)
    df2 = pd.DataFrame(data2, columns = ['SysName1','SysName2','pred'])
    
    merged_inner = pd.merge(df, df2, how='left', left_on=['SysName1','SysName2'], right_on=['SysName1','SysName2'])
    if path!=None:
        outfile=str(path)+'Data_predsMerged.txt'
        merged_inner.to_csv(outfile, sep='\t')
    else:
        merged_inner.to_csv('Data_predsMerged.txt', sep='\t')

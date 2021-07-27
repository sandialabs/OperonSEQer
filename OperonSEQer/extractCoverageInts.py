import sys
import configargparse
import csv
import random
from statistics import mean
import numpy as np
from scipy import stats


if __name__ == '__main__':
    if len (sys.argv) != 9 :
        print("Usage: python extractCoverage.py -c coverage_file -g gff_file -o outfile -p outpath\nNumber of arguements is " + str(len(sys.argv)))
        sys.exit (1)
    p = configargparse.ArgParser(description='this script takes in a coverage file and outputs a table of pertinent features for operon prediction')
    p.add('-c', required=True, help='coverage file (full path)',dest='cov')
    p.add('-g', required=True, help='gff with genes (full path)',dest='gff')
    p.add('-p', required=True, help='output path',dest='path')
    p.add('-o', required=True, help='output',dest='out')
    #print('we made it in')
    args=p.parse_args()

    with open(args.cov,"r") as f:
        reader = csv.reader(f, delimiter="\t")
        #next(reader)
        cov = list(reader)
        
    #to ask whether genes are in an operon, look at intergenic region versus. If it is less than 200 bp, check whether intergenic expression is within 2 fold of gene on either side
    #load genes

    with open(args.gff,"r") as f:
        reader = csv.reader(f, delimiter=" ")
        #next(reader)
        rna = list(reader)

    covD={}
    for a in cov:
        if a[0] in covD:
            covD[a[0]][int(a[1])]=float(a[2])
        else:
            covD[a[0]]={}
            covD[a[0]][int(a[1])]=float(a[2])
	
    data=[]
    length=[]
    ##values for gene x, gene y, mean_cov_genex, sd_cov_genex, mean_cov_geney, sd_cov_geney, mean intergene x-y, sd intergene x-y,distance, strand x, strand y
    for a in range(0,len(rna)-1):
        chr_a=rna[a][0]
        chr_b=rna[a+1][0]
        start_a=int(rna[a][3])
        stop_a=int(rna[a][4])
        start_i=int(rna[a][4])+1
        stop_i=int(rna[a+1][3])-1
        dist_genes=stop_i - start_i
        start_b=int(rna[a+1][3])
        stop_b=int(rna[a+1][4])
        strand_a=rna[a][6]
        strand_b=rna[a+1][6]
        covA=[]
        covB=[]
        #print(start_i)
        #print(stop_i)
        if chr_a != chr_b:
            continue
        else:
            x=start_a
            while x < stop_a:
                try:
                    covA.append(covD[chr_a][x])
                except:
                    print('aloha A')
                    #print(args.cov)
                    #print(args.gff)
                    print('could not do it ' + str(x)+ ' ' + str(chr_a))
                x+=1
            z=start_b
            while z < stop_b:
                try:
                    #print(covD[chr_b][z])
                    covB.append(covD[chr_b][z])
                except:
                    print('aloha B')
                    #print(covD.keys())
                    #print(args.cov)
                    print(start_b)
                    print('could not do it ' +str(stop_b) +' is end of gene and '+ str(z)+ ' is location ' + str(chr_b))
                z+=1
            ##get mean of only ends of A and B
            meanA=mean(covA[-100:])
            sdA=np.std(covA[-100:])
            meanB=mean(covB[:100])
            sdB=np.std(covB[:100])       
            if int(start_i) >= int(stop_i):
                covMid=[]
                #print('testing')
                #print(start_i)
                #print(stop_i)
                #print(dist_genes)
                #print('testing_done')
                newMid=int((int(start_b)+int(stop_a))/2)
                newStart=newMid-25
                newStop=newMid+25
                z = newStart
                while z < newStop:
                    covMid.append(covD[chr_a][z])
                    z+=1
                meanMid=mean(covMid)
                sdMid=np.std(covMid) 
                #print(covMid)
                #print(covA)
                #print(covB)
            else:
                covMid=[]
                y=start_i
                while y < stop_i:
                    covMid.append(covD[chr_a][y])
                    y+=1
                meanMid=mean(covMid)
                sdMid=np.std(covMid) 
                if len(covMid)>50:
                    mid=int(len(covMid)/2)
                    startMid=mid-25
                    stopMid=mid+25
                    meanMid=mean(covMid[startMid:stopMid])
                    sdMid=np.std(covMid[startMid:stopMid])
                    meanAend=mean(covMid[0:25])
                    sdAend=np.std(covMid[0:25])
                    meanBend=mean(covMid[-25:])
                    sdBend=np.std(covMid[-25:])
            #Kruskal-wallis test
            try:
                test=stats.kruskal(covA, covB, covMid)
            except ValueError:
                #print('error in mutli test')
                #print(a)
                test='not valid'
            try:
                testAI=stats.kruskal(covA, covMid)
            except ValueError:
                #print('error in AI test')
                #print(a)
                #print(covA)
                #print(covMid)
                testAI='not valid'
            try:
                testBI=stats.kruskal(covB, covMid)
            except ValueError:
                #print('error in BI test')
                #print(a)
                testBI='not valid'
            try:
                testAB=stats.kruskal(covB, covA)
            except ValueError:
                #print('error in AB test')
                #print(a)
                testAB='not valid'
            #print('here is ab')
            #print(testAB)
            if test != 'not valid' and testAI != 'not valid' and testBI != 'not valid' and testAB != 'not valid':
            	data.append([rna[a][8],rna[a+1][8],meanA,sdA,meanB,sdB,meanMid,sdMid,dist_genes,strand_a,strand_b,len(covA),len(covB),len(covMid),test[0],test[1],testAI[0],testAI[1],testBI[0],testBI[1],testAB[0],testAB[1]])
            #length.append([len(covA),len(covB),len(covMid)])
    outfile=str(args.path)+str(args.out)
    #print(outfile)
    header=['']
    with open(outfile, 'w') as writeFile:
        writer = csv.writer(writeFile,delimiter='\t')
        writer.writerows(header)
        writer.writerows(data)
    
    #outfile2=str(args.path)+str(args.out)+'_2'
    #with open(outfile2, 'w') as writeFile2:
    #    writer=csv.writer(writeFile2,delimiter='\t')
    #    writer.writerows(length)

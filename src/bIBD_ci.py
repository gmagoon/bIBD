import sys
import os.path
import math
import bIBD_ci
from collections import OrderedDict

if __name__ == '__main__':
    llrinpath = sys.argv[1]
    ibdinpath = sys.argv[2]
    cl = float(sys.argv[3]) #approximate confidence level, e.g. 0.95
    ibdoutpath = sys.argv[4]

    cithresh=-math.log(1.0-cl) #the "confidence interval" considered here is not statistically robust, but it should be a useful proxy for one; 0.95 cl will give approximately 3.0 here
    print 'cithresh = ' + str(cithresh)

    #read in segment info
    ibdData=[]
    ibdin=open(ibdinpath,'r')
    for line in ibdin:
        split=line.strip().split('\t')
        posL = int(split[0])
        posR = int(split[1])
        ibdData.append((posL,posR))
    ibdin.close()

    #we could either read the llr data into memory once, or read from disk twice (first pass to identify min/max cllr associated with boundaries, second pass to identify confidence intervals)
    #here we will read into memory for simplicity (higher memory use but should be less disk/CPU intensive)
    cllrData=[]
    cllrDataIndex={} #stores the corresponding array index for a position key
    cllr=0.0
    i=0
    cllrData.append((0,cllr))
    cllrDataIndex[0]=i
    llrin = open(llrinpath,'r')
    for line in llrin:
        split=line.strip().split('\t')
        pos = int(split[2])
        llr = float(split[9])
        if not math.isnan(llr): #skip nan sites
            i=i+1
            cllr=cllr+llr
            cllrData.append((pos,cllr))
            cllrDataIndex[pos]=i
    llrin.close()

    ibdout = open(ibdoutpath,'w')
    for ibd in ibdData:
        posL=ibd[0]
        posR=ibd[1]
        leftcllr=cllrData[cllrDataIndex[posL]-1] #this is a tuple
        rightcllr=cllrData[cllrDataIndex[posR]]
        #left CI of left boundary
        i=cllrDataIndex[posL]-1 #could be reused from above
        i=i-1
        a = 'chrStart' #default value if below loop fails fo find anything
        while i >= 0:
            if cllrData[i][1] - leftcllr[1] > cithresh:
                a = str(cllrData[i][0])
                break
            i=i-1
        #right CI of left boundary
        i=cllrDataIndex[posL]-1 #could be reused from above #this is the index of the minimum
        i=i+1
        b = 'chrEnd' #default value if below loop fails fo find anything #this shouldn't happen
        while i < len(cllrData):
            if cllrData[i][1] - leftcllr[1] > cithresh:
                b = str(cllrData[i][0])
                break
            i=i+1
        assert not b=='chrEnd'
        #left CI of right boundary
        i=cllrDataIndex[posR] #could be reused from above
        i=i-1
        c = 'chrStart' #default value if below loop fails fo find anything #this shouldn't happen
        while i >= 0:
            if rightcllr[1] - cllrData[i][1] > cithresh:
                c = str(cllrData[i][0])
                break
            i=i-1
        assert not c=='chrStart'
        #right CI of left boundary
        i=cllrDataIndex[posR] #could be reused from above #this is the index of the maximum
        i=i+1
        d = 'chrEnd' #default value if below loop fails fo find anything
        while i < len(cllrData):
            if rightcllr[1] - cllrData[i][1] > cithresh:
                d = str(cllrData[i][0])
                break
            i=i+1
        #determine max drawdown (here calculated as a negative or zero value), a measure of segment reliability; significantly negative values could indicate false positive, de novo mutation, genotyping error, structural variation, or undetected break
        maxcllr=cllrData[cllrDataIndex[posL]][1]
        maxddn=0.0
        maxddnstartpos=0
        potentialmaxddnstartpos=0
        maxddnendpos=0
        storePotentialStartPos=True #store potential start position as the position after a new maximum
        for i in range(cllrDataIndex[posL]+1,cllrDataIndex[posR]): #note we skip (cllrDataIndex[posL]-1) and cllrDataIndex[posR] as these are known min and max, respectively; also (cllrDataIndex[posL]) is known to be higher so this is used for initialization above and is also skipped
            if storePotentialStartPos:
                potentialmaxddnstartpos=cllrData[i][0]
                storePotentialStartPos=False
            cllr=cllrData[i][1]
            if cllr > maxcllr:
                maxcllr=cllr
                storePotentialStartPos=True
            ddn=cllr-maxcllr
            if ddn < maxddn:
                maxddn=ddn
                maxddnstartpos=potentialmaxddnstartpos
                maxddnendpos=cllrData[i][0]
        ibdout.write('\t'.join([str(posL),str(posR),str(rightcllr[1]-leftcllr[1]),str(cllrDataIndex[posR]-(cllrDataIndex[posL]-1)),str(maxddnstartpos),str(maxddnendpos),str(maxddn),a,b,c,d])+'\n')
    ibdout.close()

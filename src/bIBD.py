import sys
import os.path
import math
import bIBD

if __name__ == '__main__':
    sam1path = sys.argv[1]
    sam2path = sys.argv[2]
    afpath = sys.argv[3]
    e1 = float(sys.argv[4]) #probability of het-to-hom or hom-to-het genotyping error, e.g. 0.005
    e2 = float(sys.argv[5]) #probability of hom-to-xhom genotyping error, e.g. 0.0005
    ibdthresh = float(sys.argv[6]) #e.g. 15.0
    xibdthresh = float(sys.argv[7]) #e.g. 15.0 or 5.0
    outpath = sys.argv[8]
    ibdoutpath = sys.argv[9]

    noIndels=True #whether to exclude indels when determinining cmpflag; True setting makes cmpstr='.' for indels; ...this will enable a more apples-to-apples comparison of simple vs allele-frequency-aware (but indel-dumb or at least indel-unimplemented) approach
    noAF0obsalleles=False #we can check and fail when observed alleles have zero allele frequency; if allele frequencies are taken from up-to-date database pool frequencies, then this shouldn't be an issue; otherwise, this could happen naturally in cases of genotyping error
    assert xibdthresh <= ibdthresh
    assert e1>0.0 or e2>0.0 #with modeling of zero error, opposite homozygotes give singularity in llr; the algorithm should handle without crashing but the results could be severely impaired since the opposite homozygote sites would be effectively ignored 

    afdata = bIBD.readAFdata(afpath)

    outfile = open(outpath,'w')
    ibdout = open(ibdoutpath,'w')

    #open sam1 and sam2 and read in tandem; for now, this assumes no header and exactly matching sites, but I should ultimately generalize
    #also, doesn't handle indels at the moment, and only handles one ALT allele (biallelic)
    #also, only handles diploid

    #only single chromosome data is handled for now, but the extension is straightforward
    cllr=0.0
    cllr_max=cllr
    cllr_min_pos=0
    cllr_max_pos=0
    threshTriggered=False
    storeNextPointAsPotentialStart=True
    potentialStart=0

    sam1=open(sam1path,'r')
    sam2=open(sam2path,'r')
    line1=sam1.readline()
    line2=sam2.readline()
    while not line1=='':
        split1=line1.strip().split('\t')
        split2=line2.strip().split('\t')
        rsIDSTR1=split1[0]
        chrSTR1=split1[1]
        pos1=int(split1[2])
        gtSTR1=split1[3]
        rsIDSTR2=split2[0]
        chrSTR2=split2[1]
        pos2=int(split2[2])
        gtSTR2=split2[3]
        assert rsIDSTR1==rsIDSTR2
        assert chrSTR1==chrSTR2
        assert pos1==pos2
        #we output results in the tab-delimited format: rsID, chr, pos, gt1, gt2, cmpflag, ref, alt, alt_af, llr_halfIBD
        cmpflag='.' #A: hom-to-xhom, B: hom-to-hom, C: het-to-het, D: hom-to-het, .: at least one site is no-call
        aaf=float('NaN')
        llr=float('NaN')
        obsals=[]
        ref='.'
        alt='.'

        assert len(gtSTR1)==2 #diploid genotypes (no male chrX)
        assert len(gtSTR2)==2

        #determine cmpflag; this can treat indels
        if '-' not in gtSTR1 and '-' not in gtSTR2: #only determine cmpflag for sites with calls
            obsals=list(set(gtSTR1 + gtSTR2)) #reduce to the observed alleles; cf. https://stackoverflow.com/questions/13902805/list-of-all-unique-characters-in-a-string
            if not (('I' in obsals or 'D' in obsals) and noIndels): #exclude indels if specified (these will get cmpstr='.')
                assert len(obsals) in [1,2] # more than 2 would be inconsistent with biallelic site
                if len(obsals)==1:
                    cmpflag='B' #hom-to-hom
                    assert bIBD.gtIsHom(gtSTR1)
                    assert bIBD.gtIsHom(gtSTR2)
                    assert gtSTR1==gtSTR2
                else: #len(obsals)==2
                    if bIBD.gtIsHom(gtSTR1) and bIBD.gtIsHom(gtSTR2):
                        cmpflag='A' #hom-to-xhom
                        assert not gtSTR1==gtSTR2
                    elif gtSTR1==gtSTR2:
                        cmpflag='C' #het-to-het
                        assert not bIBD.gtIsHom(gtSTR1)
                        assert not bIBD.gtIsHom(gtSTR2)
                    else:
                        cmpflag='D' #hom-to-het
                        assert (not bIBD.gtIsHom(gtSTR1) and bIBD.gtIsHom(gtSTR2)) or (bIBD.gtIsHom(gtSTR1) and not bIBD.gtIsHom(gtSTR2))

        #determine aaf,ref,alt
        if (chrSTR1,pos1) in afdata:
            (rsID,ref,altList,aaflist)=afdata[(chrSTR1,pos1)]
            assert len(ref)==1
            assert len(altList)==1
            assert len(aaflist)==1
            alt=altList[0]
            assert len(alt)==1
            mismatch=False
            for a in obsals:
                if not a in [ref,alt]: mismatch=True
            if mismatch: print 'Allele mismatch: observed alleles: ' + str(obsals) + '; allele-frequency alleles: ' + str([ref,alt])
            else:  aaf=aaflist[0] #we want to check that the observed alleles match ref/alt before assigning aaf

        #determine llr
        if not cmpflag=='.' and not math.isnan(aaf):
            #let B = alt allele; and x=B allele frequency, y=A allele frequency
            #NOTE: to speed things when running multiple comparisons, could precompute LLRs for each possible genotype combination for each site
            x=aaf
            x2=x*x
            y=1.0-x
            y2=y*y
            A_a_a=y2*y #notation: A shared IBD, and the unshared allele is a for each sample
            B_b_b=x2*x
            A_a_b=2*x*y2
            B_a_b=2*x2*y
            ab_ab_hibd=x*y2+x2*y
            aa_bb=2*x2*y2
            aa_aa=y2*y2
            bb_bb=x2*x2
            ab_ab=2*aa_bb #4*x2*y2
            aa_ab=4*x*y2*y
            bb_ab=4*x2*x*y
            hom_hom_e_coef=1.0-2*e1-2*e2
            hom_het_e_coef=1.0-3*e1-e2
            het_het_e_coef=1.0-4*e1
            if cmpflag=='A': #hom-to-xhom, AA vs BB
                if noAF0obsalleles: assert not x==0.0 #note: this situation doesn't necessarily give singularity
                if noAF0obsalleles: assert not y==0.0 #note: this situation doesn't necessarily give singularity
                hibd=e1*(A_a_b+B_a_b) + 2*e2*(A_a_a+B_b_b) #+hom_hom_e_coef*aa_bb_hibd=0
                xibd=hom_hom_e_coef*aa_bb + e1*(aa_ab+bb_ab) + 2*e2*(aa_aa+bb_bb)
            elif cmpflag=='B': #hom-to-hom
                if ref in gtSTR1: #AA vs AA
                    if noAF0obsalleles: assert not y==0.0
                    hibd=hom_hom_e_coef*A_a_a + e1*A_a_b #+e2*aa_bb_hibd=0
                    xibd=hom_hom_e_coef*aa_aa + e1*aa_ab + e2*aa_bb
                else: #alt in gtSTR1 #BB vs BB
                    if noAF0obsalleles: assert not x==0.0
                    hibd=hom_hom_e_coef*B_b_b + e1*B_a_b #+e2*aa_bb_hibd=0
                    xibd=hom_hom_e_coef*bb_bb + e1*bb_ab + e2*aa_bb
            elif cmpflag=='C': #het-to-het, AB vs AB
                if noAF0obsalleles: assert not x==0.0
                if noAF0obsalleles: assert not y==0.0
                hibd=het_het_e_coef*ab_ab_hibd + e1*(A_a_b+B_a_b)
                xibd=het_het_e_coef*ab_ab + e1*(aa_ab+bb_ab)
            elif cmpflag=='D': #hom-to-het
                if noAF0obsalleles: assert not x==0.0 #note: this situation doesn't necessarily give singularity
                if noAF0obsalleles: assert not y==0.0 #note: this situation doesn't necessarily give singularity
                if (ref in gtSTR1) and (ref in gtSTR2): #AA vs AB
                    hibd=hom_het_e_coef*A_a_b + 2*e1*(ab_ab_hibd+A_a_a) + e2*B_a_b #+e1*aa_bb_hibd=0
                    xibd=hom_het_e_coef*aa_ab + e1*(2*ab_ab+2*aa_aa+aa_bb) + e2*bb_ab
                else: #(alt in gtSTR1) and (alt in gtSTR2) #BB vs AB
                    hibd=hom_het_e_coef*B_a_b + 2*e1*(ab_ab_hibd+B_b_b) + e2*A_a_b #+e1*aa_bb_hibd=0
                    xibd=hom_het_e_coef*bb_ab + e1*(2*ab_ab+2*bb_bb+aa_bb) + e2*aa_ab
            #hibd tally:
            #      e1 +: 2*A_a_a, 2*B_b_b, 4*ab_ab_hibd, 3*A_a_b, 3*B_a_b
            #      e1 -: 2*A_a_a, 2*B_b_b, 4*ab_ab_hibd, 3*A_a_b, 3*B_a_b
            #      e2 +: 2*A_a_a, 2*B_b_b, 1*A_a_b, 1*B_a_b
            #      e2 -: 2*A_a_a, 2*B_b_b, 1*A_a_b, 1*B_a_b
            #xibd tally:
            #      e1 +: 2*aa_bb, 2*aa_aa, 2*bb_bb, 4*ab_ab, 3*aa_ab, 3*bb_ab
            #      e1 -: 2*aa_bb, 2*aa_aa, 2*bb_bb, 4*ab_ab, 3*aa_ab, 3*bb_ab
            #      e2 +: 2*aa_bb, 2*aa_aa, 2*bb_bb, 1*aa_ab, 1*bb_ab
            #      e2 -: 2*aa_bb, 2*aa_aa, 2*bb_bb, 1*aa_ab, 1*bb_ab
            llr=math.log(hibd/xibd)

            #the below nan check is probably superfluous in most use cases but could still arise when e.g. e1=e2=0
            if not math.isnan(llr): (cllr,cllr_min_pos,cllr_max,cllr_max_pos,threshTriggered,storeNextPointAsPotentialStart,potentialStart)=bIBD.ibdanalysis(cllr,cllr_min_pos,cllr_max,cllr_max_pos,threshTriggered,storeNextPointAsPotentialStart,potentialStart,ibdthresh,xibdthresh,pos1,llr,ibdout,False)


        outfile.write('\t'.join([rsIDSTR1,chrSTR1,str(pos1),gtSTR1,gtSTR2,cmpflag,ref,alt,str(aaf),str(llr)])+'\n')
        line1=sam1.readline()
        line2=sam2.readline()

    bIBD.ibdanalysis(cllr,cllr_min_pos,cllr_max,cllr_max_pos,threshTriggered,storeNextPointAsPotentialStart,potentialStart,ibdthresh,xibdthresh,pos1,llr,ibdout,True)
    sam1.close()
    sam2.close()
    outfile.close()
    ibdout.close()

def ibdanalysis(cllr,cllr_min_pos,cllr_max,cllr_max_pos,threshTriggered,storeNextPointAsPotentialStart,potentialStart,ibdthresh,xibdthresh,pos,llr,ibdout,chrEnd):

    if not chrEnd: cllr=cllr+llr #we don't add if chrEnd=True to avoid double-counting (this function usually gets called twice on last point, first with chrEnd=False, then with chrEnd=True; the situation where it only gets called once would involve llr=NaN on last site); since NaN should only show up on the last iteration, it should be sufficient to just check for chrEnd rather than also checking for NaN which are assumed to be screened out before call

    if storeNextPointAsPotentialStart:
        potentialStart=pos
    storeNextPointAsPotentialStart=False #this should be reset to True wheneverthe cllr_min_pos is updated (we want the potentialStart to be the following marker)

    if threshTriggered:
        printAndReset=False
        if cllr >= cllr_max:
            cllr_max=cllr
            cllr_max_pos=pos
        elif llr < 0.0:
           ddn=cllr_max - cllr
           if ddn > xibdthresh:
               printAndReset=True
        if printAndReset or chrEnd:
            bIBD.printIBDseg(ibdout, potentialStart, cllr_max_pos, cllr_max) #print the segment
            cllr=0.0 #reset
            cllr_max=0.0
            cllr_max_pos=pos
            cllr_min_pos=pos
            threshTriggered=False
            storeNextPointAsPotentialStart=True
    else: # not threshTriggered
        if llr <0.0:
            if cllr<0.0 or cllr_max - cllr > xibdthresh: #reset if early reversion or new minimum; the (cllr_max-cllr>xibdthresh) check could formulated as an elseif for efficiency, if the compiler would have it evaluated even if cllr<0.0 is True
                cllr=0.0 #reset to zero for better scaling when lower precision arithmetic is used, and also allows elimination of previously-used cllr_min variable and one floating point operation
                cllr_max=0.0 #resetting cllr_max variables is not strictly necessary for cllr<0.0 situation, but it is done for consistency and simplicity
                cllr_max_pos=pos
                cllr_min_pos=pos
                storeNextPointAsPotentialStart=True
        else: # llr>= 0.0
            if cllr >= cllr_max: #this block needed to check for early reversion and also should be done if dup>ibdthresh below gets triggered
                cllr_max=cllr
                cllr_max_pos=pos
                if cllr > ibdthresh:
                    threshTriggered=True

    return (cllr,cllr_min_pos,cllr_max,cllr_max_pos,threshTriggered,storeNextPointAsPotentialStart,potentialStart)

def printIBDseg(out,start,end,segcllr):
    out.write('\t'.join([str(start),str(end),str(segcllr)])+'\n')

def readAFdata(path):
    import gzip
    myDict={}
    if path.endswith('.gz'): ts=gzip.open(path,'r')
    else: ts=open(path,'r')
    count=0
    missingcount=0
    for line in ts:
        if not line.startswith('#'):
            split=line.split('\t')
            info=split[7].split(';')
            acSTR=next(x for x in info if x.startswith('AC=')) #https://stackoverflow.com/questions/9542738/python-find-in-list
            anSTR=next(x for x in info if x.startswith('AN=')) #https://stackoverflow.com/questions/9542738/python-find-in-list
            if not anSTR=='.' and not anSTR=='0':
                count=count+1
                an=int(anSTR[3:])
                ac=[int(x) for x in acSTR[3:].split(',')]
                aaf=[float(x)/float(an) for x in ac]
                myDict[(split[0].strip(),int(split[1].strip()))]=(split[2].strip(),split[3].strip(),split[4].strip().split(','),aaf) #dictionary structured as (chrSTR,pos)=(rsID,REF,[ALT],[ALT_AF])
            else:
                missingcount=missingcount+1
    ts.close()
    print 'Records with missing AF data = ' + str(missingcount)
    print 'Records in AF file = ' + str(count)
    print 'Records in AF dictionary = ' + str(len(myDict))
    return myDict

def gtIsHom(gtstr):
    assert len(gtstr)==2
    return gtstr[0]==gtstr[1]

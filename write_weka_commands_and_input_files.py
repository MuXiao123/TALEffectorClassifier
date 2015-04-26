# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 19:57:34 2014

@author: Katie
"""

from optparse import OptionParser
import csv
import pickle
import os

##command line option parsing
usage = 'usage: %prog [options]'
parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='inFolder', type='string', default="in.txt", help='folder full of output files from target finder')
parser.add_option('-g', '--genomeInfo', dest='genomeInfo', type='string', default="genome.txt", help='filepath prefix for genomic resources')
parser.add_option('-o', '--output', dest='outFolder', type='string', default="out.txt", help='folder in which to place arff file which can be used as input for machine learning classifier')
parser.add_option('-s', '--skip', dest='toSkip', type='int', default=4, help='number of lines to skip in input file')
parser.add_option('-c', '--classifier', dest='classifierLoc', type='string', default="AllFeaturesPlusIdNB.model", help='location of classifier .model file to be used in weka commands output by this script')
parser.add_option('-r', '--resultsFile', dest='resultsFile', type='string', default="out2.txt", help='location to output results from Weka created using the commands output by this script')
parser.add_option('-w', '--wekaPath', dest='wekaPath', type='string', default="weka.jar", help='location of weka.jar')

##read in dictionaries describing various promoter features
(options, args) = parser.parse_args()
fileYPatch = options.genomeInfo+"_Y_patches.pickled"
Ypatches=pickle.load(open(fileYPatch, 'r'))
fileTATABox = options.genomeInfo+"_TATA_boxes.pickled"
TATAboxes=pickle.load(open(fileTATABox, 'r'))
filePromLen = options.genomeInfo+"_PromLen.pickled"
PromLen = pickle.load(open(filePromLen,'r'))
fileTXSLoc = options.genomeInfo+"_TXSLocations.pickled"
TXSLocations = pickle.load(open(fileTXSLoc,'r'))
fileTSSLoc = options.genomeInfo+"_TSSLocations.pickled"
TSSLocations = pickle.load(open(fileTSSLoc,'r'))
fileStrands = options.genomeInfo+"_strands.pickled"
strands = pickle.load(open(fileStrands,'r'))
fileStrands = options.genomeInfo+"_annotatedTXS.pickled"
annotatedTXS = pickle.load(open(fileStrands,'r'))

##open output file where Weka commands will be printed
outScript=open(options.outFolder+"/runWeka.txt","w")

##loop through target finder output files in input folder
for root, directories, files in os.walk(options.inFolder):
        for filename in files:
			##get input file path
            inFile=options.inFolder+"/"+filename
			##get output file path for file where EBE locations will be printed
            EBELocOutFile=options.outFolder+"/EBELocations_"+filename
			##open output files
            EBELocOut = open(EBELocOutFile,"w")
            outFile=open(options.outFolder+"/"+filename.split(".")[0]+".arff","w")
			##print header for file to be used by machine learning classifier
            outFile.write("@RELATION Tals\n\n")
            outFile.write("@ATTRIBUTE ID string\n")
            outFile.write("@ATTRIBUTE RealScore numeric\n")
            outFile.write("@ATTRIBUTE RelativeScore numeric\n")
            outFile.write("@ATTRIBUTE Rank numeric\n")
            outFile.write("@ATTRIBUTE EBEtoTLS numeric\n")
            outFile.write("@ATTRIBUTE EBEtoTXS numeric\n")
            outFile.write("@ATTRIBUTE EBEtoTATA numeric\n")
            outFile.write("@ATTRIBUTE EBEtoYPatch numeric\n")
            outFile.write("@ATTRIBUTE class {Yes,No}\n")
            outFile.write("\n")
            outFile.write("@DATA\n")
            outFile.write("\n")
			##initialize rank at 0, line at 1, and previous score at -1
            rank=0
            line=1
            prevScore=-1.0

            with open(inFile,'r') as f:
                reader=csv.reader(f,delimiter='\t')
				##skip first N lines, specified in options
                for i in range(options.toSkip):
                   s = f.readline()
				   ##grab the best possible score and RVD count from these header lines
                   if(s.find("Best Possible Score:")!=-1):
                       bestScore=float(s.split(":")[1].strip())
                   if(s.find("rvd_sequence")!=-1):
                       talLength = int(s.split("rvd_sequence = ")[1].strip().count("_"))+1
				##for other lines
                for locus,strand,score,start,seq1,seq2 in reader:
					##this assumes that transcript IDs are used in the format gene.isoformNumber
					##gets the gene name
                    gene=locus.split(".")[0]
					##and the isoform number
                    score = float(score)
					##if working with rice, chromosomes Sy and Un have different gene name formats
					## the genes are also not assigned isoform numbers
					## I give them the number two by default because they weren't included in the original
					## ranking used to train the classifier, just like all transcript isoforms not numbered 1
                    if gene == "ChrSy" or gene == "ChrUn":
                        gene = locus.replace("mRNA","gene")
                        isoform = 2
                    else:
                        isoform = int(locus.split(".")[1])
					## get the strand the EBE and gene are each on
                    strandEBE=strand
                    strand=strands[gene]
					##get the location of the EBE in the genome, done differentely depending on the strand of the gene
					##since the promoter is in a different direction depending on the gene's orientation
					## "start" is the distance from the beginning of the promoter to the 5' end (on the plus strand) of the EBE
					## the "EBELoc" being computed is the location of the beginning of the EBE (on the gene strand) relative
					## to the beginning of the chromosome on which the EBE is located
                    if strand=='+':
                        EBELoc = int(start)+TXSLocations[gene][locus]-1000-1
                    if strand=='-':
                        EBELoc = TXSLocations[gene][locus]+1000-int(start)-(talLength-2)
					##print out the gene strand and the EBE location
                    EBELocOut.write(str(locus)+"\t"+strand+"\t"+str(EBELoc)+"\n")
					##get the distance between the EBE and the TXS and TSS of the appropriate transcript
					##also different math depending on the orientation of the gene
					##if the TXS or TSS is downstream of the EBE, the distance will always be positive
					##if the EBE is downstream of the TXS (in the 5' UTR) the distance will always be negative
                    if strand=='+':
                        disToTXS = TXSLocations[gene][locus]-EBELoc
                        disToTSS = TSSLocations[gene][locus]-EBELoc
                    if strand=='-':
                        disToTXS = EBELoc-TXSLocations[gene][locus]
                        disToTSS = EBELoc-TSSLocations[gene][locus]
					##keeping track of rank and line because ties are given the same rank
					##since only EBEs in primary transcripts were counted in the original rank computations
					##used to make the classifier in Cernades et al, so only primary transcripts are counted in this ranking
					##rank is therefore equivalent to a count of EBEs in primary transcripts with a higher score than a given EBE
					##only increase line number if looking at an EBE in a primary transcripts
					##only make rank equivalent to line number when score is greater than the previous score
                    if score > prevScore and isoform == 1:
                        rank=line
                        prevScore=score
                    if isoform == 1: 
                        line+=1
                    if score > prevScore:
                        rank = line
                    relScore=float(score)/bestScore
                    ##the greatest distance a Y patch could be from the EBE is the length of the promoter, so initialize with that
                    minToYPatch=PromLen[locus]
                    ##get the Y patch closest to the EBE
                    for j in range(len(Ypatches[locus].keys())):
                        if strand=='+':
                            newYpatch=abs(int(start)-Ypatches[locus].keys()[j])
                        if strand=='-':
                            newYpatch=abs(Ypatches[locus].keys()[j]-int(start)-(talLength-2))
                        if newYpatch<abs(int(float(minToYPatch))):
                            if strand=='+':
                                minToYPatch=str(int(start)-Ypatches[locus].keys()[j])
                            if strand=='-':
                                minToYPatch=str(Ypatches[locus].keys()[j]-int(start)-(talLength-2))
                    ##if no Y patch, use a question mark which indicates a missing value
                    if len(Ypatches[locus].keys())==0:
                        minToYPatch='?'
                    ##the greatest distance a TATA box could be from the EBE is the length of the promoter, so initialize with that
                    minToTATABox=str(PromLen[locus])
                    ##get the Y patch closest to the EBE
                    for j in range(len(TATAboxes[locus].keys())):
                        if strand=='+':
                            newTATAbox=abs(int(start)-TATAboxes[locus].keys()[j])
                        if strand=='-':
                            newTATAbox=abs(TATAboxes[locus].keys()[j]-int(start)-(talLength-2))
                        if newTATAbox<abs(int(float(minToTATABox))):
                            if strand=='+':
                                minToTATABox=str(int(start)-TATAboxes[locus].keys()[j])
                            if strand=='-':
                                minToTATABox=str(TATAboxes[locus].keys()[j]-int(start)-(talLength-2))
                    if len(TATAboxes[locus].keys())==0:
                        minToTATABox ='?'
                    if annotatedTXS[locus]:    
                        outFile.write(locus+","+str(score)+","+str(relScore)+","+str(rank)+","+str(disToTSS)+","+str(disToTXS)+","+str(minToTATABox)+','+str(minToYPatch)+",?\n")
                    elif not annotatedTXS[locus]:
                        outFile.write(locus+","+str(score)+","+str(relScore)+","+str(rank)+","+str(disToTSS)+",?,"+str(minToTATABox)+','+str(minToYPatch)+",?\n")
            EBELocOut.close()        
            outFile.close()
            outScript.write("java -classpath \""+options.wekaPath+"\" weka.classifiers.bayes.NaiveBayes -l "+options.classifierLoc+" -T "+options.outFolder+"/"+filename.split(".")[0]+".arff -p 1 > \""+options.resultsFile+"/"+filename.split(".")[0]+"ClassifierResults.txt\"\n")
        
        outScript.close()

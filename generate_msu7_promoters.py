# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 09:04:43 2014

@author: Katie Wilkins modified from Nick Booher's generate_msu7_promoters.py script
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from msu_gff_parser import parse_gff
import pickle
from math import fabs
from optparse import OptionParser

##parse command line options
usage = 'usage: %prog [options]'
parser = OptionParser(usage=usage)
parser.add_option('-a', '--annot', dest='annot', type='string', default="msu7_annotation.gff3", help='filepath of genome annotation file')
parser.add_option('-g', '--genome', dest='genomeLoc', type='string', default="Osativa_204_v7.fa", help='filepath of fasta format genome sequence')
parser.add_option('-o', '--outPrefix', dest='outPrefix', type='string', default="msu7", help='filepath prefix of output files')
(options, args) = parser.parse_args()

##read in genome annotation
genome = parse_gff(options.annot)

##read in the genome
sequences = SeqIO.parse(options.genomeLoc, "fasta")
chromosome_seqs = {}
for record in sequences:
	chromosome_seqs[record.id] = record.seq

##create dictionaries to store different features
promoters = {}
TXStoTSS = {}
PromLen = {}
TXSLocations = {}
TSSLocations = {}
strands = {}
annotatedTXS = {}
transcriptNames = list()

for chromosome in genome:
    for gene in genome[chromosome]:
        strands[gene] = genome[chromosome][gene]["strand"]
        TXSLocations[gene] = {}
        TSSLocations[gene] = {}
        for mRNA in genome[chromosome][gene]["mRNA"]:
            transcriptNames.append(mRNA)
            mRNA_feature = genome[chromosome][gene]["mRNA"][mRNA]
            
			##check if a transcriptional start site is annoated
            annotatedTXS[mRNA] = len(mRNA_feature["five_prime_UTR"].keys()) != 0
            
			##could choose to remove mRNAs which are too small, but currently keeping everything
            ##if mRNA_feature["total_cds"] < 150:
                ##continue
            
			##check if the transcript is on the plus strand
            plus_strand = mRNA_feature["strand"] == "+"
            
			##check if any coding sequence is annotated for this transcript
            num_CDSs = len(mRNA_feature["CDS"].keys())
            
			##if there are annotated CDSs, 
            if num_CDSs > 0:
                
                CDSs = mRNA_feature["CDS"]
                
                if plus_strand:
                    ## if on the positive strand, the first CDS in the transcript is the one that starts first the genome
                    TSS = CDSs[min(CDSs, key=lambda x: CDSs[x]["start"])]["start"]
                    FirstCDS=CDSs[min(CDSs, key=lambda x: CDSs[x]["start"])]
                    TSSLocations[gene][mRNA_feature["name"]]=TSS
                    TXS = mRNA_feature["start"]
                    TXSLocations[gene][mRNA_feature["name"]]=TXS
					##promoter is 1000bp upstream of the transcriptional start site plus the 5' UTR
					## -1 at beginning is because python strings are zero indexed
                    promoters[mRNA_feature["name"]] = chromosome_seqs[chromosome][TXS-1000-1 : TSS]
                    TXStoTSS[mRNA_feature["name"]]=fabs(TSS-TXS)
                    PromLen[mRNA_feature["name"]]=len(promoters[mRNA_feature["name"]])
                else:
                    ## if on the positive strand, the first CDS in the transcript is the one that starts last in genome
                    ## remembering that the start is labeled the end site in the annotation when on the negative strand
                    TSS = CDSs[max(CDSs, key=lambda x: CDSs[x]["stop"])]["stop"]
                    FirstCDS=CDSs[max(CDSs, key=lambda x: CDSs[x]["stop"])]
                    TSSLocations[gene][mRNA_feature["name"]]=TSS
                    TXS = mRNA_feature["stop"]
                    TXSLocations[gene][mRNA_feature["name"]]=TXS
					##promoter is 1000bp upstream of the transcriptional start site plus the 5' UTR
					## -1 at beginning is because python strings are zero indexed
                    promoters[mRNA_feature["name"]] = chromosome_seqs[chromosome][TSS-1 : TXS+1000].reverse_complement()
                    TXStoTSS[mRNA_feature["name"]]=fabs(TSS-TXS)
                    PromLen[mRNA_feature["name"]]=len(promoters[mRNA_feature["name"]])
                if PromLen[mRNA_feature["name"]] != TXStoTSS[mRNA_feature["name"]] + 1001:
                    print mRNA_feature["name"] + " doesn't have a complete promoter! Genome assembly is probably incomplete."

##write out promoters in fasta format and write out all genome features in pickle format
SeqIO.write([SeqRecord(promoters[promoter], id=promoter, description="") for promoter in sorted(promoters)], options.outPrefix+"promoters.fasta", "fasta")
outfile = open('_promoters.pickled', 'wb')
pickle.dump(promoters, outfile)
outfile.close()
outfile = open(options.outPrefix+'_annotatedTXS.pickled', 'wb')
pickle.dump(annotatedTXS, outfile)
outfile.close()
outfile = open(options.outPrefix+'_TXStoTSS.pickled', 'wb')
pickle.dump(TXStoTSS, outfile)
outfile.close()
outfile = open(options.outPrefix+'_PromLen.pickled', 'wb')
pickle.dump(PromLen, outfile)
outfile.close()
outfile = open(options.outPrefix+'_TXSLocations.pickled', 'wb')
pickle.dump(TXSLocations, outfile)
outfile.close()
outfile = open(options.outPrefix+'_TSSLocations.pickled', 'wb')
pickle.dump(TSSLocations, outfile)
outfile.close()
outfile = open(options.outPrefix+'_strands.pickled', 'wb')
pickle.dump(strands, outfile)
outfile.close()
outfile = open(options.outPrefix+'_transcript_names.pickled', 'wb')
pickle.dump(transcriptNames, outfile)
outfile.close()

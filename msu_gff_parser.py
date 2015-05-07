# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 08:54:47 2014

@author: Nick
"""
from collections import defaultdict
from urllib import unquote

def parse_gff(gff_filepath, version=7):

    with open(gff_filepath, "rU") as gff_file:
        
        genome = defaultdict(dict)
        
        gene_id = None
        mRNA_id = None
        
        for line in gff_file:
            
            line = line.rstrip()
            
            if len(line) == 0 or line[0] == "#":
                continue
            
            chromosome, _, feature_type, start, stop, _, strand, _, annotation_str = line.rstrip().split("\t")
            
            if version == 6 and (chromosome == "ChrSy" or chromosome == "ChrUn"):
                continue
            
            annotation = {pair.split("=")[0].lower(): pair.split("=")[1] for pair in annotation_str.split(";")}
            
            if feature_type == "gene":
                
                gene_id = annotation["id"]
                
                genome[chromosome][gene_id] = {
                    "start": int(start),
                    "stop": int(stop),
                    "strand": strand,
                    "mRNA": {}
                }
                
                if version == 6:
                    genome[chromosome][gene_id]["name"] = annotation["alias"]
                    genome[chromosome][gene_id]["annotation"] = unquote(annotation["name"])
                else:
                    genome[chromosome][gene_id]["name"] = annotation["name"]
                    genome[chromosome][gene_id]["annotation"] = unquote(annotation["note"])
                
            elif feature_type == "mRNA":
                
                mRNA_id = annotation["id"]
                
                genome[chromosome][gene_id]["mRNA"][mRNA_id] = {
                    "start": int(start),
                    "stop": int(stop),
                    "strand": strand,
                    "total_cds": 0,
                    "five_prime_UTR": {},
                    "exon": {},
                    "CDS": {},
                    "three_prime_UTR": {}
                }
                
                if version == 6:
                    genome[chromosome][gene_id]["mRNA"][mRNA_id]["name"] = annotation["alias"]
                else:
                    genome[chromosome][gene_id]["mRNA"][mRNA_id]["name"] = annotation["name"]
                
            else:
                
                feature = {
                    "start": int(start),
                    "stop": int(stop),
                    "strand": strand,
                }
                
                feature_id = annotation.get("id", id(feature))
                
                genome[chromosome][gene_id]["mRNA"][mRNA_id][feature_type][feature_id] = feature
                
                if feature_type == "CDS":
                    genome[chromosome][gene_id]["mRNA"][mRNA_id]["total_cds"] += (feature["stop"] - feature["start"] + 1)
    
    return genome

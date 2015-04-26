#! usr/bin/python

#TATA_and_TLS_info_step5.py

#Description: A script to determine the distance of each TAL binding site (best sites in each promoter and 
#those below cutoff score) from a TATA box (if any) and the transcriptional start site)

#Created by Erin Doyle on Tuesday, Jan 11 2011

#Feb. 4, 2011: TATA boxes modified to reflect rice-specific frequencies in 
#Civan and Svec, Genome 53(3), pp.294-291. 2009. TATA consensus: TATAWAWA
#Feb. 4, 2011: Y-patch information added from the same paper. consensus: CYYCYYCCYC
###################################################################################################
# Section 0: import modules, scoring and cutoff data, and promoters and utr sequences
###################################################################################################
import pickle
import string

print 'Loading promoters...'
promoters = pickle.load(open('/home/kxw116/GenomicResources/OsativaMSU7/msu7_promoters.pickled', 'r'))

###################################################################################################
# Section 2: Locate TATA boxes in promoters
###################################################################################################
#TATA_sequences = ['TATAAA', 'CTATAAATAC', 'TATAAAT', 'TATTAAT', 'TATATAA', 'TTATTT']
TATA_sequences = ['TATAAAT', 'TATATAT', 'TATATAA', 'TATAAAA']

TATAs = {}
for locus in promoters.keys():
	
	TATAs[locus] = {}
	for TATA in TATA_sequences:
		##split the promoter string everywhere a TATA sequence is found
		split_promoter = string.split(promoters[locus], TATA)

		position = 0 #Track the position of the TATA boxes in the promoter
		for sequence in split_promoter:
			position += len(sequence)
			
			if position != len(promoters[locus]): #length will be length of the whole sequence if no TATA box is found
				TATAs[locus][position] = TATA
			position += len(TATA)

outfile = open('/home/kxw116/GenomicResources/OsativaMSU7/msu7_TATA_boxes.pickled', 'w')
pickle.dump(TATAs, outfile)
outfile.close()

###################################################################################################
# Section 2: Locate Y patch elements in promoters
###################################################################################################
y_patch_sequences = ['CCTCCCCCTC', 'CTTCCCCCTC', 'CCCCCCCCTC', 'CTCCCCCCTC', 'CCTCTCCCTC', 'CTTCTCCCTC', 'CCCCTCCCTC', 'CTCCTCCCTC', 'CCTCCTCCTC', 'CTTCCTCCTC', 'CCCCCTCCTC', 'CTCCCTCCTC', 'CCTCTTCCTC', 'CTTCTTCCTC', 'CCCCTTCCTC', 'CTCCTTCCTC', 'CCTCCCCCCC', 'CTTCCCCCCC', 'CCCCCCCCCC', 'CTCCCCCCCC', 'CCTCTCCCCC', 'CTTCTCCCCC', 'CCCCTCCCCC', 'CTCCTCCCCC', 'CCTCCTCCCC', 'CTTCCTCCCC', 'CCCCCTCCCC', 'CTCCCTCCCC', 'CCTCTTCCCC', 'CTTCTTCCCC', 'CCCCTTCCCC', 'CTCCTTCCCC']

y_patches = {}
for locus in promoters.keys():
	
	y_patches[locus] = {}
	for y_patch in y_patch_sequences:
		##split the promoter string everywhere a Y patch sequence is found
		split_promoter = string.split(promoters[locus], y_patch)

		position = 0 #Track the position of the y patches in the promoter
		for sequence in split_promoter:
			position += len(sequence)
			
			if position != len(promoters[locus]): #length will be length of the whole sequence if no y patch is found
				y_patches[locus][position] = y_patch
			position += len(y_patch)

outfile = open('/home/kxw116/GenomicResources/OsativaMSU7/msu7_Y_patches.pickled', 'w')
pickle.dump(y_patches, outfile)
outfile.close()

# Filename: calc_ssbonds.py
# Author: Evelyne Deplazes
# Date: May 5, 2018 

# Script to calculate the distance between the Carbon alpha (CA) in a 
# disulfide bond (CA-CA distance)

# the script relies on the definition of the disulfdie bonds by the keyword SSBOND 
# in the HEADER of the PDB file (example below) 

# SSBOND   1 CYS A    3    CYS A   20                          1555   1555  2.03  
# SSBOND   2 CYS A    7    CYS A   16                          1555   1555  2.03  

# the script reqires the standart python libraries umpy, sys, math and glob
# as well as the library  MDAnalysis https://www.mdanalysis.org 

# The script loops through all .pdb files in a the current folder.  
# For each .pdb file the script calculates the average CA-CA 
# distance in EACH SSBOND. 
# for PDB files with multiple models, the CA-CA distances are 
# averaged over all models. 
# the script produces the following output files:
# a .dat file for each .pdb file containing 
# the avg +/- stdev CA-CA distance for each SSBOND 
# a file called ssbonds_histogram.dat that contains a 
# histogram of the avg CA-CA distances from all .pdb files 

import glob
import re
import sys, os
import math as math
import numpy as np
import MDAnalysis as MDA

import MDAnalysis.analysis.distances as distances 
ssbonds_all = []
f2=str("ssbonds_histogram.dat")
outfile2=open(f2, 'w')

# loop through all .pdb files in the current folder 
for pdb_file in glob.glob('*.pdb'):
    print "processing pdb file", pdb_file, "................"
    f1=str("ssbonds_"+str(pdb_file)+".dat")
    outfile=open(f1, 'w')
    string=str(str(pdb_file)+"\n") 
    outfile.write(string)
    # find number of ss-bonds in the structure 
    # by looking for the SSBOND keyword in the .pdb file 
    # for each ssbond, get residue numbers involved 
    p = re.compile('(SSBOND)[a-z]*')
    count_ssbonds = 0
    resid1=[]
    resid2=[]
    infile=open(pdb_file,'r')  #open pdb file
    line=infile.readline()
    while (line):
        if p.match(line):
            count_ssbonds = count_ssbonds + 1
            resid1.append(int(line[19:21]))
            resid2.append(int(line[33:35]))
        line =infile.readline()
    
    # load the pdb file with MDAnalysis 
    # select protein and find number of frames/models
    u = MDA.Universe(pdb_file)
    peptide=u.select_atoms("protein").residues
    num_models = len(u.trajectory)

    # for each ssbond, loop through pdb models and
    # calculate average CA-CA distance 
    for ssbond in range(0,count_ssbonds):
        # select CA in ssbonds
        string=str(str("avg and stdev for CA-CA distance in ssbond between resids ")+str(resid1[ssbond])+str(" and ")+str(resid2[ssbond])+"\n")
        outfile.write(string)        
        CA_CA_dist = [] 
        CA_resid1 = u.select_atoms("resid "+str(resid1[ssbond])+" and name CA")
        CA_resid2 = u.select_atoms("resid "+str(resid2[ssbond])+" and name CA")
        if (len(distances.dist(CA_resid1, CA_resid2)[2]) > 1):
            print "ss bond excluded due to multiple chains", pdb_file
        else:
            for ts in u.trajectory[0:num_models:1]:
                CA_CA_dist.append(distances.dist(CA_resid1, CA_resid2)[2])
                ssbonds_all.append(distances.dist(CA_resid1, CA_resid2)[2])
            string=str(str(np.mean(CA_CA_dist))+" +/- "+str(np.std(CA_CA_dist))+"\n")
            outfile.write(string)
            if (np.mean(CA_CA_dist) < 4.20):
                print pdb_file
        #print np.mean(CA_CA_dist)
        #print np.std(CA_CA_dist)
        
print "collected data for", len(ssbonds_all), "ss-bonds from ", len(glob.glob('*.pdb')), "pdb files"

for i in range(0,len(ssbonds_all)):
    if (len(ssbonds_all[i]) > 1):
        print len(ssbonds_all[i])

hist_range_min=3.0
hist_range_max=7.0

# get histogram for ssbonds_all
hist, bin_edges = np.histogram(ssbonds_all, range=[hist_range_min, hist_range_max], bins=20, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile2.write(hist_string)	



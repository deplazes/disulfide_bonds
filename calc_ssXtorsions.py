# Filename: calc_ssXtorisons.py
# Author: Evelyne Deplazes
# Date: May 5, 2018 

# Script to calculate the five torsion angles that are defined by the 
# two residues that form a disulfide bbond (see below for defintion of 
# torsion angles)

# the script relies on the definition of the disulfdie bonds by the keyword SSBOND 
# in the HEADER of the PDB file (example below) 

# SSBOND   1 CYS A    3    CYS A   20                          1555   1555  2.03  
# SSBOND   2 CYS A    7    CYS A   16                          1555   1555  2.03  

# the script reqires the standart python libraries umpy, sys, math and glob
# as well as the library  MDAnalysis https://www.mdanalysis.org 

# The script loops through all .pdb files in a the current folder.  
# For each .pdb file the script calculates the average for each 
# of the five torsion angles. 
# for PDB files with multiple models, the CA-CA distances are 
# averaged over all models. 
# the script produces the following output files:
# histograms for each of the five torsion angles calculated from all .pdb files 

# ----- defintion of torsion angles -----

# the CYS residues that form the SSbond define 5 torsion angles 
# formed by the N CA CB and S atoms of residues Cys i and Cys j 

# tosions and atoms involved (based on schematic in Fig 1 in 
# Sirinivasan et al, J Peptide Research, 1990 

# Xss = CBi Si Sj CBj 
# Xi1 = Ni CAi CBi Si
# Xj1 = Nj CAj CBj Sj
# Xi2 = CAi CBi Si Sj
# Xj2 = CAj CBj Sj Si


import glob
import re
import numpy as np
import sys, os
import math as math
import numpy as np
import MDAnalysis as MDA

import MDAnalysis.analysis.distances as distances 
Xss_all = []
Xi1_all = []
Xi2_all = []
Xj1_all = []
Xj2_all = []

f1=str("ssXtorsions_Xss_hist.dat")
outfile1=open(f1, 'w')
f2=str("ssXtorsions_Xi1_hist.dat")
outfile2=open(f2, 'w')
f3=str("ssXtorsions_Xj1_hist.dat")
outfile3=open(f3, 'w')
f4=str("ssXtorsions_Xi2_hist.dat")
outfile4=open(f4, 'w')
f5=str("ssXtorsions_Xj2_hist.dat")
outfile5=open(f5, 'w')

for pdb_file in glob.glob('*.pdb'):

    print "processing pdb file", pdb_file, "................"                    
    multiple_chains = 0 
    # check if file has mulitiple chains, if so, skip this pdb file 
    p1 = re.compile('(COMPND)[a-z]*')
    infile=open(pdb_file,'r')  #open pdb file
    line=infile.readline()
    while (line):
        if p1.match(line):
            #print line
            if line.find('B;') > 0:
                print "found more than one chain in pdb file"
                multiple_chains=1
        line =infile.readline()
    if multiple_chains == 1:
        print "pdb file excluded due to multiple chains", pdb_file
    else:
        #print "calculating ssbond torsion angle for pdb file", pdb_file, "................"                
        #sys.exit()
        # process pdb file
        # find number of ss-bonds in the structure 
        # by looking for the SSBOND keyword in the .pdb file 
        # for each ssbond, get residue numbers involved 
        p = re.compile('(SSBOND)[a-z]*')
        count_ssbonds = 0
        resid_i=[]
        resid_j=[]
        infile=open(pdb_file,'r')  #open pdb file
        line=infile.readline()
        while (line):
            if p.match(line):
                count_ssbonds = count_ssbonds + 1
                resid_i.append(int(line[19:21]))
                resid_j.append(int(line[33:35]))
            line =infile.readline()
        #print pdb_file, count_ssbonds 
        #print resid1, resid2

        # load the pdb file with MDAnalysis 
        # select protein and find number of frames/models
        u = MDA.Universe(pdb_file)
        peptide=u.select_atoms("protein").residues
        num_models = len(u.trajectory)

        # for each ssbond, loop through pdb models and
        # calculate the torsion angles 
        for ssbond in range(0,count_ssbonds):
            # select CA in ssbonds  
            #print "------------", ssbond
            Xss = []      
            Xi1 = []
            Xj1 = []
            Xi2 = []
            Xj2 = []
            selection = 'resid {} and name {}'
            CAi = selection.format(resid_i[ssbond], 'CA')
            CAj = selection.format(resid_j[ssbond], 'CA')
            CBi = selection.format(resid_i[ssbond], 'CB')
            CBj = selection.format(resid_j[ssbond], 'CB')
            Si = selection.format(resid_i[ssbond], 'SG')
            Sj = selection.format(resid_j[ssbond], 'SG')
            Ni = selection.format(resid_i[ssbond], 'N')
            Nj = selection.format(resid_j[ssbond], 'N')
            # Xss = CBi Si Sj CBj     
            Xss_atoms = [CBi, Si, Sj, CBj]   
            print Xss_atoms
            # Xi1 = Ni CAi CBi Si
            Xi1_atoms = [Ni, CAi, CBi, Si]
            # Xj1 = Nj CAj CBj Sj
            Xj1_atoms = [Nj, CAj, CBj, Sj]
            # Xi2 = CAi CBi Si Sj
            Xi2_atoms = [CAi, CBi, Si, Sj]
            # Xj2 = CAj CBj Sj Si
            Xj2_atoms = [CAj, CBj, Sj, Si]
            Xss_angle = sum([u.select_atoms(atom) for atom in Xss_atoms])  # sum of Atoms creates an AtomGroup
            Xss_angle = Xss_angle.dihedral  # convert AtomGroup to Dihedral object
            Xi1_angle = sum([u.select_atoms(atom) for atom in Xi1_atoms])  
            Xi1_angle = Xi1_angle.dihedral         
            Xj1_angle = sum([u.select_atoms(atom) for atom in Xj1_atoms])  
            Xj1_angle = Xj1_angle.dihedral  
            Xi2_angle = sum([u.select_atoms(atom) for atom in Xi2_atoms])  
            Xi2_angle = Xi2_angle.dihedral  
            Xj2_angle = sum([u.select_atoms(atom) for atom in Xj2_atoms])  
            Xj2_angle = Xj2_angle.dihedral  
            for ts in u.trajectory[0:num_models:1]:
                Xss.append(Xss_angle.value())
                Xss_all.append(Xss_angle.value())
                Xi1.append(Xi1_angle.value())
                Xi1_all.append(Xi1_angle.value())
                Xj1.append(Xj1_angle.value())
                Xj1_all.append(Xj1_angle.value())
                Xi2.append(Xi2_angle.value())
                Xi2_all.append(Xi2_angle.value())
                Xj2.append(Xj2_angle.value())
                Xj2_all.append(Xj2_angle.value())


print "collected data on torsion angles for", len(Xss_all), "ss-bonds from ", len(glob.glob('*.pdb')), "pdb files"
# get histogram for Xss_all
hist, bin_edges = np.histogram(Xss_all, bins=50, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile1.write(hist_string)	

outfile1.close()

hist, bin_edges = np.histogram(Xi1_all, bins=50, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile2.write(hist_string)	

outfile2.close()


hist, bin_edges = np.histogram(Xj1_all, bins=50, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile3.write(hist_string)	

outfile3.close()


hist, bin_edges = np.histogram(Xi2_all, bins=50, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile4.write(hist_string)	

outfile4.close()

hist, bin_edges = np.histogram(Xj2_all, bins=50, density=True)
for i in range(0,len(hist)):
    hist_string=str(str(bin_edges[i])+" "+str(hist[i])+"\n")
    outfile5.write(hist_string)	

outfile5.close()

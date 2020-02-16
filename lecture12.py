#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: roberto.villegas-diaz
"""

# Count specific amino-acid
done=0
while (not done):
    sequence=input("Please enter a sequence: ")
    aa=input("Please enter the amino acid to look for: ")
    
    #compute the number of occurrences using for loop    
    cnt=0
    for i in sequence:
        if i == aa:
            cnt+=1
    if cnt == 1:
        print ("%s occurs in that sequence once" % aa)
    else:
        print ("%s occurs in that sequence %d times" % (aa, cnt))
    answer=input("try again? [yn]")
    if answer == "n" or answer == "N":
        done = 1


# Processing tblastn output
countHits = 0
with open("data/my_tblastn_output_nr.txt","r") as tblastn:
    for hit in tblastn.readlines():
        hit = hit.split('\t')
        countHits += 1
        print("Sequence ID: %s" % hit[4])
        print("e-value: %s " % hit[7])
        print("Sequence ID: %s \ne-value: %s \n" %(hit[4],hit[7]))
print("Number of hits: %s" % countHits)
tblastn.close()

# Read FASTA file
SeqID = ""
A = 0
C = 0
G = 0
T = 0
with open("data/exampledna.fasta","r") as fasta:
    for line in fasta.readlines():
        if(line.startswith(">")):
            SeqID = line.strip()
        else:
            line = line.upper()
            A += sum([int(i == 'A') for i in line]) 
            C += sum([int(i == 'C') for i in line]) 
            G += sum([int(i == 'G') for i in line]) 
            T += sum([int(i == 'T') for i in line])
GCcontent = (C + G)/(A + C + G + T)*100
print("Sequence ID: %s" % SeqID)
print("A: %s \nC: %s \nG: %s \nT: %s" % (A,C,G,T))
print("GC content: %.2f" % GCcontent)
fasta.close()
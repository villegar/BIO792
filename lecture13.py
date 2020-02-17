#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 16:10:22 2020

@author: roberto.villegas-diaz
"""

# Checking version of Biopythib
import Bio
print(Bio.__version__)


# Working with sequences
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
my_seq
print(my_seq)
my_seq.alphabet
my_seq.complement() 
my_seq.reverse_complement()


# Simple FASTA parsing example
from Bio import SeqIO 
for seq_record in SeqIO.parse("data/ls_orchid.fasta", "fasta"): 
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


# Simple GenBank parsing example
from Bio import SeqIO
for seq_record in SeqIO.parse("data/ls_orchid.gbk", "genbank"): 
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


# More sequences: Alphabets
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
my_seq
my_seq.alphabet


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_prot = Seq("AGTACACTGGT", IUPAC.protein) 
my_prot
my_prot.alphabet


# Sequences act like strings
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GATCG", IUPAC.unambiguous_dna)
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter)) 
print(len(my_seq))

from Bio.Seq import Seq 
"AAAA".count("AA")
Seq("AAAA").count("AA")

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna) 
len(my_seq)
my_seq.count("G")
100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna) 
GC(my_seq)


# Slicing a sequence
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna) 
my_seq[4:12] # Starting from 4 to 12

my_seq[0::3] # Starting from 0 every 3 elements
my_seq[1::3] # Starting from 1 every 3 elements
my_seq[2::3] # Starting from 2 every 3 elements

my_seq[::-1] # Reverse the sequence


# Turning Seq objects into strings
str(my_seq) 
print(my_seq)
fasta_format_string = ">Name\n%s\n" % my_seq 
print(fasta_format_string)


# Concatenating or adding sequences
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
protein_seq = Seq("EVRNAK", IUPAC.protein)
dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
#protein_seq + dna_seq
## Error expected

from Bio.Alphabet import generic_alphabet 
protein_seq.alphabet = generic_alphabet 
dna_seq.alphabet = generic_alphabet
protein_seq + dna_seq

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
list_of_seqs = [Seq("ACGT", generic_dna), Seq("AACC", generic_dna), Seq("GGTT", generic_dna)] 
sum(list_of_seqs, Seq("", generic_dna))


# Changing case
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
dna_seq = Seq("acgtACGT", generic_dna) 
dna_seq
dna_seq.upper()
dna_seq.lower()


# Transcription
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna) 
coding_dna
messenger_rna = coding_dna.transcribe()
messenger_rna


# Translation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna) 
coding_dna
coding_dna.translate()


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
    "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
    "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
    "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
    "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
    generic_dna)
gene.translate(table="Bacterial")

gene.translate(table="Bacterial", to_stop=True) 


# Comparing sequences
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
seq1 = Seq("ACGT", IUPAC.unambiguous_dna) 
seq2 = Seq("ACGT", IUPAC.ambiguous_dna) 
str(seq1) == str(seq2)
str(seq1) == str(seq1)

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
dna_seq = Seq("ACGT", generic_dna)
prot_seq = Seq("ACGT", generic_protein)
dna_seq == prot_seq
## BiopythonWarning: Incompatible alphabets DNAAlphabet() and ProteinAlphabet()


# MutableSeq objects 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)

#my_seq[5] = "G"
## Error expected


mutable_seq = my_seq.tomutable()
mutable_seq
new_seq = mutable_seq.toseq()
new_seq

from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)

mutable_seq
mutable_seq[5] = "C"
mutable_seq
mutable_seq.remove("T")
mutable_seq
mutable_seq.reverse()
mutable_seq

# UnknownSeq objects
from Bio.Seq import UnknownSeq
unk = UnknownSeq(20)
unk
print(unk)
len(unk)


from Bio.Seq import UnknownSeq
from Bio.Alphabet import IUPAC
unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna) 
unk_dna
print(unk_dna)

unk_protein = unk_dna.translate()
unk_protein

# Connecting with biological databases
﻿import math
import random
DNA_NUC = ['A', 'T', 'G', 'C']
DNA_Seq = []
RNA_Seq = []
# DNA Replication
Nucleotide = 0
while Nucleotide < 100:
    Nucleotide += 1
    DNA_Seq.append(random.choice(DNA_NUC))
print('DNA Sequence:')
print(DNA_Seq)
print(f'this sequence has {len(DNA_Seq)} bp')
# Transcription

for element in DNA_Seq:
    if element == 'A':
        RNA_Seq.append('U')
    elif element == 'T':
        RNA_Seq.append('A')
    elif element == 'G':
        RNA_Seq.append('C')
    elif element == 'C':
        RNA_Seq.append('G')
RNA_Seq=''.join(RNA_Seq)
print('RNA Sequence:')
print(RNA_Seq)
print(f'this sequence also has {len(RNA_Seq)} bp')



# Translations
mRNA_codons = []

#mRNA_codons = [RNA_Seq[i:i+3] for i in range(0, len(RNA_Seq), 3)]
mRNA_codons = [RNA_Seq[i:i+3] for i in range(0, len(RNA_Seq)-2, 3)]
print(len(mRNA_codons),"codons\n",mRNA_codons)
print(f'The mRNA_codons list has {len(mRNA_codons)} indexes')

protein_sequence = []
for element in mRNA_codons:
    if element in ['UUU', 'UUC']:
        protein_sequence.append('Phe')
    elif element in ['UUA', 'UUG']:
        protein_sequence.append('Leu')
    elif element in ['UCU', 'UCC', 'UCA', 'UCG']:
        protein_sequence.append('Ser')
    elif element in ['UAU', 'UAC']:
        protein_sequence.append('Tyr')
    elif element in ['UGU', 'UGC']:
        protein_sequence.append('Cys')
    elif element in ['UGG']:
        protein_sequence.append('Trp')
    elif element in ['CUU', 'CUC', 'CUA', 'CUG']:
        protein_sequence.append('Leu')
    elif element in ['CCU', 'CCC', 'CCA', 'CCG']:
        protein_sequence.append('Pro')
    elif element in ['CAU', 'CAC']:
        protein_sequence.append('His')
    elif element in ['CAA', 'CAG']:
        protein_sequence.append('Gln')
    elif element in ['CGU', 'CGC', 'CGA', 'CGG']:
        protein_sequence.append('Arg')
    elif element in ['AUU', 'AUC', 'AUA']:
        protein_sequence.append('Ile')
    elif element in ['AUG']:
        protein_sequence.append('Met')
    elif element in ['ACU', 'ACC', 'ACA', 'ACG']:
        protein_sequence.append('Thr')
    elif element in ['AAU', 'AAC']:
        protein_sequence.append('Asn')
    elif element in ['AAA', 'AAG']:
        protein_sequence.append('Lys')
    elif element in ['AGU', 'AGC']:
        protein_sequence.append('Ser')
    elif element in ['AGA', 'AGG']:
        protein_sequence.append('Arg')
    elif element in ['GUU', 'GUC', 'GUA', 'GUG']:
        protein_sequence.append('Val')
    elif element in ['GCU', 'GCC', 'GCA', 'GCG']:
        protein_sequence.append('Ala')
    elif element in ['GAU', 'GAC']:
        protein_sequence.append('Asp')
    elif element in ['GAA', 'GAG']:
        protein_sequence.append('Glu')
    elif element in ['GGU', 'GGC', 'GGA', 'GGG']:
        protein_sequence.append('Gly')

print('here is your protein sequence\n')
for index in range(len(protein_sequence)):
    print(f"{index+1}: {protein_sequence[index]}")

print('protein sequence:', protein_sequence)





def translate_dna_to_protein(dna_codons):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

translate_dna_to_protein(mRNA_codons)

print("new shit")






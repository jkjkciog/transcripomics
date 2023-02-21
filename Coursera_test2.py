import Bio
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
highest_record = 0
i=0



def find_orf1(sequence):                            
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            orf = codon
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    return orf
                orf += codon
            return find_orf1 
    return ""


def find_orf2(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(1, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            orf = codon
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    return orf
                orf += codon
            return find_orf2
    return ""



def find_orf3(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(2, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            orf = codon
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    return orf
                orf += codon
            return orf
    return ""



for seq_record in SeqIO.parse("dna.example (1).fasta", "fasta"):
    print("record number",i)
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    i=i+1
    if len(seq_record) > highest_record:
        highest_record = len(seq_record) 
    orf = find_orf(str(seq_record))
    if len(orf) > 0:
        print("ORF:", orf)
    else:
        print("No ORF found")
print(i,"records")
print("the highest record is",highest_record)































# Define the BLAST search parameters
#database = "nt"
#program = "blastn"
#email = "your_email@example.com"

# Perform the BLAST search
#result_handle = NCBIWWW.qblast(program, database, sequence, email=email)

# Parse the BLAST search results
#from Bio.Blast import NCBIXML
#blast_records = NCBIXML.parse(result_handle)
#blast_record = next(blast_records)

# Print the top hit's organism
#print("Organism:", blast_record.alignments[0].hit_def.split("[")[-1].split("]")[0])

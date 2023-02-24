import Bio
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
highest_record = 0
lowest_record = 1e99
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


def find_orf4(sequence):
    start = sequence.find("ATG")

    if start > -1:
        stopTaa = sequence[start:].find("TAA")
        stopTag = sequence[start:].find("TAG")
        stopTga = sequence[start:].find("TGA")

        stop = -1

        if stopTaa > -1:
            stop = stopTaa
        if stopTag > -1 and stopTag < stop:
            stop = stopTag
        if stopTga > -1 and stopTga < stop:
            stop = stopTga

        if stop > -1:
            return sequence[start:stop], stop

    return "", 0



for seq_record in SeqIO.parse("dna2.fasta", "fasta"):
    i=i+1
    print("record number",i)
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    stop = 0


    
    if len(seq_record) > highest_record:
        highest_record = len(seq_record)
        highest_record_num=i
    if len(seq_record) < lowest_record:
        lowest_record = len(seq_record)
        lowest_record_num=i
    orf, stop = find_orf4(str(seq_record)[stop:])
    
    #next seach in loop needs to pass str(seq_record)[stop:]

    if len(orf) > 0:
        print("ORF:", orf)
    else:
        print("No ORF found")
print(i,"records")
print("the longest record is number",highest_record_num,highest_record)
print("the shortest record is", lowest_record_num, lowest_record)































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

def readGenome(filename):
    genome = ''
    introns = []
    current_sequence = ''

    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()  # Remove trailing newline character
            
            if line.startswith('>'):
                if genome == '':  # If genome is not set yet
                    genome = current_sequence
                elif current_sequence != '':  # Otherwise, if there's any sequence collected
                    introns.append(current_sequence)
                current_sequence = ''  # Reset the sequence collector
            else:
                current_sequence += line

        # After the loop ends, check if there's any sequence left
        if current_sequence:
            introns.append(current_sequence)

    return genome, introns

genome, introns = readGenome('Rosalind_splc.fa.txt')

print("Genome:", genome)
print("Introns:", introns)





# %%
def remove_introns(genome, introns):
    for intron in introns:
        if intron in genome:
            genome = genome.replace(intron, "",1)
    return genome



DNA_Seq=remove_introns(genome, introns)
print(DNA_Seq)


# %%
def DNA_To_RNA(DNA_seq):
    RNA_Seq=[]
    for element in DNA_Seq:
        if element == 'A':
            RNA_Seq.append('A')
        elif element == 'T':
            RNA_Seq.append('U')
        elif element == 'G':
            RNA_Seq.append('G')
        elif element == 'C':
            RNA_Seq.append('C')
    RNA_Seq=''.join(RNA_Seq)
    return RNA_Seq

RNA_Seq=DNA_To_RNA(DNA_Seq)
print(RNA_Seq)


# %%
def translate_dna_to_protein(dna_codons):
    codon_table = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
        'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
    }

    protein_sequence = ""
    for codon in range(0, len(dna_codons), 3):
        protein = codon_table.get(dna_codons[codon:codon+3])
        if protein:
            protein_sequence += protein

    return protein_sequence




protein_seq = translate_dna_to_protein(RNA_Seq)
print("Protein sequence:")
print(protein_seq)



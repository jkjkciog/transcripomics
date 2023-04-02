with open('Homo_Sapiens_Chrom_1.fasta', 'r') as f:
    f.readline() # skip the first line (header)
    genome = f.read().replace('\n', '') # concatenate the rest of the lines

    GC_base=0
for i in range(len(genome)):
    if 'G' or 'C' in genome:
        GC_base+=1
    else:
        continue
print(GC_base)
import Bio
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq


# Define the DNA sequence
sequence = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")

# Define the BLAST search parameters
database = "nt"
program = "blastn"
email = "your_email@example.com"

# Perform the BLAST search
result_handle = NCBIWWW.qblast(program, database, sequence, email=email)

# Parse the BLAST search results
from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

# Print the top hit's organism
print("Organism:", blast_record.alignments[0].hit_def.split("[")[-1].split("]")[0])

from Bio import SeqIO

records = SeqIO.parse("THIS_IS_YOUR_INPUT_FILE.tab", "tab")
count = SeqIO.write(records, "THIS_IS_YOUR_OUTPUT_FILE.fasta", "fasta")
print("Converted %i records" % count)
from Bio import SeqIO

def find_orf(sequence, start_pos, stop_codons):
    start_codon = "ATG"
    for i in range(start_pos, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            orf = codon
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    return orf, j
                orf += codon
            return orf, len(sequence)
    return "", start_pos

def find_orf4(sequence):
    start = sequence.find("ATG")
    if start > -1:
        stop_codons = ["TAA", "TAG", "TGA"]
        stop = len(sequence)
        for stop_codon in stop_codons:
            temp_stop = sequence.find(stop_codon, start)
            if temp_stop > -1 and temp_stop < stop:
                stop = temp_stop
        return sequence[start:stop], stop
    return "", 0

def process_record(seq_record):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    orf, stop = find_orf4(str(seq_record))
    while stop < len(seq_record):
        orf, stop = find_orf(str(seq_record), stop, ["TAA", "TAG", "TGA"])
        if len(orf) > 0:
            print("ORF:", orf)
        else:
            print("No ORF found")

    print()

records = list(SeqIO.parse("dna2.fasta", "fasta"))
longest_record = max(records, key=lambda r: len(r.seq))
shortest_record = min(records, key=lambda r: len(r.seq))
print("Number of records: ", len(records))
print("ID of longest record: ", longest_record.id)
print("Length of longest record: ", len(longest_record))
print("ID of shortest record: ", shortest_record.id)
print("Length of shortest record: ", len(shortest_record))

for seq_record in records:
    process_record(seq_record)
import Bio
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO

def find_stop_codon(seq, start, stop_codons):
    """Find the next stop codon in the given sequence."""
    for i in range(start, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stop_codons:
            return i
    return None

def find_orf(sequence, start_codon="ATG", stop_codons=["TAA", "TAG", "TGA"]):
    """Find the longest ORF in the given sequence."""
    longest_orf = ""
    for i in range(len(sequence)):
        if sequence[i:i+3] == start_codon:
            stop = find_stop_codon(sequence, i+3, stop_codons)
            if stop is not None:
                orf = sequence[i:stop]
                if len(orf) > len(longest_orf):
                    longest_orf = orf
    return longest_orf


def find_orf2(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for i in range(len(sequence)):
        if sequence[i:i+3] == start_codon:
            for j in range(i+3, len(sequence), 3):
                if sequence[j:j+3] in stop_codons:
                    orfs.append(sequence[i:j+3])
                    break
    return orfs





def find_longest_shortest(records):
    """Find the longest and shortest records in the given list of records."""
    longest_record = None
    shortest_record = None
    for record in records:
        if longest_record is None or len(record.seq) > len(longest_record.seq):
            longest_record = record
        if shortest_record is None or len(record.seq) < len(shortest_record.seq):
            shortest_record = record
    return longest_record, shortest_record

def main():
    records = list(SeqIO.parse("dna2.fasta", "fasta"))
    longest_record, shortest_record = find_longest_shortest(records)
    print(f"The longest record is {longest_record.id} with length {len(longest_record.seq)}")
    print(f"The shortest record is {shortest_record.id} with length {len(shortest_record.seq)}")
    for i, record in enumerate(records):
        print(f"Record number: {i+1}")
        print(f"ID: {record.id}")
        print(f"Sequence: {record.seq}")
        orf = find_orf(str(record.seq))
        len_orf = len(orf)
        if orf:
            print(f"ORF: {orf}", "the length of this orf is:", len_orf)
        else:
            print("No ORF found")
        print()
        #for i in range(3):
            #orfs = find_orf(record.seq[i:])
            #print(f"ORFs starting at position {i}:", orfs)


if __name__ == "__main__":
    main()
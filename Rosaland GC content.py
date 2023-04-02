with open('rosalind_gc.txt', 'r') as f:
    seq = ''
    for line in f:
        if line.startswith('>'):  # if line is a header
            if seq:  # if a sequence has been read
                gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
                print(f'GC content of {header}: {gc_content:.2f}%')
                seq = ''  # reset sequence
            header = line.strip()[1:]  # remove ">" and newline characters
        else:
            seq += line.strip()

    # process the last sequence
    gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
    print(f'GC content of {header}: {gc_content:.2f}%')


with open('rosalind_gc.txt', 'r') as f:
    for line in f:
        if line.startswith('>'):  # if line is a header
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
            print(f'GC content of {header}: {gc_content:.2f}%')
            seq = ''  # reset sequence
            header = line.strip()[1:]  # remove ">" and newline characters
        else:
            seq += line.strip()

    # process the last sequence
    gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
    print(f'GC content of {header}: {gc_content:.2f}%')

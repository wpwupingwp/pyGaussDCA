import numpy as np


# pythran export load_fasta(str, float)
# pythran export load_fasta(str)
def load_fasta(fasta_file: str, max_gap_fraction=0.9):
    """
    Import from treebarcode.
    Upper letter
    """
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
               'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
               'N': 12, 'O': 21, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
               'V': 18, 'W': 19, 'Y': 20,
               'U': 21, 'Z': 21, 'X': 21, 'J': 21}
    seqs = list()
    f = open(fasta_file, 'r')
    record = []
    for line_raw in f:
        line = line_raw.strip()
        if line.startswith('>'):
            if len(record) != 0:
                join_str = ''.join(record)
                if len(join_str.strip()) != 0:
                    seqs.append(join_str.upper())
                record.clear()
        else:
            record.append(line)
    if len(record) != 0:
        join_str = ''.join(record)
        if len(join_str.strip()) != 0:
            seqs.append(join_str)
    parsed = []
    for seq in seqs:
        if seq.count('-')/len(seq) >= max_gap_fraction:
            continue
        mapped_seq = [mapping.get(i, 22) for i in seq]
        parsed.append(mapped_seq)
    assert len(parsed)>0
    f.close()
    return np.array(parsed, dtype=np.int8).T


import sys
from Bio import SeqIO, Seq
import numpy as np

if __name__ == '__main__':
    fastaPath = sys.argv[1]
    col_path = sys.argv[2]
    out_path = sys.argv[3]
    fastaFile=list(SeqIO.parse(fastaPath, 'fasta'))
    cols = list(np.int64(open(col_path, 'r').read().split(',')))

    for i in range(len(fastaFile)): 
        aln = fastaFile[i]
        fastaFile[i].seq = Seq.Seq("".join(np.array(list(aln.seq))[cols]))

    with open(out_path, "w") as output_handle:
        SeqIO.write(fastaFile, output_handle, "fasta")
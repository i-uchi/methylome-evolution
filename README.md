# Extended base substitution model for methylome evolution

This repository contains scripts to generate methylome sequences with the extended nucleotide code including methylated bases and run phylogenetic analyses of them using RAxML-NG.

## Requirements

- RAxML-NG
- ClustalOmega

## Installation

- Install scripts

```
tar xvfz conv10b.tgz
```

## Input data Preparation

Create **orig_data** directory, which contains subdirectories for individual strains to be compared. Each subdirectory should contains the following files:

- **origseq.fna**: original genomic sequences in FASTA format.
- **origseq.gff**: gene annotations in GFF format.
- **motif.data**: methylation motifs generated from the modificaiton data obtained from SMRT sequencing.

See **orig_data** directory in the above sample data.


## Exection

1. Preparation of 8-base alignment. The script contains the following steps: 1) create methylome sequences, 2) extract 4-base and 8-base CDS sequences, 3) create 4-base alignment of the extracted CDS, and 4) convert the 4-base alignments into the alignment with the extended code.

```
./conv10b_all.sh
```

2. Inference of tree topology using 4-base GTR model 


```bash
./4branch.sh
```

3. Estimation of parameters using 8-base GTR model 


```sh
./8branch.sh
```

This evaluates branch lengths and substitution parameters under the [multi state](https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model:~:text=Morphological/multistate) GTR model using a fixed tree topology.

⸻

For detailed usage of RAxML-NG, refer to the official documentation:
👉 https://github.com/amkozlov/raxml-ng/wiki

# Extended base substitution model for methylome evolution

This repository contains scripts to generate methylome sequences with the extended nucleotide code including methylated bases and run phylogenetic analyses of them using RAxML-NG.

## Requirements

- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [ClustalOmega ](http://www.clustal.org/omega/)

## Installation

- Install scripts

```
tar xvfz conv10b.tgz
```

## Input data preparation

Create `orig_data` directory containing subdirectories for each strains to be compared, with the strain identifier as the directory name. Each subdirectory should contain the following files:

- `origseq.fna`: Original genomic sequences in FASTA format.
- `origseq.gff`: Gene annotations in GFF format. Consistent sequence names should be used in `origseq.fna` and `origseq.gff`
- `motif.data`: Methylation motifs generated from modificaiton data obtained by SMRT sequencing. The file should be tab-delimited, and the first four columns should contain: motif string, position of methylated base, motification type (m6A/m4C/m5C), and methyllation fraction.


The package contains a sample dataset, which is part of the H. pylori methylome data used in the paper by Yoshida et al. See `orig_data` directory bundled with this package.


## Usage

1. Preparation of the methylome alignment using the extended code.

```
./conv10b_all.sh
```

The overall procedure consists of the following steps: 1) create methylome sequences, 2) extract 4-base and extended-base CDS sequences, 3) create 4-base alignments of the extracted CDSs in the `orig_align` directory, and 4) create the alignments with the extended code  in the `ali10b` directory.



2. Inference of the tree topology using the 4-base GTR model 


```bash
./4branch.sh
```

3. Estimation of the parameters using the 8-base GTR model 


```sh
./8branch.sh
```

This evaluates branch lengths and substitution parameters under the [multi state](https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model:~:text=Morphological/multistate) GTR model using a fixed tree topology.


For detailed usage of RAxML-NG, refer to the official documentation:
👉 https://github.com/amkozlov/raxml-ng/wiki

## References

Yoshida, S., Uchiyama, I., Fukuyo, M., Kato, M., Rao, D. N., Konno, M., Fujiwara, S., Azuma, T. Kobayashi, I., Kishino, H.
Towards molecular evolutionary epigenomics with an expanded nucleotide code involving methylated bases, submitted.

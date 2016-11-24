# signasel
This program (R package) aims to detect selection/Ne by computing the likelihood of SNP allele frequencies.

## GENERAL NOTATIONS

* '1' and '2' refer to the 'first' (initial) and 'last' (final) populations, respectively, regardless of the number of generations between them.
* ng: total number of generations between population 1 and population 2. Hence, ng = 1 if pop 2 is the direct offspring of pop 1. More generally, ng = t2 - t1 if pop 1 lives at time t1 and pop 2 at time t2.
* n1, n2: sizes of populations 1 and 2, respectively. 
**Warning**: if ng > 1, the sizes of populations at intermediates generations are all n2, so that at generation 1: pop size n1->n2.
* S1, S2: sizes of the samples taken from populations 1 and 2, respectively. **NB**: Computations are made as if the total size of population 2 was S2 (not n2). This is strictly equivalent to first considering n2 individuals then drawing S2 individuals among those n2, and obviously saves time. *Important*: as a generalization S2 > n2 is allowed. (??)
* i1, i2: number of copies of the selected allele in samples S1 and S2, respectively.
* s: coefficient of selection, such that the fitnesses are 1, 1+s, 1+2s for genotypes aa, aA, AA (in the case where A is the positively selected allele).

**IMPORTANT**: population sizes (n1, n2) and sample sizes (S1, S2) are expressed in number of (diploid) individuals, so allele frequency reads p = i / 2n

## Data format

For each sample (row) the columns are: g, i, S, N:
g1 i1 S1 N1
g2 i2 S2 N2
g3 i3 S3 N3
etc.

with g: generation, i: number of allele copies, S: sample size, N: (effective) population size.


## Objectives

1. Information about intermediate generations.
2. Combine estimation of Ne AND selection rate (on =/= markers)
3. Flexibility of recursion: a/ exacte b/ diffusion c/ simul


Exercises for 16S ribosomal RNA data analysis
-----

A common task in bacterial ecology is PCR amplification and next-gen sequencing of 16S ribosomal RNA. This is a robust and straightforward way to profile the types of bacteria in a given environment. This set of exercises starts from the point when the sequencing has been completed and the user has been sent a download link for the sequences.

First the sequence reads need to be merged and filtered for high-quality sequences. Subsequently identical sequences are collapsed. Finally the sequences are clustered into OTUs - Operational Taxonomic Units - which is bacterial concept used instead of trying to adapt the concept of species to bacteria. We will follow the [UPARSE workflow](https://www.drive5.com/usearch/manual/ex_min2.html).

Start by creating a new directory under your home folder for the exercise. Copy all the gzipped files from /wrk/bio_workshop/ to this folder. Then continue using the UPARSE workflow:

```bash
/wrk/bio_workshop/usearch -threads 3 -fastq_mergepairs *_R1_*.fastq -relabel @ -fastqout merged.fq

/wrk/bio_workshop/usearch -fastq_filter merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout filtered.fa

/wrk/bio_workshop/usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

/wrk/bio_workshop/usearch -cluster_otus uniques.fa -otus otus.fa -relabel Otu
```

Exercises:
1. How many reads are in the raw files?
2. How many unique sequences are in the data?
3. How many times is (are) the most occurring sequence(s) observed?
4. How many OTUs are in the data?


The 16S rRNA sequences can be identified by comparing them to databases of known 16S rRNA sequences. This can be done for instance by using the assign_taxonomy.py script provided by [Qiime](http://qiime.org).

```bash
module load qiime

source activate qiime-1.9.1

assign_taxonomy.py -i otus.fa -r /wrk/bio_workshop/SILVA123_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta -t /wrk/bio_workshop/SILVA123_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_all_levels.txt -o silva_bac_taxonomy

source deactivate
```

Exercises:
1. Which are the three most common phyla, classes and genera in the data set?

Hint: use two successive awk commands with different delimiters, and pipe to | sort | uniq -c | sort -rn to count the most abundant ones.


OTU tables show bacteria OTUs on rows and sampling sites on columns, permitting comparison of multiple sampling sites for their bacterial content. usearch can be used to construct OTU tables:

```bash
usearch -otutab merged.fq -otus otus.fa -otutabout otutab.txt
```


Before a phylogenetic tree of the OTUs can be constructed, the OTUs need to be made comparable by aligning their 16S rRNA sequences.

```bash
/wrk/bio_workshop/sina-1.2.11/sina -i otus.fa --intype fasta -o otus_align.fa --outtype fasta --ptdb /wrk/bio_workshop/sina-1.2.11/SSURef_119_SILVA_14_07_14_opt.arb
```


Phylogenetic tree of the aligned OTU 16S rRNA sequences can be constructed using for instance FastTree.

```bash
FastTree -nt otus_align.fa > otus_align.tre
```

Exercise:
1. How does the tree file look like?
2. Visualise the resulting tree on [iToL](http://itol.embl.de).
3. How could you present abundance data around the phylogenetic tree?


## INTRODUCTION
This is the second part of the UTU Microbial Genomics course tutorial. Before attempting this part make sure you have completed [Part I](https://github.com/slhogle/UTU_microbial_genomics/blob/master/PART_I.md) first!

## DISCOVERING NUCLEOTIDE VARIANTS IN MEEI01
In the prior section we estimated the shared nucleotide identity between MEEI01 and EC958 and they are very closely related (over 99.5% similar) For the vast majority of their genomes they have identical nucleotide sequences. However, given the unique phenotype of MEEI01 we have reason to suspect that there is something significantly different between these two strains. Here we'll perform some analyses examine the presence of nucleotide variants in the MEEI01 genome and determine which genes these variants affect.

Going from quality controlled sequencing reads to information and insights about variants is generally comprised of three steps:

1. Alignment of short sequencing reads to a reference genome and sorting/arranging those reads
2. Calling sequencing variants
3. Filtering and annotating sequencing variants.

Take a look at the below script. Create a new file called `map2var.sbatch` and paste the below contents inside.

```bash
#!/usr/bin/env bash
#SBATCH -J BWA
#SBATCH -o BWA.log
#SBATCH --time 0-02:30:00
#SBATCH -p serial
#SBATCH --mem 8GB
#SBATCH -n 8

module load biokit

mkdir alignments variants

bwa index genomes/EC958.fna

bwa mem	-t 8 genomes/EC958.fna /wrk/bio_workshop/NZ_LNAE01000014_1.fq /wrk/bio_workshop/NZ_LNAE01000014_2.fq > alignments/MEEI01.sam

samtools sort --threads 8 -l 4 -O bam -o alignments/MEEI01.bam alignments/MEEI01.sam

bcftools mpileup -Ou -f genomes/EC958.fna alignments/MEEI01.bam | bcftools call --ploidy 1 -vcO z -o variants/MEEI01.vcf.gz

bcftools filter -i'QUAL>10' variants/MEEI01.vcf.gz | bcftools filter -i'QUAL>10' | bcftools filter -i'DP>50' | bcftools filter -i'IMF > 0.8' -o variants/MEEI01.filtered.vcf
```

Now we'll submit this script to the SLURM scheduler on Taito.

```bash
sbatch map2var.sbatch
```
This will take a while to complete because it's a relatively computationally intensive process. In the meantime let's take a look at the script we just made.

What's different in the #SBATCH header from the blast script we submitted? What do you think is going on? Read about the tools [bwa](https://github.com/lh3/bwa), [samtools](http://www.htslib.org/doc/samtools.html), and [bcftools](http://www.htslib.org/doc/bcftools.html) to get a general idea of what they are for. Use `zless` to examine the header of `variants/MEEI01.vcf.gz`. From the header descriptions can you figure out what `DP>50` and `IMF>0.8` mean in the final `bcftools filter` command? Hint: `zgrep` is your friend here. How many different insertions/deletions/polymorphisms passed our quality control filters?

Now we want to take a closer look these variants and what genes they occur within. I searched around for tools to annotate VCF files using GFF files, and I even tried a few, but I couldn't get any to work. So as with most problems in bioinformatics I simply found it easier to just write my own solution instead. This solution uses [python](https://www.python.org/) which is a very powerful scripting langauge. Create a new file called `tinyVCF_parser.py` with the text below.

```python
#!/usr/bin/env python
import sys
import csv

MYGFF=sys.argv[1]
MYVCF=sys.argv[2]
MYOUT=sys.argv[3]

with open(MYGFF, 'rb') as csvfile:
	GFFreader=csv.reader((row for row in csvfile if not row.startswith('#')), delimiter='\t')
	GFFdata = [r for r in GFFreader]

with open(MYVCF, 'rb') as csvfile:
	VCFreader=csv.reader((row for row in csvfile if not row.startswith('#')), delimiter='\t')
	VCFdata = [r for r in VCFreader]

with open(MYOUT, 'w') as f_out:
	for VCFrow in VCFdata:
		for GFFrow in GFFdata:
			if (GFFrow[2] == 'CDS') and (int(GFFrow[3]) <= int(VCFrow[1]) <= int(GFFrow[4])):
				f_out.write(VCFrow[1]+"\t"+VCFrow[3]+"\t"+VCFrow[4]+"\t"+VCFrow[7].split(";")[0]+"\t"+GFFrow[8]+"\n")
```

Now execute this script by typing:

```bash
python tinyVCF_parser.py genomes/EC958.gff variants/MEEI01.filtered.vcf variants/MEEI01.filtered.annotated.vcf
```

Inspect the resulting tab delimited file `variants/MEEI01.filtered.annotated.vcf`. All of these 'high confidence' mutations have been observed in other _E coli_ strains except for __one single mutation__. The unique mutation is a frame shift indel in a gene called _igaA._ Find the entry for this gene and it's Genbank identifier (hint: it starts with WP). Go to the NCBI [Conserved Domain database](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi), enter the ID into the search box, and submit. Read about the annotation for this gene. What is the gene called in _E coli_? What is the general function of this gene family?

Read this [journal article](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1003802) (and do your best to understand what's going on). After reading this article can you figure out why this _E coli_ strain caused this rather unusual occular infection?


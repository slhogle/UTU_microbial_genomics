## INTRODUCTION
Today we'll work through a real-world scenario that uses microbial genome sequencing and sequence variant discovery to identify the genotype in an E. Coli strain partially causing a highly virulent phenotype. You'll need to use your knowledge of Unix, Bash scripting, and biology/genomics to solve this medical mystery. As always please ask for help if you become stuck!

The data and story we'll be following today comes from [this paper](https://jamanetwork.com/journals/jamaophthalmology/fullarticle/2552682): Van Tyne D, Ciolino JB, Wang J, Durand ML, Gilmore MS. Novel Phagocytosis-Resistant Extended-Spectrum β-Lactamase–Producing _Escherichia coli_ From Keratitis. JAMA Ophthalmol. 2016;134(11):1306–1309. doi:10.1001/jamaophthalmol.2016.3283

## PRETEXT
A patient was admitted to a local hospital with recurrent inflammation and pain in the right eye. A culturing biopsy revealed a corneal ulcer and a bacterial infection caused by an extended-Spectrum β-Lactamase producing (broad antibiotic resistance) strain of _E coli_. This _E coli_ strain, named MEEI01, was highly resistant to nearly all available antibiotics, but was finally successfully treated with the last-line antibiotics amikacin and polymyxin B–trimethoprim. _E coli_ is a well-known causative agent of extraintestinal infection (bacterial infection that escapes the intestine and infects other internal organs), as well as more common urinary tract, intra-abdominal, and wound infections. However, it rarely causes eye infections, which prompted medical doctors to isolate the strain and sequence its genome in order to determine the underlying virulence genotype(s). Here we'll reproduce the steps taken by the scientists in their analysis of this _E coli_ strain MEEI01 to see if it: 

1. contains known genes encoding antibiotic resistance
2. contains other genetic variants (single nucleotide polymorphisms, insertions/deletions) relevant to its virulence

## ANTIBIOTIC RESISTANCE GENES
On your computer (not the Taito cluster) go to this ftp directory - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/442/595/GCF_001442595.1_ASM144259v1 and download the amino acid fasta file (most web browsers will download after clicking the link, but you may also have to right-click "Save Link as..."). NOTE: If you are using Microsoft Windows on the UTU public computers then you are not allowed to decompress gzipped files :anguished:. The protein sequences [are here](https://github.com/slhogle/UTU_microbial_genomics/blob/master/GCF_001442595.1_ASM144259v1_protein.faa) for the convenience of Windows users.

Although this is already a full assembled and annotated genome we'll pretend that we've just assembled it without any annotations or predictions :) The first step we'll take is to predict the presence of known antibiotic resistance genes using the [CARD server](https://academic.oup.com/nar/article/45/D1/D566/2333912). Because the MEEI01 strain displayed resistance to multiple types of antibiotics we can probably assume it is carrying at least a few different antibiotic resistance genes. 

Go to this [web page](https://card.mcmaster.ca/analyze/rgi) and submit the protein amino acid file you just downloaded to the server. Make sure you choose 'Perfect and strict hits only' so you don't get too many low quality hits.

1. How many potential resistance genes does MEEI01 carry?
2. How many different resistance mechanisms do these genes encompass? Which mechanisms?
3. Genes with a 'perfect' hit are predicted to provide resistance to how many drugs and which ones? Search some of these antibiotics on wikipedia and read about their usage. Why is the presence of these resistance genes in the MEEI01 genome particularly worrisome?

Read [this article](https://www.helsinki.fi/fi/uutiset/elamantieteet/antibioottien-maailmanloppu-vai-vallankumous) and [this one](https://www.theguardian.com/society/2017/oct/08/world-faces-antibiotic-apocalypse-says-chief-medical-officer). What is one likely cause for the rise of antibiotic resistance over the past two decades?

## MUTATIONS LEADING TO INCREASED VIRULENCE
Ideally, we would have performed the above analysis on the command line on the Taito cluster, but it would have required you all to download and install some software which is beyond the scope of this class. However, we can use some software from packages already available on Taito to search for genomic variants in the MEEI01 genome. What do I mean when I talk about 'genetic variants?' Above we were looking for the presence of whole genes that encoded a specific function. However, smaller and more subtle changes like the addition or deletion of just 10-20 nucleotides to a gene or promoter region can also cause resistance to antibiotics. It's more challenging to search for these subtle differences in a newly sequenced genome compared with the presence/absence of entire genes because you need to carefully account for additional variables such as sequencing coverage and read quality. First we're going to set up our directory structure for this part of the course and then we'll do some analyses. 

### ENVIRONMENT SETUP
First thing's first... In your home directory `~/` create a directory for hosting this project. Call it something like `UTU_microbial_genomics`. 

In the `UTU_microbial_genomics` directory create a `genomes` directory, `cd` into it, and type:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/285/655/GCF_000285655.3_EC958.v1/GCF_000285655.3_EC958.v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/285/655/GCF_000285655.3_EC958.v1/GCF_000285655.3_EC958.v1_genomic.gff.gz
```

Unzip these files and rename them to `EC958.fna, and EC958.gff`. These are assembled genome (fna) and general feature format (gff) files for [_E coli_ strain EC958](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0104400). This is a model strain from the ST131 subgroup of multi-drug resistant _E coli_ with a high quality closed assembled genome. We'll be using it in this tutorial as a reference to compare with our new strain MEEI01.

Now type:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/442/595/GCF_001442595.1_ASM144259v1/GCF_001442595.1_ASM144259v1_genomic.fna.gz
```

Like above, unzip the newly downloaded file and rename it to `MEEI01.fna`. 

### PERCENT NUCLEOTIDE IDENTITY BETWEEN MEEI01 and EC958
Let's take a look and compare the similarity between the genomes of our new strain (MEEI01) and EC958. One way to do this is compare percentage of shared nucleotides between two aligned genomes. This is done by aligning two genomes pairwise and looking for the shared fraction of aligned base pairs and dividing by the length of the aligned region. Here we'll use probably the most commonly used and well-known bioinformatic algorithm named BLAST for Basic Local Alignment Search Tool. We'll align these two genomes against each other and calculate the fraction of shared nucleotides (...more-or-less because the genomes we downloaded include more highly variable plasmid sequences which usually would be excluded from an analysis like this).

Type
```bash
module load biokit
makeblastdb -help
blastn -help
```

Read the help and see if you can figure out anything about what these programs do. An introduction to algorithms behind BLAST is beyond the scope of this tutorial. There is a detailed manual [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) which may be useful. For the remainder of this tutorial you'll just be typing commands in kind of "black box" manner, but you just need to know that BLAST compares biological sequences by finding regions of	local similarity. BLAST can compute under a null statistical model the "expectation" of seeing an alignment between two unrelated sequences simply by chance and shows you the "validity" of your alignment relative to this random model. This gives you a type of statistical significance for your BLAST alignments, which you can read more about [here](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html). In practice you should probably be skeptical of BLAST alignments with E-values greater than 1e-5. If you'd like to more details about what you are actually doing when you type these commands ask your instructors!

To submit real-world tasks to Taito we don't want to just execute them from the command line because we're currently logged into the login node. You can see which login node you are using by visually checking the `@taito-loginX` after your username at the command prompt. The login nodes are shared by lots of different users and are not intended to perform intensive computational tasks. Tasks like creating directories, deleting files, renaming things and using text editors is fine for the login node because these tasks don't consume many resources. However, once you start doing tasks that use multiple processors and tens of gigabytes or RAM then this would quickly slow down and disrupt the login node. DON'T SUBMIT BIG JOBS TO THE LOGIN NODES! Doing so might get you a grumpy email from the CSC admins :)

Create the following text file in `UTU_microbial_genomics` directory and name it `blastANI.sbatch`

```bash
#!/usr/bin/env bash
#SBATCH -J blastANI
#SBATCH -o blastANI.log
#SBATCH --time 0-00:30:00
#SBATCH -p serial
#SBATCH --mem 1GB
#SBATCH -n 1

module load biokit

mkdir perc_id

makeblastdb -in genomes/EC958.fna -dbtype nucl -out genomes/EC958
blastn -query genomes/MEEI01.fna -db genomes/EC958 -outfmt 6 -out perc_id/MEEI01vEC958.tab
```

Now at the command prompt type:
```bash
sbatch blastANI.sbatch
```

You can check the progress of jobs that you've submitted by typing:
```bash
squeue -u <MYUSERNAME>
```

While you're waiting for the job to finish read about the [SLURM Scheduler](https://slurm.schedmd.com/overview.html) and in particular [sbatch](https://slurm.schedmd.com/sbatch.html). What did the script do that you just created? Can you describe what each of the lines starting with #SBATCH mean?

By this time the blast search should have finished... Take a look at the resulting output file. You can read about what the different columns mean [here](https://www.drive5.com/usearch/manual/blast6out.html). Let's calculate the mean percent identity of all alignment contigs/chromosomes using the values in the third column from the blast report.

```bash
cat MEEI01vEC958.tab | grep -E "NZ_HG941718.1" | sort -k1,1 -k5,5g -k6,6g | sort -u -k1,1 --merge | awk '{ sum += $3 } END { if (NR > 0) print sum / NR }'
```

What is the shared nucleotide percentage between these two genomes? Can you describe what the above command is doing? Note the usage of pipes "|" and how that prevented us from writing unnecessary intermediate files. Take a look at the final `awk` command. Can you figure out what is going on? What is the NR variable doing?

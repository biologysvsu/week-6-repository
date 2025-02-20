# Blast Wrap-Up
# BLAST: Basic Local Alignment Search Tool.
The software tool finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.
---

# BLASTn
## **Organizing Your Directories**
To keep our work organized, we will create separate directories for different BLAST runs.

1. **Go to your personal storage directory**:
   ```bash
   cd /ocean/projects/agr250001p/your-username

2. **Create a directory for the CDS `blastn` run**:
   ```bash
   mkdir cdsblast
   
2. **Go to the newly created directory**:
   ```bash
   cd cdsblast
---
## Obtaining Input Data
We need to download reference data from the NCBI FTP server.

1. Download the CDS (Coding DNA Sequences) file:
``` bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz
```
2. Unzip the CDS file:
``` bash
gunzip GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz
```
---
## Make the CDS blast database

2. Create a BLAST database from the CDS file:
``` bash
module load BLAST
makeblastdb -in GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna -dbtype nucl -out cds_db
```
---
# Running `blastn` searches

## *blastn*
The blastn command compares a nucleotide query sequence to a nucleotide database.
We only need a `query` sequence to run blastn. We will use the same unknown sequence you use for homework 3. You do not have this file in your current directory, so make a copy of the file in here. There are two ways to do that, a long one and a short one:

### The long way

``` bash
cp /ocean/projects/agr250001p/your-username/unkx.fasta /ocean/projects/agr250001p/your-username/cdsblast
```
Make sure to replace unkx with your actual unknown numbe, and "your-username" with your actual one.

### The short way

``` bash
cp ../unkx.fasta .
```
Make sure to replace unkx with your actual unknown number.

### You can blast now using:

``` bash
blastn -query unk8.fasta -db cds_db -out blastn_results.txt -outfmt 6
```
- Make sure you replace `
---
# BLASTx

## Organize directory

1. For the `blastx` run, make the proteinblast folder and cd into the folder
   ```bash
   mkdir proteinblast
   cd proteinblast

2. Download Proteome (protein sequences) file from ncbi ftp
   ``` bash
   cd ../proteinblast
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz
   ```
3. Unzip the proteome file:
``` bash
gunzip GCF_000001405.40_GRCh38.p14_protein.faa.gz
```
## Make the protein database
1. Create a BLAST database from the protein file:
``` bash
makeblastdb -in GCF_000001405.40_GRCh38.p14_protein.faa -dbtype prot -out protein_db
```

## *blastx*
The blastx command translates a nucleotide sequence and searches it against a protein database.
``` bash
blastx -query my_sequence.fasta -db protein_db -out blastx_results.txt -outfmt 6
```

## blastp (not running today)
The blastp command compares a protein sequence to a protein database.
``` bash
blastp -query my_protein.fasta -db protein_db -out blastp_results.txt -outfmt 6

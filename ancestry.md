We have a VCF file, it is called:
mystery.vcf
Login to the HPC
Go to your ocean folder

### Add plink to environment
```
echo 'export PATH=/ocean/projects/bio260081p/shared/software/plink2:$PATH' >> ~/.bashrc
```
and
```
source ~/.bashrc
```

NOTE: Human genotype data comes from PLINK resources: https://www.cog-genomics.org/plink/2.0/resources, pgen, pvar and psam files were downloaded using the wget command and the link of files
-Files were renamed using the mv command to all_hg38.pgen.zst  all_hg38.pvar.zst and all_hg38.psam 
- file were uncompressed using plink, interact and then:
```
interact -t 3:00:00 --ntasks-per-node=4 --mem=32G
plink2 --zst-decompress all_hg38.pgen.zst all_hg38.pgen
#verify genomes
plink2 --pfile all_hg38 vzs --write-samples --out reference_test

# Identify the sample
grep '^#CHROM' mistery.vcf

# Check the reference genome
grep '^##reference' mistery.vcf

# Examine chromosome naming
grep '^##contig' mistery.vcf | head

# Count variant records
grep -vc '^#' mistery.vcf

#checking chromosome reference names
plink2 --pfile all_hg38 vzs --write-snplist --out reference_snps
head all_hg38.psam

#And inspect the reference variant file:
zstdcat all_hg38.pvar.zst | head -n 15

#check mistery vcf has variants
grep -v '^#' mistery.vcf | head -n 5

#convert autosome names
awk 'BEGIN {FS=OFS="\t"}
/^#/ {print; next}
$1 ~ /^NC_0000[0-9][0-9]\./ {
    chr=substr($1,8,2)+0
    if (chr>=1 && chr<=22) {
        $1=chr
        print
    }
}' mistery.vcf > mystery_autosomes.vcf

#check result
grep -v '^#' mystery_autosomes.vcf | head -n 3
grep -vc '^#' mystery_autosomes.vcf
```
### TUTORIAL STARTS HERE
# STEP 1: Convert chromosome names
awk 'BEGIN {FS=OFS="\t"}
/^#/ {print; next}
$1 ~ /^NC_0000[0-9][0-9]\./ {
    chr=substr($1,8,2)+0
    if (chr>=1 && chr<=22) {
        $1=chr
        print
    }
}' mistery.vcf > mystery_autosomes.vcf

# STEP 2: Convert VCF to PLINK
plink2 \
  --vcf mystery_autosomes.vcf \
  --double-id \
  --snps-only just-acgt \
  --set-all-var-ids '@:#:$r:$a' \
  --make-pgen \
  --out mystery

# STEP 3: Prepare reference
plink2 \
  --pfile all_hg38 vzs \
  --autosome \
  --snps-only just-acgt \
  --set-all-var-ids '@:#:$r:$a' \
  --rm-dup force-first \
  --make-pgen \
  --out reference_clean

# STEP 4: LD pruning
plink2 \
  --pfile reference_clean \
  --indep-pairwise 200 50 0.2 \
  --out reference_pruned

# STEP 5: Identify shared SNPs
plink2 \
  --pfile mystery \
  --extract reference_pruned.prune.in \
  --write-snplist \
  --out mystery_shared

# STEP 6: Reference PCA and allele weights
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --pca allele-wts 10 \
  --out reference_pca

# STEP 7: Reference allele frequencies
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --nonfounders \
  --freq counts \
  --out reference_freq

# STEP 8: Mystery projection
plink2 \
  --pfile mystery \
  --extract mystery_shared.snplist \
  --read-freq reference_freq.acount \
  --score reference_pca.eigenvec.allele 2 5 header-read \
    no-mean-imputation variance-standardize \
  --score-col-nums 6-15 \
  --out mystery_projected

# STEP 9: Reference projection
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --read-freq reference_freq.acount \
  --score reference_pca.eigenvec.allele 2 5 header-read \
    no-mean-imputation variance-standardize \
  --score-col-nums 6-15 \
  --out reference_projected

# STEP 10: Plot PCA
cat > plot_ancestry.py <<'PY'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ref = pd.read_csv("reference_projected.sscore", sep=r"\s+")
mystery = pd.read_csv("mystery_projected.sscore", sep=r"\s+")
samples = pd.read_csv("reference_clean.psam", sep=r"\s+")

ref = ref.merge(samples[["#IID", "SuperPop"]], on="#IID")

fig, ax = plt.subplots(figsize=(10, 7))

for pop, group in ref.groupby("SuperPop"):
    ax.scatter(group["PC1_AVG"], group["PC2_AVG"],
               s=12, alpha=0.5, label=pop)

ax.scatter(mystery["PC1_AVG"], mystery["PC2_AVG"],
           s=180, marker="*", color="red",
           edgecolor="black", label="Mystery individual",
           zorder=10)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("Human Genetic Ancestry - Exploratory PCA")
ax.legend()
plt.tight_layout()
plt.savefig("ancestry_pca.png", dpi=300)
PY

python plot_ancestry.py




copy the files we will work on todY:
```
cp -r ../shared/ancestry/ .
```
download file
go to: https://ondemand.bridges2.psc.edu/
login with bridges credentials




plink2 \
  --pfile reference_clean \
  --extract reference_pruned.prune.in \
  --maf 0.05 \
  --pca 10 \
  --out reference_pca_clean

cat > plot_reference.py <<'PY'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pca = pd.read_csv("reference_pca_clean.eigenvec", sep=r"\s+")
samples = pd.read_csv("reference_clean.psam", sep=r"\s+")

df = pca.merge(samples[["#IID", "SuperPop"]], on="#IID")

for pop, group in df.groupby("SuperPop"):
    plt.scatter(group["PC1"], group["PC2"],
                s=10, alpha=0.5, label=pop)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("1000 Genomes — Human Population Structure")
plt.legend()
plt.tight_layout()
plt.savefig("reference_pca_clean.png", dpi=300)
PY

python plot_reference.py

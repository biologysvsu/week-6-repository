# BIOL 443 — HUMAN GENETIC ANCESTRY ANALYSIS
# PLINK 2 | 1000 Genomes + Mystery Individual

# STEP 1 — Convert chromosome names and keep autosomes
```
awk 'BEGIN {FS=OFS="\t"}
/^#/ {print; next}
$1 ~ /^NC_0000[0-9][0-9]\./ {
    chr=substr($1,8,2)+0
    if (chr>=1 && chr<=22) {
        $1=chr
        print
    }
}' mistery.vcf > mystery_autosomes.vcf
```

# STEP 2 — Convert mystery VCF to PLINK format
```
plink2 \
  --vcf mystery_autosomes.vcf \
  --double-id \
  --snps-only just-acgt \
  --set-all-var-ids '@:#:$r:$a' \
  --make-pgen \
  --out mystery
```

# STEP 3 — Prepare 1000 Genomes reference dataset
```
plink2 \
  --pfile all_hg38 vzs \
  --autosome \
  --snps-only just-acgt \
  --set-all-var-ids '@:#:$r:$a' \
  --rm-dup force-first \
  --make-pgen \
  --out reference_clean
```

# STEP 4 — LD pruning
```
plink2 \
  --pfile reference_clean \
  --indep-pairwise 200 50 0.2 \
  --out reference_pruned
```

# PART A — PCA USING ONLY 1000 GENOMES

# STEP 5 — Calculate reference PCA
```
plink2 \
  --pfile reference_clean \
  --extract reference_pruned.prune.in \
  --maf 0.05 \
  --pca 10 \
  --out reference_pca_clean
```

# STEP 6 — Plot reference PCA
```
cat > plot_reference.py <<'PY'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pca = pd.read_csv("reference_pca_clean.eigenvec", sep=r"\s+")
samples = pd.read_csv("reference_clean.psam", sep=r"\s+")

df = pca.merge(samples[["#IID", "SuperPop"]], on="#IID")

plt.figure(figsize=(10, 7))

for pop, group in df.groupby("SuperPop"):
    plt.scatter(
        group["PC1"],
        group["PC2"],
        s=10,
        alpha=0.5,
        label=pop
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("1000 Genomes - Human Population Structure")
plt.legend()
plt.tight_layout()
plt.savefig("reference_pca_clean.png", dpi=300)
PY
```
```
python plot_reference.py
```

# PART B — EXPLORATORY PCA WITH MYSTERY INDIVIDUAL

# STEP 7 — Identify shared SNPs
```
plink2 \
  --pfile mystery \
  --extract reference_pruned.prune.in \
  --write-snplist \
  --out mystery_shared
```

# STEP 8 — Calculate reference PCA and allele weights
```
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --pca allele-wts 10 \
  --out reference_pca
```

# STEP 9 — Calculate reference allele frequencies
```
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --nonfounders \
  --freq counts \
  --out reference_freq
```

# STEP 10 — Project mystery individual
```
plink2 \
  --pfile mystery \
  --extract mystery_shared.snplist \
  --read-freq reference_freq.acount \
  --score reference_pca.eigenvec.allele 2 5 header-read \
    no-mean-imputation variance-standardize \
  --score-col-nums 6-15 \
  --out mystery_projected
```

# STEP 11 — Project reference individuals
```
plink2 \
  --pfile reference_clean \
  --extract mystery_shared.snplist \
  --read-freq reference_freq.acount \
  --score reference_pca.eigenvec.allele 2 5 header-read \
    no-mean-imputation variance-standardize \
  --score-col-nums 6-15 \
  --out reference_projected
```

# STEP 12 — Plot exploratory PCA
```
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
    ax.scatter(
        group["PC1_AVG"],
        group["PC2_AVG"],
        s=12,
        alpha=0.5,
        label=pop
    )

ax.scatter(
    mystery["PC1_AVG"],
    mystery["PC2_AVG"],
    s=180,
    marker="*",
    color="red",
    edgecolor="black",
    label="Mystery individual",
    zorder=10
)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("Human Genetic Ancestry - Exploratory PCA")
ax.legend()
plt.tight_layout()
plt.savefig("ancestry_pca.png", dpi=300)
PY
```
```
python plot_ancestry.py
```

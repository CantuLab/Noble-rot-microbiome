# Noble Rot Microbiome Analysis

Code repository for:

**Botrytis cinerea infection reshapes the grape berry microbiome during noble rot**

Cantu D., Blanco-Ulate B., Allen G., Figueroa-Balderas R.

---

## Requirements

### Software
- FASTX-Toolkit v0.0.14
- QIIME 2 v2026.4 (conda environment)
- R v4.5

### R packages
```r
install.packages(c("tidyverse", "patchwork"))
```

### Python packages (QIIME 2 environment)
- scikit-bio
- scipy
- statsmodels
- networkx
- matplotlib

### QIIME 2 plugins
- cutadapt
- dada2
- feature-table
- feature-classifier
- taxa
- diversity
- phylogeny
- composition (ANCOM-BC)

---

## Directory structure

```
project/
├── raw/          # raw FASTQ files
├── data/         # processed data and QIIME2 artifacts
├── metadata/     # sample metadata mapping files
├── scripts/      # this repository
└── results/      # output figures and tables
```

All commands are run from `project/`.

---

## 1. Raw read preprocessing

```bash
RAW=path/to/raw/fastq
DATA=path/to/output

# Sequencing runs: run1 = 2015, run2 = 2012 to 2014
# (16S_run1 = 16S, 16S_run2 = 16s-2, ITS_run1 = BITS3, ITS_run2 = BITS4)
# Extract barcodes and trim reads for each sequencing run (FASTX-Toolkit)
for RUN in 16S_run1 16S_run2 ITS_run1 ITS_run2; do
    mkdir -p ${DATA}/${RUN}_import

    fastx_trimmer -f 1 -l 8 -Q33 \
        -i ${RAW}/${RUN}_R1.fastq \
        -o ${DATA}/${RUN}_import/barcodes.fastq

    fastx_trimmer -f 9 -Q33 \
        -i ${RAW}/${RUN}_R1.fastq \
        -o ${DATA}/${RUN}_import/forward.fastq

    fastx_trimmer -f 9 -Q33 \
        -i ${RAW}/${RUN}_R2.fastq \
        -o ${DATA}/${RUN}_import/reverse.fastq

    gzip ${DATA}/${RUN}_import/*.fastq
done
```

---

## 2. QIIME 2 pipeline

### 2.1 Import and demultiplex

```bash
conda activate qiime2

for RUN in 16S_run1 16S_run2 ITS_run1 ITS_run2; do
    case ${RUN} in
        16S_*) MAP=metadata/mapping_16s.txt ;;
        ITS_*) MAP=metadata/mapping_bits.txt ;;
    esac

    # Import reads
    qiime tools import \
        --type EMPPairedEndSequences \
        --input-path ${DATA}/${RUN}_import/ \
        --output-path ${DATA}/${RUN}.qza

    # Demultiplex (8-nt barcodes, no Golay error correction)
    qiime demux emp-paired \
        --i-seqs ${DATA}/${RUN}.qza \
        --m-barcodes-file ${MAP} \
        --m-barcodes-column BarcodeSequence \
        --p-no-golay-error-correction \
        --o-per-sample-sequences ${DATA}/${RUN}_demux.qza \
        --o-error-correction-details ${DATA}/${RUN}_demux_details.qza
done
```

### 2.2 Denoise with DADA2

```bash
# ITS: remove primers and 3' read-through (reverse-complemented primer, Illumina adapter)
for RUN in ITS_run1 ITS_run2; do
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${DATA}/${RUN}_demux.qza \
        --p-front-f ACCTGCGGARGGATCA \
        --p-front-r TTGYTRAAAGTT \
        --p-adapter-f AACTTTYARCAA --p-adapter-f AGATCGGAAGAGC \
        --p-adapter-r TGATCCYTCCGCAGGT --p-adapter-r AGATCGGAAGAGC \
        --p-times 3 \
        --p-match-adapter-wildcards \
        --p-discard-untrimmed \
        --o-trimmed-sequences ${DATA}/${RUN}_trimmed.qza
done

# 16S (no truncation)
for RUN in 16S_run1 16S_run2; do
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${DATA}/${RUN}_demux.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --o-table ${DATA}/${RUN}_table.qza \
        --o-representative-sequences ${DATA}/${RUN}_repseqs.qza \
        --o-denoising-stats ${DATA}/${RUN}_stats.qza \
        --o-base-transition-stats ${DATA}/${RUN}_base_transition_stats.qza
done

# ITS (no truncation)
for RUN in ITS_run1 ITS_run2; do
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${DATA}/${RUN}_trimmed.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --o-table ${DATA}/${RUN}_table.qza \
        --o-representative-sequences ${DATA}/${RUN}_repseqs.qza \
        --o-denoising-stats ${DATA}/${RUN}_stats.qza \
        --o-base-transition-stats ${DATA}/${RUN}_base_transition_stats.qza
done
```

### 2.3 Merge runs

```bash
# Each run was demultiplexed with the mapping file of all samples.
# Keep only the vintages sequenced on each run: run1 = 2015, run2 = 2012 to 2014.
for MARKER in 16S ITS; do
    case ${MARKER} in
        16S) MAP=metadata/mapping_16s.txt ;;
        ITS) MAP=metadata/mapping_bits.txt ;;
    esac

    qiime feature-table filter-samples \
        --i-table ${DATA}/${MARKER}_run1_table.qza \
        --m-metadata-file ${MAP} \
        --p-where "[Year]='2015'" \
        --o-filtered-table ${DATA}/${MARKER}_run1_2015.qza

    qiime feature-table filter-samples \
        --i-table ${DATA}/${MARKER}_run2_table.qza \
        --m-metadata-file ${MAP} \
        --p-where "[Year] IN ('2012', '2013', '2014')" \
        --o-filtered-table ${DATA}/${MARKER}_run2_2012_2014.qza

    qiime feature-table merge \
        --i-tables ${DATA}/${MARKER}_run1_2015.qza \
        --i-tables ${DATA}/${MARKER}_run2_2012_2014.qza \
        --p-overlap-method error_on_overlapping_sample \
        --o-merged-table ${DATA}/${MARKER}_merged_table.qza

    qiime feature-table merge-seqs \
        --i-data ${DATA}/${MARKER}_run1_repseqs.qza \
        --i-data ${DATA}/${MARKER}_run2_repseqs.qza \
        --o-merged-data ${DATA}/${MARKER}_merged_repseqs.qza
done
```

### 2.4 Taxonomic classification

```bash
# SILVA 138 99 % reference reads (515F/806R) and taxonomy (QIIME 2 data resources)
curl -sL -o silva-138-99-seqs-515-806.qza https://data.qiime2.org/2020.6/common/silva-138-99-seqs-515-806.qza
curl -sL -o silva-138-99-tax.qza https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza

# Train SILVA 138 classifier for 16S V4 region (515F/806R)
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138-99-seqs-515-806.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --o-classifier silva-138-99-515-806-classifier.qza

# Classify 16S ASVs
qiime feature-classifier classify-sklearn \
    --i-classifier silva-138-99-515-806-classifier.qza \
    --i-reads ${DATA}/16S_merged_repseqs.qza \
    --o-classification ${DATA}/16S_taxonomy.qza

# Classify ITS ASVs using UNITE v10.0
qiime feature-classifier classify-sklearn \
    --i-classifier unite-v10-99-classifier.qza \
    --i-reads ${DATA}/ITS_merged_repseqs.qza \
    --o-classification ${DATA}/ITS_taxonomy.qza

# Reference ITS sequences: B. cinerea CBS 261.71, B. pseudocinerea, S. sclerotiorum CBS 499.50
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MH860108,JN692379,MH856725&rettype=fasta" \
    > ${DATA}/its_reference_sequences.fasta
```

```python
# Assign Sclerotiniaceae ASVs to Botrytis: >= 97 % identity to Botrytis and closer to
# Botrytis than to Sclerotinia (semi-global alignment, ASV fully aligned)
import pandas as pd
import qiime2

DATA = "path/to/output"


def read_fasta(path):
    seqs, key = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            key = line[1:].split()[0].split(".")[0]
            seqs[key] = ""
        elif key:
            seqs[key] += line.upper()
    return seqs


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def identity(q, r, M=2, X=-3, G=-5):
    n, m = len(q), len(r)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0], T[i][0] = S[i - 1][0] + G, 1
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a = S[i - 1][j - 1] + (M if q[i - 1] == r[j - 1] else X)
            b, c = S[i - 1][j] + G, S[i][j - 1] + G
            S[i][j], T[i][j] = max((a, 0), (b, 1), (c, 2), key=lambda x: x[0])
    j = max(range(m + 1), key=lambda k: S[n][k])
    i, mat, length = n, 0, 0
    while i > 0:
        t = T[i][j]
        length += 1
        if j == 0 or t == 1:
            i -= 1
        elif t == 0:
            mat += q[i - 1] == r[j - 1]
            i, j = i - 1, j - 1
        else:
            j -= 1
    return 100 * mat / length


refs = read_fasta(f"{DATA}/its_reference_sequences.fasta")
tax = qiime2.Artifact.load(f"{DATA}/ITS_taxonomy.qza").view(pd.DataFrame)
seqs = qiime2.Artifact.load(f"{DATA}/ITS_merged_repseqs.qza").view(pd.Series)
scl = tax.index[tax["Taxon"].str.contains("f__Sclerotiniaceae")]

top = str(seqs[scl[0]])
ref = {k: (r if identity(top, r) >= identity(top, rc(r)) else rc(r)) for k, r in refs.items()}
botrytis = []
for asv in scl:
    q = str(seqs[asv])
    b = max(identity(q, ref["MH860108"]), identity(q, ref["JN692379"]))
    s = identity(q, ref["MH856725"])
    if b >= 97 and b > s:
        botrytis.append(asv)

tax.loc[botrytis, "Taxon"] = tax.loc[botrytis, "Taxon"].str.replace(
    r"f__Sclerotiniaceae.*$", "f__Sclerotiniaceae;g__Botrytis", regex=True)
qiime2.Artifact.import_data("FeatureData[Taxonomy]", tax).save(f"{DATA}/ITS_taxonomy_botrytis.qza")
```

### 2.5 Filter contaminants and low-depth samples

```bash
# Remove chloroplasts and mitochondria from 16S
qiime taxa filter-table \
    --i-table ${DATA}/16S_merged_table.qza \
    --i-taxonomy ${DATA}/16S_taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table ${DATA}/16S_merged_table_filtered.qza

# Remove samples with fewer than 1,000 reads
qiime feature-table filter-samples \
    --i-table ${DATA}/16S_merged_table_filtered.qza \
    --p-min-frequency 1000 \
    --o-filtered-table ${DATA}/16S_table_final.qza

qiime feature-table filter-samples \
    --i-table ${DATA}/ITS_merged_table.qza \
    --p-min-frequency 1000 \
    --o-filtered-table ${DATA}/ITS_table_final.qza

# Bacteria only, for genus-level analyses
qiime taxa filter-table \
    --i-table ${DATA}/16S_table_final.qza \
    --i-taxonomy ${DATA}/16S_taxonomy.qza \
    --p-include d__Bacteria \
    --o-filtered-table ${DATA}/16S_table_bacteria.qza

# Fungi other than Botrytis
qiime taxa filter-table \
    --i-table ${DATA}/ITS_table_final.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy_botrytis.qza \
    --p-exclude f__Sclerotiniaceae \
    --o-filtered-table ${DATA}/ITS_table_noBotrytis.qza
```

### 2.6 Phylogenetic tree (bacteria only)

```bash
# Single thread: MAFFT and FastTree with several threads give a different tree
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences ${DATA}/16S_merged_repseqs.qza \
    --p-n-threads 1 \
    --o-alignment ${DATA}/16S_aligned_repseqs.qza \
    --o-masked-alignment ${DATA}/16S_masked_aligned_repseqs.qza \
    --o-tree ${DATA}/16S_unrooted_tree.qza \
    --o-rooted-tree ${DATA}/16S_rooted_tree.qza
```

### 2.7 Collapse to genus and family level

```bash
# 16S genus (level 6), bacteria only
qiime taxa collapse \
    --i-table ${DATA}/16S_table_bacteria.qza \
    --i-taxonomy ${DATA}/16S_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${DATA}/16S_genus_table.qza

# ITS genus (level 6)
qiime taxa collapse \
    --i-table ${DATA}/ITS_table_final.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy_botrytis.qza \
    --p-level 6 \
    --o-collapsed-table ${DATA}/ITS_genus_table.qza

# ITS genus (level 6), Botrytis removed
qiime taxa collapse \
    --i-table ${DATA}/ITS_table_noBotrytis.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy_botrytis.qza \
    --p-level 6 \
    --o-collapsed-table ${DATA}/ITS_genus_table_noBotrytis.qza

# ITS family (level 5)
qiime taxa collapse \
    --i-table ${DATA}/ITS_table_final.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy_botrytis.qza \
    --p-level 5 \
    --o-collapsed-table ${DATA}/ITS_family_table.qza
```

---

## 3. Diversity analyses

### 3.1 Alpha diversity

```bash
for YEAR in 2012 2013 2014 2015; do
    # Bacteria - filter by year
    qiime feature-table filter-samples \
        --i-table ${DATA}/16S_table_final.qza \
        --m-metadata-file metadata/mapping_16s.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/16S_table_${YEAR}.qza

    qiime diversity alpha-phylogenetic \
        --i-table ${DATA}/16S_table_${YEAR}.qza \
        --i-phylogeny ${DATA}/16S_rooted_tree.qza \
        --p-metric faith_pd \
        --o-alpha-diversity ${DATA}/16S_faith_${YEAR}.qza

    for METRIC in shannon observed_features pielou_e chao1; do
        qiime diversity alpha \
            --i-table ${DATA}/16S_table_${YEAR}.qza \
            --p-metric ${METRIC} \
            --o-alpha-diversity ${DATA}/16S_${METRIC}_${YEAR}.qza
    done

    # Fungi - filter by year
    qiime feature-table filter-samples \
        --i-table ${DATA}/ITS_table_final.qza \
        --m-metadata-file metadata/mapping_bits.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/ITS_table_${YEAR}.qza

    for METRIC in shannon observed_features pielou_e chao1; do
        qiime diversity alpha \
            --i-table ${DATA}/ITS_table_${YEAR}.qza \
            --p-metric ${METRIC} \
            --o-alpha-diversity ${DATA}/ITS_${METRIC}_${YEAR}.qza
    done
done
```

### 3.2 Beta diversity

```bash
for YEAR in 2012 2013 2014 2015; do
    # Bacteria - UniFrac
    for METRIC in weighted_unifrac unweighted_unifrac; do
        qiime diversity beta-phylogenetic \
            --i-table ${DATA}/16S_table_${YEAR}.qza \
            --i-phylogeny ${DATA}/16S_rooted_tree.qza \
            --p-metric ${METRIC} \
            --o-distance-matrix ${DATA}/16S_${METRIC}_${YEAR}.qza

        qiime diversity pcoa \
            --i-distance-matrix ${DATA}/16S_${METRIC}_${YEAR}.qza \
            --o-pcoa ${DATA}/16S_pcoa_${METRIC}_${YEAR}.qza

        qiime diversity beta-group-significance \
            --i-distance-matrix ${DATA}/16S_${METRIC}_${YEAR}.qza \
            --m-metadata-file metadata/mapping_16s.txt \
            --m-metadata-column Stage \
            --p-method permanova \
            --p-permutations 999 \
            --o-visualization ${DATA}/16S_permanova_${METRIC}_${YEAR}.qzv
    done

    # Fungi - Bray-Curtis on genus-level relative abundance, all fungi and Botrytis removed
    for SET in genus_table genus_table_noBotrytis; do
        qiime feature-table filter-samples \
            --i-table ${DATA}/ITS_${SET}.qza \
            --m-metadata-file metadata/mapping_bits.txt \
            --p-where "[Year]='${YEAR}'" \
            --p-min-frequency 1 \
            --o-filtered-table ${DATA}/ITS_${SET}_${YEAR}.qza

        qiime feature-table relative-frequency \
            --i-table ${DATA}/ITS_${SET}_${YEAR}.qza \
            --o-relative-frequency-table ${DATA}/ITS_${SET}_${YEAR}_rel.qza

        qiime diversity beta \
            --i-table ${DATA}/ITS_${SET}_${YEAR}_rel.qza \
            --p-metric braycurtis \
            --o-distance-matrix ${DATA}/ITS_braycurtis_${SET}_${YEAR}.qza

        qiime diversity pcoa \
            --i-distance-matrix ${DATA}/ITS_braycurtis_${SET}_${YEAR}.qza \
            --o-pcoa ${DATA}/ITS_pcoa_${SET}_${YEAR}.qza

        qiime diversity beta-group-significance \
            --i-distance-matrix ${DATA}/ITS_braycurtis_${SET}_${YEAR}.qza \
            --m-metadata-file metadata/mapping_bits.txt \
            --m-metadata-column Stage \
            --p-method permanova \
            --p-permutations 999 \
            --o-visualization ${DATA}/ITS_permanova_${SET}_${YEAR}.qzv
    done
done
```

## 3.3 Alpha diversity group significance

```python
# Kruskal-Wallis test across noble rot stages, per marker, metric and year; where
# Kruskal-Wallis p < 0.05, pairwise two-sided Dunn's tests (tie-corrected, no correction
# for multiple comparisons) and compact letters (Fig. 2).
import itertools

import numpy as np
import pandas as pd
import qiime2
from scipy.stats import kruskal, norm, rankdata

DATA = "path/to/output"
STAGES = ["0_control", "1_S1", "2_S2", "3_S3"]
METRICS = {"16S": ["shannon", "observed_features", "chao1", "faith", "pielou_e"],
           "ITS": ["shannon", "observed_features", "chao1", "pielou_e"]}
MAPS = {"16S": "metadata/mapping_16s.txt", "ITS": "metadata/mapping_bits.txt"}


def dunn(d, groups):
    r = pd.Series(rankdata(d.value), index=d.index)
    n = len(d)
    _, t = np.unique(d.value, return_counts=True)
    c = 1 - (t ** 3 - t).sum() / (n ** 3 - n)
    out = {}
    for a, b in itertools.combinations(groups, 2):
        ia, ib = d.stage == a, d.stage == b
        z = (r[ia].mean() - r[ib].mean()) / np.sqrt(c * n * (n + 1) / 12 * (1 / ia.sum() + 1 / ib.sum()))
        out[(a, b)] = 2 * norm.sf(abs(z))
    return out


def letters(groups, sig):
    sets = [set(groups)]
    for a, b in sig:
        new = []
        for s in sets:
            new += [s - {a}, s - {b}] if (a in s and b in s) else [s]
        sets = [s for s in new if not any(s < o for o in new)]
        sets = [s for i, s in enumerate(sets) if s not in sets[:i]]
    sets.sort(key=lambda s: min(groups.index(g) for g in s))
    lab = {g: "" for g in groups}
    for i, s in enumerate(sets):
        for g in s:
            lab[g] += "abcdefgh"[i]
    return lab


rows, pairs_out = [], []
for marker, metrics in METRICS.items():
    meta = qiime2.Metadata.load(MAPS[marker]).to_dataframe()
    for year in [2012, 2013, 2014, 2015]:
        for metric in metrics:
            a = qiime2.Artifact.load(f"{DATA}/{marker}_{metric}_{year}.qza").view(pd.Series)
            d = pd.DataFrame({"value": a, "stage": meta.loc[a.index, "Stage"]})
            groups = [g for g in STAGES if (d.stage == g).any()]
            kw = kruskal(*[d.value[d.stage == g] for g in groups]).pvalue if len(groups) > 1 else np.nan
            ps = dunn(d, groups) if len(groups) > 1 else {}
            for (x, y), p in ps.items():
                pairs_out.append({"marker": marker, "year": year, "metric": metric,
                                  "stage_a": x, "stage_b": y, "p_dunn": p})
            lab = letters(groups, [pr for pr, p in ps.items() if p < 0.05]) if kw < 0.05 else {}
            row = {"marker": marker, "year": year, "metric": metric, "kruskal_p": kw}
            for g in groups:
                v = d.value[d.stage == g]
                row[g] = f"{v.mean():.3f} ± {v.std():.3f} (n = {len(v)})" if len(v) > 1 else f"{v.mean():.3f} (n = 1)"
                row[f"letter_{g}"] = lab.get(g, "")
            rows.append(row)

pd.DataFrame(rows).to_csv(f"{DATA}/alpha_table.csv", index=False)
pd.DataFrame(pairs_out).to_csv(f"{DATA}/alpha_pairwise.csv", index=False)
```

## 3.4 Alpha diversity summary table

Produces the per-year, per-metric tables of `mean ± SD` by stage with the Kruskal-Wallis
p-value (3.3) and the sequencing-depth checks reported in Supplementary Tables S2 and S3:
Spearman correlation between depth and each metric, and stage comparisons after rarefying
each vintage to its lowest sample depth (100 iterations).

```python
import numpy as np
import pandas as pd
import qiime2
from scipy.stats import kruskal, spearmanr
from skbio.diversity import alpha_diversity

DATA = "path/to/output"
ITER = 100


def rarefy(counts, depth, rng):
    pool = np.repeat(np.arange(len(counts)), counts)
    keep = rng.choice(pool, size=depth, replace=False)
    return np.bincount(keep, minlength=len(counts))


corr, rare = [], []
for marker in ["16S", "ITS"]:
    meta = qiime2.Metadata.load({"16S": "metadata/mapping_16s.txt", "ITS": "metadata/mapping_bits.txt"}[marker]).to_dataframe()
    for year in [2012, 2013, 2014, 2015]:
        t = qiime2.Artifact.load(f"{DATA}/{marker}_table_{year}.qza").view(pd.DataFrame)
        depth = t.sum(axis=1)
        stage = meta.loc[t.index, "Stage"]
        for metric in ["shannon", "observed_features", "chao1", "pielou_e"] + (["faith"] if marker == "16S" else []):
            a = qiime2.Artifact.load(f"{DATA}/{marker}_{metric}_{year}.qza").view(pd.Series)
            r, p = spearmanr(depth[a.index], a)
            corr.append({"marker": marker, "year": year, "metric": metric, "spearman_r": r, "p": p})

        d = int(depth.min())
        X = t.values.astype(int)
        reps = []
        for seed in range(ITER):
            rng = np.random.default_rng(seed)
            R = np.vstack([rarefy(row, d, rng) for row in X])
            reps.append(pd.DataFrame({
                "shannon": alpha_diversity("shannon", R, ids=t.index) / np.log(2),
                "observed_features": alpha_diversity("observed_features", R, ids=t.index),
                "chao1": alpha_diversity("chao1", R, ids=t.index),
                "pielou_e": alpha_diversity("pielou_e", R, ids=t.index)}))
        mean = sum(reps) / ITER
        for metric in mean.columns:
            groups = [g for g in stage.unique()]
            kw = kruskal(*[mean.loc[stage == g, metric] for g in groups]).pvalue if len(groups) > 1 else None
            rare.append({"marker": marker, "year": year, "metric": metric, "rarefied_depth": d,
                         "kruskal_p_rarefied": kw,
                         **{f"mean_{g}": mean.loc[stage == g, metric].mean() for g in groups}})

pd.DataFrame(corr).to_csv(f"{DATA}/alpha_depth_spearman.csv", index=False)
pd.DataFrame(rare).to_csv(f"{DATA}/alpha_rarefied_kruskal.csv", index=False)
```

---

## 4. Differential abundance analysis (ANCOM-BC)

```bash
# Filter genus tables by year
for YEAR in 2012 2013 2014 2015; do
    qiime feature-table filter-samples \
        --i-table ${DATA}/16S_genus_table.qza \
        --m-metadata-file metadata/mapping_16s.txt \
        --p-where "[Year]='${YEAR}'" \
        --p-min-frequency 1 \
        --o-filtered-table ${DATA}/16S_genus_table_${YEAR}.qza
done

for YEAR in 2014 2015; do
    qiime feature-table filter-samples \
        --i-table ${DATA}/ITS_genus_table.qza \
        --m-metadata-file metadata/mapping_bits.txt \
        --p-where "[Year]='${YEAR}'" \
        --p-min-frequency 1 \
        --o-filtered-table ${DATA}/ITS_genus_table_${YEAR}.qza
done

# Prevalence filter: genera present in at least three samples (3 / number of samples)
# (associative array: zsh or bash >= 4)
declare -A PRV=( [16S_2012]=0.25 [16S_2013]=0.5 [16S_2014]=0.25 [16S_2015]=0.0638
                 [ITS_2014]=0.3 [ITS_2015]=0.0714 )
declare -A MAP=( [16S]=metadata/mapping_16s.txt [ITS]=metadata/mapping_bits.txt )

# Bacteria ANCOM-BC per year; fungi 2014 and 2015 only
for KEY in 16S_2012 16S_2013 16S_2014 16S_2015 ITS_2014 ITS_2015; do
    MARKER=${KEY%_*}
    YEAR=${KEY#*_}

    qiime composition ancombc \
        --i-table ${DATA}/${MARKER}_genus_table_${YEAR}.qza \
        --m-metadata-file ${MAP[$MARKER]} \
        --p-formula Stage \
        --p-reference-levels Stage::0_control \
        --p-p-adj-method BH \
        --p-prv-cut ${PRV[$KEY]} \
        --p-lib-cut 0 \
        --p-tol 1e-5 \
        --p-max-iter 100 \
        --p-no-conserve \
        --p-alpha 0.05 \
        --o-differentials ${DATA}/${MARKER}_ancombc_${YEAR}.qza

    qiime tools export \
        --input-path ${DATA}/${MARKER}_ancombc_${YEAR}.qza \
        --output-path ${DATA}/${MARKER}_ancombc_export_${YEAR}
done
```

---

## 5. Core microbiome

```bash
for STAGE in 0_control 1_S1 2_S2 3_S3; do
    # Bacteria
    qiime feature-table filter-samples \
        --i-table ${DATA}/16S_table_bacteria.qza \
        --m-metadata-file metadata/mapping_16s.txt \
        --p-where "[Stage]='${STAGE}'" \
        --o-filtered-table ${DATA}/16S_table_stage_${STAGE}.qza

    qiime feature-table core-features \
        --i-table ${DATA}/16S_table_stage_${STAGE}.qza \
        --p-min-fraction 0.5 \
        --o-visualization ${DATA}/16S_core_${STAGE}.qzv

    qiime tools export \
        --input-path ${DATA}/16S_core_${STAGE}.qzv \
        --output-path ${DATA}/16S_core_export_${STAGE}

    # Fungi
    qiime feature-table filter-samples \
        --i-table ${DATA}/ITS_table_final.qza \
        --m-metadata-file metadata/mapping_bits.txt \
        --p-where "[Stage]='${STAGE}'" \
        --o-filtered-table ${DATA}/ITS_table_stage_${STAGE}.qza

    qiime feature-table core-features \
        --i-table ${DATA}/ITS_table_stage_${STAGE}.qza \
        --p-min-fraction 0.5 \
        --o-visualization ${DATA}/ITS_core_${STAGE}.qzv

    qiime tools export \
        --input-path ${DATA}/ITS_core_${STAGE}.qzv \
        --output-path ${DATA}/ITS_core_export_${STAGE}
done
```

---

## 6. Export QIIME 2 artifacts for R

```bash
# Genus and family tables
for TABLE in 16S_genus_table ITS_genus_table ITS_genus_table_noBotrytis ITS_family_table; do
    qiime tools export \
        --input-path ${DATA}/${TABLE}.qza \
        --output-path ${DATA}/${TABLE}_export

    biom convert \
        -i ${DATA}/${TABLE}_export/feature-table.biom \
        -o ${DATA}/${TABLE}_export/feature-table.tsv \
        --to-tsv
done

# Taxonomy
qiime tools export \
    --input-path ${DATA}/16S_taxonomy.qza \
    --output-path ${DATA}/16S_taxonomy_export

qiime tools export \
    --input-path ${DATA}/ITS_taxonomy_botrytis.qza \
    --output-path ${DATA}/ITS_taxonomy_export

# Alpha diversity
for YEAR in 2012 2013 2014 2015; do
    for METRIC in shannon observed_features pielou_e chao1 faith; do
        qiime tools export \
            --input-path ${DATA}/16S_${METRIC}_${YEAR}.qza \
            --output-path ${DATA}/16S_${METRIC}_export_${YEAR}
    done
    for METRIC in shannon observed_features pielou_e chao1; do
        qiime tools export \
            --input-path ${DATA}/ITS_${METRIC}_${YEAR}.qza \
            --output-path ${DATA}/ITS_${METRIC}_export_${YEAR}
    done
done

# PCoA
for YEAR in 2012 2013 2014 2015; do
    for METRIC in weighted_unifrac unweighted_unifrac; do
        qiime tools export \
            --input-path ${DATA}/16S_pcoa_${METRIC}_${YEAR}.qza \
            --output-path ${DATA}/16S_pcoa_${METRIC}_export_${YEAR}
    done
    for SET in genus_table genus_table_noBotrytis; do
        qiime tools export \
            --input-path ${DATA}/ITS_pcoa_${SET}_${YEAR}.qza \
            --output-path ${DATA}/ITS_pcoa_${SET}_export_${YEAR}
    done
done
```

---

## 7. Python — Correlation network analysis

```python
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import qiime2
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

DATA    = "path/to/data"
RESULTS = "path/to/results"
AAB = {"Gluconobacter", "Komagataeibacter", "Acetobacter", "Kozakia", "Gluconacetobacter", "Asaia"}


def genus(tax):
    g = tax.split(";")[-1].strip().replace("g__", "")
    return "Unclassified" if g in ("", "__") else g


# Relative abundance of genera (reads without a genus assignment and Escherichia-Shigella removed;
# fungi other than Botrytis); genera in >= 30 % of samples with mean relative abundance >= 0.01 %
def matrix(table, ids):
    df = qiime2.Artifact.load(f"{DATA}/{table}.qza").view(pd.DataFrame)
    df = df.loc[[i for i in ids if i in df.index]]
    g = df.columns.map(genus)
    keep = ~g.str.contains("uncultured|Unclassified|unassigned|Escherichia-Shigella", case=False) & ~g.str.endswith("__")
    df = df.loc[:, keep]
    rel = df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) * 100
    rel = rel.T.groupby(rel.columns.map(genus), sort=False).sum().T
    return rel.loc[:, ((rel > 0).mean() >= 0.3) & (rel.mean() >= 0.01)]


# Spearman correlations, Benjamini-Hochberg correction per network
def correlate(X, Y, same):
    rows = []
    for i, a in enumerate(X.columns):
        for j, b in enumerate(Y.columns):
            if same and j <= i:
                continue
            r, p = spearmanr(X[a], Y[b])
            rows.append({"from": a, "to": b, "r": r, "p": p})
    e = pd.DataFrame(rows).dropna()
    e["q"] = multipletests(e["p"], method="fdr_bh")[1]
    return e


# Samples with both 16S and ITS data
s16 = qiime2.Artifact.load(f"{DATA}/16S_genus_table.qza").view(pd.DataFrame).index
sits = qiime2.Artifact.load(f"{DATA}/ITS_genus_table_noBotrytis.qza").view(pd.DataFrame).index
ids = sorted(set(s16) & set(sits))

B = matrix("16S_genus_table", ids)
F = matrix("ITS_genus_table_noBotrytis", ids)
print(f"samples {len(ids)}, bacterial genera {B.shape[1]}, fungal genera {F.shape[1]}")

# Bacterial and cross-kingdom networks: |r| >= 0.3, q < 0.05
eb = correlate(B, B, True)
ec = correlate(B, F, False)
eb = eb[(eb.q < 0.05) & (eb.r.abs() >= 0.3)]
ec = ec[(ec.q < 0.05) & (ec.r.abs() >= 0.3)]
eb.to_csv(f"{RESULTS}/network_bacteria_edges.csv", index=False)
ec.to_csv(f"{RESULTS}/network_crosskingdom_edges.csv", index=False)

mean_abund = {**B.mean().to_dict(), **F.mean().to_dict()}
kingdom = {**{g: "Bacteria" for g in B.columns}, **{g: "Fungi" for g in F.columns}}


def draw(ax, edges, title):
    G = nx.Graph()
    for _, d in edges.iterrows():
        G.add_edge(d["from"], d["to"], r=d["r"])
    pos = nx.kamada_kawai_layout(G) if G.number_of_nodes() > 2 else nx.spring_layout(G, seed=1)
    for u, v, d in G.edges(data=True):
        ax.plot(*zip(pos[u], pos[v]), color="#D55E00" if d["r"] > 0 else "#0072B2",
                lw=1 + 7 * (abs(d["r"]) - 0.3), alpha=0.7, zorder=1)
    for n in G.nodes:
        ax.scatter(*pos[n], s=100 + 40 * np.sqrt(mean_abund[n]), zorder=2, edgecolor="white",
                   color="#4CAF50" if kingdom[n] == "Bacteria" else "#984EA3",
                   marker="D" if (kingdom[n] == "Bacteria" and n in AAB) else "o")
        ax.annotate(n, pos[n], xytext=(6, 6), textcoords="offset points", style="italic", fontsize=9)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.axis("off")


fig, axes = plt.subplots(2, 1, figsize=(12, 11), gridspec_kw={"height_ratios": [1, 2.4]})
draw(axes[0], eb, "A  Bacteria")
draw(axes[1], ec, "B  Bacteria and fungi")
fig.savefig(f"{RESULTS}/Fig7_correlation_networks.pdf", bbox_inches="tight")
fig.savefig(f"{RESULTS}/Fig7_correlation_networks.png", dpi=300, bbox_inches="tight")
```

---

## Citation

If you use this code, please cite:

> Cantu D., Blanco-Ulate B., Allen G., Figueroa-Balderas R. (2026).
> Botrytis cinerea infection reshapes the grape berry microbiome during noble rot. bioRxiv.
> https://doi.org/10.64898/2026.05.29.728822

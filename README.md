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
install.packages(c("tidyverse", "psych", "igraph", "ggraph"))
```

### QIIME 2 plugins
- dada2
- feature-table
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

---

## 1. Raw read preprocessing

```bash
conda activate qiime2

RAW=path/to/raw/fastq
DATA=path/to/output

# Extract barcodes and trim reads for each sequencing run
for RUN in run1 run2; do
    fastx_trimmer -f 1 -l 8 -Q33 \
        -i ${RAW}/${RUN}_R1.fastq \
        -o ${DATA}/${RUN}_barcodes.fastq

    fastx_trimmer -f 9 -Q33 \
        -i ${RAW}/${RUN}_R1.fastq \
        -o ${DATA}/${RUN}_R1_trimmed.fastq

    fastx_trimmer -f 9 -Q33 \
        -i ${RAW}/${RUN}_R2.fastq \
        -o ${DATA}/${RUN}_R2_trimmed.fastq
done
```

---

## 2. QIIME 2 pipeline

### 2.1 Import and demultiplex

```bash
# Import reads (repeat for each run)
qiime tools import \
    --type EMPPairedEndSequences \
    --input-path ${DATA}/run1_import/ \
    --output-path ${DATA}/16S_run1.qza

# Demultiplex (8-nt barcodes, no Golay error correction)
qiime demux emp-paired \
    --i-seqs ${DATA}/16S_run1.qza \
    --m-barcodes-file metadata/mapping_16S.txt \
    --m-barcodes-column BarcodeSequence \
    --p-no-golay-error-correction \
    --o-per-sample-sequences ${DATA}/16S_run1_demux.qza \
    --o-error-correction-details ${DATA}/16S_run1_demux_details.qza
```

### 2.2 Denoise with DADA2

```bash
# 16S (no truncation)
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${DATA}/16S_run1_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --o-table ${DATA}/16S_run1_table.qza \
    --o-representative-sequences ${DATA}/16S_run1_repseqs.qza \
    --o-denoising-stats ${DATA}/16S_run1_stats.qza

# ITS (no truncation)
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${DATA}/ITS_run1_demux.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --o-table ${DATA}/ITS_run1_table.qza \
    --o-representative-sequences ${DATA}/ITS_run1_repseqs.qza \
    --o-denoising-stats ${DATA}/ITS_run1_stats.qza
```

### 2.3 Merge runs

```bash
# 16S: two runs with same samples, merge by summing
qiime feature-table merge \
    --i-tables ${DATA}/16S_run1_table.qza \
    --i-tables ${DATA}/16S_run2_table.qza \
    --p-overlap-method sum \
    --o-merged-table ${DATA}/16S_merged_table.qza

qiime feature-table merge-seqs \
    --i-data ${DATA}/16S_run1_repseqs.qza \
    --i-data ${DATA}/16S_run2_repseqs.qza \
    --o-merged-data ${DATA}/16S_merged_repseqs.qza

# ITS: two runs with different samples
qiime feature-table merge \
    --i-tables ${DATA}/ITS_run1_table.qza \
    --i-tables ${DATA}/ITS_run2_table.qza \
    --p-overlap-method sum \
    --o-merged-table ${DATA}/ITS_merged_table.qza

qiime feature-table merge-seqs \
    --i-data ${DATA}/ITS_run1_repseqs.qza \
    --i-data ${DATA}/ITS_run2_repseqs.qza \
    --o-merged-data ${DATA}/ITS_merged_repseqs.qza
```

### 2.4 Taxonomic classification

```bash
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
```

### 2.5 Filter contaminants and low-depth samples

```bash
# Remove chloroplasts and mitochondria from 16S
qiime taxa filter-table \
    --i-table ${DATA}/16S_merged_table.qza \
    --i-taxonomy ${DATA}/16S_taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table ${DATA}/16S_merged_table_filtered.qza

# Remove 2015 S3 bacterial samples with fewer than 5000 reads
qiime feature-table filter-samples \
    --i-table ${DATA}/16S_merged_table_filtered.qza \
    --m-metadata-file metadata/mapping_16S.txt \
    --p-where "[Year] != '2015' OR [Stage] != '3_S3'" \
    --o-filtered-table ${DATA}/16S_table_final.qza
```

### 2.6 Phylogenetic tree (bacteria only)

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences ${DATA}/16S_merged_repseqs.qza \
    --o-alignment ${DATA}/16S_aligned_repseqs.qza \
    --o-masked-alignment ${DATA}/16S_masked_aligned_repseqs.qza \
    --o-tree ${DATA}/16S_unrooted_tree.qza \
    --o-rooted-tree ${DATA}/16S_rooted_tree.qza
```

### 2.7 Collapse to genus and family level

```bash
# 16S genus (level 6)
qiime taxa collapse \
    --i-table ${DATA}/16S_table_final.qza \
    --i-taxonomy ${DATA}/16S_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${DATA}/16S_genus_table.qza

# ITS genus (level 6)
qiime taxa collapse \
    --i-table ${DATA}/ITS_merged_table.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${DATA}/ITS_genus_table.qza

# ITS family (level 5)
qiime taxa collapse \
    --i-table ${DATA}/ITS_merged_table.qza \
    --i-taxonomy ${DATA}/ITS_taxonomy.qza \
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
        --m-metadata-file metadata/mapping_16S.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/16S_table_${YEAR}.qza

    qiime diversity alpha-phylogenetic \
        --i-table ${DATA}/16S_table_${YEAR}.qza \
        --i-phylogeny ${DATA}/16S_rooted_tree.qza \
        --p-metric faith_pd \
        --o-alpha-diversity ${DATA}/16S_faith_${YEAR}.qza

    for METRIC in shannon observed_features evenness chao1; do
        qiime diversity alpha \
            --i-table ${DATA}/16S_table_${YEAR}.qza \
            --p-metric ${METRIC} \
            --o-alpha-diversity ${DATA}/16S_${METRIC}_${YEAR}.qza
    done

    # Fungi - filter by year
    qiime feature-table filter-samples \
        --i-table ${DATA}/ITS_merged_table.qza \
        --m-metadata-file metadata/mapping_ITS.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/ITS_table_${YEAR}.qza

    for METRIC in shannon observed_features evenness chao1; do
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
            --m-metadata-file metadata/mapping_16S.txt \
            --m-metadata-column Stage \
            --p-method permanova \
            --p-permutations 999 \
            --o-visualization ${DATA}/16S_permanova_${METRIC}_${YEAR}.qzv
    done

    # Fungi - Bray-Curtis
    qiime diversity beta \
        --i-table ${DATA}/ITS_table_${YEAR}.qza \
        --p-metric braycurtis \
        --o-distance-matrix ${DATA}/ITS_braycurtis_${YEAR}.qza

    qiime diversity pcoa \
        --i-distance-matrix ${DATA}/ITS_braycurtis_${YEAR}.qza \
        --o-pcoa ${DATA}/ITS_pcoa_${YEAR}.qza

    qiime diversity beta-group-significance \
        --i-distance-matrix ${DATA}/ITS_braycurtis_${YEAR}.qza \
        --m-metadata-file metadata/mapping_ITS.txt \
        --m-metadata-column Stage \
        --p-method permanova \
        --p-permutations 999 \
        --o-visualization ${DATA}/ITS_permanova_${YEAR}.qzv
done
```

## 3.3 Alpha diversity group significance

```bash
# Kruskal-Wallis test of alpha diversity across noble rot stages, per metric and year.
# The visualization also reports pairwise comparisons; pairwise Kruskal-Wallis on two
# groups is equivalent to the Wilcoxon rank-sum test used for the figure letters.
#
# Note: 3.1 writes Faith's PD as 16S_faith_${YEAR}.qza, so "faith" is used here rather
# than the QIIME metric name faith_pd.

for YEAR in 2012 2013 2014 2015; do
    for METRIC in shannon observed_features evenness chao1 faith; do
        qiime diversity alpha-group-significance \
            --i-alpha-diversity ${DATA}/16S_${METRIC}_${YEAR}.qza \
            --m-metadata-file metadata/mapping_16S.txt \
            --o-visualization ${DATA}/16S_${METRIC}_${YEAR}_significance.qzv
    done
done

# Fungi: no phylogenetic metric
for YEAR in 2012 2013 2014 2015; do
    for METRIC in shannon observed_features evenness chao1; do
        qiime diversity alpha-group-significance \
            --i-alpha-diversity ${DATA}/ITS_${METRIC}_${YEAR}.qza \
            --m-metadata-file metadata/mapping_ITS.txt \
            --o-visualization ${DATA}/ITS_${METRIC}_${YEAR}_significance.qzv
    done
done

# Export the per-sample vectors for the summary table below
for MARKER in 16S ITS; do
    for YEAR in 2012 2013 2014 2015; do
        for METRIC in shannon observed_features evenness chao1 faith; do
            IN=${DATA}/${MARKER}_${METRIC}_${YEAR}.qza
            [ -f "${IN}" ] || continue
            qiime tools export \
                --input-path "${IN}" \
                --output-path ${DATA}/export/${MARKER}_${METRIC}_${YEAR}
        done
    done
done
```

## 3.4 Alpha diversity summary table

Produces the per-year, per-metric table of `mean ± SD` by stage with the Kruskal-Wallis
p-value — the values reported in the Results and in Supplementary Tables S2 and S3.

```r
# Alpha diversity summary: mean +/- SD per stage, with Kruskal-Wallis p-values.
DATA <- "path/to/output"

read_mapping <- function(path) {
  hdr <- strsplit(readLines(path, n = 1), "\t")[[1]]
  hdr[1] <- sub("^#", "", hdr[1])            # "#SampleID" -> "SampleID"
  m <- read.delim(path, sep = "\t", header = FALSE, skip = 2,   # skip header + #q2:types
                  col.names = hdr, check.names = FALSE, colClasses = "character")
  names(m)[1] <- "SampleID"
  m
}

alpha_table <- function(prefix, metrics, mapping,
                        years = c(2012, 2013, 2014, 2015)) {
  map  <- read_mapping(mapping)
  rows <- list()
  for (y in years) for (met in metrics) {
    f <- file.path(DATA, "export", sprintf("%s_%s_%s", prefix, met, y),
                   "alpha-diversity.tsv")
    if (!file.exists(f)) next                 # stage/year combinations not sampled
    a <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE)
    names(a)[1:2] <- c("SampleID", "value")
    d <- merge(a, map, by = "SampleID")
    d$Stage <- as.factor(d$Stage)

    stats <- tapply(d$value, d$Stage, function(v)
               sprintf("%s ± %s", round(mean(v), 3), round(sd(v), 3)))
    p <- if (nlevels(droplevels(d$Stage)) > 1)
           round(kruskal.test(value ~ Stage, data = d)$p.value, 4) else NA

    row <- list(year = y, metric = met,
                n = round(mean(table(d$Stage)[table(d$Stage) > 0])))
    for (s in c("0_control", "1_S1", "2_S2", "3_S3"))
      row[[s]] <- if (s %in% names(stats)) stats[[s]] else NA
    row$kruskal_p <- p
    rows[[length(rows) + 1]] <- as.data.frame(row, check.names = FALSE)
  }
  do.call(rbind, rows)
}

bact <- alpha_table("16S", c("shannon", "observed", "chao1", "faith", "pielou"),
                    "metadata/mapping_16S.txt")
fung <- alpha_table("ITS", c("shannon", "observed", "chao1", "pielou"),
                    "metadata/mapping_ITS.txt")

write.csv(bact, file.path(DATA, "16S_alpha_table.csv"), row.names = FALSE)
write.csv(fung, file.path(DATA, "ITS_alpha_table.csv"), row.names = FALSE)
```

---

## 4. Differential abundance analysis (ANCOM-BC)

```bash
# Filter genus tables by year
for YEAR in 2012 2013 2014 2015; do
    qiime feature-table filter-samples \
        --i-table ${DATA}/16S_genus_table.qza \
        --m-metadata-file metadata/mapping_16S.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/16S_genus_table_${YEAR}.qza
done

for YEAR in 2014 2015; do
    qiime feature-table filter-samples \
        --i-table ${DATA}/ITS_genus_table.qza \
        --m-metadata-file metadata/mapping_ITS.txt \
        --p-where "[Year]='${YEAR}'" \
        --o-filtered-table ${DATA}/ITS_genus_table_${YEAR}.qza
done

# Bacteria ANCOM-BC per year
for YEAR in 2012 2013 2014 2015; do
    qiime composition ancombc \
        --i-table ${DATA}/16S_genus_table_${YEAR}.qza \
        --m-metadata-file metadata/mapping_16S.txt \
        --p-formula Stage \
        --p-reference-levels Stage::0_control \
        --o-differentials ${DATA}/16S_ancombc_${YEAR}.qza

    qiime tools export \
        --input-path ${DATA}/16S_ancombc_${YEAR}.qza \
        --output-path ${DATA}/16S_ancombc_export_${YEAR}
done

# Fungi ANCOM-BC (2014 and 2015 only)
for YEAR in 2014 2015; do
    qiime composition ancombc \
        --i-table ${DATA}/ITS_genus_table_${YEAR}.qza \
        --m-metadata-file metadata/mapping_ITS.txt \
        --p-formula Stage \
        --p-reference-levels Stage::0_control \
        --o-differentials ${DATA}/ITS_ancombc_${YEAR}.qza

    qiime tools export \
        --input-path ${DATA}/ITS_ancombc_${YEAR}.qza \
        --output-path ${DATA}/ITS_ancombc_export_${YEAR}
done
```

---

## 5. Core microbiome

```bash
for STAGE in 0_control 1_S1 2_S2 3_S3; do
    # Bacteria
    qiime feature-table filter-samples \
        --i-table ${DATA}/16S_table_final.qza \
        --m-metadata-file metadata/mapping_16S.txt \
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
        --i-table ${DATA}/ITS_merged_table.qza \
        --m-metadata-file metadata/mapping_ITS.txt \
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
for TABLE in 16S_genus_table ITS_genus_table ITS_family_table; do
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
    --input-path ${DATA}/ITS_taxonomy.qza \
    --output-path ${DATA}/ITS_taxonomy_export

# Alpha diversity
for YEAR in 2012 2013 2014 2015; do
    for METRIC in shannon observed_features evenness chao1 faith_pd; do
        qiime tools export \
            --input-path ${DATA}/16S_${METRIC}_${YEAR}.qza \
            --output-path ${DATA}/16S_${METRIC}_export_${YEAR}
    done
    for METRIC in shannon observed_features evenness chao1; do
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
    qiime tools export \
        --input-path ${DATA}/ITS_pcoa_${YEAR}.qza \
        --output-path ${DATA}/ITS_pcoa_export_${YEAR}
done
```

---

## 7. R — Co-occurrence network analysis

```r
library(tidyverse)
library(psych)
library(igraph)
library(ggraph)

DATA    <- "path/to/data"
RESULTS <- "path/to/results"

# Load metadata
metadata <- read.table(file.path(DATA, "mapping_16S.txt"),
                       header=TRUE, sep="\t", comment.char="",
                       check.names=FALSE) %>%
  rename(SampleID=`#SampleID`) %>%
  filter(!grepl("#q2:types", SampleID))

metadata_its <- read.table(file.path(DATA, "mapping_ITS.txt"),
                           header=TRUE, sep="\t", comment.char="",
                           check.names=FALSE) %>%
  rename(SampleID=`#SampleID`) %>%
  filter(!grepl("#q2:types", SampleID))

# Helper: parse genus from taxonomy string
extract_genus <- function(tax) {
  parts <- strsplit(as.character(tax), ";")[[1]]
  g <- trimws(tail(parts, 1))
  g <- gsub("g__", "", g)
  if (g == "" | g == "__") "Unclassified" else g
}

# Load genus tables
tax_16s <- read.delim(
  file.path(DATA, "16S_genus_table_export", "feature-table.tsv"),
  skip=1, check.names=FALSE) %>%
  rename(taxonomy=`#OTU ID`)

tax_its <- read.delim(
  file.path(DATA, "ITS_genus_table_export", "feature-table.tsv"),
  skip=1, check.names=FALSE) %>%
  rename(taxonomy=`#OTU ID`)

# Shared samples between 16S and ITS datasets
shared_samples <- intersect(metadata$SampleID, metadata_its$SampleID)
cat("Shared samples:", length(shared_samples), "\n")

# Prepare relative abundance matrix
# Filters: minimum prevalence >= 30%, minimum mean abundance >= 0.01%
prepare_matrix <- function(tax_df, meta,
                            min_prevalence=0.3, min_mean_abund=0.01) {
  long <- tax_df %>%
    pivot_longer(-taxonomy, names_to="SampleID", values_to="count") %>%
    inner_join(meta, by="SampleID") %>%
    filter(!is.na(Stage)) %>%
    mutate(genus=sapply(taxonomy, extract_genus)) %>%
    filter(!grepl("uncultured|Unclassified|unassigned|__$",
                  genus, ignore.case=TRUE))

  long <- long %>%
    group_by(SampleID) %>%
    mutate(total     = sum(count),
           rel_abund = ifelse(total > 0, count / total * 100, 0)) %>%
    ungroup() %>%
    group_by(SampleID, genus) %>%
    summarize(rel_abund=sum(rel_abund), .groups="drop")

  prevalent <- long %>%
    group_by(genus) %>%
    summarize(prevalence=mean(rel_abund > 0, na.rm=TRUE),
              mean_abund=mean(rel_abund,     na.rm=TRUE), .groups="drop") %>%
    filter(prevalence >= min_prevalence, mean_abund >= min_mean_abund) %>%
    pull(genus)

  mat <- long %>%
    filter(genus %in% prevalent) %>%
    pivot_wider(names_from=genus, values_from=rel_abund, values_fill=0)
  rnames <- mat$SampleID
  mat    <- as.data.frame(mat[, -1])
  rownames(mat) <- rnames
  mat
}

# Build matrices for shared samples only
mat_16s <- prepare_matrix(tax_16s,
                          metadata     %>% filter(SampleID %in% shared_samples))
mat_its <- prepare_matrix(tax_its,
                          metadata_its %>% filter(SampleID %in% shared_samples))

# Align rows
common_rows <- intersect(rownames(mat_16s), rownames(mat_its))
mat_16s     <- mat_16s[common_rows, ]
mat_its     <- mat_its[common_rows, ]

cat("Bacterial genera:", ncol(mat_16s), "\n")
cat("Fungal genera:",    ncol(mat_its), "\n")
cat("Samples used:",     length(common_rows), "\n")

# Spearman correlations with Benjamini-Hochberg correction
cor_16s_test   <- corr.test(mat_16s,          method="spearman", adjust="BH")
cor_cross_test <- corr.test(mat_16s, mat_its,  method="spearman", adjust="BH")

cat("Bacterial significant pairs (q<0.05, |r|>=0.5):",
    sum(cor_16s_test$p < 0.05 &
        abs(cor_16s_test$r) >= 0.5, na.rm=TRUE), "\n")
cat("Cross-kingdom significant pairs (q<0.05, |r|>=0.3):",
    sum(cor_cross_test$p < 0.05, na.rm=TRUE), "\n")

# Mean abundance for node sizing
mean_abund_16s <- colMeans(mat_16s)
mean_abund_its <- colMeans(mat_its)
all_abund      <- c(mean_abund_16s, mean_abund_its)

# AAB genera for node shape coding
aab_genera <- c("Gluconobacter", "Komagataeibacter", "Acetobacter",
                "Kozakia", "Gluconacetobacter", "Asaia")

# Build igraph network from correlation matrix
build_network <- function(cor_mat, p_mat, min_cor=0.3, max_p=0.05) {
  cor_long <- as.data.frame(as.table(cor_mat)) %>%
    rename(from=Var1, to=Var2, correlation=Freq) %>%
    filter(as.character(from) < as.character(to)) %>%
    left_join(as.data.frame(as.table(p_mat)) %>%
                rename(from=Var1, to=Var2, pval=Freq),
              by=c("from","to")) %>%
    filter(pval < max_p, abs(correlation) >= min_cor)

  if (nrow(cor_long) == 0) return(NULL)

  g <- graph_from_data_frame(cor_long, directed=FALSE)
  E(g)$weight     <- cor_long$correlation
  E(g)$positive   <- cor_long$correlation > 0
  E(g)$abs_weight <- abs(cor_long$correlation)
  g
}

# Bacterial network: |r| >= 0.5, q < 0.05
g_16s   <- build_network(cor_16s_test$r,   cor_16s_test$p,   min_cor=0.5)

# Cross-kingdom network: |r| >= 0.3, q < 0.05
g_cross <- build_network(cor_cross_test$r, cor_cross_test$p, min_cor=0.3)

# Add node attributes - bacteria
V(g_16s)$kingdom    <- "Bacteria"
V(g_16s)$is_aab     <- V(g_16s)$name %in% aab_genera
V(g_16s)$mean_abund <- mean_abund_16s[V(g_16s)$name]

# Add node attributes - cross-kingdom
V(g_cross)$kingdom    <- ifelse(V(g_cross)$name %in% colnames(mat_16s),
                                "Bacteria", "Fungi")
V(g_cross)$is_aab     <- V(g_cross)$name %in% aab_genera
V(g_cross)$mean_abund <- ifelse(V(g_cross)$name %in% names(all_abund),
                                all_abund[V(g_cross)$name], 1)

# Network plot function
plot_network <- function(g, title) {
  ggraph(g, layout="fr", weights=E(g)$abs_weight) +
    geom_edge_link(aes(color=positive, width=abs(weight)), alpha=0.7) +
    geom_node_point(aes(color=kingdom, shape=is_aab, size=mean_abund)) +
    geom_node_text(aes(label=name), repel=TRUE, size=3, fontface="italic") +
    scale_edge_color_manual(
      values=c("TRUE"="#D55E00", "FALSE"="#0072B2"),
      labels=c("TRUE"="Positive", "FALSE"="Negative"),
      name="Correlation") +
    scale_edge_width_continuous(range=c(0.3, 2), name="|r|") +
    scale_color_manual(
      values=c("Bacteria"="#4CAF50", "Fungi"="#984EA3"),
      name="Kingdom") +
    scale_shape_manual(
      values=c("TRUE"=18, "FALSE"=16),
      labels=c("TRUE"="AAB", "FALSE"="Other"),
      name="Group") +
    scale_size_continuous(range=c(3, 8), name="Mean\nabundance (%)") +
    labs(title=title) +
    theme_graph(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
}

p_net_16s   <- plot_network(g_16s,   "A  Bacterial co-occurrence network")
p_net_cross <- plot_network(g_cross, "B  Cross-kingdom co-occurrence network")

fig7 <- p_net_16s / p_net_cross

pdf(file.path(RESULTS, "Fig7_cooccurrence_networks.pdf"), width=12, height=12)
print(fig7)
dev.off()

ggsave(file.path(RESULTS, "Fig7_cooccurrence_networks.png"),
       fig7, width=12, height=12, dpi=300)

# Save correlation matrices
write.csv(as.data.frame(cor_16s_test$r),
          file.path(RESULTS, "cooccurrence_16S_rho.csv"))
write.csv(as.data.frame(cor_16s_test$p),
          file.path(RESULTS, "cooccurrence_16S_qval.csv"))
write.csv(as.data.frame(cor_cross_test$r),
          file.path(RESULTS, "cooccurrence_cross_rho.csv"))
write.csv(as.data.frame(cor_cross_test$p),
          file.path(RESULTS, "cooccurrence_cross_qval.csv"))

cat("Co-occurrence analysis complete.\n")
```

---

## Citation

If you use this code, please cite:

> Cantu D., Blanco-Ulate B., Allen G., Figueroa-Balderas R. (2026).
> Botrytis cinerea infection reshapes the grape berry microbiome during noble rot. bioRxiv - to complete with DOI

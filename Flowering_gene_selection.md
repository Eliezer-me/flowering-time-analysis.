# Flowering Gene  Selection in Chamaecrista fasciculata

This repository documents a reproducible, conservative, and biologically grounded pipeline for identifying putative orthologous genes in Chamaecrista fasciculata using reciprocal best hit (RBH) comparisons against Arabidopsis thaliana and Glycine max (soybean), and for integrating these orthologs with selective sweep signals detected using SweeD. It also contains PCA analyses of Chamaecrista populations (by state and by latitude), a US map showing sample coverage, and space reserved for figures.

The workflow was developed through iterative analysis, validation, and debugging, and is suitable for:

- Evolutionary genomics
- Functional gene annotation in non-model species
- Detecting flowering-time genes under recent positive selection. 

---

## Introduction & Biological Motivation

Chamaecrista fasciculata is adapted to diverse environments, with photoperiod sensitivity varying across latitudes.

- Annual legume, native to North America
- Locally adapted to regional photoperiods and thermal environments
- The survival and reproductive success of C. fasciculata depends on its ability to flower at the right time.

Understanding its adaptive strategies is essential for conservation in a changing climate. Flowering time affects:
- reproductive success
- gene flow
- adaptation
- survival and fitness

Challenges and caveats:
- Phenological mismatches are a key concern: populations adapted to a given photoperiod may suffer if climate shifts flowering schedules.
- Gene flow between populations can blur adaptive differences.
- Flowering-time gene divergence can lead to mismatched flowering schedules with fitness consequences.

---

## Repository overview

This repo contains:

- RBH-based ortholog detection code and notes
- SweeD integration and post-processing to map CLR peaks to genes
- PCA scripts and plots for Chamaecrista populations (by state and by latitude)
- A US map showing states with sampled populations
- Final tables linking Chamaecrista genes to Arabidopsis and soybean orthologs with selection evidence
- Space for images/figures (see `docs/images/`)
- Guidance for reproducibility and environment recording

---

## BLAST database creation & searches (R script — use placeholder paths)

Below is the R script used to create BLAST databases and run `blastp`. Paths are placeholders: replace them with your local paths or variables loaded from a gitignored config.

```r
# Example R script used to create blast DBs and run blastp (paths are placeholders)
# Replace /path/to/... with your actual paths.

ath_proteins <- "/path/to/athaliana_proteins.fa"
chamaecrista_proteins <- "/path/to/Cfasciculatavar_ISC494698HAP2_855_v1.1.protein.fa"
gmax_proteins <- "/path/to/Gmax.protein.fa"

blast_db_cha <- "/path/to/chamaecrista_db"
blast_db_ath <- "/path/to/athaliana_db"
blast_db_gmax <- "/path/to/gmax_db"

ath2_out <- "/path/to/output/ath2_vs_chamaecrista.tsv"
cham2_ath_out <- "/path/to/output/cham2_vs_athaliana.tsv"
gmax2_out <- "/path/to/output/gmax2_vs_chamaecrista.tsv"
cham2_gmax_out <- "/path/to/output/cham2_vs_gmax.tsv"

blastp_bin <- "/path/to/blastp"
makeblastdb_bin <- "/path/to/makeblastdb"

# Create BLAST databases
system(paste(makeblastdb_bin, "-in", chamaecrista_proteins, "-dbtype prot -out", blast_db_cha))
system(paste(makeblastdb_bin, "-in", ath_proteins, "-dbtype prot -out", blast_db_ath))
system(paste(makeblastdb_bin, "-in", gmax_proteins, "-dbtype prot -out", blast_db_gmax))

# Arabidopsis proteins vs Chamaecrista proteins
cmd1 <- paste(
  blastp_bin, "-query", ath_proteins,
  "-db", blast_db_cha,
  "-evalue 1e-5",
  "-num_threads 4",
  "-out", ath2_out,
  "-outfmt '6 qseqid sseqid pident length evalue bitscore'"
)
system(cmd1)

# Chamaecrista proteins vs Arabidopsis proteins
cmd2 <- paste(
  blastp_bin, "-query", chamaecrista_proteins,
  "-db", blast_db_ath,
  "-evalue 1e-5",
  "-num_threads 4,
  "-out", cham2_ath_out,
  "-outfmt '6 qseqid sseqid pident length evalue bitscore'"
)
system(cmd2)

# Glycine max proteins vs Chamaecrista proteins
cmd3 <- paste(
  blastp_bin, "-query", gmax_proteins,
  "-db", blast_db_cha,
  "-evalue 1e-5",
  "-num_threads 4",
  "-out", gmax2_out,
  "-outfmt '6 qseqid sseqid pident length evalue bitscore'"
)
system(cmd3)

# Chamaecrista proteins vs Glycine max proteins
cmd4 <- paste(
  blastp_bin, "-query", chamaecrista_proteins,
  "-db", blast_db_gmax,
  "-evalue 1e-5",
  "-num_threads 4,
  "-out", cham2_gmax_out,
  "-outfmt '6 qseqid sseqid pident length evalue bitscore'"
)
system(cmd4)
```
---

## Isoform handling (critical)

To avoid inflated or misleading matches caused by transcript isoforms, collapse identifiers to gene-level IDs before RBH:

```r
clean_id <- function(ids) {
  gsub("\\.\\d+(\\.p)?$", "", ids)
}
```

This removes isoform suffixes like `.1`, `.2`, `.3`, or `.p`.

---

## Post-BLAST processing & RBH filtering (conceptual R steps)

- Read BLAST tabular results with columns: qseqid, sseqid, pident, length, evalue, bitscore
- Collapse isoform IDs using `clean_id()` for both query and subject
- For each query, keep the best hit (by bitscore or best evalue)
- Compute reciprocal best hits: keep pairs where query's best subject and subject's best query match each other
- Generate RBH sets for:
  - Chamaecrista ↔ Arabidopsis
  - Chamaecrista ↔ Glycine max
- Merge RBH results by Chamaecrista gene ID to produce the master ortholog table

Conceptual pseudocode:

```r
# read results, collapse ids, keep top hit per query, compute RBH by joining forward+reverse top-hits
```

---

## Selective sweeps (SweeD) processing and mapping to genes

SweeD produces per-chromosome `.sweed` files that are merged to `sweed_all` with columns: `Chr`, `Position`, `Likelihood`.

To focus on strong signals we used the top 10% CLR threshold:

```r
threshold <- quantile(sweed_all$Likelihood, 0.90, na.rm = TRUE)
```

Map sweeps to gene models (GFF3) using GenomicRanges and keep the maximum CLR per gene:
<img width="1259" height="778" alt="image" src="https://github.com/user-attachments/assets/05828bfd-ea55-4827-a28f-691a13999b3f" />

```r
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

sweed_gr <- GRanges(
  seqnames = sweed_all$Chr,
  ranges = IRanges(start = sweed_all$Position, end = sweed_all$Position),
  Likelihood = sweed_all$Likelihood
)

sweed_top <- sweed_gr[sweed_gr$Likelihood > threshold]

gff <- import.gff3("/path/to/Cfasciculatavar_ISC494698HAP2_855_v1.1.gene.gff3")
genes <- gff[gff$type == "gene"]

hits <- findOverlaps(sweed_top, genes)

overlap_info <- data.frame(
  gene_index = subjectHits(hits),
  Likelihood = mcols(sweed_top)[queryHits(hits), "Likelihood"]
)

max_likelihood_by_gene <- overlap_info %>%
  group_by(gene_index) %>%
  summarise(Max_Likelihood = max(Likelihood, na.rm = TRUE))

```

Key fix: if multiple sweep positions overlap the same gene, assign only the maximum CLR value to that gene.

---

## Integration with flowering-time genes

- Intersect genes with Max_CLR > threshold with curated flowering-time gene sets (Arabidopsis-based lists, literature lists).
- Annotate these Chamaecrista genes with RBH-derived Arabidopsis and soybean orthologs.
- Produce split, biologically interpretable tables:
  - `orthologs_ath_only_under_sweep.csv`
  - `orthologs_gmax_only_under_sweep.csv`
  - `orthologs_shared_under_sweep.csv`
- All orthologs reported are putative (RBH-based), not experimentally validated.

---

## PCA of Chamaecrista populations.

Two PCA plotting scripts are included (color by State and color by Latitude). The `.eigenvec` and environmental file paths below are placeholders — update to your local path or read from a config.

PCA colored by State:

```r
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Reload PCA data and extract State from FamilyID
pca_data <- read.table("path/to/chamaecrista_allStates_PCA0.5.eigenvec", header = FALSE, sep = "")
colnames(pca_data) <- c("FamilyID", "IndividualID", "PC1", "PC2", "PC3", "PC4", "PC5","PC6", "PC7", "PC8", "PC9", "PC10")
pca_data$State <- substr(pca_data$FamilyID, 1, 2)
pca_data$State <- as.factor(pca_data$State)

# Pick a distinct color palette (at least as many colors as states)
palette <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(pca_data$State)))

# Plot with custom colors
ggplot(pca_data, aes(x = PC1, y = PC2, color = State)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = palette) +
  theme_minimal(base_size = 10) +
  labs(
    title = "PCA of Chamaecrista_0.05%",
    x = "PC1",
    y = "PC2",
    color = "State"
  )
<img width="1400" height="865" alt="image" src="https://github.com/user-attachments/assets/2c953332-370c-4087-bc06-799e20e735ec" />
```
#scree plot to show what percentage each Pca contributes to the variations
```r
library(ggplot2)

# Load eigenvalues
eigenvals <- scan("learning/chamaecrista_pca.eigenval")

# Calculate proportions and cumulative variance
prop_var <- eigenvals / sum(eigenvals)
cum_var <- cumsum(prop_var)
pc_labels <- paste0("PC", 1:length(eigenvals))

# Create data frame
df <- data.frame(
  PC = factor(pc_labels, levels = pc_labels),
  Variance = prop_var,
  Cumulative = cum_var
)

# Scree plot with cumulative line and elbow annotation
ggplot(df, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "#377eb8", alpha = 0.8) +
  geom_line(aes(y = Cumulative), group = 1, color = "#e41a1c", size = 1.2) +
  geom_point(aes(y = Cumulative), color = "#e41a1c", size = 2) +
  geom_text(aes(y = Cumulative, label = paste0(round(Cumulative * 100, 1), "%")),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "gray40") +  # elbow around PC4
  annotate("text", x = 4.5, y = 0.8, label = "Elbow", angle = 90, vjust = -0.5, color = "purple") +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance"
    
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



```
<img width="1400" height="865" alt="image" src="https://github.com/user-attachments/assets/47b56479-fa4b-4ba3-a655-8415f656b340" />

PCA colored by Latitude:

```r
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)

# Step 1: Load PCA data
pca_data <- read.table("path/to/chamaecrista_allStates_PCA0.5.eigenvec", header = FALSE, sep = "")
colnames(pca_data) <- c("FamilyID", "IndividualID", "PC1", "PC2", "PC3", "PC4", "PC5","PC6", "PC7", "PC8", "PC9", "PC10")

# Step 2: Extract clean state codes (2-letter prefix from FamilyID)
pca_data$State <- sub("^([A-Z]{2})_.*", "\\1", pca_data$FamilyID)

# Optional: check for malformed IDs
# table(pca_data$State, useNA = "always")

# Step 3: Create a data frame of states and approximate latitudes
state_latitudes <- data.frame(
  State = c("AL", "AR", "FL", "GA", "IA", "IL", "IN", "KS", "LA", 
            "MN", "MO", "MS", "NC", "OH", "OK", "PA", "SC", "TN", "TX", "VA"),
  Latitude = c(32.3, 35.2, 27.7, 32.2, 41.9, 40.6, 39.8, 38.5, 30.9,
               46.7, 38.6, 32.4, 35.5, 40.4, 35.5, 41.2, 33.8, 35.9, 31.5, 37.5)
)

# Step 4: Merge latitude into PCA data
pca_data <- left_join(pca_data, state_latitudes, by = "State")

# Step 5: Plot PC1 vs PC2 colored by latitude
ggplot(pca_data, aes(x = PC1, y = PC2, color = Latitude)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "C") +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCA by Latitude",
    x = "PC1",
    y = "PC2",
    color = "Latitude"
  )
```
<img width="1400" height="865" alt="image" src="https://github.com/user-attachments/assets/0a5a7710-2f1c-43c9-a3ec-f61be318e20f" />

Notes from PCA:
- Populations cluster partly by state and latitude, consistent with regional adaptation and photoperiod differences.
- These patterns support downstream analyses of flowering-time divergence along latitudinal gradients.

---

## US map of sampling coverage

Below is the R code used to create a US map highlighting states with at least one sample. Update the file paths (`path/to/...`) to point to your `.eigenvec` and environmental tables (or keep relative `learning/...` if that reflects your project layout).

```r
library(ggplot2)
library(dplyr)
library(maps)
library(viridis)

# 1. Load PCA data (change path as needed)
pca_raw <- read.table("path/to/chamaecrista_allStates_PCA0.5.eigenvec", header = FALSE)
colnames(pca_raw) <- c("FamilyID", "IndividualID", paste0("PC", 1:10))

# 2. Load environmental data (with Latitude/Longitude columns)
env <- read.table("path/to/chamae_BioClim.tsv", header = TRUE, sep = "\t")

# 3. Join datasets
df <- left_join(pca_raw, env, by = c("IndividualID" = "IID"))

# 4. Add State column from FamilyID AFTER the join, so it definitely exists in df
df <- df %>% mutate(State = substr(FamilyID, 1, 2))

# 5. Filter rows with valid Lat/Long
df_geo <- df %>% filter(!is.na(Latitude), !is.na(Longitude))

# 6. Load US map data
us_map <- map_data("state")

# 7. Create state abbreviation to full lowercase name mapping
state_abbs <- data.frame(
  abb = state.abb,
  full = tolower(state.name)
)

# 8. Extract sample states and map to full names
sample_states <- df_geo %>%
  distinct(State) %>%
  left_join(state_abbs, by = c("State" = "abb")) %>%
  filter(!is.na(full)) %>%
  pull(full)

# 9. Flag states in map data with samples
us_map <- us_map %>%
  mutate(has_sample = region %in% sample_states)

# 10. Plot
ggplot(us_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = has_sample), color = "white") +
  scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "purple"),
                    name = "Sample Present",
                    labels = c("No", "Yes")) +
  coord_fixed(1.3) +
  theme_minimal(base_size = 11) +
  labs(
       subtitle = "At least one sample",
       fill = "Samples used")
```
<img width="1400" height="865" alt="image" src="https://github.com/user-attachments/assets/c4a22e8c-5969-432a-9124-50b6617aeb33" />
```r


library(ggplot2)
library(dplyr)
library(maps)
library(viridis)

# === 1. Load PCA and metadata ===
pca_raw <- read.table("learning/chamaecrista_allStates_PCA0.5.eigenvec", header = FALSE)
colnames(pca_raw) <- c("FamilyID", "IndividualID", paste0("PC", 1:10))
pca_raw$State <- substr(pca_raw$FamilyID, 1, 2)

# === 2. Load environmental data (assumed to have Latitude, Longitude) ===
env <- read.table("learning/chamae_BioClim.tsv", header = TRUE, sep = "\t")
df <- left_join(pca_raw, env, by = c("IndividualID" = "IID"))  # Adjust join key if needed

# === 3. Filter complete cases with coordinates ===
df_geo <- df %>% filter(!is.na(Latitude), !is.na(Longitude))

# === 4. Plot PC1 on US map ===
us_map <- map_data("state")

ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "white") +
  geom_point(data = df_geo, aes(x = Longitude, y = Latitude, color = PC1),
             size = 3, alpha = 0.8) +
  scale_color_viridis(option = "C", name = "PC1 Score") +
  coord_fixed(1.3) +
  theme_minimal(base_size = 11) 


```
<img width="1259" height="778" alt="image" src="https://github.com/user-attachments/assets/93a0b3e7-4f2f-4935-bde7-0454fe04f706" />

Map note:
- The map highlights U.S. states with at least one sampled individual in our dataset. It is intended as a quick visual summary of geographic sampling coverage and should be complemented with point maps (latitude/longitude) for precise sample locations and with overlays of environmental gradients when relevant.
- Ensure the `env` table contains columns named `IID`, `Latitude`, and `Longitude` (or adjust join keys accordingly).

---

## Evolutionary notes (summary)

This study examines how C. fasciculata populations vary across latitude in traits related to flowering time. Key points:

- Photoperiod sensitivity varies across latitudes and is an important axis of adaptation.
- Flowering time influences reproductive success and gene flow; mismatches in phenology can reduce fitness.
- Detecting selection on flowering-time genes helps pinpoint loci potentially involved in local adaptation.

---

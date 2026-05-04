# RNA-seq Processing & Downstream Analysis Workflow

This document describes the complete bioinformatics workflow used in:

> **Pandey RS, Kotredes KP, Sasner M, Howell GR, Carter GW.** Differential splicing of neuronal genes in a Trem2\*R47H mouse model mimics alterations associated with Alzheimer's disease. *BMC Genomics*. 2023;24:172. [https://doi.org/10.1186/s12864-023-09280-x](https://doi.org/10.1186/s12864-023-09280-x)

The workflow covers upstream RNA-seq read processing through downstream statistical analysis of whole-brain transcriptomes from aging mouse models carrying humanized *APOE4* and/or the *Trem2\*R47H* variant on a C57BL/6J background. It is provided to ensure full reproducibility of published results.

---

## Table of Contents

- [Part 1 — Upstream RNA-seq Processing](#part-1--upstream-rna-seq-processing)
  - [Step 1: Quality Control (FastQC)](#step-1-quality-control)
  - [Step 2: Adapter & Quality Trimming (Trimmomatic)](#step-2-trimming)
  - [Step 3: Alignment (STAR)](#step-3-alignment)
  - [Step 4: Gene-level Quantification (HTSeq)](#step-4-gene-level-quantification)
  - [Step 5: Isoform Quantification (RSEM)](#step-5-isoform-quantification)
- [Part 2 — Downstream Analysis](#part-2--downstream-analysis)
  - [Step 6: Differential Gene Expression (DESeq2)](#step-6-differential-gene-expression)
  - [Step 7: Differential Exon Usage (DEXSeq)](#step-7-differential-exon-usage)
  - [Step 8: Isoform Switch Analysis (IsoformSwitchAnalyzeR)](#step-8-isoform-switch-analysis)
  - [Step 9: Functional Annotation](#step-9-functional-annotation)
  - [Step 10: Cell Type Enrichment](#step-10-cell-type-enrichment)
  - [Step 11: RNA-Binding Protein Site Prediction (RBPmap)](#step-11-rna-binding-protein-site-prediction)
  - [Step 12: Overlap with Human AD Splicing Studies](#step-12-overlap-with-human-ad-splicing-studies)
- [Software Versions](#software-versions)
- [References](#references)

---

## Part 1 — Upstream RNA-seq Processing

All upstream steps were run in a Unix/Linux environment on paired-end FASTQ files. Generic placeholders (`R1.fq`, `R2.fq`, `/path/to/`) should be replaced with actual file paths.

---

### Step 1: Quality Control

**Tool:** FastQC v0.11.3

```bash
fastqc R1.fq R2.fq -o outputdir/
```

FastQC generates per-base quality scores, GC content, adapter content, and duplication rate reports for each FASTQ file. Results should be reviewed before proceeding to trimming.

---

### Step 2: Trimming

**Tool:** Trimmomatic v0.33 | Java v1.7.0_79

```bash
java -jar trimmomatic-0.33.jar PE \
  R1.fq R2.fq \
  R1_paired.fq R1_unpaired.fq \
  R2_paired.fq R2_unpaired.fq \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36
```

**Parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `LEADING` | 3 | Remove leading bases below quality 3 |
| `TRAILING` | 3 | Remove trailing bases below quality 3 |
| `SLIDINGWINDOW` | 4:15 | Trim when 4-base sliding window quality drops below 15 |
| `MINLEN` | 36 | Drop reads shorter than 36 bp after trimming |

---

### Step 3: Alignment

**Tools:** STAR v2.5.3 | Picard v1.95 | SAMtools v1.10

#### 3.1 Build STAR Genome Index

```bash
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /path/to/genome/index/ \
  --genomeFastaFiles /path/to/genome/referenceGenomeFile \
  --sjdbGTFfile /path/to/gtfFile/referencegtfFile \
  --sjdbOverhang 100
```

#### 3.2 Align Reads to Reference Genome

```bash
STAR \
  --genomeDir /path/to/genome/index/ \
  --readFilesIn R1_paired.fq R2_paired.fq \
  --runThreadN 12 \
  --outFileNamePrefix filePrefix \
  --quantMode GeneCounts \
  --sjdbGTFfile /path/to/gtfFILE
```

#### 3.3 Sort SAM Files by Coordinate

```bash
java -Xms1g -Xmx4g -jar picard-v1.95/SortSam.jar \
  INPUT=alignment_sample1.sam \
  OUTPUT=sorted_sample1.sam \
  SO=coordinate
```

#### 3.4 Convert SAM to BAM

```bash
samtools view -bS sorted_sample1.sam -o sample1.out.bam
samtools index sample1.out.bam
```

---

### Step 4: Gene-level Quantification

**Tool:** HTSeq v0.8.0

```bash
htseq-count \
  -s no \
  -f sam \
  sorted_sample1.sam \
  /path/to/gtfFILE \
  > rawcount_sample1.txt
```

The `-s no` flag specifies unstranded library. Output is a per-gene raw count file used as input to DESeq2.

---

### Step 5: Isoform Quantification

**Tool:** RSEM v1.3.3 (with STAR v2.5.3 backend)

#### 5.1 Build RSEM-STAR Index

```bash
rsem-prepare-reference \
  --gtf /path/to/genome/referencegtfFile \
  --star \
  --star-path /path/to/STAR/2.5.3/bin \
  -p 8 \
  /path/to/genome/referenceGenomeFile \
  /path/to/RSEM_star_index/
```

#### 5.2 Alignment and Quantification

```bash
rsem-calculate-expression \
  --paired-end \
  --star \
  --star-path /path/to/STAR/2.5.3/bin \
  -p 8 \
  R1_paired.fq.gz R2_paired.fq.gz \
  /path/to/RSEM_star_index/rsem_index_prefix
```

> **Outputs used downstream:**
> - Raw gene count matrix (HTSeq) → input to **DESeq2** (Step 6)
> - SAM files (STAR) → input to **DEXSeq** (Step 7)
> - Isoform abundance matrix (RSEM) → input to **IsoformSwitchAnalyzeR** (Step 8)

---

## Part 2 — Downstream Analysis

All downstream analyses were performed in R. Steps 6–12 assume that the raw count matrix, isoform abundance matrix, and metadata file are available from Part 1.

---

### Step 6: Differential Gene Expression

**Tool:** DESeq2 v1.16.1

Differentially expressed genes (DEGs) were identified for each mouse model compared to sex and age-matched control mice.

#### 6.1 Load Libraries

```r
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(xlsx)
```

#### 6.2 Load Metadata and Count Data

```r
metadata <- read.csv("/path/to/metadataFile", row.names = 1)
metadata <- metadata[sort(rownames(metadata)), ]

rawdata <- read.csv("/path/to/htseq_countdata_File",
                    check.names = FALSE, row.names = 1)
rawdata <- rawdata[, sort(colnames(rawdata))]
```

#### 6.3 DEG Analysis — Example

The analysis was run for each mouse model × sex × age combination. The example below identifies DEGs in 12-month-old male mice of a given model compared to age and sex-matched controls:

```r
# Subset metadata for the comparison of interest
selected.samples <- metadata[
  metadata$Sex == "Male" &
  metadata$Age == 12 &
  metadata$Genotype %in% c("ModelName", "C57BL/6J"),
]

# Run DEG analysis
DEG(rawdata = rawdata, meta = selected.samples, include.batch = FALSE)

degs_model_M12 <- dseq_res    # significant DEGs only (padj < 0.05)
all_model_M12  <- All_res     # all genes with statistics
```

#### 6.4 Core DEG Function

```r
#' Differential Expression Analysis using DESeq2
#'
#' @param rawdata      Raw count matrix (genes × samples)
#' @param meta         Metadata data frame (rownames = sample IDs)
#' @param include.batch Logical. Whether to include Batch as a covariate.
#' @param ref          Reference genotype level for relevel()
#'
#' Results are assigned to dseq_res (significant DEGs) and All_res (all genes)
#' in the global environment.

DEG <- function(rawdata, meta, include.batch = FALSE, ref = "C57BL/6J") {

  dseq_res <- data.frame()
  All_res  <- data.frame()

  design_formula <- if (include.batch) ~ Batch + Genotype else ~ Genotype

  dat2 <- as.matrix(rawdata[, colnames(rawdata) %in% rownames(meta)])

  ddsHTSeq <- DESeqDataSetFromMatrix(
    countData = dat2,
    colData   = meta,
    design    = design_formula
  )

  ddsHTSeq               <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 10, ]
  ddsHTSeq$Genotype      <- relevel(ddsHTSeq$Genotype, ref = ref)
  dds                    <- DESeq(ddsHTSeq, parallel = TRUE)
  res                    <- results(dds, alpha = 0.05)

  res_sub        <- subset(res[order(res$padj), ], padj < 0.05)
  res_sub$symbol <- map_gene_ids(res_sub, "ENSEMBL", "SYMBOL")
  res$symbol     <- map_gene_ids(res,     "ENSEMBL", "SYMBOL")

  res_sub$EntrezGene <- map_gene_ids(res_sub, "ENSEMBL", "ENTREZID")
  res$EntrezGene     <- map_gene_ids(res,     "ENSEMBL", "ENTREZID")

  dseq_res <<- as.data.frame(res_sub[, c(7:8, 1:6)])
  All_res  <<- as.data.frame(res[,     c(7:8, 1:6)])
}


#' Map Ensembl IDs to Gene Symbols or Entrez IDs
#'
#' @param x          Data frame with Ensembl IDs as rownames
#' @param inputtype  Key type of input (e.g., "ENSEMBL")
#' @param outputtype Key type of output (e.g., "SYMBOL", "ENTREZID")

map_gene_ids <- function(x, inputtype, outputtype) {
  mapIds(
    org.Mm.eg.db,
    keys      = rownames(x),
    column    = outputtype,
    keytype   = inputtype,
    multiVals = "first"
  )
}
```

Step 6.3 was repeated for each mouse model across all sex × age combinations.

---

### Step 7: Differential Exon Usage

**Tool:** DEXSeq v1.40.0

Differential exon usage (DEU) was identified using reads per exon counting bin, following the standard DEXSeq workflow.

#### 7.1 Prepare Annotation (GFF)

```bash
python dexseq_prepare_annotation.py referencegtfFile flatt.mm.gff
```

#### 7.2 Count Reads per Exon Bin

```bash
# Convert BAM to SAM
samtools view /path/to/sample1.out.bam > sample1.sam

# Count exon-level reads
python dexseq_count.py \
  -p yes \
  -r pos \
  flatt.mm.gff \
  sorted_sample1.sam \
  Sample1.txt
```

#### 7.3 Load Libraries and Data in R

```r
library(DEXSeq)
library(BiocParallel)
library(dplyr)
library(tidyverse)

inDir      <- "/path/to/DEXCount/"
countFiles <- list.files(inDir, pattern = ".txt$", full.names = TRUE)
flatFile   <- list.files(inDir, pattern = "gff$",  full.names = TRUE)
```

#### 7.4 Define Sample Table

```r
# Example: Control vs model comparison
sampleTable <- data.frame(
  row.names = c(
    "control1", "control2", "control3",
    "control4", "control5", "control6",
    "model1",   "model2",   "model3",
    "model4",   "model5",   "model6"
  ),
  condition = c(
    rep("control", 6),
    rep("model",   6)
  )
)
```

#### 7.5 Build DEXSeqDataSet

```r
dxd <- DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData   = sampleTable,
  design       = ~ sample + exon + condition:exon,
  flattenedfile = flatFile
)
```

#### 7.6 Run Differential Exon Usage Analysis

```r
dxr <- DEXSeq(dxd, BPPARAM = MulticoreParam(workers = 4))
```

#### 7.7 Filter Significant Results

```r
# Keep exonic regions at FDR < 10%
dxr_significant <- dxr[dxr$padj < 0.1, ]
```

#### 7.8 Save Results

```r
write.xlsx(as.data.frame(dxr_significant),
           file      = "/path/to/Supplementary_Table1.xlsx",
           sheetName = "Model_Age_Sex",
           append    = TRUE)
```

#### 7.9 Overlap Between DEGs and DEU Genes

```r
# Genes with both differential expression and exon usage
overlap_deg_deu <- intersect(degs_model$symbol, dxr_significant$symbol)
# Results saved to Supplementary Table 3, sheet "Overlap_DEGs_DEU"
```

Steps 7.4 to 7.8 were repeated for each mouse model × sex × age comparison.

---

### Step 8: Isoform Switch Analysis

**Tool:** IsoformSwitchAnalyzeR (ISAR) v1.20.0

Isoform switch analysis uses RSEM isoform-level quantification as input. All steps use default parameter settings as described in the [ISAR documentation](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html).

#### 8.1 Load Libraries and Paths

```r
library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)

filedir  <- "/path/to/isoform/countdata/"
metadata <- read.csv("/path/to/metadata.csv",
                     row.names = 1, stringsAsFactors = FALSE)
```

#### 8.2 Select Samples for Analysis

```r
# Example: 12-month-old male model vs control
sub.meta <- metadata[
  metadata$Sex      == "Male" &
  metadata$Age      == 12 &
  metadata$Genotype %in% c("ModelName", "C57BL/6J"),
]

myQuant <- importIsoformExpression(
  sampleVector = paste0(filedir, rownames(sub.meta), ".isoforms.results")
)
```

#### 8.3 Design Matrix

```r
myDesign <- data.frame(
  sampleID  = rownames(sub.meta),
  condition = sub.meta$Genotype
)
```

#### 8.4 Create SwitchAnalyzeRlist Object

```r
mySwitchList <- importRdata(
  isoformCountMatrix   = myQuant$counts,
  isoformRepExpression = myQuant$abundance,
  designMatrix         = myDesign,
  ignoreAfterPeriod    = TRUE,
  isoformExonAnnoation = "/path/to/gtfFILE/mm10.GRCm38.97.gtf.gz",
  isoformNtFasta       = c(
    "/path/to/cdnaFile/Mus_musculus.GRCm38.cdna.all.fa.gz",
    "/path/to/ncrnaFile/Mus_musculus.GRCm38.ncrna.fa.gz"
  ),
  showProgress = TRUE
)
```

#### 8.5 Identify Isoform Switches (Part 1)

```r
mySwitchListPart1 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = mySwitchList,
  dIFcutoff            = 0.05,
  pathToOutput         = "/path/to/output/",
  genomeObject         = Mmusculus,
  outputSequences      = TRUE,
  prepareForWebServers = FALSE
)
```

Two FASTA files are written to the output directory: nucleotide sequences and amino acid sequences. These were submitted to the following external tools at default settings:

| Tool | Purpose |
|------|---------|
| [Pfam](http://pfam.xfam.org) | Protein domain prediction |
| [CPAT](http://rna-tools.online/tools/cpat/) | Coding potential prediction |
| [IUPred2A](https://iupred2a.elte.hu) | Intrinsically disordered region prediction |
| [SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/) | Signal peptide prediction |

#### 8.6 Annotate and Plot Isoform Switches (Part 2)

```r
mySwitchListPart2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist      = mySwitchListPart1,
  dIFcutoff               = 0.05,
  n                       = Inf,
  removeNoncodinORFs      = FALSE,
  consequencesToAnalyze   = c(
    "intron_retention", "coding_potential",
    "ORF_seq_similarity", "NMD_status",
    "domains_identified", "signal_peptide_identified"
  ),
  pathToCPATresultFile    = "/path/to/cpat_result.txt",
  codingCutoff            = 0.721,
  pathToPFAMresultFile    = "/path/to/pfam_results.txt",
  pathToIUPred2AresultFile = "/path/to/IUpred2a_isoform_AA.result",
  pathToSignalPresultFile  = "/path/to/SignalP_output_protein_type.txt",
  pathToOutput            = "/path/to/outputdir/",
  outputPlots             = TRUE
)

# Save results — repeat for each comparison
write.xlsx(mySwitchListPart2$isoformFeatures,
           file      = "/path/to/Supplementary_Table5.xlsx",
           sheetName = "Model_Age_Sex",
           append    = TRUE)
```

#### 8.7 Overlap Between DEGs and Isoform Switch Genes

```r
# Genes with both differential expression and isoform usage
overlap_deg_dtu <- intersect(degs_model$symbol,
                              mySwitchListPart2$isoformFeatures$gene_name)
# Results saved to Supplementary Table 3, sheet "Overlap_DEGs_DTU"
```

Steps 8.2 to 8.6 were repeated for each mouse model × sex × age comparison.

---

### Step 9: Functional Annotation

**Tool:** clusterProfiler (R Bioconductor)

#### 9.1 KEGG Pathway Enrichment — DEU Genes

```r
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Build named list of Entrez IDs per model
geneList <- list(
  Model1_4M  = dxr_model1_4M$EntrezGene,
  Model1_12M = dxr_model1_12M$EntrezGene,
  Model2_4M  = dxr_model2_4M$EntrezGene
  # ... add all comparisons
)

# Run KEGG enrichment across all models simultaneously
kegg_result <- compareCluster(
  geneCluster = geneList,
  fun         = "enrichKEGG",
  organism    = "mmu"
)

# Visualize
dotplot(kegg_result,
        font.size    = 13,
        showCategory = 15,
        title        = "KEGG Pathway Enrichment")

# Save
write.xlsx(as.data.frame(kegg_result),
           file      = "/path/to/Supplementary_Table2.xlsx",
           sheetName = "LOAD_Models",
           append    = TRUE)
```

#### 9.2 GO Biological Process Enrichment — DTU Genes

```r
#' GO Enrichment for Differential Isoform Usage Genes
#'
#' @param x        Character vector of Entrez gene IDs
#' @param universe Background gene universe
#' @param n        Minimum gene count threshold for GO terms

ENRICHGO <- function(x, universe, n = 2) {
  GO.term <- enrichGO(
    gene          = x,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  GO.term@result <- GO.term@result[GO.term@result$Count > n, ]
  GO.term
}

# Example: 12-month-old male mice
dat_male <- results_isar$Model_12M_Male

# Convert Ensembl IDs to Entrez IDs
entrez_male <- as.character(
  as.data.frame(
    mapIds(org.Mm.eg.db,
           keys     = unique(dat_male$isoformFeatures$gene_id),
           column   = "ENTREZID",
           keytype  = "ENSEMBL",
           multiVals = "first")
  )[, 1]
)

go_male <- ENRICHGO(entrez_male, universe, n = 2)

# Save to Supplementary Table 6
write.xlsx(as.data.frame(go_male),
           file      = "/path/to/Supplementary_Table6.xlsx",
           sheetName = "Male")
```

---

### Step 10: Cell Type Enrichment

Brain cell type gene signatures (neurons, endothelial cells, astrocytes, microglia, oligodendrocyte precursor cells, newly formed oligodendrocytes, myelinating oligodendrocytes) were sourced from Zhang et al. \[11\].

Fisher's exact test was used to assess enrichment of cell type signatures in DEGs, DEU genes, and DTU genes:

```r
# Fisher's exact test for cell type enrichment
# A = genes of cell type X in the gene set of interest (e.g., DEU genes in model Y)
# B = genes of cell type X in the database minus A
# C = other genes in the gene set of interest (not cell type X)
# D = total genes in the database minus B and C

pval <- fisher.test(
  matrix(c(A, B, C, D), nrow = 2, ncol = 2),
  alternative = "greater"
)$p.value

OR <- fisher.test(
  matrix(c(A, B, C, D), nrow = 2, ncol = 2),
  alternative = "greater"
)$estimate
```

---

### Step 11: RNA-Binding Protein Site Prediction

**Tool:** RBPmap web server ([rbpmap.technion.ac.il](https://rbpmap.technion.ac.il))

Binding site prediction for RNA-binding proteins (RBM25, HNRNPM, CELF5) was performed using the RBPmap web server at default settings against the mouse genome (GRCm38/mm10 assembly).

**Input:** Genomic coordinates of differentially spliced exonic regions from DEU and DTU analyses (Supplementary Tables 7 and 8).

**Output:** Predicted binding sites were downloaded and appended to Supplementary Tables 7 and 8.

---

### Step 12: Overlap with Human AD Splicing Studies

Hypergeometric tests were used to assess the significance of overlap between differentially spliced genes in mouse models and those reported in human AD splicing studies.

```r
# Compute overlap
Hmgenes <- intersect(Hsgenes, mmgenes)

# Hypergeometric test
pvalue <- phyper(
  q           = length(Hmgenes) - 1,
  m           = length(Hsgenes),
  n           = total_mouse_genes_with_human_orthologs - length(Hsgenes),
  k           = length(mmgenes),
  lower.tail  = FALSE
)
```

**Variable definitions:**

| Variable | Definition |
|----------|-----------|
| `Hsgenes` | All genes from human AD splicing studies with mouse orthologs |
| `mmgenes` | All differentially spliced genes in the given mouse model |
| `Hmgenes` | Overlapping genes between the mouse model and human study |
| `total_mouse_genes_with_human_orthologs` | Background gene universe |

---

## Software Versions

### Upstream Processing

| Tool | Version |
|------|---------|
| FastQC | v0.11.3 |
| Trimmomatic | v0.33 |
| Java | v1.7.0_79 |
| STAR | v2.5.3 |
| Picard | v1.95 |
| SAMtools | v1.10 |
| HTSeq | v0.8.0 |
| RSEM | v1.3.3 |

### Downstream Analysis

| Tool | Version |
|------|---------|
| R | — |
| DESeq2 | v1.16.1 |
| DEXSeq | v1.40.0 |
| IsoformSwitchAnalyzeR | v1.20.0 |
| clusterProfiler | — |
| org.Mm.eg.db | — |
| AnnotationDbi | — |
| BSgenome.Mmusculus.UCSC.mm10 | — |
| BiocParallel | — |

---

## References

1. Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*. 2014;30(15):2114–2120.

2. Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2012;29(1):15–21.

3. Anders S, Pyl PT, Huber W. HTSeq — a Python framework to work with high-throughput sequencing data. *Bioinformatics*. 2015;31(2):166–169.

4. Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. *BMC Bioinformatics*. 2011;12(1):323.

5. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15(12):550.

6. Pandey RS, Kotredes KP, Sasner M, Howell GR, Carter GW. Differential splicing of neuronal genes in a Trem2\*R47H mouse model mimics alterations associated with Alzheimer's disease. *BMC Genomics*. 2023;24:172. https://doi.org/10.1186/s12864-023-09280-x (**This study**)

7. Anders S, Reyes A, Huber W. Detecting differential usage of exons from RNA-seq data. *Genome Research*. 2012;22(10):2008–2017.

8. Carbajosa G, et al. Loss of Trem2 in microglia leads to widespread disruption of cell coexpression networks in mouse brain. *Neurobiology of Aging*. 2018;69:151–166.

9. Vitting-Seerup K, Sandelin A. IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. *Bioinformatics*. 2019;35(21):4469–4471.

10. Yu G, et al. clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*. 2012;16(5):284–287.

11. Zhang Y, et al. An RNA-Sequencing Transcriptome and Splicing Database of Glia, Neurons, and Vascular Cells of the Cerebral Cortex. *Journal of Neuroscience*. 2014;34(36):11929.

12. Paz I, et al. RBPmap: a web server for mapping binding sites of RNA-binding proteins. *Nucleic Acids Research*. 2014;42(W1):W361–W367.

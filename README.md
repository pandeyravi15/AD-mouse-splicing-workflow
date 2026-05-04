# Trem2*R47H Differential Splicing Analysis Workflow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs12864--023--09280--x-blue)](https://doi.org/10.1186/s12864-023-09280-x)
[![BMC Genomics](https://img.shields.io/badge/Published-BMC%20Genomics%202023-green)](https://link.springer.com/article/10.1186/s12864-023-09280-x)

This repository contains the complete bioinformatics workflow used in:

> **Pandey RS, Kotredes KP, Sasner M, Howell GR, Carter GW.**
> Differential splicing of neuronal genes in a Trem2\*R47H mouse model mimics alterations associated with Alzheimer's disease.
> *BMC Genomics*. 2023;24:172.
> [https://doi.org/10.1186/s12864-023-09280-x](https://doi.org/10.1186/s12864-023-09280-x)

---

## Overview

We performed differential gene expression and differential splicing analyses on whole-brain transcriptomes from aging mouse models carrying humanized *APOE4* and/or the *Trem2\*R47H* variant on a C57BL/6J background. This workflow documents every computational step from raw FASTQ processing through functional annotation, enabling full reproduction of published results.

**Key findings:**
- Differentially expressed genes in *Trem2\*R47H* mice were enriched in immune and metabolic pathways
- Differentially spliced genes were enriched in neuronal functions (GABAergic and glutamatergic synapse)
- Significant overlap was observed between spliced genes in *Trem2\*R47H* mice and human AD subjects
- These effects were absent in *APOE4* mice and suppressed in *APOE4·Trem2\*R47H* double mutant mice

---

## Workflow Summary

The full step-by-step workflow is documented in [WORKFLOW.md](WORKFLOW.md).

### Part 1 — Upstream RNA-seq Processing (Unix)

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| 1 | FastQC | v0.11.3 | Read quality assessment |
| 2 | Trimmomatic | v0.33 | Adapter and quality trimming |
| 3 | STAR | v2.5.3 | Alignment to reference genome |
| 4 | HTSeq | v0.8.0 | Gene-level read counting |
| 5 | RSEM | v1.3.3 | Isoform-level quantification |

### Part 2 — Downstream Analysis (R)

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| 6 | DESeq2 | v1.16.1 | Differential gene expression |
| 7 | DEXSeq | v1.40.0 | Differential exon usage |
| 8 | IsoformSwitchAnalyzeR | v1.20.0 | Isoform switch analysis |
| 9 | clusterProfiler | — | KEGG and GO enrichment |
| 10 | — | — | Cell type enrichment (Fisher's exact test) |
| 11 | RBPmap | — | RNA-binding protein site prediction |
| 12 | — | — | Overlap with human AD splicing studies |

---

## Data Availability

Raw sequencing data and processed count matrices are available via the AD Knowledge Portal:

> [https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn66318364](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn66318364)

---

## Repository Structure

```
Trem2-RNAseq-splicing-workflow/
│
├── README.md          # This file
├── WORKFLOW.md        # Complete 12-step documented workflow
├── LICENSE            # MIT license
│
└── scripts/
    ├── DEG_function.R              # DESeq2 wrapper function
    ├── ENRICHGO_function.R         # GO enrichment function
    └── cell_type_enrichment.R      # Fisher's exact test for cell types
```

---

## Requirements

### Unix Environment (Part 1)

| Tool | Version | Installation |
|------|---------|-------------|
| FastQC | v0.11.3 | [bioinformatics.babraham.ac.uk](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| Trimmomatic | v0.33 | [usadellab.org](http://www.usadellab.org/cms/?page=trimmomatic) |
| STAR | v2.5.3 | [github.com/alexdobin/STAR](https://github.com/alexdobin/STAR) |
| Picard | v1.95 | [broadinstitute.github.io/picard](https://broadinstitute.github.io/picard/) |
| SAMtools | v1.10 | [samtools.sourceforge.net](http://samtools.sourceforge.net) |
| HTSeq | v0.8.0 | [htseq.readthedocs.io](https://htseq.readthedocs.io) |
| RSEM | v1.3.3 | [github.com/deweylab/RSEM](https://github.com/deweylab/RSEM) |

### R Environment (Part 2)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "DEXSeq",
  "IsoformSwitchAnalyzeR",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "BiocParallel",
  "BSgenome.Mmusculus.UCSC.mm10"
))

install.packages(c("tidyverse", "xlsx"))
```

---

## Reference Genome

All analyses used the mouse reference genome **GRCm38/mm10** with Ensembl annotation (release 97).

---

## Citation

If you use this workflow, please cite:

> Pandey RS, Kotredes KP, Sasner M, Howell GR, Carter GW. Differential splicing of neuronal genes in a Trem2\*R47H mouse model mimics alterations associated with Alzheimer's disease. *BMC Genomics*. 2023;24:172. https://doi.org/10.1186/s12864-023-09280-x

---

## License

MIT — see [LICENSE](LICENSE) for details.

---

## Contact

**Ravi S. Pandey**
GitHub: [@pandeyravi15](https://github.com/pandeyravi15)

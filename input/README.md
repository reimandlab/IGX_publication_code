# IGX: Cancer genomic alterations and microenvironmental features encode synergistic interactions with disease outcomes

## Detailes and instructions for input files

List of input files/folders:
- [NIHMS948705-supplement-8.xlsx](#1)
- [mc3.v0.2.8.PUBLIC.maf](#2)
- [tcga_cna_gistic/](#3)
- [tcga_code_tables/](#4)
- [TCGA.Kallisto.fullIDs.cibersort.relative.tsv](#5)
- [TCGA-CDR-SupplementalTableS1.xlsx](#6)
- [PanCancerAtlas_subtypes.rds](#7)
- [expr_tpm_20_genes_in_tcga.rds](#8)
- [expr_counts_tcga_brca_lumA.rds](#9)
- [NIHMS958212-supplement-2.xlsx](#10)
- [brca_metabric/](#11)
- [NIHMS45243-supplement-7/table_S43/](#12)
- [LM22.txt](#13)
- [hsapiens.REAC_GOBP.name.gmt](#14)
- [Cosmic_CancerGeneCensus_v102_GRCh37.tsv](#15)
- [OncoKB_cancerGeneList.tsv](#16)

---

### <a id="1"></a>NIHMS948705-supplement-8.xlsx
List of 299 cancer driver genes, downloaded from https://doi.org/10.1016/j.cell.2018.02.060.

---

### <a id="2"></a>mc3.v0.2.8.PUBLIC.maf
Download this file from https://gdc.cancer.gov/about-data/publications/pancanatlas.

---

### <a id="3"></a>tcga_cna_gistic/
Results of GISTIC2 analysis of TCGA tumors. With the current directory being the `input` folder, run the following bash script to generate the `tcga_cna_gistic` folder. The data is downloaded from _FireBrowse_.
```bash
#!/usr/bin/env bash

mkdir tcga_cna_gistic
cd tcga_cna_gistic

base="https://gdac.broadinstitute.org/runs/analyses__latest/data/"
folders=$(curl -s "$base" \
  | grep -oP '(?<=href=")[^"]+/' \
  | sed 's:/$::' \
  | grep '-' |
  | grep -v 'COADREAD' \
  | grep -v 'GBMLGG' \
  | grep -v 'KIPAN' \
  | grep -v 'STES')


for f in $folders; do
  url="${base}/${f}/20160128/gdac.broadinstitute.org_${f}.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz"
  echo "Downloading ${f}..."
  wget -c "$url" || echo "Warning: failed to download $f"
done

for f in *.tar.gz; do
  tar -xvzf "$f" && rm "$f"
done
```

<hr style="border: 0; height: 1px; background: #eee;" />

### <a id="4"></a>tcga_code_tables/
TCGA code tables downloaded from https://gdc.cancer.gov/system/files/public/file/tcga_code_tables.zip

---

### <a id="5"></a>TCGA.Kallisto.fullIDs.cibersort.relative.tsv
Results of CIBERSORT analysis of TCGA tumors, downloaded from https://gdc.cancer.gov/about-data/publications/panimmune.

---

### <a id="6"></a>TCGA-CDR-SupplementalTableS1.xlsx
TCGA clinical information of patients, downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas.

---

### <a id="7"></a>PanCancerAtlas_subtypes.rds
Cancer subtype annotation of TCGA tumors. With the current directory being the `input` folder, run the following R code to generate this file.
```R
library(TCGAbiolinks)
subtypes <- PanCancerAtlas_subtypes()
saveRDS(subtypes, "PanCancerAtlas_subtypes.rds")
```

---

### <a id="8"></a>expr_tpm_20_genes_in_tcga.rds
Follow the steps below to generate this file:
1) download the mRNA expression data from TCGAbiolinks with the following arguments: {project = "_TCGA-***_", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts"}. The "_TCGA-***_" should be set separately to each TCGA study used in the anlaysis.
2) extract the "TPM" data and turn it into a dataframe with columns {"gene_name", "gene_type", "TCGA barcodes..."}.
3) Filter to only include these 20 genes: {_CTLA4_, _LAG3_, _CD274_, _PDCD1_, _GZMA_, _PRF1_, _MEN1_, _JUN_, _ARID1A_, _VHL_, _FBXW7_, _CASP3_, _NPM1_, _EGFR_, _MYC_, _COX6C_, _HIF1A_, _TP53_, _ERBB2_, _AXIN2_}.
4) Filter TCGA barcodes to primary tumor samples.
5) Convert TCGA barcodes to patient-id, i.e. _TCGA-TSS-Participant_.
6) Remove patient-ids with duplicated samples (if there is any).
7) Filter patient-ids to those included in the IGX analysis (refer to the file "intermediate/tcga_all_datasets.rds").
8) Save this dataframe as an RDS object in the `input` folder.

---

### <a id="9"></a>expr_counts_tcga_brca_lumA.rds
Follow the steps below to generate this file:
1) download the mRNA expression data from TCGAbiolinks with the following arguments: {project = "TCGA-BRCA", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts" }.
2) extract the "raw counts" and turn the data into a dataframe with columns {"gene_name", "gene_type", "TCGA barcodes..."}.
3) Filter to protein-coding genes.
4) Filter TCGA barcodes to primary tumor samples
5) Convert TCGA barcodes to patient-id, i.e. _TCGA-TSS-Participant_.
6) Remove patient-ids with duplicated samples (if there is any).
7) Filter patient-ids to those included in the IGX analysis (refer to the file "intermediate/tcga_all_datasets.rds").
8) Save this dataframe as an RDS object in the `input` folder.

---

### <a id="10"></a>NIHMS958212-supplement-2.xlsx
Collection of Immuno-genomic measures for TCGA tumors downloaded from https://doi.org/10.1016/j.immuni.2018.03.023.

---

### <a id="11"></a>brca_metabric/
The folder containing the METABRIC dataset, including clinical information, copy number alterations and transcriptomic data, downloaded from cBioPortal: https://www.cbioportal.org/study/summary?id=brca_metabric.

---

### <a id="12"></a>NIHMS45243-supplement-7/table_S43/
The list of patient IDs in _discovery_ and _validation_ cohorts of METABRIC were obtained from the two tables in this folder. Downloaded from https://doi.org/10.1038/nature10983.

---

### <a id="13"></a>LM22.txt
The signature matrix for distinguishing 22 human hematopoietic cell subsets in bulk tissues, downloaded from https://doi.org/10.1038/nmeth.3337.

---

### <a id="14"></a>hsapiens.REAC_GOBP.name.gmt
GMT file containing GO:BP (biological process) and Reactome terms, downloaded from _g:Profiler_ on August 5th, 2025.

---

### <a id="15"></a>Cosmic_CancerGeneCensus_v102_GRCh37.tsv
Cancer Gene Census v102 GRCh37, downloaded from https://cancer.sanger.ac.uk/cosmic/download/cosmic/v102/cancergenecensus.

---

### <a id="16"></a>OncoKB_cancerGeneList.tsv
OncoKB cancer gene list, downloaded from https://www.oncokb.org/cancer-genes, on August 5th, 2025.

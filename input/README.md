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

---

### <a id="1"></a>NIHMS948705-supplement-8.xlsx
List of 299 cancer driver genes, downloaded from https://doi.org/10.1016/j.cell.2018.02.060.

---

### <a id="2"></a>mc3.v0.2.8.PUBLIC.maf
Download this file from https://gdc.cancer.gov/about-data/publications/pancanatlas.

---

### <a id="3"></a>tcga_cna_gistic/
Results of GISTIC2 analysis of TCGA tumors. With the current directory being the `input` folder, run the following bash script to generate the `tcga_cna_gistic` folder. The data is downloaded from FireBrowse.
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
TCGA code tables downloaded from .

---

### <a id="5"></a>TCGA.Kallisto.fullIDs.cibersort.relative.tsv
Results of CIBERSORT analysis of TCGA tumors, downloaded from .

---

### <a id="6"></a>TCGA-CDR-SupplementalTableS1.xlsx
TCGA clinical information of patients, downloaded from .

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
1) download the mRNA expression data from TCGAbiolinks with the following arguments: {project = "TCGA-BRCA", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts"}.
2) extract the "TPM" data and turn it into a dataframe with columns {"gene_name", "gene_type", "TCGA barcodes..."}.
3) Filter to only include these 20 genes: {_CTLA4_, _LAG3_, _CD274_, _PDCD1_, _GZMA_, _PRF1_, _MEN1_, _JUN_, _ARID1A_, _VHL_, _FBXW7_, _CASP3_, _NPM1_, _EGFR_, _MYC_, _COX6C_, _HIF1A_, _TP53_, _ERBB2_, _AXIN2_}.
4) Filter TCGA barcodes to primary tumor samples.
5) Convert TCGA barcodes to patient-id, i.e. _TCGA-TSS-Participant_.
6) Remove patient-ids with duplicated samples (if there is any).
7) Filter patient-ids to those included in the IGX analysis (refer to the file "intermediate/tcga_all_datasets.rds").
8) Save this dataframe as an RDS object in the `input` folder.










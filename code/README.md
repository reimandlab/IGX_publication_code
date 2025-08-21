# IGX: Cancer genomic alterations and microenvironmental features encode synergistic interactions with disease outcomes

## Source code scripts
The scripts run in numeric order of the filenames. The intermediate files are stored in the `intermediate` folder. The figure panels and the supplementary table are stored in the `output` folder.

**Upstream scripts**
- 001: Compile the mutation features of TCGA tumors.
- 002: Compile the copy number alteration (CNA) features of TCGA tumors.
- 003: Compile the immune features of TCGA tumors.
- 004: Compile the clinical features of TCGA tumors.
- 005: Create the object holding all features of tumors in each TCGA cohort (cancer type) separately.
- 006: Run PACIFIC-step1. Accumulate the required number of feature selection iterations in the designated intermediate folder.
- 007: Run PACIFIC-step2. Aggregate the feature selection iterations and store the _candidate IGXs_ for each cohort.

**Downstream scripts**
- Catalog of 34 IGXs:
  - 101: Create the final _catalog of IGXs_ by combining the PACIFIC results across the cohorts and applying FDR correction.
  - 102: Create the supplementary table containing general information about the IGXs in the catalog.
  - 103: Implement functions used for making Kaplan-Meier (KM) plots.
  - 104: Make KM plots for all IGXs in the catalog (used in panels 2E,F,G, 5A, S2, S6).
- Figure 2:
  - 201: Grid plot of the catalog (panel 2A).
  - 202: Plot P-value and logHR values of IGXs vs genomic or immune features (panels 2B,C).
  - 203: Plot fraction of cancer samples in each IGX (panel 2D).
  - 204: Plot expression of CNA-affected genes for the CNA features involved in IGXs (panel 2H).
- Figure 3:
  - 301: Barplots of logHRs per I or G featrures in each IGX (panel 3A).
  - 302: Plot the average values of immune features per each IGX group (panel 3B).
- Figure 4:
  - 401: Test for and plot the associations of IGXs with additional genomic and immunogenic characteristics (panels 4A,B,C,D, 5H).
- Figure 5:
  - 501: Compile the clinical features of METABRIC tumors.
  - 502: Compile the CNA features of METABRIC tumors.
  - 503: Prepare the gene expression profiles of METABRIC tumors used for CIBERSORT analysis.
  - 504: Compile the immune features of METABRIC tumors.
  - 505: Create the object holding all features of tumors in METABRIC luminal A and luminal B subtypes separately.
  - 506: Run differenrial expression analysis for the "luminal A IGX" applied to METABRIC-LumA data.
  - 507: Run differenrial expression analysis for the "luminal A IGX" applied to TCGA-LumA data.
  - 508: Make a combined object of the two differential expression analyses above.
  - 509: Run pathway enrichment analysis (files used for panel 5F stored in `output/enrichment_map/`).
  - 510: Make KM plot for the luminal A IGX in METABRIC-LumA and for the luminal B IGX in METABRIC-LumB (panels 5B, S6).
  - 511: Plot _MEN1_ expression with respect to 11q13.1 loss in the two datasets (panel 5C).
  - 512: Plot the Neutrophil levels in the two datasets with respect to the luminal A IGX (panel 5D).
  - 513: Plot the differentially expressed genes with respect to the luminal A IGX in the two datasets (panel 5E).
  - 514: Plot the genes contributing to the enrichment of selected pathways (panel 5G).
- Remaining supplementary figures:
  - 601: plot expression of _EGFR_ with respect to the IGX found in TCGA-CESC (panel S4).
  - 602: Plot genes down-regulated by 11q13.1 deletion in the two luminal A breast cancer datasets (panel S5).
  - 602: Plot the distributions of immune cell levels from CIBERSORT analyses, and the median threshold used for defining the binary immune features in cohorts of individual cancer types (panel S3).
 


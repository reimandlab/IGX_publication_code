# IGX: Cancer genomic alterations and microenvironmental features encode synergistic interactions with disease outcomes

## Source code scripts
The scripts run in numeric order of the filenames. The intermediate files are stored in the `intermediate` folder. The figure panels and the supplementary table are stored in the `output` folder.

**_Upstream scripts_**
- 001: Compile the mutation features of TCGA tumors.
- 002: Compile the copy number alteration (CNA) features of TCGA tumors.
- 003: Compile the immune features of TCGA tumors.
- 004: Compile the clinical features of TCGA tumors.
- 005: Create the object holding all features of tumors for each cohort (cancer type) separately.
- 006: Run PACIFIC-step1. Accumulates the required number of feature selection iterations in the designated intermediate folder.
- 007: Run PACIFIC-step2. Aggregates the feature selection iterations and stores the _candidate IGXs_ for each cohort.

**_Downstream scripts_**
- Catalog of 34 IGXs:
  - 101: Create the final _catalog of IGXs_ by combining the PACIFIC results across the cohorts and applying FDR correction. Store the catalog in a specially structured object.
  - 102: Create the supplementary table containing general information about the IGXs in the catalog.
  - 103: Implement functions used for making Kaplan-Meier (KM) plots.
  - 104: Make KM plots for all IGXs in the catalog.
- Figure 2:
  - 201: Grid plot of the catalog (panel A).
  - 202: Plot P-value and logHR values of IGXs vs G or I features (panels B & C).

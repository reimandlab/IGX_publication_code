# IGX: Cancer genomic alterations and microenvironmental features encode synergistic interactions with disease outcomes

## Source code and input files

This code repository contains the source code and the majority of input and intermediate files to generate results and figures of the manuscript: 

"Cancer genomic alterations and microenvironmental features encode synergistic interactions with disease outcomes" by Bayati et al (2025). 

To understand the structure and flow of the pipeline, see the README in the `code` folder. For instructions and details about the input files—especially those that are initially missing and must be downloaded from publicly available resources—refer to the README in the `input` folder.

### Environment Setup

The source code is in R (version 4.3.3). The following R package versions were used:

#### For upstream scripts \^
``` bash
data.table v1.15.4 
maftools v2.18.0 
survival v3.7.0 
readxl v1.4.3
PACIFIC v1.0.0 (install from https://github.com/reimandlab/PACIFIC)
```
#### For downstream scripts \^
``` bash
data.table v1.15.4
survival v3.7.0
ggplot2 v3.5.1 
ggrepel v0.9.6
patchwork v1.3.0
stringr v1.5.1
coin v1.4.3
limma v3.58.1
DESeq2 v1.42.1
ActivePathways v2.0.5
ggplotify v0.1.2
egg v0.4.5
```
\^ See the README in the `code` folder.

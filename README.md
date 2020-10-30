# Repository with workflow for competitive read mapping and allele-specific expression counts


## Description of the repository

In the main directory you will find `bash` scripts named from 1 to 4. They should be run in that order (see [Requirements](#Requirements)).

### Script 1

It will download the raw reads directly from SRA with `fasterq-dump` using the accession numbers stored in `metadata-7938743-processed-ok.tsv` and rename them. Next step is to filter and trim reads with `trimmomatic`.

### Script 2

Generates genome indices for `star` based on files stored in [references](https://github.com/santiagosnchez/competitive_mapping_workflow/tree/master/references).

### Script 3

Maps reads using `star` and generates BAM files as well as their indices.

### Script 4

Counts reads for parent species alignments using `featureCounts` and for hybrid data and allele-specific expression it uses `CompMap`. It will also generate tables with raw counts for each sample.

### Simulations

This directory includes `bash` and `r` scripts running RNA-seq data simulations for allele-specific expression using the package `polyester`.

### Counts

All the raw count tables are included here.

### Analyses

This contains two subdirectories [scripts](https://github.com/santiagosnchez/competitive_mapping_workflow/tree/master/analyses/scripts) and [tables](https://github.com/santiagosnchez/competitive_mapping_workflow/tree/master/analyses/tables).

`scripts` has further subdirectories with `r` scripts to generate summary tables from differential expression (and other) analyses.
`tables` has all the raw summary tables used for this work.

Running some of the `r` scripts will overwrite the files already in that directory. So do it carefully!

## Requirements

### Unix command-line programs

* `samtools`
* [`CompMap`](https://github.com/santiagosnchez/CompMap)
* `GNU-parallel`
* `star`
* `BWA`
* `trimmomatic`
* `sratools`
* `featureCounts`

### R packages

* `polyester`
* `glimma`
* `edgeR`
* `DESeq2`
* `MASS`
* `RColorBrewer`
* `cowplot`
* `dplyr`
* `edgeR`
* `ggVennDiagram`
* `ggforce`
* `ggplot2`
* `ggplotify`
* `ggpointdensity`
* `ggrepel`
* `ggridges`
* `ggstance`
* `ggtext`
* `grid`
* `gridExtra`
* `gtable`
* `lemon`
* `parallel`
* `scales`
* `tibble`
* `tidyr`
* `venneuler`
* `viridis`

### Kallisto Alignment and Gene Annotation

Raw reads are aligned to the *S. papillosus* reference transcriptome
(PRJEB525.WBPS14.mRNA\_transcripts, downloaded from [WormBase
ParaSite](https://parasite.wormbase.org/Strongyloides_papillosus_prjeb525/Info/Index)
on 17 August 2020), using Kallisto. Kallisto alignments are imported
into the R environment using `Tximport`. Counts are generated from
abundance files using the `lengthScaledTPM` option.  

Count data is annotated with information imported via the WormBase
ParaSite BioMaRT. Annotation information includes:

-   *C. elegans* homologs/percent homology
-   UniProtKB number
-   Interpro terms
-   GO terms
-   general Description information

[Hunt *et al* 2016](https://www.nature.com/articles/ng.3495) establishes
two distinct subclades from the four sequenced *Strongyloides* species:
*S. venezuelensis-S. papillosus* and *S. ratti-S. stercoralis*. Thus, we
also include annotation information for the appropriate in-group and both members of the out-group, here:

-   In-group: *S. venezuelensis* homologs/percent homology
-   Out-group: *S. stercoralis* and *S. ratti* homologs/percent homology

### Filtering and Normalization Steps

Raw reads were quantified as log2 counts per million (CPM) using the `EdgeR` package,
then filtered to remove transcripts with low counts (less than 1 log2CPM in at least 1 sample). Non-discarded gene values are
normalized using the trimmed mean of M-values method [(TMM, Robinson and
Oshlack)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)
to permit between-samples comparisons. The mean-variance relationship
was modeled using a precision weights approach [(Law *et al*
2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29).  

These steps produced a variance-stabilized digital gene expression list (vDGEList), which serves as the primary input to the Shiny App.

### Saved Output Files

During pre-processing, the following required data files were saved and
are imported into Shiny application upon species initialization:

1.  a gene annotation R object (`Sp_geneAnnotations`)
2.  the variance-stabilized vDGEList, saved as an R object
    (`Sp_vDGEList`)
3.  a matrix of discarded genes and their raw counts
    (`SpRNAseq_discardedGene_counts.csv`) - this data is downloadable
    from within the Browser App
4.  a matrix of variance-stabilized gene expression data, extracted from
    the vDGEList (`SpRNAseq_log2cpm_filtered_norm_voom.csv`) - this data
    is downloadable from within the Browser App

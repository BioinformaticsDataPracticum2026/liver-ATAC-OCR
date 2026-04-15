# Motif Enrichment Analysis of Liver ATAC-seq Data of Mouse and Human Tissues

## Aim:

The aim of this task is to identify transcription factor binding motifs (TF motifs) that are enriched in the conserved and species-specific open chromatin regions identified in the mouse and human liver ATAC-seq data. This analysis will provide insights into the transcription factors that may be driving regulatory activity in the liver and how this may differ between human and mouse.

## Input Data:

The input data for this analysis consists of the conserved and species-specific open chromatin regions (OCRs) identified in the mouse and human liver ATAC-seq data. These OCRs were classified into candidate promoters and candidate enhancers based on their overlap with annotated transcription start sites (TSS) using BEDTools. The OCRs were further categorized into conserved (open in both species) and species-specific (open in one species but closed in the other) groups. All the input files for Motif Analysis are present in the [PE_classification/output](liver-ATAC-OCR/tree/main/PE_classification/output) directory.

## Methodology:


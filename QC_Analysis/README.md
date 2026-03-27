# Quality Control Analysis of ATAC-seq Data of Mouse and Human Tissues

## Aim:

The aim of this analysis is to evaluate the quality of ATAC-seq data obtained from mouse and human liver and adrenal gland tissues by assessing various quality control (QC) metrics. This will help determine the reliability and reproducibility of the data for downstream analyses. Out of the 2 organs (liver and adrenal gland [alternative] tissues) among 2 species (mouse and human),the organ with best quality data will be identified based on the QC metrics for both organisms.

## Metrics Identified for QC Analysis:

1. **% Mapped Reads** - Should be higher (~90% or more) for good quality data.
2. **% Properly Paired Reads** - Should be higher (~90% or more) for good quality data.
3. **Periodicity plots** - Should show three clear, distinct humps for good quality data, indicating nucleosome structure.
4. **TSS-Enrichment Score** - Should be higher (>=10) for good quality data.
5. **NRF (Non-Redundant Fraction)** - Should be higher (>0.8) for good quality data.
6. **IDR (Irreproducible Discovery Rate) Plots** - Should contain peaks that are reproducible across replicates for good quality data. Reproducible peaks = Strong peaks that are present in both replicates, and not just in one replicate.
7. **Rescue Ratio** - Should be lower (<2) for good quality data.
8. **Self-consistency Ratio** - Should be lower (<2) for good quality data.

## QC Analysis Results:

| Metric | Mouse Liver Replicate 1 | Mouse Liver Replicate 2 | Human Liver Replicate 1 | Human Liver Replicate 2 | Mouse Adrenal Gland Replicate 1 | Mouse Adrenal Gland Replicate 2 | Human Adrenal Gland Replicate 1 | Human Adrenal Gland Replicate 2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| % Mapped Reads | $\color{green}{97.8\%}$ | 97.89999999999999% | 98.9% | 98.6% | 98.5% | 98.6% | 98.5% | 98.9% |	
| % Properly Paired Reads | 95.7% | 95.6% | 97.89999999999999% | 97.5% | ９4.6% | ９4.6９９９９９９９９９９９９% | ９７.０% | ９７.５% |
| Periodicity Plots | 3 Clear, distinct humps | 3 Clear, distinct humps | 3 Clear, distinct humps | 3 Clear, distinct humps | 2 Clear, distinct humps | 2 Clear, distinct humps | 2 Clear and distinct humps | 2 Clear and distinct humps |
| TSS-Enrichment Score | 7.692382151686117 | 7.354598361825624 | 23.294805448649445 | 20.962438980086663 | 18.365117628803656 | 18.997852749462098 | 26.084253634636482 | 14.011771178374902 |
| NRF (Non-Redundant Fraction) | 0.935887 | 0.94404 | 0.884593 | 0.899531 | 0.301612 | 0.270462 | 0.778014 | 0.969527 |
| Rescue Ratio | 1.0043626038838573 | 1.0122929315643505 | 1.0906310179587084 | 1.24175568252391 | 1.013263779368353 | 1.0098552058921995 | 2.0842882741797046 | 4.283663165487085|
| Self-consistency Ratio | 1.0756982081825852 | 1.1043673083661616 | 1.2192105681999046 | 1.4284339019600292 | 1.2220589050973432 | 1.2891362811933205 | 3.6038942470169237 | 6.26995057660626 |
| Reproducibility Test | Pass | Pass | Pass | Pass | Pass | Pass | Fail | Fail |

## Conclusion:

Based on the QC metrics evaluated, the mouse liver and human liver samples show good quality data with high % mapped reads, high % properly paired reads, clear periodicity plots, high TSS-enrichment scores, high NRF values, low rescue ratios, and low self-consistency ratios. Both replicates for mouse liver and human liver passed the reproducibility test.

In contrast, the mouse adrenal gland samples show moderate quality data with slightly lower % mapped reads and % properly paired reads, less clear periodicity plots, moderate TSS-enrichment scores, low NRF values, and higher rescue ratios and self-consistency ratios. Both replicates for mouse adrenal gland passed the reproducibility test.

The human adrenal gland samples show poor quality data with high % mapped reads and % properly paired reads, but less clear periodicity plots, moderate TSS-enrichment scores, variable NRF values, and very high rescue ratios and self-consistency ratios. Both replicates for human adrenal gland failed the reproducibility test.

Overall, the mouse liver and human liver samples have the best quality data among the four tissue types analyzed, while the human adrenal gland samples have the poorest quality data.
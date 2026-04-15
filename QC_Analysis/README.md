# Quality Control Analysis of ATAC-seq Data of Mouse and Human Tissues

## Aim:

The aim of this task is to evaluate the quality of ATAC-seq data obtained from mouse and human liver and adrenal gland tissues by assessing various quality control (QC) metrics. This will help determine the reliability and reproducibility of the data for downstream analyses. Out of the 2 organs (liver and adrenal gland [alternative] tissues) among 2 species (mouse and human),the organ with best quality data will be identified based on the QC metrics for both organisms.

## Metrics Identified for QC Analysis:

1. **% Mapped Reads** - Should be higher (~90% or more) for good quality data.
2. **% Properly Paired Reads** - Should be higher (~90% or more) for good quality data.
3. **Periodicity plots** - Should show three clear, distinct humps for good quality data, indicating nucleosome structure.
4. **TSS-Enrichment Score** - Should be higher (>=10) for good quality data.
5. **NRF (Non-Redundant Fraction)** - Should be higher (>0.8) for good quality data.
6. **IDR (Irreproducible Discovery Rate)** - Should contain peaks that are reproducible across replicates for good quality data. Reproducible peaks = Strong peaks that are present in both replicates, and not just in one replicate. The reproducibility is measured by;
    1. **Rescue Ratio** - Should be lower (<2) for good quality data.
    2. **Self-consistency Ratio** - Should be lower (<2) for good quality data.

## QC Analysis Results:

| Metric | Mouse Liver Replicate 1 | Mouse Liver Replicate 2 | Human Liver Replicate 1 | Human Liver Replicate 2 | Mouse Adrenal Gland Replicate 1 | Mouse Adrenal Gland Replicate 2 | Human Adrenal Gland Replicate 1 | Human Adrenal Gland Replicate 2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| % Mapped Reads (%) | $\color{green}97.8$ | $\color{green}97.89999999999999$ | $\color{green}98.9$ | $\color{green}98.6$ | $\color{green}98.5$ | $\color{green}98.6$ | $\color{green}98.5$ | $\color{green}98.9$ |	
| % Properly Paired Reads (%) | $\color{green}95.7$ | $\color{green}95.6$ | $\color{green}97.89999999999999$ | $\color{green}97.5$ | $\color{green}94.6$ | $\color{green}94.6999999999999$ | $\color{green}97.0$ | $\color{green}97.5$ |
| Periodicity Plots | $\color{green}3 \text{ Clear, distinct}$<br>$\color{green} \text{humps}$ | $\color{green}3 \text{ Clear, distinct}$<br>$\color{green} \text{humps}$ | $\color{green}3 \text{ Clear, distinct}$<br>$\color{green} \text{humps}$ | $\color{green}3 \text{ Clear, distinct}$<br>$\color{green} \text{humps}$ | $\color{yellow}2 \text{ Clear, distinct}$<br>$\color{yellow} \text{humps}$ | $\color{yellow}2 \text{ Clear, distinct}$<br>$\color{yellow} \text{humps}$ | $\color{yellow}2 \text{ Clear and distinct}$<br>$\color{yellow} \text{humps}$ | $\color{yellow}2 \text{ Clear and distinct}$<br>$\color{yellow} \text{humps}$ |
| TSS-Enrichment Score | $\color{yellow}7.692382151686117$ | $\color{yellow}7.354598361825624$ | $\color{green}23.294805448649445$ | $\color{green}20.962438980086663$ | $\color{green}18.365117628803656$ | $\color{green}18.997852749462098$ | $\color{green}26.084253634636482$ | $\color{green}14.011771178374902$ |
| NRF (Non-Redundant Fraction) | $\color{green}0.935887$ | $\color{green}0.94404$ | $\color{green}0.884593$ | $\color{green}0.899531$ | $\color{red}0.301612$ | $\color{red}0.270462$ | $\color{yellow}0.778014$ | $\color{green}0.969527$ |
| Rescue Ratio | $\color{green}1.0043626038838573$ | $\color{green}1.0122929315643505$ | $\color{green}1.0906310179587084$ | $\color{green}1.24175568252391$ | $\color{green}1.013263779368353$ | $\color{green}1.0098552058921995$ | $\color{red}2.0842882741797046$ | $\color{red}4.283663165487085$ |
| Self-consistency Ratio | $\color{green}1.0756982081825852$ | $\color{green}1.1043673083661616$ | $\color{green}1.2192105681999046$ | $\color{green}1.4284339019600292$ | $\color{green}1.2220589050973432$ | $\color{green}1.2891362811933205$ | $\color{red}3.6038942470169237$ | $\color{red}6.26995057660626$ |
| Reproducibility Test | $\color{green}\text{Pass}$ | $\color{green}\text{Pass}$ | $\color{green}\text{Pass}$ | $\color{green}\text{Pass}$ | $\color{green}\text{Pass}$ | $\color{green}\text{Pass}$ | $\color{red}\text{Fail}$ | $\color{red}\text{Fail}$ |

## Conclusion:

Based on the QC metrics evaluated, the mouse liver and human liver samples show good quality data with high % mapped reads, high % properly paired reads, clear periodicity plots, high TSS-enrichment scores, high NRF values, low rescue ratios, and low self-consistency ratios. Both replicates for mouse liver and human liver passed the reproducibility test.

In contrast, the mouse adrenal gland samples show moderate quality data with slightly lower % mapped reads and % properly paired reads, less clear periodicity plots, moderate TSS-enrichment scores, low NRF values, and higher rescue ratios and self-consistency ratios. Both replicates for mouse adrenal gland passed the reproducibility test.

The human adrenal gland samples show poor quality data with high % mapped reads and % properly paired reads, but less clear periodicity plots, moderate TSS-enrichment scores, variable NRF values, and very high rescue ratios and self-consistency ratios. Both replicates for human adrenal gland failed the reproducibility test.

Overall, the mouse liver and human liver samples have the best quality data among the four tissue types analyzed, while the human adrenal gland samples have the poorest quality data.
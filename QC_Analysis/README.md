# Quality Control Analysis of ATAC-seq Data of Mouse and Human Tissues

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
| % Mapped Reads | 95% | 94% | 92% |  | 91% | 96% | 95% | 93% | 92% |
| % Properly Paired Reads | 93% | 92% | 90 % |  | 89% | 94% | 92% | 90% | 91% |
| Periodicity Plots | Clear, distinct humps | Clear, distinct humps | Clear, distinct humps |  | Clear, distinct humps | Clear, distinct humps | Clear, distinct humps | Clear, distinct humps | Clear, distinct humps |
| TSS-Enrichment Score | 12 | 11 | 10 |  | 9 | 13 | 11 | 10 | 9 |
| NRF (Non-Redundant Fraction) | 0.85 | 0.83 | 0.80 |  | 0.78 | 0.88 | 0.84 | 0.81 | 0.79 |
| IDR Plots | Reproducible peaks present | Reproducible peaks present | Reproducible peaks present |  | Reproducible peaks present | Reproducible peaks present | Reproducible peaks present | Reproducible peaks present | Reproducible peaks present |
| Rescue Ratio | 1.5 | 1.4 | | 1.6 |  | 1.7 | 1.3 | 1.5 | 1.4 |
| Self-consistency Ratio | 1.8 | | 1.7 | 1.6 |  | 1.5 | 1.4 | 1.6 | 1.5 |

## Conclusion:


# ***T***ranscription ***R***egulatory ***A***nalysis of ***C***onserved ***E***lements (***TRACE***)

![Python](https://img.shields.io/badge/python-3.10+-blue.svg)
![NumPy](https://img.shields.io/badge/numpy-1.21+-blue.svg)
![Matplotlib](https://img.shields.io/badge/matplotlib-3.4+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

## Aim:
The project involves comparing promoters and enhancers between Human and Mouse data for a specific tissue.

## Tissues Used:
*Adrenal Gland [Alternative]*

***Liver***

## Tasks:
1. **Evaluation of data quality for human and mouse open chromatin data**
2. **Mapping and comparing open chromatin regions across all species**
3. **Identifying biological processes likely to be regulated by open chromatin regions**
4. **Comparing candidate enhancers to candidate promoters**
5. **Find transcription factors that tend to bind open chromatin regions**

### 1. Quality Control Analysis
The [QC_Analysis](QC_Analysis/README.md) provides a comprehensive evaluation of the quality of the human and mouse open chromatin data. It includes assessments of read quality, mapping efficiency, and peak calling performance. The analysis ensures that the mouse and human liver ATAC-seq data are of high quality and suitable for downstream analyses.

### 2. Cross-Species mapping using HALPER
HALPER is used to identify orthologous open chromatin regions between the human and mouse genome. The pipeline takes ATAC-seq peak files from each species as input and performs liftover of the genomes across species using the HAL alignment format. HALPER then post-processes the liftover results to identify high confidence regulatory regions. The mapping enables comparison of conserved and species-specific regions for downstream analysis.

### 3. GO enrichment analysis using rGREAT
Regions are grouped into conserved (open in both species) and species-specific (open in one species but closed in the other). rGREAT is run on each group to identify enriched GO Biological Processes. It assigns regions to genes and performs enrichment based on genome coverage.

### 4. Promoter and enhancer classification using BEDTools
OCRs are classified into candidate promoters and candidate enhancers using BEDTools. Promoter regions are generated from annotated transcription start sites, and OCRs are grouped based on whether they overlap promoter regions or fall outside them. This analysis helps compare promoter-like and enhancer-like regulatory regions across species and conservation categories.

### 5. Motif enrichment analysis using MEME-ChIP
MEME-ChIP is used to identify transcription factor binding motifs that are enriched in the conserved and species-specific open chromatin regions. The analysis provides insights into the transcription factors that may be driving regulatory activity in the liver and how this may differ between human and mouse.

## Dependencies:
Python version *3.10 or newer* (https://www.python.org/downloads/release/python-3100/)

R version *4.0 or newer* (https://cran.r-project.org/bin/macosx/base/)

**Python modules**
- numpy (https://numpy.org/)
- matplotlib (https://pypi.org/project/matplotlib/)

**R modules**
- Bioconductor (https://bioconductor.org/install/)

## Citations:
1. ***Liu et al. An ATAC-seq atlas of chromatin accessibility in mouse tissues. Sci Data 6, 65 (2019). (https://doi.org/10.1038/s41597-019-0071-0)***

2. ***Currin et al. Genetic effects on liver chromatin accessibility identify disease regulatory variants. The American Journal of Human Genetics, 2021; 108, 1169-1189. (https://doi.org/10.1016/j.ajhg.2021.05.001)***

3. ***ENCODE portal (https://www.encodeproject.org/)***

4. ***HALPER (https://github.com/pfenninglab/halLiftover-postprocessing)***

5. ***rGREAT (https://github.com/jokergoo/rGREAT)***

6. ***BEDTools (https://bedtools.readthedocs.io/en/latest/)***

7. ***MEME-ChIP (https://meme-suite.org/meme/doc/meme-chip.html)***

8. ***Philip Machanick and Timothy L. Bailey, "MEME-ChIP: motif analysis of large DNA datasets", Bioinformatics 27(12):1696-1697, 2011. (https://doi.org/10.1093/bioinformatics/btr189)***

9. ***Timothy L. Bailey and Charles Elkan, "Fitting a mixture model by expectation maximization to discover motifs in biopolymers", Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology, pp. 28-36, AAAI Press, Menlo Park, California, 1994. (http://www.aaai.org/Papers/ISMB/1994/ISMB94-004.pdf)***

10. ***Timothy L. Bailey, "STREME: accurate and versatile sequence motif discovery", Bioinformatics, Mar. 24, 2021. (https://doi.org/10.1093/bioinformatics/btab203)***

## Contributors:

***Evan Lin (evanlin@andrew.cmu.edu)***

***Arun Sujatha Bharath Raj (asujatha@andrew.cmu.edu)***

***Nikita Rajesh (nrajesh@andrew.cmu.edu)***

***Suratha Sriram (surathas@andrew.cmu.edu)***

## License:
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

# ***T***ranscription ***R***egulatory ***A***nalysis of ***C***onserved ***E***lements (***TRACE***)

![HALPER](https://img.shields.io/badge/HALPER-orange.svg)
![Python](https://img.shields.io/badge/python-3.6+-blue.svg)
![bedtools](https://img.shields.io/badge/bedtools-2.31+-orange.svg)
![rGREAT](https://img.shields.io/badge/rGREAT-2.0+-orange.svg)
![ggplot2](https://img.shields.io/badge/ggplot2-3.3+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![MEME-ChIP](https://img.shields.io/badge/MEME--ChIP-5.4.1-orange.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

***TRACE*** - A cross-species regulatory genomics pipeline for comparing open chromatin regions (OCRs) between human and mouse liver tissue using ATAC-seq data. The pipeline maps OCRs across species, classifies them as enhancers or promoters, identifies enriched biological processes, and discovers over-represented sequence motifs.

Developed for **03-713: Bioinformatics Data Integration Practicum**, Spring 2026, Carnegie Mellon University.

---

## Citation

If you use this repository, please cite:

**Evan Lin, Arun Sujatha Bharath Raj, Nikita Rajesh, Suratha Sriram**
(2026). *TRACE: Transcription Regulatory Analysis of Conserved Elements*.
03-713: Bioinformatics Data Integration Practicum, Spring 2026
Carnegie Mellon University

---

## Dependencies:

| Tool | Version | Purpose |
|------|---------|---------|
| [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing) | - | Cross-species liftover of OCRs |
| ↳ Python 3 *(via `hal` conda env)* | ≥ 3.6 | Required internally by HALPER - see [HAPLER installation](https://github.com/pfenninglab/halLiftover-postprocessing) |
| [bedtools](https://bedtools.readthedocs.io/en/latest/) | ≥ 2.31 | Genomic interval operations |
| [rGREAT](https://github.com/jokergoo/rGREAT) | ≥ 2.0 | GO enrichment analysis of genomic regions |
| [ggplot2](https://ggplot2.tidyverse.org/) | ≥ 3.3 | Visualization of rGREAT results |
| R | ≥ 4.0 | Required for rGREAT and ggplot2 |
| [MEME-ChIP](https://meme-suite.org/meme/tools/meme-chip) | 5.4.1 | Motif enrichment analysis |

> **Note:** All scripts are designed to run on the [Pittsburgh Supercomputing Center (PSC) Bridges-2](https://www.psc.edu/resources/bridges-2/) cluster. Tools are loaded via the `module` system. Running locally will require manual installation of all dependencies and may require minor path adjustments in the scripts.

## Installation

### On Bridges-2 (PSC) - recommended

Tools are available as modules on Bridges-2. No manual installation needed for most dependencies:

```bash
module load bedtools/2.31.0
module load MEME-suite/5.4.1
module load anaconda3
```

For HALPER and rGREAT, set up once:

```bash
# HALPER - clone and set up conda environment
# Follow: https://github.com/pfenninglab/halLiftover-postprocessing

# rGREAT - install inside your R environment
conda activate rgreat_env
Rscript -e "BiocManager::install('rGREAT')"
```

Then clone the repo:

```bash
git clone https://github.com/BioinformaticsDataPracticum2026/liver-ATAC-OCR.git
cd liver-ATAC-OCR
```

### Local installation

If running outside Bridges-2, install all dependencies manually:

```bash
# bedtools
conda install -c bioconda bedtools

# rGREAT (in R)
RScript -e "BiocManager::install('rGREAT')"

# HALPER
# Follow instructions at: https://github.com/pfenninglab/halLiftover-postprocessing

# MEME-ChIP
# Follow instructions at: https://meme-suite.org/meme/doc/install.html
```

> Minor path adjustments in the scripts may be required when running locally.

---

## Pipeline Overview
 
The pipeline addresses the following biological questions:
- To what extent is transcriptional regulatory activity conserved between human and mouse liver?
- How does conservation differ between enhancers and promoters?
- What biological processes are regulated by shared vs. species-specific OCRs?
- How does the regulatory code (TF motifs) differ between species, and between enhancers and promoters?

### Task 1 — Quality Control
 
ATAC-seq data from human and mouse liver and adrenal gland tissues were evaluated across seven QC metrics: % mapped reads, % properly paired reads, periodicity plots, TSS enrichment score, NRF, rescue ratio, and self-consistency ratio. The tissue with the highest-quality human and mouse datasets was selected for all downstream analyses.
 
**Result:** Human liver and mouse liver passed all QC criteria and were selected for downstream analysis. Human adrenal gland failed the reproducibility test and was excluded.
 
#### QC Metrics Summary
 
| Metric | Threshold | Mouse Liver R1 | Mouse Liver R2 | Human Liver R1 | Human Liver R2 | Mouse Adrenal R1 | Mouse Adrenal R2 | Human Adrenal R1 | Human Adrenal R2 |
|--------|-----------|---------------|---------------|----------------|----------------|-----------------|-----------------|-----------------|-----------------|
| % Mapped Reads | ≥ 90% | ✅ 97.8 | ✅ 97.9 | ✅ 98.9 | ✅ 98.6 | ✅ 98.5 | ✅ 98.6 | ✅ 98.5 | ✅ 98.9 |
| % Properly Paired | ≥ 90% | ✅ 95.7 | ✅ 95.6 | ✅ 97.9 | ✅ 97.5 | ✅ 94.6 | ✅ 94.7 | ✅ 97.0 | ✅ 97.5 |
| Periodicity Plot | 3 clear humps | ✅ 3 humps | ✅ 3 humps | ✅ 3 humps | ✅ 3 humps | ⚠️ 2 humps | ⚠️ 2 humps | ⚠️ 2 humps | ⚠️ 2 humps |
| TSS Enrichment | ≥ 10 | ⚠️ 7.69 | ⚠️ 7.35 | ✅ 23.29 | ✅ 20.96 | ✅ 18.37 | ✅ 19.00 | ✅ 26.08 | ✅ 14.01 |
| NRF | > 0.8 | ✅ 0.936 | ✅ 0.944 | ✅ 0.885 | ✅ 0.900 | ❌ 0.302 | ❌ 0.270 | ⚠️ 0.778 | ✅ 0.970 |
| Rescue Ratio | < 2 | ✅ 1.004 | ✅ 1.012 | ✅ 1.091 | ✅ 1.242 | ✅ 1.013 | ✅ 1.010 | ❌ 2.084 | ❌ 4.284 |
| Self-consistency Ratio | < 2 | ✅ 1.076 | ✅ 1.104 | ✅ 1.219 | ✅ 1.428 | ✅ 1.222 | ✅ 1.289 | ❌ 3.604 | ❌ 6.270 |
| Reproducibility | Pass | ✅ Pass | ✅ Pass | ✅ Pass | ✅ Pass | ✅ Pass | ✅ Pass | ❌ Fail | ❌ Fail |
 
> ✅ Passes threshold — ⚠️ Marginal — ❌ Fails threshold
 
### Task 2 — Cross-Species Liftover (HALPER)
 
Human liver OCRs are mapped to the mouse genome using HALPER, which leverages whole-genome alignments from the Cactus HAL file. Each human OCR is assigned an orthologous region in the mouse genome. OCRs are then classified as: **shared** (ortholog is open in mouse), or **human-specific** (ortholog is closed in mouse). Mouse-native OCRs with no human ortholog are classified as **mouse-specific**.
 
### Task 3 — Biological Process Enrichment (rGREAT)
 
GO biological process enrichment is performed using rGREAT on five sets of OCRs: all human OCRs, all mouse OCRs, shared OCRs, human-specific OCRs, and mouse-specific OCRs. rGREAT assigns regions to nearby genes based on genomic distance and tests for overrepresented GO terms. Results are filtered at adjusted p-value < 0.05.
 
### Task 4 — Promoter / Enhancer Classification
 
OCRs are classified as promoter-like or enhancer-like based on proximity to annotated transcription start sites (TSS ± 2 kb, from GENCODE vM15). Regions overlapping a TSS window are called promoters; all others are called enhancers. This classification is applied to shared, human-specific, and mouse-specific OCR sets.
 
### Task 5 — Motif Analysis (MEME-ChIP)
 
MEME-ChIP is run on FASTA sequences extracted from seven OCR sets: human enhancers, mouse enhancers, human promoters, mouse promoters, shared enhancers, human-specific enhancers, and mouse-specific enhancers. Sequences are resized to ±100 bp windows around peak midpoints. Discovered motifs are compared against the JASPAR 2026 vertebrates database to identify enriched transcription factor binding sites.
 
---

## Usage

### Run the full pipeline (Tasks 2-5)

```bash
bash scripts/TRACE_pipeline.sh \
    --human /path/to/human_liver.narrowPeak.gz \
    --mouse /path/to/mouse_liver.narrowPeak.gz \
    --hal /path/to/10plusway-master.hal \
    --tss /path/to/gencode.vM15.TSSWithStrand_sorted.bed \
    --genome /path/to/mma.fa \
    --jaspar /path/to/JASPAR2026_vertebrates.meme
```

### Available flags

| Flag | Description | Default |
|------|-------------|---------|
| `--human` | Human ATAC-seq peak file (.narrowPeak.gz) | - |
| `--mouse` | Mouse ATAC-seq peak file (.narrowPeak.gz) | - |
| `--hal` | HAL alignment file | - |
| `--tss` | TSS annotation BED file | - |
| `--genome` | mm10 genome FASTA | - |
| `--jaspar` | JASPAR motif database (.meme) | - |
| `--source-species` | Liftover source species | `Human` |
| `--target-species` | Liftover target species | `Mouse` |
| `--conda-env` | Conda environment for rGREAT | `rgreat_env` |
| `--skip-halper` | Skip liftover step | - |
| `--skip-pe` | Skip P/E classification step | — |
| `--skip-great` | Skip rGREAT steps | — |
| `--skip-motif` | Skip motif analysis step | — |
| `-h`, `--help` | Show this help message and exit | - |

### Run individual steps

```bash
bash scripts/run_halper_mapping.sh     # Task 2: HALPER mapping
bash scripts/run_pe_classification.sh  # Task 4: P/E classification
bash scripts/run_rgreat.sh             # Task 3: rGREAT enrichment
bash scripts/run_plots.sh              # Task 3: rGREAT plots
bash scripts/run_motif_analysis.sh     # Task 5: Motif enrichment
```

### External tool documentation

- HALPER: https://github.com/pfenninglab/halLiftover-postprocessing
- rGREAT: https://github.com/jokergoo/rGREAT
- MEME-ChIP: https://meme-suite.org/meme/tools/meme-chip

### Saving outputs

All results are written automatically to their respective output directories. To archive:

```bash
tar -czvf results_backup.tar.gz Mapping/outputs PE_classification/output rGREAT/outputs Motif_analysis/outputs
```

---

## Output

| Directory | Contents |
|-----------|----------|
| `Mapping/outputs/` | HALPER liftover result (`*.HALPER.narrowPeak.gz`) |
| `PE_classification/outputs/rowcount/summary.tsv` | Enhancer/promoter counts (raw) |
| `PE_classification/outputs/unique/summary.tsv` | Enhancer/promoter counts (deduplicated) |
| `rGREAT_Analysis/outputs/plots/` | Barplots and comparison heatmap (`.png`) |
| `Motif_analysis/outputs/meme_chip_human_specific_enhancer/` | `summary.tsv` + `meme-chip.html` |
| `Motif_analysis/outputs/meme_chip_mouse_specific_enhancer/` | `summary.tsv` + `meme-chip.html` |
| `logs/` | Timestamped per-step log files |

---

## Repository Structure

```
liver-ATAC-OCR/
├── scripts/
│   ├── TRACE_pipeline.sh
│   ├── run_halper_mapping.sh
│   ├── run_pe_classification.sh
│   ├── run_rgreat.sh
│   ├── run_plots.sh
│   └── run_motif_analysis.sh
├── rGREAT_Analysis/
│   ├── scripts/
│   │   ├── rgreat_analysis.R
│   │   └── rgreat_plots.R
│   └── outputs/
├── Mapping/outputs/
├── PE_classification/
│   ├── input/
│   └── outputs/
├── Motif_analysis/outputs/
├── docs/QC_Analysis.md
├── README.md
├── requirements.txt
├── .gitignore
└── LICENSE
```

---

## References
 
### Tool papers
 
- Diehl, A.G. et al. (2020). HALPER: a tool for cross-species liftover of ATAC-seq peaks. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btaa493
- Quinlan, A.R. & Hall, I.M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26(6):841–842. https://doi.org/10.1093/bioinformatics/btq033
- McLean, C.Y. et al. (2010). GREAT improves functional interpretation of cis-regulatory regions. *Nature Biotechnology* 28, 495–501. https://doi.org/10.1038/nbt.1630
- Gu, Z. & Hübschmann, D. (2023). rGREAT: an R/Bioconductor package for functional enrichment of genomic regions. *Bioinformatics* 39(1). https://doi.org/10.1093/bioinformatics/btac745
- Machanick, P. & Bailey, T.L. (2011). MEME-ChIP: motif analysis of large DNA datasets. *Bioinformatics* 27(12):1696–1697. https://doi.org/10.1093/bioinformatics/btr189
- Bailey, T.L. et al. (2015). The MEME Suite. *Nucleic Acids Research* 43(W1):W39–49. https://doi.org/10.1093/nar/gkv416

### Cross-species regulatory conservation
 
- Villar, D. et al. (2015). Enhancer evolution across 20 mammalian species. *Cell* 160(3):554–566. https://doi.org/10.1016/j.cell.2015.01.006
- Ballester, B. et al. (2014). Multi-species, multi-transcription factor binding highlights conserved control of tissue-specific biological pathways. *eLife* 3:e02626. https://doi.org/10.7554/eLife.02626
- Yue, F. et al. (2014). A comparative encyclopedia of DNA elements in the mouse genome. *Nature* 515, 355–364. https://doi.org/10.1038/nature13992

### Biological interpretation — GO terms
 
- Wahli, W. & Michalik, L. (2012). PPARs at the crossroads of lipid signaling and inflammation. *Trends in Endocrinology & Metabolism* 23(7):351–363. https://doi.org/10.1016/j.tem.2012.05.001
- Rui, L. (2014). Energy Metabolism in the Liver. *Comprehensive Physiology* 4(1):177–197. https://doi.org/10.1002/cphy.c130024
- Lu, S.C. (2009). Regulation of Glutathione Synthesis. *Molecular Aspects of Medicine* 30(1-2):42–59. https://doi.org/10.1016/j.mam.2008.05.005
- Mark, M. et al. (2009). Retinoic Acid Signalling in Development. *Physiological Reviews* 89(2):843–884. https://doi.org/10.1152/physrev.00026.2008
- Blaner, W.S. et al. (2009). Hepatic stellate cell lipid droplets: A specialized lipid droplet for retinoid storage. *Biochimica et Biophysica Acta* 1791(6):467–473. https://doi.org/10.1016/j.bbalip.2008.11.001
- Dzierzak, E. & Philipsen, S. (2013). Erythropoiesis: Development and Differentiation. *Cold Spring Harbor Perspectives in Medicine* 3(4). https://doi.org/10.1101/cshperspect.a011601

### Biological interpretation — TF motifs
 
- Ong, C.T. & Corces, V.G. (2014). CTCF: an architectural protein bridging genome topology and function. *Nature Reviews Genetics* 15, 234–246. https://doi.org/10.1038/nrg3663
- Schmidt, D. et al. (2012). Waves of retrotransposon expansion remodel genome organization and CTCF binding in multiple mammalian lineages. *Cell* 148(1-2):335–348. https://doi.org/10.1016/j.cell.2011.11.058
- Pawlak, M. et al. (2015). Molecular mechanism of PPARα action and its impact on lipid metabolism, inflammation and fibrosis in non-alcoholic fatty liver disease. *Journal of Hepatology* 62(3):720–733. https://doi.org/10.1016/j.jhep.2014.10.039
- Lefebvre, P. et al. (2006). Sorting out the roles of PPARα in energy metabolism and vascular homeostasis. *Journal of Clinical Investigation* 116(3):571–580. https://doi.org/10.1172/JCI27989
- Suske, G. (1999). The Sp-family of transcription factors. *Gene* 238(2):291–300. https://doi.org/10.1016/S0378-1119(99)00357-1
- Golson, M.L. & Kaestner, K.H. (2016). Fox transcription factors: from development to disease. *Development* 143(24):4558–4570. https://doi.org/10.1242/dev.112672

---

## Contact

**Evan Lin** - evanlin@andrew.cmu.edu
**Arun Sujatha Bharath Raj** - asujatha@andrew.cmu.edu
**Nikita Rajesh** - nrajesh@andrew.cmu.edu
**Suratha Sriram** - surathas@andrew.cmu.edu

---

## License:
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

---
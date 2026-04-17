library(rGREAT) #loads rGREAT package

base <- "/ocean/projects/bio230007p/ssriram6/liver-ATAC-OCR" #project directory path
outdir <- file.path(base, "rGREAT_Analyis/outputs") #path for output directory
dir.create(outdir, showWarnings = FALSE) #creates the output directory, won't complain if the folder already exists

load_gr <- function(f) { #helper function that takes a file path f
    d <- read.table(f, sep = "\t", header = FALSE) #reads the narrowpeak file as a tab separated table with no header row, d is a dataframe where col 1 is chromosome, col 2 is start position and col 3 is end position
    GRanges(seqnames=d[,1], ranges = IRanges(d[,2],d[,3])) #converts the dataframe into a GRanges object(bioconductors standard format for genomic regions), seqnames is the chromosome name, Iranges defines the coordinate range for each region
} 

inputs <- list(
  list("human_all",      file.path(base, "Mapping/outputs/human_liver.narrowPeak"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_all",      file.path(base, "Mapping/outputs/mouse_liver.narrowPeak"), "TxDb.Mmusculus.UCSC.mm10.knownGene"),
  list("shared",         file.path(base, "PE_classification/output/unique/shared_open.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("human_specific", file.path(base, "PE_classification/output/unique/human_open_mouse_closed.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_specific", file.path(base, "PE_classification/output/unique/mouse_open_human_closed.bed"), "TxDb.Mmusculus.UCSC.mm10.knownGene")
) #creates a list of 4 sub-lists containing a label name, file path and genome assembly

for(inp in inputs) { #iterates through the 4 input sets
    nm <- inp[[1]]; fp <- inp[[2]];gn <-inp[[3]] #extracts 3 elements - name label, file path and genome assembly
    cat("Running:", nm, "\n") #to know which analysis is currently running
    gr <- load_gr(fp) #calls the helper function to convert the BED file and convert it to a GRanges object
    job <- great(gr, gene_sets="GO:BP", tss_source=gn) #runs the rGREAT algorithm 
    res <- getEnrichmentTable(job) #extracts full results from the job
    sig <- res[res$p_adjust < 0.05, ] #filters the results to keep only rows where the adjusted p-value is below 0.05
    cat(" ", nrow(sig), "significant terms\n") #prints number of significant terms
    write.csv(res, file.path(outdir, paste0(nm, "_GOBP.csv")), row.names=FALSE) #saves the results as csv
    write.csv(sig, file.path(outdir, paste0(nm, "_GOBP_sig.csv")), row.names=FALSE) #saves only significant results
}

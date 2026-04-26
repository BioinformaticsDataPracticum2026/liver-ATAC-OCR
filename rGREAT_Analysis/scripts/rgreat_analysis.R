library(rGREAT) #loads rGREAT package

args <- commandArgs(trailingOnly = TRUE) #captures command-line arguments passed to the script, trailingOnly=TRUE means only the arguments after --args are captured (excludes R's own flags)
base <- args[1] #project directory path
data_root <- args[2] #root path to where the ATAC-seq peak files are stored
outdir <- file.path(base, "rGREAT_Analysis/outputs") #path for output directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE) #creates the output directory, won't complain if the folder already exists

inputs <- list(
  list("human_all",      file.path(data_root, "HumanAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_all",      file.path(data_root, "MouseAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"), "TxDb.Mmusculus.UCSC.mm10.knownGene"),
  list("shared",         file.path(base, "PE_classification/outputs/unique/shared_open.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("human_specific", file.path(base, "PE_classification/outputs/unique/human_open_mouse_closed.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_specific", file.path(base, "PE_classification/outputs/unique/mouse_open_human_closed.bed"), "TxDb.Mmusculus.UCSC.mm10.knownGene")
)  #creates a list of 4 sub-lists containing a label name, file path and genome assembly

load_gr <- function(f) { #helper function that takes a file path f
    con <- if (grepl("\\.gz$", f)) gzfile(f) else f #checks if file ends in .gz, if so opens as gzip connection, otherwise reads as plain text
    d <- read.table(con, sep = "\t", header = FALSE) #reads the narrowpeak file as a tab separated table with no header row
    GRanges(seqnames=d[,1], ranges = IRanges(d[,2],d[,3])) #converts the dataframe into a GRanges object
}

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

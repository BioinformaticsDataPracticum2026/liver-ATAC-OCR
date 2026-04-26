library(ggplot2)  
#loads the ggplot2 package for creating plots

args <- commandArgs(trailingOnly = TRUE) #captures command-line arguments passed to the script
base <- args[1] #project directory path

indir  <- file.path(base, "rGREAT_Analysis/outputs")
#path to folder containing GO enrichment CSV results
outdir <- file.path(base, "rGREAT_Analysis/outputs/plots")
#path where generated plots will be saved
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
#creates the plots folder if it doesn't already exist
#showWarnings suppresses the warning if the folder already exists
#recursive creates creates any parent directories if needed

conditions <- c("human_all","mouse_all","shared","human_specific","mouse_specific") 
#vector for 5 analysis conditions
results <- lapply(conditions, function(c) read.csv(file.path(indir, paste0(c, "_GOBP_sig.csv")))) 
#lappy loops over each condition name
#for each one, it builds the filename
#file.path joins the directory and filename into a full path
#read.csv reads that csv into a data frame

names(results) <- conditions
#names each element in the list so it can be accessed like results[["human_all"]]

safe_log <- function(p) {
  p[p == 0] <- 1e-20
  -log10(p)
}
#Helper function that takes p-values and returns -log10 of them
#If any p-value is exactly 0, it replaces it with 1e-20 first to avoid -log10(0)=Infinity

for (cond in conditions) {
  #loops through each of the 5 conditions to make a bar plot for each
  df <- head(results[[cond]][order(results[[cond]]$p_adjust), ], 15)
  #takes the dataframe for this condition, order sorts rows by p_adjust, head keeps only the top 15 most significant terms
  df$description <- factor(df$description, levels=rev(df$description))
  #converts the GO term names to a factor with reversed order, this controls the y-axis order so the most significant term appears at the top
  df$neg_log_p <- safe_log(df$p_adjust)
  #adds a new column with -log10(p_adjust) for use as the fill colour
  p <- ggplot(df, aes(fold_enrichment, description, fill=neg_log_p)) +
  #Initializes the plot, x-axis = fold enrichment, y-axis=GO term name, fill colour = statistical significance
    geom_col() +
    #draws hortizontal bars
    scale_fill_gradient(low="peachpuff", high="darkred", name="-log10\n(p_adjust)") +
    #sets the colour gradient : light peachpuff for less significant, dark red for more significant, name sets the legend title
    labs(title=cond, x="Fold Enrichment", y=NULL) +
    #title = condition name, x = x-axis label, y=NULL removes the y-axis label since terms are self-explanatory
    theme_minimal() +
    #applies a clean minimal theme
    theme(
      plot.title = element_text(face="bold", size=14),
      #makes the title bold and size 14
      axis.text.y = element_text(size=9)
      #sets the GO term label size to 9
    )
  ggsave(file.path(outdir, paste0("barplot_", cond, ".png")), p, width=9, height=6, dpi=150)
  #saves the plot as a PNG file, width = 9 inches, height= 6 inches, resolution = 150 dots per inch
}
library(ggplot2)

indir  <- "/ocean/projects/bio230007p/ssriram6/liver-ATAC-OCR/rGREAT_Analysis/outputs"
outdir <- "/ocean/projects/bio230007p/ssriram6/liver-ATAC-OCR/rGREAT_Analysis/outputs/plots"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

conditions <- c("human_all","mouse_all","shared","human_specific","mouse_specific")
results <- lapply(conditions, function(c) read.csv(file.path(indir, paste0(c, "_GOBP_sig.csv"))))
names(results) <- conditions

safe_log <- function(p) {
  p[p == 0] <- 1e-20
  -log10(p)
}

for (cond in conditions) {
  df <- head(results[[cond]][order(results[[cond]]$p_adjust), ], 15)
  df$description <- factor(df$description, levels=rev(df$description))
  df$neg_log_p <- safe_log(df$p_adjust)
  p <- ggplot(df, aes(fold_enrichment, description, fill=neg_log_p)) +
    geom_col() +
    scale_fill_gradient(low="peachpuff", high="darkred", name="-log10\n(p_adjust)") +
    labs(title=cond, x="Fold Enrichment", y=NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face="bold", size=14),
      axis.text.y = element_text(size=9)
    )
  ggsave(file.path(outdir, paste0("barplot_", cond, ".png")), p, width=9, height=6, dpi=150)
}

top_terms <- unique(unlist(lapply(results, function(d) head(d[order(d$p_adjust), "description"], 10))))
hm <- do.call(rbind, lapply(conditions, function(c) {
  df <- results[[c]]
  data.frame(description=top_terms, category=c,
    value=sapply(top_terms, function(t) {
      v <- df$fold_enrichment[df$description == t]
      if (length(v)==0) NA else v[1]
    }))
}))
hm$category <- factor(hm$category, levels=conditions)

p <- ggplot(hm, aes(category, description, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient(low="white", high="darkred", na.value="grey90", name="Fold\nEnrichment") +
  labs(title="GO:BP comparison", x=NULL, y=NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face="bold", size=14),
    axis.text.x = element_text(angle=30, hjust=1, size=10),
    axis.text.y = element_text(size=7),
    panel.grid = element_blank()
  )
ggsave(file.path(outdir, "comparison_heatmap.png"), p, width=12, height=14, dpi=150)



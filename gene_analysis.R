library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(dplyr)

edb <- EnsDb.Hsapiens.v86
gene_id <- "ENSG00000091513"

cat("Extracting transcripts...\n")
tx <- transcripts(edb, filter = GeneIdFilter(gene_id))
tx_df <- as.data.frame(tx)
write.csv(tx_df, "transcripts.csv", row.names = FALSE)

#protein
cat("Extracting proteins...\n")
prot <- proteins(edb, filter = GeneIdFilter(gene_id))
prot_df <- as.data.frame(prot)
write.csv(prot_df, "proteins.csv", row.names = FALSE)


#protein → genome mapping
if (nrow(prot_df) > 0) {
  cat("Mapping proteins to genome...\n")
  
  pr <- IRanges(start = rep(1, nrow(prot_df)), end = prot_df$protein_length)
  
  names(pr) <- prot_df$protein_id

  res <- tryCatch({
    proteinToGenome(pr, edb)
  }, error = function(e) {
    message("Mapping error: ", e$message)
    return(NULL)
  })

  if (!is.null(res)) {
    res_grl <- GRangesList(res)
    map_df <- as.data.frame(res_grl)

    map_df_clean <- map_df %>%
      select(-group, -group_name) 

    write.csv(map_df_clean, "protein_genome_map.csv", row.names = FALSE)

    cat("\nPreview of mapping data:\n")
    print(head(map_df_clean))
  }

} else {
  print("No protein-coding transcripts found.")
}

cat("\n--- Summary ---\n")
cat("Transcripts extracted:", nrow(tx_df), "\n")
cat("Proteins extracted:", nrow(prot_df), "\n")
cat("Success! Files saved to your working directory.\n")
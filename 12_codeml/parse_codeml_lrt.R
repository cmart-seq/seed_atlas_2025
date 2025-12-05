library(stringr)
library(dplyr)
library(purrr)
library(readr)

#model pairs and degrees of freedom
model_pairs <- list(
  M2a_vs_M1a = list(alt = 3, null = 2, df = 2),
  M8_vs_M7 = list(alt = 5, null = 4, df = 2)
)

root_dir <- "orthofinder/clean/orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees"
#A function to extract LNLs and omega values from a codeml.out file
extract_codeml_info <- function(file_path) {
  cat(paste0("Processing ", file_path, "\n"))
  text_lines <- readLines(file_path)
  #log-likelihoods
  lnl_matches <- str_split(grep("lnL", text_lines, value = T), "\\s+", simplify = T)
  #catch anything that didn't finish
  
  if(nrow(lnl_matches) == 0 || ncol(lnl_matches) < 5){
    lnls <- rep(NA, 5)
    omegas <- list(
      p0 = NA,
      p1 = NA,
      p2 = NA,
      w0 = NA,
      w1 = NA, 
      w2 = NA
    )
    return(list(lnls = lnls, omegas = omegas))
  }else{
    # [1] "Model 0: one-ratio"                       
    # [2] "Model 1: NearlyNeutral (2 categories)"    
    # [3] "Model 2: PositiveSelection (3 categories)"
    # [4] "Model 7: beta (10 categories)"            
    # [5] "Model 8: beta&w>1 (11 categories)" 
    lnls <- as.numeric(unlist(lnl_matches[,5]))
    
    #omega values from model2
    omega_idx <- grep("for site classes \\(K=3\\)", text_lines)
    all_idx <- sort(unique(unlist(lapply(omega_idx, function(i) i:(i+3)))))
    omega_match <- text_lines[all_idx[c(3, 4)]]
    
    omegas <- list(
      p0 = as.numeric(str_split(omega_match, "\\s+")[[1]][2]),
      p1 = as.numeric(str_split(omega_match, "\\s+")[[1]][3]),
      p2 = as.numeric(str_split(omega_match, "\\s+")[[1]][4]),
      
      w0 = as.numeric(str_split(omega_match, "\\s+")[[2]][2]),
      w1 = as.numeric(str_split(omega_match, "\\s+")[[2]][3]),
      w2 = as.numeric(str_split(omega_match, "\\s+")[[2]][4])
    )
  return(list(lnls = lnls, omegas = omegas))
  }
}
    
#A function to perform LRT and return test statistic and p-value
run_lrt <- function(lnl_alt, lnl_null, df) {
  stat <- 2 * (lnl_alt - lnl_null)
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  return(list(stat = round(stat, 3), pval = round(pval, 6)))
}

#recursively read .codeml.out files and extract data
all_dirs <- list.dirs(root_dir, recursive = FALSE)
results <- list()
raw_pvals <- list()


for (subdir in all_dirs) {
    files <- list.files(subdir, pattern = "\\.codeml\\.out$", full.names = TRUE)
    for (file in files) {
      info <- extract_codeml_info(file)
      lnls <- info$lnls
      omegas <- info$omegas

      # Skip if any NA in lnls
      if (any(is.na(lnls))) next
      
      row <- list(Gene = basename(subdir))
      
      for (name in names(model_pairs)) {
        alt <- model_pairs[[name]]$alt
        null <- model_pairs[[name]]$null
        df <- model_pairs[[name]]$df
        
        if (length(lnls) >= max(alt, null)) {
          res <- run_lrt(lnls[alt], lnls[null], df)
          row[[paste0(name, "_null_raw")]] <- lnls[null]
          row[[paste0(name, "_alt_raw")]] <- lnls[alt]
          row[[paste0(name, "_LRT")]] <- res$stat
          row[[paste0(name, "_pval")]] <- res$pval
          raw_pvals <- append(raw_pvals, list(list(gene = row$Gene, model = name, pval = res$pval)))
        } else {
          row[[paste0(name, "_null_raw")]] <- lnls[null]
          row[[paste0(name, "_alt_raw")]] <- lnls[alt]
          row[[paste0(name, "_LRT")]] <- res$stat
          row[[paste0(name, "_LRT")]] <- NA
          row[[paste0(name, "_pval")]] <- NA
        }
      }
      
      #Add omega values
      row <- c(row, omegas)
      results <- append(results, list(row))
    }
}

# Convert to dataframe
df <- bind_rows(results)

# Adjust p-values by model
for (model_name in names(model_pairs)) {
  pvals <- df[[paste0(model_name, "_pval")]]
  padj <- rep(NA, length(pvals))
  valid <- !is.na(pvals)
  padj[valid] <- p.adjust(pvals[valid], method = "fdr")
  df[[paste0(model_name, "_padj")]] <- round(padj, 6)
}

# Save results
write_csv(df, "orthofinder_inputs4/codeml_LRT_summary_full.csv")

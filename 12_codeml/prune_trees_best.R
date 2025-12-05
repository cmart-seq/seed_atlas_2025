#preserve only the one best ortholog from each species to a single Arabidopsis thaliana gene

library(ape)

#getting all of the tree files from OrthoFinder outputs
tree_files <- list.files(pattern = "_tree\\.txt$")

for (file in tree_files) {
  tree <- read.tree(file)

  #Get species from tip labels
  species <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)

  #Identify Arabidopsis tips (starts with Athaliana_)
  arabidopsis_tips <- tree$tip.label[species == "Cleanthaliana"]
  other_species <- unique(species[species != "Cleanthaliana"])

  if (length(arabidopsis_tips) == 0) next

  #getting a distance matrix for tree members
  dist_matrix <- cophenetic(tree)

  #iterating through tips to find the one with the minimum distance
  for (ath_tip in arabidopsis_tips) {
    best_hits <- c()
    for (spp in other_species) {
      spp_tips <- tree$tip.label[species == spp]
      if (length(spp_tips) > 0) {
        # Get distances to the Arabidopsis tip
        dists <- dist_matrix[ath_tip, spp_tips]
        closest <- spp_tips[which.min(dists)]
        best_hits <- c(best_hits, closest)
      }
    }

    keep_tips <- c(ath_tip, best_hits)
    
    #not saving if the tree has less than three members
    if (length(keep_tips) >= 3) {
      pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))
      gene_id <- sub("Cleanthaliana_", "", ath_tip)
      output_name <- gsub("_tree.txt", paste0("_", gene_id, "_pruned.tree"), file)
      write.tree(pruned_tree, file=output_name)
    }
  }
}

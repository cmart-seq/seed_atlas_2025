#OrthoFinder appends prefixes to tree tip genes that do not match protein and CDS names in the FASTA files, 
#this script removes all prefixes from gene trees.

import os

#List of prefixes to remove
prefixes = ["Cleanlyrata_Proteins_","Cleanbrassica_Proteins_", "Cleancapsella_Proteins_", "Cleanthaliana_Proteins_"]

#Function to remove prefixes from tree file, iterates over each prefix and remove it if present
def remove_prefixes_from_tree(treefile):
    with open(treefile, 'r') as f:
        tree = f.read()
    for prefix in prefixes:
        tree = tree.replace(prefix, "")
    with open(treefile, 'w') as f:
        f.write(tree)

tree_dir = 'OrthoFinder/Results_Apr17/Gene_Trees'

#Looping over all tree files
for filename in os.listdir(tree_dir):
    if filename.endswith(".tree"):
        treefile_path = os.path.join(tree_dir, filename)
        remove_prefixes_from_tree(treefile_path)
        print(f"Processed: {filename}")

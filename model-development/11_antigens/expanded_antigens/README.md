This folder contains the antigenicity tag details after expanding the antigenicity for both proteins (using BLAST) and peptides (using kmer expansion). While these files can be created from scratch using the data and code provided in *model-development*, we put this information here for ease of use for anyone interested specifically in the antigenicity.

### **Proteins:**
We created similarity clusters between the proteins using BLAST. After doing an all vs all BLAST we kept matches with a pident >= 0.75, an evalue <= 1 x 10-12 and a match length of at least half of the length of the shortest protein in the match. Using these matches, we created a distance matrix where distance = 1 – pident and applied a single-linkage hierarchical clustering method. We then cut this tree using a cutoff of 0.25, resulting in a set of clusters of similar proteins.

Inside each folder there is a file called *protein_antigenic_tag.tsv*. These are its columns:

- **protein**: id of the protein
- **is_antigenic**: if we had bibliography indicating that the protein is antigenic (0 or 1)
- **similarity_cluster**: the ID of the cluster where the protein was found after the BLAST. These clusters are shared between all 15 species
- **in_antigenic_cluster**: if at least 1 of the proteins in the cluster had information in bibliography indicating that the protein is antigenic (0 or 1)
- **selfSimilarity_cluster**: the ID of the cluster where the protein was found after the BLAST, using only proteins from this specific species
- **in_selfAntigenic_cluster**: if at least 1 of the proteins in the self-similarity cluster had information in bibliography indicating that the protein is antigenic (0 or 1)

### **Peptides:**
Using a method we called "kmer expansion", we expanded the antigenicity from peptides found in bibliography to others in that same protein. This method works by marking as antigenic any peptide that shared a kmer of at least 8 amino acids with a curated antigenic sequence for that same protein.

Inside most folders there is a file called peptide_antigenic_tag.tsv. Some folders don't have this file because we didn't have antigenicity information at peptide level for that species. Due to size constraints, we only put the antigenic peptides in these files. All other 15mers were considered to be negative when training the models.

These are peptide_antigenic_tag.tsv columns:

- **protein** : id of the protein
- **start** : the position of the first amino acid of the peptide inside the protein (the first peptide has a start of 1)
- **peptide** : the sequence of the peptide
- **is_antigenic** : if we had bibliography indicating that the peptide is antigenic (0 or 1, here is all 1 due to the filtering)

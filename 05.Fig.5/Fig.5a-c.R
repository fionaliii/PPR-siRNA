library(tidyverse)
library(readr)
library(Biostrings)
library(readxl)
library(ggtree)



## 1. Import the taxonomy metadata and reference sequences:
RefSeq_metadata <- read_xlsx("/path/to/Data_SX1.csv")

Species_meta <- RefSeq_metadata[,c(1,3,40,42,47)]



## 2. Import the reference and outgroup PPRs:
reference_PPRs <- readAAStringSet("/path/to/reference_PPRs.fasta")

AtPPR <- readAAStringSet("/path/to/AtPPR_ref.fasta")

outgroup_PPR <- readAAStringSet("/path/to/Outgroup_PPR.fasta")



## 3. Import the HMMsearch outputs for different profiles:
# PPR
PPR_tbl <- read_table("/path/to/PF01535_tbl.out",
                      col_names = FALSE, col_types = cols(X19 = col_skip()), 
                      skip = 3) %>% head(-10)

PPR_tbl <- filter(PPR_tbl, PPR_tbl$X5 <= 0.01) # keep the significant ones
PPR <- PPR_tbl$X1 %>% as.data.frame() # just the ids


# PPR_1
PPR_1_tbl <- read_table("/path/to/PF12854_tbl.out",
                        col_names = FALSE, col_types = cols(X19 = col_skip()), 
                        skip = 3) %>% head(-10)

PPR_1_tbl <- filter(PPR_1_tbl, PPR_1_tbl$X5 <= 0.01) # keep the significant ones
PPR_1 <- PPR_1_tbl$X1 %>% as.data.frame() # just the ids


# PPR_2
PPR_2_tbl <- read_table("/path/to/PF13041_tbl.out",
                        col_names = FALSE, col_types = cols(X19 = col_skip()), 
                        skip = 3) %>% head(-10)

PPR_2_tbl <- filter(PPR_2_tbl, PPR_2_tbl$X5 <= 0.01) # keep the significant ones
PPR_2 <- PPR_2_tbl$X1 %>% as.data.frame() # just the ids


# PPR_3
PPR_3_tbl <- read_table("/path/to/PF13812_tbl.out",
                        col_names = FALSE, col_types = cols(X19 = col_skip()), 
                        skip = 3) %>% head(-10)

PPR_3_tbl <- filter(PPR_3_tbl, PPR_3_tbl$X5 <= 0.01) # keep the significant ones
PPR_3 <- PPR_3_tbl$X1 %>% as.data.frame() # just the ids



## 4. Merge the outputs and extract the sequences:
PPR_all <- Reduce(function(x,y) merge(x,y, by = ".", all = TRUE), list(PPR,PPR_1,PPR_2,PPR_3)) # merge all PPR containing proteins

# export the IDs and extract the from RefSeq proteomes
write(PPR_all$., "/path/to/PPR_ids.txt")



## 5. Extract the sequences and import the in R:
#### * When extracting the sequences we added the assembly accession to the sequences as part of the sequence name.
#### * We will use that to connect sequences to genome metadata
PPR_seq <- readAAStringSet("/path/to/PPR.fasta")

PPR_seq_ids <- PPR_seq@ranges@NAMES %>% as.data.frame() %>% setNames("seqname")

PPR_seq_ids <- PPR_seq_ids %>% 
  separate(seqname, into = c("seqname", "refseq"), sep = " ", extra = "merge", fill = "right")

PPR_seq_ids <- PPR_seq_ids %>%
  mutate(refseq = str_remove(refseq, "\\[folder=")) %>%
  mutate(refseq = str_remove(refseq, "\\]"))

PPR_seq_meta <- PPR_seq_ids %>% left_join(Species_meta, join_by(refseq == `Assembly Accession`))

PPR_seq_meta <- PPR_seq_meta %>% 
  mutate(ID = paste(`Order name`, `Organism Name`, seqname, sep = "_")) %>%
  mutate(ID = gsub(" ", "_", ID))

PPR_seq@ranges@NAMES <- PPR_seq_meta$ID



## 6. Clean up and filter out unwated and duplicated sequences:
PPR_seq_filtered <- PPR_seq[PPR_seq@ranges@width > 150]
PPR_seq_filtered <- PPR_seq_filtered[PPR_seq_filtered@ranges@width < 1500]

# only keep the non-redundant sequences
PPR_seq_filtered_unique <- unique(PPR_seq_filtered)

# export
writeXStringSet(PPR_seq, "/path/to/PPR.fasta")
writeXStringSet(PPR_seq_filtered, "/path/to/PPR_filtered.fasta")
writeXStringSet(PPR_seq_filtered_unique, "/path/to/PPR_filtered_unique.fasta")

# add the references and export for phylogenetic analysis
all_seq_ref_out <- c(PPR_seq_filtered, reference_PPRs, outgroup_PPR)
writeXStringSet(all_seq_ref_out, "/path/to/PPR_filtered_ref_out.fasta")

all_seq_ref_out_unique <- c(PPR_seq_filtered_unique, reference_PPRs, outgroup_PPR)
writeXStringSet(all_seq_ref_out_unique, "/path/to/PPR_filtered_unique_ref_out.fasta")


# export the supplementary sequence and metadata
PPR_seq_meta_seq <- PPR_seq_meta
PPR_seq_meta_seq$sequence <- PPR_seq
write_csv(PPR_seq_meta_seq, "/path/to/Data_SX_PPR.csv")


PPR_seq_meta_filtered <- PPR_seq_meta[PPR_seq_meta$ID %in% PPR_seq_filtered@ranges@NAMES,]
PPR_seq_meta_filtered_seq <- PPR_seq_meta_filtered
PPR_seq_meta_filtered_seq$sequence <- PPR_seq_filtered
write_csv(PPR_seq_meta_filtered_seq, "/path/to/Data_SX_PPR_filtered.csv")


PPR_seq_meta_filtered_unique <- PPR_seq_meta[PPR_seq_meta$ID %in% PPR_seq_filtered_unique@ranges@NAMES,]
PPR_seq_meta_filtered_unique_seq <- PPR_seq_meta_filtered_unique
PPR_seq_meta_filtered_unique_seq$sequence <- PPR_seq_filtered_unique
write_csv(PPR_seq_meta_filtered_unique_seq, "/path/to/Data_SX_PPR_filtered_unique.csv")




## 5. Extract the PPR-rich clade and remake a new tree:
#### * anything with "out" in the name includes the paraphyletic monocot clade
#### * anything with "para_clade" in the name is only the paraphyletic clade
# import the trees
PPR_rich_out_clade_tree <- read.tree("/path/to/sRNA_rich_clade_out.tree")
PPR_rich_clade_tree <- read.tree("/path/to/sRNA_rich_clade.tree")

# extract the sequences for each clade
PPR_rich_out_clade_seq <- PPR_seq_filtered_unique[PPR_seq_filtered_unique@ranges@NAMES %in% PPR_rich_out_clade_tree$tip.label]
PPR_rich_clade_seq <- PPR_seq_filtered_unique[PPR_seq_filtered_unique@ranges@NAMES %in% PPR_rich_clade_tree$tip.label]
PPR_rich_para_clade_seq <- PPR_rich_out_clade_seq[!PPR_rich_out_clade_seq@ranges@NAMES %in% PPR_rich_clade_tree$tip.label]

PPR_rich_clade_seq_ref <- all_seq_ref_out_unique[PPR_rich_clade_tree$tip.label]


# extract the metadata for each clade
PPR_rich_out_clade_meta <- PPR_seq_meta[PPR_seq_meta$ID %in% PPR_rich_out_clade_tree$tip.label,]
PPR_rich_clade_meta <- PPR_seq_meta[PPR_seq_meta$ID %in% PPR_rich_clade_tree$tip.label,]
PPR_rich_para_clade_meta <- PPR_rich_out_clade_meta[!PPR_rich_out_clade_meta$ID %in% PPR_rich_clade_tree$tip.label,]

PPR_rich_clade_meta$sequence <- PPR_rich_clade_seq
PPR_rich_para_clade_meta$sequence <- PPR_rich_para_clade_seq

# export the sequence and metadata
write_csv(PPR_rich_clade_meta, "/path/to/PPR_sRNA_clade_meta.csv")
write_csv(PPR_rich_para_clade_meta, "/path/to/PPR_sRNA_para_clade_meta.csv")

# Extract the siRNA-rich reference PPRs
reference_PPRs_rich_metadata <- reference_PPRs_metadata[reference_PPRs_metadata$Pep_ID %in% PPR_rich_out_clade_tree$tip.label,]

PPR_rich_reference_seq <- reference_PPRs[reference_PPRs@ranges@NAMES %in% reference_PPRs_rich_metadata$Pep_ID]


# export the sRNA-rich clade + paraphyletic clade + outgroups
c(PPR_rich_out_clade_seq, PPR_rich_reference_seq, LPPRC) %>% 
  writeXStringSet("/path/to/PPR_rich_ref_out.fasta")

c(PPR_rich_out_clade_seq, LPPRC) %>% 
  writeXStringSet("/path/to/PPR_rich_out.fasta")

c(PPR_rich_out_clade_seq, PPR_rich_reference_seq) %>% 
  writeXStringSet("/path/to/PPR_rich_ref.fasta")



## 6. Perpare the iTOL annotation. To be used in the iTOL excel annotator:

# to annotate the all PPR tree with reference sequences
PPR_reference_annot <- reference_PPRs_metadata %>% mutate(color = case_when(
  Species == "Arabidopsis_thaliana-Mouse_ear_cress" ~ "#885693",
  Species == "Brassica_napus-Rape" ~ "#604F98",
  Species == "Citrus_sinensis-Sweet_orange" ~ "#E790AC",
  Species == "Coffea_canephora-Robusta_coffee" ~ "#A3C9EB",
  Species == "Glycine_max-Soybean" ~ "#E25A26",
  Species == "Gossypium_raimondii-New_world_cotton" ~ "#B4446D",
  Species == "Medicago_truncatula-Barrel_medic" ~ "#F18421",
  Species == "Populus_euphratica-Euphrates_poplar" ~ "#F2C318",
  Species == "Populus_trichocarpa-Western_balsam_poplar" ~ "#F5A71D",
  Species == "Prunus_persica-Peach" ~ "#BF1F36",
  Species == "Solanum_lycopersicum-Tomato" ~ "#0067A6",
  Species == "Vitis_vinifera-Grape" ~ "#F6937A",
  Species == "Musa_acuminata-Banana" ~ "#C2B280",
  Species == "Zea_mays-Maize" ~ "#8DB73F",
  Species == "Oryza_sativa-Japonica_rice" ~ "#008856"))

write_csv(PPR_reference_annot, "/path/to/PPR_reference_annot.csv")


# to annotate the sRNA-rich clade
PPR_rich_out_clade_annot <- PPR_rich_out_clade_meta %>% mutate(color = case_when(
  `Organism Name` == "Arabidopsis thaliana" ~ "#885693",
  `Organism Name` == "Brassica napus" ~ "#604F98",
  `Organism Name` == "Citrus sinensis" ~ "#E790AC",
  `Organism Name` == "Coffea arabica" ~ "#A3C9EB",
  `Organism Name` == "Cucumis sativus" ~ "#882F1C",
  `Organism Name` == "Glycine max" ~ "#E25A26",
  `Organism Name` == "Gossypium raimondii" ~ "#B4446D",
  `Organism Name` == "Medicago truncatula" ~ "#F18421",
  `Organism Name` == "Populus euphratica" ~ "#F2C318",
  `Organism Name` == "Populus trichocarpa" ~ "#F5A71D",
  `Organism Name` == "Prunus persica" ~ "#BF1F36",
  `Organism Name` == "Solanum lycopersicum" ~ "#0067A6",
  `Organism Name` == "Vitis vinifera" ~ "#F6937A",
  `Organism Name` == "Musa acuminata AAA Group" ~ "#C2B280",
  `Organism Name` == "Zea mays" ~ "#8DB73F",
  `Organism Name` == "Oryza sativa Japonica Group" ~ "#008856"))

write_csv(PPR_rich_out_clade_annot, "/path/to/PPR_clade_annot.csv")


PPR_rich_reference_annot <- reference_PPRs_rich_metadata %>% mutate(color = case_when(
  Species == "Arabidopsis_thaliana-Mouse_ear_cress" ~ "#885693",
  Species == "Brassica_napus-Rape" ~ "#604F98",
  Species == "Citrus_sinensis-Sweet_orange" ~ "#E790AC",
  Species == "Coffea_canephora-Robusta_coffee" ~ "#A3C9EB",
  Species == "Glycine_max-Soybean" ~ "#E25A26",
  Species == "Gossypium_raimondii-New_world_cotton" ~ "#B4446D",
  Species == "Medicago_truncatula-Barrel_medic" ~ "#F18421",
  Species == "Populus_euphratica-Euphrates_poplar" ~ "#F2C318",
  Species == "Populus_trichocarpa-Western_balsam_poplar" ~ "#F5A71D",
  Species == "Prunus_persica-Peach" ~ "#BF1F36",
  Species == "Solanum_lycopersicum-Tomato" ~ "#0067A6",
  Species == "Vitis_vinifera-Grape" ~ "#F6937A",
  Species == "Musa_acuminata-Banana" ~ "#C2B280",
  Species == "Zea_mays-Maize" ~ "#8DB73F",
  Species == "Oryza_sativa-Japonica_rice" ~ "#008856"))

write_csv(PPR_rich_reference_annot, "/path/to/PPR_rich_reference_annot.csv")


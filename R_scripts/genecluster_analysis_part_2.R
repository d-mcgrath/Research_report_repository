# load packages -----------------------------------------------------------
packs <- c('tidyverse', 'caret', 'ranger', 'RColorBrewer', 'devtools', 'reticulate', 'reshape2',
           'pheatmap', 'gplots', 'viridis', 'gridExtra', 'ggrepel', 'grid', 'formattable', 'data.table')
lapply(packs, library, character.only = TRUE)
library('DESeq2')

 # create a key to connect genes to their gene clusters --------------------
# create a dataframe which links secondary metabolites to their MAG, geneclusters, contigs, contig location, and SM class
gbk_map <- gbk_geneclusters %>%
  dplyr::rename('geneclust' = 'V1', 'mag' = 'V2', 'sm_class' = 'V3', 'ctg' = 'V4') %>%
  select(-c(V5)) %>%
  separate(mag, into = c('magdraft', 'magnum', 'ex'), sep = '_') %>%
  unite('mag', c('magdraft', 'magnum'), sep = '_') %>%
  select(-c(ex)) %>%
  select(mag, geneclust, sm_class, ctg)

# collapse secondary metabolite class factor into bins within the gbk_map dataframe
gbk_map$sm_class <- gbk_map$sm_class %>%
  fct_collapse('Polyketide (PK)' = c('t1pks', 't1pks-PUFA', 
                                     't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 'transatpks-t1pks',
                                     't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks',
                                     't2pks-t1pks', 't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 
                                     'nucleoside-otherks', 'ladderane, t2pks',
                                     't2pks-ladderane', 't2pks', 't2pks-arylpolyene', 't2pks-arylpolyene-resorcinol', 
                                     't3pks', 'phosphonate-t3pks')) %>%
  fct_collapse('NRP-PK hybrid' = c('t1pks-nrps', 'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps')) %>%
  fct_collapse('Non-ribosomal peptide (NRP)' = c('arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps', 'nrps, terpene')) %>%
  fct_collapse('Aryl polyene' = c('arylpolyene', 'arylpolyene-resorcinol', 'arylpolyene-ladderane')) %>%
  fct_collapse('Other secondary metabolite*' = c('acyl_amino_acids', 'nucleoside', 'PUFA', 'pbde', 'indole', 
                                                 'amglyccycl', 'siderophore', 'resorcinol', 'resorcinol-ladderane',
                                                 'ectoine')) %>%
  fct_collapse('Unclassified' = 'other') %>%
  fct_collapse('Ladderane' = 'ladderane') %>%
  fct_collapse('Phosphonate' = c('phosphonate', 'phosphonate-terpene', 'phosphonate-bacteriocin', 'bacteriocin, phosphonate, other',
                                 'phosphonate, other')) %>%
  fct_collapse('Terpene' = c('terpene', 'bacteriocin-lantipeptide-terpene', 'other, terpene', 'terpene, phosphonate')) %>%
  fct_collapse('Lactone' = c('butyrolactone', 'hserlactone', 'hserlactone-acyl_amino_acids')) %>%
  fct_collapse('RiPP' = c('lassopeptide', 'cyanobactin', 'thiopeptide', 'head_to_tail', 'microviridin', 
                          'sactipeptide', 'proteusin', 'lantipeptide',
                          'bacteriocin-lassopeptide', 'bacteriocin', 'bacteriocin-terpene', 
                          'proteusin-bacteriocin', 'bacteriocin-lantipeptide', 'thiopeptide-bacteriocin',
                          'bacteriocin-thiopeptide', 'bacteriocin-proteusin'))

# order the SM classes within the SM class factor
gbk_map$sm_class <- factor(gbk_map$sm_class, levels = c('Aryl polyene', 'RiPP', 'Phosphonate', 'Lactone', 'Ladderane',
                                                          'Other secondary metabolite*', 'Unclassified', 'Terpene', 
                                                          'NRP-PK hybrid', 'Non-ribosomal peptide (NRP)', 'Polyketide (PK)'))

# rename 16 geneclusters that do not have unique identifiers (16 out of 1035 geneclusters have duplicate names, 
# likely from being located on the contig of a MAG)
length(unique(gbk_map$geneclust))
head(table(sort(fct_infreq(gbk_map$geneclust))), 17)

#set geneclust col to character class to allow for manual editing then edit the duplicates, add MAG information,
   # and add phylum column for matching 
gbk_map$geneclust <- as.character(gbk_map$geneclust)
gbk_map$geneclust[c(646,647,648)] <- c('MAG_46_1A', 'MAG_46_1B', 'MAG_46_1C')
gbk_map$geneclust[c(909,910,911)] <- c('MAG_815_2A', 'MAG_815_2B', 'MAG_815_2C')
gbk_map$geneclust[c(223,224)] <- c('MAG_1305_1A', 'MAG_1305_1B')
gbk_map$geneclust[c(316,317)] <- c('MAG_1507_1A', 'MAG_1507_1B')
gbk_map$geneclust[c(326,327)] <- c('MAG_1507_12A', 'MAG_1507_12B')
gbk_map$geneclust[c(448,449)] <- c('MAG_180_3A', 'MAG_180_3B')
gbk_map$geneclust[c(466,467)] <- c('MAG_196_3A', 'MAG_196_3B')
gbk_map$geneclust[c(501,502)] <- c('MAG_216_19A', 'MAG_216_19B')
gbk_map$geneclust[c(519,520)] <- c('MAG_24_5A', 'MAG_24_5B')
gbk_map$geneclust[c(606,607)] <- c('MAG_379_1A', 'MAG_379_1B')
gbk_map$geneclust[c(771,772)] <- c('MAG_646_33A', 'MAG_646_33B')
gbk_map$geneclust[c(799,800)] <- c('MAG_699_2A', 'MAG_699_2B')
gbk_map$geneclust[c(854,855)] <- c('MAG_771_18A', 'MAG_771_18B')
gbk_map$geneclust[c(1008,1009)] <- c('MAG_957_7A', 'MAG_957_7B')

gbk_map <- gbk_map %>%
  mutate(ctg = strsplit(as.character(ctg), ';')) %>%
  unnest(ctg)

gbk_map <- gbk_map %>%
  mutate(smBGC = gbk_map$geneclust) %>%
  unite('gene', c('smBGC', 'ctg'), sep = '-')

gbk_map <- gbk_map %>%
  mutate(ctg = gene) %>%
  separate(ctg, into = c('mg', 'ctg'), sep = '-') %>%
  select(-c(mg))

gbk_gene <- gbk_gene %>%
  mutate(ctg = smBGC) %>%
  separate(ctg, into = c('ctg', 'mg'), sep = '-') %>%
  select(-c(mg)) %>%
  add_column(geneclust = gbk_map$geneclust) %>%
  add_column(phylum = NA)

# add phylum associated with each MAG into the gbk_gene dataframe
for (i in seq_along(gbk_gene$mag)) {
  for (k in seq_along(mtb2_to_mag_key$bin_id)) {
    if (as.character(gbk_gene$mag)[i] == as.character(mtb2_to_mag_key$bin_id)[k]) {
      gbk_gene$phylum[i] <- as.character(mtb2_to_mag_key$Phylum)[k]
      next
    }
  }
}

# Import combined May & Nov data for DESeq2 and Wald Tests ------------------------
prodFinal_counts <- as.matrix(read.csv('~/Documents/cariaco/R_analysis/final-gbk-read-count-table.tsv', sep = '\t', row.names = 1))
colnames(prodFinal_counts) <- sub('.counts.txt', '', colnames(prodFinal_counts))

# create the 'coldata' metadata - as a matrix - this will be metadata associated with the RNA-seq samples
coldata <- as.data.frame(colnames(prodFinal_counts))
names(coldata) <- 'sample'

coldata$season <- ifelse(str_detect(coldata$sample, 'M') == TRUE, 
                         coldata$season <- 'May', coldata$season <- 'Nov')

coldata$depth <- coldata$sample
coldata$depth <- sub('[[:upper:]]FL', '', coldata$depth)
coldata$depth <- sub('[[:upper:]]PA', '', coldata$depth)
coldata$depth <- sub('_\\d', '', coldata$depth)
coldata <- column_to_rownames(coldata, var = 'sample')

all(rownames(coldata) == colnames(cts_matrix)) # needs to be true to build DESeq2 object

coldata$season <- as.factor(coldata$season)
coldata$depth <- as.factor(coldata$depth)

#re-order the columns of counts to match the order of rows in treatment
final_coldata <- final_coldata %>%
  rownames_to_column(var = 'rnms') %>%
  slice(-29) %>%
  column_to_rownames(var = 'rnms')

prodFinal_counts <- prodFinal_counts[, rownames(final_coldata)] # re-order the counts matrix columns to match rowname 
  # order from the coldata
all(colnames(prodFinal_counts) == rownames(final_coldata))


# define coldata by O2 concentration, not by depth in meters - this allows for comparison across seasons, and layer-based
 # comparison (oxycline, shallow anoxic, euxinic layers) vs. depth-wise - gives a big picture look at the data
oxy_coldata <- prodFinal_coldata
oxy_coldata$depth <- as.numeric(as.character(oxy_coldata$depth))
oxy_coldata$layer <- ifelse(oxy_coldata$depth <= 237, oxy_coldata$layer <- 'oxycline',
                            ifelse(oxy_coldata$depth > 237 & oxy_coldata$depth < 900, 
                                   oxy_coldata$layer <- 'shallow_anoxic', oxy_coldata$layer <- 'euxinic'))

# convert depth and layer columns to factor class
oxy_coldata$depth <- as.factor(oxy_coldata$depth)
oxy_coldata$layer <- as.factor(oxy_coldata$layer)

# the DESeq2 model will be run with the fraction (FL vs PA) and layer (water layer) factors combined in unique
 # pairs that will give the same result as if they were left uncombined, but is easier for pulling out DESeq2 Wald Test results
 # (as suggested by Micheal Love - one of the creators of DESeq2)
all(colnames(prodFinal_counts) == rownames(oxy_coldata))
final_coldata <- oxy_coldata[c(3, 4)]
final_coldata <- final_coldata %>%
  unite(col = 'group', c('fraction', 'layer'), sep = '_')

final_coldata$group <- factor(final_coldata$group)
levels(final_coldata$group)

prodFinal_counts <- prodFinal_counts %>%
  data.frame() %>%
  rownames_to_column(var= 'gene')

# subsetted dataframe of all the counts for all predicted genes to include just the counts of biosynthetic gene clusters
prodFinal_counts <- subset(prodFinal_counts, !(gene %in% final_gene$smBGC))
rownames(prodFinal_counts) <- c()
prodFinal_counts <- prodFinal_counts %>%
  column_to_rownames(var = 'gene') %>%
  as.matrix()
  
#load the DESeq2 object
prodFinal_dds <- DESeqDataSetFromMatrix(countData = prodFinal_counts, colData = final_coldata,
                                    design = ~ group)
prod_vst <- vst(prodFinal_dds, blind = FALSE) # take vst normalized count data
prod_vst_mat <- assay(prod_vst)
prod_vst_cor <- cor(prod_vst_mat)
pheatmap(prod_vst_cor) # visualize correlations of the normalized count data, look at how the samples are different or similar

fin_dds <- estimateSizeFactors(prodFinal_dds)
normProd_counts <- counts(fin_dds, normalized = TRUE) # save the normalized count data

gbkFinal_dds <- DESeq(prodFinal_dds) # run the DESeq2 statistical model which assumes a negative binomial distribution for
    # the RNA-seq count matrix data
resultsNames(gbkFinal_dds)

# pull Wald test results comparing the oxycline PA to FL fractions 
gbk_oxyRes <- results(gbkFinal_dds, contrast = c('group',
                                          'PA_oxycline', 'FL_oxycline'), alpha = 0.05)
summary(gbk_oxyRes)

# pull Wald test results comparing the euxinic PA to FL fractions 
gbk_euxRes <- results(gbkFinal_dds, contrast = c('group',
                                          'PA_euxinic', 'FL_euxinic'), alpha = 0.05)
summary(gbk_euxRes)

# pull Wald test results comparing the shallow anoxic PA to FL fractions 
gbk_anoxRes <- results(gbkFinal_dds, contrast = c('group', 'PA_shallow_anoxic',
                                           'FL_shallow_anoxic'), alpha = 0.05)
summary(gbk_anoxRes)

# make a copy of the coldata
copy_fc <- final_coldata %>%
  rownames_to_column(var = 'sample')

# create another key connecting MAG names to taxonomy
mtb2_to_mag_key <- taxonomy_draftMags[, c(1, 3)]
mtb2_to_mag_key$Phylum <- ifelse(mtb2_to_mag_key$Phylum == 'p__', 
                                 sub('p__', 'Unclassified', mtb2_to_mag_key$Phylum), 
                                 sub('p__', "", mtb2_to_mag_key$Phylum))
mtb2_to_mag_key$bin_id <- sub('CarAnox_mtb2-', 'MAG_', mtb2_to_mag_key$bin_id)

mtb2_to_mag_key <- mtb2_to_mag_key %>%
  add_column(mag = mtb2_to_mag_key$bin_id) %>%
  add_column(phy = mtb2_to_mag_key$Phylum)

mtb2_to_mag_key <- mtb2_to_mag_key %>%
  separate(mag, into = c('mg', 'magnum'), sep = '_') %>%
  unite(taxmag, c('phy', 'magnum'), sep = '_') %>%
  select(-c(mg))

# function to subset the normalized RNA-seq counts matrix for significant genes from a given water layer
get_normCounts <- function(results, columns, overall_normCounts) {
  res_tb <- results %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as_tibble()
  sig <- res_tb %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 2)
  norm_sig <- overall_normCounts[, columns] %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    filter(gene %in% sig$gene) %>%
    filter(gene %in% gbk_gene$smBGC) %>%
    column_to_rownames(var = 'gene') %>%
    as.matrix()
  return(norm_sig)
}

# function to pull out geneclusters and their associated taxonomies
de_bar <- function(layer_particle, title, size = 9) {
  layer_particle <- layer_particle %>%
    rownames_to_column(var = 'gene')
  layer_particle$sm_type <- NA
  layer_particle$geneclust <- NA
  
  for (i in seq_along(layer_particle$gene)) {
    for (k in seq_along(gbk_gene$smBGC)) {
      if (as.character(layer_particle$gene)[i] == as.character(gbk_gene$smBGC)[k]) {
        layer_particle$sm_type[i] <- as.character(gbk_gene$sm_class)[k]
        layer_particle$geneclust[i] <- as.character(gbk_gene$geneclust)[k]
        next
      }
    }
  }
  
  layer_particle <- layer_particle %>% select(geneclust, sm_type, gene, position)
  layer_particle <- aggregate(. ~ layer_particle$geneclust, layer_particle, function(x) toString(unique(x)))
  layer_particle <- layer_particle %>% select(-`layer_particle$geneclust`)
  layer_particle <- layer_particle %>% add_column(mag = layer_particle$geneclust, .after = 2)
  layer_particle <- layer_particle %>%
    separate(mag, into = c('mag', 'magnum', 'cluster'), sep = '_') %>%
    unite(col = mag, c('mag', 'magnum'), sep = '_') %>%
    select(-c(cluster)) %>%
    add_column(taxonomy = NA, .after = 3)
  
  for (i in seq_along(layer_particle$mag)) {
    for (k in seq_along(mtb2_to_mag_key$bin_id)) {
      if (as.character(layer_particle$mag)[i] == as.character(mtb2_to_mag_key$bin_id)[k]) {
        layer_particle$taxonomy[i] <- as.character(mtb2_to_mag_key$Phylum)[k]
        next
      }
    }
  }
  
  layer_particle$sm_type <- layer_particle$sm_type %>%
    fct_collapse('Polyketide (PK)' = c('t1pks', 't1pks-PUFA', 
                                       't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 'transatpks-t1pks',
                                       't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks',
                                       't2pks-t1pks', 't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 
                                       'nucleoside-otherks', 'ladderane, t2pks',
                                       't2pks-ladderane', 't2pks', 't2pks-arylpolyene', 't2pks-arylpolyene-resorcinol', 
                                       't3pks', 'phosphonate-t3pks')) %>%
    fct_collapse('NRP-PK hybrid' = c('t1pks-nrps', 'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps')) %>%
    fct_collapse('Non-ribosomal peptide (NRP)' = c('arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps', 'nrps, terpene')) %>%
    fct_collapse('Aryl polyene' = c('arylpolyene', 'arylpolyene-resorcinol', 'arylpolyene-ladderane')) %>%
    fct_collapse('Other secondary metabolite*' = c('acyl_amino_acids', 'nucleoside', 'PUFA', 'pbde', 'indole', 
                                                   'amglyccycl', 'siderophore', 'resorcinol', 'resorcinol-ladderane',
                                                   'ectoine')) %>%
    fct_collapse('Unclassified' = 'other') %>%
    fct_collapse('Ladderane' = 'ladderane') %>%
    fct_collapse('Phosphonate' = c('phosphonate', 'phosphonate-terpene', 'phosphonate-bacteriocin', 'bacteriocin, phosphonate, other',
                                   'phosphonate, other')) %>%
    fct_collapse('Terpene' = c('terpene', 'bacteriocin-lantipeptide-terpene', 'other, terpene', 'terpene, phosphonate')) %>%
    fct_collapse('Lactone' = c('butyrolactone', 'hserlactone', 'hserlactone-acyl_amino_acids')) %>%
    fct_collapse('RiPP' = c('lassopeptide', 'cyanobactin', 'thiopeptide', 'head_to_tail', 'microviridin', 
                            'sactipeptide', 'proteusin', 'lantipeptide',
                            'bacteriocin-lassopeptide', 'bacteriocin', 'bacteriocin-terpene', 
                            'proteusin-bacteriocin', 'bacteriocin-lantipeptide', 'thiopeptide-bacteriocin',
                            'bacteriocin-thiopeptide', 'bacteriocin-proteusin'))
  layer_particle <- within(layer_particle, taxonomy <- factor(taxonomy, levels = names(sort(table(taxonomy), decreasing = TRUE))))
  
  layer_particle$sm_type <- factor(layer_particle$sm_type, levels = c('Aryl polyene', 'RiPP', 'Phosphonate', 'Lactone', 'Ladderane',
                                                                      'Other secondary metabolite*', 'Unclassified', 'Terpene', 
                                                                      'NRP-PK hybrid', 'Non-ribosomal peptide (NRP)', 'Polyketide (PK)'))
  
  layer_palette <- c("khaki", "#2A7FB7", "#99CD91", "#52AF43", "#B89B74", "#ED4F50", 
                     "#A6CEE3", "#FDA440", "#FFCCCB", "#B294C7", "#825D99")
  names(layer_palette) <- levels(layer_particle$sm_type) 
  
  # plot the geneclusters by taxa, with each SM class color coded within a barplot
  smbgc_plot <- 
    ggplot(arrange(layer_particle, sm_type), aes(taxonomy, fill = sm_type)) +
    geom_bar() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1),
          axis.title = (element_text(size = 10)),
          legend.position = c(0.96, 0.96), legend.justification = c(0.96, 0.96),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.background = element_blank(), panel.border = 
            element_rect(color = 'black', fill = NA),
          legend.box.background = element_rect(color = 'black'),
          plot.title = element_text(size = 10)) +
    labs(x = 'Phylum of MAG expressing cluster', fill = 'Secondary metabolite class', y = 'Number of biosynthetic gene clusters', title = title) +
    scale_fill_manual(values = layer_palette)
  
  output <- list(layer_particle, smbgc_plot)
  return(output)
}

# function to extract annotation columns
get_colAnno <- function(norm_layerSig) {
  colAnno <- as.data.frame(colnames(norm_layerSig))
  names(colAnno)[1] <- 'rnms'
  colAnno <- colAnno %>%
    add_column('Fraction' = NA) %>%
    column_to_rownames(var = 'rnms')
  colAnno$Fraction <- ifelse(str_detect(rownames(colAnno), 'PA'), 'PA', 'FL')
  return(colAnno)
}

#function to extraction annotation rows
get_rowAnno <- function(norm_layerSig) {
  rowAnno <- gbk_gene %>%
    filter(smBGC %in% rownames(norm_layerSig)) %>%
    select(-c(mag)) %>%
    column_to_rownames(var = 'smBGC')
  return(rowAnno)
}

# create heatmaps, pull out signifcant genes ------------------------------

# pull out significantly expressed (padj < 0.05 & absolute value of log2FoldChange > 2) for oxycline data
gbk_oxy_SigGenes <- gbk_oxyRes %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

# pull out significantly expressed (padj < 0.05 & absolute value of log2FoldChange > 2) for shallow anoxic data
gbk_anox_SigGenes <- gbk_anoxRes %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

# pull out significantly expressed (padj < 0.05 & absolute value of log2FoldChange > 2) for euxinic data
gbk_eux_SigGenes <- gbk_euxRes %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

# run function to subset counts matrix to pull out count data for significantly expressed genes in each water layer
gbk_norm_oxySig <- get_normCounts(gbk_oxyRes, c(1:5, 12:17, 24:28, 36:40), normProd_counts)
gbk_norm_anoxSig <- get_normCounts(gbk_anoxRes, c(6:9, 18:21, 29:32, 41:44), normProd_counts)
gbk_norm_euxSig <- get_normCounts(gbk_euxRes, c(10, 11, 22, 23, 33, 34, 45, 46), normProd_counts)

# create column annotation dataframes for heatmap creation
gbk_oxy_colAnno <- get_colAnno(gbk_norm_oxySig)
gbk_anox_colAnno <- get_colAnno(gbk_norm_anoxSig)
gbk_eux_colAnno <- get_colAnno(gbk_norm_euxSig)

# create row annotation dataframes for heatmap creation
gbk_oxy_rowAnno <- get_rowAnno(gbk_norm_oxySig)
gbk_anox_rowAnno <- get_rowAnno(gbk_norm_anoxSig)
gbk_eux_rowAnno <- get_rowAnno(gbk_norm_euxSig)

# rename sample names to make them more readable
copy_gbkNormOxySig <- gbk_norm_oxySig
copy_gbkOxyColAnno <- gbk_oxy_colAnno

# rename columns to make them more readable
colnames(copy_gbkNormOxySig) <- c("May_FL_103_2", "May_FL_198_1", "May_FL_198_2", "May_FL_234_1", "May_FL_234_2", "May_PA_103_1", 
                                  "May_PA_103_2", "May_PA_198_1", "May_PA_198_2", "May_PA_234_1", "May_PA_234_2",
                                  "Nov_FL_148_1", "Nov_FL_148_2", "Nov_FL_200_1", "Nov_FL_200_2", "Nov_FL_237_1",
                                  "Nov_PA_200_1", "Nov_PA_200_2", "Nov_PA_237_1", "Nov_PA_237_2")
                              
# apply the naming scheme to the column annotation
rownames(copy_gbkOxyColAnno) <- colnames(copy_gbkNormOxySig)

#create a heatmap of the normalized counts of significantly expressed genes in the oxycline samples
gbk_oxyHeat <- pheatmap(copy_gbkNormOxySig, scale = 'row', show_rownames = F,
         clustering_method = 'ward.D', annotation_col = copy_gbkOxyColAnno,
         annotation_colors = relAbunPalette, treeheight_row = 25, treeheight_col = 0)

# create a dataframe which shows which genes belong to which cluster (4 clusters total)
gbk_norm_oxySig_test <- rownames(gbk_norm_oxySig[gbk_oxyHeat$tree_row[['order']], ]) %>%
  data.frame() %>%
  add_column(clust = cutree(gbk_oxyHeat$tree_row, 4)[gbk_oxyHeat$tree_row[['order']]])
gbk_oxy_smBGC <- cutree(gbk_oxyHeat$tree_row, 4)

# dataframes for genes belonging to PA fraction and FL fraction clusters from the heatmap
gbk_oxyFL <- as.data.frame(which(gbk_oxy_smBGC == 1)); colnames(gbk_oxyFL) <- 'position'
gbk_oxyPA <- as.data.frame(which(gbk_oxy_smBGC == 2 | gbk_oxy_smBGC == 3 |
                                   gbk_oxy_smBGC == 4)); colnames(gbk_oxyPA) <- 'position'

# create barplots for the genes which represent the SM geneclusters they belong to (for each of their rep phylums)
gbk_oxyFL <- de_bar(gbk_oxyFL, 'Differentially expressed smBGCs in the oxycline FL fraction')
gbk_oxyPA <- de_bar(gbk_oxyPA, 'Differentially expressed smBGCs in the oxycline PA fraction')


#rername the sample names to make them more readable
copy_gbkNormAnoxSig <- gbk_norm_anoxSig
copy_gbkAnoxColAnno <- gbk_anox_colAnno

colnames(copy_gbkNormAnoxSig) <- c("May_FL_295_1", "May_FL_95_2", "May_FL_314_1", "May_FL_314_2", "May_PA_295_1", 
                                   "May_PA_295_2", "May_PA_314_1", "May_PA_314_2", "Nov_FL_247_1", "Nov_FL_247_2", 
                                   "Nov_FL_267_1", "Nov_FL_267_2", "Nov_PA_247_1", "Nov_PA_247_2", "Nov_PA_267_1", 
                                   "Nov_PA_267_2")
rownames(copy_gbkAnoxColAnno) <- colnames(copy_gbkNormAnoxSig)

#create a heatmap of the normalized counts of significantly expressed genes in the shallow anoxic samples
gbk_anoxHeat <- pheatmap(copy_gbkNormAnoxSig, scale = 'row', show_rownames = FALSE,
         clustering_method = 'ward.D', annotation_col = copy_gbkAnoxColAnno,
         annotation_colors = relAbunPalette, treeheight_row = 25, treeheight_col = 0)

# create a dataframe which shows which genes belong to which cluster (2 clusters total)
gbk_norm_anoxSig_test <- rownames(gbk_norm_anoxSig[gbk_anoxHeat$tree_row[['order']], ]) %>%
  data.frame() %>%
  add_column(clust = cutree(gbk_anoxHeat$tree_row, 2)[gbk_anoxHeat$tree_row[['order']]])

gbk_anox_smBGC <- cutree(gbk_anoxHeat$tree_row, 2)

# dataframes for genes belonging to PA fraction and FL fraction clusters from the heatmap
gbk_anoxFL <- as.data.frame(which(gbk_anox_smBGC == 2)); colnames(gbk_anoxFL) <- 'position'
gbk_anoxPA <- as.data.frame(which(gbk_anox_smBGC == 1)); colnames(gbk_anoxPA) <- 'position'

# create barplots for the genes which represent the SM geneclusters they belong to (for each of their rep phylums)
gbk_anoxFL <- de_bar(gbk_anoxFL, 'Differentially expressed smBGCs in the shallow anoxic FL fraction')
gbk_anoxPA <- de_bar(gbk_anoxPA, 'Differentially expressed smBGCs in the shallow anoxic PA fraction')


#rename samples to make them more readable
copy_gbkNormEuxSig <- gbk_norm_euxSig
copy_gbkEuxColAnno <- gbk_eux_colAnno

colnames(copy_gbkNormEuxSig) <- c("May_FL_900_1", "May_FL_900_2", "May_PA_900_1", "May_PA_900_2", "Nov_FL_900_1",
                                  "Nov_FL_900_2", "Nov_PA_900_1", "Nov_PA_900_2")
rownames(copy_gbkEuxColAnno) <- colnames(copy_gbkNormEuxSig)

#create the heatmap for euxinic depths
gbk_euxHeat <- pheatmap(copy_gbkNormEuxSig, scale = 'row',
         clustering_method = 'ward.D', annotation_col = copy_gbkEuxColAnno,
         annotation_colors = relAbunPalette, show_rownames = FALSE,
         treeheight_row = 25, treeheight_col = 0)

# create a dataframe which shows which genes belong to which cluster (4 clusters total)
gbk_norm_euxSig_test <- rownames(gbk_norm_euxSig[gbk_euxHeat$tree_row[['order']], ]) %>%
  data.frame() %>%
  add_column(clust = cutree(gbk_euxHeat$tree_row, 4)[gbk_euxHeat$tree_row[['order']]])

gbk_eux_smBGC <- cutree(gbk_euxHeat$tree_row, 4)

# dataframes for genes belonging to PA fraction and FL fraction clusters from the heatmap
gbk_euxFL <- as.data.frame(which(gbk_eux_smBGC == 4)); colnames(gbk_euxFL) <- 'position'
gbk_euxPA <- as.data.frame(which(gbk_eux_smBGC == 1 | gbk_eux_smBGC == 2 |
                                   gbk_eux_smBGC == 3)); colnames(gbk_euxPA) <- 'position'

# create barplots for the genes which represent the SM geneclusters they belong to (for each of their rep phylums)
gbk_euxFL <- de_bar(gbk_euxFL, 'Differentially expressed smBGCs in the euxinic FL fraction')
gbk_euxPA <- de_bar(gbk_euxPA, 'Differentially expressed smBGCs in the euxinic PA fraction')

# save heatmaps arranged with PA and FL barplots for the oxycline samples
png(filename = 'oxycline_de_FINAL.png', units = 'in', width = 11, height = 16.3, res = 500)
grid.arrange(gbk_oxyHeat[[4]], gbk_oxyPA[[2]], gbk_oxyFL[[2]],
             layout_matrix = rbind(c(1,1,1,2,2), c(1,1,1,3,3)))
dev.off()

# save heatmaps arranged with PA and FL barplots for the shallow anoxic samples
png(filename = 'anox_de_FINAL.png', units = 'in', width = 12, height = 16.3, res = 500)
grid.arrange(gbk_anoxHeat[[4]], gbk_anoxFL[[2]], gbk_anoxPA[[2]],
             layout_matrix = rbind(c(1,1,1,2,2), c(1,1,1,3,3)))
dev.off()

# save heatmaps arranged with PA and FL barplots for the euxinic samples
png(filename = 'eux_de_FINAL.png', units = 'in', width = 11, height = 16.3, res = 500)
grid.arrange(gbk_euxHeat[[4]], gbk_euxFL[[2]], gbk_euxPA[[2]],
             layout_matrix = rbind(c(1,1,1,2,2), c(1,1,1,3,3)))
dev.off()

# function to generate inputs for volanco plot function
getVolcInput <- function(layer_results, copy_gbkNorm_layerSig) { #VolcInput = DESeq2 results for a given water layer
  VolcInput <- layer_results
  VolcInput <- VolcInput %>%
    as.data.frame() %>%
    mutate(gene = rownames(VolcInput)) %>%
    as_tibble() %>%
    filter(gene %in% rownames(copy_gbkNorm_layerSig))
  clusterCount <- VolcInput %>%
    select(gene, log2FoldChange) %>%
    add_column(geneclust = NA)
  
  for (i in seq_along(clusterCount$gene)) {
    for (k in seq_along(gbk_gene$smBGC)) {
      if (as.character(clusterCount$gene)[i] == as.character(gbk_gene$smBGC)[k]) {
        clusterCount$geneclust[i] <- as.character(gbk_gene$geneclust)[k]
        next
      }
    }
  }
  
  clusterCount$log2FoldChange <- ifelse(clusterCount$log2FoldChange > 0, 'up', 'down')
  clusterCount <- aggregate(. ~ geneclust, clusterCount, function(x) toString(unique(x)))
  
  ## Obtain logical vector where TRUE values denote 
  ## padj values < 0.05 and fold change > 1.5 in either direction
  VolcInput <- VolcInput %>% 
    #mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)
    mutate(threshold = padj < 0.05)
  VolcInput$threshold[is.na(VolcInput$threshold)] <- FALSE
  
  ## Create a column to indicate which genes to label
  VolcInput <- VolcInput %>% arrange(padj)
  VolcInput <- VolcInput %>%
    add_column(sm_type = "")
  
  for (i in seq_along(VolcInput$gene[VolcInput$threshold == TRUE])) {
    for (k in seq_along(gbk_gene$smBGC)) {
      if (as.character(VolcInput$gene)[i] == as.character(gbk_gene$smBGC)[k]) {
        VolcInput$sm_type[i] <- as.character(gbk_gene$sm_class)[k]
        next
      }
    }
  }
  
  VolcInput$sm_type <- VolcInput$sm_type %>%
    fct_collapse('Polyketide (PK)' = c('t1pks', 't1pks-PUFA', 
                                       't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 'transatpks-t1pks',
                                       't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks',
                                       't2pks-t1pks', 't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 
                                       'nucleoside-otherks', 'ladderane, t2pks',
                                       't2pks-ladderane', 't2pks', 't2pks-arylpolyene', 't2pks-arylpolyene-resorcinol', 
                                       't3pks', 'phosphonate-t3pks')) %>%
    fct_collapse('NRP-PK hybrid' = c('t1pks-nrps', 'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps')) %>%
    fct_collapse('Non-ribosomal peptide (NRP)' = c('arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps', 'nrps, terpene')) %>%
    fct_collapse('Aryl polyene' = c('arylpolyene', 'arylpolyene-resorcinol', 'arylpolyene-ladderane')) %>%
    fct_collapse('Other secondary metabolite*' = c('acyl_amino_acids', 'nucleoside', 'PUFA', 'pbde', 'indole', 
                                                   'amglyccycl', 'siderophore', 'resorcinol', 'resorcinol-ladderane',
                                                   'ectoine')) %>%
    fct_collapse('Unclassified' = 'other') %>%
    fct_collapse('Ladderane' = 'ladderane') %>%
    fct_collapse('Phosphonate' = c('phosphonate', 'phosphonate-terpene', 'phosphonate-bacteriocin', 'bacteriocin, phosphonate, other',
                                   'phosphonate, other')) %>%
    fct_collapse('Terpene' = c('terpene', 'bacteriocin-lantipeptide-terpene', 'other, terpene', 'terpene, phosphonate')) %>%
    fct_collapse('Lactone' = c('butyrolactone', 'hserlactone', 'hserlactone-acyl_amino_acids')) %>%
    fct_collapse('RiPP' = c('lassopeptide', 'cyanobactin', 'thiopeptide', 'head_to_tail', 'microviridin', 
                            'sactipeptide', 'proteusin', 'lantipeptide',
                            'bacteriocin-lassopeptide', 'bacteriocin', 'bacteriocin-terpene', 
                            'proteusin-bacteriocin', 'bacteriocin-lantipeptide', 'thiopeptide-bacteriocin',
                            'bacteriocin-thiopeptide', 'bacteriocin-proteusin'))
  
  VolcInput$sm_type <- factor(VolcInput$sm_type, levels = c('Aryl polyene', 'RiPP', 'Phosphonate', 'Lactone', 'Ladderane',
                                                            'Other secondary metabolite*', 'Unclassified', 'Terpene', 
                                                            'NRP-PK hybrid', 'Non-ribosomal peptide (NRP)', 'Polyketide (PK)'))
  VolcInput <- VolcInput %>% add_column('gc_match_final' = VolcInput$gene)
  VolcInput <- VolcInput %>%
    separate(gc_match_final, into = c('contig', 'mag_final'), sep = '-') %>%
    select(-c(contig))
  
  for (i in seq_along(VolcInput$mag_final[VolcInput$threshold == TRUE])) {
    for (k in seq_along(mp_key$mag_id)) {
      if (as.character(VolcInput$mag_final)[i] == as.character(mp_key$mag_id)[k]) {
        VolcInput$mag_final[i] <- as.character(mp_key$full_phylum)[k]
        next
      }
    }
  }
  
  VolcInput$mag_final <- paste0(VolcInput$mag_final, '_', VolcInput$sm_type)
  VolcInput$mag_final[VolcInput$threshold == FALSE] <- ''
  VolcInput$sm_type <- NULL
  
  VolcInput <- VolcInput %>%
    add_column(sm = VolcInput$mag_final) %>%
    separate(sm, into = c('tax', 'tax_num', 'sm_class'), sep = '_') %>%
    unite(genelabels, c(tax, tax_num), sep = '_')
  
  VolcInput$genelabels[VolcInput$threshold == FALSE] <- ''
  VolcInput$genelabels_final <- ''
  VolcInput$genelabels_final[1:20] <- VolcInput$genelabels[1:20]
  
  VolcInput$sm_class <- factor(VolcInput$sm_class, levels = c('Aryl polyene', 'RiPP', 'Phosphonate', 'Lactone', 'Ladderane',
                                                              'Other secondary metabolite*', 'Unclassified', 'Terpene', 
                                                              'NRP-PK hybrid', 'Non-ribosomal peptide (NRP)', 'Polyketide (PK)'))
  VolcInput$reg <- ''
  VolcInput$reg[VolcInput$threshold == TRUE] <- ifelse(VolcInput$log2FoldChange > 0, 'up', 'down')
  
  output <- list(VolcInput, clusterCount)
  return(output)
}

# volcano plot function
getVolcano <- function(res_tb, layer, nudgex = -5, nudgey = 3, tag = NULL) {
  
  sm_palette <- c("khaki", "#2A7FB7", "#99CD91", "#52AF43", "#B89B74", "#ED4F50", 
                  "#A6CEE3", "#FDA440", "#FFCCCB", "#B294C7", "#825D99")
  names(sm_palette) <- levels(res_tb$sm_class)
  
  ggplot(res_tb, aes(x = log2FoldChange, y = -log10(padj), color = sm_class)) +
    theme_bw() +
    geom_text_repel(label = res_tb$genelabels_final, segment.size = 0.2, segment.color = 'grey50', direction = 'y',
                    nudge_x = nudgex, nudge_y = nudgey, show.legend = F, size = 3) +
    labs(x = 'Log2 fold change', y = '-Log10 adjusted p-value', title = paste(layer, 'differential expression across fraction (PA vs FL)'),
         color = 'Biosynthetic gene cluster class', tag = tag) +
    theme(plot.title = element_text(size = 9),
          axis.title = element_text(size = 9),
          #legend.position = c(0.04, 0.92), legend.justification = c(0.04, 0.92),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.background = element_blank(), panel.border = 
            element_rect(color = 'black', fill = NA),
          legend.box.background = element_rect(color = 'black')) +
    scale_color_manual(values = sm_palette, na.translate = FALSE) +
    geom_point() +
    xlim(-10, 10) +
    ylim(0, 30)
}

# generate input for volcanco plots for each water layer
oxy_volcInput <- getVolcInput(gbk_oxyRes, copy_gbkNormOxySig)
oxy_volc <- getVolcano(oxy_volcInput[[1]], 'Oxycline', -5, 6, tag = '7a')

anox_volcInput <- getVolcInput(gbk_anoxRes, copy_gbkNormAnoxSig)
anox_volc <- getVolcano(anox_volcInput[[1]], 'Shallow anoxic', tag = '7b')

eux_volcInput <- getVolcInput(gbk_euxRes, copy_gbkNormEuxSig)
eux_volc <- getVolcano(eux_volcInput[[1]], 'Euxinic', -5, 6, tag = '7c')

png(filename = 'volcano_plots_FINAL_final.png', units = 'in', width = 8, height = 11, res = 500)
grid.arrange(oxy_volc, anox_volc, eux_volc,
             layout_matrix = rbind(c(1,1,1), c(2,2,2), c(3,3,3)))
dev.off()

# create supplementary table 1 with sample naming scheme information --------
suppTab1 <- data.table('Sample Name' = rownames(coldata), 'Time point' = c(rep('May', 24), rep('November', 24)),
                       'Year' = 2014, 'Sampling depth' = paste(coldata$depth, 'm', sep = ''),
                       'Sample fraction' = c(rep('Free-living', 12), rep('Particle-associated', 12),
                                             rep('Free-living', 12), rep('Particle-associated', 12)),
                       'Replicate' = rep(c(1,2), 24))

suppTab1 <- suppTab1 %>%
  add_column(samp = suppTab1$`Sampling depth`, .after = 4)
suppTab1 <- suppTab1 %>%
  mutate(samp = sub('m', '', suppTab1$samp)) %>%
  mutate(samp = as.numeric(samp)) %>%
  mutate(`Water layer` = ifelse(suppTab1$samp < 238, 'Oxycline', ifelse(suppTab1$samp > 237 & suppTab1$samp < 900,
                                                                'Shallow anoxic', 'Euxinic'))) %>%
  select(-c(samp))

suppTab1$`Sample Name` <- c("May_FL_103_1", "May_FL_103_2", "May_FL_198_1", "May_FL_198_2", "May_FL_234_1",
                            "May_FL_234_2", "May_FL_295_1", "May_FL_295_2", "May_FL_314_1", "May_FL_314_2", 
                            "May_FL_900_1", "May_FL_900_2", "May_PA_103_1", "May_PA_103_2", "May_PA_198_1", 
                            "May_PA_198_2", "May_PA_234_1", "May_PA_234_2", "May_PA_295_1", "May_PA_295_2", 
                            "May_PA_314_1", "May_PA_314_2", "May_PA_900_1", "May_PA_900_2", "Nov_FL_148_1", 
                            "Nov_FL_148_2", "Nov_FL_200_1", "Nov_FL_200_2", "Nov_FL_237_1", "Nov_FL_237_2",
                            "Nov_FL_247_1", "Nov_FL_247_2", "Nov_FL_267_1", "Nov_FL_267_2", "Nov_FL_900_1",
                            "Nov_FL_900_2", "Nov_PA_148_1", "Nov_PA_148_2", "Nov_PA_200_1", "Nov_PA_200_2",
                            "Nov_PA_237_1", "Nov_PA_237_2", "Nov_PA_247_1", "Nov_PA_247_2", "Nov_PA_267_1",
                            "Nov_PA_267_2", "Nov_PA_900_1", "Nov_PA_900_2")

write.table(suppTab1, file = '~/Documents/cariaco/figures/supp-table-1-sample-info.txt', quote = FALSE,
            row.names = FALSE, sep = '\t')

#suppTab1 <- formattable(suppTab1)
#function to save supp table 1
#library("htmltools")
#library("webshot")    
#export_formattable <- function(f, file, width = "100%", height = NULL, 
#                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
# save supp table 1
#export_formattable(suppTab1, 'supp_table_1.png')

df_endOxyRes <- as.data.frame(end_oxyRes)
df_endOxyRes <- df_endOxyRes %>%
  rownames_to_column(var = 'gene') %>%
  filter(padj < 0.05)


get_normCounts <- function(results, columns, overall_normCounts) {
  res_tb <- results %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as_tibble()
  sig <- res_tb %>%
    filter(padj < 0.05)
  norm_sig <- overall_normCounts[, columns] %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    filter(gene %in% sig$gene) %>%
    filter(str_detect(gene, 'cluster')) %>%
    column_to_rownames(var = 'gene') %>%
    as.matrix()
  return(norm_sig)
}


# set wd to location of the annotation files
mypath = "~/Documents/cariaco/draftMag-gene-annotations/"
setwd(mypath)

# index vector of the .tsv file names
fileIndex <- list.files(path = mypath, pattern = "*.tsv")

# import all genecluster text files
geneAnno <- lapply(fileIndex, function(x) {read.csv(x, sep = '\t', fill = TRUE,
            col.names = c('protein_acc', 'MD5_digest', 'seq_len', 'database', 'sig_acc', 'sig_descrip',
                          'start_loc', 'stop_loc', 'score', 'status', 'date', 'interpro_anno', 'interpro_descrip'), 
            na.strings = c('', '-'), header = FALSE)})

# combine them
geneAnno <- do.call("rbind", lapply(geneAnno, as_tibble))

mags_with_smBGCs <- tibble(as.character(unique(gbk_gene$mag))) %>%
  rename('mag' = 1)

#
geneAnno <- geneAnno %>%
  select(-c(MD5_digest, status, date)) %>%
  add_column(dmag = geneAnno$protein_acc, .after = 1) %>%
  separate(dmag, into = c('ctg', 'mag'), sep = '-') %>%
  select(-c('ctg')) %>%
  filter(mag %in% mags_with_smBGCs$mag) %>%
  rename('gene' = 'protein_acc')

# a subset of the all unique gene annotation discriptors
sig_descrip_uniq <- tibble(unique(geneAnno$sig_descrip)) %>%
  rename('gene_description' = 1)

# a subset of the all unique interpro gene annotation discriptors
intropro_descrip_uniq <- tibble(unique(geneAnno$interpro_descrip)) %>%
  rename('gene_description' = 1)


# subsetting the interproscan gene annotations for amo and rubis --------
rubis_amo_genes <- geneAnno %>%
  filter(sig_descrip == 'RuBisCO' | sig_descrip == 'RuBisCO_small' | sig_descrip == 'RuBisCO_large_II' | 
           sig_descrip == 'rubisco_III: ribulose bisphosphate carboxylase, type III' |
           sig_descrip == 'RuBisCO_IV_RLP' | sig_descrip == 'RuBisCO small subunit signature' | 
           interpro_descrip == 'RuBisCO large subunit, N-terminal domain superfamily' |  interpro_descrip == 'RuBisCO' |
           interpro_descrip == 'Ammonia monooxygenase/particulate methane monooxygenase, subunit C' |
           interpro_descrip == 'Ammonia monooxygenase/particulate methane monooxygenase, subunit C domain superfamily' |
           interpro_descrip == 'Ammonia/methane monooxygenase, subunitB, N-terminal' |
           interpro_descrip == 'Ammonia/methane monooxygenase, subunit B, C-terminal' |
           interpro_descrip == 'Ammonia/methane monooxygenase, subunit B, hairpin domain superfamily' |
           interpro_descrip == 'Ammonia monooxygenase/particulate methane monooxygenase, subunit B' |
           interpro_descrip == 'Ammonia monooxygenase/particulate methane monooxygenase, subunit A' |
           interpro_descrip == 'Ammonia/particulate methane monooxygenase, subunit A superfamily' |
           sig_descrip == 'Ammonia monooxygenase' | 
           sig_descrip == 'Ammonia monooxygenase/methane monooxygenase, subunit C' | 
           sig_descrip == 'CH4_NH3mon_ox_C: methane monooxygenase/ammonia monooxygenase, subunit C' | 
           sig_descrip == 'CH4_NH3mon_ox_B: methane monooxygenase/ammonia monooxygenase, subunit B' | 
           sig_descrip == 'CH4_NH3mon_ox_A: methane monooxygenase/ammonia monooxygenase, subunit A')
rubis_amo_genes$sig_descrip <- as.character(rubis_amo_genes$sig_descrip); rubis_amo_genes$interpro_descrip <- as.character(rubis_amo_genes$interpro_descrip)
rubis_amo_genes$interpro_descrip[is.na(rubis_amo_genes$interpro_descrip)] <- 
  as.character(rubis_amo_genes$sig_descrip[is.na(rubis_amo_genes$interpro_descrip)])

rubis_amo_genes <- rubis_amo_genes %>%
  add_column(car_mag = rubis_amo_genes$mag, .before = 2)
rubis_amo_genes <- rubis_amo_genes %>%  
  mutate(car_mag = sub('MAG', 'CarAnox_mtb2', rubis_amo_genes$car_mag)) %>%
  select(-c(seq_len, sig_acc, sig_descrip, start_loc, stop_loc, database)) %>%
  rename('qseqid' = 'interpro_anno')
rubis_amo_genes <- rubis_amo_genes %>%
  select(5, 1, 2, 3, 4, 6)
rubis_amo_genes <- rubis_amo_genes %>%
  add_column('phylum' = NA, .after = 3) %>%
  add_column('class' = NA, .after = 3)

for (i in seq_along(rubis_amo_genes$car_mag)) {
  for (k in seq_along(taxonomy_draftMags$bin_id)) {
    if (as.character(rubis_amo_genes$car_mag)[i] == as.character(taxonomy_draftMags$bin_id)[k]) {
      rubis_amo_genes$class[i] <- as.character(taxonomy_draftMags$Class)[k]
      rubis_amo_genes$phylum[i] <- as.character(taxonomy_draftMags$Phylum)[k]
      next
    }
  }
}

#rename columns to match blast table (made next)
rubis_amo_genes <- rubis_amo_genes %>% rename('CarAnox_mag' = 'car_mag', 'annotation' = 'interpro_descrip')

#pick the gene annotation with the highest bitscore, and remove the rest
rubis_amo_genes <- rubis_amo_genes %>% 
  arrange(gene) %>% 
  group_by(gene) %>% 
  top_n(1, score)

# manually check and get rid of ties
rubis_amo_genes <- as.data.frame(rubis_amo_genes)
rubis_amo_genes <- rubis_amo_genes %>% slice(-c(45, 46))
rubis_amo_genes <- rubis_amo_genes %>% select(-c(score))
rubis_amo_genes <- rubis_amo_genes %>% select(1,2,3,6,4,5,7)

#  
blast_index <- list.files(path = '~/Documents/cariaco/blast_results/', pattern = '*.txt')
setwd('~/Documents/cariaco/blast_results/')
blast_res <- lapply(blast_index, function (x) {read.csv(x, sep = '\t', fill = TRUE,
                    col.names = c('qseqid', 'gene', 'pident', 'length', 'mismatch', 'gapopen', 
                                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'),
                    na.strings = '', header = FALSE, stringsAsFactors = FALSE)})

blast_res <- do.call('rbind', lapply(blast_res, as_tibble))

blast_res$qseqid[blast_res$qseqid == '2264262959'] <- '2264262959_A288L16DRAFT_00600'
blast_res$qseqid[blast_res$qseqid == '2682017948'] <- '2682017948_Ga0131072_101230'

#
queryKey <- read.csv('~/Documents/cariaco/queryseq_key.txt', 
                     sep = '\t', col.names = c('gene_name', 'gene_descrip'), header = FALSE,
                     stringsAsFactors = FALSE)

#pick the gene annotation with the highest bitscore, and remove the rest
blast_finalTable <- blast_res %>% 
  select(1,2, 11, 12) %>%
  arrange(gene) %>% 
  group_by(gene) %>% 
  top_n(1, bitscore) %>%
  add_column(taxonomy = NA) %>%
  add_column(annotation = NA)

#there were some ties for highest bitscore, this command removed some ties based on lowest e-value
blast_finalTable <- blast_finalTable %>%
  top_n(1, evalue) %>%
  ungroup() %>%
  as.data.frame()

# manually removed one of the elements from each of the last few ties
# these genes are the same in terms of annotation, just from diff representative microbes
blast_finalTable <- blast_finalTable %>% slice(-40)
blast_finalTable <- blast_finalTable %>% slice(-75)
blast_finalTable <- blast_finalTable %>% slice(-127)
blast_finalTable <- blast_finalTable %>% slice(-222)
blast_finalTable <- blast_finalTable %>% slice(-279)
blast_finalTable <- blast_finalTable %>% slice(-354)
blast_finalTable <- blast_finalTable %>% slice(-370)
blast_finalTable <- blast_finalTable %>% slice(-427)

#verify there is only one representation of each gene
View(table(sort(fct_infreq(blast_finalTable$gene))))

# add gene annotations
for (i in seq_along(blast_finalTable$qseqid)) {
  for (k in seq_along(queryKey$gene_name)) {
    if (as.character(blast_finalTable$qseqid)[i] == as.character(queryKey$gene_name)[k]) {
      blast_finalTable$annotation[i] <- as.character(queryKey$gene_descrip)[k]
      next
    }
  }
}

blast_finalTable <- blast_finalTable %>%
  add_column(c_mag = blast_finalTable$gene, .before = 2) %>%
  separate(c_mag, into = c('ctg', 'mag'), sep = '-') %>%
  select(-c('ctg')) %>%
  add_column('CarAnox_mag' = NA, .before = 2)

for (i in seq_along(blast_finalTable$mag)) {
  for (k in seq_along(mp_key$mag_id)) {
    if (as.character(blast_finalTable$mag)[i] == as.character(mp_key$mag_id)[k]) {
      blast_finalTable$CarAnox_mag[i] <- as.character(mp_key$bin_id)[k]
      blast_finalTable$taxonomy[i] <- as.character(mp_key$full_phylum)[k]
      next
    }
  }
}

blast_finalTable <- blast_finalTable %>%
  separate(taxonomy, into = c('taxonomy', 'num'), sep = '_') %>%
  select(-c('num')) %>%
  rename('phylum' = 'taxonomy') %>%
  add_column('class' = NA, .before = 'phylum')

for (i in seq_along(blast_finalTable$CarAnox_mag)) {
  for (k in seq_along(taxonomy_draftMags$bin_id)) {
    if (as.character(blast_finalTable$CarAnox_mag)[i] == as.character(taxonomy_draftMags$bin_id)[k]) {
      blast_finalTable$class[i] <- as.character(taxonomy_draftMags$Class)[k]
      next
    }
  }
}


# merge the rubis_amo and blast_final tables into one and tidy t --------
a <- as.character(unique(sort(fct_infreq(blast_finalTable$class))))
blast_finalTable$class <- factor(blast_finalTable$class, levels = a)
blast_finalTable$numeric_class <- as.numeric(blast_finalTable$class)
blast_finalTable <- blast_finalTable %>%
  arrange(numeric_class, CarAnox_mag) %>%
  select(-c(numeric_class, bitscore, evalue))

blast_finalTable <- blast_finalTable %>%
  select(1,4,2,3,5,6,7)

write.table(blast_finalTable, file = '~/Documents/cariaco/blast-results-table.tsv', quote = FALSE,
            row.names = FALSE, sep = '\t')

merged_blast <- bind_rows(blast_finalTable, rubis_amo_genes)

c_merged_blast <- merged_blast
c_merged_blast$annotation <- c_merged_blast$annotation %>%
  fct_collapse('Ammonia/methane monooxygenase' = 
                 c("Ammonia/methane monooxygenase, subunitB, N-terminal"  ,
                   "Ammonia monooxygenase/particulate methane monooxygenase, subunit A",                                     
                   "Ammonia monooxygenase/particulate methane monooxygenase, subunit C",
                   "Ammonia/particulate methane monooxygenase, subunit A superfamily",                                       
                   "Ammonia/methane monooxygenase, subunit B, C-terminal"))

c_merged_blast$annotation <- sub('anaerobic carbon-monoxide dehydrogenase catalytic subunit.*',
                                 'anaerobic carbon-monoxide dehydrogenase catalytic subunit',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('dissimilatory-type sulfite reductase subunit alpha.*',
                                 'dissimilatory-type sulfite reductase subunit alpha',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('.ulfide.quinone reductase.*',
                                 'sulfide-quinone reductase', c_merged_blast$annotation)
c_merged_blast$annotation <- sub('nitrate reductase alpha subunit.*', 'nitrate reductase alpha subunit',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('ATP citrate synthase.*', 'ATP citrate synthase',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('ATP.citrate lyase', 'ATP citrate lyase', c_merged_blast$annotation)
c_merged_blast$annotation <- sub('sulfite reductase alpha subunit',
                                 'dissimilatory-type sulfite reductase subunit alpha',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('CO dehydrogenase/acetyl-CoA synthase complex subunit epsilon.*',
                                 'CO dehydrogenase/acetyl-CoA synthase complex subunit epsilon',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- sub('Succinyl-CoA synthetase, alpha subunit.*',
                                 'Succinyl-CoA synthetase, alpha subunit',
                                 c_merged_blast$annotation)
c_merged_blast$annotation <- c_merged_blast$annotation %>%
  fct_collapse('RuBisCO' = c('RuBisCO large subunit, N-terminal domain superfamily',
                             'RuBisCO_small', 'Ribulose bisphosphate carboxylase, small chain'))
c_merged_blast$annotation <- c_merged_blast$annotation %>%
  fct_collapse('Ammonia/methane monooxygenase' = c('Ammonia monooxygenase/particulate methane monooxygenase, subunit C',
                                                   'Ammonia monooxygenase/particulate methane monooxygenase, subunit A',
                                                   'Ammonia/methane monooxygenase, subunit B, C-terminal',
                                                   'Ammonia/methane monooxygenase, subunitB, N-terminal',
                                                   'Ammonia/particulate methane monooxygenase, subunit A superfamily'))
c_merged_blast <- as_tibble(c_merged_blast) %>%
  mutate(annotation = as.factor(annotation)) %>%
  mutate(mag = as.factor(mag)) %>%
  mutate(gene = as.factor(gene)) %>%
  select(4,1,2,3,5,6,7,8)
blast_concat <- c_merged_blast %>%
  select(-c(qseqid, gene)) %>%
  pivot_wider(names_from = annotation, values_from = count, values_fn = list(count = sum))

a <- as.character(unique(sort(fct_infreq(blast_concat$class))))
blast_concat$class <- factor(blast_concat$class, levels = a)
blast_concat$numeric_class <- as.numeric(blast_concat$class)
blast_concat <- blast_concat %>%
  arrange(numeric_class, mag) %>%
  select(-c(numeric_class))

write.table(blast_concat, file = '~/Documents/cariaco/blast-results-concatenated.tsv', quote = FALSE,
            row.names = FALSE, sep = '\t')



write.table(c_merged_blast$gene, file = '~/Documents/cariaco/blast-genes-list.txt', quote = FALSE,
            row.names = FALSE, col.names = FALSE) #eol = ','


#import mag sizes (in terms of base pairs) -----------------------------

mag_bp <- read.csv('~/Documents/cariaco/final-fasta-base-counts.txt', sep = '\n', stringsAsFactors = FALSE,
                   header = FALSE)
mag_bp <- mag_bp %>% rename('mag' = 'V1')
mag_bp$mag <- sub('draft-mags-subset-[0-9]/', '', mag_bp$mag)
mag_bp$mag <- sub('draft-mags-subset-[0-9][0-9]/', '', mag_bp$mag)
mag_bp$mag <- sub('.fa', '', mag_bp$mag)
mag_bp_names <- mag_bp %>%
  filter(str_detect(mag, 'CarAnox'))
mag_bp <- mag_bp %>%
  filter(str_detect(mag, '^[0-9]'))
colnames(mag_bp) <- 'bp'

mag_bp <- bind_cols(mag_bp_names, mag_bp)
mag_bp <- mag_bp %>% add_column(mbp = mag_bp$bp)
mag_bp <- mag_bp %>% 
  mutate(mbp = as.numeric(mbp)) %>%
  mutate(mbp = mbp / 1e+06) %>%
  mutate(mbp = round(mbp, digits = 2))
mag_bp$mag <- sub('CarAnox_mtb2', 'MAG', mag_bp$mag)


# characterize clusters based on their known cluster hits within MiBIG DB --------
knownClust <- as.data.frame(read.csv('~/Documents/cariaco/knownclusterblast-results-allcontigs.txt', sep = '\n', 
                                     header = FALSE, stringsAsFactors = FALSE))
knownClust <- unique(knownClust); rownames(knownClust) <- c()

knownClust <- knownClust %>%
  rename('gene' = 'V1') %>%
  add_column('geneclust' = NA)

#
for (i in seq_along(knownClust$gene)) {
  for (k in seq_along(gbk_gene$sm_class)) {
    if (as.character(knownClust$gene)[i] == as.character(gbk_gene$smBGC)[k]) {
      knownClust$geneclust[i] <- as.character(gbk_gene$geneclust)[k]
      next
    }
  }
}

#
all_clusters <- as.data.frame(unique(gbk_gene$geneclust)); names(all_clusters) <- 'geneclust'
clusts_withHits <- as.data.frame(unique(knownClust$geneclust)); names(clusts_withHits) <- 'geneclust'
clusts_noKnownHits <- subset(all_clusters, !(all_clusters$geneclust %in% clusts_withHits$geneclust))

clusts_noKnownHits$sm_class <- NA
clusts_noKnownHits$phylum <- NA
#
for (i in seq_along(clusts_noKnownHits$geneclust)) {
  for (k in seq_along(gbk_gene$geneclust)) {
    if (as.character(clusts_noKnownHits$geneclust)[i] == as.character(gbk_gene$geneclust)[k]) {
      clusts_noKnownHits$sm_class[i] <- as.character(gbk_gene$sm_class)[k]
      next
    }
  }
}

for (i in seq_along(clusts_noKnownHits$geneclust)) {
  for (k in seq_along(gbk_gene$geneclust)) {
    if (as.character(clusts_noKnownHits$geneclust)[i] == as.character(gbk_gene$geneclust)[k]) {
      clusts_noKnownHits$phylum[i] <- as.character(gbk_gene$phylum)[k]
      next
    }
  }
}

clusts_noKnownHits <- within(clusts_noKnownHits, sm_class <- factor(sm_class, levels = names(sort(table(sm_class), decreasing = TRUE))))

clusts_noKnownHits %>% ggplot(aes(clusts_noKnownHits$sm_class)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

#
clusts_withHits$sm_class <- NA

for (i in seq_along(clusts_withHits$geneclust)) {
  for (k in seq_along(gbk_gene$geneclust)) {
    if (as.character(clusts_withHits$geneclust)[i] == as.character(gbk_gene$geneclust)[k]) {
      clusts_withHits$sm_class[i] <- as.character(gbk_gene$sm_class)[k]
      next
    }
  }
}

clusts_withHits <- clusts_withHits %>%
  mutate('MAG' = geneclust) %>%
  separate('MAG', into = c('mag', 'magnum', 'clustnum'), sep = '_') %>%
  unite('mag', c('mag', 'magnum'), sep = '_') %>%
  select(-c('clustnum'))

#
unk <- unique(clusts_withHits$mag[clusts_withHits$sm_class == 'Unclassified'])
write.table(unk, file = '~/Documents/cariaco/magNames-with-unk-smBGC-class', quote = FALSE,
            row.names = FALSE, col.names = FALSE)
                            

# parse hmm antibiotic resistance hits ------------------------------------
hmmRes <- read.csv('~/Documents/cariaco/R_analysis/csv-hmm-final-results.txt', header = FALSE, fill = TRUE)

genes_with_hmmHits <- as.data.frame(hmmRes$V1) %>%
  rename('gene' = 'hmmRes$V1') %>%
  mutate(gene = sub('_1.cluster.*', '', genes_with_hmmHits$gene)) %>%
  add_column('geneclust' = NA) %>%
  add_column('phylum' = NA) %>%
  add_column('sm_class' = NA)

for (i in seq_along(genes_with_hmmHits$gene)) {
  for (k in seq_along(gbk_gene$smBGC)) {
    if (as.character(genes_with_hmmHits$gene)[i] == as.character(gbk_gene$smBGC)[k]) {
      genes_with_hmmHits$geneclust[i] <- as.character(gbk_gene$geneclust)[k]
      genes_with_hmmHits$phylum[i] <- as.character(gbk_gene$phylum)[k]
      genes_with_hmmHits$sm_class[i] <- as.character(gbk_gene$sm_class)[k]
      next
    }
  }
}

hmmHits_onGeneClusters <- genes_with_hmmHits %>%
  filter(gene != '-') %>%
  select(-c(gene))

hmmHits_onGeneClusters <- aggregate(. ~ hmmHits_onGeneClusters$geneclust, hmmHits_onGeneClusters, 
                                    function(x) toString(unique(x)))
hmmHits_onGeneClusters <- hmmHits_onGeneClusters %>% select(-c(`hmmHits_onGeneClusters$geneclust`))

hmmHits_onGeneClusters <- within(hmmHits_onGeneClusters, phylum <- factor(phylum, 
                                                                    levels = names(sort(table(phylum), decreasing = TRUE))))

hmmHits_onGeneClusters <- within(hmmHits_onGeneClusters, sm_class <- factor(sm_class, 
                                                                          levels = names(sort(table(sm_class), decreasing = TRUE))))

hmmHits_onGeneClusters %>% ggplot(aes(sm_class, fill = phylum)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

hmmHits_onGCs_PKS_NRPS <- hmmHits_onGeneClusters %>%
  filter(sm_class == 'Polyketide (PK)' | sm_class == 'Non-ribosomal peptide (NRP)' |
           sm_class == 'NRP-PK hybrid')

hmmHits_onGCs_PKS_NRPS %>% ggplot(aes(sm_class, fill = phylum)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

#checking what kinds of smBGCs are within MAGs with rhodopsin genes - NOTE: this analysis was not included in the report - 
                                    # it was exploratory
rhodop <- geneAnno %>%
  filter(str_detect(sig_descrip, 'rhodopsin') | str_detect(interpro_descrip, 'rhodopsin'))
View(rhodop)

rhodop_mags <- gbk_gene %>%
  filter(mag %in% unique(rhodop$mag))


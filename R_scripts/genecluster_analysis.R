
# Load packages, import data into R, format the data ----------------------

# load potentially useful packages for this analysis
packs <- c('tidyverse', 'caret', 'ranger', 'RColorBrewer', 'devtools', 'reticulate', 'reshape2',
           'pheatmap', 'gplots', 'viridis', 'gridExtra', 'data.table')
lapply(packs, library, character.only = TRUE)


# set wd to location of the genecluster files
mypath = "~/Documents/cariaco/gbk-final-geneclusters/"
setwd(mypath)

# index vector of the .txt file names
fileIndex <- list.files(path = mypath, pattern = "*.txt")

# import all genecluster text files
gbk_geneclusters <- lapply(fileIndex, FUN = read.table, sep ="\t")

# combine them
gbk_geneclusters <- do.call("rbind", lapply(gbk_geneclusters, as.data.frame)) 

# copy geneclusters to edited, final version
gbk_gene <- gbk_geneclusters
gbk_gene <- gbk_gene %>%
  dplyr::rename('smBGC' = 'V1', 'tax' = 'V2', 'sm_class' = 'V3', 'ctg' = 'V4') %>%
  select(-c(V5)) %>%
  separate(tax, into = c('magdraft', 'magnum', 'ex'), sep = '_') %>%
  unite('mag', c('magdraft', 'magnum'), sep = '_') %>%
  select(-c(ex))

gbk_gene$sm_class <- gbk_gene$sm_class %>%
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

gbk_gene$sm_class <- factor(gbk_gene$sm_class, levels = c('Aryl polyene', 'RiPP', 'Phosphonate', 'Lactone', 'Ladderane',
                                                                      'Other secondary metabolite*', 'Unclassified', 'Terpene', 
                                                                      'NRP-PK hybrid', 'Non-ribosomal peptide (NRP)', 'Polyketide (PK)'))

gbk_gene <- gbk_gene %>%
  mutate(ctg = strsplit(as.character(ctg), ';')) %>%
  unnest(ctg)

gbk_gene <- gbk_gene %>%
  mutate(smBGC = gbk_gene$mag) %>%
  unite('smBGC', c('ctg', 'smBGC'), sep = '-')

gbk_gene <- gbk_gene %>%
  add_column(tax = gbk_gene$mag) %>%
  add_column(phy = gbk_gene$phylum)

gbk_gene <- gbk_gene %>%
  separate(tax, into = c('mg', 'magnum'), sep = '_') %>%
  unite(taxmag, c('phy', 'magnum'), sep = '_') %>%
  select(-c(mg))

# a function to visualize categorical columns
plot_bars <- function(df, title = NULL, y_label = NULL, tag = NULL) {
  cols <- colnames(df)
  for (col in cols) {
    if(is.factor(df[[col]])) {
      p = ggplot(df, aes_string(col)) +
        geom_bar() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
        labs(title = title, y = y_label, tag = tag)
      return(p)
    }
  }
}


# SM and Taxonomic analysis -----------------------------------------------

#import table with information about all MAGs
draft_mags <- read.csv("~/Documents/cariaco/tables_magInfo/mags_75compl_5cont.csv")
taxonomy_draftMags <- draft_mags[, c(1,7)] # keep only relevant columns from the table

# separate the concatenated taxonomic column into indiviudal columns by tax. rank
taxonomy_draftMags <- taxonomy_draftMags %>%
  separate(classification, into = c('Domain', 'Phylum', 'Class', 'Order', 'Family',
                                    'Genus', 'Species'), sep = ';')

draft_phyla <- as.data.frame(taxonomy_draftMags$Phylum); names(draft_phyla) <- "Phylum"

# edit the phyla to remove unneccesary text, or classify unknown phyla as 'Unclassified'
draft_phyla$Phylum <- ifelse(draft_phyla$Phylum == 'p__', sub('p__', 'Unclassified', draft_phyla$Phylum), 
                             sub('p__', "", draft_phyla$Phylum)) 
draft_phyla$Phylum <- as.factor(draft_phyla$Phylum) # set the Phylum to class factor for plotting

mp_key <- taxonomy_draftMags[c(1, 3)]
mp_key$Phylum <- ifelse(mp_key$Phylum == 'p__', sub('p__', 'Unclassified', mp_key$Phylum),
                        sub('p__', '', mp_key$Phylum))

taxonomy_draftMags$Class <- ifelse(taxonomy_draftMags$Class == 'c__', sub('c__', 'Unclassified', taxonomy_draftMags$Class),
                                   sub('c__', '', taxonomy_draftMags$Class))
taxonomy_draftMags$bin_id <- sub('-', '_', taxonomy_draftMags$bin_id)
taxonomy_draftMags$Phylum <- ifelse(taxonomy_draftMags$Phylum == 'p__', sub('p__', 'Unclassified', taxonomy_draftMags$Phylum),
                                   sub('p__', '', taxonomy_draftMags$Phylum))

# set working directory to save plots
setwd('~/Documents/cariaco/figures/')

# archaeal phyla
arch_phyla <- subset(draft_phyla, Phylum %in% c('Nanoarchaeota', 'Thermoplasmatota', 
                                   'Iainarchaeota', 'Altiarchaeota', 'UAP2', 'Aenigmarchaeota',
                                   'Crenarchaeota', 'Halobacterota'))
rownames(arch_phyla) <- c()
arch_phyla <- within(arch_phyla, Phylum <- 
                        factor(Phylum, levels = names(sort(table(Phylum), decreasing = TRUE))))

# bacterial phyla ############################
bact_phyla <- subset(draft_phyla, !(Phylum %in% c('Nanoarchaeota', 'Thermoplasmatota', 
                                                 'Iainarchaeota', 'Altiarchaeota', 'UAP2', 'Aenigmarchaeota',
                                                 'Crenarchaeota', 'Halobacterota')))
rownames(bact_phyla) <- c()
bact_phyla <- within(bact_phyla, Phylum <- 
                       factor(Phylum, levels = names(sort(table(Phylum), decreasing = TRUE))))

# final bacterial phyla (cut phyla with only 1 representative MAG)
final_bp <- subset(bact_phyla, !(Phylum %in% c('Aerophobota', 'Chloroflexota_B', 'Cyanobacteria',
                                               'Dadabacteria', 'Delongbacteria', 'Dependentiae',
                                               'FEN-1099', 'Fermentibacterota', 'Fibrobacterota',
                                               'Gemmatimonadota_A', 'KSB1', 'OLB16', 'Poribacteria',
                                               'Ratteibacteria', 'UBP17', 'UBP3', 'Verrucomicrobiota_A')))

# plot bacterial MAGs
bp <- plot_bars(final_bp, y_label = 'Number of MAGs per phylum', tag = '1a')

# plot archaeal MAGs
ap <- plot_bars(arch_phyla, y_label = 'Number of MAGs per phylum', tag = '1b')

# group the plots together with grid.arrange, set quality and size and save plot grid
png(filename = 'arch_bact_mags_distr.png', units = 'in', width = 8, height = 8, res = 500)
grid.arrange(bp, ap, layout_matrix = rbind(c(1,1,1,1,1), c(2,2,NA,NA,NA)))
dev.off()

# convert to columns to class character for matching
taxonomy_draftMags$Phylum <- as.character(taxonomy_draftMags$Phylum)
taxonomy_draftMags$bin_id <- as.character(taxonomy_draftMags$bin_id)

draftContigs_geneclusters <- 
  allContigs_geneclusters[allContigs_geneclusters$Genome %in% taxonomy_draftMags$bin_id, ]
rownames(draftContigs_geneclusters) <- c()

draftContigs_geneclusters$Genome <- as.character(draftContigs_geneclusters$Genome)
draftContigs_geneclusters <- add_column(draftContigs_geneclusters, Phylum = NA, 
                                      .after = 'Genome')

# For each element of the allContigs_geneclusters$Genome
# column, scanned the entire taxonomy_allMags$bin_id column for a match, then added the match to the 
# allContigs_geneclusters$Phylum column for that iteration of the loop, building the Phylum column up
# one by one until complete
for (i in seq_along(draftContigs_geneclusters$Genome)) {
  for (k in seq_along(taxonomy_draftMags$bin_id)) {
    if (draftContigs_geneclusters$Genome[i] == taxonomy_draftMags$bin_id[k]) {
      draftContigs_geneclusters$Phylum[i] <- taxonomy_draftMags$Phylum[k]
      next
    }
  }
}

# edit the phyla to remove unneccesary text, or classify unknown phyla as 'Unclassified'
draftContigs_geneclusters$Phylum <- ifelse(draftContigs_geneclusters$Phylum == 'p__', 
                                           sub('p__', 'Unclassified', draftContigs_geneclusters$Phylum), 
                             sub('p__', "", draftContigs_geneclusters$Phylum))

# set the Phylum to class factor for plotting
draftContigs_geneclusters$Phylum <- as.factor(draftContigs_geneclusters$Phylum) 


# plot

#creating custom palette to fit the amount of Phyla represented in the plot
colorCount <- length(unique(draftContigs_geneclusters$Phylum))
getPalette <- colorRampPalette(brewer.pal(n = 8, name = 'Dark2'))(colorCount)

# duplicate the dataframe just in case it gets messed up
copy_draftContigs_geneclusters <- draftContigs_geneclusters

# rename NA's to 'Unclassified'
draftContigs_geneclusters$Phylum <- as.character(draftContigs_geneclusters$Phylum)
draftContigs_geneclusters$Phylum[is.na(draftContigs_geneclusters$Phylum)] <- 'Unclassified'
draftContigs_geneclusters$Phylum <- as.factor(draftContigs_geneclusters$Phylum)

draftContigs_geneclusters <- within(draftContigs_geneclusters, Secondary_metabolite_class <- 
              factor(Secondary_metabolite_class, levels = 
                       names(sort(table(Secondary_metabolite_class), 
                                  decreasing = TRUE))))

##### ALL CONTIGS - COLLAPSING CATEGORIES OF SECONDARY METABOLITES ########## 

# collapse categories of secondary metabolites
draftContigs_sm_collapsed <- draftContigs_geneclusters

# the collapse step
draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Polyketide (PK)' = c('t1pks', 't1pks-PUFA', 
    't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 'transatpks-t1pks',
    't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks',
    't2pks-t1pks', 't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 'nucleoside-otherks',
    't2pks-ladderane', 't2pks', 't2pks-arylpolyene', 't2pks-arylpolyene-resorcinol', 't3pks', 'phosphonate-t3pks'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('NRP-PK hybrid' = c('t1pks-nrps', 'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Non-ribosomal peptide (NRP)' = c('arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Aryl polyene' = c('arylpolyene', 'arylpolyene-resorcinol', 'arylpolyene-ladderane'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Other secondary metabolite*' = c('acyl_amino_acids', 'nucleoside', 'PUFA', 'pbde', 'indole', 
                                                'amglyccycl', 'siderophore', 'resorcinol', 'resorcinol-ladderane',
                                                'ectoine'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Unclassified' = 'other')

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Ladderane' = 'ladderane')

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Phosphonate' = c('phosphonate', 'phosphonate-terpene', 'phosphonate-bacteriocin'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Terpene' = c('terpene', 'bacteriocin-lantipeptide-terpene'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Lactone' = c('butyrolactone', 'hserlactone', 'hserlactone-acyl_amino_acids'))

draftContigs_sm_collapsed$Secondary_metabolite_class <- draftContigs_sm_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('RiPP' = c('lassopeptide', 'cyanobactin', 'thiopeptide', 'head_to_tail', 'microviridin', 
                          'sactipeptide', 'proteusin', 'lantipeptide',
                          'bacteriocin-lassopeptide', 'bacteriocin', 'bacteriocin-terpene', 
                          'proteusin-bacteriocin', 'bacteriocin-lantipeptide', 'thiopeptide-bacteriocin',
                          'bacteriocin-thiopeptide', 'bacteriocin-proteusin'))

# drop unused factor levels
draftContigs_sm_collapsed$Secondary_metabolite_class <- factor(draftContigs_sm_collapsed$Secondary_metabolite_class)

#pie chart of all secondary metabolite classes
draft_pie <- as.data.frame(table(draftContigs_sm_collapsed$Secondary_metabolite_class))
draft_pie$percent <- draft_pie$Freq / sum(draft_pie$Freq) * 100
draft_pie$percent <- round(draft_pie$percent, digits = 1)
names(draft_pie)[1] <- 'Secondary_metabolite_class'
draft_pie <- draft_pie[order(draft_pie$percent, decreasing = TRUE), ]
rownames(draft_pie) <- c()
draft_pie$Secondary_metabolite_class <-
  factor(draft_pie$Secondary_metabolite_class, levels =
           draft_pie$Secondary_metabolite_class[order(draft_pie$percent, decreasing = TRUE)])

pie_colorCount <- length(unique(draft_pie$Secondary_metabolite_class))
pie_Palette <- colorRampPalette(brewer.pal(n = 12, name = 'Dark2'))(pie_colorCount)

#png(filename = 'all_smbgc.png', units = 'in', width = 7, height = 5, res = 1000)
#draft_pie %>% ggplot(aes(x = "", y = percent,
#                         fill = Secondary_metabolite_class)) +
#  geom_bar(width = 1, stat = 'identity') +
#  coord_polar(theta = 'y') +
#  theme(axis.text.x = element_blank(),
#        legend.text = element_text(size = 8),
#        legend.title = element_text(size = 10),
#        legend.background = element_blank(), panel.border = 
#          element_rect(color = 'black', fill = NA),
#        legend.box.background = element_rect(color = 'black')) +
#  theme_void() +
#  geom_text(aes(label = percent(percent/100)), position = position_stack(vjust = 0.5), size = 2) +
#  scale_fill_manual(values = pie_Palette) +
#  labs(fill = 'smBGC class')
#dev.off()

library(plotly)

#png(filename = 'pie_chart.png', units = 'in', width = 6, height = 5, res = 500)
plot_ly(draft_pie, labels = ~Secondary_metabolite_class, values = ~Freq, type = 'pie',
        textposition = 'inside',
        textinfo = 'label+percent+value',
        insidetextfont = list(color = '#FFFFFF'),
        text = ~Secondary_metabolite_class,
        marker = list(colors = pie_Palette,
                      line = list(color = '#FFFFFF', width = 1)),
        showlegend = FALSE) %>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
#dev.off()
htmlwidgets::saveWidget(pie, file = 'pie_chart_Q.html')


# PKS/NRPS pie collapse step
pknrp_collapsed <- subset(draftContigs_geneclusters, Secondary_metabolite_class %in% c('t1pks', 
                      't1pks-PUFA', 't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 
                      'transatpks-t1pks', 't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks', 
                      't2pks-t1pks', 't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 
                      'nucleoside-otherks', 't2pks-ladderane', 't2pks', 't2pks-arylpolyene', 
                      't2pks-arylpolyene-resorcinol', 't3pks', 'phosphonate-t3pks', 't1pks-nrps', 
                      'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps',
                      'arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps'))

pknrp_collapsed$Secondary_metabolite_class <- pknrp_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Type 1 PKS' = c('t1pks', 't1pks-PUFA', 
                                't1pks-PUFA-otherks', 't1pks-otherks', 'otherks', 'transatpks', 'transatpks-t1pks',
                                't1pks-terpene', 'bacteriocin-t1pks', 't3pks-t1pks', 't2pks-t1pks', 
                                't1pks-arylpolyene', 'thiopeptide-otherks', 'arylpolyene-otherks', 
                                'nucleoside-otherks'))

pknrp_collapsed$Secondary_metabolite_class <- pknrp_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Type 2 PKS' = c('t2pks-ladderane', 't2pks', 't2pks-arylpolyene', 't2pks-arylpolyene-resorcinol'))

pknrp_collapsed$Secondary_metabolite_class <- pknrp_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('Type 3 PKS' = c('t3pks', 'phosphonate-t3pks'))

pknrp_collapsed$Secondary_metabolite_class <- pknrp_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('NRPS-PKS hybrid' = c('t1pks-nrps', 'transatpks-t1pks-nrps', 'terpene-t1pks-nrps', 'bacteriocin-t1pks-nrps'))

pknrp_collapsed$Secondary_metabolite_class <- pknrp_collapsed$Secondary_metabolite_class %>% 
  fct_collapse('NRPS' = c('arylpolyene-nrps', 'nrps', 'bacteriocin-nrps', 'transatpks-nrps'))

# drop unused factor levels
pknrp_collapsed$Secondary_metabolite_class <- factor(pknrp_collapsed$Secondary_metabolite_class)

#pie chart of PK/NRP secondary metabolite classes
pknrp_pie <- as.data.frame(table(pknrp_collapsed$Secondary_metabolite_class))
pknrp_pie$percent <- pknrp_pie$Freq / sum(pknrp_pie$Freq) * 100
pknrp_pie$percent <- round(pknrp_pie$percent, digits = 1)
names(pknrp_pie)[1] <- 'Secondary_metabolite_class'
pknrp_pie <- pknrp_pie[order(pknrp_pie$percent, decreasing = TRUE), ]
rownames(pknrp_pie) <- c()
pknrp_pie$Secondary_metabolite_class <-
  factor(pknrp_pie$Secondary_metabolite_class, levels =
           pknrp_pie$Secondary_metabolite_class[order(pknrp_pie$percent, decreasing = TRUE)])

library(scales)
pie_colorCount <- length(unique(pknrp_pie$Secondary_metabolite_class))
pie_Palette <- colorRampPalette(brewer.pal(n = 12, name = 'Dark2'))(pie_colorCount)

png(filename = 'pknrp_smbgc.png', units = 'in', width = 6.5, height = 4, res = 300)
pknrp_pie %>% ggplot(aes(x = "", y = percent,
                         fill = Secondary_metabolite_class)) +
  theme_void() +
  geom_bar(width = 1, stat = 'identity') +
  coord_polar(theta = 'y') +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = percent(percent/100)), position = position_stack(vjust = 0.5), size=4) +
  scale_fill_manual(values = pie_Palette) +
  labs(fill = 'PK/NRP class', title = 'Distribution of PK/NRP smBCGs')
dev.off()




# plot secondary metabolite distributions by phyla, normalized by MAG counts per phylum
dp_key <- as.data.frame(table(draft_phyla))
draftContigs_sm_collapsed$cat <- paste(draftContigs_sm_collapsed$Phylum, 
                                       draftContigs_sm_collapsed$Secondary_metabolite_class, sep = '`')
dp_w_sm_table <- as.data.frame(table(draftContigs_sm_collapsed$cat))
dp_w_sm_table <- dp_w_sm_table %>%
  separate(Var1, into = c('Phylum', 'Secondary_metabolite_class'), sep = '`')
dp_w_sm_table$mag_abun <- NA
dp_w_sm_table$Phylum <- as.character(dp_w_sm_table$Phylum); dp_key$draft_phyla <- as.character(dp_key$draft_phyla)

for (i in seq_along(dp_w_sm_table$Phylum)) {
  for (k in seq_along(dp_key$draft_phyla)) {
    if (dp_w_sm_table$Phylum[i] == dp_key$draft_phyla[k]) {
      dp_w_sm_table$mag_abun[i] <- dp_key$Freq[k]
      next
    }
  }
}

dp_w_sm_table$Phylum <- as.factor(dp_w_sm_table$Phylum); dp_key$draft_phyla <- as.factor(dp_key$draft_phyla)

# normalize the SM distributions by MAG abundance per phyla
dp_w_sm_table$norm_sm_freq <- dp_w_sm_table$Freq / dp_w_sm_table$mag_abun


#creating custom palette to fit the amount of Phyla represented in the plot
smc_colorCount <- length(unique(dp_w_sm_table$Secondary_metabolite_class))
smc_getPalette <- colorRampPalette(brewer.pal(n = 12, name = 'Paired'))(smc_colorCount)
smc_getPalette <- c("khaki", "#2A7FB7", "#99CD91", "#52AF43", "#B89B74", "#ED4F50", 
                    "#A6CEE3", "#FDA440", "#FFCCCB", "#B294C7", "#825D99","#A6CEE3", 
                    "#B15928")

#sorting factors by abundance order, create a key with highest to lowest total normalized count
dpsm_agg <- dp_w_sm_table[c(1, 5)]
dpsm_agg <- aggregate(. ~ dpsm_agg$Phylum, dpsm_agg, sum)
dpsm_agg <- dpsm_agg %>%
  select(-Phylum) %>%
  arrange(desc(norm_sm_freq)) %>%
  dplyr::rename(phylum = `dpsm_agg$Phylum`)

# order levels of phyla for plot by normalized count total, in descending order
dp_w_sm_table$Phylum <- factor(dp_w_sm_table$Phylum, levels = dpsm_agg$phylum)

# add an asterix (*) next to underrepresented phyla (phyla w/ only 1 MAG representative)
dp_w_sm_table$Phylum <- dp_w_sm_table$Phylum %>% 
  fct_collapse(paste(expression(bold('*')), 'Aerophobota', sep = '') = '*Aerophobota') %>%
  fct_collapse('*Cyanobacteria' = 'Cyanobacteria*') %>%
  fct_collapse('*Delongbacteria' = 'Delongbacteria*') %>%
  fct_collapse('*FEN-1099' = 'FEN-1099*') %>%
  fct_collapse('*Fibrobacterota' = 'Fibrobacterota*') %>%
  fct_collapse('*Gemmatimonadota_A' = 'Gemmatimonadota_A*') %>%
  fct_collapse('*KSB1' = 'KSB1*') %>%
  fct_collapse('*OLB16' = 'OLB16*') %>%
  fct_collapse('*Poribacteria' = 'Poribacteria*') %>%
  fct_collapse('*UBP17' = 'UBP17*') %>%
  fct_collapse('*UBP3' = 'UBP3*')
paste

# plot the normalized SM distributions for bacterial MAGs
bact_sm <- subset(dp_w_sm_table, !(Phylum %in% c('Nanoarchaeota', 'Thermoplasmatota', 
                                                 'Iainarchaeota', 'Altiarchaeota', 'UAP2', 'Aenigmarchaeota',
                                                 'Crenarchaeota', 'Halobacterota')))

bact_sm$Secondary_metabolite_class <- factor(bact_sm$Secondary_metabolite_class,
                                             levels = c('Aryl polyene', 'Lactone', 'Ladderane', 'Non-ribosomal peptide (NRP)',
                                                        'NRP-PK hybrid', 'Other secondary metabolite*', 'Phosphonate', 
                                                        'Polyketide (PK)', 'RiPP', 'Terpene', 'Unclassified'))

bact <- ggplot(bact_sm, 
  aes(x = bact_sm$Phylum, y = bact_sm$norm_sm_freq,
  fill = bact_sm$Secondary_metabolite_class)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10, 
                                   face = c('bold','plain','bold','bold','plain','plain',
                                   'plain','plain','plain','plain','bold','bold',
                                   'bold','bold','plain','plain','plain','plain',
                                   'bold','bold','plain','plain','plain','plain',
                                   'plain','plain','plain','plain','plain','bold',
                                   'bold','plain','plain','plain','plain','plain',
                                   'plain','plain')),
        axis.title = (element_text(size = 8)),
        legend.position = c(0.96, 0.96), legend.justification = c(0.96, 0.96),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.background = element_blank(), panel.border = 
          element_rect(color = 'black', fill = NA),
        legend.box.background = element_rect(color = 'black')) +
  labs(x = 'Phylum', y = 'Normalized biosynthetic cluster count per MAG',
       fill = 'Secondary metabolite class', tag = '1a') +
  scale_fill_manual(values = smc_getPalette)

# plot the normalized SM distributions for archaeal MAGs
arch_sm <- subset(dp_w_sm_table, Phylum %in% c('Nanoarchaeota', 'Thermoplasmatota', 
                                                 'Iainarchaeota', 'Altiarchaeota', 'UAP2', 'Aenigmarchaeota',
                                                 'Crenarchaeota', 'Halobacterota'))

arch_sm$Secondary_metabolite_class <- factor(arch_sm$Secondary_metabolite_class,
                                             levels = c('Aryl polyene', 'Lactone', 'Ladderane', 'Non-ribosomal peptide (NRP)',
                                                        'NRP-PK hybrid', 'Phosphonate', 'Polyketide (PK)', 'RiPP', 'Terpene', 
                                                        'Other secondary metabolite*', 'Unclassified'))

arch <- ggplot(arch_sm, 
  aes(x = arch_sm$Phylum, y = arch_sm$norm_sm_freq,
  fill = arch_sm$Secondary_metabolite_class)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10), 
        axis.title = (element_text(size = 8)), legend.position = 'none') +   #, legend.title = element_text(size = 8),
        #legend.position = c(0.96, 0.96), legend.justification = c(0.96, 0.96),
        #legend.text = element_text(size = 7),
        #legend.background = element_blank(), panel.border = 
          #element_rect(color = 'black', fill = NA),
        #legend.box.background = element_rect(color = 'black')) +
  labs(x = 'Phylum', y = 'Normalized biosynthetic cluster count per MAG', tag = '1b') +
  scale_fill_manual(values = smc_getPalette[c(7:10, 6, 11)])

png(filename = 'bact_arch_distr_FINAL_1.png', units = 'in', width = 8, height = 7, res = 500)
grid.arrange(bact, arch, layout_matrix = rbind(c(1,1,1,1,2)))
dev.off()

# Working with meganomes - the reads --------------------------------------
metagenome_names <- read.table('~/Documents/cariaco/R_analysis/metagenome_names.txt', sep = '\n')

metagenome_names$V1 <- sub('_R.*', '', metagenome_names$V1)
metagenome_names <- distinct(metagenome_names)

R1_readNames <- as.vector(read.table('~/Documents/cariaco/R_analysis/R1_readNames.txt', sep = '\n'))
R1_readNames <- as.character(R1_readNames$V1)
R2_readNames <- read.table('~/Documents/cariaco/R_analysis/R2_readNames.txt', sep = '\n')
R2_readNames <- as.character(R2_readNames$V1)

for (i in seq_along(metagenome_names$V1)) {
  for (k in seq_along(R1_readNames)) {
    if (str_detect(R1_readNames[k], metagenome_names$V1[i]) == TRUE) {
      metagenome_names$R1[i] <- R1_readNames[k]
      next
    }
  }
}

for (i in seq_along(metagenome_names$V1)) {
  for (k in seq_along(R2_readNames)) {
    if (str_detect(R2_readNames[k], metagenome_names$V1[i]) == TRUE) {
      metagenome_names$R2[i] <- R2_readNames[k]
      next
    }
  }
}

names(metagenome_names) <- c('sample', 'r1', 'r2')
write.table(metagenome_names,"samples.txt",sep="\t",row.names=FALSE, quote = FALSE)

final_metaNames <- read.delim('~/Documents/cariaco/R_analysis/2nd_final_samples.txt'); View(final_metaNames)
summary(str_detect(as.character(final_metaNames$r1), as.character(final_metaNames$sample)))
summary(str_detect(as.character(final_metaNames$r2), as.character(final_metaNames$sample)))

# Generating Pearson correlations between the draft mags ------------------
library(reshape2)

options(max.print=1000000000)

mean_coverage <- t(read.delim(file="~/Documents/cariaco/R_analysis/mean_coverage.txt",
                              header = TRUE,
                              stringsAsFactors = FALSE,check.names = F, row.names = 1))

correlations <- melt(data = cor(mean_coverage, use="complete.obs", method="pearson"))

names(correlations) <- c('MAG_1', 'MAG_2', 'correlation')

write.table(correlations, "~/Documents/cariaco/R_analysis/REDUNDANT-MAGs-PEARSON.txt", sep="\t", quote=FALSE,  col.names=NA)

redMAGs_aff <- data.frame(draft_mags$bin_id); names(redMAGs_aff) <- 'mag_id'

for (i in seq_along(redMAGs_aff$mag_id)) {
  for (k in seq_along(taxonomy_draftMags$bin_id)) {
    if (str_detect(as.character(taxonomy_draftMags$bin_id)[k], as.character(redMAGs_aff$mag_id)[i]) == TRUE) {
      redMAGs_aff$phylum[i] <- taxonomy_draftMags$Phylum[k]
    }
  }
}

redMAGs_aff$mag_id <- sub('-', "_", redMAGs_aff$mag_id)
redMAGs_aff$phylum <- ifelse(redMAGs_aff$phylum == 'p__', sub('p__', 'Unclassified', redMAGs_aff$phylum), 
                             sub('p__', "", redMAGs_aff$phylum))

write.table(redMAGs_aff, "~/Documents/cariaco/R_analysis/REDUNDANT-MAGs-AFFILIATIONS.txt", sep = "\t", quote = FALSE, col.names = FALSE,
            row.names = FALSE)

stats_redMAGS <- read.delim("~/Documents/cariaco/R_analysis/REDUNDANT-MAGs-STATS.txt", header = TRUE, stringsAsFactors = FALSE,
                            check.names = FALSE)

taxonomy_draftMags$bin_id <- sub('-', "_", taxonomy_draftMags$bin_id)

for (i in seq_along(stats_redMAGS$MAG)) {
  for (k in seq_along(taxonomy_draftMags$bin_id)) {
    if (str_detect(as.character(taxonomy_draftMags$bin_id)[k], as.character(stats_redMAGS$MAG)[i]) == TRUE) {
      stats_redMAGS$domain[i] <- taxonomy_draftMags$Domain[k]
    }
  }
}

stats_redMAGS$domain <- sub('d__', '', stats_redMAGS$domain)
write.table(stats_redMAGS, "~/Documents/cariaco/R_analysis/REDUNDANT-MAGs-STATS.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

# mag relative abundance workflow -----------------------------------------

samples_summary <- read.delim('~/Documents/cariaco/R_analysis/samples_summary.txt')
samples_totalReadCounts <- read.csv('~/Documents/cariaco/R_analysis/metagenome-read-counts.txt', sep = ':', header = FALSE)

names(samples_summary)[1] <- 'sample'
names(samples_totalReadCounts) <- c('sample', 'total_read_count')
samples_totalReadCounts$sample <- sub('paired.fastq.gz', '', samples_totalReadCounts$sample)
samples_totalReadCounts$sample <- sub('_R1_', '', samples_totalReadCounts$sample);samples_totalReadCounts$sample <- sub('_R2_', '', samples_totalReadCounts$sample)
samples_totalReadCounts <- distinct(samples_totalReadCounts)

samples_summary$sample <- sub('_IN_RMAGS', '', samples_summary$sample)

for (i in seq_along(samples_summary$sample)) {
  for (k in seq_along(samples_totalReadCounts$sample)) {
    if (str_detect(samples_totalReadCounts$sample[k], samples_summary$sample[i]) == TRUE) {
      samples_summary$total_read_count[i] <- samples_totalReadCounts$total_read_count[k]
    }
  }
}

# realized that samples_summary$total_read_count needs to be doubled to account for both R1 and R2 read sets, not just 
# the value of one of them
samples_summary$total_read_count <- samples_summary$total_read_count * 2

relAbun <- read.delim('~/Documents/cariaco/R_analysis/bins_percent_recruitment.txt')
names(relAbun)[1] <- 'sample'
relAbun$sample <- sub('_IN_RMAGS', '', relAbun$sample)

relAbun <- column_to_rownames(relAbun, var = 'sample')
relAbun <- t(relAbun)
relAbun <- as.data.frame(relAbun)

# multiplying each value from bins_percent_recruitment.txt by the total reads mapped, then dividing by the total read 
# counts (for each metagenome)
for (i in seq_along(colnames(relAbun))) {
  for (k in seq_along(samples_summary$sample)) {
    if (str_detect(samples_summary$sample[k], colnames(relAbun)[i]) == TRUE) {
      relAbun[i] <- relAbun[i] * samples_summary$total_reads_mapped[k]
      relAbun[i] <- relAbun[i] / samples_summary$total_read_count[k]
    }
  }
}

#save a table of the relative abundances of the MAGs across all samples
write.table(relAbun, "~/Documents/cariaco/R_analysis/Real_relative_abundance.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = TRUE)

#reform sample summary, find tot percent (%) of metagenome reads which mapped to samples
samples_summary <- samples_summary %>%
  mutate(percent_reads_mapped = total_reads_mapped / total_read_count * 100) %>%
  mutate(sample = sub('_D[0-9]', '', samples_summary$sample))
samples_summary <- samples_summary %>%
  mutate(sample = sub('M', 'May_', samples_summary$sample)) %>%
  mutate(sample = sub('N', 'Nov_', samples_summary$sample)) %>%
  mutate(sample = sub('PA', 'PA_', samples_summary$sample)) %>%
  mutate(sample = sub('FL', 'FL_', samples_summary$sample))

write.table(samples_summary, '~/Documents/cariaco/figures/metagen-reads-info-supp-table.txt', sep = '\t', quote = FALSE,
            row.names = FALSE)

#transverse matrix to put MAGS in columns, samples in rows
t_relAbund <- t(relAbun)

#coerce into dataframe to allow for column-wise calculations
t_relAbund <- as.data.frame(t_relAbund)

#get max of each MAG column associated with its MAG
max_AbunMAGS <- as.data.frame(colnames(t_relAbund))
max_AbunMAGS$max_abundance <- apply(t_relAbund, 2, max)
names(max_AbunMAGS)[1] <- 'mag'

#get mean of MAG max vector
mean_maxAbun <- mean(max_AbunMAGS$max_abundance)

#how many MAGS will be leftover with cutoffs equal to 20% and 50% of the average max abundance
length(max_AbunMAGS$max_abundance); length(max_AbunMAGS$max_abundance[max_AbunMAGS$max_abundance > (.2 * mean_maxAbun)]) #330 left over
length(max_AbunMAGS$max_abundance); length(max_AbunMAGS$max_abundance[max_AbunMAGS$max_abundance > (.5 * mean_maxAbun)]) #180 left over

summary(colnames(t_relAbund) == max_AbunMAGS$mag)

#index t_relAbun to just take MAGs with a max abundance larger than 50% of the mean max abundance
columns_180 <- t_relAbund[max_AbunMAGS$max_abundance > (.5* mean_maxAbun)]
maxColumns_180 <- as.data.frame(colnames(columns_180))
maxColumns_180$max_abundance <- apply(columns_180, 2, max)
summary(maxColumns_180$max_abundance > (.5 * mean_maxAbun))


# Relative abundance - heatmap creation -----------------------------------

palette_ttra <- colorRampPalette(brewer.pal(n = 8, name = 'Dark2')); View(palette_ttra)

refPallete <- colorRampPalette(color = brewer.pal(n = 8, name = 'Dark2'))(8)

#create depth and sample metadata dataframe for follow creation of annotation column for heatmap
tra_depth <- data.frame(sample = colnames(trimmed_relAbun), depth = colnames(trimmed_relAbun))
tra_depth$depth <- sub('MFL', '', tra_depth$depth); tra_depth$depth <- sub('MPA', '', tra_depth$depth)
tra_depth$depth <- sub('NFL', '', tra_depth$depth); tra_depth$depth <- sub('NPA', '', tra_depth$depth)
tra_depth$depth <- sub('_D1_\\d.*', '', tra_depth$depth); tra_depth$depth <- sub('_D2_\\d.*', '', tra_depth$depth)
tra_depth$depth <- sub('_D3_\\d.*', '', tra_depth$depth)
tra_depth$sample <- sub('_D1_\\d.*', '', tra_depth$sample); tra_depth$sample <- sub('_D2_\\d.*', '', tra_depth$sample)
tra_depth$sample <- sub('_D3_\\d.*', '', tra_depth$sample)
tra_depth <- distinct(tra_depth)

t_relAbun_colAnno <- data.frame(sample = colnames(trimmed_relAbun))

for (i in seq_along(t_relAbun_colAnno$sample)) {
  for (k in seq_along(tra_depth$sample)) {
    if (str_detect(as.character(t_relAbun_colAnno$sample)[i], 
                   as.character(tra_depth$sample)[k]) == TRUE) {
      t_relAbun_colAnno$depth[i] <- tra_depth$depth[k]
    }
  }
}

str_detect(as.character(t_relAbun_colAnno$sample), as.character(t_relAbun_colAnno$depth))
t_relAbun_colAnno <- column_to_rownames(t_relAbun_colAnno, var = 'sample')
names(t_relAbun_colAnno) <- 'Depth'

t_tra <- data.matrix(t(trimmed_relAbun))

trimmed_relAbun <- as.data.frame(t(columns_180))
trimmed_relAbun <- round(trimmed_relAbun, digits = 3)
trimmed_relAbun <- log1p(trimmed_relAbun)

#trimmed_relAbun <- log10(t(trimmed_relAbun))
#remove -infinity values
#is.na(trimmed_relAbun)<-sapply(trimmed_relAbun, is.infinite)
#convert -infity values to 0
#trimmed_relAbun[is.na(trimmed_relAbun)] <- 0

t_tra <- data.matrix(t(trimmed_relAbun))

breaks2 <- c(seq(0, 0.01, by = 0.01), seq(0.02, max(t_tra), by = 0.01))

colors <- c(colorRampPalette(colors = 'white')(length(breaks2)/100),
            colorRampPalette(colors = c('#ABD9E9', '#FEE090', '#D6604D', '#D73027'))(length(breaks2)))

#copy and then change t_tra column names to representative phyla
copy_ttra <- t_tra

mp_key <- mp_key %>%
  add_column(bin_copy = mp_key$bin_id, .after = 'Phylum') %>%
  separate(bin_copy, into = c('mag_prefix', 'mag_num'), sep = '-') %>%
  unite(full_phylum, c('Phylum', 'mag_num'), sep = '_') %>%
  select(-mag_prefix)

mp_key$full_phylum[182] <- 'Bdellovibrionota-B_723'; mp_key$full_phylum[208] <- 'Bdellovibrionota-B_1534'

mp_key$bin_id <- sub('-', '_', mp_key$bin_id)

for (i in seq_along(colnames(copy_ttra))) {
  for (k in seq_along(mp_key$bin_id)) {
    if (colnames(copy_ttra)[i] == as.character(mp_key$bin_id)[k]) {
      colnames(copy_ttra)[i] <- mp_key$full_phylum[k]
      next
    }
  }
}

mp_key <- mp_key %>%
  add_column(mag = 'MAG') %>%
  add_column(mg = mp_key$full_phylum) %>%
  separate(mg, into = c('taxa', 'taxa_num'), sep = '_') %>%
  unite(mag_id, c('mag', 'taxa_num'), sep = '_') %>%
  select(-c(taxa))

final_gene$phylum <- NA

for (i in seq_along(final_gene$mag)) {
  for (k in seq_along(mp_key$mag_id)) {
    if (as.character(final_gene$mag)[i] == as.character(mp_key$mag_id)[k]) {
      final_gene$phylum[i] <- mp_key$full_phylum[k]
      next
    }
  }
}

# create water column layer annotation column
wl_colAnno <- t_relAbun_colAnno
wl_colAnno$Layer <- ifelse(wl_colAnno$Depth <= 237, wl_colAnno$Layer <- 'Oxycline',
                            ifelse(wl_colAnno$Depth > 237 & wl_colAnno$Depth < 900, 
                                   wl_colAnno$Layer <- 'Shallow anoxic', wl_colAnno$Layer <- 'Euxinic'))
wl_colAnno <- wl_colAnno %>%
  select(-Depth)

# edit rownames for matrix and annotation col
rownames(copy_ttra) <- sub('_D\\d_', '_', rownames(copy_ttra))
rownames(wl_colAnno) <- sub('_D\\d_', '_', rownames(wl_colAnno))

# set working directory to figure directory
setwd('~/Documents/cariaco/figures/')

#reorder the matrix columns by MAG phylum frequency, create annotation column
colAnno <- as.data.frame(colnames(copy_ttra)); colnames(colAnno) <- 'mag'; colAnno$mag <- as.character(colAnno$mag)
colAnno$mag[32] <- 'Gemmatimonadota-A_1254'
colAnno <- colAnno %>%
  add_column(phy = colAnno$mag) %>%
  separate(phy, into = c('Phylum', 'num')) %>%
  select(-c('num'))

colAnno$Phylum <- as.character(colAnno$Phylum); colAnno$Phylum[colAnno$Phylum == 'AABM5'] <- 'AABM5-125-24'; colAnno$Phylum <- as.factor(colAnno$Phylum)
colAnno$Phylum <- as.character(colAnno$Phylum); colAnno$Phylum[colAnno$Phylum == 'FEN'] <- 'FEN-1099'; colAnno$Phylum <- as.factor(colAnno$Phylum)
colAnno$Phylum <- as.character(colAnno$Phylum); colAnno$Phylum[c(26,28,30,31,34,36,37,40:43,
                                                               46:50)] <- 'Gammaproteobacteria'; colAnno$Phylum <- as.factor(colAnno$Phylum)
colAnno$Phylum <- as.character(colAnno$Phylum); colAnno$Phylum[c(27,29,32,33,35,38,39,44,45)] <- 'Alphaproteobacteria'; colAnno$Phylum <- as.factor(colAnno$Phylum)

a <- as.character(unique(sort(fct_infreq(colAnno$Phylum))))
colAnno$Phylum <- factor(colAnno$Phylum, levels = a)
colAnno$numeric_phylum <- as.numeric(colAnno$Phylum)
colAnno <- colAnno %>%
  arrange(numeric_phylum) %>%
  select(-c(numeric_phylum))

colAnno <- colAnno %>%
  arrange(-row_number())

copy2_ttra <- copy_ttra
copy2_ttra <- as.data.frame(copy2_ttra)
colnames(copy2_ttra)[32] <- 'Gemmatimonadota-A_1254'
copy2_ttra <- copy2_ttra[, colAnno$mag]

colAnno <- colAnno %>%
  column_to_rownames(var = 'mag')

rel_rowAnno <- wl_colAnno %>%
  t() %>%
  as.data.frame() %>%
  select(11,12,35,36,31:34,7:10,1,2,25,26,3,4,27,28,5,6,29,30,
         23,24,47,48,43:46,19:22,13,14,37,38,15,16,39,40,17,18,41,42) %>%
  t() %>%
  as.data.frame() %>%
  add_column('Fraction' = c(rep('FL', 24), rep('PA', 24)))

copy3ttra <- copy2_ttra
copy3ttra <- copy3ttra[rownames(rel_rowAnno), ]
copy3ttra <- copy3ttra %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag') %>%
  arrange(-row_number()) %>%
  column_to_rownames(var = 'mag') %>%
  t()

relPalette <- list(Layer = c('Oxycline' = 'steelblue2', 'Shallow anoxic' = 'lightcoral', 'Euxinic' = 'springgreen2'),
                       Fraction = c('FL' = 'plum3', 'PA' = 'gold'),
                       Phylum = c('Planctomycetota' = 'khaki', 'Desulfobacterota' = 'darkslategray1', 'Myxococcota' = 'orange',
                                  'Gammaproteobacteria' = 'orchid1', 'Chloroflexota' = 'slateblue1', 'Marinisomatota' = 'aquamarine', 
                                  'Alphaproteobacteria' = 'gold', 'Krumholzibacteriota' = 'deeppink', 'Omnitrophota' = 'darkviolet',
                                  'Verrucomicrobiota' = 'khaki', 'Thermoplasmatota' = 'firebrick1',
                                  'AABM5-125-24' = 'gainsboro', 'Bacteroidota' = 'springgreen',
                                  'SAR324' = 'chocolate4', 'Cloacimonadota' = 'grey28',
                                  'Crenarchaeota' = 'yellowgreen', 'Acidobacteriota' = 'deepskyblue',
                                  'Actinobacteriota' = '#666666', 'Eisenbacteria' = '#B89B74', 'Gemmatimonadota' = 'lightsteelblue1',
                                  'Nanoarchaeota' = 'ivory', 'Nitrospinota' = '#2A7FB7', 
                                  'FEN-1099' = '#9B7329', 'Fermentibacterota' = '#1B9E77',
                                  'Firmicutes' = 'darksalmon', 'KSB1' = 'wheat', 'Latescibacterota' = '#AD4C9E',
                                  'Patescibacteria' = 'thistle1', 'UAP2' = 'turquoise', 'Unclassified' = '#5A8950'))

oddsPalette <- list(Layer = c('May oxycline' = 'plum3', 'May shallow anoxic' = 'paleturquoise1', 'May euxinic' = 'darkseagreen1',
                                 'November oxycline' = 'springgreen3', 'November shallow anoxic' = 'gold', 
                                 'November euxinic' = 'darkorchid3'),
                       Phylum = c('Planctomycetota' = 'khaki', 'Desulfobacterota' = 'darkslategray1', 'Myxococcota' = 'orange',
                                  'Gammaproteobacteria' = 'orchid1', 'Chloroflexota' = 'slateblue1', 'Marinisomatota' = 'aquamarine', 
                                  'Alphaproteobacteria' = 'gold', 'Krumholzibacteriota' = 'deeppink', 'Omnitrophota' = 'darkviolet',
                                  'Verrucomicrobiota' = 'khaki', 'Thermoplasmatota' = 'firebrick1',
                                  'AABM5-125-24' = 'gainsboro', 'Bacteroidota' = 'springgreen',
                                  'SAR324' = 'chocolate4', 'Cloacimonadota' = 'grey28',
                                  'Crenarchaeota' = 'yellowgreen', 'Acidobacteriota' = 'deepskyblue',
                                  'Actinobacteriota' = '#666666', 'Eisenbacteria' = '#B89B74', 'Gemmatimonadota' = 'lightsteelblue1',
                                  'Nanoarchaeota' = 'ivory', 'Nitrospinota' = '#2A7FB7', 
                                  'FEN-1099' = '#9B7329', 'Fermentibacterota' = '#1B9E77',
                                  'Firmicutes' = 'darksalmon', 'KSB1' = 'wheat', 'Latescibacterota' = '#AD4C9E',
                                  'Patescibacteria' = 'thistle1', 'UAP2' = 'turquoise', 'Unclassified' = '#5A8950'))

#make sample names more readable for the figure
rownames(copy3ttra) <- c("May_FL_900_1", "May_FL_900_2", "Nov_FL_900_1", "Nov_FL_900_2", "Nov_FL_247_1", "Nov_FL_247_2", 
                         "Nov_FL_267_1", "Nov_FL_267_2", "May_FL_295_1", "May_FL_295_2", 
                         "May_FL_314_1", "May_FL_314_2", "May_FL_103_1", "May_FL_103_2", "Nov_FL_143_1", "Nov_FL_143_2", 
                         "May_FL_198_1", "May_FL_198_2", "Nov_FL_200_1", "Nov_FL_200_2",
                         "May_FL_234_1", "May_FL_234_2", "Nov_FL_237_1", "Nov_FL_237_2", "May_PA_900_1", "May_PA_900_2", 
                         "Nov_PA_900_1", "Nov_PA_900_2", "Nov_PA_247_1", "Nov_PA_247_2", 
                         "Nov_PA_267_1", "Nov_PA_267_2", "May_PA_295_1", "May_PA_295_2", "May_PA_314_1", "May_PA_314_2", 
                         "May_PA_103_1", "May_PA_103_2", "Nov_PA_143_1", "Nov_PA_143_2", 
                         "May_PA_198_1", "May_PA_198_2", "Nov_PA_200_1", "Nov_PA_200_2", "May_PA_234_1", "May_PA_234_2", 
                         "Nov_PA_237_1", "Nov_PA_237_2")
rownames(rel_rowAnno) <- rownames(copy3ttra)

# plot the relative abundances for the most abundant 180 MAGs, across all water layers
#png(filename = 'rel_abun_abs_FINAL.png', units = 'in', width = 16.4, height = 10, res = 500)
#pheatmap(copy3ttra, annotation_row = rel_rowAnno, annotation_col = colAnno, cluster_cols = FALSE,
#         border_color = 'grey60', annotation_colors = relPalette, cluster_rows = FALSE,
#         cellwidth = 6, color = colors, treeheight_row = 13, fontsize_row = 9, fontsize_col = 7,
#         treeheight_col = 0, show_colnames = FALSE,
#         annotation_legend = FALSE, legend = FALSE)
#dev.off()

#
test_ttra <- copy2_ttra %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag') %>%
  arrange(-row_number()) %>%
  column_to_rownames(var = 'mag')
View(test_ttra)

#
odds_ratio <- t(trimmed_relAbun)
odds_ratio <- as.data.frame(odds_ratio)
odds_ratio <- rbind((odds_ratio[13:24, ] / odds_ratio[1:12, ]), (odds_ratio[37:48, ] / odds_ratio[25:36, ]))
odds_ratio <- log10(odds_ratio)
is.na(odds_ratio) <- sapply(odds_ratio, is.infinite)
is.na(odds_ratio) <- sapply(odds_ratio, is.nan)
odds_ratio[is.na(odds_ratio)] <- 0
odds_ratio <- t(odds_ratio)
colnames(odds_ratio) <- sub('_D[0-9]', '', colnames(odds_ratio)); colnames(odds_ratio) <- sub('PA', '', colnames(odds_ratio))
odds_ratio <- odds_ratio %>%
  as.data.frame() %>%
  rownames_to_column(var = 'car_mag') %>%
  add_column(mag = NA, .after = 1)

#
for (i in seq_along(odds_ratio$car_mag)) {
  for (k in seq_along(mp_key$bin_id)) {
    if (as.character(odds_ratio$car_mag)[i] == as.character(mp_key$bin_id)[k]) {
      odds_ratio$mag[i] <- as.character(mp_key$full_phylum)[k]
      next
    }
  }
}

#
odds_ratio$mag[32] <- 'Gemmatimonadota-A_1254'
odds_ratio <- odds_ratio %>%
  select(-c(car_mag)) %>%
  column_to_rownames(var = 'mag')

colAnno <- colAnno %>% rownames_to_column(var = 'mag')
odds_ratio <- odds_ratio[colAnno$mag, ]
colAnno <- colAnno %>% column_to_rownames(var = 'mag')

#odds_ratio <- odds_ratio %>%
#  gather(sample, abun, M103_1:N900_2)
#odds_ratio %>% ggplot(aes(sample, sam, size = abun)) +
#  geom_point()


# plot the odds ratio heatmap
OR_color <- colorRampPalette(c("royalblue1", "white", "firebrick1"))(800)

# use floor and ceiling to deal with even/odd length pallettelengths
OR_breaks <- c(seq(min(odds_ratio), 0, length.out = ceiling(800/2) + 1), 
               seq(max(odds_ratio)/800, max(odds_ratio), length.out = floor(800/2)))

# finish creating odds ratio figure
copy_odds_ratio <- odds_ratio
c_colAnno <- colAnno
c_colAnno <- c_colAnno %>% add_column('mag' = rownames(c_colAnno))
c_colAnno <- c_colAnno[rownames(copy_odds_ratio), ]
c_colAnno <- c_colAnno %>% select(-c(mag))

c_colAnno <- c_colAnno %>%
  rownames_to_column(var = 'mag') %>%
  mutate(c_mag = mag) %>%
  mutate(copy_phylum = Phylum) %>%
  separate(c_mag, into = c('mag', 'num'), sep = '_') %>%
  unite(rn, c(copy_phylum, num), sep = '_') %>%
  select(-c(mag)) %>%
  column_to_rownames(var = 'rn')

rownames(copy_odds_ratio) <- rownames(c_colAnno)

copy_odds_ratio <- copy_odds_ratio %>%
  rownames_to_column(var = 'mag') %>%
  mutate(sum = rowSums(copy_odds_ratio[, 2:length(copy_odds_ratio)])) %>%
  mutate(phylum = mag) %>%
  separate(phylum, into = c('phylum', 'num'), sep = '_') %>%
  select(-c(num))

a <- as.character(unique(sort(fct_infreq(copy_odds_ratio$phylum))))
copy_odds_ratio$phylum <- factor(copy_odds_ratio$phylum, levels = a)
copy_odds_ratio$phy_numeric <- as.numeric(copy_odds_ratio$phylum)

copy_odds_ratio <- copy_odds_ratio %>%
  group_by(mag) %>%
  arrange(phy_numeric, desc(sum)) %>%
  ungroup() %>%
  select(-c(phy_numeric, phylum)) %>%
  #arrange(-row_number()) %>%
  column_to_rownames(var = 'mag')
fracPref <- copy_odds_ratio
copy_odds_ratio <- copy_odds_ratio %>% select(-c(sum))

c_colAnno <- c_colAnno %>% add_column('mag' = rownames(c_colAnno))
c_colAnno <- c_colAnno[rownames(copy_odds_ratio), ]
c_colAnno <- c_colAnno %>% select(-c(mag))

c_rowAnno <- wl_colAnno %>%
  rownames_to_column(var = 'sample')
c_rowAnno$sample <- sub('FL', '', c_rowAnno$sample)
c_rowAnno$sample <- sub('PA', '', c_rowAnno$sample)
c_rowAnno <-  aggregate(. ~ sample, c_rowAnno, function(x) toString((unique(x)))) %>%
  column_to_rownames(var = 'sample')
c_rowAnno$Layer[1:6] <- 'May oxycline'
c_rowAnno$Layer[7:10] <- 'May shallow anoxic'
c_rowAnno$Layer[11:12] <- 'May euxinic'
c_rowAnno$Layer[13:18] <- 'November oxycline'
c_rowAnno$Layer[19:22] <- 'November shallow anoxic'
c_rowAnno$Layer[23:24] <- 'November euxinic'

#change odds ratio matrix column names and rowAnno to make them more readable in the figure
colnames(copy_odds_ratio) <- c("May_103_1", "May_103_2", "May_198_1", "May_198_2", "May_234_1", "May_234_2", 
                               "May_295_1", "May_295_2", "May_314_1", "May_314_2", "May_900_1", "May_900_2", 
                               "Nov_143_1", "Nov_143_2", "Nov_200_1", "Nov_200_2", "Nov_237_1", "Nov_237_2", 
                               "Nov_247_1", "Nov_247_2", "Nov_267_1", "Nov_267_2", "Nov_900_1", "Nov_900_2")

rownames(c_rowAnno) <- colnames(copy_odds_ratio)

# call the function to produce the heatmap, with rows and columns annotated  
png(filename = 'odds_ratio_legend_abs_final.png', units = 'in', width = 22, height = 12, res = 500)
pheatmap(t(copy_odds_ratio), annotation_col = c_colAnno, cluster_cols = FALSE, cluster_rows = FALSE,
         border_color = 'grey60', annotation_colors = oddsPalette, annotation_row = c_rowAnno,
         cellwidth = 6, treeheight_row = 13, fontsize_row = 9, fontsize_col = 7,
         treeheight_col = 0, show_colnames = FALSE, color = OR_color, breaks = OR_breaks,
         legend = T, annotation_legend = T)
dev.off()

#editing the rel abun heatmap so that MAGs are in the same order as in the odds ratio heatmap
copy4ttra <- as.data.frame(copy3ttra)

#Rename a few Mags that don't match - proteos to alpha/gamma proteos, etc.
bad <- subset(colnames(copy4ttra), !(colnames(copy4ttra) %in% rownames(copy_odds_ratio)))
colnames(copy4ttra)[62:77] <- sub('Proteobacteria', 'Gammaproteobacteria', colnames(copy4ttra)[62:77])
colnames(copy4ttra)[104:112] <- sub('Proteobacteria', 'Alphaproteobacteria', colnames(copy4ttra)[104:112])
colnames(copy4ttra)[167] <- sub('Gemmatimonadota-A_1254', 'Gemmatimonadota_1254', colnames(copy4ttra)[167])

#put dataframe into same order as odds ratio object
copy4ttra <- copy4ttra[, rownames(copy_odds_ratio)]

png(filename = 'rel_abun_SORTED_FINAL.png', units = 'in', width = 16.4, height = 10, res = 500)
pheatmap(copy4ttra, annotation_row = rel_rowAnno, annotation_col = c_colAnno, cluster_cols = FALSE,
         border_color = 'grey60', annotation_colors = relPalette, cluster_rows = FALSE,
         cellwidth = 6, color = colors, treeheight_row = 13, fontsize_row = 9, fontsize_col = 7,
         treeheight_col = 0, show_colnames = FALSE,
         annotation_legend = FALSE, legend = FALSE)
dev.off()


# format the draft_mags table to use as supplementary table ---------------
draft_mags <- draft_mags %>%
  add_column('genome' = draft_mags$bin_id, .before = 1)
draft_mags <- draft_mags %>%
  mutate(genome = sub('CarAnox_mtb2-', 'MAG_', draft_mags$genome)) %>%
  select(-c(bin_id))

draft_mags <- draft_mags %>%
  add_column('taxmag' = NA, .before = 1)

for (i in seq_along(draft_mags$genome)) {
  for (k in seq_along(mtb2_to_mag_key$bin_id)) {
    if (as.character(draft_mags$genome)[i] == as.character(mtb2_to_mag_key$bin_id)[k]) {
      draft_mags$taxmag[i] <- as.character(mtb2_to_mag_key$taxmag)[k]
      next
    }
  }
}

draft_mags <- draft_mags %>%
  select(-c(genome)) %>%
  rename('genome' = taxmag)

#subset the full table for the mags from the rel abun heatmap
mags_relAbun_key <- data.frame(taxmag = colnames(copy4ttra))
mags_relAbun_key <- mags_relAbun_key %>%
  mutate(mag = taxmag) %>%
  separate(mag, into = c('taxonomy', 'magnum'), sep = '_')
mags_relAbun_key <- mags_relAbun_key %>%
  mutate(mag = paste('MAG', mags_relAbun_key$magnum, sep = '_')) %>%
  select(-c(magnum))
suppl_relAbun_mags <- draft_mags %>%
  filter(genome %in% mags_relAbun_key$mag)

suppl_relAbun_mags <- suppl_relAbun_mags %>%
  column_to_rownames(var = 'genome')
suppl_relAbun_mags <- suppl_relAbun_mags[mags_relAbun_key$mag, ]
suppl_relAbun_mags <- suppl_relAbun_mags %>%
  rownames_to_column(var = 'genome') %>%
  mutate(genome = mags_relAbun_key$taxmag)

write.table(suppl_relAbun_mags, "~/Documents/cariaco/figures/relabun-subsetted-mag-table.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

write.table(draft_mags, "~/Documents/cariaco/figures/full-mag-table.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)


# create table showing preferred fraction preference based on rel --------
rel_ab <- round(relAbun, digits = 3)
rel_ab <- rel_ab %>% rownames_to_column(var = 'mag')
rel_ab$mag <- sub('CarAnox_mtb2', 'MAG', rel_ab$mag)
rel_ab <- rel_ab %>% column_to_rownames(var = 'mag')
rel_ab <- t(rel_ab)

rel_ab <- rbind((rel_ab[13:24, ] / rel_ab[1:12, ]), (rel_ab[37:48, ] / rel_ab[25:36, ]))
rel_ab <- log10(rel_ab)
is.na(rel_ab) <- sapply(rel_ab, is.infinite)
is.na(rel_ab) <- sapply(rel_ab, is.nan)
rel_ab[is.na(rel_ab)] <- 0
rel_ab <- t(rel_ab)
colnames(rel_ab) <- sub('_D[0-9]', '', colnames(rel_ab))
colnames(rel_ab) <- sub('PA', '', colnames(rel_ab))
rel_ab <- as.data.frame(rel_ab)

rel_ab <- rel_ab %>%
  rownames_to_column(var = 'mag') %>%
  mutate(sum = rowSums(rel_ab[, 2:length(rel_ab)])) %>%
  select(mag, sum) %>%
  mutate(preference = ifelse(sum > 0, 'PA', ifelse(sum < 0, 'FL', 'none'))) %>%
  column_to_rownames(var = 'mag')

rel_ab <- rel_ab %>%
  rownames_to_column(var = 'mag') %>%
  select(-c(sum))


write.table(rel_ab, file = '~/Documents/cariaco/fracPref.txt', quote = FALSE, sep = '\t',
            row.names = FALSE)

#convert data frame from a "wide" format to a "long" format

#odds_ratio <- reshape2::melt(odds_ratio, id = c('sample'))

#odds_ratio %>% ggplot(aes(x = sample, y = colnames(odds_ratio))) + 
#  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
#  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
#  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
#  theme(legend.key=element_blank(), 
#        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
#        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
#        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
#        legend.title = element_text(size = 12, face = "bold"), 
#        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#        legend.position = "right") +  
#  scale_y_discrete(limits = rev(levels(odds_ratio$variable)))






# loading libraries -------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(circlize)
  library(ComplexHeatmap)
  library(gridBase)
  library(tidyverse)
  library(limma)
  library(magrittr)
  source("../../src/interactome_functions_v4_updated.R")
})

# loading data ------------------------------------------------------------

ob <- readRDS('../../results/rpca_integrated_annotated.rds')

# preprocessing the data --------------------------------------------------

# check groups
table(ob$DiseaseState)
meta_label <- 'cell_types'

# read in reference LR
LR_database <- read_csv("../../data/Literature Supported Receptor Ligand Pairs Ramilowski.csv")
LR_database <- LR_database[,1:2]
LR_database <- apply(LR_database, 2, alias2SymbolTable, species = "Hs") %>% data.frame()

# list of wanted genes
ligand_genes <- as.character(LR_database$Ligand.ApprovedSymbol...1)
receptor_genes <- as.character(LR_database$Receptor.ApprovedSymbol...2)

# choose cell types that you want to include in interactome
table(ob$cell_types, ob$DiseaseState) # all groups have the same cell types present, so I will use all
cell_types <- levels(as.factor(ob$cell_types))
cell_types <- cell_types[-c(15:17)]

# get genes from object
merge_genes <- rownames(ob)

# split object
Idents(ob) <- 'DiseaseState'
control <- subset(ob, idents = "Healthy")
stim <- subset(ob, idents = "Tumor")

# change metadata column
Idents(control) <- meta_label
Idents(stim) <- meta_label

rm(ob) # just to make more memory space


# Calculations ------------------------------------------------------------

# Calculations
not_found_ligands <- check_genes(ligand_genes, LR_database, merge_genes)
not_found_receptors <- check_genes(receptor_genes, LR_database, merge_genes)

# separate genes into ligands and receptors and get rid of genes not found
not_in_lig <- unique(not_found_ligands$object)
wanted_ligands <- unique(ligand_genes[!ligand_genes %in% not_in_lig])

not_in_rec <- unique(not_found_receptors$object)
wanted_receptors <- unique(receptor_genes[!receptor_genes %in% not_in_rec])

genes <- unique(c(wanted_ligands, wanted_receptors))

# get need info from object into table
# change second vars to what meta.data refers to labels
control_data <- FetchData(control, vars = c(genes, meta_label))
stim_data <- FetchData(stim, vars = c(genes, meta_label))

# Make all possible LR pairs
LRs <- make_LR_pairs(wanted_ligands, wanted_receptors, LR_database)

# Calculate average expression for all cell types for each group
avg_0 <- AverageExpression(object = control, features = genes, verbose = F, assays = 'CorrectedCounts')[[1]]
avg_1 <- AverageExpression(object = stim, features = genes, verbose = F, assays = 'CorrectedCounts')[[1]]

# Create table w/ average expressions for every ligand/receptor and source/target cell combination
LR_table <- create_LR_table(wanted_ligands, wanted_receptors, cell_types, LRs, avg_0, avg_1)
LR_table[,1:4] <- data.frame(apply(LR_table[,1:4], 2, as.character), stringsAsFactors = FALSE)
write_csv(LR_table, '../../results/11-30-2022_interactome_analysis/pasca_panc_lr_table_full.csv')

# Choose threshold --> here I'm using median
summary(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))
threshold <- median(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))

# Filter LR pairs based on average expressions
LR_table <- avg_LR_filt(LR_table, threshold)

# Calculate p-value between groups (can add alpha = )
LR_table <- LR_diff(LR_table, control_data, stim_data, genes, meta_label)

# Filter LR pairs where expression is constant between groups
LR_table <- LR_table[which(LR_table$lig_diff != 0 | LR_table$rec_diff != 0),]

# Sort by ligand
LR_table <- arrange(LR_table, desc(lig_diff != 0 & rec_diff != 0), lig_diff_p_adj)

# Add additional annotation metadata
## continuous: expression (avg log2(FC)) (two columns)
LR_table$log2_avg_FC_lig <- log2(LR_table$avg_lig_1 + 1) - log2(LR_table$avg_lig_0 + 1)
LR_table$log2_avg_FC_rec <- log2(LR_table$avg_rec_1 + 1) - log2(LR_table$avg_rec_0 + 1)

# Save table
write_csv(LR_table, file = '../../results/panc_lr_table_filt.csv')

# Plotting ----------------------------------------------------------------

# > Processing data -------------------------------------------------------
lr <- LR_table
## importing data
# lr <- read.csv('../../results/panc_lr_table_filt.csv')
## Filtering out celltypes we don't need
lr %<>% filter(source %in% c("Fibroblast", "Macrophage", "Epithelial", "Endothelial", "Pericytes"))
lr %<>% filter(target %in% c("Epithelial", "Granulocyte", "Fibroblast", "Macrophage", "Acinar", "Endothelial", "Pericytes", "CD4" ,"CD8"))

## can change this (change source sector)
source_cell <- "Endothelial"

## can change this
## decrease alpha level until get about 200 interactions (~20 ligands)
alpha <- .05 # needs to be < .05
# lr_sub <- lr %>% filter(lig_diff_p_adj < alpha & source == source_cell) %>% arrange(target)
lr_sub <- lr %>% filter(lig_diff_p_adj < alpha, source == source_cell)
lr_sub %>% distinct(ligands,log2_avg_FC_lig) %>% top_n(10, log2_avg_FC_lig) %>% pull(ligands) -> top_ligands
lr_sub <- lr_sub %>% filter(ligands %in% top_ligands) %>% arrange(target)
lr_sub <- lr_sub %>% mutate(source_split_labels = paste(source, "L", sep = "_"),
                            target_split_labels = ifelse(grepl(paste0("^", source_cell), target), paste(target, "R", sep = "_"), target))

### check how many ligands and receptors
length(unique(lr_sub$ligands))
length(unique(lr_sub$receptors))
#check ligands singled out as significant by the alpha inputed above
# unique(lr_sub$ligands)
# "CD274" %in% lr$ligands
# "PDCD1" %in% lr$receptors
# x <- lr[lr$ligands %in% "CD274" & lr$receptors %in% "PDCD1",]

### filtering out interactions by ligands
# lr_sub <- lr_sub %>% filter(ligands != 'ANXA1') # interactions based on ligands
# lr_sub <- lr_sub %>% filter(receptors != 'ABCA1') # interactions based on receptors
# lr_sub <- lr_sub %>% filter(ligands != 'HLA-C' & receptors != 'KIR2DL1') # specific interaction


# > Plot 1 ----------------------------------------------------------------

# Plot #1: include ligands and receptors from same source and target cell (no split)
## make gene df
lr_express <- rbind(lr_sub %>% select(ligands, source, log2_avg_FC_lig) %>%
                      rename('genes'=ligands, 'cell_type'=source, "log2_avg_FC" = log2_avg_FC_lig) %>%
                      unique() %>% mutate(LR = 'L'),
                    lr_sub %>% select(receptors, target, log2_avg_FC_rec) %>%
                      rename('genes'=receptors, 'cell_type'=target, "log2_avg_FC" = log2_avg_FC_rec) %>%
                      unique() %>% mutate(LR = 'R'))
lr_express <- lr_express %>% mutate(cell_type = factor(cell_type))

### refactor cell type levels
x <- levels(lr_express$cell_type)[grep(paste0("^", source_cell), levels(lr_express$cell_type))]
lr_express$cell_type <- fct_relevel(lr_express$cell_type, x)

## can change this (change track 1 colors; make sure there the same number of colors as cell types)
track_1_colors <- c("#d90016","#8200d9","#0700d9","#0062d9","#00d9d2","#d9cb00","#d95a00",
                    "#d90000","#008b8b")

## making matrix list to put into track 2
### can change this (change cell types; make sure source cell is first)
z <- c( "Fibroblast", "Epithelial", "Granulocyte", "Macrophage", "Acinar", "Endothelial", "Pericytes", "CD4" ,"CD8")
mat_list <- sapply(z,
                   function(x) {
                     y <- lr_express$log2_avg_FC[lr_express$cell_type == x]
                     names(y) <- lr_express$genes[lr_express$cell_type == x]
                     y
                   })
names(mat_list) <- z

## making annotations customize table
x <- c()
x <- unlist(sapply(names(mat_list), function(y) x <- c(x,names(mat_list[[y]]))))
annote_label_df <- data.frame(sector = rep(names(mat_list), times = lengths(mat_list)),
                              name = x,
                              show = T,
                              label = x,
                              size = .5,
                              font = 2)
### can change this (change colors and scale)
col_fun <- colorRamp2(seq(-2,2,length.out = 100), colorRampPalette(colors = c('blue', 'white', 'red'))(100))

## making interactions customize table
### can change this (change width of lines)
width <- .3

sector_widths <- table(lr_express$cell_type, lr_express$LR)
temp_sectors <- unique(lr_express$cell_type)
gene_sector_start <- data.frame()
for (i in 1:nrow(sector_widths)) {

  temp_df <- data.frame(LR = rep(c(colnames(sector_widths)[1], colnames(sector_widths)[2]),
                                 times = c(sector_widths[i,1], sector_widths[i,2])),
                        start = seq(0, (sector_widths[i,1] + sector_widths[i,2] - 1)),
                        cell_type = rownames(sector_widths)[i])

  gene_sector_start <- rbind(gene_sector_start, temp_df)
}

gene_sector_start <- gene_sector_start %>% arrange(LR, cell_type)
lr_express$start <- gene_sector_start$start

### can change this (change line color and transparency)
line_color <- '#289c9c'
line_transparent <- .2
line_length <- .675

## plotting circos plot
### initialize plot
circos.par("gap.after" = 2, "cell.padding" = c(0,0,0,0))
circos.initialize(lr_express$cell_type, xlim = cbind(c(0,0), table(lr_express$cell_type)))

### make source and targets track (track 1)
circos.track(lr_express$cell_type,
             track.index = 1,
             ylim = c(0,1),
             bg.col = track_1_colors,
             cell.padding = c(0, 0, 0, 0),
             track.height = .03,
             panel.fun = function(x,y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(2),
                           CELL_META$sector.index, cex = 1, niceFacing = T)
               circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
                           CELL_META$cell.xlim[2], CELL_META$cell.ylim[2])
             })

## set gap b/w source and targets track and annotations track
set_track_gap(gap = cm_h(1.5))

## make annotations track (track 2)
circos.track(ylim = c(0, 1), bg.border = NA, track.height = .025, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index

  m <- mat_list[[sector.index]]

  annote <- annote_label_df[annote_label_df$sector == sector.index,]

  col_mat <- col_fun(m)
  nr <- length(m)
  nc <- 1

  circos.rect(1:nr - 1, rep(nc, nr),
              1:nr, rep(nc + 1, nr),
              border = 'black', col = col_mat)

  circos.text((1:nr - 1) + cm_h(2), 2.5, annote$label, facing = 'clockwise', niceFacing = T, adj = c(0, .5), cex = annote$size, font = annote$font)

})

# interactions
for (i in 1:nrow(lr_sub)) {

  current_source <- lr_sub$source[i]
  current_target <- lr_sub$target[i]

  current_ligand <- lr_sub$ligands[i]
  current_rec <- lr_sub$receptors[i]

  source_start <- lr_express$start[lr_express$cell_type == current_source & lr_express$genes == current_ligand]
  target_start <- lr_express$start[lr_express$cell_type == current_target & lr_express$genes == current_rec]

  source_lim <- c((.5-(width/2))+source_start, (.5+(width/2))+source_start)
  target_lim <- c((.5-(width/2))+target_start, (.5+(width/2))+target_start)

  circos.link(lr_sub$source[i], point1 = source_lim, lr_sub$target[i], point2 = target_lim, col = adjustcolor(line_color, alpha.f = line_transparent),
              rou = line_length)
}

# clear if want to redo plot
circos.clear()



# > Plot 2 ------------------------------------------------------------------
# Plot #2: don't include ligands and receptors from same source and target cell

## filter out same source and target cell type rows
lr_sub <- lr_sub %>% filter(target_split_labels != paste0(source_cell, "_R"))

## make gene df
lr_express <- rbind(lr_sub %>% select(ligands, source, log2_avg_FC_lig) %>%
                      rename('genes'=ligands, 'cell_type'=source, "log2_avg_FC" = log2_avg_FC_lig) %>%
                      unique() %>% mutate(LR = 'L'),
                    lr_sub %>% select(receptors, target, log2_avg_FC_rec) %>%
                      rename('genes'=receptors, 'cell_type'=target, "log2_avg_FC" = log2_avg_FC_rec) %>%
                      unique() %>% mutate(LR = 'R'))
lr_express <- lr_express %>% mutate(cell_type = factor(cell_type))

### refactor cell type levels
x <- levels(lr_express$cell_type)[grep(paste0("^", source_cell), levels(lr_express$cell_type))]
lr_express$cell_type <- fct_relevel(lr_express$cell_type, x)

## can change this (change track 1 colors; make sure there the same number of colors as cell types)
track_1_colors <- c("#d90016","#8200d9","#0700d9","#0062d9","#00d9d2","#d9cb00","#d95a00",
                    "#d90000","#008b8b", "#ffd700", "#32cd32", "#006400")

## making matrix list to put into track 2
### can change this (change cell types; make sure source cell is first)
z <- c( "Endothelial", "Pericytes", "Epithelial", "Macrophage", "Fibroblast",  "Granulocyte",  "Acinar","CD4" ,"CD8")

mat_list <- sapply(z,
                   function(x) {
                     y <- lr_express$log2_avg_FC[lr_express$cell_type == x]
                     names(y) <- lr_express$genes[lr_express$cell_type == x]
                     y
                   })
names(mat_list) <- z

## making annotations customize table
x <- c()
x <- unlist(sapply(names(mat_list), function(y) x <- c(x,names(mat_list[[y]]))))
annote_label_df <- data.frame(sector = rep(names(mat_list), times = lengths(mat_list)),
                              name = x,
                              show = T,
                              label = x,
                              size = .5,
                              font = 2)
### can change this (change colors and scale)
col_fun <- colorRamp2(seq(-2,2,length.out = 100), colorRampPalette(colors = c('blue', 'white', 'red'))(100))

## making interactions customize table
### can change this (change width of lines)
width <- .3

sector_widths <- table(lr_express$cell_type, lr_express$LR)
temp_sectors <- unique(lr_express$cell_type)
gene_sector_start <- data.frame()
for (i in 1:nrow(sector_widths)) {

  temp_df <- data.frame(LR = rep(c(colnames(sector_widths)[1], colnames(sector_widths)[2]),
                                 times = c(sector_widths[i,1], sector_widths[i,2])),
                        start = seq(0, (sector_widths[i,1] + sector_widths[i,2] - 1)),
                        cell_type = rownames(sector_widths)[i])

  gene_sector_start <- rbind(gene_sector_start, temp_df)
}

gene_sector_start <- gene_sector_start %>% arrange(LR, cell_type)
lr_express$start <- gene_sector_start$start

### can change this (change line color and transparency)
line_color <- '#289c9c'
line_transparent <- .2
line_length <- 0.6

## plotting circos plot
### initialize plot
circos.par("gap.after" = 2, "cell.padding" = c(0,0,0,0))
circos.initialize(lr_express$cell_type, xlim = cbind(c(0,0), table(lr_express$cell_type)))

### make source and targets track (track 1)
circos.track(lr_express$cell_type,
             track.index = 1,
             ylim = c(0,1),
             bg.col = track_1_colors,
             cell.padding = c(0, 0, 0, 0),
             track.height = .03,
             panel.fun = function(x,y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(2),
                           CELL_META$sector.index, cex = 1, niceFacing = T)
               circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
                           CELL_META$cell.xlim[2], CELL_META$cell.ylim[2])
             })

## set gap b/w source and targets track and annotations track
set_track_gap(gap = cm_h(1.5))

## make annotations track (track 2)
circos.track(ylim = c(0, 1), bg.border = NA, track.height = .025, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index

  m <- mat_list[[sector.index]]

  annote <- annote_label_df[annote_label_df$sector == sector.index,]

  col_mat <- col_fun(m)
  nr <- length(m)
  nc <- 1

  circos.rect(1:nr - 1, rep(nc, nr),
              1:nr, rep(nc + 1, nr),
              border = 'black', col = col_mat)

  # circos.text((1:nr - 1) + cm_h(2), 2.5, annote$label, facing = 'clockwise', niceFacing = T, adj = c(0, .5), cex = annote$size, font = annote$font)

})

# set_track_gap(gap = cm_h(1.5))

# interactions
for (i in 1:nrow(lr_sub)) {

  current_source <- lr_sub$source[i]
  current_target <- lr_sub$target[i]

  current_ligand <- lr_sub$ligands[i]
  current_rec <- lr_sub$receptors[i]

  source_start <- lr_express$start[lr_express$cell_type == current_source & lr_express$genes == current_ligand]
  target_start <- lr_express$start[lr_express$cell_type == current_target & lr_express$genes == current_rec]

  source_lim <- c((.5-(width/2))+source_start, (.5+(width/2))+source_start)
  target_lim <- c((.5-(width/2))+target_start, (.5+(width/2))+target_start)

  circos.link(lr_sub$source[i], point1 = source_lim, lr_sub$target[i], point2 = target_lim, col = adjustcolor(line_color, alpha.f = line_transparent),
              rou = line_length)
}

# clear if want to redo plot
circos.clear()

# > Plot 3 ----------------------------------------------------------------
# Plot 3: Highlight interactions (without include same source and target cell genes)
## filter out same source and target cell type rows
lr_sub <- lr_sub %>% filter(target_split_labels != paste0(source_cell, "_R"))

## make gene df
lr_express <- rbind(lr_sub %>% select(ligands, source, log2_avg_FC_lig) %>%
                      rename('genes'=ligands, 'cell_type'=source, "log2_avg_FC" = log2_avg_FC_lig) %>%
                      unique() %>% mutate(LR = 'L'),
                    lr_sub %>% select(receptors, target, log2_avg_FC_rec) %>%
                      rename('genes'=receptors, 'cell_type'=target, "log2_avg_FC" = log2_avg_FC_rec) %>%
                      unique() %>% mutate(LR = 'R'))
lr_express <- lr_express %>% mutate(cell_type = factor(cell_type))

### refactor cell type levels
x <- levels(lr_express$cell_type)[grep(paste0("^", source_cell), levels(lr_express$cell_type))]
lr_express$cell_type <- fct_relevel(lr_express$cell_type, x)

## can change this (change track 1 colors; make sure there the same number of colors as cell types)
track_1_colors <- c("#d90016","#8200d9","#0700d9","#0062d9","#00d9d2","#d9cb00","#d95a00",
                    "#d90000","#008b8b", "#ffd700", "#32cd32", "#006400")

## making matrix list to put into track 2
### can change this (change cell types; make sure source cell is first)
z <- c("Fibroblast","Pericytes", "Endothelial", "Macrophage", "Epithelial",   "Granulocyte", "Acinar", "CD4" ,"CD8")
mat_list <- sapply(z,
                   function(x) {
                     y <- lr_express$log2_avg_FC[lr_express$cell_type == x]
                     names(y) <- lr_express$genes[lr_express$cell_type == x]
                     y
                   })
names(mat_list) <- z

## making annotations customize table
x <- c()
x <- unlist(sapply(names(mat_list), function(y) x <- c(x,names(mat_list[[y]]))))
annote_label_df <- data.frame(sector = rep(names(mat_list), times = lengths(mat_list)),
                              name = x,
                              show = T,
                              label = x,
                              size = .5,
                              font = 2)
### can change this (change colors and scale)
col_fun <- colorRamp2(seq(-2,2,length.out = 100), colorRampPalette(colors = c('blue', 'white', 'red'))(100))

## making interactions customize table
### can change this (change width of lines)
width <- .3

sector_widths <- table(lr_express$cell_type, lr_express$LR)
temp_sectors <- unique(lr_express$cell_type)
gene_sector_start <- data.frame()
for (i in 1:nrow(sector_widths)) {

  temp_df <- data.frame(LR = rep(c(colnames(sector_widths)[1], colnames(sector_widths)[2]),
                                 times = c(sector_widths[i,1], sector_widths[i,2])),
                        start = seq(0, (sector_widths[i,1] + sector_widths[i,2] - 1)),
                        cell_type = rownames(sector_widths)[i])

  gene_sector_start <- rbind(gene_sector_start, temp_df)
}

gene_sector_start <- gene_sector_start %>% arrange(LR, cell_type)
lr_express$start <- gene_sector_start$start

### can change this (change line color and transparency)
line_color <- '#289c9c'
line_transparent <- .2
line_length <- .673

### this is what's different from plot 3
interact_pars_df <- data.frame(color = rep(line_color, times = nrow(lr_sub)),
                               transparency = line_transparent,
                               width = width,
                               rou = line_length)

#### highlight different interactions
##### by ligand (and paired receptors)
x <- lr_sub$ligands == 'FN1'
interact_pars_df$color[x] <- 'orange'
interact_pars_df$transparency[x] <- .5

# ##### by receptor (and paired ligands)
# x <- lr_sub$receptors == 'ITGAM'
# interact_pars_df$color[x] <- 'blue'
# interact_pars_df$transparency[x] <- .5

# ##### a certain interaction (ligand, receptor)
# x <- lr_sub$ligands == 'GNAS' & lr_sub$receptors == 'ADORA1'
# interact_pars_df$color[x] <- 'green'
# interact_pars_df$transparency[x] <- .5

##### gray out other interactions
# x <- !(lr_sub$ligands =='AREG' | lr_sub$receptors == 'ITGAM' | (lr_sub$ligands == 'GNAS' & lr_sub$receptors == 'ADORA1'))
x <- !(lr_sub$ligands =='AREG' | lr_sub$receptors == 'ITGAM' | (lr_sub$ligands == 'GNAS' & lr_sub$receptors == 'ADORA1'))
interact_pars_df$color[x] <- '#f7f5f5'
interact_pars_df$transparency[x] <- 1

##### order interactions
interact_order <- rownames(interact_pars_df)[order(factor(interact_pars_df$color,levels = c('#f7f5f5', 'orange', 'blue', 'green')))]

## plotting circos plot
### initialize plot
circos.par("gap.after" = 1, "cell.padding" = c(0,0,0,0))
circos.initialize(lr_express$cell_type, xlim = cbind(c(0,0), table(lr_express$cell_type)))

### make source and targets track (track 1)
circos.track(lr_express$cell_type,
             track.index = 1,
             ylim = c(0,1),
             bg.col = track_1_colors,
             cell.padding = c(0, 0, 0, 0),
             track.height = .03,
             panel.fun = function(x,y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(2),
                           CELL_META$sector.index, cex = 1, niceFacing = T)
               circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
                           CELL_META$cell.xlim[2], CELL_META$cell.ylim[2])
             })

## set gap b/w source and targets track and annotations track
set_track_gap(gap = cm_h(1))

## make annotations track (track 2)
circos.track(ylim = c(0, 1), bg.border = NA, track.height = .025, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index

  m <- mat_list[[sector.index]]

  annote <- annote_label_df[annote_label_df$sector == sector.index,]

  col_mat <- col_fun(m)
  nr <- length(m)
  nc <- 1

  circos.rect(1:nr - 1, rep(nc, nr),
              1:nr, rep(nc + 1, nr),
              border = 'black', col = col_mat)

  circos.text((1:nr - 1) + cm_h(2), 2.5, annote$label, facing = 'clockwise', niceFacing = T, adj = c(0, .5), cex = annote$size, font = annote$font)

})

set_track_gap(gap = cm_h(1))

## interactions
for (i in as.numeric(interact_order)) {

  current_source <- lr_sub$source[i]
  current_target <- lr_sub$target[i]

  current_ligand <- lr_sub$ligands[i]
  current_rec <- lr_sub$receptors[i]

  source_start <- lr_express$start[lr_express$cell_type == current_source & lr_express$genes == current_ligand]
  target_start <- lr_express$start[lr_express$cell_type == current_target & lr_express$genes == current_rec]

  source_lim <- c((.5-(interact_pars_df$width[i]/2))+source_start, (.5+(interact_pars_df$width[i]/2))+source_start)
  target_lim <- c((.5-(interact_pars_df$width[i]/2))+target_start, (.5+(interact_pars_df$width[i]/2))+target_start)

  circos.link(lr_sub$source[i], point1 = source_lim, lr_sub$target[i], point2 = target_lim,
              col = adjustcolor(interact_pars_df$color[i], alpha.f = interact_pars_df$transparency[i]),
              rou = interact_pars_df$rou[i])
}

# clear if want to redo plot
circos.clear()
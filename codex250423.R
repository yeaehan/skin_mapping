library(readxl)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)

### 4i lung

setwd('//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping')

load('Data/241114.RData')


save.image('241114.RData')

####

smallercrop2_done <- read_excel("LB06 BP 4i LaminB1 green BCL2 red 27mar25_R 1_Merged_ch00.tif clean up_no zeros.xlsx", 
                                sheet = "nucleiProps", col_types = c("text", "numeric",
                                                                     "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                                     "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                                     "numeric", "numeric"))


####


####
smallercrop2 <- smallercrop2_done

colnames(smallercrop2)[2:3] <- c("umap_2", "umap_1")
library(tidyr)

smallercrop2 <- smallercrop2[, -1]

smallercrop2 <- drop_na(smallercrop2, 1)

smallercrop2[which(duplicated(smallercrop2$Label)), ] #check if any duplicated label
#dup <- key[which(!duplicated(key$Genotype)), ]

mat1 <- (as.data.frame(smallercrop2))

rownames(mat1) <- smallercrop2$Label
mat1 <- mat1[, colnames(mat1)!="Label"]
mat1 <- mat1[, colnames(mat1)!="umap_2"]
mat1 <- mat1[, colnames(mat1)!="umap_1"]
mat1 <- mat1[, colnames(mat1)!="area"]
mat1 <- mat1[, colnames(mat1)!="area(um)"]
mat1 <- mat1[, colnames(mat1)!="Y(um)"]
mat1 <- mat1[, colnames(mat1)!="X(um)"]
mat1 <- mat1[, colnames(mat1)!="EquivDiameter"]
mat1 <- mat1[, colnames(mat1)!="...1"]
mat1 <- mat1[, colnames(mat1)!="filePath"]
mat1 <- mat1[, colnames(mat1)!="Channel"]
mat1 <- mat1[, colnames(mat1)!="MeanIntensity"]
mat1 <- mat1[, colnames(mat1)!="p16 layer 6"]
mat1 <- mat1[, colnames(mat1)!="p16 Thermal"]
mat1 <- mat1[, colnames(mat1)!="p16 thermo bad"]
mat1 <- mat1[, colnames(mat1)!="p16 bad"]
mat1 <- mat1[, colnames(mat1)!="EquivDiameter(um)"]


table(rowSums(is.na(mat1)) > 0)
colSums(is.na(mat1)) > 0

which(rowSums(is.na(mat1)) > 0)


gene_list <- c(
  'Pan-cytokeratin', 
  'SMA', 
  'MelanA', 
  'Vimentin', 
  'Keratin 5', 
  'p21', 
  'Ki67', 
  'LaminB1', 
  'HMGB1', 
  'BCL2', 
  'CD4', 
  'CD3e', 
  'CD8', 
  'FOXP3', 
  'CD11c', 
  'HLA-DR'
)

# Subset the columns from the matrix
subset_mat <- mat1[, gene_list, drop = FALSE]
####

#### making seurat object
save.image('241114.RData')

smallercrop2_seurat <- CreateSeuratObject(t(mat1))

smallercrop2_seurat <- FindVariableFeatures(smallercrop2_seurat)
smallercrop2_seurat <- NormalizeData(smallercrop2_seurat)
smallercrop2_seurat <- ScaleData(smallercrop2_seurat)
smallercrop2_seurat <- RunPCA(smallercrop2_seurat,features =VariableFeatures(smallercrop2_seurat), approx=FALSE)

smallercrop2_seurat <- FindNeighbors(smallercrop2_seurat, dims = 1:5)
smallercrop2_seurat <- FindClusters(smallercrop2_seurat, resolution = 0.5)

#UMAP calculation itself will just be used as a place holder here
smallercrop2_seurat <- RunUMAP(smallercrop2_seurat, reduction = "pca", n.neighbors=5, dims=1:5)

#get "real" UMAP-values
cell_embeddings <- smallercrop2[, 1:2]
cell_embeddings <- cell_embeddings/5000 ##what does the 5000 here mean?? was 5000 but I will change to 1000
rownames(cell_embeddings) <- rownames(smallercrop2_seurat@meta.data) #check that the rows match!! subsetting the seurat in preprocessing can alter things
smallercrop2_seurat@reductions$umap@cell.embeddings <- as.matrix(cell_embeddings)

##Yeaeun: look for clusters
smallercrop2_seurat <- FindNeighbors(smallercrop2_seurat, dims = 1:10)
smallercrop2_seurat <- FindClusters(smallercrop2_seurat, resolution = 0.5)

#### groups

pc <- smallercrop2_seurat@assays$RNA$counts["Pancytokeratin",] > 0.0
sma <- smallercrop2_seurat@assays$RNA$counts["aSMA",] > 0.0
vimentin <- smallercrop2_seurat@assays$RNA$counts["vimentin",] > 0.0
CD45 <- smallercrop2_seurat@assays$RNA$counts["CD45",] > 0.0

smallercrop2_seurat@meta.data$skincelltype <- case_when(
  sma  ~ 'SMA', 
  vimentin & !CD45 ~ 'vimentin',
  pc ~ 'Pancytokeratin ',
  CD45 ~ 'CD45',
  TRUE ~ 'none'
)


p16 <- smallercrop2_seurat@assays$RNA$counts["p16",] > 0.0
p21 <- smallercrop2_seurat@assays$RNA$counts["p21",] > 0.0

smallercrop2_seurat@meta.data$sencelltype <- case_when(
  p16 ~ "p16",
  p21 ~ "p21",
  TRUE ~ 'none'
)


#save as 800 by 500
Idents(smallercrop2_seurat) <- 'sencelltype'
DimPlot(smallercrop2_seurat,reduction='umap', pt.size=1, order=c('p16', 'p21', 'none'), cols = c("gray", "red", "blue")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "lung_sendim.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

Idents(smallercrop2_seurat) <- 'skincelltype'
DimPlot(smallercrop2_seurat, pt.size=1, order=c('vimentin', 'Pancytokeratin', 'SMA', 'CD45', 'none'), cols = c("orange", "grey", "blue", "#CD4014", "#68B02B")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "lung_celldim.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

####


#### dotplots
Idents(smallercrop2_seurat) <- 'skincelltype'

f <- c("rna_p21", "rna_p16")

DotPlot(smallercrop2_seurat, features = f, dot.scale=10, dot.min = 0.01 , scale=TRUE)+scale_colour_gradient2(low="#006699", mid="#D9D9D9", high="#CC0000")+theme(axis.text.y=element_text(size=14), axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.line = element_blank())+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+scale_y_discrete(limits=rev)
ggsave(filename = "lung_dot2.svg", device = "svg", width=8, height=4, dpi=1000)  # Works fine

####

### pie charts
smallercrop2_seurat@meta.data$p16 <- case_when(
  p16 ~ "pos",
  TRUE ~ 'neg'
)

smallercrop2_seurat@meta.data$p21 <- case_when(
  p21 ~ "pos",
  TRUE ~ 'neg'
)


library(plotly)
library(ggplot2) 

mycols <- c("#3156a6", "#cd8729")

skinsen <- as.data.frame(table(smallercrop2_seurat$skincelltype, smallercrop2_seurat$p21))
skinsen <- as.data.frame(table(smallercrop2_seurat$skincelltype, smallercrop2_seurat$p16))

your_data <- as.data.frame(skinsen %>%
                             group_by(Var1) %>%
                             mutate(Freq = Freq / sum(Freq)))


colnames(your_data) <- c('cell_types', 'p21', 'Freq')
your_data <- your_data[your_data$cell_types != "none", ]
png("lung_cellp21_2.png", res = 1200, width = 10, height = 8, units='in', pointsize=4)

ggplot(data=your_data, aes(x=" ", y=Freq, group=p21, fill=p21)) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +  
  facet_grid(.~ cell_types) +theme_void() +
  scale_fill_manual(values = mycols) 
dev.off()


###

### pie charts per p16/p21


skinsen <- as.data.frame(table(seurat$immcelltype2, seurat$p21))

your_data <- as.data.frame(skinsen %>%
                             group_by(Var1) %>%
                             mutate(Freq = Freq / sum(Freq)))

colnames(your_data) <- c('cell_types', 'p21', 'Freq')

your_data <- your_data[your_data$cell_types != "none", ]

your_data$Freq <- round(your_data$Freq * 100, 2)

your_data <- your_data %>%
  arrange(desc(cell_types)) %>%
  group_by(cell_types) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5 * Freq) %>%
  ungroup()

png("combined_young_skinp21.png", res = 1200, width = 10, height = 8, units='in', pointsize=4)
png("combined_old_imm2p21.png", res = 1200, width = 10, height = 8, units='in', pointsize=4)


ggplot(your_data, aes(x = "", y = Freq, group=p21, color=p21, fill = p21)) +
  geom_bar(width = 2, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Freq), color = "black")+
  scale_fill_manual(values = mycols) +
  facet_grid(.~ cell_types) +theme_void()
dev.off()


###

###
library(ggplot2)
library(patchwork)
library(dplyr)

# Subset for young and old
combined$age <- factor(combined$age, levels = c("young", "old"))

skinsen3 <- as.data.frame(table(combined$age, combined$skincelltype))
immsen3 <- as.data.frame(table(combined$age, combined$immcelltype4))

cells <- as.data.frame(immsen3 %>%
                           group_by(Var1) %>%
                           mutate(Freq = Freq / sum(Freq)))
cells <- cells[cells$Var2 != 'none',]
cells <- cells[cells$Var2 != 'CD45+',]

cells <- cells %>%
  group_by(Var1) %>%
  mutate(Percentage = round(Freq * 100, 1))  # Convert proportion to percentage

library(ggplot2)
stacked_barplot <- ggplot(cells, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Proportion of Skin Cell Types in Young vs. Old",
    x = "Age Group",
    y = "Proportion",
    fill = "Skin Cell Types"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Display the plot
print(stacked_barplot)


stacked_barplot <- ggplot(cells, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  # 1. Change the color palette to a professional, muted scheme
  scale_fill_brewer(palette = "RdGy") + # "Set2" is a good choice for discrete, muted data
  
  # 2. Use theme_classic() or theme_minimal() and remove grid lines
  theme_classic() + # theme_classic is a cleaner default than theme_minimal
  
  # 3. Add labels and titles
  labs(
    title = "Proportion of Skin Cell Types in Young vs. Old",
    x = "Age Group",
    y = "Proportion",
    fill = "Skin Cell Types"
  ) +
  
  # 4. Apply custom text formatting
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # Explicitly remove gridlines and background elements if using theme_minimal
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Display the plot
print(stacked_barplot)


### THE BLACK BORDER & NUMBERS
# Create stacked bar plot with black borders and percentage labels
stacked_barplot <- ggplot(cells, aes(x = Var1, y = Freq, fill = Var2))+
  # 1. Change the color palette to a professional, muted scheme
  scale_fill_brewer(palette = "RdGy") + # "Set2" is a good choice for discrete, muted data +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black borders
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") +  # Percentage labels inside bars
  theme_minimal() +
  labs(
    title = "Proportion of Immune Cell Types in Young vs. Old",
    x = "Age Group",
    y = "Proportion",
    fill = "Immune Cell Types"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    # Explicitly remove gridlines and background elements if using theme_minimal
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Display the plot
print(stacked_barplot)
####

ggsave("Results/codex_immune_cell_types_old_vs_young(4).png", plot = stacked_barplot, width = 5, height = 5, dpi = 150)

# Save the plot
ggsave("immune_cell_types_old_vs_young.png", plot = stacked_barplot, width = 10, height = 6, dpi = 300)
ggsave("skin_cell_types_old_vs_young.png", plot = stacked_barplot, width = 10, height = 6, dpi = 300)


###


### violin plots & heatmap

f <-c("hIL-1alpha", "CXCR3", "IP10", "IL6", "IL8")
f <- c("p21", "p16", "BCL2", "H2A.X")


young_samples <- c("201", "202", "207", "4614", "7503")
old_samples <- c("4687", "4774", "5014", "CAMS18a", "subj1")

# Add a new column categorizing samples as 'young' or 'old'
combined$aged_sample <- ifelse(combined$sample %in% young_samples, 
                          paste0("young_", combined$sample), 
                          paste0("old_", combined$sample))

# Print the modified table
t<- table(combined$sasp, combined$aged_sample)

t<- table(combined$sencelltype, combined$aged_sample)

# Calculate column-wise percentages (including "none")
pct <- sweep(t, 2, colSums(t), FUN = "/") * 100

# Remove the "none" row after calculating percentages
pct_none <- pct[rownames(pct) != "none", ]

# Print the final percentage table
pct_none_df <- as.data.frame.matrix(pct_none)
pct_none_df$Gene <- rownames(pct_none_df)

pct_long <- melt(pct_none_df, id.vars = "Gene", variable.name = "Sample", value.name = "Percentage")

# Categorize as "young" or "old" based on sample names
pct_long$Group <- factor(ifelse(grepl("young", pct_long$Sample), "young", "old"), levels = c("young", "old"))

# Perform t-tests for each gene to calculate p-values
p_values <- pct_long %>%
  group_by(Gene) %>%
  summarise(p_value = t.test(Percentage[Group == "young"], Percentage[Group == "old"], var.equal = FALSE, na.rm = TRUE)$p.value)

# Define significance levels
p_values$significance <- case_when(
  p_values$p_value < 0.001 ~ "***",
  p_values$p_value < 0.01 ~ "**",
  p_values$p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Merge significance back to main data
pct_long <- merge(pct_long, p_values, by = "Gene")

# Function to generate dot plot with box and error bars
plot_gene_expression <- function(gene_name) {
  gene_data <- pct_long %>% filter(Gene == gene_name)
  
  gene_data$Group <- factor(gene_data$Group, levels = c("young", "old"))  # Ensure correct order
  
  ggplot(gene_data, aes(x = Group, y = Percentage, color = Group)) +
    geom_jitter(position = position_jitter(width = 0.2), size = 3) +  # Dots
    scale_color_manual(values = c("young" = "red", "old" = "darkgreen")) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  # Mean line
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +  # Error bars
    geom_text(aes(x = 1.5, y = max(Percentage) + 5, label = significance, color = "black"), size = 6) +  # Significance stars
    labs(title = gene_name, x = NULL, y = paste(gene_name,"positive cells (% total cells)")) +  
    guides(color = "none") +  # Removes legend entirely
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
}

# Generate individual plots
plot_hIL1alpha <- plot_gene_expression("hIL-1alpha")
plot_CXCR3 <- plot_gene_expression("CXCR3")
plot_IL6 <- plot_gene_expression("IL6")

plot_BCL2 <- plot_gene_expression("BCL2")
plot_HMGB1 <- plot_gene_expression("HMGB1")
plot_LaminB1 <- plot_gene_expression("LaminB1")
plot_p16 <- plot_gene_expression("p16")
plot_p21 <- plot_gene_expression("p21")
plot_PCNA <- plot_gene_expression("PCNA")



# Combine plots into one figure
plot_grid(plot_hIL1alpha, plot_CXCR3, plot_IL6, ncol = 3)

plot_grid(plot_BCL2, plot_HMGB1, plot_LaminB1, plot_p16, plot_p21, plot_PCNA, ncol = 6)

# Display final combined plot
print(final_plot)



f <-c("hIL-1alpha", "CXCR3", "IP10", "IL6", "IL8")
f <- c("p21", "p16", "BCL2", "HMGB1", "LaminB1", "PCNA")

### this is for the average expression of genes per sample, useful code!
a <- AverageExpression(combined, assays = "RNA", slot = "data", group.by=c("age", "sample"), features = f)

# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(Matrix)

# Convert sparse matrix to a data frame
df <- as.data.frame(as.matrix(a$RNA))

# Add gene names as a column
df$Gene <- rownames(df)

rna_long <- melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")

# Categorize as "young" or "old" based on sample names
rna_long$Group <- factor(ifelse(grepl("young", rna_long$Sample), "young", "old"), levels = c("young", "old"))

# Remove NA or zero values if needed
rna_long <- rna_long %>% filter(!is.na(Expression), Expression > 0)

# Perform t-tests for each gene to calculate p-values
p_values <- rna_long %>%
  group_by(Gene) %>%
  summarise(p_value = t.test(Expression[Group == "young"], Expression[Group == "old"], var.equal = FALSE, na.rm = TRUE)$p.value)

# Define significance levels
p_values$significance <- case_when(
  p_values$p_value < 0.001 ~ "***",
  p_values$p_value < 0.01 ~ "**",
  p_values$p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Merge significance back to main data
rna_long <- merge(rna_long, p_values, by = "Gene")

# Function to generate dot plot with box and error bars
plot_gene_expression <- function(gene_name) {
  gene_data <- rna_long %>% filter(Gene == gene_name)
  
  gene_data$Group <- factor(gene_data$Group, levels = c("young", "old"))  # Ensure correct order
  
  ggplot(gene_data, aes(x = Group, y = Expression, color = Group)) +
    geom_jitter(position = position_jitter(width = 0.2), size = 3) +  # Dots
    scale_color_manual(values = c("young" = "red", "old" = "darkgreen")) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  # Mean line
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") +  # Error bars
    geom_text(aes(x = 1.5, y = max(Expression) + 5, label = significance), size = 6) +  # Significance stars
    labs(title = gene_name, x = "Group", y = "Expression Level") +  
    guides(color = "none") +  # Removes legend entirely
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
}

# Generate individual plots
plot_hIL1alpha <- plot_gene_expression("hIL-1alpha")
plot_CXCR3 <- plot_gene_expression("CXCR3")
plot_IL6 <- plot_gene_expression("IL6")

plot_BCL2 <- plot_gene_expression("BCL2")
plot_HMGB1 <- plot_gene_expression("HMGB1")
plot_LaminB1 <- plot_gene_expression("LaminB1")
plot_p16 <- plot_gene_expression("p16")
plot_p21 <- plot_gene_expression("p21")
plot_PCNA <- plot_gene_expression("PCNA")



# Combine plots into one figure
plot_grid(plot_hIL1alpha, plot_CXCR3, plot_IL6, ncol = 3)

plot_grid(plot_BCL2, plot_HMGB1, plot_LaminB1, plot_p16, plot_p21, plot_PCNA, ncol = 6)


# Display final combined plot
print(final_plot)




VlnPlot(object = combined, layer='scale.data', features = f, group.by = "age",  pt.size = 0.5) & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)) & labs(x= NULL)
VlnPlot(object = combined, layer='scale.data', features = f, group.by = "age",  pt.size = 0) & stat_compare_means(method = "wilcox.test") & stat_summary(fun = mean, geom='crossbar', width = 0.3, position=position_dodge(0.9)) & labs(x= NULL)
ggsave(filename = "combined_saspvln_numbered.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine

ggsave(filename = "combined_senvln.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine
DoHeatmap(combined, slot='scale.data', features = f, group.by='age') + NoLegend()
ggsave(filename = "combined_senheatmap.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine

####

#### combined dotplots
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

GetMatrixFromSeuratByGroupMulti<- function(obj, features, group1, group2){
  exp_mat<- obj@assays$RNA$data[features, ,drop=FALSE]
  count_mat<- obj@assays$RNA$counts[features,,drop=FALSE ]
  
  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  
  # get the average expression matrix        
  exp_df<- as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene,{{group1}}, {{group2}}) %>%
    summarise(average_expression = mean(expression)) %>%
    # the trick is to make the data wider in columns: cell_type|group
    tidyr::pivot_wider(names_from = c({{group1}}, {{group2}}), 
                       values_from= average_expression,
                       names_sep="|") 
  
  # convert to a matrix
  exp_mat<- exp_df[, -1] %>% as.matrix()
  rownames(exp_mat)<- exp_df$gene
  
  # get percentage of positive cells matrix
  count_df<- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, {{group1}}, {{group2}}) %>%
    summarise(percentage = mean(count >0)) %>%
    tidyr::pivot_wider(names_from = c({{group1}}, {{group2}}), 
                       values_from= percentage,
                       names_sep="|")
  
  percent_mat<- count_df[, -1] %>% as.matrix()
  rownames(percent_mat)<- count_df$gene
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  
  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}

layer_fun2 = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = "gray", fill = NA))
  grid.circle(x=x,y=y,
              r= sqrt(pindex(mat2$percent_mat, i, j)) * unit(4, "mm"),
              gp = gpar(fill = col_fun(pindex(mat2$exp_mat, i, j)), col = NA))}

col_fun<-  circlize::colorRamp2(c(0, 2, 1), c("#FDE725FF", "#238A8DFF", "#440154FF"))
col_fun<-  circlize::colorRamp2(c(0, 2, 5),c("#006699", "#D9D9D9", "#CC0000"))
col_fun<-  circlize::colorRamp2(c(0, 2, 1),c("#006699", "#D9D9D9", "#CC0000"))


lgd_list = list(
  Legend( labels = c(0, 10, 25, 50, 75), title = "percentage",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0  * unit(4, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1)  * unit(4, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(4, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(4, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(4, "mm"),
                                             gp = gpar(fill = "black"))),
          row_gap = unit(2.5, "mm")
  ))

levels(seurat_object$cell_type)

# Reorder the levels
combined$age <- factor(combined$age, levels = c("young", "old"))


f <- c("hIL-1alpha", "CXCR3", "IL6", "IL8", "IP10")
f <- c("p21", "p16", "BCL2", "H2A.X")

mat2<- GetMatrixFromSeuratByGroupMulti(obj= combined, features =f, skincelltype, age)
hp2<- Heatmap(mat2$exp_mat,
              heatmap_legend_param=list(title= "expression"),
              column_title = "genes", 
              col=col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun2,
              row_names_gp = gpar(fontsize = 8),
              border = "black",
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              row_names_side  = "left")


draw(hp2, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm"),
     annotation_legend_side = "right")

ggsave(filename = "combined_oldyoung_skinsaspdot.svg", device = "svg", width=8, height=5, dpi=1000)  # Works fine
ggsave(filename = "combined_oldyoung_imm2saspdot.svg", device = "svg", width=8, height=5, dpi=1000)  # Works fine
800 by 400

####
seurat <- subset(combined, subset=age== 'old')

Idents(seurat) <- 'skincelltype'
Idents(seurat) <- 'immcelltype2'


f <- c("hIL-1alpha", "CXCR3", "IP10", "IL6", "IL8")

f <- c("p21", "p16", "BCL2", "H2A.X")

f <- c("PCNA", "Ki67")

DotPlot(seurat, features = f, dot.scale=10, dot.min = 0.01 , scale=TRUE)+scale_colour_gradient2(low="#006699", mid="#D9D9D9", high="#CC0000")+theme(axis.text.y=element_text(size=14), axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.line = element_blank())+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+scale_y_discrete(limits=rev)

ggsave(filename = "combined_old_imm2saspdot.svg", device = "svg", width=8, height=5, dpi=1000)  # Works fine


####

Idents(o4687) <- 'seurat_clusters'
DoHeatmap(o4687, slot='scale.data', group.by='age') + NoLegend()


###
library(cowplot)
remotes::install_github("genomics-kl/seurathelpeR")
library(seurathelpeR)
oseurat <- subset(combined, subset=age== 'old')
yseurat <- subset(combined, subset=age== 'young')

+ labs(title = 'Expression of proliferation markers in p16+ cells in old samples')

f <- c("PCNA", "LaminB1") 

add_title_ggplot(VlnPlot(object = oseurat, layer='scale.data', features = f, group.by = "p16",  pt.size = 0)  & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)),
                 'Expression of proliferation markers in p16+ cells in old samples')
ggsave(filename = "old_prolifvln_p16.svg", device = "svg", width=10, height=5, dpi=1000)  # Works fine

add_title_ggplot(VlnPlot(object = combined_1, layer='scale.data', features = f, group.by = "p16",  pt.size = 0)  & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)),
                 'Expression of proliferation markers in p16+ cells in combined samples')
ggsave(filename = "combined_prolifvln_p16.svg", device = "svg", width=10, height=5, dpi=1000)

add_title_ggplot(VlnPlot(object = oseurat, layer='scale.data', features = f, group.by = "p21",  pt.size = 0)  & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)),
                 'Expression of proliferation markers in p21+ cells in old samples')
ggsave(filename = "old_prolifvln_p21.svg", device = "svg", width=10, height=5, dpi=1000)  # Works fine

add_title_ggplot(VlnPlot(object = combined_1, layer='scale.data', features = f, group.by = "p21",  pt.size = 0)  & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)),
                 'Expression of proliferation markers in p21+ cells in combined samples')
ggsave(filename = "combined_prolifvln_p21.svg", device = "svg", width=10, height=5, dpi=1000)

add_title_ggplot(VlnPlot(object = combined_1, layer='scale.data', features = f, group.by = "age",  pt.size = 0)  & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)),
                 'Expression of proliferation markers in young and old cells in combined samples')
ggsave(filename = "combined_prolifvln_age.svg", device = "svg", width=10, height=5, dpi=1000)

f <- c("PCNA","HMGB1", "CXCR3", "IP10", "IL6", "IL8")

add_title_ggplot(VlnPlot(object = combined, layer='scale.data', features = f, group.by = "p16",  pt.size = 0)  & stat_compare_means(method = "t.test", label = "p.signif", label.y.npc = 'top', label.x.npc = 'center') & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)) & labs(y = "Intensity", x = element_blank()),
                 'Intensity of SASP markers in p16+ cells in combined samples') 
ggsave(filename = "combined_saspvln_p16.svg", device = "svg", width=10, height=8, dpi=500)

add_title_ggplot(VlnPlot(object = combined, layer='scale.data', features = f, group.by = "p21",  pt.size = 0)  & stat_compare_means(method = "t.test", label = "p.signif", label.y.npc = 'top', label.x.npc = 'center') & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)) & labs(y = "Intensity", x = element_blank()),
                 'Intensity of SASP markers in p21+ cells in combined samples')
ggsave(filename = "combined_saspvln_p21.svg", device = "svg", width=10, height=8, dpi=500)


VlnPlot(object = combined_1, layer='data', features = f, group.by = "p16",  pt.size = 0) & stat_compare_means(method = "wilcox.test", label = "p.signif") & stat_summary(fun = mean, geom='crossbar', width = 0.5, position=position_dodge(0.9)) & labs(x= NULL)
VlnPlot(object = combined, layer='scale.data', features = f, group.by = "age",  pt.size = 0) & stat_compare_means(method = "wilcox.test") & stat_summary(fun = mean, geom='crossbar', width = 0.3, position=position_dodge(0.9)) & labs(x= NULL)
ggsave(filename = "combined_saspvln_numbered.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine

ggsave(filename = "combined_senvln.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine
DoHeatmap(combined, slot='scale.data', features = f, group.by='age') + NoLegend()
ggsave(filename = "combined_senheatmap.svg", device = "svg", width=10, height=10, dpi=1000)  # Works fine


### histograms for p16, p21 expression in young vs old
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
sen <- c('rna_p21', 'rna_p16')
Idents(combined) <- 'age'
RidgePlot(combined, features = sen, slot = 'scale.data')
RidgePlot(combined, features = 'rna_p21', slot = 'data', group.by = "age", log = TRUE, same.y.lims = 3,sort = TRUE) + 
  scale_x_continuous(trans = "log1p")


library(ggplot2)
library(ggridges)
library(ggpubr)

# Extract p21 expression and metadata
df <- FetchData(combined, vars = c("rna_p21", "rna_p16", "age"))
colnames(df) <- c("p21", "p16", "Age")

# Ridge plot with manual bandwidth
ggplot(df, aes(x = p16, y = Age, fill = Age)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.1) + 
  scale_x_continuous(trans = "log1p") + 
  theme_minimal()


### Cell type-specific expression of p16/p21 by age
library(scCustomize)
Idents(combined) <- 'skincelltype'
f <- c('rna_p16', 'rna_p21')
VlnPlot(object = combined, features = f, split.by = 'age', layer='data')

Stacked_VlnPlot(seurat_object = combined, features = "rna_p16", split.by = 'age',
                x_lab_rotate = TRUE)

# Set cell type as the identity class
Idents(combined) <- 'skincelltype'

# Define genes of interest and clean names
features <- c('rna_p16', 'rna_p21')
features_clean <- gsub("^rna_", "", features)  # Remove 'rna_' prefix

BiocManager::install("dittoSeq")

library(dittoSeq)
dittoBarPlot(combined, 'p16', group.by = "age")

# Generate dittoBarPlots for each cell type for 'p16'
plots <- lapply(setdiff(levels(combined), "none"), function(celltype) {
  dittoBarPlot(
    subset(combined, idents = celltype),
    var = "p21",
    group.by = "age",
    main = paste0(celltype, " - p16"),
    ylab = "Intensity",
    color.panel = c("young" = "#006699", "old" = "#CC0000")
  ) +
    theme(axis.title.x = element_blank())

  })



# Create violin plots for each cell type and gene
plots <- lapply(setdiff(levels(combined), "none"), function(celltype) {
  lapply(seq_along(features), function(i) {
    feature <- features[i]
    feature_clean <- features_clean[i]  # Corresponding clean name
    
    VlnPlot(object = subset(combined, idents = celltype),
            pt.size = 0.5,
            features = feature, 
            group.by = 'age', 
            layer = 'scale.data') &
      
      theme(axis.title.x = element_blank()) &
      stat_compare_means(method = "t.test", 
                         label = "p.signif", 
                         label.y.npc = 'top',
                         label.x.npc = 'center', na.rm = TRUE) &  # Adjust height
      stat_summary(fun = mean, geom = 'crossbar', width = 0.5, position = position_dodge(0.9)) &
      labs(fill = "age") &
      scale_fill_manual(values = c("young" = "#006699", "old" = "#CC0000")) &
      ggtitle(paste0(celltype, " - ", feature_clean)) &  # Title is just "p16" or "p21"
      ylab("Intensity")  # Change y-axis title
  })
})

# Flatten the nested list and combine all plots
plots <- do.call(c, plots)  # Flatten list to get individual plots
combined_plot <- plot_grid(plotlist = plots, ncol = 4)  # Arrange in grid

# Display the plot
print(combined_plot)


# Save the plot
ggsave("combined_violin_plot.png", plot = combined_plot, width = 11, height = 10, dpi = 150)

ggsave("celltype_p16_p21.png", combined_plot, width = 30, height = 30, dpi = 300)


### Stacked violin plot

# Function to modify individual violin plots
modify_vlnplot_by_celltype <- function(obj, 
                                       feature, 
                                       celltype,
                                       pt.size = 0, 
                                       group.by = "age",
                                       plot.margin = unit(c(0, 0, -0.75, 0), "cm"),
                                       ...) {
  
  obj_subset <- subset(obj, subset = skincelltype == celltype)
  
  p <- VlnPlot(obj_subset, 
               features = feature, 
               pt.size = pt.size, 
               group.by = group.by,
               ...) +
    xlab("") + ylab(celltype) + ggtitle("") +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = rel(1), angle = 0),
      axis.text.y = element_text(size = rel(1)),
      plot.margin = plot.margin
    ) +
    stat_compare_means(method = "t.test", 
                       label = "p.signif", 
                       label.y.npc = 'top',
                       label.x.npc = 'center', 
                       na.rm = TRUE) +
    stat_summary(fun = mean, 
                 geom = 'crossbar', 
                 width = 0.5, 
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = c("young" = "#006699", "old" = "#CC0000"))
  
  return(p)
}

# Extract y-axis max
extract_max <- function(p) {
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

# Stacked violin plot across skin cell types
StackedVlnPlotBySkinCellType <- function(obj, 
                                         feature, 
                                         pt.size = 0.5, 
                                         group.by = "age",
                                         plot.margin = unit(c(-100, 0, -5, 0), "cm"),
                                         ...) {
  title <- gsub("^rna_", "", feature)
  
  skin_types <- unique(obj$skincelltype)
  skin_types <- setdiff(skin_types, "none")
  
  plot_list <- purrr::map(skin_types, function(ct) 
    modify_vlnplot_by_celltype(obj = obj, feature = feature, celltype = ct, pt.size = pt.size, group.by = group.by, plot.margin =plot.margin, ...)
  )
  
  # Restore x-axis in bottom plot
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  
  # Adjust y-axis to only show max value
  # ymaxs <- purrr::map_dbl(plot_list, extract_max)
  # plot_list <- purrr::map2(plot_list, ymaxs, function(x, y) 
  #   x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y)
  # )
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)  +
    patchwork::plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ))
  return(p)
}

combined_plot <- plot_grid(plotlist = plots, ncol = 4)  # Arrange in grid

###

combined_plot <- (StackedVlnPlotBySkinCellType(obj = combined, feature = "rna_p21"))
ggsave("combined_p21_codex.png", plot = combined_plot, width = 5, height = 20, dpi = 150)

### Proportions of senescent immune cells

CD45 <- combined@assays$RNA$counts["CD45",] > 0.0

p21 <- combined@assays$RNA$counts["p21",] > 0.0
p16 <- combined@assays$RNA$counts["p16",] > 0.0

PCNA <- combined@assays$RNA$counts["PCNA",] > 0.0

IL6 <- combined@assays$RNA$counts["IL6",] > 0.0
IL8 <- combined@assays$RNA$counts["IL8",] > 0.0
hILa <- combined@assays$RNA$counts["hIL-1alpha",] > 0.0
IP10 <- combined@assays$RNA$counts["IP10",] > 0.0
CXCR3 <- combined@assays$RNA$counts["CXCR3",] > 0.0

combined@meta.data$senimmcells <- case_when(
  CD45 & (p21 | p16) &(!PCNA) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+ CD45+",
  CD45 ~ "CD45+",
  (p21 | p16) &(!PCNA) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+",
  TRUE ~ 'none'
)

# Subset for young and old
combined$age <- factor(combined$age, levels = c("young", "old"))

immsen3 <- as.data.frame(table(combined$age, combined$senimmcells))

cells <- as.data.frame(immsen3 %>%
                         group_by(Var1) %>%
                         mutate(Freq = Freq / sum(Freq)))
#cells <- cells[cells$Var2 == 'SnC+ CD45+',]
cells <- cells[cells$Var2 != 'none',]


cells <- cells %>%
  group_by(Var1) %>%
  mutate(Percentage = round(Freq * 100, 1))  # Convert proportion to percentage

library(ggplot2)




### THE BLACK BORDER & NUMBERS
# Create stacked bar plot with black borders and percentage labels
stacked_barplot <- ggplot(cells, aes(x = Var1, y = Freq, fill = Var2))+
  # 1. Change the color palette to a professional, muted scheme
  scale_fill_brewer(palette = "RdGy") + # "Set2" is a good choice for discrete, muted data +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black borders
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") +  # Percentage labels inside bars
  theme_minimal() +
  labs(
    title = "Proportion of Immune Cell Types in Young vs. Old",
    x = "Age Group",
    y = "Proportion",
    fill = "Cell Types"
  ) +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    # Explicitly remove gridlines and background elements if using theme_minimal
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Display the plot
print(stacked_barplot)
ggsave("combined_sncimm.png", plot = stacked_barplot, width = 5, height = 5, dpi = 150)

###

### spatial plots of senescent immune cells 

o5104
osubj1
ocams
o4687
o4774

CD45 <- o5104@assays$RNA$counts["CD45",] > 0.0

p21 <- o5104@assays$RNA$counts["p21",] > 0.0
p16 <- o5104@assays$RNA$counts["p16",] > 0.0

PCNA <- o5104@assays$RNA$counts["PCNA",] > 0.0
Ki67 <- o5104@assays$RNA$counts["Ki67",] > 0.0

IL6 <- o5104@assays$RNA$counts["IL6",] > 0.0
IL8 <- o5104@assays$RNA$counts["IL8",] > 0.0
hILa <- o5104@assays$RNA$counts["hIL-1alpha",] > 0.0
IP10 <- o5104@assays$RNA$counts["IP10",] > 0.0
CXCR3 <- o5104@assays$RNA$counts["CXCR3",] > 0.0

o5104@meta.data$senimmcells <- case_when(
  #CD45 & (p21 | p16) &(!PCNA) & (!Ki67) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+ CD45+",
  
  (p21 | p16) &(!PCNA) & (!Ki67) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+",
  CD45 ~ "SnC- CD45+",
  TRUE ~ 'none'
)


# 1. Access the original UMAP coordinate matrix
umap_coords <- o5104@reductions$umap@cell.embeddings

# 2. Store the original X and Y columns temporarily
original_X <- umap_coords[, 1]
original_Y <- umap_coords[, 2]

# 3. Perform the simplified transformation (Y becomes new X, X becomes new Y)
# The first column ([ , 1]) is typically X, and the second ([ , 2]) is Y.
umap_coords[, 1] <- original_Y # New X-coordinate is the original Y
umap_coords[, 2] <- original_X  # New Y-coordinate is the original X

umap_coords[, 1] <- umap_coords[, 1] * 5000
umap_coords[, 2] <- umap_coords[, 2] * 20

# 4. Update the object with the new, transformed coordinates
o5104@reductions$umap@cell.embeddings <- umap_coords

o5104@reductions$umap@cell.embeddings <- as.matrix(umap_coords)

Idents(o5104) <- 'senimmcells'
DimPlot(o5104,reduction='umap', pt.size=0.5, order=c("SnC+ CD45+", 'SnC+','SnC- CD45+', 'none'), cols = c("#add8e6", "#31476e", "red", "purple")) + coord_fixed(ratio = 1/200)  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "Results/o5104_senimmcells_2.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

DimPlot(o5104,reduction='umap', pt.size=0.5, order=c('CD45+', 'SnC+','none'), cols = c("white", "white", "#31476e")) + coord_fixed(ratio = 1/100)  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "Results/o5104_senimmcells_cd45.svg", device = "svg", width=8, height=8, dpi=1000)  # Works fine

Idents(o5104) <- 'senimmcells'
DimPlot(o5104,reduction='umap', pt.size=0.5, order=c('SnC+','SnC- CD45+', 'none'), cols = c("#add8e6", "black", "red")) + coord_fixed(ratio = 1/200)  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "Results/o5104_senimmcells_2.png", device = "png", width=8, height=8, dpi=1000)  # Works fine



CD45 <- osubj1@assays$RNA$counts["CD45",] > 0.0

p21 <- osubj1@assays$RNA$counts["p21",] > 0.0
p16 <- osubj1@assays$RNA$counts["p16",] > 0.0

PCNA <- osubj1@assays$RNA$counts["PCNA",] > 0.0

IL6 <- osubj1@assays$RNA$counts["IL6",] > 0.0
IL8 <- osubj1@assays$RNA$counts["IL8",] > 0.0
hILa <- osubj1@assays$RNA$counts["hIL-1alpha",] > 0.0
IP10 <- osubj1@assays$RNA$counts["IP10",] > 0.0
CXCR3 <- osubj1@assays$RNA$counts["CXCR3",] > 0.0

osubj1@meta.data$senimmcells <- case_when(
  CD45 ~ "CD45+",
  (p21 | p16) &(!PCNA) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+",
  TRUE ~ 'none'
)

Idents(osubj1) <- 'senimmcells'
DimPlot(osubj1,reduction='umap', pt.size=1, order=c('SnC+','CD45+', 'none'), cols = c("gray", "blue", "red")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "osubj1_senimmcells.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

CD45 <- ocams@assays$RNA$counts["CD45",] > 0.0

p21 <- ocams@assays$RNA$counts["p21",] > 0.0
p16 <- ocams@assays$RNA$counts["p16",] > 0.0

PCNA <- ocams@assays$RNA$counts["PCNA",] > 0.0

IL6 <- ocams@assays$RNA$counts["IL6",] > 0.0
#IL8 <- ocams@assays$RNA$counts["IL8",] > 0.0
hILa <- ocams@assays$RNA$counts["hIL-1alpha",] > 0.0
CXCR3 <- ocams@assays$RNA$counts["CXCR3",] > 0.0

ocams@meta.data$senimmcells <- case_when(
  CD45 ~ "CD45+",
  (p21 | p16) &(!PCNA) & (IL6|hILa|CXCR3) ~ "SnC+",
  TRUE ~ 'none'
)

Idents(ocams) <- 'senimmcells'
DimPlot(ocams,reduction='umap', pt.size=1, order=c('SnC+','CD45+', 'none'), cols = c("gray", "blue", "red")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "ocams_senimmcells.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

CD45 <- o4687@assays$RNA$counts["CD45",] > 0.0

p21 <- o4687@assays$RNA$counts["p21",] > 0.0
p16 <- o4687@assays$RNA$counts["p16",] > 0.0

PCNA <- o4687@assays$RNA$counts["PCNA",] > 0.0

IL6 <- o4687@assays$RNA$counts["IL6",] > 0.0
IL8 <- o4687@assays$RNA$counts["IL8",] > 0.0
hILa <- o4687@assays$RNA$counts["hIL-1alpha",] > 0.0
IP10 <- o4687@assays$RNA$counts["IP10",] > 0.0
CXCR3 <- o4687@assays$RNA$counts["CXCR3",] > 0.0

o4687@meta.data$senimmcells <- case_when(
  CD45 ~ "CD45+",
  (p21 | p16) &(!PCNA) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+",
  TRUE ~ 'none'
)
Idents(o4687) <- 'senimmcells'
DimPlot(o4687,reduction='umap', pt.size=1, order=c('SnC+','CD45+', 'none'), cols = c("gray", "blue", "red")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "o4687_senimmcells.png", device = "png", width=8, height=8, dpi=1000)  # Works fine


CD45 <- o4774@assays$RNA$counts["CD45",] > 0.0

p21 <- o4774@assays$RNA$counts["p21",] > 0.0
p16 <- o4774@assays$RNA$counts["p16",] > 0.0

PCNA <- o4774@assays$RNA$counts["PCNA",] > 0.0

IL6 <- o4774@assays$RNA$counts["IL6",] > 0.0
IL8 <- o4774@assays$RNA$counts["IL8",] > 0.0
hILa <- o4774@assays$RNA$counts["hIL-1alpha",] > 0.0
IP10 <- o4774@assays$RNA$counts["IP10",] > 0.0
CXCR3 <- o4774@assays$RNA$counts["CXCR3",] > 0.0

o4774@meta.data$senimmcells <- case_when(
  CD45 ~ "CD45+",
  (p21 | p16) &(!PCNA) & (IL6|IL8|hILa|IP10|CXCR3) ~ "SnC+",
  TRUE ~ 'none'
)

Idents(o4774) <- 'senimmcells'
DimPlot(o4774,reduction='umap', pt.size=1, order=c('SnC+','CD45+', 'none'), cols = c("#add8e6", "#31476e", "red")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename = "Data/o4774_senimmcells_1.png", device = "png", width=8, height=8, dpi=1000)  # Works fine

###

### pie charts
p16 <- combined@assays$RNA$counts["p16",] > 0.0
p21 <- combined@assays$RNA$counts["p21",] > 0.0

combined@meta.data$sencelltype3 <- case_when(
  p16  & p21 ~ "p16/p21",
  p16 ~ "p16",
  #p16 ~ "p16",
  p21 ~ "p21",
  TRUE ~ 'none'
)

CD4 <- combined@assays$RNA$counts["CD4",] > 0.0
CD68 <- combined@assays$RNA$counts["CD68",] > 0.0
CD8 <- combined@assays$RNA$counts["CD8",] > 0.0
FOXP3 <- combined@assays$RNA$counts["FOXP3",] > 0.0
#CD11c <- combined@assays$RNA$counts["CD11c",] > 0.0
CD14 <- combined@assays$RNA$counts["CD14",] > 0.0
CXCR3 <- combined@assays$RNA$counts["CXCR3",] > 0.0
CD45 <- combined@assays$RNA$counts["CD45",] > 0.0


combined@meta.data$immcelltype2 <- case_when(
  CD8 ~ "CD8",
  CD4  & FOXP3 ~ "CD4+ FOXP3+",
  CD4  & CXCR3 ~ "CD4+ CXCR3+",
  CD14 ~ "CD14+",
 # CD11c ~ "CD11c+",
  CD45 ~ "CD45+",
  CD68 ~ "CD68+",
  CD4 ~ "CD4+",
  TRUE ~ 'none'
)


library(plotly)
library(ggplot2)
library(dplyr)
library(scales)

smallercrop2_seurat <- subset(combined, subset=age== 'old')
smallercrop2_seurat <- subset(combined, subset=age== 'young')
smallercrop2_seurat <- subset(smallercrop2_seurat, subset=sample != '7503')

immsen <- as.data.frame(table(smallercrop2_seurat$immcelltype2, smallercrop2_seurat$sencelltype3))
your_data <- as.data.frame(immsen %>%
                             group_by(Var1) %>%
                             mutate(Freq = Freq / sum(Freq)))
colnames(your_data) <- c('cell_types', 'Markers', 'Freq')
your_data <- your_data[your_data$cell_types != "none", ]

my_manual_colors <- c(
  "p16/p21" = "#3156a6",
  "p21" = "#cd8729",
  "p16" = "#CC0000",
  "none" = "gray"
)

# 1. Calculate the position for the labels and create the formatted text.
#    This version only creates a label for slices larger than 5% to avoid clutter.
your_data <- your_data %>%
  group_by(cell_types) %>%
  arrange(Markers) %>% # Arrange to ensure consistent order for cumsum
  mutate(
    label_position = cumsum(Freq) - 0.5 * Freq,
    percent_label = if_else(Freq > 0.05, percent(Freq, accuracy = 1), "")
  )

# Create the formatted percentage label. The manual position calculation is no longer needed.
your_data <- your_data %>%
  mutate(
    percent_label = if_else(Freq > 0.05, percent(Freq, accuracy = 1), "")
  )


png("young_immsen.png", res = 1200, width = 10, height = 8, units='in', pointsize=4)

ggplot(data=your_data, aes(x=" ", y=Freq, group=Markers, fill=Markers)) +
  geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.5) + # Added an outline
  coord_polar("y", start=0) +
  
  # --- ADD THIS LINE TO PLOT THE LABELS ---
  geom_text(
    aes(label = percent_label), 
    color = "black", 
    size = 2.5,
    position = position_stack(vjust = 0.5) # <-- This is the key fix
  ) +
  
  facet_grid(.~ cell_types) +
  theme_void() +
  scale_fill_manual(values = my_manual_colors)

dev.off()


###

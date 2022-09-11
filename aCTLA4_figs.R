##############Libraries####
library(Matrix)
library(reshape)
library(readxl)
library(ggplot2)
library(devtools)
library(circlize)
library(ggpubr)
library(plyr)
library(dplyr)
library(scales)
library(viridis)
library(grid)
library(Seurat)
library(mgcv)
library(MASS)
library(colorRamps)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpmisc)
library(cowplot)
library(M3C)
library(ConsensusClusterPlus)
library(clusterProfiler)
library(rstatix)
library(fgsea)
library(org.Mm.eg.db)
library(GO.db)
library(limma)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


##############load data fig 1-4####
projdir = "/home/labs/amit/tomerlan/UCL_project"
id = "final_1"
workdir = paste0("/home/labs/amit/tomerlan/UCL_project/", id)
setwd(workdir)
seurat = readRDS(file = "seurat.rds")
meta_plate = read.delim("metadata.txt")
# setwd(paste0(projdir, "/figures"))


assay = "integrated2"
DefaultAssay(seurat) = assay
seurat@active.ident = as.factor(seurat$integrated2_nn_res.2)
cluster_ids = read.csv(paste0(workdir, "/", assay, "_annot.csv"), row.names = 1); 
cluster_annot = cluster_ids[,1]; names(cluster_annot) = rownames(cluster_ids)
cluster_color = cluster_ids[,2]; names(cluster_color) = as.character(cluster_annot)
# seurat@active.ident = cluster_annot[as.character(seurat@active.ident)]
seurat = RenameIdents(seurat, cluster_annot)

unique(seurat@active.ident)
clusters_lymphoid =  c("Treg_Foxp3", "CD4_Cxcr6",  "CD4_Ifit1", "CD8_Xcl1", "CD8_Pdcd1", "CD8_Ccl4", "CD8_Ccl5", "CD8_Lef1", "NK_Nrc1") # "B_cells
clusters_myeloid = c("Mreg_Hmox1", "TAM_C1qa", "TAM_Arg1", "TAM_Ifit1", "TAM_Ccl7", "Mon_Arg1", "Mon_Ifit1", "Mon_Ccl5", "Mon_Lyz2",  "Mon_H2-Ab1", "cDC2_Cd209a", "migDC_Ccr7", "Neut_S100a8") #, "Basophil", "Mon_UI"
cluster_monocyte = c("Mon_Arg1", "Mon_Ifit1","Mon_Lyz2", "Mon_H2-Ab1", "Mon_Ccl5")
cluster_macrophage = c("TAM_C1qa", "TAM_Arg1",  "TAM_Ifit1", "TAM_Ccl7")

stress_genes= c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes =  unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))), "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"))
IFN_genes = unique(c(grep("Irf", rownames(seurat@assays$RNA@counts), v = T), grep("Ifi", rownames(rownames(seurat@assays$RNA@counts)), v = T), "Cd2","Cd3d", "Cmpk2","Cxcl10", "Isg15","Isg20","Oasl1" ,"Phf11b" ,"Plac8", "Rsad2","Rtp4","Sdf4", "Slfn4","Tnfrsf18" ,"Tnfrsf4","Tnfrsf9","Trex1","Usp18","Xaf1"))
ccl_genes = grep("Ccl", rownames(seurat@assays$RNA@counts), v = T)
MHC_genes =  c(grep("^H2-", rownames(seurat@assays$RNA@counts), v = T), "Cd74", "B2m")
comp_genes =  c(grep("^C1q", rownames(seurat@assays$RNA@counts), v = T))
ig_genes = c(grep("^Igj|^Igh|^Igk|^Igl", rownames(seurat@assays$RNA@counts), v = T))
hb_genes = c(grep("^Hba|^Hbb", rownames(seurat@assays$RNA@counts), v = T))
markers = read.csv(paste0(workdir, "/", assay, "_markers.csv"))[, 1:2]
LR_pairs = read.table(paste0(projdir, "/mouse_lr_pair.txt"), header = T)


##############fig 1b####
stats = as.data.frame(read_xlsx("/home/labs/amit/tomerlan/UCL_project/tumor.xlsx", sheet = 3, col_names = T))
class(stats$X) = "numeric"
class(stats$Y) = "numeric"
class(stats$vol) = "numeric"
class(stats$vol_fixed) = "numeric"
stats$treatment = factor(stats$treatment, levels = c("Ctrl", "aCTLA4-m1", "aCTLA4-m2a"))

stats2 = as.data.frame(read_xlsx("/home/labs/amit/tomerlan/UCL_project/tumor.xlsx", sheet = 10, col_names = T))
class(stats2$day) = "numeric"
class(stats2$vol) = "numeric"
stats2$treatment = mapvalues(stats2$treatment, from = c("UT", "m1", "m2a"), to = c("Ctrl", "aCTLA4-m1", "aCTLA4-m2a"))
stats2$rescued = as.character(stats2$rescued)

stats = cbind(rbind(stats[,1:4], stats2[,1:4]), Model = c(rep("MCA205", nrow(stats)), rep("MC38", nrow(stats2))))
stats$Model = factor(stats$Model, levels = c("MCA205", "MC38"))
ggplot(stats, mapping = aes(x = day, y = vol, group = paste(mouse, treatment))) +
  facet_grid(.~ treatment) +
  geom_line(aes(linetype =Model)) +
  # geom_path(arrow = arrow()) +
  # geom_point(aes(color = treatment), size = 2) +
  scale_color_manual(values = c("grey90", "grey60", "grey20")) +
  xlab("Days after tumor injection") +
  ylab("Tumor volume [mm^3]") +
  scale_x_continuous(breaks = unique(stats$day)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size= 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))
ggsave(file = "fig1_tumor_growth4.png", width = 300, height = 120, units = "mm", dpi = 150, bg = "transparent")
ggsave(file = "fig1_tumor_growth4.pdf", width = 300, height = 120, units = "mm", dpi = 150, bg = "transparent")


summary(aov(vol ~ day*treatment, data = stats[stats$treatment %in% c("Ctrl", "aCTLA4-m2a"),]))
summary(aov(vol ~ day*treatment, data = stats[stats$treatment %in% c("aCTLA4-m1", "aCTLA4-m2a"),]))
summary(aov(vol ~ day*treatment, data = stats[stats$treatment %in% c("Ctrl", "aCTLA4-m1"),]))


##############fig 1c#####
dot_days = c("d14")
dot_exp = c(1, 2)
dot_treat = c("UT","m1", "m2a")
dot_cells =  c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_lymphoid & seurat$Staining.panel.name %in% "T" & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Treatment %in% dot_treat]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid & seurat$Staining.panel.name %in% c("Innate", "Broad") & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Treatment %in% dot_treat])); length(dot_cells)

bad_features_pca = c(stress_genes, cc_genes, ig_genes, MHC_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
seurat.reduction1 = RunPCA(seurat[, dot_cells], assay = assay, features = setdiff(seurat@assays[[assay]]@var.features, bad_features_pca))
seurat.reduction1 = RunUMAP(seurat.reduction1[, dot_cells], assay = assay, dims = 1:30, a = 0.5, b = 2)

cluster_num = mapvalues(c(clusters_lymphoid, clusters_myeloid), from = c(clusters_lymphoid, clusters_myeloid), to = paste0(c(clusters_lymphoid, clusters_myeloid), " (", 1:length(c(clusters_lymphoid, clusters_myeloid)),")"))
cluster_num_color = cluster_color[c(clusters_lymphoid, clusters_myeloid)]; names(cluster_num_color) = cluster_num
ggumap = data.frame(seurat.reduction1@reductions$umap@cell.embeddings, Clust = seurat.reduction1@active.ident, Clust_no = mapvalues(seurat.reduction1@active.ident, from = c(clusters_lymphoid, clusters_myeloid), to = 1:length(c(clusters_lymphoid, clusters_myeloid))), Treatment = seurat.reduction1$Treatment, Day = seurat.reduction1$Day.of.tumor.harvest)
ggcenters = aggregate(ggumap[, 1:2], by = list(ggumap$Clust_no), FUN = median)
ggumap$Treatment = mapvalues(factor(ggumap$Treatment , levels = c("UT", "m1", "m2a")), from = c("UT", "m1", "m2a"), to = c("Control", "aCTLA4-m1", "aCTLA4-m2a"))
ggumap$Clust = factor(ggumap$Clust, levels = c(clusters_lymphoid, clusters_myeloid))

#clusters
ggplot(ggumap) +
  theme(legend.key.size = unit(0.4, "cm")) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = factor(paste0(Clust, " (", Clust_no,")"), levels = cluster_num)), shape = 21, stroke = 0.01, size = 1, alpha = 1) +
  scale_fill_manual(values = cluster_num_color) +
  geom_text(data = ggcenters, aes(x = UMAP_1, y = UMAP_2, label = Group.1), size = 3, fontface = "bold") +
  guides(fill = guide_legend(title = "Cell population", override.aes = list(size = 3))) +
  theme_void() +
  # lims(x = c( min(ggumap[,1:2]), max(ggumap[,1:2])), y = c(min(ggumap[,1:2]), max(ggumap[,1:2]))) +
  theme(panel.background = element_rect(fill = "transparent", size = 0), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom")
ggsave(file = paste0("fig1_umap_cluster.png"), width = 100, height = 60, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig1_umap_cluster.pdf"), width = 100, height = 60, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')

##treatment cell density
umap_x_limit = -5

dot_panel = "T"
dot_clusters = clusters_lymphoid
dot_cells_lymphoid =  intersect(names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]), rownames(ggumap[ggumap$UMAP_1 < umap_x_limit,])); length(dot_cells_lymphoid)

dot_panel = c("Innate", "Broad")
dot_clusters = clusters_myeloid
dot_cells_myeloid =  intersect(names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]), rownames(ggumap[ggumap$UMAP_1 > umap_x_limit,])); length(dot_cells_myeloid)

ggumap = data.frame(ggumap, density = rep(0, nrow(ggumap)))

ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "UT"))), "density"] =  get_density(ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "UT"))),"UMAP_1"], ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "UT"))),"UMAP_2"], n = 100)
ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m1"))), "density"] =  get_density(ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m1"))),"UMAP_1"], ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m1"))),"UMAP_2"], n = 100)
ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m2a"))), "density"] =  0.5*get_density(ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m2a"))),"UMAP_1"], ggumap[intersect(dot_cells_lymphoid, names(which(seurat$Treatment == "m2a"))),"UMAP_2"], n = 100)
ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "UT"))), "density"] =  get_density(ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "UT"))),"UMAP_1"], ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "UT"))),"UMAP_2"], n = 100)
ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m1"))), "density"] =  get_density(ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m1"))),"UMAP_1"], ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m1"))),"UMAP_2"], n = 100)
ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m2a"))), "density"] =  get_density(ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m2a"))),"UMAP_1"], ggumap[intersect(dot_cells_myeloid, names(which(seurat$Treatment == "m2a"))),"UMAP_2"], n = 100)
# ggumap[dot_cells_lymphoid, "density_l"] =  get_density(ggumap[dot_cells_lymphoid,"UMAP_1"], ggumap[dot_cells_lymphoid,"UMAP_2"], n = 100)
# ggumap[dot_cells_myeloid, "density_m"] =  get_density(ggumap[dot_cells_myeloid,"UMAP_1"], ggumap[dot_cells_myeloid,"UMAP_2"], n = 100)

ggplot(ggumap) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = density), alpha = 1, size = 1.2,  shape = 21, stroke = 0) +
  # scale_fill_gradientn(colors = matlab.like(100)) +
  # scale_fill_gradientn(colors = heat.colors(100)) +
  scale_fill_gradientn(colors = rev(inferno(100))) +
  # scale_fill_gradientn(colors = wes_palette("Zissou1", 10, type = "continuous")) +
  # scale_fill_viridis(option = "F") +
  facet_wrap(. ~ Treatment, scales = "free", ncol = 1,strip.position = "top") +
  theme_void() +
  # lims(x = c( min(ggumap[,1:2]), max(ggumap[,1:2])), y = c(min(ggumap[,1:2]), max(ggumap[,1:2]))) +
  theme(panel.background = element_rect(fill = "transparent", size = 0), 
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 12),
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent", size =0),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom") +
  guides(colour = F, alpha = F)
ggsave(file = paste0("fig1_umap_treatment_density_cell.png"), width = 100, height = 170, units = "mm", dpi = 300, device = "png", bg = "transparent")
ggsave(file = paste0("fig1_umap_treatment_density_cell.pdf"), width = 100, height = 170, units = "mm", dpi = 300, device = "pdf", bg = "transparent")


##CTLA4 expression
FeaturePlot(seurat.reduction1, features = "Ctla4",reduction = "umap", cols = c("white", "red"), pt.size = 0.1, label = T, repel = T, label.size = 2, slot = "data")

umap_assay = "integrated2"
ggumap$"CTLA4" = CTLA4 = seurat.reduction1@assays[[umap_assay]]@scale.data["Ctla4",]

#clusters
ggplot(ggumap) +
  theme(legend.key.size = unit(0.4, "cm")) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = CTLA4), color = "grey", shape = 21, stroke = 0, size = 1.2, alpha = 1) +
  scale_fill_gradient2(low = "white", high =  "red") +
  geom_text(data = ggcenters, aes(x = UMAP_1, y = UMAP_2, label = Group.1), size = 3, fontface = "bold") +
  theme_void() +
  # lims(x = c( min(ggumap[,1:2]), max(ggumap[,1:2])), y = c(min(ggumap[,1:2]), max(ggumap[,1:2]))) +
  theme(panel.background = element_rect(fill = "transparent", size = 0), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom")
ggsave(file = "fig1_ctla4_umap.png", width = 100, height = 50, units = "mm", dpi = 250, device = "png")
ggsave(file = "fig1_ctla4_umap.pdf", width = 100, height = 50, units = "mm", dpi = 250, device = "pdf")




##############fig 1d#####
unique(seurat@active.ident)
dot_logFC_cut = log2(1.5)
droplet_bad_features = c(stress_genes, cc_genes, ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dotplot_lin = "Lymphoid"
dot_days = c("d14")
dot_exp = c(1, 2)
dotplot_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_lymphoid & seurat$Staining.panel.name %in% "T" & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(dotplot_cells)
pct = 0.1

dotplot_markers = FindAllMarkers(seurat[,dotplot_cells], assay = assay, only.pos = T, min.pct = pct, logfc.threshold = dot_logFC_cut, base = 2, test.use = "wilcox", pseudocount.use = mean(seurat@assays[[assay]]@data), features = setdiff(seurat@assays[[assay]]@var.features, droplet_bad_features), max.cells.per.ident = 1000, random.seed = 1, min.cells.feature = 10)
dotplot_genes = intersect(markers[,dotplot_lin], rownames(seurat@assays[[assay]]@scale.data)); length(dotplot_genes)
dotplot_means = apply(as.matrix(clusters_lymphoid), 1, function(x) rowMeans(seurat@assays[[assay]]@scale.data[dotplot_genes, names(seurat@active.ident[seurat@active.ident %in% x])])); colnames(dotplot_means) = clusters_lymphoid
dotplot_clusts_hc = hclust(dist(t(dotplot_means)))
dotplot_clusts_order = rev(dotplot_clusts_hc$labels[dotplot_clusts_hc$order])
dotplot_genes_order = unique(c(unlist(apply(as.matrix(dotplot_clusts_order), 1, function(x) intersect(head(dotplot_markers[dotplot_markers$cluster == x, "gene"], 100), dotplot_genes))))); length(dotplot_genes_order)
seurat@active.ident = factor(seurat@active.ident, levels = dotplot_clusts_order)
droplet_markers_top = t(aaply(as.matrix(clusters_lymphoid), 1, function(ct) head(dotplot_markers[dotplot_markers$cluster == ct, "gene"], 12))); colnames(droplet_markers_top) = clusters_lymphoid

DotPlot(seurat[,dotplot_cells], features = rev(dotplot_genes_order), assay = assay, cols = c("white", "red"), cluster.idents = F) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "none")
# scale_color_gradientn(values = c(-0.1,0,2), colors = c("blue", "white", "red")) 
ggsave(file = paste0("fig1_dotmap_lymphoid.png"), width = length(clusters_lymphoid)*10, height = 200, units = "mm", dpi = 150, device = "png")
ggsave(file = paste0("fig1_dotmap_lymphoid.pdf"), width = length(clusters_lymphoid)*10, height = 200, units = "mm", dpi = 150, device = "pdf", bg = "transparent")



##############fig 1e#####
unique(seurat@active.ident)
droplet_bad_features = c(stress_genes, cc_genes, ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dot_days = c("d14")
dot_exp = c(1, 2)
dot_assay = "integrated2"
dotplot_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid & seurat$Staining.panel.name %in% "Innate" & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(dotplot_cells)

dotplot_markers = FindAllMarkers(seurat[,dotplot_cells], assay = dot_assay, only.pos = T, min.pct = 0.05, logfc.threshold = 0.25, base = 2, test.use = "wilcox", pseudocount.use = mean(seurat@assays[[assay]]@data), features = setdiff(seurat@assays[[assay]]@var.features, droplet_bad_features), max.cells.per.ident = 1000, random.seed = 1, min.cells.feature = 10)
dotplot_genes = intersect(markers[,"Myeloid"], rownames(seurat@assays[[dot_assay]]@data)); length(dotplot_genes)
dotplot_means = apply(as.matrix(clusters_myeloid), 1, function(x) rowMeans(seurat@assays[[dot_assay]]@data[dotplot_genes, names(seurat@active.ident[seurat@active.ident %in% x])])); colnames(dotplot_means) = clusters_myeloid
dotplot_clusts_hc = hclust(dist(t(dotplot_means)))
dotplot_clusts_order = rev(dotplot_clusts_hc$labels[dotplot_clusts_hc$order])
dotplot_genes_order = unique(c(unlist(apply(as.matrix(dotplot_clusts_order), 1, function(x) intersect(head(dotplot_markers[dotplot_markers$cluster == x, "gene"], 300), dotplot_genes))))); length(dotplot_genes_order)
seurat@active.ident = factor(seurat@active.ident, levels = dotplot_clusts_order)

DotPlot(seurat[,dotplot_cells], features = rev(dotplot_genes_order), scale = T, assay = dot_assay, cols = c("white", "red"), cluster.idents = F) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal") 
ggsave(file = paste0("fig1_dotmap_myeloid.png"), width = length(clusters_myeloid)*9, height = 200, units = "mm", dpi = 150, device = "png")
ggsave(file = paste0("fig1_dotmap_myeloid.pdf"), width = length(clusters_myeloid)*9, height = 200, units = "mm", dpi = 250, device = "pdf",  bg = "transparent")



##############fig 1g#####
dot_panel = "T"
dot_clusters = clusters_lymphoid
dot_days = c("d14")
dot_exp = c(1, 2)
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_plot = melt(dot_sum)
dot_plot_mean = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = mean)
dot_plot_sd = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = sd)
dot_plot = cbind(dot_plot_mean, sd = dot_plot_sd$x, Name = paste(dot_plot_mean$variable))
dot_plot$Treatment = mapvalues(dot_plot$Treatment, from = c("UT", "m1","m2a"), to = c("Ctrl", "aCTLA4-m1","aCTLA4-m2a"))
temp = data.frame(Name = dot_plot$Name, Treatment = dot_plot$Treatment, x = dot_plot$x)
ref = "Ctrl"
temp = temp[temp$Treatment == ref,]
dot_plot = merge(dot_plot, temp, by = "Name")
dot_plot = dot_plot[-which(dot_plot$Treatment.x == ref) ,]

dot_order_temp = dot_plot[dot_plot$Treatment.x %in% c("aCTLA4-m2a"),]
dot_order_values = log2(dot_order_temp$x.x/dot_order_temp$x.y)
dot_order_lymphoid = as.character(dot_order_temp$variable[order(dot_order_values)])

dot_pvalues_m2a = apply(as.matrix(clusters_lymphoid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m2a",clust])$p.value); names(dot_pvalues_m2a) = clusters_lymphoid
dot_pvalues_m1 = apply(as.matrix(clusters_lymphoid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m1",clust])$p.value); names(dot_pvalues_m1) = clusters_lymphoid
dot_plot = arrange(dot_plot, by_group = Treatment.x)
dot_plot = data.frame(dot_plot, pvalue = c(dot_pvalues_m1[as.character(unique(dot_plot$variable))], dot_pvalues_m2a[as.character(unique(dot_plot$variable))]))
dot_plot_lymph = dot_plot
ggplot(dot_plot, aes(x = factor(variable, levels = dot_order_lymphoid), y =  log2(x.x/x.y))) +
  facet_grid(.~Treatment.x) +
  # geom_bar(aes(fill = variable), stat = "identity", position = "dodge2") + #, width = x.y + x.x
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = pvalue, size = (x.y + x.x)/2), shape = 21) +
  scale_fill_gradientn(colours = c("red", "white"), trans = "log10") +
  # scale_fill_manual(values =  cluster_color) +
  theme_classic() +
  coord_flip() +
  ylim(-3,3)+
  labs(y = paste0("log2(X/Ctrl)")) +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", size = 0), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent", size =0),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(size = guide_legend(title  = "%"))
ggsave(file = paste0("fig1_ratios_d14_lymphoid.png"), width = 120, height = length(dot_clusters)*7, units = "mm", dpi = 350, device = "png", bg = 'transparent')
ggsave(file = paste0("fig1_ratios_d14_lymphoid.pdf"), width = 120, height = length(dot_clusters)*8, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')




##############fig 1h-k####
dot_clusters = clusters_lymphoid
dot_days = c("d14")
dot_exp = c(1, 2)
num_cells =  setdiff(names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% "T" & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days])); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)
dot_sum_lymph = dot_sum
dot_show = c("Treg_Foxp3", "CD4_Cxcr6")
ggplot(dot_plot[dot_plot$variable %in% dot_show,], aes(x = Treatment, y = value)) +
  facet_wrap(. ~ variable,  strip.position = 'top', nrow = 1, scales = "free") +
  stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  stat_summary(fun.data=mean_sdl, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + #c("cornflowerblue", "orange", "red")
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("% of lymphoid") 
ggsave(file = paste0("fig1_makeup_l.png"), width = 50, height = 60, units = "mm", dpi = 150, device = "png", bg = 'transparent')
ggsave(file = paste0("fig1_makeup_l.pdf"), width = 50, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/(CD4_Cxcr6 + CD4_Ifit1)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/CD4 ratio") 
ggsave(file = paste0("fig1_CD8-4.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_CD8-4.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2((CD4_Ifit1)/c(CD4_Cxcr6)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD4 Ifn/CD4 Cxcr6 ratio")
ggsave(file = paste0("fig1_CD4_Ifn-Cxcr6.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_CD4_Ifn-Cxcr6.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/c(Treg_Foxp3)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, test.args = list(alternative = "less", var.equal = F, paired=F)) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/Treg ratio")
ggsave(file = paste0("fig1_CD8-Treg.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_CD8-Treg.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/c(Treg_Foxp3)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "wilcox.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, test.args = list(alternative = "less", var.equal = F, paired=F)) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/Treg ratio")
ggsave(file = paste0("figrev1_CD8-Treg.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev1_CD8-Treg.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/(CD4_Cxcr6 + CD4_Ifit1)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "wilcox.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/CD4 ratio") 
ggsave(file = paste0("figrev1_CD8-4.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev1_CD8-4.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


##############fig 1i#####
dot_panel = ("Innate")
dot_clusters = clusters_myeloid
dot_days = c("d14")
dot_exp = c(1, 2)
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_plot = melt(dot_sum)
dot_plot_mean = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = mean)
dot_plot_sd = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = sd)
dot_plot = cbind(dot_plot_mean, sd = dot_plot_sd$x, Name = paste(dot_plot_mean$variable))
dot_plot$Treatment = mapvalues(dot_plot$Treatment, from = c("UT", "m1","m2a"), to = c("Ctrl", "aCTLA4-m1","aCTLA4-m2a"))
temp = data.frame(Name = dot_plot$Name, Treatment = dot_plot$Treatment, x = dot_plot$x)
ref = "Ctrl"
temp = temp[temp$Treatment == ref,]
dot_plot = merge(dot_plot, temp, by = "Name")
dot_plot = dot_plot[-which(dot_plot$Treatment.x == ref) ,]

dot_order_temp = dot_plot[dot_plot$Treatment.x %in% c("aCTLA4-m2a"),]
dot_order_values = log2(dot_order_temp$x.x/dot_order_temp$x.y)
dot_order_myeloid = as.character(dot_order_temp$variable[order(dot_order_values)])

dot_pvalues_m2a = apply(as.matrix(clusters_myeloid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m2a",clust])$p.value); names(dot_pvalues_m2a) = clusters_myeloid
dot_pvalues_m1 = apply(as.matrix(clusters_myeloid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m1",clust])$p.value); names(dot_pvalues_m1) = clusters_myeloid
dot_plot = arrange(dot_plot, by_group = Treatment.x)
dot_plot = data.frame(dot_plot, pvalue = c(dot_pvalues_m1[as.character(unique(dot_plot$variable))], dot_pvalues_m2a[as.character(unique(dot_plot$variable))]))
dot_plot_mye = dot_plot
ggplot(dot_plot, aes(x = factor(variable, levels = dot_order_myeloid), y =  log2(x.x/x.y))) +
  facet_grid(.~Treatment.x) +
  # geom_bar(aes(fill = variable), stat = "identity", position = "dodge2") + #, width = x.y + x.x
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = pvalue, size = (x.y + x.x)/2), shape = 21) +
  scale_fill_gradientn(colours = c("red", "white"), trans = "log10") +
  # scale_fill_manual(values =  cluster_color) +
  theme_classic() +
  coord_flip() +
  ylim(-5, 7)+
  labs(y = paste0("log2(X/Ctrl)")) +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", size = 0), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent", size =0),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(size = guide_legend(title  = "%"))
ggsave(file = paste0("fig1_enrichment_m.png"), width = 140, height = 80, units = "mm", dpi = 350, device = "png", bg = 'transparent')
ggsave(file = paste0("fig1_enrichment_m.pdf"), width = 140, height = 80, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')



##############fig 1m-p####
dot_panel = c("Innate", "Broad")
dot_clusters = clusters_myeloid
dot_days = c("d14")
dot_exp = c(1, 2)
num_cells =  setdiff(names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]), names(seurat$amp_batch_id[seurat$amp_batch_id == "AB3636"])); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)

ggplot(dot_sum, aes(x = Treatment, y = log2((Mon_Ifit1 + 0.5)/(Mon_Arg1 + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Mon Ifn/Mon Arg1 ratio") 
ggsave(file = paste0("fig1_mon_Ifn-Arg1.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_mon_Ifn-Arg1.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2((Mon_Ccl5 + 0.5)/(Mon_Arg1 + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Mon Ccl5/Mon Arg1 ratio") 
ggsave(file = paste0("fig1_mon_Ccl5-Arg1.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_mon_Ccl5-Arg1.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2((TAM_Ifit1 + 0.5)/(TAM_Arg1 + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Mac Ifn/Mac Arg1 ratio)") 
ggsave(file = paste0("fig1_mac_Ifn-Arg1.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("fig1_mac_Ifn-Arg1.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


ggplot(dot_sum, aes(x = Treatment, y = log2(c(Mon_Ifit1 + Mon_Lyz2 + Mon_Arg1 + Mon_H2-Ab1 + Mon_Ccl5)/c(TAM_C1qa + TAM_Arg1 + TAM_Ifit1 + TAM_Ccl7)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab(paste0("log2 Mons/TAMs ratio")) 
ggsave(file = paste0("fig1_mon-mac.png"), width =  30, height = 60, units = "mm", dpi = 350)
ggsave(file = paste0("fig1_mon-mac.pdf"), width =  30, height = 60, units = "mm", dpi = 350, device = "pdf", bg = 'transparent')


dot_show = c("Neut_S100a8")
ggplot(dot_plot[dot_plot$variable %in% dot_show,], aes(x = Treatment, y = value)) +
  facet_wrap(. ~ variable,  strip.position = 'top', nrow = 1, scales = "free") +
  stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + #c("cornflowerblue", "orange", "red")
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("% of myeloid") 
ggsave(file = paste0("figrev_neut.png"), width = 30, height = 60, units = "mm", dpi = 200)
ggsave(file = paste0("figrev_neut.pdf"), width = 30, height = 60, units = "mm", dpi = 200, device = "pdf", bg = 'transparent')



##############fig 2e#####
dot_panel = c("T")
dot_clusters = clusters_lymphoid
dot_days = c("d07")
dot_exp = c(3, 6)
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_plot = melt(dot_sum)
dot_plot_mean = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = mean)
dot_plot_sd = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = sd)
dot_plot = cbind(dot_plot_mean, sd = dot_plot_sd$x, Name = paste(dot_plot_mean$variable))
dot_plot$Treatment = mapvalues(dot_plot$Treatment, from = c("UT", "m1","m2a"), to = c("Ctrl", "aCTLA4-m1","aCTLA4-m2a"))
temp = data.frame(Name = dot_plot$Name, Treatment = dot_plot$Treatment, x = dot_plot$x)
ref = "Ctrl"
temp = temp[temp$Treatment == ref,]
dot_plot = merge(dot_plot, temp, by = "Name")
dot_plot = dot_plot[-which(dot_plot$Treatment.x == ref) ,]

dot_order_temp = dot_plot[dot_plot$Treatment.x %in% c("aCTLA4-m2a"),]
dot_order_values = log2(dot_order_temp$x.x/dot_order_temp$x.y)
dot_order_lymphoid = as.character(dot_order_temp$variable[order(dot_order_values)])

dot_pvalues_m2a = apply(as.matrix(clusters_lymphoid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m2a",clust])$p.value); names(dot_pvalues_m2a) = clusters_lymphoid
dot_pvalues_m1 = apply(as.matrix(clusters_lymphoid),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m1",clust])$p.value); names(dot_pvalues_m1) = clusters_lymphoid
dot_plot = arrange(dot_plot, by_group = Treatment.x)
dot_plot = data.frame(dot_plot, pvalue = c(dot_pvalues_m1[as.character(unique(dot_plot$variable))], dot_pvalues_m2a[as.character(unique(dot_plot$variable))]))
dot_plot_lymph = dot_plot
ggplot(dot_plot, aes(x = factor(variable, levels = dot_order_lymphoid), y =  log2(x.x/x.y))) +
  facet_grid(.~Treatment.x) +
  # geom_bar(aes(fill = variable), stat = "identity", position = "dodge2") + #, width = x.y + x.x
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = pvalue, size = (x.y + x.x)/2), shape = 21) +
  scale_fill_gradientn(colours = c("red", "white"), trans = "log10") +
  # scale_fill_manual(values =  cluster_color) +
  theme_classic() +
  coord_flip() +
  ylim(-1.5, 1.5)+
  labs(y = paste0("log2(X/Ctrl)")) +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 14),
        panel.background = element_rect(fill = "transparent", size = 0), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent", size =0),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  guides(size = guide_legend(title  = "%"))
ggsave(file = paste0("fig2_enrichment_l.png"), width = 120, height = 80, units = "mm", dpi = 350, device = "png", bg = 'transparent')
ggsave(file = paste0("fig2_enrichment_l.pdf"), width = 120, height = 80, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')



##############fig 2f-i####
dot_panel = c("T")
dot_clusters = clusters_lymphoid
dot_days = c("d07")
dot_exp = c(3, 6)
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% "T" & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment", "Experiment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)
dot_sum_lymph = dot_sum
ggplot(dot_plot[dot_plot$variable %in% "CD4_Cxcr6",], aes(x = Treatment, y = value)) +
  facet_wrap(. ~ variable,  strip.position = 'top', nrow = 1, scales = "free") +
  stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + #c("cornflowerblue", "orange", "red")
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("% of lymphoid") 
ggsave(file = paste0("fig2_makeup_l.png"), width = 30, height = 60, units = "mm", dpi = 150, device = "png", bg = 'transparent')
ggsave(file = paste0("fig2_makeup_l.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

dat = dot_plot[dot_plot$variable %in% "CD4_Cxcr6" & dot_plot$Treatment %in% c("m1", "m2a"),]
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_plot[dot_plot$variable %in% "CD4_Cxcr6" & dot_plot$Treatment %in% c("Ctrl", "m2a"),]
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)



ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/Treg_Foxp3))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, test.args = list(alternative = "two.sided", var.equal = F, paired=F)) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/Treg ratio")
ggsave(file = paste0("figrev2_CD8-Treg.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev2_CD8-Treg.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2((dat$CD8_Ccl4 + dat$CD8_Ccl5 + dat$CD8_Lef1 + dat$CD8_Pdcd1 + dat$CD8_Xcl1 )/(dat$Treg_Foxp3)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2((dat$CD8_Ccl4 + dat$CD8_Ccl5 + dat$CD8_Lef1 + dat$CD8_Pdcd1 + dat$CD8_Xcl1 )/(dat$Treg_Foxp3)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/(CD4_Cxcr6 + CD4_Ifit1)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + 
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/CD4 T ratio") 
ggsave(file = paste0("figrev2_CD8-4.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev2_CD8-4.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2((dat$CD8_Ccl4 + dat$CD8_Ccl5 + dat$CD8_Lef1 + dat$CD8_Pdcd1 + dat$CD8_Xcl1 )/(dat$CD4_Cxcr6 + dat$CD4_Ifit1)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2((dat$CD8_Ccl4 + dat$CD8_Ccl5 + dat$CD8_Lef1 + dat$CD8_Pdcd1 + dat$CD8_Xcl1 )/(dat$CD4_Cxcr6 + dat$CD4_Ifit1)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


##############fig 2i#####
dot_panel = c("T")
dge_assay = "RNA"
dot_clusters = setdiff(clusters_lymphoid, "NK")
dot_days = c("d07")
dot_exp = c(3, 6)
dot_treat = c("m2a", "UT")
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days  & seurat$Treatment %in% dot_treat]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_plot = melt(dot_sum)
dot_plot_mean = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = mean)
dot_plot_sd = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = sd)
dot_plot = cbind(dot_plot_mean, sd = dot_plot_sd$x, Name = paste(dot_plot_mean$variable))
dot_plot$Treatment = mapvalues(dot_plot$Treatment, from = c("UT", "m2a"), to = c("Ctrl", "aCTLA4-m2a"))
temp = data.frame(Name = dot_plot$Name, Treatment = dot_plot$Treatment, x = dot_plot$x)
ref = "Ctrl"
temp = temp[temp$Treatment == ref,]
dot_plot = merge(dot_plot, temp, by = "Name")
dot_plot = dot_plot[-which(dot_plot$Treatment.x == ref) ,]

num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days  & seurat$Treatment %in% "UT"]); length(num_cells)
pct_cluster = dot_plot$x.y; names(pct_cluster) = dot_plot$Name
logfc = log2(dot_plot$x.x/dot_plot$x.y); names(logfc) = dot_plot$Name
dot_pvalues_m2a = apply(as.matrix(names(logfc)),1, function(clust) t.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m2a",clust])$p.value); names(dot_pvalues_m2a) = names(logfc)
ctla4_expr =  as.data.frame(AverageExpression(seurat["Ctla4", num_cells], verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); ctla4_expr$gene = rownames(ctla4_expr)
ggdepletion = data.frame(Cluster = names(logfc), log2FC = logfc[names(logfc)], Ctla4 = as.numeric(ctla4_expr[names(logfc)]), Pct = pct_cluster)

ggplot(ggdepletion, aes(x = log1p(Ctla4), y = log2FC)) +
  geom_point(aes(fill = Cluster, size = Pct), shape = 21) +
  geom_text_repel(aes(label = Cluster)) +
  scale_fill_manual(values = cluster_color[names(logfc)]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "log-normalized Ctla4 expression", y = "Log2 m2a/Ctrl ratio") +
  guides(fill = F)
ggsave(file = paste0("fig2_ctla4.png"), width = 120, height = 100, units = "mm", dpi = 350, device = "png", bg = 'transparent')
ggsave(file = paste0("fig2_ctla4.pdf"), width = 120, height = 100, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')




##############fig 2j#####
dot_clusters = clusters_myeloid
dot_days = c("d07")
dot_exp = c(3, 6)
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% c("Innate", "Broad") & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_plot = melt(dot_sum)
dot_plot_mean = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = mean)
dot_plot_sd = aggregate(dot_plot$value, by = list(Treatment = dot_plot$Treatment, variable = dot_plot$variable), FUN = sd)
dot_plot = cbind(dot_plot_mean, sd = dot_plot_sd$x, Name = paste(dot_plot_mean$variable))
dot_plot$Treatment = mapvalues(dot_plot$Treatment, from = c("UT", "m1","m2a"), to = c("Ctrl", "aCTLA4-m1","aCTLA4-m2a"))
temp = data.frame(Name = dot_plot$Name, Treatment = dot_plot$Treatment, x = dot_plot$x)
ref = "Ctrl"
temp = temp[temp$Treatment == ref,]
dot_plot = merge(dot_plot, temp, by = "Name")
dot_plot = dot_plot[-which(dot_plot$Treatment.x == ref) ,]

dot_order_temp = dot_plot[dot_plot$Treatment.x %in% c("aCTLA4-m2a"),]
dot_order_values = log2(dot_order_temp$x.x/dot_order_temp$x.y)
dot_order_myeloid = as.character(dot_order_temp$variable[order(dot_order_values)])

dot_pvalues_m2a = apply(as.matrix(clusters_myeloid),1, function(clust) wilcox.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m2a",clust])$p.value); names(dot_pvalues_m2a) = clusters_myeloid
dot_pvalues_m1 = apply(as.matrix(clusters_myeloid),1, function(clust) wilcox.test(dot_sum[dot_sum$Treatment == "UT",clust], dot_sum[dot_sum$Treatment == "m1",clust])$p.value); names(dot_pvalues_m1) = clusters_myeloid
dot_plot = arrange(dot_plot, by_group = Treatment.x)
dot_plot = data.frame(dot_plot, pvalue = c(dot_pvalues_m1[as.character(unique(dot_plot$variable))], dot_pvalues_m2a[as.character(unique(dot_plot$variable))]))
dot_plot_mye = dot_plot
ggplot(dot_plot, aes(x = factor(variable, levels = dot_order_myeloid), y =  log2(x.x/x.y))) +
  facet_grid(.~Treatment.x) +
  # geom_bar(aes(fill = variable), stat = "identity", position = "dodge2") + #, width = x.y + x.x
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = pvalue, size = (x.y + x.x)/2), shape = 21) +
  scale_fill_gradientn(colours = c("red", "white"), trans = "log10") +
  # scale_fill_manual(values =  cluster_color) +
  theme_classic() +
  coord_flip() +
  ylim(-2, 2)+
  labs(y = paste0("log2(X/Ctrl)")) +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", size = 0), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent", size =0),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(size = guide_legend(title  = "%"))
ggsave(file = paste0("fig2_enrichment_m.png"), width = 120, height = 80, units = "mm", dpi = 350, device = "png", bg = 'transparent')
ggsave(file = paste0("fig2_enrichment_m.pdf"), width = 120, height = 80, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')


##############fig k-o####
dot_clusters = clusters_myeloid
dot_days = c("d07")
dot_exp = c(3, 6)
mice_out = NA #c("3.07.11", "3.07.15")
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% c("Innate", "Broad") & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment", "Experiment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)
dot_sum_mye = dot_sum

dat = dot_plot[dot_plot$variable %in% "Mreg_Hmox1" & dot_plot$Treatment %in% c("m1", "m2a"),]
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_plot[dot_plot$variable %in% "Mreg_Hmox1" & dot_plot$Treatment %in% c("Ctrl", "m2a"),]
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2(dat$Mon_Ifit1/dat$Mon_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2(dat$Mon_Ifit1/dat$Mon_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2(dat$Mon_Ccl5/dat$Mon_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2(dat$Mon_Ccl5/dat$Mon_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2(dat$TAM_Ifit1/dat$TAM_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2(dat$TAM_Ifit1/dat$TAM_Arg1), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


colnames(dot_sum)[5] = "Mon_MHCII"
ggplot(dot_sum, aes(x = Treatment, y = log2(c(Mon_Ifit1 + Mon_Lyz2 + Mon_Arg1 + Mon_MHCII + Mon_Ccl5)/c(TAM_C1qa + TAM_Arg1 + TAM_Ifit1 + TAM_Ccl7)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a"), c("Ctrl", "m1")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab(paste0("log2 Mons/TAms ratio")) 
ggsave(file = paste0("fig2_mon-mac.png"), width =  30, height = 60, units = "mm", dpi = 350)
ggsave(file = paste0("fig2_mon-mac.pdf"), width =  30, height = 60, units = "mm", dpi = 350, device = "pdf", bg = 'transparent')

dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2(c(dat$Mon_Lyz2 + dat$Mon_Ifit1 + dat$Mon_Ccl5+ dat$`Mon_H2-Ab1` + dat$Mon_Arg1)/(dat$TAM_Arg1 + dat$TAM_C1qa + dat$TAM_Ifit1 + dat$TAM_Ccl7)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment, type = 1)

dat = dot_sum[dot_sum$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2(c(dat$Mon_Lyz2 + dat$Mon_Ifit1 + dat$Mon_Ccl5+ dat$`Mon_H2-Ab1` + dat$Mon_Arg1)/(dat$TAM_Arg1 + dat$TAM_C1qa + dat$TAM_Ifit1 + dat$TAM_Ccl7)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment, type = 1)


dot_clusters = c(clusters_myeloid, "NK")
dot_days = c("d07")
dot_exp = c(3, 6)
mice_out = NA #c("3.07.11", "3.07.15")
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Staining.panel.name %in% c("Innate", "Broad") & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment", "Experiment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)


ggplot(dot_plot[dot_plot$variable %in% "Mreg_Hmox1",], aes(x = Treatment, y = value)) +
  # facet_wrap(. ~ factor(variable, levels = dot_show),  strip.position = 'top', nrow = 1, scales = "free") +
  stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  # geom_text(aes(label = Mouse))+
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  labs(title =  "Mreg_Hmox1", y = "% of myeloid") 
ggsave(file = paste0("figrev2_Mreg_Hmox1.png"), width = 27, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev2_Mreg_Hmox1.pdf"), width = 27, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

colnames(dot_sum)[5] = "Mon_MHCII"
ggplot(dot_sum, aes(x = Treatment, y = log2(c(Mon_Ifit1 + Mon_Lyz2 + Mon_Arg1 + Mon_MHCII + Mon_Ccl5)/c(TAM_C1qa + TAM_Arg1 + TAM_Ifit1 + TAM_Ccl7)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab(paste0("log2 Mons/TAMs ratio")) 
ggsave(file = paste0("figrev2_mon-mac.png"), width =  30, height = 60, units = "mm", dpi = 350)
ggsave(file = paste0("figrev2_mon-mac.pdf"), width =  30, height = 60, units = "mm", dpi = 350, device = "pdf", bg = 'transparent')


ggplot(dot_sum, aes(x = Treatment, y = log2(TAM_Ifit1/TAM_Arg1))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 TAM Ifn/TAM Arg1 ratio")
ggsave(file = paste0("figrev2_TAM_Ifn-Arg1.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev2_TAM_Ifn-Arg1.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

ggplot(dot_sum, aes(x = Treatment, y = log2((TAM_Ifit1 + Mon_Ifit1)/(TAM_Arg1 + Mon_Arg1)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 TAM + Mon Ifn / TAM + Mon Arg1 ratio")
ggsave(file = paste0("figrev2_Ifn-Arg1.png"), width = 30, height = 60, units = "mm", dpi = 150)
ggsave(file = paste0("figrev2_Ifn-Arg1.pdf"), width = 30, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


##############fig 2p-q####
dot_days = c("d07")
dot_exp = c(3, 6)
dot_treat = c("m2a")
clusters_m = c("TAM_Arg1", "TAM_C1qa", "TAM_Ifit1", "TAM_Ccl7",  "Mon_Ifit1","Mon_Arg1", "Mon_Lyz2", "Mon_Ccl5","Mreg_Hmox1",  "Mon_H2-Ab1", "Neut_S100a8") #, "Basophil","cDC2", "migDC"
clusters_l =  c("CD8_Xcl1", "CD8_Pdcd1", "CD8_Ccl4", "CD8_Ccl5", "CD8_Lef1", "CD4_Cxcr6",  "CD4_Ifit1", "Treg_Foxp3") # "B_cells

bad_cells = names(seurat$Biological.replicate[seurat$Biological.replicate %in% "1.07.03"])
num_cells = setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)

consensus_clust = ConsensusClusterPlus(dot_sum1, maxK = 5, reps = 500, pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "spearman", seed = 9, plot = "png")
# icl = calcICL(consensus_clust, plot = "jpg")
# M3C = M3C(dot_sum1, clusteralg= "hc")
# M3C$realdataresults[[4]]
cluster_no = 2
cluster_cons = consensus_clust[[cluster_no]]
cluster_order = cluster_cons$consensusTree$order
cluster_class = cluster_cons$consensusClass[cluster_order]
cluster_mat = cluster_cons$consensusMatrix[cluster_order, cluster_order]
rownames(cluster_mat) = names(cluster_class)
colnames(cluster_mat) = names(cluster_class)
cluster_cor = cor(dot_sum1, method = "spearman")[names(cluster_class),names(cluster_class)]

png(paste0("fig2_consensus_matrix.png"), width = 150, height = 120, units = "mm", res = 250)
pdf(paste0("fig2_consensus_matrix.pdf"), width = 7, height = 6, pointsize = 1)
Heatmap(cluster_mat,
        
        name = "Consensus", 
        col = colorRamp2(breaks = c(0, 1), colors = c("white", "red3")),
        
        cluster_rows = F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        
        cluster_columns = F,
        show_column_names = T,
        column_names_gp = gpar(fontsize = 10, rot = 45),
        column_names_rot = 45,
        
        row_split = cluster_class, 
        column_split = cluster_class,
        
        # right_annotation =  rowAnnotation(cell_type = cell_type[mice_cor_clusts], col =list(cell_type = cluster_color[mice_cor_clusts]), show_legend = F), 
        # bottom_annotation = HeatmapAnnotation(cell_type = cell_type[mice_cor_clusts], col =list(cell_type = cluster_color[mice_cor_clusts]), show_legend = F),
        rect_gp = gpar(col = "grey", lwd = 1),
        border = T,
        use_raster = F, 
        
)
dev.off()


png(paste0("fig2_cor_cell_module.png"), width = 150, height = 120, units = "mm", res = 250)
pdf(paste0("fig2_cor_cell_module.pdf"), width = 7, height = 6, pointsize = 1)
Heatmap(cluster_cor,
        
        name = "Correlation", 
        col = colorRamp2(breaks = c(min(cluster_cor), min(cluster_cor)/2, 0, max(cluster_cor)/2, max(cluster_cor)), colors = c("slateblue4", "skyblue", "white", "tomato", "red3")),
        
        cluster_rows = F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        
        cluster_columns = F,
        show_column_names = T,
        column_names_gp = gpar(fontsize = 10, rot = 45),
        column_names_rot = 45,
        
        row_split = cluster_class, 
        column_split = cluster_class,
        
        # right_annotation =  rowAnnotation(cell_type = cell_type[mice_cor_clusts], col =list(cell_type = cluster_color[mice_cor_clusts]), show_legend = F), 
        # bottom_annotation = HeatmapAnnotation(cell_type = cell_type[mice_cor_clusts], col =list(cell_type = cluster_color[mice_cor_clusts]), show_legend = F),
        rect_gp = gpar(col = "grey", lwd = 1),
        border = T,
        use_raster = F, 
        
)
dev.off()


dot_days = c("d07")
dot_exp = c(3, 6)
dot_treat = c("UT","m1","m2a")
num_cells =   setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
##
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)
##
attributes(dot_sum1)$class <- "matrix"

Mouse_experiment = c("3","6")[apply(as.matrix(table(Mouse, Experiment)[substr(rownames(dot_sum1), 1,9),] > 0),MARGIN = 1, FUN = which)]
mice_mat_meta = data.frame(dot_sum1, Treatment = substr(rownames(dot_sum1), 1,2), Mouse =  substr(rownames(dot_sum1), 1,9), Experiment = Mouse_experiment)
ggmice = melt(mice_mat_meta, id.vars = c("Mouse", "Treatment", "Experiment")); 
ggmice = data.frame(ggmice, Module = factor(cluster_class[as.character(ggmice$variable)]))
ggmice = ddply(ggmice, .(Module, Mouse, Treatment, Experiment), summarize, score = sum(value))
ggmice = ddply(ggmice, .(Module, Treatment, Mouse, Experiment), summarize, mean = mean(score))
ggmice$mean = ggmice$mean/2
ggmice$Treatment = factor(mapvalues(ggmice$Treatment, from = c("1.", "2.", "3."), to = c("Ctrl", "m1", "m2a")), levels = c("Ctrl", "m1", "m2a"))

ggmice_cast = cast(ggmice, Treatment + Mouse  + Experiment ~ Module, mean, value.var = "mean")
colnames(ggmice_cast) = c("Treatment", "Mouse", "Experiment", paste0("Mod", 1:cluster_no))
ggplot(ggmice_cast[ggmice_cast$Treatment %in% c("Ctrl", "m1", "m2a"),], aes(x = Treatment, y = log2((Mod1)/(Mod2)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0, size = 0.3) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20", "navajowhite2", "navajowhite4")) +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Anti-tumor/immune-repressor module")
ggsave(file = "fig2_module.png", width = 30, height = 80, units = "mm", dpi = 250, device = "png", bg = "transparent")
ggsave(file = "fig2_module.pdf", width = 30, height = 80, units = "mm", dpi = 250, device = "pdf", bg = "transparent")

dat = ggmice_cast[ggmice_cast$Treatment %in% c("Ctrl", "m2a"),]
dat = cbind(value = log2(c(dat$Mod1)/(dat$Mod2)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment, type = 2)



##############fig 3f####
alpha = 0.05
dge_assay = "RNA"
dge_logFC_cut = log2(1.25)
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
cluster_coarse = c("Treg", "CD4", "CD8", "Mon", "TAM")  #unique(gsub("_.*","",c(clusters_myeloid, clusters_lymphoid)))

dge_treat =  c("UT","m1", "m2a")
dge_day = c("d07")
dge_panel = c(rep("T", 3), rep("Innate", 2)); names(dge_panel) = cluster_coarse 
dge_exp = c(3)
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_list_m1 = list()
dge_activity_m1 = matrix(0, nrow = length(cluster_coarse), ncol = 2); rownames(dge_activity_m1) = cluster_coarse; colnames(dge_activity_m1) = c("aCTLA4-m1", "aCTLA4-m2a")
dge_list_UT = list()
dge_activity_UT = matrix(0, nrow = length(cluster_coarse), ncol = 2); rownames(dge_activity_UT) = cluster_coarse; colnames(dge_activity_UT) = c("Ctrl", "aCTLA4-m2a")
dge_list_ref = list()
dge_activity_ref = matrix(0, nrow = length(cluster_coarse), ncol = 2); rownames(dge_activity_ref) = cluster_coarse; colnames(dge_activity_ref) = c("Ctrl", "aCTLA4-m1")

for(dge_clusts in cluster_coarse){
  dge_cells_UT = names(seurat@active.ident[gsub("_.*","",seurat@active.ident) %in% dge_clusts & seurat$Treatment %in% c("UT") & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel[dge_clusts]])
  dge_cells_m1 = names(seurat@active.ident[gsub("_.*","",seurat@active.ident) %in% dge_clusts & seurat$Treatment %in% c("m1") & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel[dge_clusts]])
  dge_cells_m2a = names(seurat@active.ident[gsub("_.*","",seurat@active.ident) %in% dge_clusts & seurat$Treatment %in% c("m2a") & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel[dge_clusts]])
  
  dge_sample = min(length(dge_cells_m1), length(dge_cells_m2a))
  dge_cells_m1 = sample(dge_cells_m1, dge_sample, replace = F)
  dge_cells_m2a = sample(dge_cells_m2a, dge_sample, replace = F)
  
  dge_cells = c(dge_cells_UT, dge_cells_m2a); length(dge_cells)
  dge_genes = setdiff( names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
  dge_seurat = seurat[dge_genes, dge_cells]
  Idents(dge_seurat) = "Treatment"
  dge_seurat.wilcox_UT = FindMarkers(dge_seurat, ident.1 = "UT", ident.2 = "m2a", test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, logfc.threshold = 0, min.pct = 0.01, min.cells.group = dge_min_cell, slot = dge_values)
  dge_seurat.wilcox_UT$p_val_adj = p.adjust(dge_seurat.wilcox_UT$p_val, method = "fdr")
  dge_activity_UT[dge_clusts,] = c(length(which(dge_seurat.wilcox_UT$avg_log2FC > dge_logFC_cut & dge_seurat.wilcox_UT$p_val_adj < alpha)), -length(which(dge_seurat.wilcox_UT$avg_log2FC < -dge_logFC_cut & dge_seurat.wilcox_UT$p_val_adj < alpha)))
  dge_list_UT[[dge_clusts]] = dge_seurat.wilcox_UT
  
  dge_cells = c(dge_cells_UT, dge_cells_m1); length(dge_cells)
  dge_genes = setdiff( names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
  dge_seurat = seurat[dge_genes, dge_cells]
  Idents(dge_seurat) = "Treatment"
  dge_seurat.wilcox_ref = FindMarkers(dge_seurat, ident.1 = "UT", ident.2 = "m1", test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, logfc.threshold = 0, min.pct = 0.01, min.cells.group = dge_min_cell, slot = dge_values)
  dge_seurat.wilcox_ref$p_val_adj = p.adjust(dge_seurat.wilcox_ref$p_val, method = "fdr")
  dge_activity_ref[dge_clusts,] = c(length(which(dge_seurat.wilcox_ref$avg_log2FC > dge_logFC_cut & dge_seurat.wilcox_ref$p_val_adj < alpha)), -length(which(dge_seurat.wilcox_ref$avg_log2FC < -dge_logFC_cut & dge_seurat.wilcox_ref$p_val_adj < alpha)))
  dge_list_ref[[dge_clusts]] = dge_seurat.wilcox_ref
}

ggactivity = melt(data.frame(Clust = rownames(dge_activity_ref), -dge_activity_ref, -dge_activity_UT))
ggactivity = data.frame(ggactivity, Treatment = c(rep("Ctrl vs m1", nrow(ggactivity)/2), rep("Ctrl vs m2a", nrow(ggactivity)/2)))
clust_order = names(dge_activity_ref[,2][order(dge_activity_ref[,2], decreasing = F)])
ggplot(ggactivity) + 
  # geom_bar(aes(x = factor(Clust, levels = rev(clust_order)), y = sign(value)*log10(abs(value)), fill = variable), stat = "identity") +
  geom_bar(aes(x = factor(Clust, levels = cluster_coarse), y = value, fill = variable), stat = "identity", width = 0.5) +
  facet_grid(.~ Treatment) +
  scale_fill_manual(values = c("grey90", "grey60", "grey90","grey20")) +
  # scale_x_continuous(trans = "exp") +
  theme_classic() +
  labs(y = "#DGEs") +
  theme(title = element_text(size = 10),
        line = element_line(size = 0.2),
        text = element_text(size = 10),
        strip.background = element_blank(),
        # axis.text.x = element_text(size = 5, angle = 30, vjust = 1, hjust = 1),
        axis.title.y =  element_blank(),
        legend.key.size = unit(0.2, "cm")) +
  coord_flip() 
ggsave(filename = paste0("fig3_reactivity_coarse_",dge_assay,".png"), width = 80, height = 70, units = "mm", dpi = 250)
ggsave(filename = paste0("fig3_reactivity_coarse_",dge_assay,".pdf"), width = 80, height = 70, units = "mm", dpi = 250, device = "pdf", bg = "transparent")



##############fig 3a####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon_Ifit1"
dge_treat =  c("UT","m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c("Mon_Ifit1", "Mon_Arg1") 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
# dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
# dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "Mon_Arg1"
dge_g2_name = "Mon_Ifit1"
# Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 15)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 15)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  scale_color_gradientn(colors = rev(inferno(100))) +
  # scale_fill_viridis() +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 1, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 1, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".png"), width = 150, height = 150, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".pdf"), width = 150, height = 150, units = "mm", dpi = 200, device = "pdf", bg = "transparent")


##############fig 3b####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon_Ccl5"
dge_treat =  c("UT","m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c("Mon_Ccl5", "Mon_Arg1") 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
# dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
# dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "Mon_Arg1"
dge_g2_name = "Mon_Ccl5"
# Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 15)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 15)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  # scale_fill_viridis() +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 1, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 1, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".png"), width = 150, height = 150, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".pdf"), width = 150, height = 150, units = "mm", dpi = 200, device = "pdf", bg = "transparent")


##############fig 3c####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon_Lyz2"
dge_treat =  c("UT","m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c("Mon_Lyz2", "Mon_Arg1") 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
# dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
# dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "Mon_Arg1"
dge_g2_name = "Mon_Lyz2"
# Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 15)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 15)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 1, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 1, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".png"), width = 150, height = 150, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".pdf"), width = 150, height = 150, units = "mm", dpi = 200, device = "pdf", bg = "transparent")



##############fig 3d####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "TAM_Ifit1"
dge_treat =  c("UT","m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c("TAM_Ifit1", "TAM_Arg1") 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
# dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
# dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "TAM_Arg1"
dge_g2_name = "TAM_Ifit1"
# Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 15)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 15)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 1, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 1, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".png"), width = 150, height = 150, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".pdf"), width = 150, height = 150, units = "mm", dpi = 200, device = "pdf", bg = "transparent")

z

##############fig 3e####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "TAM_Ccl7"
dge_treat =  c("UT","m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c("TAM_Ccl7", "TAM_Arg1") 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
# dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
# dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "TAM_Arg1"
dge_g2_name = "TAM_Ccl7"
# Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 15)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 15)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 1, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 1, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".png"), width = 150, height = 150, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,".pdf"), width = 150, height = 150, units = "mm", dpi = 200, device = "pdf", bg = "transparent")



##############fig 3g####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon"
dge_treat =  c("UT","m2a")
dge_day = c("d07")
dge_exp = c(3)
dge_panel = c("Innate")
dge_clusts = c(cluster_monocyte) 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"
dge_pct = 0.1
pseudocount = 0.1
dge_logFC_cut = log2(1.25)

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "UT"
dge_g2_name = "m2a"
Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = pseudocount, verbose = T, min.pct = dge_pct, logfc.threshold = dge_logFC_cut, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 20)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 20)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 2, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 2, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"),
       title = dge_name)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.png"), width = 150, height = 160, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.pdf"), width = 150, height = 160, units = "mm", dpi = 200, device = "pdf", bg = "transparent")


##############fig 3h####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "TAM"
dge_treat =  c("UT","m2a")
dge_day = c("d07")
dge_exp = c(3)
dge_panel = c("Innate")
dge_clusts = c(cluster_macrophage) 
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"
dge_pct = 0.1
pseudocount = 0.1
dge_logFC_cut = log2(1.25)

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "UT"
dge_g2_name = "m2a"
Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = pseudocount, verbose = T, min.pct = dge_pct, logfc.threshold = dge_logFC_cut, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 20)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 20)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 2, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 2, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"),
       title = dge_name)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.png"), width = 150, height = 160, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.pdf"), width = 150, height = 160, units = "mm", dpi = 200, device = "pdf", bg = "transparent")



##############fig 3i####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon-TAM"
dge_treat =  c("UT","m2a")
dge_day = c("d07")
dge_exp = c(3)
dge_panel = c("Innate")
dge_clusts = c(cluster_macrophage, cluster_monocyte) # c("CD8_Ccl5", "CD8_Ccl4", "CD8_Gzma", "CD8_Lef1", "CD8_Xcl1")
dge_min_UMI = 2
dge_min_cell = 5
dge_stat = "wilcox"
dge_values = "data"
dge_pct = 0.1
pseudocount = 0.1
dge_logFC_cut = log2(1.25)


dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "UT"
dge_g2_name = "m2a"
Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = pseudocount, verbose = T, min.pct = dge_pct, logfc.threshold = dge_logFC_cut, slot = dge_values)
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg =  as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay, slot = "data")[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

positive_genes = unique(c(head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 18), "Apoe", "Arg1","Tgm2", "Aldoa", "Eno1"))
negative_genes = unique(c(head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 18), "Ccl8", "Ccl2", "Ccl12", "Cxcl10", "Ccl5"))

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 2, nudge_x = 0, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 6,  nudge_y = 0, nudge_x = 2, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"),
       title = dge_name)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.png"), width = 150, height = 160, units = "mm", dpi = 200)
ggsave(filename = paste0("fig3_DGE_", dge_name, "_", dge_assay,"_rev.pdf"), width = 150, height = 160, units = "mm", dpi = 200, device = "pdf", bg = "transparent")




##############fig 3j####
##prepare ranked gene list 
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon-Mac"
dge_treat =  c("UT","m2a")
dge_day = c("d07")
dge_exp = c(3)
dge_panel = c("Innate")
dge_clusts = c(cluster_macrophage, cluster_monocyte)
dge_min_UMI = 2
dge_min_cell = 5

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[dge_genes, dge_cells]
dge_g1_name = "m2a"
dge_g2_name = "UT"
Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, latent.vars = "Biological.replicate", slot = "data")
logFC = dge_seurat.wilcox$avg_log2FC; names(logFC) = rownames(dge_seurat.wilcox)
gsea_list = sort(logFC, decreasing = T)
egSYMBOL = AnnotationDbi::toTable(org.Mm.egSYMBOL)
egSYMBOL = egSYMBOL[-which(duplicated(egSYMBOL$symbol)),] 
rownames(egSYMBOL) = egSYMBOL$symbol
names(gsea_list) = egSYMBOL[names(gsea_list),1]
gsea_list = gsea_list[-which(is.na(names(gsea_list)))]


##fgsea
gsea_pathways = get(load("/home/labs/amit/tomerlan/analysis/data_files/mouse_c5_v5p2.RData"))

fgsea = fgsea(gsea_pathways, gsea_list, minSize = 10, maxSize = 500)

gene_key = egSYMBOL
rownames(gene_key) = gene_key$gene_id

for(i in 1:nrow(fgsea)){
  fgsea[i, 8] = gene_key[fgsea$leadingEdge[[i]], "symbol"]
}


gse_select = read.csv("fig3_fgsea_select.csv")
gse_order = gse_select$pathway

ggplot(gse_select) +
  geom_point(aes(y = factor(substr(tolower(pathway), 4, 100), levels = rev(substr(tolower(gse_order), 4, 100))), color = NES, size = -log10(padj)), x = 0) +
  scale_color_gradient2(low = "royalblue4", mid = "white", high = "red") +
  scale_y_discrete(position = "right") +
  theme(axis.title = element_blank(),
        legend.position = "bottom") 
ggsave(file = paste0("fig3_GSEA.png"), width = 120, height = 130, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_GSEA.pdf"), width = 250, height = 130, units = "mm", dpi = 250, device = "pdf", bg = "transparent")



##############fig 3k-i####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "Mon-Mac"
dge_treat =  c("UT", "m2a")
dge_day = c("d07")
dge_exp = c(3)
dge_panel = c("Innate")
dge_clusts = c(cluster_macrophage, cluster_monocyte)
dge_logFC_cut = log2(1.25)
dge_stat = "wilcox"
dge_values = "data"
dge_pct = 0.1
pseudocount = 0.1

dge_lig_genes = intersect(unique(LR_pairs$ligand_gene_symbol), rownames(seurat))
dge_rec_genes = intersect(unique(LR_pairs$receptor_gene_symbol), rownames(seurat))
dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells])); dge_min_mouse
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells)
dge_g1_name = "UT"
dge_g2_name = "m2a"

lr_clusters = c(clusters_lymphoid, clusters_myeloid)
lr_treat = c("m2a", "UT")
lr_days = c("d07")
lr_experiments = c(3)
lr_assay = "RNA"
lr_mcells = names(seurat$Treatment[seurat$Day.of.tumor.harvest %in% lr_days & seurat$Experiment %in% lr_experiments & seurat$Treatment %in% lr_treat & seurat@active.ident %in% intersect(lr_clusters, clusters_myeloid) & seurat$Staining.panel.name %in% "Innate"]); length(lr_mcells)
# lr_min_mouse = min(table(seurat$Biological.replicate[lr_mcells])); lr_min_mouse
# lr_mcells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[lr_mcells])), 1, function(x) sample(x = names(seurat$Biological.replicate[lr_mcells])[seurat$Biological.replicate[lr_mcells] %in% x], size = lr_min_mouse, replace = F))); length(lr_mcells)
lr_lcells = names(seurat$Treatment[seurat$Day.of.tumor.harvest %in% lr_days & seurat$Experiment %in% lr_experiments & seurat$Treatment %in% lr_treat & seurat@active.ident %in% intersect(lr_clusters, clusters_lymphoid) & seurat$Staining.panel.name %in% "T"]); length(lr_lcells)
# lr_min_mouse = min(table(seurat$Biological.replicate[lr_lcells])); lr_min_mouse
# lr_lcells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[lr_lcells])), 1, function(x) sample(x = names(seurat$Biological.replicate[lr_lcells])[seurat$Biological.replicate[lr_lcells] %in% x], size = lr_min_mouse, replace = F))); length(lr_lcells)
lr_cells = c(lr_mcells, lr_lcells)

##ligand DGE
dge_seurat = seurat[dge_lig_genes, dge_cells]
Idents(dge_seurat) = "Treatment"
dge_lig = FindMarkers(dge_seurat, features = dge_lig_genes, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = dge_stat, assay = dge_assay, pseudocount.use = pseudocount, verbose = T, min.pct = dge_pct, logfc.threshold = dge_logFC_cut, slot = dge_values)
dge_lig$p_val_adj = p.adjust(dge_lig$p_val, method = "fdr")

##ligand-receptor UT
dge_lig_down = rownames(dge_lig[dge_lig$p_val_adj < alpha & dge_lig$avg_log2FC > dge_logFC_cut,]); dge_lig_down
dge_rec_down = unique(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_down, "receptor_gene_symbol"]); dge_rec_down
sort(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_down, "lr_pair"])
# lr_table = cast(lr_table, formula =  receptor_gene_symbol ~  ligand_gene_symbol)

rec_cells = intersect(names(seurat$Treatment[seurat$Treatment == "UT"]), lr_cells)
rec_genes = dge_rec_down[rowSums(seurat@assays$RNA@counts[dge_rec_down, rec_cells]) > 30]; rec_genes
lr_table_UT = as.data.frame(t(table(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_down & LR_pairs$receptor_gene_symbol %in% rec_genes, c("ligand_gene_symbol", "receptor_gene_symbol")])))
rec_order = hclust(dist(cast(lr_table_UT, formula =  receptor_gene_symbol ~  ligand_gene_symbol)))$labels

ggplot(lr_table_UT) +
  geom_point(aes(x = ligand_gene_symbol, y = factor(receptor_gene_symbol, levels = rec_order), size = Freq), show.legend = F, shape = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank()) 
ggsave(file = paste0("fig3_lig_down.png"), width = 70, height = 120, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_lig_down.pdf"), width = 70, height = 120, units = "mm", dpi = 250, device = "pdf")

DotPlot(seurat[, rec_cells], features = rec_order, assay = dge_assay, cols = c("white", "red"), cluster.idents = T) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))
ggsave(file = paste0("fig3_rec_down.png"), width = 150, height = 100, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_rec_down.pdf"), width = 150, height = 100, units = "mm", dpi = 250, device = "pdf", bg = "transparent")


##ligand-receptor m2a
dge_lig_up = rownames(dge_lig[dge_lig$p_val_adj < alpha & dge_lig$avg_log2FC < -dge_logFC_cut,]); dge_lig_up
dge_rec_up = intersect(unique(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_up, "receptor_gene_symbol"]), rownames(seurat@assays$RNA@counts)); dge_rec_up
sort(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_up, "lr_pair"])

rec_cells = intersect(names(seurat$Treatment[seurat$Treatment == "m2a"]), lr_cells)
rec_genes = dge_rec_up[rowSums(seurat@assays$RNA@counts[dge_rec_up, rec_cells]) > 20]; rec_genes
lr_table_m2a = as.data.frame(t(table(LR_pairs[LR_pairs$ligand_gene_symbol %in% dge_lig_up & LR_pairs$receptor_gene_symbol %in% rec_genes, c("ligand_gene_symbol", "receptor_gene_symbol")])))
rec_order = hclust(dist(cast(lr_table_m2a, formula =  receptor_gene_symbol ~  ligand_gene_symbol)))$labels

ggplot(lr_table_m2a) +
  geom_point(aes(x = ligand_gene_symbol, y = factor(receptor_gene_symbol, levels = rec_order), size = Freq), show.legend = F, shape = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank()) 
ggsave(file = paste0("fig3_lig_up.png"), width = 60, height = 90, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_lig_up.pdf"), width = 60, height = 90, units = "mm", dpi = 250, device = "pdf")

DotPlot(seurat[, rec_cells], features = rec_order, assay = dge_assay, cols = c("white", "red"), cluster.idents = T) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))
ggsave(file = paste0("fig3_rec_up.png"), width = 150, height = 90, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_rec_up.pdf"), width = 150, height = 90, units = "mm", dpi = 250, device = "pdf", bg = "transparent")


##############fig 3m#####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_treat =  c("UT", "m2a", "m1")
dge_day = c("d07")
dge_exp = c(3, 6)
dge_panel = c("Innate")
dge_clusts = c(cluster_macrophage, cluster_monocyte)
dge_logFC_cut = log2(1.25)
dge_values = "data"
dge_pct = 0.1

dge_lig_genes = c("Fcgr1", "Fcgr4","Fcgr2b", "Fcgr3")
dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)

DotPlot(seurat[, dge_cells], features = dge_lig_genes, assay = dge_assay, cols = c("white", "red"), cluster.idents = T) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"))
ggsave(file = paste0("fig3_fcgr_",dge_assay,".png"), width = 80, height = 55, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig3_fcgr_",dge_assay,".pdf"), width = 80, height = 55, units = "mm", dpi = 250, device = "pdf", bg = "transparent")


##############fig 4a####
stats = as.data.frame(read_xlsx("/home/labs/amit/tomerlan/UCL_project/experiments_summary.xlsx", sheet = 4, col_names = T))

ggplot(stats[stats$Treatment != "m2a" & stats$Tissue != "LN",], aes(x = factor(paste(Treatment, Dosage), levels = c("PBS 0", "m2a 0", "DT 1", "DT 2", "DT 5", "DT 10")), y = Freq)) +
  geom_boxplot(show.legend = F, width = 0.6) +
  geom_jitter(show.legend = F, width = 0.2, size = 1) +
  facet_grid(factor(Tissue, levels = c("Tumor", "LN")) ~., scales = "free") +
  ylab("% Treg of CD4+") +
  xlab("DT dosage [ug/gr]") +
  theme_classic() +
  expand_limits(y = 0)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size= 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))
ggsave(file = "fig6_DT_titration.png", width = 100, height = 60, units = "mm", dpi = 150, bg = "transparent")
ggsave(file = "fig6_DT_titration.pdf", width = 100, height = 60, units = "mm", dpi = 150, bg = "transparent")



##############fig 4b####
stats = as.data.frame(read_xlsx("/home/labs/amit/tomerlan/UCL_project/tumor.xlsx", sheet = 11, col_names = T))
class(stats$day) = "numeric"
class(stats$vol) = "numeric"
stats$rescued = as.character(stats$rescued)

ggplot(stats, mapping = aes(x = day, y = vol, group = paste(mouse, treatment))) +
  facet_grid(.~ factor(treatment, levels = c("UT", "m2a", "DT"))) +
  geom_line() +
  # geom_path(arrow = arrow()) +
  # geom_point(aes(color = treatment), size = 2) +
  scale_color_manual(values = c("firebrick1", "dodgerblue")) +
  xlab("Days after tumor injection") +
  ylab("Tumor volume [mm^3]") +
  scale_x_continuous(breaks = unique(stats$day)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size= 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))
ggsave(file = "fig6_tumor_growth_DTR.png", width = 420, height = 120, units = "mm", dpi = 150, bg = "transparent")
ggsave(file = "fig6_tumor_growth_DTR.pdf", width = 420, height = 120, units = "mm", dpi = 150, bg = "transparent")



##############fig 4d####
dot_days = c("d07")
dot_exp = c(15, 17)
dot_treat = c("DTR_DT", "DTR_UT")
dot_panel = c("T")
num_cells =  names(seurat$orig.ident[seurat$Treatment %in% dot_treat & as.character(seurat@active.ident) %in% clusters_lymphoid & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("DTR_UT", "DTR_DT")), from = c("DTR_UT", "DTR_DT"), to = c("Ctrl", "DT"))
dot_plot = melt(dot_sum)
dot_sum_lymph = dot_sum

ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/c(Treg_Foxp3)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("navajowhite2", "navajowhite4")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("Ctrl", "DT")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, test.args = list(alternative = "less", var.equal = F, paired=F)) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/Treg ratio")
ggsave(file = paste0("fig4_CD8-Treg.png"), width = 30, height = 80, units = "mm", dpi = 150)
ggsave(file = paste0("fig4_CD8-Treg.pdf"), width = 30, height = 80, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')



##############fig 4e####
dot_days = c("d07")
dot_exp = c(15, 17)
dot_treat = c("DTR_UT", "DTR_DT")
dot_panel = c("Innate","Broad")
num_cells =  setdiff(names(seurat$orig.ident[seurat$Treatment %in% dot_treat & as.character(seurat@active.ident) %in% clusters_myeloid & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]), names(seurat$amp_batch_id[seurat$amp_batch_id == "AB3636"])); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells]); #ct[ct == "Mac_Pf4"] = "Mac_Arg1"
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_sum$Treatment = factor(mapvalues(dot_sum$Treatment, from = c("DTR_UT", "DTR_DT"), to = c("Ctrl", "DT")), levels = c("Ctrl", "DT"))
dot_sum_mye = dot_sum
colnames(dot_sum)[5] = "Mon_MHCII"
ggplot(dot_sum, aes(x = Treatment, y = log2(c(Mon_Ifit1 + Mon_Lyz2 + Mon_Arg1 + Mon_MHCII + Mon_Ccl5)/c(TAM_C1qa + TAM_Arg1 + TAM_Ifit1 + TAM_Ccl7)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0, size = 0.3) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("navajowhite2", "navajowhite4")) +
  geom_signif(comparisons = list(c("DT", "Ctrl")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  ylab(paste0("log2 Monocytes/Macrophages ratio)"))
ggsave(file = paste0("fig4_DTR_mon-mac.png"), width = 30, height = 80, units = "mm", dpi = 350)
ggsave(file = paste0("fig4_DTR_mon-mac.pdf"), width = 30, height = 80, units = "mm", dpi = 350)





##############fig 4f###########
dot_days = c("d07")
dot_exp = c(3, 6)
dot_treat = c("m2a")
clusters_m = c(cluster_macrophage, cluster_monocyte, "Mreg_Hmox1", "Neut_S100a8") #, "Basophil","cDC2", "migDC"
clusters_l =  c("CD8_Xcl1", "CD8_Pdcd1", "CD8_Ccl4", "CD8_Ccl5", "CD8_Lef1", "CD4_Cxcr6",  "CD4_Ifit1", "Treg_Foxp3") # "B_cells
bad_cells = c(names(seurat$Biological.replicate[seurat$Biological.replicate %in% "1.07.03"]))
num_cells =   setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)

consensus_clust = ConsensusClusterPlus(dot_sum1, maxK = 5, reps = 50, pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "spearman", seed = 19, plot = "png")

clust_no = 2
mice_clust_tree = consensus_clust[[clust_no]]$consensusClass;mice_clust_tree
mice_cor_clusts = names(sort(mice_clust_tree))


dot_days = c("d07")
dot_exp = c(15, 17)
dot_treat = c("DTR_UT", "DTR_DT")
num_cells =   setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
##
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)
##
attributes(dot_sum1)$class <- "matrix"

mice_mat_meta = data.frame(dot_sum1, Treatment = substr(rownames(dot_sum1), 1,2), Mouse =  substr(rownames(dot_sum1), 1,9))
ggmice = melt(mice_mat_meta, id.vars = c("Mouse", "Treatment")); 
ggmice = data.frame(ggmice, Module = factor(mice_clust_tree[ggmice$variable]))
ggmice = ddply(ggmice, .(Module, Mouse, Treatment), summarize, score = sum(value))
ggmice = ddply(ggmice, .(Module, Treatment, Mouse), summarize, mean = mean(score))
ggmice$mean = ggmice$mean/2
ggmice$Treatment = factor(mapvalues(ggmice$Treatment, from = c("11","13"), to = c( "Ctrl", "DT")), levels = c("Ctrl", "DT"))

ggmice_cast = cast(ggmice, Treatment + Mouse  ~ Module, mean, value.var = "mean")
colnames(ggmice_cast) = c("Treatment", "Mouse", paste0("Mod", 1:clust_no))
mice_clust_tree
ggplot(ggmice_cast, aes(x = Treatment, y = log2((Mod1)/(Mod2)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0, size = 0.3) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c( "navajowhite2", "navajowhite4")) +
  geom_signif(comparisons = list(c("Ctrl", "DT")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Anti-tumor/immune-repressor module")
ggsave(file = "fig4_DTR_atmodule.png", width = 30, height = 80, units = "mm", dpi = 250, device = "png", bg = "transparent")
ggsave(file = "fig4_DTR_atmodule.pdf", width = 30, height = 80, units = "mm", dpi = 250, device = "pdf", bg = "transparent")



##############fig 4g######
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, stress_genes, cc_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dge_name = "Mon-mac"
dge_min_UMI = 2
dge_min_cell = 5
dge_clusters = c(cluster_monocyte, cluster_macrophage)
dge_panel =  c("Innate", "Broad")

unique(seurat@active.ident)
g1_clust =  dge_clusters
g1_treat =  "m2a"
g1_day = c("d07")
g1_exp = c(3)
g1_panel = dge_panel
g1_name = "WT_m2a"

g2_clust =  dge_clusters
g2_treat =  "UT"
g2_day = c("d07")
g2_exp = c(3)
g2_panel = dge_panel
g2_name = "WT_UT"

unique(seurat@active.ident)
g3_clust = dge_clusters
g3_treat = "DTR_DT"
g3_day = c("d07")
g3_exp = c(15,17)
g3_panel = dge_panel
g3_name = "DTR_DT"

g4_clust = dge_clusters
g4_treat =  "DTR_UT"
g4_day = c("d07")
g4_exp = c(15,17)
g4_panel = dge_panel
g4_name = "DTR_UT"

g1 = names(seurat@active.ident[seurat@active.ident %in% g1_clust & seurat$Treatment %in% g1_treat & seurat$Day.of.tumor.harvest %in% g1_day & seurat$Experiment %in% g1_exp & seurat$Staining.panel.name %in% g1_panel]); length(g1)
g2 = names(seurat@active.ident[seurat@active.ident %in% g2_clust & seurat$Treatment %in% g2_treat & seurat$Day.of.tumor.harvest %in% g2_day & seurat$Experiment %in% g2_exp & seurat$Staining.panel.name %in% g2_panel]); length(g2)
g3 = names(seurat@active.ident[seurat@active.ident %in% g3_clust & seurat$Treatment %in% g3_treat & seurat$Day.of.tumor.harvest %in% g3_day & seurat$Experiment %in% g3_exp & seurat$Staining.panel.name %in% g3_panel]); length(g3)
g4 = names(seurat@active.ident[seurat@active.ident %in% g4_clust & seurat$Treatment %in% g4_treat & seurat$Day.of.tumor.harvest %in% g4_day & seurat$Experiment %in% g4_exp & seurat$Staining.panel.name %in% g4_panel]); length(g4)

dge_cells12 = c(g1, g2)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells12])); dge_min_mouse
dge_cells12 = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells12])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells12])[seurat$Biological.replicate[dge_cells12] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells12)
high_genes12 =  names(which(apply(seurat@assays$RNA@counts[, dge_cells12], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell)))
dge_genes12 = setdiff(high_genes12, dge_bad_features); length(dge_genes12)
dge_seurat12 = seurat[dge_genes12, dge_cells12]
# dge_seurat12 = NormalizeData(dge_seurat12, assay = dge_assay, normalization.method = "LogNormalize", scale.factor = 10000)
dge_g1_name = "UT"
dge_g2_name = "m2a"
Idents(dge_seurat12) = "Treatment"
dge_seurat.avg12 =  as.data.frame(AverageExpression(dge_seurat12, verbose = T, assays = dge_assay)[[dge_assay]]); dge_seurat.avg12$gene = rownames(dge_seurat.avg12)
dge_seurat.avg12$density =  get_density(dge_seurat.avg12[[dge_g1_name]], dge_seurat.avg12[[dge_g2_name]], n = 1000)
dge_seurat.wilcox12 = FindMarkers(dge_seurat12, ident.1 = dge_g2_name, ident.2 = dge_g1_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = "data")
dge_seurat.wilcox12$p_val_adj = p.adjust(dge_seurat.wilcox12$p_val, method = "fdr")

dge_cells34 = c(g3, g4)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells34])); dge_min_mouse
dge_cells34 = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells34])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells34])[seurat$Biological.replicate[dge_cells34] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells34)
high_genes34 =  names(which(apply(seurat@assays$RNA@counts[, dge_cells34], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell)))
dge_genes34 = setdiff(high_genes34, dge_bad_features); length(dge_genes34)
dge_seurat34 = seurat[dge_genes34, dge_cells34]
# dge_seurat34 = NormalizeData(dge_seurat34, assay = dge_assay, normalization.method = "LogNormalize", scale.factor = 10000)
dge_g3_name = "DTR_UT"
dge_g4_name = "DTR_DT"
Idents(dge_seurat34) = "Treatment"
dge_seurat.avg34 =  as.data.frame(AverageExpression(dge_seurat34, verbose = T, assays = dge_assay)[[dge_assay]]); dge_seurat.avg34$gene = rownames(dge_seurat.avg34)
dge_seurat.avg34$density =  get_density(dge_seurat.avg34[[dge_g3_name]], dge_seurat.avg34[[dge_g4_name]], n = 1000)
dge_seurat.wilcox34 = FindMarkers(dge_seurat34, ident.1 = dge_g4_name, ident.2 = dge_g3_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = "data")
dge_seurat.wilcox34$p_val_adj = p.adjust(dge_seurat.wilcox34$p_val, method = "fdr")

dge_logFC_cut = log2(1.5)
dge_seurat.wilcox12.sig.pos = dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox12.sig.pos)
dge_seurat.wilcox12.sig.neg = dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox12.sig.neg)
dge_seurat.wilcox34.sig.pos = dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox34.sig.pos)
dge_seurat.wilcox34.sig.neg = dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox34.sig.neg)

ngenes = 15
positive_genes12 = head(rownames(dge_seurat.wilcox12.sig.pos)[order(dge_seurat.wilcox12.sig.pos$avg_log2FC, decreasing = T)], ngenes)
negative_genes12 = head(rownames(dge_seurat.wilcox12.sig.neg)[order(dge_seurat.wilcox12.sig.neg$avg_log2FC, decreasing = F)], ngenes)
positive_genes34 = head(rownames(dge_seurat.wilcox34.sig.pos)[order(dge_seurat.wilcox34.sig.pos$avg_log2FC, decreasing = T)], ngenes)
negative_genes34 = head(rownames(dge_seurat.wilcox34.sig.neg)[order(dge_seurat.wilcox34.sig.neg$avg_log2FC, decreasing = F)], ngenes)
shared_genes_pos = intersect(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]))
shared_genes_neg = intersect(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]))
shared_genes_posneg = intersect(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]))
shared_genes_negpos = intersect(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]))
shared_genes = intersect(c(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,])), c(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,])))

fc_genes = intersect(rownames(dge_seurat.wilcox12), rownames(dge_seurat.wilcox34)); length(fc_genes)

ggde_table = data.frame(gene = fc_genes, wt = dge_seurat.wilcox12[fc_genes,]$avg_log2FC, ko = dge_seurat.wilcox34[fc_genes,]$avg_log2FC)
ggde_table$density =  get_density(ggde_table$wt, ggde_table$ko, n = 100)
nudge = 0.4
force = 30
size = 6

ggplot(ggde_table, aes(x = wt, y = ko)) +
  geom_hline(yintercept = c(-dge_logFC_cut, dge_logFC_cut), linetype = "dashed") +
  geom_vline(xintercept = c(-dge_logFC_cut, dge_logFC_cut), linetype = "dashed") +
  geom_point(aes(color = density), size = 2, shape = 19, stroke = 0) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes12, shared_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes12, shared_genes),], color = "red", aes(label = gene), segment.size  = 0.15, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes12, shared_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes12, shared_genes),], color = "red", aes(label = gene), segment.size  = 0.15, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes34, shared_genes),], color = "blue", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes34, shared_genes),], color = "blue", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes34, shared_genes),], color = "blue", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes34, shared_genes),], color = "blue", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_pos,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_pos,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_neg,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_neg,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_posneg,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_posneg,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_negpos,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_negpos,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  theme(text = element_text(size = 25),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0)) +
  expand_limits(x = c(min(ggde_table$wt, ggde_table$ko), max(ggde_table$wt, ggde_table$ko)), y =  c(min(ggde_table$wt, ggde_table$ko), max(ggde_table$wt, ggde_table$ko))) +
  labs(x = "log2FC m2a/Ctrl WT", y = "log2FC DT/Ctrl DTR", title = dge_name) 
ggsave(filename = paste0("fig4_DTR_fcfc_mon-mac_",dge_assay,"_combined.png"), width = 260, height = 230, units = "mm", dpi = 200)
ggsave(filename = paste0("fig4_DTR_fcfc_mon-mac_",dge_assay,"_combined.pdf"), width = 260, height = 230, units = "mm", dpi = 200, bg = "transparent")


##############fig 4i####
unique(seurat@active.ident)
dot_days = c("d07")
dot_exp = c(11, 12, 21)
dot_clusters = clusters_lymphoid
dot_panel = c("T")
dot_treat = c("FCGR_KO+m1", "FCGR_KO+m2a")
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Treatment %in% dot_treat & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("FCGR_KO+m1", "FCGR_KO+m2a")), from = c("FCGR_KO+m1", "FCGR_KO+m2a"), to = c("m1", "m2a"))
dot_plot = melt(dot_sum)
dot_sum_lymph2 = dot_sum

ggplot(dot_sum, aes(x = Treatment, y = log2(c(CD8_Ccl5 + CD8_Ccl4 + CD8_Pdcd1 + CD8_Lef1 + CD8_Xcl1)/c(Treg_Foxp3)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("lightblue1", "lightblue4")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, test.args = list(alternative = "less", var.equal = F, paired=F)) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 CD8/Treg ratio")
ggsave(file = paste0("fig4_fcgr_CD8-Treg.png"), width = 30, height = 80, units = "mm", dpi = 150)
ggsave(file = paste0("fig4_fcgr_CD8-Treg.pdf"), width = 30, height = 80, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')


##############fig 4j####
unique(seurat@active.ident)
dot_days = c("d07")
dot_exp = c(11, 12, 21)
dot_clusters = clusters_myeloid
dot_panel = c("Innate", "Broad")
dot_treat = c("FCGR_KO+m1", "FCGR_KO+m2a")
num_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dot_clusters & seurat$Treatment %in% dot_treat & seurat$Staining.panel.name %in% dot_panel & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells]); #ct[ct == "Mac_Pf4"] = "Mac_Arg1"
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse); dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_sum$Treatment = factor(mapvalues(dot_sum$Treatment, from = c("FCGR_KO+m1", "FCGR_KO+m2a"), to = c("m1", "m2a")), levels = c("m1", "m2a"))
dot_sum_mye2 = dot_sum
colnames(dot_sum)[5] = "Mon_MHCII"
ggplot(dot_sum, aes(x = Treatment, y = log2(c(Mon_Ifit1 + Mon_Lyz2 + Mon_Arg1 + Mon_MHCII + Mon_Ccl5)/c(TAM_C1qa + TAM_Arg1 + TAM_Ifit1 + TAM_Ccl7)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0, size = 0.3) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c("lightblue1", "lightblue4")) +
  geom_signif(comparisons = list(c("m1", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15, size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  ylab("log2 Monocytes/Macrophages ratio")
ggsave(file = paste0("fig4_FCGRKO_mon-mac.png"), width = 30, height = 80, units = "mm", dpi = 350)
ggsave(file = paste0("fig4_FCGRKO_mon-mac.pdf"), width = 30, height = 80, units = "mm", dpi = 350)



##############fig 4k###########
dot_days = c("d07")
dot_exp = c(3, 6)
dot_treat = c("m2a")
clusters_m = c(cluster_macrophage, cluster_monocyte, "Mreg_Hmox1", "Neut_S100a8") #, "Basophil","cDC2", "migDC"
clusters_l =  c("CD8_Xcl1", "CD8_Pdcd1", "CD8_Ccl4", "CD8_Ccl5", "CD8_Lef1", "CD4_Cxcr6",  "CD4_Ifit1", "Treg_Foxp3") # "B_cells
bad_cells = names(seurat$Biological.replicate[seurat$Biological.replicate %in% "1.07.03"])
num_cells =   setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)

consensus_clust = ConsensusClusterPlus(dot_sum1, maxK = 5, reps = 50, pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "spearman", seed = 9, plot = "png")

clust_no = 2
mice_clust_tree = consensus_clust[[clust_no]]$consensusClass
mice_cor_clusts = names(sort(mice_clust_tree))

dot_days = c("d07")
dot_exp = c(11, 12, 21)
dot_treat = c("FCGR_KO+m1", "FCGR_KO+m2a")
num_cells =   setdiff(c(names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_l &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% "T"]), names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_m &  seurat$Treatment %in% dot_treat & as.vector(seurat$Experiment) %in% dot_exp & seurat$Day.of.tumor.harvest %in% dot_days & seurat$Staining.panel.name %in% c("Innate", "Broad")])), bad_cells); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, ct = ct)
dot_sum1 = table(dot_data)/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
##
dot_sum1_lymph = table(dot_data)[,clusters_l]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_l])*50
dot_sum1_myel = table(dot_data)[,clusters_m]/rowSums(table(dot_data[,c("Mouse", "ct")])[,clusters_m])*50
dot_sum1 = cbind(dot_sum1_lymph, dot_sum1_myel)
##
attributes(dot_sum1)$class <- "matrix"

Mouse_experiment = dot_exp[apply(as.matrix(table(Mouse, Experiment)[substr(rownames(dot_sum1), 1,9),]>0),MARGIN = 1, FUN = which)]
mice_mat_meta = data.frame(dot_sum1, Treatment = substr(rownames(dot_sum1), 1,2), Mouse =  substr(rownames(dot_sum1), 1,9), Experiment = Mouse_experiment)
ggmice = melt(mice_mat_meta, id.vars = c("Mouse", "Treatment", "Experiment")); 
ggmice = data.frame(ggmice, Module = factor(mice_clust_tree[ggmice$variable]))
ggmice = ddply(ggmice, .(Module, Mouse, Treatment, Experiment), summarize, score = sum(value))
ggmice = ddply(ggmice, .(Module, Treatment, Mouse, Experiment), summarize, mean = mean(score))
ggmice$mean = ggmice$mean/2
ggmice$Treatment = factor(mapvalues(ggmice$Treatment, from = c("7.","8."), to = c ("m1", "m2a")), levels = c("m1", "m2a"))

ggmice_cast = cast(ggmice, Treatment + Mouse  + Experiment ~ Module, mean, value.var = "mean")
colnames(ggmice_cast) = c("Treatment", "Mouse", "Experiment", paste0("Mod", 1:clust_no))
mice_clust_tree
ggplot(ggmice_cast, aes(x = Treatment, y = log2((Mod1)/(Mod2)))) +
  # stat_summary(color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.4, fill = NA) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0, size = 0.3) +
  # geom_text_repel(aes(label = Mouse), size = 2) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  scale_fill_manual(values = c( "lightblue1", "lightblue4")) +
  geom_signif(comparisons = list(c("m1", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 Anti-tumor/immune-repressor module")
ggsave(file = "fig4_FCGRKO_atmodule.png", width = 30, height = 80, units = "mm", dpi = 250, device = "png", bg = "transparent")
ggsave(file = "fig4_FCGRKO_atmodule.pdf", width = 30, height = 80, units = "mm", dpi = 250, device = "pdf", bg = "transparent")

dat = ggmice_cast[ggmice_cast$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2(c(dat$Mod1)/(dat$Mod2)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment, type = 1)




##############fig 4i######
alpha = 0.05
dge_min_UMI = 2
dge_min_cells = 5
dge_assay = "RNA"
bad_features_DGE = c(ig_genes, hb_genes, stress_genes, cc_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dge_name= "Mon-mac"

unique(seurat@active.ident)
g1_clust =  c(cluster_monocyte, cluster_macrophage)
g1_treat =  "m2a"
g1_day = c("d07")
g1_exp = c(3)
g1_panel = c("Innate")
g1_name = "WT_m2a"

g2_clust =  c(cluster_monocyte, cluster_macrophage)
g2_treat =  "m1"
g2_day = c("d07")
g2_exp = c(3)
g2_panel = c("Innate")
g2_name = "WT_m1"

unique(seurat@active.ident)
g3_clust =  c(cluster_monocyte, cluster_macrophage)
g3_treat = "FCGR_KO+m2a"
g3_day = c("d07")
g3_exp = c(11, 12, 21)
g3_panel = c("Broad")
g3_name = "FCGRKO_m2a"

g4_clust =  c(cluster_monocyte, cluster_macrophage)
g4_treat =  "FCGR_KO+m1"
g4_day = c("d07")
g4_exp = c(11, 12, 21)
g4_panel = c("Broad")
g4_name = "FCGRKO_m1"

g1 = names(seurat@active.ident[seurat@active.ident %in% g1_clust & seurat$Treatment %in% g1_treat & seurat$Day.of.tumor.harvest %in% g1_day & seurat$Experiment %in% g1_exp & seurat$Staining.panel.name %in% g1_panel]); length(g1)
g2 = names(seurat@active.ident[seurat@active.ident %in% g2_clust & seurat$Treatment %in% g2_treat & seurat$Day.of.tumor.harvest %in% g2_day & seurat$Experiment %in% g2_exp & seurat$Staining.panel.name %in% g2_panel]); length(g2)
g3 = names(seurat@active.ident[seurat@active.ident %in% g3_clust & seurat$Treatment %in% g3_treat & seurat$Day.of.tumor.harvest %in% g3_day & seurat$Experiment %in% g3_exp & seurat$Staining.panel.name %in% g3_panel]); length(g3)
g4 = names(seurat@active.ident[seurat@active.ident %in% g4_clust & seurat$Treatment %in% g4_treat & seurat$Day.of.tumor.harvest %in% g4_day & seurat$Experiment %in% g4_exp & seurat$Staining.panel.name %in% g4_panel]); length(g4)

dge_cells12 = c(g1, g2)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells12])); dge_min_mouse
dge_cells12 = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells12])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells12])[seurat$Biological.replicate[dge_cells12] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells12)
high_genes12 =  names(which(apply(seurat@assays$RNA@counts[, dge_cells12], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cells)))
dge_genes12 = setdiff(high_genes12, bad_features_DGE); length(dge_genes12)
dge_seurat12 = seurat[dge_genes12, dge_cells12]
# dge_seurat12 = NormalizeData(dge_seurat12, assay = dge_assay, normalization.method = "LogNormalize", scale.factor = 10000)
dge_g1_name = "m1"
dge_g2_name = "m2a"
Idents(dge_seurat12) = "Treatment"
dge_seurat.avg12 =  as.data.frame(AverageExpression(dge_seurat12, verbose = T, assays = dge_assay)[[dge_assay]]); dge_seurat.avg12$gene = rownames(dge_seurat.avg12)
dge_seurat.avg12$density =  get_density(dge_seurat.avg12[[dge_g1_name]], dge_seurat.avg12[[dge_g2_name]], n = 1000)
dge_seurat.wilcox12 = FindMarkers(dge_seurat12, ident.1 = dge_g2_name, ident.2 = dge_g1_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0)
dge_seurat.wilcox12$p_val_adj = p.adjust(dge_seurat.wilcox12$p_val, method = "fdr")

dge_cells34 = c(g3, g4)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells34])); dge_min_mouse
dge_cells34 = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells34])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells34])[seurat$Biological.replicate[dge_cells34] %in% x], size = dge_min_mouse, replace = F))); length(dge_cells34)
high_genes34 =  names(which(apply(seurat@assays$RNA@counts[, dge_cells34], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cells)))
dge_genes34 = setdiff(high_genes34, bad_features_DGE); length(dge_genes34)
dge_seurat34 = seurat[dge_genes34, dge_cells34]
# dge_seurat34 = NormalizeData(dge_seurat34, assay = dge_assay, normalization.method = "LogNormalize", scale.factor = 10000)
dge_g3_name = "FCGR_KO+m1"
dge_g4_name = "FCGR_KO+m2a"
Idents(dge_seurat34) = "Treatment"
dge_seurat.avg34 =  as.data.frame(AverageExpression(dge_seurat34, verbose = T, assays = dge_assay)[[dge_assay]]); dge_seurat.avg34$gene = rownames(dge_seurat.avg34)
dge_seurat.avg34$density =  get_density(dge_seurat.avg34[[dge_g3_name]], dge_seurat.avg34[[dge_g4_name]], n = 1000)
dge_seurat.wilcox34 = FindMarkers(dge_seurat34, ident.1 = dge_g4_name, ident.2 = dge_g3_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0)
dge_seurat.wilcox34$p_val_adj = p.adjust(dge_seurat.wilcox34$p_val, method = "fdr")

dge_logFC_cut = log2(1.5)
dge_seurat.wilcox12.sig.pos = dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox12.sig.pos)
dge_seurat.wilcox12.sig.neg = dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox12.sig.neg)
dge_seurat.wilcox34.sig.pos = dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox34.sig.pos)
dge_seurat.wilcox34.sig.neg = dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox34.sig.neg)

ngenes = 15
positive_genes12 = head(rownames(dge_seurat.wilcox12.sig.pos)[order(dge_seurat.wilcox12.sig.pos$avg_log2FC, decreasing = T)], ngenes)
negative_genes12 = head(rownames(dge_seurat.wilcox12.sig.neg)[order(dge_seurat.wilcox12.sig.neg$avg_log2FC, decreasing = F)], ngenes)
positive_genes34 = head(rownames(dge_seurat.wilcox34.sig.pos)[order(dge_seurat.wilcox34.sig.pos$avg_log2FC, decreasing = T)], ngenes)
negative_genes34 = head(rownames(dge_seurat.wilcox34.sig.neg)[order(dge_seurat.wilcox34.sig.neg$avg_log2FC, decreasing = F)], ngenes)
shared_genes_pos = intersect(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]))
shared_genes_neg = intersect(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]))
shared_genes_posneg = intersect(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]))
shared_genes_negpos = intersect(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]))
shared_genes = intersect(c(rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox12[dge_seurat.wilcox12$p_val_adj < alpha & dge_seurat.wilcox12$avg_log2FC < -dge_logFC_cut,])), c(rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC > dge_logFC_cut,]), rownames(dge_seurat.wilcox34[dge_seurat.wilcox34$p_val_adj < alpha & dge_seurat.wilcox34$avg_log2FC < -dge_logFC_cut,])))

fc_genes = intersect(rownames(dge_seurat.wilcox12), rownames(dge_seurat.wilcox34)); length(fc_genes)

ggde_table = data.frame(gene = fc_genes, wt = dge_seurat.wilcox12[fc_genes,]$avg_log2FC, ko = dge_seurat.wilcox34[fc_genes,]$avg_log2FC)
ggde_table$density =  get_density(ggde_table$wt, ggde_table$ko, n = 100)
nudge = 0.4
force = 30
size = 6

ggplot(ggde_table, aes(x = wt, y = ko)) +
  geom_hline(yintercept = c(-dge_logFC_cut, dge_logFC_cut), linetype = "dashed") +
  geom_vline(xintercept = c(-dge_logFC_cut, dge_logFC_cut), linetype = "dashed") +
  geom_point(aes(color = density), size = 2, shape = 19, stroke = 0) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes12, shared_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes12, shared_genes),], color = "red", aes(label = gene), segment.size  = 0.15, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes12, shared_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes12, shared_genes),], color = "red", aes(label = gene), segment.size  = 0.15, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes34, shared_genes),], color = "blue", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(positive_genes34, shared_genes),], color = "blue", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes34, shared_genes),], color = "blue", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% setdiff(negative_genes34, shared_genes),], color = "blue", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_pos,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_pos,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_neg,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_neg,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_posneg,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_posneg,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = -nudge, nudge_x = nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  geom_point(data = ggde_table[ggde_table$gene %in% shared_genes_negpos,], color = "black", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% shared_genes_negpos,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_y = nudge, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  
  theme(text = element_text(size = 25),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size =0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0)) +
  expand_limits(x = c(min(ggde_table$wt, ggde_table$ko), max(ggde_table$wt, ggde_table$ko)), y =  c(min(ggde_table$wt, ggde_table$ko), max(ggde_table$wt, ggde_table$ko))) +
  labs(x = "log2FC m2a/m1 WT", y = "log2FC m2a/m1 FCGRKO", title = dge_name) 
ggsave(filename = paste0("fig4_FCGRKO_fcfc_mon-mac_",dge_assay,"_combined.png"), width = 260, height = 230, units = "mm", dpi = 200)
ggsave(filename = paste0("fig4_FCGRKO_fcfc_mon-mac_",dge_assay,"_combined.pdf"), width = 260, height = 230, units = "mm", dpi = 200, bg = "transparent")

##############load data fig 5####
projdir = "/home/labs/amit/tomerlan/UCL_project"
id = "Seurat_in_vitro_all_split"
workdir = paste0("/home/labs/amit/tomerlan/UCL_project/", id)
setwd(workdir)
seurat = readRDS(file = "seurat.rds")
meta_plate = read.delim("metadata.txt")
markers = read.csv(paste0(projdir, "/markers_manuscript.csv"))[,1:2]
markers_all = read.csv(paste0(projdir, "/markers.csv"))[,1]

assay = "integrated2"
res = "0.7"
DefaultAssay(seurat) = assay
Idents(seurat) = paste0(assay,"_nn_res.", res)
cluster_ids = read.csv(paste0(workdir, "/", assay, "_annot.csv"), row.names = 1); 
cluster_annot = cluster_ids[,1]; names(cluster_annot) = rownames(cluster_ids)
cluster_color = cluster_ids[,2]; names(cluster_color) = as.character(cluster_annot)
seurat = RenameIdents(seurat, cluster_annot)

clusters_myeloid = setdiff(unique(seurat@active.ident), c("Treg", "Neutrophil"))

stress_genes= c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes =  unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))), "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"))
IFN_genes = unique(c(grep("Irf", rownames(seurat@assays$RNA@counts), v = T), grep("Ifi", rownames(rownames(seurat@assays$RNA@counts)), v = T), "Cd2","Cd3d", "Cmpk2","Cxcl10", "Isg15","Isg20","Oasl1" ,"Phf11b" ,"Plac8", "Rsad2","Rtp4","Sdf4", "Slfn4","Tnfrsf18" ,"Tnfrsf4","Tnfrsf9","Trex1","Usp18","Xaf1"))
ccl_genes = grep("Ccl", rownames(seurat@assays$RNA@counts), v = T)
MHC_genes =  c(grep("^H2-", rownames(seurat@assays$RNA@counts), v = T), "Cd74", "B2m")
hist_genes = grep("Hist", rownames(seurat@assays$RNA@counts), v = T)
comp_genes =  c(grep("^C1q", rownames(seurat@assays$RNA@counts), v = T))
ig_genes = c(grep("^Igj|^Igh|^Igk|^Igl", rownames(seurat@assays$RNA@counts), v = T))
hb_genes = c(grep("^Hba|^Hbb", rownames(seurat@assays$RNA@counts), v = T))
cc_gene_module = get(load("cc_gene_module.RData"))
markers = read.csv(paste0(workdir, "/", assay, "_markers_short.csv"), header = F)[, 1]



##############fig 5b#####
dot_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid]); length(dot_cells)

bad_features_pca = unique(c(cc_gene_module, hist_genes, cc_genes, stress_genes, ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T)))
seurat.reduction5 = RunPCA(seurat[, dot_cells], assay = assay, features = setdiff(seurat@assays[[assay]]@var.features, bad_features_pca))
seurat.reduction5 = RunUMAP(seurat.reduction5, assay = assay, dims = 1:50, a = 0.5, b = 2)

cluster_num = mapvalues(c(clusters_myeloid), from = c(clusters_myeloid), to = paste0(c(clusters_myeloid), " (", 1:length(c(clusters_myeloid)),")"))
cluster_num_color = cluster_color[c(clusters_myeloid)]; names(cluster_num_color) = cluster_num
ggumap = data.frame(seurat.reduction5@reductions$umap@cell.embeddings, Clust = seurat.reduction5@active.ident, Clust_no = mapvalues(seurat.reduction5@active.ident, from = c(clusters_myeloid), to = 1:length(c(clusters_myeloid))), Treatment = seurat.reduction5$Treatment, Day = seurat.reduction5$Day.of.tumor.harvest)
ggcenters = aggregate(ggumap[, 1:2], by = list(ggumap$Clust_no), FUN = median)
ggumap$Treatment = factor(ggumap$Treatment , levels = c("UT", "m1", "m2a","FCGRKO_UT","FCGRKO_m1", "FCGRKO_m2a" , "IFNARKO_UT" , "IFNARKO_m1" , "IFNARKO_m2a"))
ggumap$Clust = factor(ggumap$Clust, levels = c(clusters_myeloid))

#clusters
ggplot(ggumap) +
  theme(legend.key.size = unit(0.4, "cm")) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = factor(paste0(Clust, " (", Clust_no,")"), levels = cluster_num)), shape = 21, stroke = 0.05, size = 1, alpha = 1) +
  scale_fill_manual(values = cluster_num_color) +
  geom_text(data = ggcenters, aes(x = UMAP_1, y = UMAP_2, label = Group.1), size = 3, fontface = "bold") +
  guides(fill = guide_legend(title = "Cell population", override.aes = list(size = 3))) +
  theme_void() +
  # lims(x = c( min(ggumap[,1:2]), max(ggumap[,1:2])), y = c(min(ggumap[,1:2]), max(ggumap[,1:2]))) +
  theme(panel.background = element_rect(fill = "transparent", size = 0), 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom")
ggsave(file = paste0("fig5_umap_cluster.png"), width = 80, height = 90, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig5_umap_cluster.pdf"), width = 80, height = 90, units = "mm", dpi = 250, device = "pdf", bg = 'transparent')






##############fig 5c#####
dot_assay = "integrated2"
dotplot_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid]); length(dotplot_cells)

dotplot_markers = FindAllMarkers(seurat[,dotplot_cells], assay = dot_assay, only.pos = T, min.pct = 0.05, logfc.threshold = 0.25, base = 2, test.use = "wilcox", pseudocount.use = 1, features = seurat@assays[[dot_assay]]@var.features, max.cells.per.ident = 1000, random.seed = 1, min.cells.feature = 10)
dotplot_genes = markers; length(dotplot_genes)
dotplot_means = apply(as.matrix(clusters_myeloid), 1, function(x) rowMeans(seurat@assays[[dot_assay]]@data[dotplot_genes, names(seurat@active.ident[seurat@active.ident %in% x])])); colnames(dotplot_means) = clusters_myeloid
dotplot_clusts_hc = hclust(dist(t(dotplot_means)))
dotplot_clusts_order = rev(dotplot_clusts_hc$labels[dotplot_clusts_hc$order])
dotplot_genes_order = unique(c(unlist(apply(as.matrix(dotplot_clusts_order), 1, function(x) intersect(head(dotplot_markers[dotplot_markers$cluster == x, "gene"], 1000), dotplot_genes))))); length(dotplot_genes_order)
seurat@active.ident = factor(seurat@active.ident, levels = dotplot_clusts_order)

DotPlot(seurat[,dotplot_cells], features = rev(dotplot_genes_order), scale = T, assay = dot_assay, cols = c("white", "red"), cluster.idents = F) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
ggsave(file = paste0("fig5_dotmap_myeloid.png"), width = 120, height = 120, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig5_dotmap_myeloid.pdf"), width = 120, height = 120, units = "mm", dpi = 250, device = "pdf",  bg = "transparent")



##############fig 5d####
dge_alpha = 0.05
dge_logFC_cut = log2(1.25)
dge_min_UMI = 2
dge_min_cell = 5
dge_assay = "RNA"
dge_ngenes = 7

dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dge_cells = names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid]); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_wilcox_in_vitro = FindAllMarkers(seurat[,dge_cells], assay = dge_assay, only.pos = F, min.pct = 0.05, base = 2, test.use = "wilcox", pseudocount.use = 1, max.cells.per.ident = 1000, random.seed = 1, min.cells.feature = 10)
dge_wilcox_in_vitro$p_val_adj = p.adjust(dge_wilcox_in_vitro$p_val, method = "fdr")

#load in vivo
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
dge_clusters = c(cluster_macrophage, cluster_monocyte, "Mreg")
dge_lin = "Myeloid"
dge_days = c("d07")
dge_exp = c(3, 6)
dge_assay = "RNA"
dge_panel = c("Innate", "Broad")
dge_treat = c("m2a", "m1", "UT")
dge_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% dge_clusters & seurat$Treatment %in% dge_treat & seurat$Staining.panel.name %in% dge_panel & as.vector(seurat$Experiment) %in% dge_exp & seurat$Day.of.tumor.harvest %in% dge_days]); length(dge_cells)
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_wilcox_in_vivo = FindAllMarkers(seurat[,dge_cells], assay = dge_assay, only.pos = F, min.pct = 0.05, base = 2, test.use = "wilcox", pseudocount.use = 1, max.cells.per.ident = 1000, random.seed = 1, min.cells.feature = 10)
dge_wilcox_in_vivo$p_val_adj = p.adjust(dge_wilcox_in_vivo$p_val, method = "fdr")

##Ifn
dge_in_vivo = dge_wilcox_in_vivo[dge_wilcox_in_vivo$cluster == "Mac_Ifn" & abs(dge_wilcox_in_vivo$avg_log2FC) > dge_logFC_cut,] ; rownames(dge_in_vivo) = dge_in_vivo$gene
dge_in_vitro = dge_wilcox_in_vitro[dge_wilcox_in_vitro$cluster == "Mac_Ifn" & abs(dge_wilcox_in_vitro$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vitro) = dge_in_vitro$gene

dge_genes_joint = intersect(dge_in_vitro$gene, dge_in_vivo$gene)
dge_joint = data.frame(cbind(in_vitro = dge_in_vitro[dge_genes_joint, "avg_log2FC"], in_vivo = dge_in_vivo[dge_genes_joint,  "avg_log2FC"])); rownames(dge_joint) = dge_genes_joint
dge_joint$gene = rownames(dge_joint)

ggjoint = melt(dge_joint)
ggorder = ggjoint[order(ggjoint[ggjoint$variable == "in_vitro", "value"], decreasing = T),]
genes_pos = head(ggorder[ggorder$value > 0,"gene"], dge_ngenes)
genes_neg = tail(ggorder[ggorder$value < 0,"gene"], dge_ngenes)

genes_order_show = c(genes_pos, genes_neg)

ggplot(ggjoint[ggjoint$gene %in% genes_order_show,]) +
  geom_bar(aes(x = factor(gene, levels = genes_order_show), y = value, fill = variable), stat = "identity", position='dodge') +
  labs(y = "log2FC", title = "Mac Ifn - Mac Ifn") +
  scale_fill_manual(values = c("grey20", "grey80")) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.position = "none")
ggsave(filename = "fig5_vv_genes_Ifn.png", height = 60, width = 30, units = "mm", dpi = 250)
ggsave(filename = "fig5_vv_genes_Ifn.pdf", height = 60, width = 30, units = "mm", dpi = 250)


##Ccl5
dge_in_vivo = dge_wilcox_in_vivo[dge_wilcox_in_vivo$cluster == "Mon_Ccl5"  & abs(dge_wilcox_in_vivo$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vivo) = dge_in_vivo$gene
dge_in_vitro = dge_wilcox_in_vitro[dge_wilcox_in_vitro$cluster == "Mac_Ccl5" & abs(dge_wilcox_in_vitro$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vitro) = dge_in_vitro$gene

dge_genes_joint = intersect(dge_in_vitro$gene, dge_in_vivo$gene)
dge_joint = data.frame(cbind(in_vitro = dge_in_vitro[dge_genes_joint, "avg_log2FC"], in_vivo = dge_in_vivo[dge_genes_joint,  "avg_log2FC"])); rownames(dge_joint) = dge_genes_joint
dge_joint$gene = rownames(dge_joint)

ggjoint = melt(dge_joint)
ggorder = ggjoint[order(ggjoint[ggjoint$variable == "in_vitro", "value"], decreasing = T),]
genes_pos = head(ggorder[ggorder$value > 0,"gene"], dge_ngenes)
genes_neg = tail(ggorder[ggorder$value < 0,"gene"], dge_ngenes)

genes_order_show = c(genes_pos, genes_neg)
ggplot(ggjoint[ggjoint$gene %in% genes_order_show,]) +
  geom_bar(aes(x = factor(gene, levels = genes_order_show), y = value, fill = variable), stat = "identity", position='dodge') +
  labs(y = "log2FC", title = "Mac Ccl5 - Mon Ccl5") +
  scale_fill_manual(values = c("grey20", "grey80")) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.position = "none")
ggsave(filename = "fig5_vv_genes_Ccl5.png", height = 60, width = 30, units = "mm", dpi = 250)
ggsave(filename = "fig5_vv_genes_Ccl5.pdf", height = 60, width = 30, units = "mm", dpi = 250)


##Ccl7
dge_in_vivo = dge_wilcox_in_vivo[dge_wilcox_in_vivo$cluster == "Mac_Ccl7"  & abs(dge_wilcox_in_vivo$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vivo) = dge_in_vivo$gene
dge_in_vitro = dge_wilcox_in_vitro[dge_wilcox_in_vitro$cluster == "Mac_Ccl7" & abs(dge_wilcox_in_vitro$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vitro) = dge_in_vitro$gene

dge_genes_joint = intersect(dge_in_vitro$gene, dge_in_vivo$gene)
dge_joint = data.frame(cbind(in_vitro = dge_in_vitro[dge_genes_joint, "avg_log2FC"], in_vivo = dge_in_vivo[dge_genes_joint,  "avg_log2FC"])); rownames(dge_joint) = dge_genes_joint
dge_joint$gene = rownames(dge_joint)

ggjoint = melt(dge_joint)
ggorder = ggjoint[order(ggjoint[ggjoint$variable == "in_vitro", "value"], decreasing = T),]
genes_pos = head(ggorder[ggorder$value > 0,"gene"], dge_ngenes)
genes_neg = tail(ggorder[ggorder$value < 0,"gene"], dge_ngenes)

genes_order_show = c(genes_pos, genes_neg)
ggplot(ggjoint[ggjoint$gene %in% genes_order_show,]) +
  geom_bar(aes(x = factor(gene, levels = genes_order_show), y = value, fill = variable), stat = "identity", position='dodge') +
  labs(y = "log2FC", title = "Mac Ccl7 - Mac Ccl7") +
  scale_fill_manual(values = c("grey20", "grey80")) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.position = "none")
ggsave(filename = "fig5_vv_genes_Ccl7.png", height = 60, width = 30, units = "mm", dpi = 250)
ggsave(filename = "fig5_vv_genes_Ccl7.pdf", height = 60, width = 30, units = "mm", dpi = 250)

##Apoe
dge_in_vivo = dge_wilcox_in_vivo[dge_wilcox_in_vivo$cluster == "Mac_Arg1"  & abs(dge_wilcox_in_vivo$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vivo) = dge_in_vivo$gene
dge_in_vitro = dge_wilcox_in_vitro[dge_wilcox_in_vitro$cluster == "Mac_Apoe" & abs(dge_wilcox_in_vitro$avg_log2FC) > dge_logFC_cut,]; rownames(dge_in_vitro) = dge_in_vitro$gene

dge_genes_joint = intersect(dge_in_vitro$gene, dge_in_vivo$gene)
dge_joint = data.frame(cbind(in_vitro = dge_in_vitro[dge_genes_joint, "avg_log2FC"], in_vivo = dge_in_vivo[dge_genes_joint,  "avg_log2FC"])); rownames(dge_joint) = dge_genes_joint
dge_joint$gene = rownames(dge_joint)

ggjoint = melt(dge_joint)
ggorder = ggjoint[order(ggjoint[ggjoint$variable == "in_vitro", "value"], decreasing = T),]
genes_pos = head(ggorder[ggorder$value > 0,"gene"], dge_ngenes)
genes_neg = tail(ggorder[ggorder$value < 0,"gene"], dge_ngenes)

genes_order_show = c(genes_pos, genes_neg)
ggplot(ggjoint[ggjoint$gene %in% genes_order_show,]) +
  geom_bar(aes(x = factor(gene, levels = genes_order_show), y = value, fill = variable), stat = "identity", position='dodge') +
  labs(y = "log2FC", title = "Mac Apoe - Mac Arg1") +
  scale_fill_manual(values = c("grey20", "grey80")) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.position = "none")
ggsave(filename = "fig5_vv_genes_Apoe.png", height = 60, width = 30,  units = "mm", dpi = 250)
ggsave(filename = "fig5_vv_genes_Apoe.pdf", height = 60, width = 30,  units = "mm", dpi = 250)




##############fig 5e-g####
dot_panel =  c("Myeloid")
num_cells = names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid]); length(num_cells)

Mouse = as.vector(seurat$Biological.replicate[num_cells])
Experiment = as.vector(seurat$Experiment[num_cells])
Treatment = as.vector(seurat$Treatment[num_cells])
Day = as.vector(seurat$Day.of.tumor.harvest[num_cells])
Panel = as.vector(seurat$Staining.panel.name[num_cells])
ct = as.vector(seurat@active.ident[num_cells])
dot_data = data.frame(Mouse = Mouse, Treatment = Treatment, ct = ct, Experiment = Experiment, Day = Day)
dot_factors = c("Mouse", "Treatment", "Experiment")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$Mouse);  dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
dot_sum$Treatment = factor(dot_sum$Treatment, levels = c("UT", "m1", "m2a","FCGRKO_UT","FCGRKO_m1", "FCGRKO_m2a","IFNARKO_UT", "IFNARKO_m1", "IFNARKO_m2a"))
dot_plot = melt(dot_sum)


ggplot(dot_sum[dot_sum$Treatment %in% c("m1",  "UT",  "m2a"),], aes(x = Treatment, y = log2((TAM_Isg15 + 0.5)/(TAM_Apoe + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values =  c("grey100", "grey60", "grey10","azure", "azure3", "azure4","lightyellow", "lightyellow3", "lightyellow4")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1",  "m2a"), c("UT",  "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2(TAM Isg15/TAM Apoe fractions)") 
ggsave(file = paste0("fig5_mac_Ifn-Apoe_wt.png"), width = 30, height = 80, units = "mm", dpi = 150)
ggsave(file = paste0("fig5_mac_Ifn-Apoe_wt.pdf"), width = 30, height = 80, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2((dat$TAM_Isg15)/(dat$TAM_Apoe)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("UT", "m2a"),]
dat = cbind(value = log2((dat$TAM_Isg15)/(dat$TAM_Apoe)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

ggplot(dot_sum[dot_sum$Treatment %in% c("m1",  "UT",  "m2a"),], aes(x = Treatment, y = log2((TAM_Ccl5 + 0.5)/(TAM_Apoe + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values =  c("grey100", "grey60", "grey10","azure", "azure3", "azure4","lightyellow", "lightyellow3", "lightyellow4")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1",  "m2a"), c("UT",  "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 TAM Ccl5/TAM Apoe ratio") 
ggsave(file = paste0("fig5_mac_Ccl5-Apoe_wt.png"), width = 30, height = 80, units = "mm", dpi = 150)
ggsave(file = paste0("fig5_mac_Ccl5-Apoe_wt.pdf"), width = 30, height = 80, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2((dat$TAM_Ccl5)/(dat$TAM_Apoe)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("UT", "m2a"),]
dat = cbind(value = log2((dat$TAM_Ccl5)/(dat$TAM_Apoe)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)



ggplot(dot_sum[dot_sum$Treatment %in% c("m1",  "UT",  "m2a"),], aes(x = Treatment, y = log2((TAM_Ccl7 + 0.5)/(TAM_Apoe + 0.5)))) +
  geom_boxplot(aes(fill = Treatment), outlier.alpha = 0) +
  # stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  # stat_summary(fun.data=mean_se, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values =  c("grey100", "grey60", "grey10","azure", "azure3", "azure4","lightyellow", "lightyellow3", "lightyellow4")) +
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1",  "m2a"), c("UT",  "m2a")), test = "t.test", map_signif_level = T, family = "serif", textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("log2 TAM Ccl7/TAM Apoe ratio") 
ggsave(file = paste0("fig5_mac_Ccl7-Apoe_wt.png"), width = 30, height = 80, units = "mm", dpi = 150)
ggsave(file = paste0("fig5_mac_Ccl7-Apoe_wt.pdf"), width = 30, height = 80, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')

dat = dot_sum[dot_sum$Treatment %in% c("m1", "m2a"),]
dat = cbind(value = log2((dat$TAM_Ccl7+0.5)/(dat$TAM_Apoe+0.5)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)

dat = dot_sum[dot_sum$Treatment %in% c("UT", "m2a"),]
dat = cbind(value = log2((dat$TAM_Ccl7+0.5)/(dat$TAM_Apoe+0.5)), dat[,c("Mouse", "Treatment", "Experiment")])
anova = aov(data = dat, value ~ Treatment*Experiment)
summary(anova)
anova_test(data = dat, value ~ Treatment*Experiment)


r
##############fig 5h####
alpha = 0.05
dge_assay = "RNA"
dge_bad_features = c(ig_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))

unique(seurat@active.ident)
dge_name = "All"
dge_treat =  c("UT","m2a")
dge_day = c("12h")
dge_exp = unique(seurat$Experiment)
dge_panel = c("Myeloid")
dge_clusts = clusters_myeloid
dge_min_UMI = 2
dge_min_cell = 5

dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)
dge_min_mouse = min(table(seurat$Biological.replicate[dge_cells]))
dge_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[dge_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[dge_cells])[seurat$Biological.replicate[dge_cells] %in% x], size = dge_min_mouse, replace = F)))
dge_genes = setdiff(names(which(apply(seurat@assays$RNA@counts[, dge_cells], 1, function(x) length(which(x > dge_min_UMI)) > dge_min_cell))), dge_bad_features); length(dge_genes)
dge_seurat = seurat[, dge_cells]
# dge_seurat = NormalizeData(dge_seurat, assay = dge_assay, normalization.method = "RC", scale.factor = 1000000)
dge_seurat = dge_seurat[dge_genes,]
dge_g1_name = "UT"
dge_g2_name = "m2a"
Idents(dge_seurat) = "Treatment"
dge_seurat.wilcox = FindMarkers(dge_seurat, ident.1 = dge_g1_name, ident.2 = dge_g2_name, test.use = "wilcox", assay = dge_assay, pseudocount.use = 1, verbose = T, min.pct = 0.01, logfc.threshold = 0, slot = "data")
dge_seurat.wilcox$p_val_adj = p.adjust(dge_seurat.wilcox$p_val, method = "fdr")

dge_seurat.avg = as.data.frame(AverageExpression(dge_seurat, verbose = T, assays = dge_assay)[[dge_assay]]); dge_seurat.avg$gene = rownames(dge_seurat.avg)
dge_seurat.avg$density =  get_density(dge_seurat.avg[[dge_g1_name]], dge_seurat.avg[[dge_g2_name]], n = 1000)

dge_logFC_cut = log2(1.25)
dge_seurat.wilcox.sig.pos = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC > dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.pos)
dge_seurat.wilcox.sig.neg = dge_seurat.wilcox[dge_seurat.wilcox$p_val_adj < alpha & dge_seurat.wilcox$avg_log2FC < -dge_logFC_cut,]; nrow(dge_seurat.wilcox.sig.neg)

ngenes = 28
positive_genes = head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], ngenes)
negative_genes = head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], ngenes)

ggde_table = data.frame(gene = dge_seurat.avg$gene, x = dge_seurat.avg[[dge_g1_name]], y = dge_seurat.avg[[dge_g2_name]])
ggde_table$density =  get_density(ggde_table$x, ggde_table$y, n = 1000)

ggplot(ggde_table, aes(x = log1p(x), y = log1p(y), fill = density)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 3, shape = 19, stroke = 0, color = "grey") +
  geom_point(data = ggde_table[ggde_table$gene %in% c(positive_genes, negative_genes),], color = "red", size = 2, shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(negative_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 4,  nudge_y = 1.5, nudge_x = 0.5, max.overlaps = 100, force = 70) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% c(positive_genes),], aes(label = gene), color = "black", segment.color = "grey", segment.size  = 0.15, size = 4,  nudge_y = 0, nudge_x = 1.5, max.overlaps = 100, force = 70) +
  theme(axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", size = 0), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", size = 0), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", size = 0),
        legend.position = "none") +
  labs(x = paste(dge_g1_name, " log2 expression"),
       y = paste(dge_g2_name, " log2 expression"))
ggsave(filename = paste0("fig5_DGE_", dge_g1_name, "_", dge_g2_name, "_", dge_assay,".png"), width = 150, height = 160, units = "mm", dpi = 200)
ggsave(filename = paste0("fig5_DGE_", dge_g1_name, "_", dge_g2_name, "_", dge_assay,".pdf"), width = 150, height = 160, units = "mm", dpi = 200)


##############fig 5i####
heat_assay = "SCT"
heat_cells =  names(seurat$orig.ident[as.character(seurat@active.ident) %in% clusters_myeloid]); length(dotplot_cells)
heat_min_mouse = min(table(seurat$Biological.replicate[heat_cells]))
heat_cells = as.vector(apply(as.matrix(unique(seurat$Biological.replicate[heat_cells])), 1, function(x) sample(x = names(seurat$Biological.replicate[heat_cells])[seurat$Biological.replicate[heat_cells] %in% x], size = heat_min_mouse, replace = F)))
heat_genes = c(head(rownames(dge_seurat.wilcox.sig.pos)[order(dge_seurat.wilcox.sig.pos$avg_log2FC, decreasing = T)], 10), head(rownames(dge_seurat.wilcox.sig.neg)[order(dge_seurat.wilcox.sig.neg$avg_log2FC, decreasing = F)], 20))

heat_seurat = seurat[heat_genes, heat_cells]
Idents(heat_seurat) = "Treatment"
heat_seurat.avg =  as.data.frame(AverageExpression(heat_seurat, verbose = T, assays = heat_assay)[[heat_assay]]); heat_seurat.avg$gene = rownames(heat_seurat.avg); heat_seurat.avg = heat_seurat.avg[,-ncol(heat_seurat.avg)]
heat_column_order = c("UT", "m1", "m2a","FCGRKO_UT","FCGRKO_m1", "FCGRKO_m2a","IFNARKO_UT","IFNARKO_m1","IFNARKO_m2a")
heat_seurat.avg.scaled = t(apply(heat_seurat.avg, 1, scale)); 
colnames(heat_seurat.avg.scaled) = colnames(heat_seurat.avg)


pdf(paste0("fig5_heatmap_ifn.pdf"), width = 2, height = 4, pointsize = 1)
png(paste0("fig5_heatmap_ifn.png"), width = 100, height = 100, units = "mm", res = 200)
hm= Heatmap(heat_seurat.avg.scaled[,heat_column_order],
            
            name = "Normalized expression",
            col = colorRamp2(breaks = c(min(heat_seurat.avg.scaled), min(heat_seurat.avg.scaled)/2, 0, max(heat_seurat.avg.scaled)/2, max(heat_seurat.avg.scaled)), colors = c("slateblue4", "skyblue", "white", "tomato", "red3")),
            
            cluster_rows = T,
            show_row_names = T,
            row_names_gp = gpar(fontsize = 7),
            show_row_dend = F,
            
            cluster_columns = F,
            show_column_names = T,
            column_names_gp = gpar(fontsize = 7, rot = 45),
            column_names_rot = 45,
            show_column_dend = F,
            column_title_gp = gpar(fontsize = 7),
            
            heatmap_legend_param = list(direction = "horizontal"),
            
            cluster_column_slices = F,
            column_split = factor(c("WT", "WT", "WT", "FCGRKO", "FCGRKO", "FCGRKO", "IFNARKO", "IFNARKO","IFNARKO"), levels = c("WT", "FCGRKO", "IFNARKO")),
            column_labels = c("UT",  "m1","m2a","UT",  "m1","m2a","UT",  "m1","m2a"),
            
            rect_gp = gpar(col = "grey", lwd = 1),
            border = T,
            use_raster = F
            
)
draw(hm, heatmap_legend_side = "bottom")
dev.off()




##############fig 5j#####
cor_bad_features = c(cc_genes, hb_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T))
unique(seurat$Treatment)
dge_treat = c("UT","m1", "m2a") #c("FCGRKO_m2a" , "FCGRKO_UT"  ,"FCGRKO_m1")# 
dge_day = c("12h")
dge_exp = unique(seurat$Experiment)
dge_panel = c("Myeloid")
dge_clusts = clusters_myeloid
dge_assay = "RNA"
cor_genes = c("Fcgr1", "Fcgr4", "Fcgr2b" ,"Fcgr3")
dge_values = "data"
dge_pct = 0.1

dge_lig_genes = c("Fcgr1", "Fcgr4","Fcgr2b", "Fcgr3")
dge_cells = names(seurat@active.ident[seurat@active.ident %in% dge_clusts & seurat$Treatment %in% dge_treat & seurat$Day.of.tumor.harvest %in% dge_day & seurat$Experiment %in% dge_exp & seurat$Staining.panel.name %in% dge_panel]); length(dge_cells)

DotPlot(seurat[, dge_cells], features = dge_lig_genes, assay = dge_assay, cols = c("white", "red"), cluster.idents = T) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"))
ggsave(file = paste0("fig5_fcgr_bubble_cluster.png"), width = 80, height = 55, units = "mm", dpi = 250, device = "png")
ggsave(file = paste0("fig5_fcgr_bubble_cluster.pdf"), width = 80, height = 55, units = "mm", dpi = 250, device = "pdf", bg = "transparent")




##############fig 5i####
stats = as.data.frame(read_xlsx("/home/labs/amit/tomerlan/UCL_project/tumor.xlsx", sheet = 8, col_names = T))
class(stats$day) = "numeric"
class(stats$vol) = "numeric"
stats$treatment = factor(stats$treatment, levels = c("WT", "IfnarKO"))
stats$rescued = as.character(stats$rescued)

ggplot(stats, mapping = aes(x = day, y = vol, color =  rescued, group = paste(mouse, treatment))) +
  facet_grid(.~ treatment) +
  geom_line() +
  scale_color_manual(values = c("firebrick1", "dodgerblue")) +
  xlab("Days after tumor injection") +
  ylab("Tumor volume [mm^3]") +
  scale_x_continuous(breaks = unique(stats$day)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size= 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))
ggsave(file = "fig6_tumor_growth.png", width = 300, height = 120, units = "mm", dpi = 150, bg = "transparent")
ggsave(file = "fig6_tumor_growth.pdf", width = 300, height = 120, units = "mm", dpi = 150, bg = "transparent")

summary(aov(vol ~ day*treatment, data = stats[stats$treatment %in% c("WT", "IfnarKO"),]))

response = data.frame(response = c(0,0,0,0,1,0,1,1,1,1,1,1,1,1,1), Treatment = c(rep("IfnarKO", 5), rep("WT", 10)))
response$response = mapvalues(response$response, from =c(0,1), to = c("non-responsive", "responsive"))
ggplot(response) +
  geom_dotplot(binaxis = "y", stackdir = "center",stackratio = 1, aes(x = Treatment, y = response)) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_blank())
ggsave(file = "fig6_tumor_response.png", width = 70, height = 60, units = "mm", dpi = 150, bg = "transparent")
ggsave(file = "fig6_tumor_response.pdf", width = 70, height = 60, units = "mm", dpi = 150, bg = "transparent")



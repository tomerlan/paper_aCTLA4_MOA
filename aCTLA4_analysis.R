##############Libraries####
library(reshape)
library(ggplot2)
library(devtools)
library(plyr)
library(dplyr)
library(Seurat)

##############QC########
seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^mt-|^Mtmr")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = T, pt.size = 0.1)
ggsave(file = "QC.png", width = 300, height = 300, units = "mm", dpi = 150, device = "png")
seurat = subset(seurat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)

##############Scoring########
seurat[["total.UMI"]] = colSums(seurat@assays$RNA@counts)
seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^mt-|^Mtmr")
seurat[["ribo.genes"]] = PercentageFeatureSet(seurat, pattern = "^Rps|^Rpl")
s.genes =  unique(c(paste0(toupper(substr(tolower(as.character(cc.genes.updated.2019$s.genes)), 1, 1)), substr(tolower(as.character(cc.genes.updated.2019$s.genes)), 2, nchar(tolower(as.character(cc.genes.updated.2019$s.genes)))))))
g2m.genes = unique(c(paste0(toupper(substr(tolower(as.character(cc.genes.updated.2019$g2m.genes)), 1, 1)), substr(tolower(as.character(cc.genes.updated.2019$g2m.genes)), 2, nchar(tolower(as.character(cc.genes.updated.2019$g2m.genes))))))) 
seurat = CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
seurat[["cc.difference"]] =  seurat$S.Score - seurat$G2M.Score

##############Standard normalization########
assay = "RNA"
seurat = NormalizeData(seurat, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000)
seurat = FindVariableFeatures(seurat, assay = assay, selection.method = "vst", nfeatures = 3000)

vars_to_regress = c("S.Score", "G2M.Score", "percent.mt", "ribo.genes", "total.UMI")
seurat = ScaleData(seurat, assay = assay, features = rownames(seurat@assays[[assay]]@counts), vars.to.regress = vars_to_regress)


##############SCTtransform########
DefaultAssay(seurat) = "RNA"
vars_to_regress = c("S.Score", "G2M.Score", "percent.mt", "ribo.genes") #automatically regresses out depth
seurat = SCTransform(seurat, new.assay.name = "SCT", vars.to.regress = vars_to_regress, variable.features.n = 3000, do.scale = F, do.center = T, verbose = T, return.only.var.genes = F)


##############Data integration with SCTransform + cca + reference-based########
DefaultAssay(seurat) = "RNA"
integration_factor = "Experiment"
vars_to_regress = c("S.Score", "G2M.Score", "percent.mt", "ribo.genes")
seurat.list = SplitObject(seurat, split.by = integration_factor)
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] = SCTransform(seurat.list[[i]], verbose = T, do.correct.umi = T, return.only.var.genes = F, vars.to.regress = vars_to_regress, assay = "RNA", new.assay.name = "SCT", variable.features.n = 5000)
}
saveRDS(seurat.list, file = "seurat_list.rds")

integration_assay = "SCT"
seurat.features = SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000, assay = rep(integration_assay, length(seurat.list)))
seurat.list = PrepSCTIntegration(object.list = seurat.list, anchor.features = seurat.features, verbose = T, assay = integration_assay)

seurat.reduction = "cca"
seurat.dims = 1:50
# seurat.ref = 3
# seurat.list = lapply(X = seurat.list, FUN = RunPCA, features = seurat.features) # when using rpca
seurat.anchors = FindIntegrationAnchors(object.list = seurat.list, reduction = seurat.reduction, dims = seurat.dims, normalization.method = integration_assay, anchor.features = seurat.features, verbose = T) # , reference = which(names(seurat.list) == seurat.ref)
saveRDS(seurat.anchors, file = "seurat_anchors.rds")
seurat = IntegrateData(anchorset = seurat_anchors, dims = seurat.dims, normalization.method = integration_assay, new.assay.name = "integrated", verbose = T)


##############Clustering########
assay = "integrated"

stress_genes= c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes =  unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))))) #, "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"
IFN_genes = unique(c(grep("Irf", rownames(seurat@assays$RNA@counts), v = T), grep("Ifi", rownames(rownames(seurat@assays$RNA@counts)), v = T), "Cd2","Cd3d", "Cmpk2","Cxcl10", "Isg15","Isg20","Oasl1" ,"Phf11b" ,"Plac8", "Rsad2","Rtp4","Sdf4", "Slfn4","Tnfrsf18" ,"Tnfrsf4","Tnfrsf9","Trex1","Usp18","Xaf1"))
ccl_genes = grep("Ccl", rownames(seurat@assays$RNA@counts), v = T)
MHC_genes =  c(grep("^H2-", rownames(seurat@assays$RNA@counts), v = T), "Cd74", "B2m")
hist_genes = grep("Hist", rownames(seurat@assays$RNA@counts), v = T)
comp_genes =  c(grep("^C1q", rownames(seurat@assays$RNA@counts), v = T))
ig_genes = c(grep("^Igj|^Igh|^Igk|^Igl", rownames(seurat@assays$RNA@counts), v = T))
hb_genes = c(grep("^Hba|^Hbb", rownames(seurat@assays$RNA@counts), v = T))
cc_gene_module = get(load("cc_gene_module.RData"))
bad_features = unique(c(hist_genes, cc_genes, stress_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T)))

seurat = RunPCA(seurat, assay = assay, features = setdiff(seurat@assays[[assay]]@var.features, bad_features))

ElbowPlot(seurat, ndims = 50) + ggsave(file = paste0("elbow_",assay, ".png"), width = 200, height = 100, units = "mm", dpi = 150, device = "png")

seurat = FindNeighbors(seurat, assay = assay, dims = 1:50, do.plot = T, k.param = 50)

res = 2
seurat = FindClusters(seurat, graph.name = paste0(assay, "_nn"), resolution = res)


##############UMAP########
seurat = RunPCA(seurat, assay = assay, features = setdiff(seurat@assays[[assay]]@var.features, bad_features))
seurat = RunUMAP(seurat, assay = assay, dims = 1:50)



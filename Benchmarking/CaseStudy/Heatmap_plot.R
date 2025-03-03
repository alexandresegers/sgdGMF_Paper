#' ---
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(BiocParallel)
library(DelayedMatrixStats)
library(pheatmap)
nworkers <- 1

#' 
## -------------------------------------------------------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "TENxBrainData_preprocessed")
rownames(sce) <- rowData(sce)$Ensembl

dec.sce <- readRDS(file = "Output_files/dec.sce.RDS")
rowData(sce)$bio.var <- dec.sce$bio

hvg_rownames <- rownames(sce)[order(rowData(sce)$bio.var, decreasing = TRUE)][1:500]
hvg <- rowData(sce)$Symbol[order(rowData(sce)$bio.var, decreasing = TRUE)][1:500]


sgd_fit <- readRDS(file = "Output_files/CaseStudy_modelfit_full.RDS")
reducedDim(sce, type = "SGD_rescaled_10", withDimnames=F) <- prcomp(sgd_fit$U %*% t(sgd_fit$V), rank = 10)$x

sce <- sce[hvg_rownames,]


#' 
## -------------------------------------------------------------------------------------------------
set.seed(100)

graph <- buildSNNGraph(x = sce, 
                       use.dimred = "SGD_rescaled_10") 

#saveRDS(graph, file = "graph_10LF.RDS")


#' 
#' 
## -------------------------------------------------------------------------------------------------
set.seed(100)
leiden_clustering_00001 <- igraph::cluster_leiden(graph,resolution_parameter = 0.0002)


#' 
#' 
## -------------------------------------------------------------------------------------------------
gene_ids <- c("Mdk","Mki67", "Nde1", "Vim", "Fabp7", "Sox9", "Dbi", "Gad1", "Gad2", "Foxp2", "Syt6", 
  "Bcl11b", "Sst", "Lhx6","Neurod1", "Synpr", "Prox1","Sox4","Sox11", "Meis2",
  "Neurod6", "Satb2", "Crym", "Tnr", "Cck", "Ntf3", "Nr4a2", "Sstr2", "Mdk", "Eomes", "Pdgfra",
  "Vim", "Dbi", "Rgs5", "Igfbp7", "Hbb-bs", "Hba-a1", "Hba-a2", "Reln", "Pcp4", "Tbr1")

gene_ids <- intersect(gene_ids, hvg)

hvg_rownames_subset <- hvg_rownames[hvg %in% gene_ids]

lcounts <- as.matrix(logcounts(sce[hvg_rownames_subset, ]))


sce_pseudo_leiden_1 <- summarizeAssayByGroup(x = lcounts, 
                                                 ids = DataFrame(labels =
                                                              leiden_clustering_00001$membership),
                                                 statistics = "mean",
                                                 BPPARAM = MulticoreParam(nworkers))

    
    

#' 
## -------------------------------------------------------------------------------------------------
clusterMeans_leiden_1 <- assay(sce_pseudo_leiden_1)
colnames(clusterMeans_leiden_1) <- sce_pseudo_leiden_1$labels
rownames(clusterMeans_leiden_1) <- hvg[hvg %in% gene_ids]


#' 
#' 
## -------------------------------------------------------------------------------------------------
markers <- data.frame(
  gene = c("Meg3", "Mbp", "Pdgfra", "Dock2", "Rgs5", "Col1a2", "Aldh1l1", "Dnah11", "Mybpc1", "Deptor", "Rarb", "Satb2", "Tfap2d", "Fign", "Arap1", "Pax3", "Ntn1", "Pax2", "Slc6a3", "Fn1", "Tspan18", "Pde11a", "Dlx6os1", "Ntf3", "Itpr1", "Pvrl3", "Rorb", "Thsd7a", "Kcnk2", "Etv1", "Grik3", "Tle4", "Syt6", "Nr4a2", "Mki67", "Prox1", "Dock10", "Spock1", "Meis2", "Aoah", "Emr1", "Dab2", "Fyb", "Rbm47"),
  celltype = c("Neuronal", "Oligo", "OPC", "Macrophage/Microglia", "Endothelia/SMC", "VLMC", "Astro", "Ependyma", "OEC", "Mitra", "Medium Spiny Neuron", "CTX Pyr", "MTt Glut", "THAL Glut/Int", "Purkinje", "CB Int Progenitor", "CB Int Stellate/Golgi", "MD Glyc Int", "Nigral Dopaminergic", "HIPP/SUB Pyr", "HIPP Pyr/Granule", "SC Glut", "Migrating Int", "Pyr L2/3", "Pyr", "Pyr L2/3/4", "Pyr L4/5", "Pyr L4/5", "Pyr L5", "Pyr L5", "Pyr L5/6", "Pyr L5/6", "Pyr L6", "Pyr", "Granule Progenitors", "Granule", "Granule", "Pyr", "Pyr Progenitors", rep("Macrophage", 5)),
  study = rep("Rosenberg", 44)
)

markers <- rbind(markers,
                 data.frame(
                   gene = c("Neurod6", "Eomes", "Mki67", "Reln", "Satb2", "Fezf2", "Crym", "Bcl11b", "Sst", "Lhx6", "Adarb2", "Gad2", "Isl1", "Tcf7l2", "Hes5", "Aldh1l1", "Olig2", "Otx2", "Trem2", "Igfbp7"),
                   celltype = c("Excitatory Neurons", "Neuronal Progenitors", "Proliferating and Glia", "L1", "L2/4", "L5/6", "L5/6", "L5/6", "Interneurons", "Interneurons", "Interneurons",  "Interneurons/Striatal", "Striatal", "Thalamic", "Astrocytes", "Proliferating", "Oligodendrocytes", "Choroid Plexus", "Microglia", "Endothelial"),
                   study = rep("Loo", 20)
                 ))

markers <- rbind(markers,
                 data.frame(
                   gene = c("Rbfox3", "Slc17a6", "Slc17a7", "Slc17a8", "Gad1", "Gad2", "Reln", "C1qb", "P2ry12", "Aqp4", "Gja1", "Mbp", "Trf", "Plp1", "Tnr", "Cspg4", "Flt1", "Dcn", "Igfbpl1", "Rgs5", "Acta2", "Sox4", "Sox11", "Tgfbi", "Coch", "Ccdc153"),
                   celltype = c("Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Microglia/Macrophage", "Microglia/Macrophage", "Astrocyte", "Astrocyte", "Oligodendrocyte", "Choroid Plexus", "Oligodendrocyte", "Polydendrocyte", "Polydendrocyte", "Endothelial", "Fibroblast", "Fibroblast", "Mural", "Mural/Ependyma", "Neurogenesis & Mitosis", "Neurogenesis & Mitosis", "Choroid Plexus", "Choroid Plexus", "Ependyma"),
                   study = rep("Saunders", 26)
                 ))

markers <- rbind(markers,
                 data.frame(
                   gene = c("Tbr1", "Rasgrf2", "Pvrl3", "Cux2", "Rorb", "Plcxd2", "Thsd7a", "Kcnk2", "Cplx3", "Sulf2", "Foxp2", "Syt6", "Rprm", "Nr4a2", "Synpr", "Pcp4", "Gad1", "Pvalb", "Sst", "Htr3a", "Vip", "Reln", "Cck", "Npy", "Lhx6", "Calb2", "Pde1a", "Lphn2", "Kcnip2", "Rgs10", "Nov", "Cpne5", "Slc5a7", "Crh", "Pax6", "Cxcl14", "Gda", "Sema3e", "Aldh1l1", "Gfap", "Aqp4", "Serpinf1", "Gja1", "Mfge8", "Slco1c1", "Rnf122", "9630013A20Rik", "Itpr2", "Cnksr3", "Rhob", "Omg", "Klk6"),
                   celltype = c("Pyr", rep("Pyr L2/3", 3), rep("Pyr L5", 4), "Pyr L5/6", "Pyr L5/6", rep("Pyr L6", 3), "ClauPyr", "ClauPyr", "Pyr L5/6", rep("Interneurons", 22), rep("Astro", 7), rep("Oligo", 7)),
                   study = rep("Zeisel", 52)
                 ))

markers <- rbind(markers,
                 data.frame(
                   gene = c("Gad1", "Gad2", "Sst", "Chodl", "Rorb", "Fezf2", "Neurod6", "Akap7", "Htr3a", "Foxp2", "Mki67", "Top2a", "Fkbp5", "Tsc22d3", "Axl", "Snx6", "Snx2", "Dab2", "Ap1b1"),
                   celltype = c(rep("Inhibitory neurons", 4), rep("Excitatory Neurons", 2), rep("Inhibitory neurons", 4), rep("Macrophages", 2), rep("Adult Macrophages", 3), rep("Postnatal Macrophages", 4)),
                   study = c(rep("Allen", 6), rep("10X", 4), rep("Elizabeth", 2), rep("Li", 7))
                 ))

markers <- rbind(markers,
                 data.frame(
                   gene = c("Fabp7", "Vim", "Aldoc", "Dbi", "Ttyh1", "Slc1a3", "Ednrb", "Ngn2", "Gadd45g", "Neurod1", "Sstr2", "Vcam1", "Pax6", "Sox9", "Sox2", "Tlx", "Mdk", "Hopx", "Hes1", "Hes5", "Sox21", "Id4", "Nde1", "Nes"),
                   celltype = c(rep("Radial precursors", 7), rep("intermediate progenitor", 4), rep("Radial precursors", 13)),
                   study = "Yuzwa"
                 ))

markers <- rbind(markers,
                 data.frame(
                   gene = c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt", "Alas2", "Ube2l6"),
                   celltype = c(rep("DE mbk14", 6)),
                   study = "DE"
                 ))


markers <- markers[markers$gene %in% gene_ids,1:2]


markers[markers$gene %in% c("Mdk","Mki67", "Nde1", "Vim", "Fabp7", "Sox9", "Dbi"),"celltype"] <- 
  "Radial precursor"

markers[markers$gene %in% c("Gad1", "Gad2", "Foxp2", "Syt6", 
                                           "Bcl11b", "Sst", "Lhx6"),"celltype"] <- "Interneurons"

markers[markers$gene %in% c("Neurod1", "Synpr", "Prox1"),"celltype"] <- "Granule Cells"

markers[markers$gene %in% c("Sox4","Sox11", "Meis2", "Neurod6", "Satb2"),"celltype"] <- "Pyramidal progenitors"

markers[markers$gene %in% c("Crym", "Tnr", "Cck", "Ntf3", "Nr4a2"),"celltype"] <- "Pyramidal neurons"

markers[markers$gene %in% c("Sstr2", "Mdk", "Eomes"),"celltype"] <- "Neuronal progenitors"

markers[markers$gene %in% c("Pdgfra", "Vim", "Dbi"),"celltype"] <- "OPC"

markers[markers$gene %in% c("Rgs5", "Igfbp7"),"celltype"] <- "Endothelial/SMC/Mural"


markers[markers$gene %in% c("Hbb-bs", "Hba-a1", "Hba-a2"),"celltype"] <- "Blood"

markers[markers$gene %in% c("Reln", "Pcp4", "Tbr1"),"celltype"] <- "Cajal-Retzius cells"

markers[markers$gene %in% c("Mdk"), "celltype"] <- "Neuronal progenitors"
markers[markers$gene %in% c("Vim"), "celltype"] <- "Radial precursor"
markers[markers$gene %in% c("Vim", "Dbi"), "celltype"] <- "Radial precursor"



markers <- markers[!duplicated(markers),]
rownames(markers) <- markers$gene
markers <- markers[gene_ids,]

#' 
## -------------------------------------------------------------------------------------------------
metadata <- data.frame(loadings = (markers$celltype))
rownames(metadata) <- paste0(markers$gene,": ", markers$celltype)


colors_loadings <- clusterExperiment::bigPalette[1:11]
names(colors_loadings) <- (unique(markers$celltype))



#' 
## -------------------------------------------------------------------------------------------------

clusterMeans_leiden_1_restructured <- clusterMeans_leiden_1[gene_ids,]
rownames(clusterMeans_leiden_1_restructured) <- paste0(markers$gene,": ", markers$celltype)

p_clusterMeans_leiden_1 <- pheatmap(clusterMeans_leiden_1_restructured, 
                               scale = "row", 
                               annotation_legend = FALSE,
                               # annotation_col = data.frame(
                               #   cluster = sce_pseudo_leiden_1$labels, 
                                 #           ncells = sce_pseudo$ncells,
                                # row.names = sce_pseudo_leiden_1$labels),
                               annotation_row = metadata,
                               annotation_colors = list(
                                                        loadings = colors_loadings),
                               color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                               labels_row = paste0(markers$gene,": ", markers$celltype),
                               labels_col = paste0(sce_pseudo_leiden_1$labels),
                               annotation_names_row = F, annotation_names_col = F)

ggsave(p_clusterMeans_leiden_1, filename = "Figures/heatmap_clusters_10LF.pdf")

saveRDS(p_clusterMeans_leiden_1, file = "Figures/heatmap_clusters_10LF.RDS")

p_clusterMeans_leiden_2 <- pheatmap(clusterMeans_leiden_1_restructured, 
                               scale = "row", 
                               annotation_legend = FALSE,
                               #annotation_col = data.frame(
                                 #cluster = sce_pseudo_leiden_1$labels, 
                                 #           ncells = sce_pseudo$ncells,
                                # row.names = sce_pseudo_leiden_1$labels),
                               annotation_row = metadata,
                               annotation_colors = list(
                                                        loadings = colors_loadings),
                               color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                               labels_row = paste0(markers$gene),
                               labels_col = paste0(sce_pseudo_leiden_1$labels),
                               annotation_names_row = F, annotation_names_col = F)

ggsave(p_clusterMeans_leiden_2, filename = "Figures/heatmap_clusters_10LF_2.pdf")
saveRDS(p_clusterMeans_leiden_2, file = "Figures/heatmap_clusters_10LF_2.RDS")


#' 

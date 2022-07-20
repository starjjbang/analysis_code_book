library(Seurat)
library(scater)
library(SingleR)
library(dplyr)
library(magrittr)
#library(celldex)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(miQC)
library(SeuratWrappers)
library(AUCell)
library(pheatmap)
library(GSEABase )
library(corrplot)
library(RColorBrewer)

source(paste0(q_dir, 'code/single_function.R'))
single_cell_col = c("#F2B035","#879BC0","#32917F","#DF9644","#47458C",
                    "#1C4E90","#881B2B","#EB7A5B","#CF8336","#ECB1C8",
                    "#711F2D","#D43B3E","#187C65","#F89A8D","#CCD654",
                    "#D9CED2","#E9B65F","#7B727D","#797977","#B2C096",
                    "#A27FB6","#BACD74","#95CBBD","#807B85","#1B7B37","#EBF4F8","#163D5C")


label_extract<-function(f_in_dat){
  lab = c()
  for(i in 1:length(f_in_dat)){
    ch_line = grepl('sharm',f_in_dat[i])
    if(ch_line == TRUE){
      lab = c(lab, paste0(unlist(strsplit(f_in_dat[i],'_'))[1:2],collapse ='-'))
    }else{
      lab = c(lab, paste(unlist(strsplit(f_in_dat[i],'_'))[1:3],collapse = '-'))
    }
  }  
  return(lab)
}


sample_composition<-function(in_dat){
  plot_list <- list()
  for(i in 1:length(clu_id)){
    harmony.obj.plot = in_dat
    clu_pos = which(harmony.obj.plot@meta.data$seurat_clusters==clu_id[i] &  harmony.obj.plot@meta.data$gcell_type =="Oligodendrocyte "  )
    plot_01 = rep("Other",nrow(harmony.obj.plot@meta.data))
    pos_tmp1 = intersect(clu_pos,which(harmony.obj.plot@meta.data$orig.ident=="sharm"))
    plot_01[pos_tmp1] = "Sharm"
    pos_tmp2 = intersect(clu_pos,which(harmony.obj.plot@meta.data$orig.ident=="T2_nega"))
    plot_01[pos_tmp2] = "T2_nega"
    pos_tmp3 = intersect(clu_pos,which(harmony.obj.plot@meta.data$orig.ident=="T2_posi"))
    plot_01[pos_tmp3] = "T2_posi"
    harmony.obj.plot = AddMetaData(object = harmony.obj.plot, metadata = plot_01, col.name = 'plot_01')
    Idents(object = harmony.obj.plot) <- harmony.obj.plot@meta.data$plot_01
    plot_list[[i]] <-DimPlot(harmony.obj.plot,reduction = "umap") +ggtitle(paste0("Cluster ",clu_id[i],"(",length(clu_pos),")" ))+
      scale_color_manual(values=c("Sharm" = "#E08C66","T2_nega" = "#458181","T2_posi" = "#E6D6A4","Other" = "#F4F4EA"))+ theme(legend.position="none")
    
  }
  return(plot_list)
}

heat_per_cluster_compare<-function(clu_id,gn){
  plot_list <- list()
  for(j in 1:length(gn)){
    for_tmp = c()
    for(i in 1:length(clu_id)){
      re_f_each = paste0(re_f_dir1,clu_id[i])
      tmp = as.data.frame(read.table(file =paste0(re_f_each,'/','Culster_', clu_id[i],'_0_total_log2FC.csv'),sep='\t',header = TRUE))
      tmp1 = tmp[which(rownames(tmp)==gn[j]),]
      if(nrow(tmp1)==0){
        for_tmp = rbind(for_tmp,c(0,0,0)) 
      }else{
        for_tmp = rbind(for_tmp,tmp1)  
      }
    }
    rownames(for_tmp) = clu_id
    g_dat = cbind(for_tmp[,1],for_tmp[,3],for_tmp[,2])
    colnames(g_dat) = c(colnames(for_tmp)[1],colnames(for_tmp)[3],colnames(for_tmp)[2])
    rownames(g_dat) = clu_id
    
    library(pheatmap)
    paletteLength <- 50
    myColor <- colorRampPalette(c("#3B577C", "white", "#E5959A"))(paletteLength)
    myBreaks <- c(seq(min(g_dat), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(g_dat)/paletteLength, max(g_dat), length.out=floor(paletteLength/2)))
    x = pheatmap(g_dat,fontsize = 10,color=myColor,breaks=myBreaks,cluster_cols = FALSE,cluster_rows = FALSE,show_rownames=TRUE,main=gn[j])
    plot_list[[j]] <- x[[4]]
  }
  return(plot_list)
}

gene_cut_plot1<-function(gn){
  gene1<- FetchData(loigo.obj, vars = gn, slot = "data")
  range_exp = seq(floor(min(gene1)),ceiling(max(gene1)))
  clu_id = sort(unique(loigo.obj@meta.data$seurat_clusters))
  gr_id = unique(loigo.obj@meta.data$orig.ident)
  sub_clu_id = sort(unique(loigo.obj@meta.data$RNA_snn_res.0.6))
  for_clu = c()
  for_gr = c()
  for_sub_clu = c()
  for(i in 2:length(range_exp)){
    cell_id = rownames(gene1)[which(gene1[,1] <= range_exp[i])]
    tmp = loigo.obj@meta.data[which(rownames(loigo.obj@meta.data) %in% cell_id),]
    gr_tmp = table(tmp$orig.ident)
    gr_sort = sapply(gr_id, function(x) {ifelse( length(gr_tmp[which(names(gr_tmp)==x)]) > 0, gr_tmp[which(names(gr_tmp)==x)], 0 )})
    for_gr = rbind(for_gr,gr_sort)
    
    sub_clu_tmp= table(tmp$RNA_snn_res.0.6)
    sub_clu_sort = sapply(sub_clu_id, function(x) {ifelse( length(sub_clu_tmp[which(names(sub_clu_tmp)==x)]) > 0, sub_clu_tmp[which(names(sub_clu_tmp)==x)], 0 )})
    for_sub_clu = rbind(for_sub_clu,sub_clu_sort)
    
    
    clu_tmp = table(tmp$seurat_clusters)
    clu_sort = sapply(clu_id, function(x) {ifelse( length(clu_tmp[which(names(clu_tmp)==x)]) > 0, clu_tmp[which(names(clu_tmp)==x)], 0 )})
    for_clu = rbind(for_clu,clu_sort)
    
  }
  rownames(for_clu) = paste0("< ",range_exp[2:length(range_exp)])
  colnames(for_clu) = clu_id
  rownames(for_gr) = paste0("< ",range_exp[2:length(range_exp)])
  colnames(for_gr) = gr_id
  rownames(for_sub_clu) = paste0("< ",range_exp[2:length(range_exp)])
  colnames(for_sub_clu) = sub_clu_id
  
  for_clu1 = apply(for_clu, 1, function(x) x/mean(x))
  for_gr1 = apply(for_gr, 1, function(x) x/mean(x))
  for_sub_clu1 = apply(for_sub_clu, 1, function(x) x/mean(x))
  
  png(filename = paste0(fig_dir,gn,"_marker_fc_window.png"), width = 3500, height = 3800,
      units = "px", bg = "white", res = 300)
  col_p  = c("#B3E2CD","#FDCDAC","#CBD5E8","#F3CAE4","#E6F5C8","#FFF2AE","#F1E2CC","#CCCCCC")
  par(mfrow = c(4,2))
  par(mar = c(2.1, 4.1, 4.1, 2.1))
  plot(density(gene1[,1]), main =gn ,xlab="log-normed version of counts")
  barplot(for_clu1,col =col_p, ylab="Rate of cells in each cutoff",xlab="log2value", main = paste0(gn,'-',"cluster"))
  par(mar=c(0, 0, 0, 0))
  plot.new()
  plot.new()
  legend("top", legend=paste0("Cluster",clu_id), col=col_p,cex=0.8, box.lty=0,pch = 16,ncol = 3)
  par(mar = c(2.1, 4.1, 2.1, 2.1))
  barplot(for_gr1,col =col_p, ylab="Rate of cells in each cutoff",xlab="log2value",main = paste0(gn,'-',"Group"))
  barplot(for_sub_clu1,col =col_p, ylab="Rate of cells in each cutoff",xlab="log2value", main = paste0(gn,'-',"Sub-cluster"))
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("top", legend=gr_id, col=col_p,cex=0.8, box.lty=0,pch = 16,ncol = 3)
  plot.new()
  legend("top", legend=paste0("Cluster",sub_clu_id), col=col_p,cex=0.8, box.lty=0,pch = 16,ncol = 3)
  dev.off()
}

val_check<-function(in_dat,g_val){
  tmp = in_dat[which(rownames(in_dat)== g_val),2]
  if(length(tmp) >0){
    return(tmp)
  }else{
    return(0)
  }
}

cal_sep<-function(mer_in_dat,sample_na){
  q_dir = '/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/'
  rds_dir = '/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/rds/1_each_snslaysis/'
  fig_dir = '/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/plot/4_each_analysis/'
  dir.create(file.path(rds_dir, sample_na), showWarnings = FALSE)
  dir.create(file.path(fig_dir, sample_na), showWarnings = FALSE)
  rds_dir = paste0(rds_dir, sample_na,'/')
  fig_dir = paste0(fig_dir,sample_na,'/')
  #QC
  data.fit = mer_in_dat[["percent.mt"]] = PercentageFeatureSet(object = mer_in_dat, pattern = "^Mt-")
  data.fit = RunMiQC(mer_in_dat, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, model.slot = "flexmix_model")
  data.fit.qc.obj = subset(data.fit, miQC.keep == "keep")
  saveRDS(data.fit.qc.obj, paste0(rds_dir,paste0("1_",sample_na,"_after_QC.rds")))
  #Filter
  counts <- GetAssayData(data.fit.qc.obj, assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% rownames(data.fit.qc.obj)[which(grepl("^Mt-", rownames(data.fit.qc.obj)))] )),]
  data.fit.qc.obj1 <- subset(data.fit.qc.obj, features = rownames(counts))
  saveRDS(data.fit.qc.obj1, paste0(rds_dir,"2_",sample_na,"_filter_QC.rds"))
  #Cell cycle
  data.fit.qc.obj <- CellCycleScoring(object = data.fit.qc.obj1, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"3_",sample_na,"_after_cell_cycle.rds"))
  #Doublet Proportion Estimate
  #data.fit.qc.obj = readRDS(paste0(rds_dir,"3_",sample_na,"_after_cell_cycle.rds"))
  data.fit.qc.obj = NormalizeData(object = data.fit.qc.obj, normalization.method = "LogNormalize")
  data.fit.qc.obj = FindVariableFeatures(object = data.fit.qc.obj, selection.method = "vst")
  data.fit.qc.obj= ScaleData(object = data.fit.qc.obj, features = rownames(data.fit.qc.obj))
  data.fit.qc.obj = RunPCA(data.fit.qc.obj, verbose = F, npcs = 20)
  
  sweep.res <- paramSweep_v3(data.fit.qc.obj)
  sweep.stats <-summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  pk_val = as.numeric(as.character(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),2]))
  print(pk_val)#I notice that in the method part you mention large pK(>0.1)
  annotations <- data.fit.qc.obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp <- round(ncol(data.fit.qc.obj) * pk_val)  # expect 4% doublets
  data.fit.qc.obj <- doubletFinder_v3(data.fit.qc.obj, pN = 0.25, pK = pk_val, nExp = nExp, PCs = 1:10)
  DF.name = colnames(data.fit.qc.obj@meta.data)[grepl("DF.classification", colnames(data.fit.qc.obj@meta.data))]
  Idents(object = data.fit.qc.obj) <- DF.name
  data.fit.qc.obj1 =  subset(data.fit.qc.obj, idents = "Singlet")
  saveRDS(data.fit.qc.obj1, paste0(rds_dir,"4_",sample_na,"_after_Doublet.rds"))
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"4_1_",sample_na,"_non_filter_Doublet.rds"))
  print("QC END")
  data.fit.qc.obj = JackStraw(object = data.fit.qc.obj1, num.replicate = 100, dims = 50)
  data.fit.qc.obj = ScoreJackStraw(object = data.fit.qc.obj, dims = 1:20)
  
  data.fit.qc.obj = readRDS(paste0(rds_dir,"4_",sample_na,"_after_Doublet.rds"))
  data.fit.qc.obj <- FindNeighbors(object = data.fit.qc.obj, dims = 1:20, force.recalc = T, reduction = "pca")
  data.fit.qc.obj <- FindClusters(object = data.fit.qc.obj, resolution = 0.6)
  data.fit.qc.obj <- RunUMAP(object = data.fit.qc.obj, reduction = "pca", dims = 1:12)
  data.fit.qc.obj <- RunTSNE(object = data.fit.qc.obj, reduction = "pca", dims = 1:12)
  #DimPlot(data.fit.qc.obj,reduction = "umap")
  
  
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"5_",sample_na,"_after_reduction.rds"))
  print("After_reduction END")
  
  cell_marker_lst = as.data.frame(read.table( paste0(q_dir,'1_raw_data/marker_gene_lst_1.txt' ),sep='\t'))
  ct_id = unique(cell_marker_lst[,2])
  g_lst = list()
  for(i in 1:length(ct_id)){
    g_lst[[i]] = cell_marker_lst[which(cell_marker_lst[,2]==ct_id[i]),1]
  }
  names(g_lst) = ct_id
  
  all.sets <- lapply(names(g_lst), function(x) {
    GeneSet(g_lst[[x]], setName=x)        
  })
  
  
  source('/home/starjjbang/project/1_prof_kim_lab_deg/1_single_cell_project/code/single_function_ver02_0412.R')
  
  all.sets <- GeneSetCollection(all.sets)
  rankings <- AUCell_buildRankings(data.fit.qc.obj@assays$RNA@counts, plotStats=FALSE, verbose=FALSE)
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank = nrow(rankings) * 0.05)
  cell_assignment <- AUCell::AUCell_exploreThresholds(
    cell.aucs, plotHist = TRUE, 
    nCores = 1, assignCells = TRUE
  )
  
  results <- t(getAUC(cell.aucs))
  cellsAssigned <- lapply(cell_assignment, function(x) x$assignment)
  assignmentTable <- reshape2::melt(cellsAssigned, value.name = "cell")
  colnames(assignmentTable)[2] <- "geneSet"
  new.labels <- colnames(results)[max.col(results)]
  data.fit.qc.obj <- AddMetaData(object = data.fit.qc.obj,metadata = new.labels, col.name = 'gcell_type')
  tab <- table(new.labels, data.fit.qc.obj@meta.data$gcell_type)
  assignmentMat <- table(assignmentTable[, "geneSet"], assignmentTable[, "cell"])
  miniAssigMat <- assignmentMat[, sample(1:ncol(assignmentMat), min(500, ncol(assignmentMat)))]
  saveRDS(data.fit.qc.obj, paste0(rds_dir,paste0("3_",sample_na,"_QC_RD_cell_an.rds")))
  print("Cell Annotation END")
  
  library(RColorBrewer)
  n <- 20
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  p1<- VizDimLoadings(data.fit.qc.obj, dims = 1:2, reduction = "pca")
  p2<- DimPlot(data.fit.qc.obj, reduction = "pca",group.by = "orig.ident")
  p3 = DimPlot(object = data.fit.qc.obj,label=TRUE)+ ggplot2::theme(legend.position = "none")
  p4 = DimPlot(object = data.fit.qc.obj, reduction = "umap", label = T,label.size = 4,cols = col_vector,group.by ="gcell_type" )+
    guides(guide_legend(ncol=3,nrow=3,byrow=TRUE)) +theme(legend.title = element_text( size=6), legend.text=element_text(size=6))+ ggtitle("")
  cell_names = rownames(data.fit.qc.obj@meta.data[which(data.fit.qc.obj@meta.data$gcell_type=="Oligodendrocyte "),])
  p5 = DimPlot(object = data.fit.qc.obj, cells.highlight = cell_names, order = TRUE, 
               split.by = "orig.ident",reduction = "umap",group.by  = "gcell_type")+
    scale_color_manual(values = c("grey", "#5880B4"),labels = c("Other","Oligodendrocyte"))+ ggtitle("")
  olg_lst = cell_marker_lst[which(cell_marker_lst[,2] == "Oligodendrocyte "),1]
  p6 = FeaturePlot(data.fit.qc.obj, features = olg_lst,ncol =3,cols =c("#F1ECDA","#305388"))
  
  png(filename = paste0(fig_dir,paste0("1_",sample_na,"_ct_map.png")), width = 5500, 
      height = 6000,units = "px", bg = "white", res = 300)
  grid.arrange(arrangeGrob(p1[[1]], p1[[2]],p2,ncol=3,nrow=1),
               arrangeGrob(p3,p4,p5,ncol=3,nrow=1),
               arrangeGrob(p6[[1]],p6[[2]],p6[[3]],p6[[4]],p6[[5]],ncol=3,nrow=2),
               nrow=3,heights = c(2,2,2.5))
  dev.off()
  
  png(filename = paste0(fig_dir,paste0("2_",sample_na,"_ct_cal.png")), width = 5500, 
      height = 6000,units = "px", bg = "white", res = 300)
  par(mfrow = c(3, 3))
  cell_assignment <- AUCell::AUCell_exploreThresholds(cell.aucs, plotHist = TRUE, nCores = 1, assignCells = TRUE)
  dev.off()
  print("END")
}

cluster_num_sel <-function(in_dat){
  ct_tmp = c()
  c_id = sort(unique(in_dat@meta.data$seurat_clusters))
  for(i in 1:length(c_id)){
    tmp = in_dat@meta.data[which(in_dat@meta.data$seurat_clusters==c_id[i]),]
    ct = table(tmp$gcell_type)/nrow(tmp)
    ct1 = length(which(names(ct)=="Oligodendrocyte "))
    if(ct1 > 0){
      ct_tmp = c(ct_tmp,ct[which(names(ct)=="Oligodendrocyte ")])
    }else{
      ct_tmp = c(ct_tmp, 0) 
    }
  }
  names(ct_tmp) = c_id
  return(as.numeric(names(ct_tmp)[which(ct_tmp>=0.5)])  )
}
pie_plot1<-function(in_dat,main_title){
  mycols <- c("#D6E6E3","#F4ECD6")
  g_df = data.frame(val = in_dat,class=names(in_dat),
                    y_pos =cumsum(in_dat) - 0.5*in_dat )
  ggplot(g_df, aes(x = 2, y = val.Freq, fill = class)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0)+
    geom_text(aes(y = y_pos.Freq, label = paste0(round(in_dat/sum(in_dat)*100,2),'%')), color = "navy")+
    scale_fill_brewer(palette=3) +
    theme_void() + xlim(0.5, 2.5) + ggtitle(main_title) + theme(plot.margin = unit(c(0,0,0.5,0), "cm"))
}


############################################################################
# 3. QC
############################################################################ 
qc_fun <-function(in_rds_dir){
  mer_in_dat = readRDS( in_rds_dir)
  data.fit = mer_in_dat[["percent.mt"]] = PercentageFeatureSet(object = mer_in_dat, pattern = "^Mt-") 
  data.fit = RunMiQC(mer_in_dat, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",  posterior.cutoff = 0.75, model.slot = "flexmix_model")
  saveRDS(data.fit, paste0(rds_dir,"1_miQC.rds"))
}
############################################################################
# 4. FILTER
############################################################################ 

filter_fn <-function(in_rds_dir){
  data.fit  = readRDS(in_rds_dir)
  data.fit.qc.obj = subset(data.fit, miQC.keep == "keep")
  counts <- GetAssayData(data.fit.qc.obj, assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% rownames(data.fit.qc.obj)[which(grepl("^Mt-", rownames(data.fit.qc.obj)))] )),]
  data.fit.qc.obj1 <- subset(data.fit.qc.obj, features = rownames(counts))
  saveRDS(data.fit.qc.obj1, paste0(rds_dir,"2_after_QC.rds"))
  rm(list=ls())
  
}

############################################################################
# 5. Cell cycle
############################################################################ 
cell_cy_fn<-function(in_rds_dir){
  data.fit.qc.obj <- readRDS( in_rds_dir)
  data.fit.qc.obj <- CellCycleScoring(object = data.fit.qc.obj, 
                                      g2m.features = cc.genes$g2m.genes, 
                                      s.features = cc.genes$s.genes)
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"3_after_cell_cycle.rds"))
}

############################################################################
# 6. PCA analysis
############################################################################ 
pca_an_fn<-function(in_rds_dir){
  data.fit.qc.obj = readRDS(in_rds_dir)
  data.fit.qc.obj = NormalizeData(object = data.fit.qc.obj, normalization.method = "LogNormalize")
  data.fit.qc.obj = FindVariableFeatures(object = data.fit.qc.obj, selection.method = "vst")
  data.fit.qc.obj= ScaleData(object = data.fit.qc.obj, features = rownames(data.fit.qc.obj))
  data.fit.qc.obj = RunPCA(data.fit.qc.obj, verbose = F, npcs = 20)
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"4_after_pca.rds"))
}


############################################################################
# 7. Doublet Proportion Estimate
############################################################################
Doublet_an_fn <-function(in_rds_dir){
  library(DoubletFinder)
  data.fit.qc.obj = readRDS(in_rds_dir)
  sweep.res <- paramSweep_v3(data.fit.qc.obj)
  sweep.stats <-summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  saveRDS(bcmvn, paste0(rds_dir,"5_bcmvn.rds"))
  pk_val = as.numeric(as.character(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),2]))
  print(pk_val)#I notice that in the method part you mention large pK(>0.1)
  annotations <- data.fit.qc.obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp <- round(ncol(data.fit.qc.obj) * pk_val)  # expect 4% doublets
  data.fit.qc.obj <- doubletFinder_v3(data.fit.qc.obj, pN = 0.25, pK = pk_val, nExp = nExp, PCs = 1:10)
  DF.name = colnames(data.fit.qc.obj@meta.data)[grepl("DF.classification", colnames(data.fit.qc.obj@meta.data))]
  colnames(data.fit.qc.obj@meta.data)[which(names(data.fit.qc.obj@meta.data) == DF.name)] <- "df_class"
  data.fit.qc.obj1 =  subset(data.fit.qc.obj, df_class == "Singlet")
  saveRDS(data.fit.qc.obj1, paste0(rds_dir,"5_after_Doublet.rds"))
  saveRDS(data.fit.qc.obj, paste0(rds_dir,"5_1_non_filter_Doublet.rds"))  
}
############################################################################
# 7. harmony batch correction
###########################################################################
harmony_an_fn <-function(in_rds_dir){
  library(harmony)
  d.data.fit = readRDS(in_rds_dir)
  d.data.fit = JackStraw(object = d.data.fit, num.replicate = 100, dims = 50)
  d.data.fit = ScoreJackStraw(object = d.data.fit, dims = 1:20)
  harmony.obj <- RunHarmony(object = d.data.fit, group.by.vars = "orig.ident")
  harmony.obj <- FindNeighbors(object = harmony.obj, dims = 1:20, force.recalc = T, reduction = "harmony")
  harmony.obj <- FindClusters(object = harmony.obj, resolution = 0.6)
  harmony.obj <- RunUMAP(object = harmony.obj, reduction = "harmony", dims = 1:12)
  harmony.obj <- RunTSNE(object = harmony.obj, reduction = "harmony", dims = 1:12)
  saveRDS(harmony.obj, paste0(rds_dir,"6_after_harmony.rds"))
  
}
############################################################################
# 8. cell annotation
############################################################################ 

cell_an_fn <-function(in_rds_dir){
  harmony.obj = readRDS(in_rds_dir)
  cell_marker_lst = as.data.frame(read.table( paste0(q_dir,'1_raw_data/marker_gene_lst_1.txt' ),sep='\t'))
  ct_id = unique(cell_marker_lst[,2])
  g_lst = list()
  for(i in 1:length(ct_id)){
    g_lst[[i]] = cell_marker_lst[which(cell_marker_lst[,2]==ct_id[i]),1]
  }
  names(g_lst) = ct_id
  all.sets <- lapply(names(g_lst), function(x) {
    GeneSet(g_lst[[x]], setName=x)        
  })
  all.sets <- GeneSetCollection(all.sets)
  rankings <- AUCell_buildRankings(harmony.obj@assays$RNA@counts, plotStats=FALSE, verbose=FALSE,splitByBlocks=TRUE)
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank = nrow(rankings) * 0.05)
  results <- t(getAUC(cell.aucs))
  head(results)
  cell_assignment <- AUCell::AUCell_exploreThresholds(cell.aucs, plotHist = TRUE, nCores = 1, assignCells = TRUE)
  cellsAssigned <- lapply(cell_assignment, function(x) x$assignment)
  assignmentTable <- reshape2::melt(cellsAssigned, value.name = "cell")
  colnames(assignmentTable)[2] <- "geneSet"
  head(assignmentTable)
  new.labels <- colnames(results)[max.col(results)]
  harmony.obj <- AddMetaData(object = harmony.obj,metadata = new.labels, col.name = 'gcell_type')
  tab <- table(new.labels, harmony.obj@meta.data$gcell_type)
  saveRDS(harmony.obj, paste0(rds_dir,"7_cell_type_gene_obj.rds"))
  assignmentMat <- table(assignmentTable[, "geneSet"], assignmentTable[, "cell"])
  miniAssigMat <- assignmentMat[, sample(1:ncol(assignmentMat), min(500, ncol(assignmentMat)))]
  png(filename = paste0(fig_dir,"6_cell_type.png"),width = 2400, height = 2000,units = "px", bg = "white", res = 300)
  pheatmap(miniAssigMat, cluster_rows = FALSE, show_colnames = FALSE)
  dev.off()
  png(filename = paste0(fig_dir,"7_ct_map.png"), width = 2400, height = 2000,units = "px", bg = "white", res = 300)
  par(mfrow=c(3,3))
  AUCell_exploreThresholds(cell.aucs, plotHist=TRUE, assign=TRUE) 
  #fullheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)
  dev.off()
}
############################################################################
# stat(#cell plot)
############################################################################ 
stat_plot<-function(in_rds_dir){
  mer_in_dat = readRDS(in_rds_dir)
  g_dat = table(mer_in_dat$orig.ident)
  tab_tmp = data.frame(Statistics = c("#Cells"),sharm = c(g_dat[1] ),T2_nega = c(g_dat[2] ),T2_pos = c(g_dat[3] ))
  tab<- ggtexttable(tab_tmp, rows = NULL, theme = ttheme("light"))
  lab_ex = label_extract(names(mer_in_dat$orig.ident))
  tab_tmp1 = data.frame(Statistics = c("#Cells"),sharm = c(g_dat[1] ),T2_nega = c(g_dat[2] ),T2_pos = c(g_dat[3] ))
  tab_tmp1= t(data.frame(c(Statistics = c("#Cells"),data.frame(table(lab_ex))$Freq)))
  colnames(tab_tmp1) = c("Statistics",names(table(lab_ex)))
  tab1<- ggtexttable(tab_tmp1, rows = NULL, theme = ttheme("light"))
  png(filename = paste0(fig_dir,"0_st_graph.png"),width = 2300, height = 500,units = "px", bg = "white", res = 300)
  grid.arrange(tab,tab1,nrow=2)
  dev.off()
}
############################################################################
# QC plot
############################################################################ 
qc_plot<-function(in_rds_dir){
  data.fit = readRDS(in_rds_dir)
  p = VlnPlot(object = data.fit, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p1 = PlotMiQC(data.fit, color.by = "miQC.probability") + scale_color_gradient(low = "grey", high = "purple")
  p2 = PlotMiQC(data.fit, color.by = "miQC.keep")
  p3 = VlnPlot(object = data.fit, features = c("nFeature_RNA", "percent.mt","nCount_RNA"), group.by = "miQC.keep")
  p5 = pie_plot1(table(data.fit@meta.data$miQC.keep), "Filter data using MiQC")
  qc_dat = table(data.fit@meta.data$miQC.keep)
  tab_tmp = data.frame(Attribute = c("QC by MiQC"),Fail = c(qc_dat[1]), Pass = c(qc_dat[2]),Total = c(sum(qc_dat)))
  tab<- ggtexttable(tab_tmp, rows = NULL, theme = ttheme("light"))
  tab<- tab %>%
    tab_add_title(text = paste0("Analysis tools : Seurat with MiQC
  posterior.cutoff = 0.75
  model.slot: flexmix_model"),
                  face = "plain", size = 10) %>%
    tab_add_title(text = "Statistics of Filtering Condition", face = "bold", padding = unit(0.1, "line")) %>%
    tab_add_footnote(text = "*AKB-AI", size = 10, face = "italic")
  
  png(filename = paste0(fig_dir,"1_QC_miqc_graph.png"), 
      width = 3500, height = 3800,units = "px", bg = "white", res = 300)
  grid.arrange(arrangeGrob(p[[1]],p[[2]],p[[3]],ncol=3,nrow=1),
               arrangeGrob(p1, p2, ncol=2, nrow=1),
               arrangeGrob(p3[[1]],p3[[2]],p3[[3]],ncol=3,nrow=1),
               arrangeGrob(p5, tab,ncol=2,nrow=1),
               nrow=4,heights = c(2,1.5,2,1.5) )
  dev.off()
}

############################################################################
#Cell cycle plot
############################################################################ 
#we assign each cell a score, based on its expression of G2/M and S phase markers. 
#ylab = "average expression of this gene list "
cc_plot<-function(in_rds_dir){
  data.fit.qc.obj <- readRDS(in_rds_dir)
  p_tmp = VlnPlot(data.fit.qc.obj, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", ncol = 2, pt.size = 0.1)+ 
    labs(y = "cell cycle score")
  p1 = p_tmp[[1]] + theme(plot.title = element_text(size=10),axis.text.x=element_text(size=rel(0.8)),
                          axis.text.y=element_text(size=rel(0.8)),axis.title.x=element_blank())+ 
    labs(y = "cell cycle score")
  p2 = p_tmp[[2]] + theme(plot.title = element_text(size=10),axis.text.x=element_text(size=rel(0.8)),
                          axis.text.y=element_text(size=rel(0.8)),axis.title.x=element_blank())+ 
    labs(y = "cell cycle score")
  in_dat = table(data.fit.qc.obj@meta.data$Phase)
  g_df = data.frame(val = as.numeric(in_dat),class=names(in_dat),
                    y_pos =cumsum(as.numeric(in_dat)) - 0.5*as.numeric(in_dat) )
  g_df$class = paste0(names(in_dat),'(',round(in_dat/sum(in_dat)*100,2),'%)')
  p3 <- ggplot(g_df, aes(x = 2, y = val, fill = class)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0)+
    scale_fill_brewer(palette=3) +
    theme_void() + xlim(0.5, 2.5) + ggtitle("Cell Cycle")
  
  png(filename = paste0(fig_dir,"2_cellcycle_graph.png"), 
      width = 3500, height = 2000,units = "px", bg = "white", res = 300)
  grid.arrange(arrangeGrob(p_tmp[[1]],p_tmp[[2]],ncol=2,nrow=1),
               arrangeGrob(p1,p2,p3,ncol=3,nrow=1),nrow=2)
  dev.off()
}
############################################################################
#Doublet plot
############################################################################

in_rds_dir = paste0(rds_dir,"5_bcmvn.rds")
in_rds_dir1 = paste0(rds_dir,"5_1_non_filter_Doublet.rds")

Doublet_plot<-function(in_rds_dir,in_rds_dir1){
  bcmvn = readRDS(in_rds_dir)
  pk_val = as.numeric(as.character(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),2]))
  df2 = data.frame(pK = as.numeric(as.character(bcmvn$pK)),BCmetric= bcmvn$BCmetric)
  p0<-ggplot(df2, aes(x=pK, y=BCmetric, group = 1)) + geom_line(color = '#97A8C8') +geom_point(color = '#97A8C8')
  p0<- p0+ theme_classic() + geom_vline(xintercept=pk_val, linetype = 'dotted', color='red', size = 1)
  p0<- p0 + annotate("text", x= pk_val, y=max(bcmvn$BCmetric)+100, 
                     label=paste0("pK : ",pk_val), family="serif", fontface="italic", colour="red", size=3)
  df3 = data.frame(pK = bcmvn$pK,BCmetric= bcmvn$BCmetric)
  p1<-ggplot(data=df3, aes(x=pK, y=BCmetric)) + geom_bar(stat="identity") 
  p1<- p1 + geom_bar(stat="identity", fill="steelblue")+ theme_classic() 
  p1<-  p1 + theme(axis.text.x=element_text(size=rel(0.35)))
  d.data.fit = readRDS(in_rds_dir1)
  
  DF.name = colnames(d.data.fit@meta.data)[grepl("df_class", colnames(d.data.fit@meta.data))]
  
  p2 = DimPlot(d.data.fit, group.by = "orig.ident",cols = c("#A7C9B9","#009BB1","#E6927A"))+
    NoAxes() +  labs(title = "Sample Types")+ theme(plot.title = element_text(size=10))
  p3 = DimPlot(d.data.fit, group.by = DF.name,cols = c("#F7D3BB","#D6E6E3")) + NoAxes() +  
    labs(title = "Doublet composition")+ theme(plot.title = element_text(size=10))
  p4 = VlnPlot(d.data.fit, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)+
    theme(plot.title = element_text(size=10),axis.text.x=element_text(size=rel(0.8)),
          axis.text.y=element_text(size=rel(0.8)))
  df = c(length(which(d.data.fit@meta.data[, DF.name] == "Singlet")),
         length(which(d.data.fit@meta.data[, DF.name] == "Doublet")))
  df = data.frame(Total_Cell  = ncol(d.data.fit),Total_genes=nrow(d.data.fit),
                  Single=df[1],Doublet=df[2])
  tab<- ggtexttable(df, rows = NULL, theme = ttheme("light"))
  tab<- tab %>%
    tab_add_title(text = paste0("Analysis tools : DoubletFinder, pN = 0.25, ","pK = ",pk_val ),
                  face = "plain", size = 10) %>%
    tab_add_title(text = "Statistics of Doublet Filtering", face = "bold", padding = unit(0.1, "line")) %>%
    tab_add_footnote(text = "*AKB-AI", size = 10, face = "italic")
  png(filename = paste0(fig_dir,"3_doublet.png"), 
      width = 3500, height = 2500,units = "px", bg = "white", res = 300)
  grid.arrange(arrangeGrob(p0,p1,p2,p3,ncol=4,nrow=1),
               arrangeGrob(p4,tab,ncol=2,nrow=1),
               nrow=2,heights = c(1.5,2))
  dev.off()
}
############################################################################
#normalization and correlation plot 
############################################################################
normal_plot <- function(in_rds_dir){
  dat_nor = readRDS(in_rds_dir)
  before_nor = data.frame(value=colSums(as.data.frame(GetAssayData(object = dat_nor, slot = "counts"))))
  after_nor = data.frame(value=colSums(as.data.frame(GetAssayData(object = dat_nor, slot = "data"))))
  scale_dat = data.frame(value=colSums(as.data.frame(GetAssayData(object = dat_nor, slot = "scale"))))
  p<-ggplot(before_nor, aes(x=value)) + geom_histogram(color="darkblue", fill="lightblue",bins=100)+
    theme_classic()+  labs(title = "Before Normalization")
  p1<-ggplot(after_nor, aes(x=value)) + geom_histogram(color="darkblue", fill="lightblue",bins=100)+
    theme_classic()+  labs(title = "After Normalization")
  p2<-ggplot(scale_dat, aes(x=value)) + geom_histogram(color="darkblue", fill="lightblue",bins=100)+
    theme_classic()+  labs(title = "scaled data")
  
  png(filename = paste0(fig_dir,"4_normalization.png"), width = 3500, height = 1000,units = "px", bg = "white", res = 300)
  grid.arrange(arrangeGrob(p,p1,p2,ncol=3,nrow=1))
  dev.off()
}


cor_plot <-function(in_rds_dir){
  harmony.obj = readRDS(in_rds_dir)
  harmony.obj@meta.data <-cbind(harmony.obj@meta.data,label_extract(rownames(harmony.obj@meta.data)))
  colnames(harmony.obj@meta.data)[which(names(harmony.obj@meta.data) == "label_extract(rownames(harmony.obj@meta.data))")] <- "sample.name"
  Idents(object = harmony.obj) <- "sample.name"
  av.exp = AverageExpression(harmony.obj,group.by = "sample.name",slot = "counts")
  av.exp1 = apply(av.exp$RNA, 2, function(x) as.numeric(x))
  res <- cor(av.exp1)
  colfunc <- colorRampPalette(c("white","navy"))
  png(filename = paste0(fig_dir,"5_sample_correlation.png"), width = 3500, height = 1500,units = "px", bg = "white", res = 300)
  corrplot(res, type = "upper",t1.col = "white",is.corr = FALSE,diag = TRUE,
           col = colorRampPalette(c("white", "deepskyblue", "blue4"))(100),addCoef.col = "grey",tl.col="black")
  dev.off()  
}

############################################################################
#PCA and UMAP
############################################################################

in_rds_dir = paste0(rds_dir,"6_after_harmony.rds")

pca_unap_plot<-function(in_rds_dir){
  harmony.obj = readRDS(in_rds_dir)
  p0 = JackStrawPlot(harmony.obj, dims = 1:20) + guides(guide_legend(ncol=2,nrow=2,byrow=TRUE))
  
  p = DimPlot(harmony.obj, reduction = "umap",group.by ="seurat_clusters" ,label = TRUE,cols = single_cell_col) + NoLegend() +
    ggtitle("Total Cell types")
  
  harmony.obj@meta.data <-cbind(harmony.obj@meta.data,label_extract(rownames(harmony.obj@meta.data)))
  colnames(harmony.obj@meta.data)[which(names(harmony.obj@meta.data) == "label_extract(rownames(harmony.obj@meta.data))")] <- "sample.name"
  
  p1 = DimPlot(harmony.obj, reduction = "umap",group.by ="seurat_clusters" ,label = TRUE,cols = single_cell_col,split.by = "sample.name") + 
    NoLegend() +ggtitle("per samples")
  
  p2 = DimPlot(harmony.obj, reduction = "umap",split.by = "orig.ident",group.by ="seurat_clusters" ,label = TRUE,
               cols = single_cell_col) + NoLegend() + ggtitle("per sample group")
  pt = ggarrange(ggarrange(p0,p,ncol = 2), 
                 ggarrange(p1,ncol=1), 
                 ggarrange(p2), 
                 nrow = 3,heights = c(1.5,1,1))  
  png(filename = paste0(fig_dir,"8_umap.png"), width = 3500, height = 4200,units = "px", bg = "white", res = 300)
  pt
  dev.off()
  
}
############################################################################
#cell annotation
############################################################################

cell_an_plot<-function(in_rds_dir){
  harmony.obj = readRDS(in_rds_dir)
  p1 = DimPlot(object = harmony.obj, reduction = "umap", label = T,label.size = 4,cols = single_cell_col,group.by ="gcell_type" )+
    guides(guide_legend(ncol=3,nrow=3,byrow=TRUE)) +
    theme(legend.title = element_text( size=6), legend.text=element_text(size=6))+ 
    ggtitle("")
  cell_names = rownames(harmony.obj@meta.data[which(harmony.obj@meta.data$gcell_type=="Oligodendrocyte "),])
  p2 = DimPlot(object = harmony.obj, cells.highlight = cell_names, order = TRUE, split.by = "orig.ident",reduction = "umap",group.by  = "gcell_type")+
    scale_color_manual(values = c("grey", "#5880B4"),labels = c("Other","Oligodendrocyte"))+ ggtitle("")
  cell_marker_lst = as.data.frame(read.table( paste0(q_dir,'1_raw_data/marker_gene_lst_1.txt' ),sep='\t'))
  olg_lst = cell_marker_lst[which(cell_marker_lst[,2] == "Oligodendrocyte "),1]
  p3 = FeaturePlot(harmony.obj, features = olg_lst,ncol =3,cols =c("#F1ECDA","#305388"))
  pt = ggarrange(ggarrange(p1,ncol = 1), 
                 ggarrange(p2,ncol=1), 
                 ggarrange(p3), 
                 nrow = 3)  
  png(filename = paste0(fig_dir,"9_cell_type_annotation_umap.png"), width = 3500, height = 4200,units = "px", bg = "white", res = 300)
  pt
  dev.off()  
}


single_cell_run<-function(in_dat_vec,in_dat_lab_vec){
  mer_in_dat <- merge(sharm_01, y = in_dat_vec, add.cell.ids = in_dat_lab_vec, project = "hyperintensities")
  saveRDS(mer_in_dat, paste0(rds_dir,"0_merge_dat.rds"))  
  
  qc_fun( paste0(rds_dir,"0_merge_dat.rds"))   #QC
  filter_fn(paste0(rds_dir,"1_miQC.rds"))   #filter
  cell_cy_fn(paste0(rds_dir,"2_after_QC.rds"))   #Cell cycle
  pca_an_fn(paste0(rds_dir,"3_after_cell_cycle.rds"))   #PCA analysis
  Doublet_an_fn(paste0(rds_dir,"4_after_pca.rds"))   # Doublet Proportion Estimate
  harmony_an_fn(paste0(rds_dir,"5_after_Doublet.rds"))   #harmony batch correction
  cell_an_fn(paste0(rds_dir,"6_after_harmony.rds"))    # 8. cell annotation
  
  #plot
  stat_plot(paste0(rds_dir,"0_merge_dat.rds"))    # stat(#cell plot)
  qc_plot(paste0(rds_dir,"1_miQC.rds"))    # QC plot
  cc_plot(paste0(rds_dir,"3_after_cell_cycle.rds"))     #Cell cycle plot
  Doublet_plot(paste0(rds_dir,"5_bcmvn.rds"),paste0(rds_dir,"5_1_non_filter_Doublet.rds"))   #Doublet plot
  normal_plot(paste0(rds_dir,"4_after_pca.rds"))    #normalization plot
  cor_plot(paste0(rds_dir,"6_after_harmony.rds"))  #correlation plot 
  pca_unap_plot(paste0(rds_dir,"6_after_harmony.rds"))     #PCA and UMAP
  cell_an_plot(paste0(rds_dir,"7_cell_type_gene_obj.rds"))       #cell annotation
  
}
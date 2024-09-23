convertgenes <- function(genes, inputspecies) {
  ref <- read.csv('./data/human_mouse_ortholog.csv')
  if (inputspecies == 'mouse') {
    ref_filtered <- ref %>%
      filter(Common.Organism.Name == "mouse, laboratory" & Symbol %in% genes)
    class_keys <- ref_filtered$DB.Class.Key
    output <- ref %>%
      filter(DB.Class.Key %in% class_keys & Common.Organism.Name == "human") %>%
      pull(Symbol)
    not_found <- setdiff(genes, ref_filtered$Symbol)
  } else if (inputspecies == 'human') {
    ref_filtered <- ref %>%
      filter(Common.Organism.Name == "human" & Symbol %in% genes)
    class_keys <- ref_filtered$DB.Class.Key
    output <- ref %>%
      filter(DB.Class.Key %in% class_keys & Common.Organism.Name == "mouse, laboratory") %>%
      pull(Symbol)
    not_found <- setdiff(genes, ref_filtered$Symbol)
  } else {
    stop("Invalid input species. Please specify either 'mouse' or 'human'.")
  }
  if (length(not_found) > 0) {
    warning(paste("Ortholog not found for the following genes:", paste(not_found, collapse = ", ")))
  }
  return(output)
}


preprocess_10x_seurat_create<-function(filepath,pipeline='starsolo'){
  if (pipeline=='starsolo'){
    #matrix.dir=paste0(filepath,"/Gene/filtered/")
    matrix.dir<-filepath
    barcode.path <- paste0(matrix.dir,"/barcodes.tsv")
    features.path <- paste0(matrix.dir,"/features.tsv")
    matrix.path <- paste0(matrix.dir, "/matrix.mtx")
    STARmatrix <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(STARmatrix) = barcode.names$V1
    rownames(STARmatrix) = feature.names$V2
    STARmatrix<-as.matrix(STARmatrix)
    seurat<- CreateSeuratObject(STARmatrix)
  }
  else {
    matrix.dir<-paste0(filepath,"/outs/filtered_feature_bc_matrix")
    counts <- Read10X(data.dir = matrix.dir)
    seurat<- CreateSeuratObject(counts)
  }
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-|^MT-")
  seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl]|RP[SL]")
  return(seurat)
}

readmtx<-function(path){
  matrix_path <- path
  mat <- readMM(paste0(matrix_path, "matrix.mtx.gz"))
  # Read in the features and barcodes
  features <- read.delim(gzfile(paste0(matrix_path, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(paste0(matrix_path, "barcodes.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  table_of_counts<-mat
  rownames(table_of_counts) <- features$V2  # Assuming gene IDs are in the first column of features.tsv
  colnames(table_of_counts) <- barcodes$V1
  return(table_of_counts)
}

doublet_detect<-function(seurat,sct=TRUE) {
  #detect doublet
  sweep.res.list <- paramSweep_v3(seurat, PCs = 1:10, sct = sct)
  for(i in 1:length(sweep.res.list)){
    if(length(sweep.res.list[[i]]$pANN[is.nan(sweep.res.list[[i]]$pANN)]) != 0){
      if(i != 1){
        sweep.res.list[[i]] <- sweep.res.list[[i - 1]]
      }else{
        sweep.res.list[[i]] <- sweep.res.list[[i + 1]]
      }
    }
  }
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_v <- as.numeric(as.character(bcmvn$pK))
  pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
  nExp_poi <- round(0.1*length(colnames(seurat)))
  seurat <- doubletFinder_v3(seurat, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  colnames(seurat@meta.data)[ncol(seurat@meta.data)]="DoubletFinder"
  return(seurat)
} 

mad_outlier <- function(sobj, metric, nmads){
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}

preprocess_10x_seurat_filter_low_high<-function(seurat,min.feature=500,min.cell=3,perc_mt=20,perc_ribo=5,max.feature=8000){
  selected_c<-WhichCells(seurat,expression=nFeature_RNA >min.feature & percent.mt< perc_mt & percent.ribo>perc_ribo & nFeature_RNA < max.feature)
  selected_f<-rownames(seurat)[Matrix::rowSums(seurat)>min.cell]
  seurat<-subset(seurat,features=selected_f,cells=selected_c)
  return(seurat)
}


preprocess_10x_seurat_filter_genes<-function(seurat){
  #seurat<-seurat[!grepl('MALAT1',rownames(seurat),ignore.case = TRUE),]
  #seurat<-seurat[!grepl('^MT-',rownames(seurat),ignore.case = TRUE),]
  #seurat<-seurat[!rownames(seurat) %in% gene_XY,]
  #seurat <- seurat[!grepl('^(RPS|RPL)', rownames(seurat), ignore.case = TRUE), ]
}


seurat_process <- function(seurat) {
  # Check if 'percent.mt' and 'percent.ribo' metadata exist
  if(!"percent.mt" %in% colnames(seurat@meta.data)) {
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-|^MT-")
  }
  if(!"percent.ribo" %in% colnames(seurat@meta.data)) {
    seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl]|RP[SL]")
  }
  seurat <- preprocess_10x_seurat_filter_genes(seurat)
  seurat <-preprocess_10x_seurat_filter_low_high(seurat = seurat, min.feature = 400, min.cell = 50, perc_mt = 5, perc_ribo = -1)
  seurat$log1p_total_counts <- log1p(seurat@meta.data$nCount_RNA)
  seurat$log1p_n_genes_by_counts <- log1p(seurat@meta.data$nFeature_RNA)
  bool_vector <- !mad_outlier(seurat, 'log1p_total_counts', 3) & !mad_outlier(seurat, 'log1p_n_genes_by_counts', 3) & !mad_outlier(seurat, 'percent.mt', 3)
  seurat <- subset(seurat, cells = which(bool_vector))
  seurat <- preprocess_10x_seurat_filter_genes(seurat)
  seurat <- NormalizeData(seurat, verbose = FALSE)
  seurat <- FindVariableFeatures(seurat, nfeatures = 3000, selection.method = 'vst') 
  seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  seurat <- ScaleData(seurat, vars.to.regress = c("percent.mt", 'nCount_RNA','S.Score','G2M.Score'))
  seurat <- RunPCA(seurat, verbose = FALSE, npcs = 50)
  seurat <- RunUMAP(seurat, dims = 1:50, verbose = FALSE, n.neighbors = 30, metric = "cosine")
  seurat <- FindNeighbors(seurat, dims = 1:50, verbose = FALSE, k.param = 30)
  seurat <- FindClusters(seurat, verbose = FALSE, resolution = 1.0)
  return(seurat)
}


seurat_norm_dimreduc<-function(seurat,reduction='pca'){
  #seurat <- NormalizeData(seurat, verbose = FALSE)
  seurat<-FindVariableFeatures(seurat,nfeatures=3000,selection.method = 'vst') 
  seurat<- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F,layer='data')
  seurat<-ScaleData(seurat,vars.to.regress = c("percent.mt",'nCount_RNA','S.Score','G2M.Score'))
  seurat <- RunPCA(seurat, verbose = FALSE,npcs=50,)
  seurat <- RunUMAP(seurat, dims =1: 50, verbose = FALSE,n.neighbors = 50,metric = "cosine")
  seurat<- FindNeighbors(seurat, dims = 1:50, verbose = FALSE,k.param = 50)
  seurat <- FindClusters(seurat, verbose = FALSE,resolution = 1)
}

make_soup<-function(path,seurat){
  table_of_droplets<-readmtx(path)
  # Remove duplicate rows based on row names, keeping the first occurrence
  table_of_droplets <- table_of_droplets[!duplicated(rownames(table_of_droplets)), ]
  all_droplets <- colnames(table_of_droplets)
  genes<- rownames(table_of_droplets) %in% rownames(seurat)
  true_cells <- all_droplets %in% colnames(seurat)
  table_of_droplets <- table_of_droplets[genes,]
  table_of_counts <- table_of_droplets[,true_cells]
  sc = SoupChannel(tod=table_of_droplets,toc=table_of_counts)
  cluster_labels<- seurat$seurat_clusters
  names(cluster_labels)<-colnames(seurat)
  sc = setClusters(sc,cluster_labels)
  sc = autoEstCont(sc)
  out = adjustCounts(sc)
}


label_transfer<-function(seurat,ref){
  DefaultAssay(seurat)<-'RNA'
  g.anchors <- FindTransferAnchors(reference = ref, query = seurat )
  predictions <- TransferData(anchorset = g.anchors, refdata = ref$celltype_main2 )
  star_combined <- AddMetaData(object = seurat, metadata = predictions)
}


seurat_plots <- function(seurat, features, colors = NULL, plot = c('vln', 'dot', 'heatmap')) {
  # Ensure plot is one of the allowed types
  plot <- match.arg(plot)
  if (plot == 'vln') {
    return(
      MySeuratWrappers::VlnPlot(seurat, features = features,  
                                stacked = TRUE,
                                pt.size = 0,
                                direction = "vertical",
                                x.lab = NULL, y.lab = NULL) +
        theme(axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 4))
    )
  }
  if (plot == 'dot') {
    return(
      DotPlot(seurat, features = rev(features)) +
        coord_flip() +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 4)) +
        scale_color_gradientn(values = seq(0, 1, 0.2),
                              colours = c('#330066', '#336699', '#66CC66', '#FFCC33')) +
        labs(x = NULL, y = NULL) +
        guides(size = guide_legend(order = 3))
    )
  }
  if (plot == 'heatmap') {
    markerdata <- ScaleData(seurat, features = features, assay = 'RNA')
    return(
      DoHeatmap(markerdata,
                features = features,
                group.by = "celltype",
                assay = 'RNA',
                group.colors = colors,
                size = 3, angle = 90) +
        scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) +
        theme(legend.position = 'left')
    )
  }
}
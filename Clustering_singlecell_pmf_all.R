## Single cell transcripotme analysis of all hematopoietic progenitors

### import data single transcriptome with their metadata 
data<-read.table("data.txt",h=T)
meta<-read.csv("allmeta.csv",h=T,sep="\t",row.names=1)

### create single cell seurat object

library(Seurat)
npm <- CreateSeuratObject(counts = data, project = "npm", min.cells = 3, min.features = 200)
npm <- AddMetaData(npm, metadata = meta)
head(npm[[]])


### evaluate ERCC amount by cells
npm[["percent.ERCC"]] <- PercentageFeatureSet(npm, pattern = "^ERCC-")
VlnPlot(npm, features = c("nFeature_RNA", "percent.ERCC"), ncol = 2)


### preprocessed data
npm <- NormalizeData(npm)
npm <- FindVariableFeatures(npm, selection.method = "vst", nfeatures = 5000)
npm<- ScaleData(npm)

### run PCA for dimension reduction 
npm <- RunPCA(npm, features = VariableFeatures(object = npm),npcs = 50)
ElbowPlot(npm,ndims = 50)


### Run UMAP and TSNE
npm <- RunUMAP(npm, reduction = "pca", dims = 1:20)
DimPlot(npm, reduction = "umap",group.by = "patients",pt.size=1.5)
npm <- RunTSNE(npm, reduction = "pca", dims = 1:20)
library(pals)
DimPlot(npm, reduction = "tsne",group.by = "patients",pt.size=1.5,cols=cols25())
DimPlot(npm, reduction = "tsne",group.by = "patients",pt.size=1.5, cols=cols25(10))
DimPlot(sub, reduction = "tsne",group.by = "patients",pt.size=2.5,cols=cols25(8))



### set identities
Idents(object = npm) <- 'mutations'
Idents(object = npm)

## find clusters
npm <- FindNeighbors(npm, dims = 1:30)
npm <- FindClusters(npm, resolution = 1,algorithm=2)

### subseting
wt<-WhichCells(npm, idents = "WT" , slot ="data")
length(wt)
wt_jak2<-WhichCells(npm, idents = "JAK2_WT" , slot ="data")
length(wt_jak2)
mix<-c(wt,wt_jak2)
length(mix)
sub <- SubsetData(object = npm, cells = mix)

### export subsetting data and metadata as tables
write.table(sub@assays[["RNA"]]@counts, file='sub_count.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(sub@active.ident, file='sub_pheno.tsv', quote=FALSE, sep='\t', col.names = TRUE)

### save Seurat object
save(npm,file="npm.rda")

# work on subset object

DimPlot(sub, reduction = "tsne",group.by = "patients",pt.size=2.5,cols=cols25(8))
DimPlot(sub, reduction = "umap",group.by = "patients",pt.size=2.5,cols=cols25(8))
DimPlot(sub, reduction = "umap",group.by = "diagnosis",pt.size=2.5,cols=cols25(8))
Idents(object = sub) <- 'diagnosis'
Idents(object = sub)
markers <- FindMarkers(sub, ident.1 = "myelofibrosis", ident.2 = "healthy control")
write.table(markers,file="markersmfwt.tsv")
VlnPlot(sub, features = c("GIMAP7"), slot = "data", log = TRUE,split.by= "orig.ident",pt.size=1)
DoHeatmap(object=sub,features=row.names(head(markers,n=200)))
FeaturePlot(sub, features = c("VIM"),min.cutoff = "q9",cols=c("#CCFFFF","darkred"),split.by= "diagnosis",reduction="tsne",pt.size=1.5)


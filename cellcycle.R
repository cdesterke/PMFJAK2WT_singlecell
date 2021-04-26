
### load library
library(Seurat)

### load data
data<-read.table("pmfwt_count.tsv",h=T)
meta<-read.table("phenotypeID1ID2.txt",h=T)
dim(meta)
dim(data)

### create Seurat single cell object
ccobject <- CreateSeuratObject(counts = data, project = "npm", min.cells = 3, min.features = 200)
ccobject[["percent.ERCC"]] <- PercentageFeatureSet(ccobject, pattern = "^ERCC-")
VlnPlot(ccobject, features = c("nFeature_RNA", "percent.ERCC"), ncol = 2)
ccobject <- AddMetaData(ccobject, metadata = meta)


### load cell cycle molecular components

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


### preprocess single cell object
ccobject <- NormalizeData(ccobject)
ccobject<- ScaleData(ccobject)
ccobject <- FindVariableFeatures(ccobject, selection.method = "vst")
ccobject <- ScaleData(ccobject, features = rownames(ccobject))

### regress single cell object on cell cycle molecular components
ccobject <- CellCycleScoring(ccobject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(ccobject[[]])

### barplot of cell numbers cell groups
library(ggplot2)
p<-ggplot(data=res,mapping = aes(x=factor(CellType),fill=CellType)) +
  geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=-20))+
  theme_classic(base_size = 14)+
  geom_bar(aes(label=..count..),stat='count')+coord_flip()
p


### barplot and summary exploring cell cycle phases 
res<-ccobject[[]]

p=ggplot(res,aes(x=factor(CellType),fill=factor(Phase)))+
  geom_bar(position="fill")+
  geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=8,colour="white")+
  scale_fill_brewer(palette="Set2")+
  labs(fill = "Cells")+  
  ggtitle("Cell Cycle WT PMF") +
  xlab("Groups of cells") + ylab("relative proportions ")+theme_classic(base_size = 16)+coord_flip()


p   



table(res$CellType,res$Phase)
prop.table(table(res$CellType,res$Phase))


### cell cycle PCA analysis
ccobject <- RunPCA(ccobject, features = c(s.genes, g2m.genes))
DimPlot(ccobject)

ccobject <- ScaleData(ccobject, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ccobject))
ccobject <- RunPCA(ccobject, features = c(s.genes, g2m.genes))
DimPlot(ccobject,pt.size = 2)


### change metadata with cell cycle phases information 
res$conca <- paste(res$Phase, "-", res$CellType)
ccobject<<-AddMetaData(object=ccobject,metadata=res$conca,col.name='conca')

head(ccobject[[]])

### barplot cell cycle phases versus groups of cells
p<-ggplot(data=res,mapping = aes(x=factor(conca),fill=CellType)) +
  geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=-20))+
  theme_classic(base_size = 14)+
  geom_bar(aes(label=..count..),stat='count')+coord_flip()
p


### Ridegplot cell cycle phases versus groups of cells for 3 markers
RidgePlot(ccobject, features = c("CCND2","CDK6","ID1"), ncol = 3, group.by="conca")

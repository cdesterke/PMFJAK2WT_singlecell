
### load library

library(monocle)


### load data 

data<-read.table("pmfwt_count.tsv",h=T)
meta<-read.table("pmfwt_pheno.tsv",h=T,sep="\t")
dim(meta)
dim(data)
mat<-as.matrix(data)

### build monocle single cell object
pd <- new("AnnotatedDataFrame", data = meta)

genes<-as.data.frame(row.names(data))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)

cds <- newCellDataSet(mat, phenoData = pd,featureData = fd,expressionFamily=tobit())



### definition of single cell hierarchy
cth <- newCellTypeHierarchy()

ID1_id <- row.names(subset(fData(cds), gene_short_name == "ID1"))
ID2_id <- row.names(subset(fData(cds), gene_short_name == "ID2"))

cth <- addCellType(cth, "ID1pos_ID2pos", classify_func = function(x) { x[ID1_id,] > 0 & x[ID2_id,] > 0})
cth <- addCellType(cth, "ID1neg_ID2pos", classify_func = function(x) { x[ID1_id,] <= 0 & x[ID2_id,] > 0})
cth <- addCellType(cth, "ID1neg_ID2neg", classify_func = function(x) { x[ID1_id,] <= 0 & x[ID2_id,] <= 0})
cth <- addCellType(cth, "ID1pos_ID2neg", classify_func = function(x) { x[ID1_id,] > 0 & x[ID2_id,] <= 0})
cds <- classifyCells(cds, cth, 0.1)

table(pData(cds)$CellType)



my_feat <- fData(cds)
my_feat$id<-my_feat$gene_short_name
head(my_feat)


### reduce dimmensions
plot_pc_variance_explained(cds, return_all = FALSE)
cds <- reduceDimension(cds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = TRUE)

### cluster cells
cds <- clusterCells(cds)


#ggplot graph cell type
library(ggplot2)

pie <- ggplot(pData(cds),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())



### distinct feature plots
plot_cell_clusters(cds, 1, 2, color = "cd34") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "CellType",cell_size=1 )


### minimum of expression
cds <- detectGenes(cds, min_expr = 0)
print(head(fData(cds)))
print(head(pData(cds)))
expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))

### search DEG on cell hierarchy (long process time) 
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")


### filtering significant genes
sig_gene_names <- row.names(subset(diff_test_res, pval < 0.05))

sig_gene_names  <- sig_gene_names[!grepl('^ERCC-', sig_gene_names)]

diff_test_res$id<-row.names(diff_test_res)


dim(diff_test_res) 
sig_gene_names  

str(sig_gene_names)

my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:832]


### subseting monocle object with significant genes
cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)

### dimensional reduction DDRTree algorithm and and order cells on the pseudotime
cds2 <- reduceDimension(cds2, method = 'DDRTree')

cds2 <- orderCells(cds2)

### build graphs on trajectory 
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="ID1",cell_size=2 )+  facet_wrap(~CellType)
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="ID2",cell_size=2 )+  facet_wrap(~diseasesub)

plot_cell_trajectory(cds2, color_by = "Pseudotime",cell_size=1.5)+ theme(legend.position = "right")+scale_color_viridis_c()
plot_cell_trajectory(cds2, color_by = "cd123",cell_size=1)

plot_cell_trajectory(cds2, color_by = "diseasesub",markers="ID1",markers_linear = F,show_branch_points=T)+
   theme(legend.position = "right")
plot_cell_trajectory(cds2, color_by = "Pseudotime",markers="EGR1",markers_linear = F) +  facet_wrap(~group)


plot_cell_trajectory(cds2, color_by = "blasts",markers="ID2",markers_linear = F,show_branch_points=T)+
  geom_point(alpha=0.1)+ theme(legend.position = "right")+ scale_color_manual(breaks = waiver(),values=c("lightblue","darkred"))

plot_cell_trajectory(cds2, color_by = "State",markers="SRSF9",markers_linear = T)+ theme(legend.position = "right")

plot_cell_trajectory(cds2, color_by = "CellType",cell_size=2 ) +  facet_wrap(~group)
plot_cell_trajectory(cds2, color_by = "CellType",cell_size=2 ) +  facet_wrap(~group)


### gene_to_cluster pseudotime heatmap
gene_to_cluster <-my_ordering_genes[1:100]
gene_to_cluster  <- gene_to_cluster[!grepl('^ERCC-', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^ZNF790-AS1.1', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^LOC', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^HLA-', gene_to_cluster)]
conca<-c("ID1","ID2","JAK2",gene_to_cluster)



my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[conca,],cores = 8,show_rownames = TRUE,return_heatmap = TRUE,cluster_rows = TRUE)


### pseudotime expression plot
plot_genes_in_pseudotime(cds2[c("ID1","LMO2","FOS","PMAIP1","ZFP36","UNK","TSPYL2"),],cell_size = 1, color_by = "cd45ra",ncol = 1)

plot_genes_in_pseudotime(cds2[c("ID2","ID3","HES1","RELA","MEIS1","JUNB"),],cell_size = 2, color_by = "State",ncol = 1)
plot_genes_in_pseudotime(cds2[c("ID1","JUN","NFKBIA","SOCS2","BEX1","FOS","ZFP36","LMO2"),],cell_size = 2, color_by = "State",ncol = 1)
plot_genes_in_pseudotime(cds2[c("ID1","ID2","EID2","JAK2","HSPB1","CAT"),],cell_size = 2, color_by = "State",ncol = 1)


### save the data
write.table(pData(cds2),file="phenotypeID1ID2.txt")
write.table(diff_test_res,file="difMoncoletrajectoryID1ID2.txt")
save(cds,file="monoclecdsID1ID2.rda")
save(cds2,file="monoclecdssuperID1ID2.rda")


### THE END ###








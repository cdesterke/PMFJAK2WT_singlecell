library(monocle)
data<-read.table("pmfwt_count.tsv",h=T)
meta<-read.table("pmfwt_pheno.tsv",h=T,sep="\t")
dim(meta)
dim(data)
mat<-as.matrix(data)

pd <- new("AnnotatedDataFrame", data = meta)

genes<-as.data.frame(row.names(data))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)

cds <- newCellDataSet(mat, phenoData = pd,featureData = fd,expressionFamily=tobit())



cth <- newCellTypeHierarchy()

ID1_id <- row.names(subset(fData(cds), gene_short_name == "ID1"))
ID2_id <- row.names(subset(fData(cds), gene_short_name == "ID2"))

cth <- addCellType(cth, "ID1pos_ID2pos", classify_func = function(x) { x[ID1_id,] > 0 & x[ID2_id,] > 0})
cth <- addCellType(cth, "ID1neg_ID2pos", classify_func = function(x) { x[ID1_id,] <= 0 & x[ID2_id,] > 0})
cth <- addCellType(cth, "ID1neg_ID2neg", classify_func = function(x) { x[ID1_id,] <= 0 & x[ID2_id,] <= 0})
cth <- addCellType(cth, "ID1pos_ID2neg", classify_func = function(x) { x[ID1_id,] > 0 & x[ID2_id,] <= 0})
cds <- classifyCells(cds, cth, 0.1)

table(pData(cds)$CellType)

#ID1neg_ID2neg ID1neg_ID2pos ID1pos_ID2neg ID1pos_ID2pos 
#149            89            66            84



my_feat <- fData(cds)
my_feat$id<-my_feat$gene_short_name
head(my_feat)

plot_pc_variance_explained(cds, return_all = FALSE)
cds <- reduceDimension(cds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = TRUE)


cds <- clusterCells(cds)

table(pData(cds)$CellType)
head(pData(cds))


#ggplot graph cell type
library(ggplot2)

pie <- ggplot(pData(cds),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


plot_cell_clusters(cds, 1, 2, color = "cd34") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "CellType",cell_size=1 )


cds <- detectGenes(cds, min_expr = 0)
print(head(fData(cds)))
print(head(pData(cds)))


expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")

sig_gene_names <- row.names(subset(diff_test_res, pval < 0.05))

sig_gene_names  <- sig_gene_names[!grepl('^ERCC-', sig_gene_names)]

diff_test_res$id<-row.names(diff_test_res)


dim(diff_test_res) 
sig_gene_names  

str(sig_gene_names)

my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:832]
cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)

#dimensional reduction and gene selection
cds2 <- reduceDimension(cds2, method = 'DDRTree')

cds2 <- orderCells(cds2)

#build graphs on trajectory 
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


#gene_to_cluster <- row.names(diff_test_res)[order(diff_test_res$qval)][1:75] 
gene_to_cluster <-my_ordering_genes[1:100]
gene_to_cluster  <- gene_to_cluster[!grepl('^ERCC-', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^ZNF790-AS1.1', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^LOC', gene_to_cluster)]
gene_to_cluster  <- gene_to_cluster[!grepl('^HLA-', gene_to_cluster)]
conca<-c("ID1","ID2","JAK2",gene_to_cluster)



my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[conca,],cores = 8,show_rownames = TRUE,return_heatmap = TRUE,cluster_rows = TRUE)
plot_genes_in_pseudotime(cds2[c("ID1","LMO2","FOS","PMAIP1","ZFP36","UNK","TSPYL2"),],cell_size = 1, color_by = "cd45ra",ncol = 1)

plot_genes_in_pseudotime(cds2[c("ID2","ID3","HES1","RELA","MEIS1","JUNB"),],cell_size = 2, color_by = "State",ncol = 1)
plot_genes_in_pseudotime(cds2[c("ID1","JUN","NFKBIA","SOCS2","BEX1","FOS","ZFP36","LMO2"),],cell_size = 2, color_by = "State",ncol = 1)
plot_genes_in_pseudotime(cds2[c("ID1","ID2","EID2","JAK2","HSPB1","CAT"),],cell_size = 2, color_by = "State",ncol = 1)
#save data
write.table(pData(cds2),file="phenotypeID1ID2.txt")
write.table(diff_test_res,file="difMoncoletrajectoryID1ID2.txt")
save(cds,file="monoclecdsID1ID2.rda")
save(cds2,file="monoclecdssuperID1ID2.rda")
###


## beam analysis B2
BEAM_res <- BEAM(cds2, branch_point = 2, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

sig_gene_names2 <- row.names(subset(BEAM_res, pval < 0.05))
my_ordering_genes <- row.names(BEAM_res)[order(BEAM_res$pval)][1:2735]
gene_to_cluster2 <- row.names(BEAM_res)[order(BEAM_res$qval)][1:90] 
conca2<-c("ID1","ID2","JAK2",gene_to_cluster2)


plot_genes_branched_heatmap(cds2[conca2,],branch_point = 2, num_clusters = 3,
                            cores = 1, use_gene_short_name = T, show_rownames = T)


genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("ID2", "ID3", "HES1","MEIS1","RELA","JUNB")))
plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 3,
                               color_by = "State",
                               ncol = 1)





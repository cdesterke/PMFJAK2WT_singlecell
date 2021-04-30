# R script to perform down sampling with equilbrate group of cells

## load libraries
library(Seurat)
library(dplyr)


## work on count and metadata of Seurat object nammed sub

### extract matrix of count
mat_sub<-sub@assays[["RNA"]]@counts
mat_sub<-as.matrix(mat_sub)

### extract identities
ident_sub<-sub@active.ident
ident_sub<-as.data.frame(ident_sub)
tmat_sub<-t(mat_sub)
tmat_sub<-as.data.frame(tmat_sub)

### merge count and identities
all<-merge(ident_sub,tmat_sub,by="row.names")


### select equilibrate groups
library(dplyr)
df <- all %>% group_by(ident_sub) %>% sample_n(80)
head(df[1:5,1:5])
table(df$ident_sub)
tdf<-t(df)
write.table(tdf,file="GSEA8080.tsv")

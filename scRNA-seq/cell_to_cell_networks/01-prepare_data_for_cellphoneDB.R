################################################
###### FOR Github TONSIL atlas
################################################
OBJ<-readRDS(“path/to/SeuratObjectWithSelectedCellTypes.rds”)
#create INPUT data
MD<-data.frame(cell=rownames(OBJ@meta.data),celltype=as.character(OBJ$CellType))
write.table(MD, file=“path/to/CPDB/MDcelllabels.txt”,sep=“\t”,quote=F, row.names = F)
writeMM(OBJ@assays$RNA@data, file =“path/to/CPDB/input/data/matrix.mtx”)
write(x = rownames(OBJ@assays$RNA@data), file = “path/to/CPDB/input/data/features.tsv”)
write(x = colnames(OBJ@assays$RNA@data), file = “path/to/CPDB/input/data/barcodes.tsv”)
Idents(OBJ)<-OBJ$CellType
DEGtb <- FindAllMarkers(OBJ,  verbose = T, only.pos = F, random.seed = 1, logfc.threshold = 0, min.pct = 0.1,
                        return.thresh = 1,slot = “data”)
write.table(DEGtb,file=“path/to/CPDB/input/data/TableDEG_allvsRest_Wilcox.txt”,quote=F,sep=“\t”,row.names = F)
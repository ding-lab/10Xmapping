##########################################################################################
### Author: Lijun Yao
### Date: April 16, 2019
### Usage: sc_mutation_mapping.R -p <project_name> -s <sample_id> -r <raw_path> -o <output_directory> -b <object_path> -e <seurat_dir> -v <verbose> -h <help>
##########################################################################################

##########################################################################################
### Load libraries and set working directory
##########################################################################################

#system("source activate r_env")

library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tibble)
library(getopt)

spec = matrix(c(
      'verbose', 'v', 2, "integer",
      'help'   , 'h', 0, "logical",
      'project'  , 'p', 1, "character",
      'sample'   , 's', 1, "character",
      'raw_path' , 'r', 1, "character",
      'output_path', 'o', 1, "character",
      'object_path','b',1,"character",
      'seurat_dir','e',1,"character"
      ), byrow=TRUE, ncol=4)

opts <- getopt(spec)

# Print help and usage

if ( !is.null(opts$help) ) { cat(getopt(spec, usage=TRUE)); q(status=1) }

# Assign defaults if value not provided

if ( is.null(opts$verbose ) ) { opts$verbose = FALSE }

# Check if all required arguments were provided

if ( is.null(opts$project) ) {print("Project name is a required argument"); q(status=1)}
if ( is.null(opts$sample) ) {print("Sample name is a required argument"); q(status=1)}
if ( is.null(opts$raw_path) ) {print("Raw matrix path is a required argument"); q(status=1)}
if ( is.null(opts$output_path) ) {print("Output directory is a required argument"); q(status=1)}

# Assign to variables

project_name <- opts$project
sample <- opts$sample
input_dir <- opts$raw_path #dir contains 10xmapping heatmap .txt
out_dir <- opts$output_path 
obj_path <- opts$object_path #.rds path
seurat_dir <- opts$seurat_dir



mapping_files <- Sys.glob(file.path(input_dir, "*_heatmap_*.txt"))
#############################################################################
###Cell barcodes with var/ref of mutations respectively
#############################################################################
ct_ref_bc <- list() #directory: cell type as keys, cell bc as values
ct_var_bc <- list()
for (file in mapping_files) {
  ct <- strsplit(strsplit(strsplit(file,"/")[[1]][12],"_")[[1]][5],".t")[[1]][1]
  ref_bc <- list() #directory: mutation as keys, cell bc as values
  var_bc <- list()
  heatmap <- read.table(file,sep=",",header=TRUE)
  rownames(heatmap) <- heatmap[,1]
  mutations <- intersect(gsub('-Var', '', heatmap[,1]),gsub('-Ref', '', heatmap[,1]))
  for (mut in mutations) {
    mut_ref <- paste0(mut,"-Ref")
    mut_var <- paste0(mut,"-Var")
    mut_rows <- heatmap[c(mut_ref, mut_var),]
    mut_rows_NA_rm <- mut_rows[ ,colSums(is.na(mut_rows)) == 0]
    mut_rows_NA_rm[,1] <- NULL
    mut_ref_vaf_above0 <- as.data.frame(mut_rows_NA_rm[,mut_rows_NA_rm[1,]>0])
    mut_var_vaf_above0 <- as.data.frame(mut_rows_NA_rm[,mut_rows_NA_rm[2,]>0])
    ref_bc[[mut]]=colnames(mut_ref_vaf_above0)
    var_bc[[mut]]=colnames(mut_var_vaf_above0)
  }
  ct_ref_bc[[ct]] <- unlist(ref_bc,use.names = FALSE)# unlist a dataframe in R without column names
  ct_var_bc[[ct]] <- unlist(var_bc,use.names = FALSE)# unlist a dataframe in R without column names
}

#unlist dir and save the cell barcodes for REF ALLELES
all_ct_ref_bc <- unlist(ct_ref_bc,use.names = FALSE)
all_ct_ref_bc <- all_ct_ref_bc[!duplicated(all_ct_ref_bc)]

#unlist dir and save the cell barcodes for VAR ALLELES
all_ct_var_bc <- unlist(ct_var_bc,use.names = FALSE)
all_ct_var_bc <- all_ct_var_bc[!duplicated(all_ct_var_bc)]

bc_with_both <- intersect(all_ct_ref_bc,all_ct_var_bc)
bc_with_both <- gsub("[.]","-", bc_with_both)
all_ct_ref_bc_unique <- all_ct_ref_bc[!(all_ct_ref_bc %in% bc_with_both)]
all_ct_ref_bc_unique <- gsub("[.]","-", all_ct_ref_bc_unique)
all_ct_var_bc_unique <- all_ct_var_bc[!(all_ct_var_bc %in% bc_with_both)]
all_ct_var_bc_unique <- gsub("[.]","-", all_ct_var_bc_unique)

#################################################################################################
###save the coordinate file
#################################################################################################

Sobj <- readRDS(obj_path)
tsne_coord <- Embeddings(object = Sobj[["umap"]]) %>% as.data.frame
write.table(tsne_coord,paste0(seurat_dir,"/",sample,"_umap_coord.txt"),sep="\t",row.names = TRUE,col.names = TRUE,quote = FALSE)

###########################################################################
###read the coordinate file
########################################################################
tsne_coord <- read.table(paste0(seurat_dir,"/",sample,"_umap_coord.txt"),header = TRUE,stringsAsFactors = FALSE)
tsne_coord_tmp <- tsne_coord %>% rownames_to_column("BC")
tsne_coord_tmp$MUT <- "NA"
tsne_coord_tmp$MUT[tsne_coord_tmp$BC %in% all_ct_ref_bc_unique] <- "Reference"
tsne_coord_tmp$MUT[tsne_coord_tmp$BC %in% all_ct_var_bc_unique] <- "Variant"
tsne_coord_tmp$MUT[tsne_coord_tmp$BC %in% bc_with_both] <- "Variant"
rownames(tsne_coord_tmp) <- tsne_coord_tmp$BC

#####################################################################
###mapping mutations to the tsne
###################################################################
#pdf(paste0(out_dir,"/coding_mutations_mapping_to_",sample,".pdf"))
pdf(paste0(out_dir,"/coding_mutations_var_ref_mapping_to_",sample,".pdf"))
p<-ggplot()
p<-p+geom_point(data=tsne_coord_tmp,aes(x=UMAP_1,y=UMAP_2),color="grey",alpha=0.3,shape=16)
cols <- c("Variant" = "#EF3B2C", "Reference" = "#4292C6", "Both" = "black")
tsne_coord_tmp_rm_NA <- tsne_coord_tmp[tsne_coord_tmp$MUT!="NA",]
p<-p+geom_point(data=tsne_coord_tmp_rm_NA,aes(x=UMAP_1,y=UMAP_2,color=MUT),alpha=0.75)+scale_color_manual(values =cols)
print(p)
dev.off()


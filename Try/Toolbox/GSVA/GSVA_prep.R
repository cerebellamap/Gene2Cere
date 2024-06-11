## Code to prepare for GSVA

## Required packages
library(ggpubr)
library(stringr)
library(GSVA)
library(GSEABase)
library(ComplexHeatmap)
library(limma)
library(msigdbr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)
library(cowplot)
library(patchwork)
library(cluster)
# ---------------------------




## GMT
output_dir = "/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/GSVA/GSVA_mgt/"
setwd(output_dir)
# ---------------------------
# GO
GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>%
  dplyr::select(gene_symbol,gs_exact_source,gs_name,gs_subcat)
GO_df = GO_df[GO_df$gs_subcat!="HPO",]
GO_df$gs_name_short <- 0
for (i in 1:length(GO_df$gs_name)){
  print(i)
  GO_df$gs_name_short[i] = str_to_title(gsub('[_]', ' ',str_split_fixed(GO_df$gs_name[i], '_',2)[2]))
}
save(GO_df, file = "GO_df.Rda")


# KEGG
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>%
  dplyr::select(gs_exact_source,gene_symbol,gs_name)
KEGG_df$gs_name_short <- 0
for (i in 1:length(KEGG_df$gs_name)){
  KEGG_df$gs_name_short[i] = str_to_title(gsub('[_]', ' ',str_split_fixed(KEGG_df$gs_name[i], '_',2)[2]))
}
save(KEGG_df, file = "KEGG_df.Rda")


# Cell Type
# 11 cluster has download the from web
library(biomaRt)
library(httr)
set_config(use_proxy(url="http://127.0.0.1:1212", port=211)) # 需要 global
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
CT_path = '/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/GSVA/GSVA_mgt/Cerebellum_culster_11/'
file_list = list.files(path = CT_path, pattern = "cluster")
for (file in file_list){
  print(file)
  ms_ct <- read.csv(paste0(CT_path,file),header = T)$gene
  m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                  values = ms_ct, mart = mouse,
                  attributesL = c("hgnc_symbol","entrezgene_id","entrezgene_description"), #,"chromosome_name","start_position"
                  martL = human,uniqueRows = T)
  write.csv(x=m2h.g, file = paste0(str_split_fixed(file, '.csv',2)[1],'_2hgnc.csv'))
}
# Construct gmt
f=("CellType_Cerebellum_11.gmt")
file_list = list.files(path = output_dir, pattern = "cluster")
sink(f)
for (file in file_list){
  # print(file)
  geneset <- data.frame(unique(read.csv(paste0(output_dir,file),header = T)$HGNC.symbol))
  gs_name = str_split_fixed(str_split_fixed(file, '_',2)[2], '.csv',2)[1]
  names(geneset) = str_split_fixed(gs_name,'_2hgnc',2)[1]
  lapply(names(geneset),function(i){
  cat(paste(c(i,'tmp',geneset[[i]]),collapse='\t'))
  cat('\n') 
  })
}
sink()


# Disease
library(devtools)
install_bitbucket("ibi_group/disgenet2r")
gmturl = file.path("http://www.disgenet.org", "static/disgenet_ap1/files/downloads/gmt_files", "disgenet.curated.v7.symbols.gmt")
download.file(url = gmturl, destfile = "DisGeNET.gmt")
## ---------------------------


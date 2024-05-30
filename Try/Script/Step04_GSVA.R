## Code to run GSVA
# source("./Toolbox/GSVA/GSVA_prep.R")  # Need to run in first time
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


## Read the normalized expression matrix
# ---------------------------
data_dir <- '/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/'
gmt_dir = "/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/GSVA/GSVA_mgt/"
output_dir = '/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/GSVA'
setwd(output_dir)
gene_FG1=read.csv(paste0(data_dir,'Step03_GCIsig.csv'), row.names = 1, header=TRUE, check.names=FALSE)
gene_FG1_sig = gene_FG1[gene_FG1[['p_perm_fdr0.05']]=='True',] # 0.05 1024
expression = read.csv(paste0(data_dir,'Step01_Gene_expression.csv'), row.names = 1, header=TRUE, check.names=FALSE)
report = read.csv(paste0(data_dir,'Step01_Sample_info.csv'), row.names = 1, header=TRUE, check.names=FALSE)
a = split(report, 20)
expression_sig = expression[, colnames(expression) %in% rownames(gene_FG1_sig)]
report$Net=0
for (x in rownames(report)){
  if (report[x,'Y']<0.5){
    report[x,'Net'] = 'Motor'
  }
  else{
    report[x,'Net'] = 'Nonmotor'
  }
}
group_list <- report$Net

# color
RdBu = read.csv('/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/RdBu_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
PiYG = read.csv('/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/PiYG_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
YlGn = read.csv('/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/YlGn_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
RdYlGn = read.csv('/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/RdYlGn_Hex_200.csv', row.names = 1, header=TRUE, check.names=FALSE)
set.seed(123)
# ---------------------------


## GMT
# ---------------------------
load(file = paste0(gmt_dir, "GO_df.Rda"))
GO_list = split(GO_df$gene_symbol,GO_df$gs_name_short)

load(file =  paste0(gmt_dir,"KEGG_df.Rda"))
KEGG_list = split(KEGG_df$gene_symbol,KEGG_df$gs_name_short)

CT_list <- getGmt( paste0(gmt_dir,"CellType_Cerebellum_11.gmt"))

DGN_list <- getGmt( paste0(gmt_dir,"DisGeNET.gmt"))
for (x in 1:length(DGN_list)){
  names(DGN_list) <- str_split(DGN_list@.Data[[x]]@shortDescription,'DISGENET_')[[1]][2]
  DGN_list@.Data[[x]]@setName <- str_split(DGN_list@.Data[[x]]@shortDescription,'DISGENET_')[[1]][2]
}
# ---------------------------


## bin-GSVA
# ---------------------------
rep_sort = report[order(report$Y),]
exp_sort = expression_sig[order(report$Y),]

Exp_sp <- NULL
Rep_sp <- NULL
GSVA <- NULL
DEG <- NULL
for(i in seq(20,80,20)){
  print(i)
  # data
  Exp_sp[[i]] = sapply(split(exp_sort,factor(sort(rank(row.names(exp_sort))%%i))),
                FUN=function(x) colMeans(x))
  colnames(Exp_sp[[i]]) = seq(1,i,1)
  Rep_sp[[i]] = sapply(split(rep_sort$Y,factor(sort(rank(row.names(rep_sort))%%i))),
                FUN=function(x) mean(x))
  Rep_sp[[i]]=as.data.frame(Rep_sp[[i]])
  colnames(Rep_sp[[i]]) = 'Y'
  rownames(Rep_sp[[i]]) = seq(1,i,1)
  Rep_sp[[i]]$Net=0
  for (x in rownames(Rep_sp[[i]])){
    if (Rep_sp[[i]][x,'Y']<0.5){Rep_sp[[i]][x,'Net'] = 'Motor'}
    else{Rep_sp[[i]][x,'Net'] = 'Nonmotor'}
  }
  group_list <- Rep_sp[[i]]$Net
  save(Rep_sp, file = 'Rep_sp.Rda')

  # GSVA
  for (j in c('GO', 'DGN')) { # 'CT', 'KEGG' # j = 'DGN'
    print(j)
    GSVA <- gsva(expr=as.matrix(Exp_sp[[i]]), gset.idx.list=get(paste0(j,'_list')),
                        mx.diff=FALSE, verbose=FALSE, min.sz=2, parallel.sz=32)
    write.csv(x=GSVA, file=paste0('GSVA_',j,'_', i,'bin','.csv'),row.names = T)
    save(GSVA, file = paste0('GSVA_',j,'_',i,'bin','.Rda'))
    # limma
    design <- model.matrix(~0+factor(group_list))
    colnames(design) <- levels(factor(group_list))
    rownames(design) <- colnames(GSVA)
    contrast.matrix<-makeContrasts("Nonmotor-Motor",levels = design)
    fit <- lmFit(GSVA,design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    DEG <- topTable(fit2, number=Inf)
    write.csv(x=DEG, file=paste0('GSVA_',j,'_limma_',i,'bin','.csv'),row.names = T)
    save(DEG, file = paste0('GSVA_',j,'_limma_',i,'bin','.Rda'))
  }

  # Visulization
  for (j in c('GO','DGN')) { # , 'KEGG', 'CT',
    load(file = 'Rep_sp.Rda')
    load(file = paste0('GSVA_',j,'_',i,'bin','.Rda'))
    load(file = paste0('GSVA_',j,'_limma_',i,'bin','.Rda'))
    # plot
    df = data.frame(label = row.names(DEG),x = DEG$logFC, y = -log10(DEG$adj.P.Val))
    row.names(df) = df$label
    df$Group = 'no'
    for (x in df$label){
      if (df[x, 'y']> -log10(0.05) && df[x,'x']>0){
          df[x,'Group'] = 'Association'
          }
      else if (df[x, 'y']> -log10(0.05) && df[x,'x']<0){
          df[x,'Group'] = 'Sensorimotor'
      }
      else{df[x,'Group'] = 'no'}
    }

    top5l <- filter(df, Group != "no"& x<0) %>%
              distinct(label, .keep_all = T) %>%
              top_n(5, abs(x))
    top5r <- filter(df, Group != "no"& x>0) %>%
              distinct(label, .keep_all = T) %>%
              top_n(5, abs(x))
    limmaplot = ggplot(data = df, aes(x = x, y =y)) +
        geom_point(data = filter(df, x<0 & y>-log10(0.05)), aes(x = x, y = y, color = x),
        size = 2.5, show.legend = T,alpha=1) +
        scale_color_gradientn(colours = RdBu[0:99,1], name = 'log2 (count + 1)') +
        new_scale_color() +
        geom_point(data = filter(df, x>0 & y>-log10(0.05)), aes(x = x, y = y, color = -x),
        size = 2.5, show.legend = T,alpha=1) +
        scale_color_gradientn(colours = rev(RdBu[100:199,1])) + # rev(RdBu[100:199,1])
        new_scale_color() +
        geom_point(data = filter(df,Group == "no"), aes(x = x, y = y), color ="#bbbbbb", size = 2.5, show.legend = T,alpha=0.7) +
        labs(x = base::expression(paste(log[2], 'FC')), y = base::expression(paste('-',log[10], 'FDR')))+
        geom_hline(yintercept = -log10(0.05), linetype = 'dotdash', color = 'grey20') +
        geom_vline(xintercept = 0, linetype = 'dotdash', color = 'grey20') +
        geom_text_repel(data = top5l, aes(x = x, y = y, label = label),
            force = 10, color = 'black', size = 2.5, point.padding = 0.5, hjust = 0.1,
            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
            segment.color="black", segment.size = 0.2, nudge_x = 0.01, nudge_y = 0.1)+
        geom_text_repel(data = top5r, aes(x = x, y = y, label = label),
            force = 10, color = 'black', size = 2.5, point.padding = 0.5, hjust = 0.1,
            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
            segment.color="black", segment.size = 0.2, nudge_x = -0.01, nudge_y = 0.1)+# 
        # ggtitle('Geneset differential analysis') +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "grey20"),
              axis.title.x = element_text(size = 11,  face ="plain"),
              axis.title.y = element_text(size = 11,  face ="plain"),
              axis.text.x  = element_text(size = 9,  face ="plain"),
              axis.text.y  = element_text(size = 9,  face ="plain"),
              plot.margin = margin(t = 25, r = 0.02, b = 5, l = 0.0002, unit = "pt"),
              plot.title = element_text(size = 15,  hjust = 0.5, vjust = 2),
              axis.line.x  = element_line(size = 0.8),
              axis.line.y  = element_line(size = 0.8),
              legend.background = element_blank(), legend.key = element_blank(),
              legend.position = 'none', legend.title = element_blank(),
              # plot.margin=unit(rep(0.1,4),'lines')
        )


    if (j=='GO'){
      DEG_sig <- DEG[DEG$adj.P.Val<0.005 & abs(DEG$logFC) > 0,]
    }else {
      DEG_sig <- DEG[DEG$adj.P.Val<0.01 & abs(DEG$logFC) > 0,]
    }
    if(dim(DEG_sig)[1]!=0){
      dat <- GSVA[match(rownames(DEG_sig),rownames(GSVA)),]
      
      annt_top = HeatmapAnnotation(FG1 = Rep_sp[[i]]$Y,
                value = anno_points(Rep_sp[[i]]$Y, size=unit(1, "mm"),ylim=c(0,1),
                                    axis_param = list(at = c(0,0.5,1))), #
                annotation_name_gp=gpar(fontsize = c(7,7)),
                                  annotation_name_side =c('left','left'),
                col = list(FG1 = circlize::colorRamp2(c(0,0.5,1),hcl_palette = 'RdBu', reverse = T)),
                annotation_height = c(0.7, 3.5),show_legend = c(F, F),
                annotation_label = c("", ""),border = c(F, T))
      km <- kmeans(t(dat), 2)
      km.order <- as.numeric(names(sort(km$centers[,1])))
      names(km.order) <- toupper(letters)[1:2]
      km.order <- sort(km.order)
      # clus.order <- factor(names(km.order[km$cluster]))
      split =  factor(names(km.order[km$cluster]), levels=c("A","B"))

      
      heatmap = Heatmap(dat, name='Exp', 
                width = ncol(dat)*unit(7, "mm"), height = nrow(dat)*unit(7, "mm"),
                cluster_rows = T, row_dend_reorder = F, row_title = "Gene set", 
                row_title_gp = gpar(fontsize = 12), show_row_dend = F, 
                row_labels = str_wrap(rownames(dat), width=30), row_names_gp = gpar(fontsize = 9),
                column_title = "Percentile along FG", column_title_side = "bottom", 
                column_title_gp = gpar(fontsize = 12),
                column_split=split, cluster_column_slices =  T, column_names_gp = gpar(fontsize = 9),#cluster_columns = dend,# column_split = rep(c("Non-motor", "Motor"), 10), #column_dend_height = 4,# col = circlize::colorRamp2(c(-1,0,1), hcl_palette = 'YlGn', reverse = T),
                col = colorRampPalette(c("#BCBCBC","#FFFFFF","#5188ba"))(50),
                top_annotation = annt_top, show_heatmap_legend = F) 
      heatmap = draw(heatmap) # ,column_title="Geneset expression",column_title_gp=grid::gpar(fontsize=15)



      png(file = paste0('GSVA_',j,'_',i,'bin','.png'), width = 4000, height = 3500, units = 'px', res = 300)
      grob = grid.grabExpr(draw(heatmap)) 
      patchwork = limmaplot +  grob + # + plot_spacer()
        plot_layout(ncol=1, guides = "collect",
                    width=c( 0.35,0.7), # heights = c(0.3,0.7),# 
                    design = 'AB')
        # + plot_layout(widths = c(4, -1.1 ,4.5),guides = "collect")& theme(legend.position = "top")
      patchwork = patchwork # + plot_annotation( title = 'Gene ontology', theme = theme(plot.title = element_text(size = 20,hjust = 0.01, vjust = 0)))
      print(patchwork) # heights = c(1.7), width=c( 0.1,2)
      dev.off()
    }else {
      print(paste0(j, ' has no sig term'))
    }

  }
}
# ---------------------------



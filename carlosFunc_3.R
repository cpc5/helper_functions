

gg_geneLabel4 <- function(df, 
gene_col='gene_name', 
dotsize=2,
dotcolor='black',
text_lbl_size=4,
colorText='black',
plotGenes=c() ) {

lfc <- grep('log2FoldChange|logFC|avg_log2FC', colnames(df))
p.col <- grep('padj|FDR|p_val_adj', colnames(df))

df$logp <-  -log10( df[[p.col]] )
df$lfc <-  df[[lfc]]
df$gene_col <-  df[[gene_col]]

require(ggplot2)
require(ggrastr)
require(ggrepel)

genes <- subset(df, gene_col %in% plotGenes ) 

list(
geom_point_rast(data=genes,
aes(x=lfc, y= logp),
size=dotsize,color=dotcolor),
geom_label_repel(fill = "white",data=genes,
aes(x=lfc,y= logp,
label=genes$gene_col),
size=text_lbl_size,color=colorText)
)
}

############################################

gg_volcano3 <- function(df,
    abs_lfc= 0 , 
    p.cutoff = 0.05,
    up.color="steelblue3",
    down.color="steelblue3",
    unSig.color='grey88', 
    dotSize=0.5, 
    alpha=1                   
    ) {

require(dplyr)
require(ggplot2)
require(ggrastr)
require(ggrepel)


lfc <- grep('log2FoldChange|logFC|lfc_|_lfc', colnames(df),ignore.case=T )
p.col <- grep('padj|FDR', colnames(df))

df$logp <-  -log10( df[[p.col]] )
df$lfc <-  df[[lfc]]

maxlfc <- max(abs(df$lfc), na.rm = TRUE)

lfc_breaks <- seq(-maxlfc, maxlfc, length.out=8) %>% round_any(accuracy=2, x=.)

maxpadj <- range(df$logp , finite=T , na.rm=T)[2]
padj_breaks <- seq(0, maxpadj,  length.out=5) %>% round_any(accuracy=10, x=.)

Sig <- subset(df, 
df[p.col] < p.cutoff  & 
abs(lfc) > abs_lfc)

unSig <- subset(df ,  
!rownames(df) %in% rownames(Sig)  )

u <- subset(Sig, lfc > 0 )
d <- subset(Sig, lfc < 0 )


ggplot()+
geom_point_rast(data=unSig,
    aes(x=lfc,y= logp), size=dotSize, colour=unSig.color) + 
geom_point_rast(data= u,
    aes(x=lfc ,y= logp ),col = up.color, alpha = alpha, size=dotSize, inherit.aes = T) +
geom_point_rast(data=d,
    aes(x=lfc ,y= logp ),col = down.color, alpha = alpha, size=dotSize, inherit.aes = T) +
xlab("log2 Fold Change")+
ylab(" -log10(FDR)")+
scale_y_continuous(trans='asinh',
    limits=c(0, maxpadj),
    breaks=c(padj_breaks))+
scale_x_continuous(trans='asinh',
    limits=c(-maxlfc,maxlfc),
    breaks=c(lfc_breaks))

}


gg_geneLabel3 <- function(df, 
gene_col='gene_id', 
dotsize=4,
text_lbl_size=4,
colorText='black',
dotcolor='black',
plotGenes=c() ) {

lfc <- grep('log2FoldChange|logFC|lfc_|_lfc', colnames(df),ignore.case=T )
p.col <- grep('padj|FDR', colnames(df))
df$logp <-  -log10( df[[p.col]] )
df$lfc <-  df[[lfc]]
df$gene_col <-  df[[gene_col]]

require(ggplot2)
require(ggrastr)
require(ggrepel)

genes <- subset(df, gene_col %in% plotGenes ) 

list(
geom_point_rast(data=genes,
      aes(x=lfc, y= logp),
      size=dotsize,color=dotcolor)
#geom_text_repel(data=genes,
#      aes(x=lfc,y= logp,
#          label=genes$gene_col),
#      size=text_lbl_size,color=colorText)
)
}

#################################
#################################

writeTable2 <- function(x,y) {

write.table(x,
file=paste0(y),
sep='\t',
quote=F,
row.names=F, 
col.names=T
)
}

#################################
#################################

library(scales)
asinh_trans<-function(){
trans_new(name='asinh', transform=function(x) asinh(x),
inverse=function(x) sinh(x))}

#################################
#################################

gg_pca <- function(df, textsize=3) {
require(ggrepel)

percentVar <- round(100 * attr(df, "percentVar"))

ggplot(df,
  aes(PC1, PC2,color=group, label=df$name))+
  geom_point(size=3)+
  ylim( min(df$PC2)-3 , max(df$PC2)+3 )+
  xlim( min(df$PC1)-3 , max(df$PC1)+3 )+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", size=textsize)+
  theme(legend.title = element_blank(), legend.position = c(1,1)) 

}


#################################
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
#################################

gg_geneLabel2 <-function(df, 
gene_col='gene_name', 
dotsize=4,
dotcolor='black',
text_lbl_size=4,
colorText='black',
plotGenes=c() ) {

lfc <- grep('log2FoldChange|logFC', colnames(df))
p.col <- grep('padj|FDR', colnames(df))

df$logp <-  -log10( df[[p.col]] )
df$lfc <-  df[[lfc]]
df$gene_col <-  df[[gene_col]]

require(ggplot2)
require(ggrastr)
require(ggrepel)

genes <- subset(df, gene_col %in% plotGenes ) 

list(
geom_point_rast(data=genes,
aes(x=lfc, y= logp),
size=dotsize,color=dotcolor),
geom_text_repel(data=genes,
aes(x=lfc,y= logp,
label=genes$gene_col),
size=text_lbl_size,color=colorText)
)
}

#################################


gg_volcano2 <- function(df,
abs_lfc= 0 , 
p.cutoff = 0.05,
sig.color="steelblue3",
unSig.color='grey88', 
dotSize=1, 
alpha=0.5
                       
) {
require(dplyr)
require(ggplot2)
require(ggrastr)
require(ggrepel)

lfc <- grep('log2FoldChange|logFC|lfc|lfc', colnames(df),ignore.case=T )
p.col <- grep('padj|FDR', colnames(df))

df$logp <-  -log10( df[[p.col]] )
df$lfc <-  df[[lfc]]

maxlfc <- range( df$lfc , finite=T , na.rm=T) %>% max()
lfc_breaks <- seq(-maxlfc, maxlfc, length.out=8) %>% round_any(accuracy=2, x=.)

maxpadj <- range(df$logp , finite=T , na.rm=T)[2]
padj_breaks <- seq(0, maxpadj,  length.out=5) %>% round_any(accuracy=10, x=.)

Sig <- subset(df, 
df[p.col] < p.cutoff  & 
abs(lfc) > abs_lfc)

unSig <- subset(df ,  
!rownames(df) %in% rownames(Sig)  )

ggplot()+
geom_point_rast(data=unSig,
aes(x=lfc,y= logp),size=dotSize, colour=unSig.color) + 
geom_point_rast(data =Sig, 
aes(x=lfc ,
y= logp ),
col = sig.color,
alpha = alpha, 
size=dotSize, 
inherit.aes = T) +
xlab("log2 Fold Change") + 
ylab(" -log10(FDR)")+
scale_y_continuous(trans='asinh',
limits=c(0, maxpadj),
breaks=c(padj_breaks))+
scale_x_continuous(trans='asinh',
limits=c(-maxlfc,maxlfc),
breaks=c(lfc_breaks))

}

#################################
txi_fromtximport_salmon <- function(salmonFolder, 
metaTable) {

require(DESeq2)
require(tximport)
require(dplyr)
require(stringr)
require(readr)

star_salmon_path = salmonFolder
tx2gene = read_tsv( paste0(star_salmon_path, "/salmon_tx2gene.tsv"))
metad = metaTable
rownames(metad) <- metad$Sample
files <- file.path(star_salmon_path, metad$Sample, "quant.sf")
toremove <- c(star_salmon_path,"/", "quant.sf")
names(files) <- str_remove_all(files, paste(toremove, collapse = "|"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
return(txi)
}

##################################################################

Deseq2_fromMatrix <- function(df,
    condition1,
    condition2,
    condition1reps,
    condition2reps,
    lfcshrinkType= 'normal') {

res <- data.frame()

require(DESeq2, quietly = TRUE)

arg1=condition1
arg2=condition2

condition1reps=condition1reps
condition2reps=condition2reps

conditionLevels <- c(arg1, 
arg2)

sampleCon <- c(
rep(arg1, condition1reps), 
rep(arg2, condition2reps)
)

ExpDesign <- data.frame(row.names=colnames(df), 
condition=sampleCon)

print(ExpDesign)

ddsMtx <- DESeqDataSetFromMatrix(countData = df, 
colData=ExpDesign, 
design=~condition)

dds <- DESeq(ddsMtx)
rldHT <- vst(dds, blind=FALSE)

res <- results(dds,
alpha = 0.5,
contrast = c("condition",
condition2,
condition1) 
) 

res <- lfcShrink(
dds=dds,
res=res,
coef = 2,
type = lfcshrinkType)

res <- data.frame(res) %>% 
rownames_to_column("gene_id") %>%
.[complete.cases(.$padj),]

normCounts <- counts(dds, normalized=T) %>%
data.frame() %>% 
rownames_to_column("gene_id")

res <- merge(res, normCounts ,   
by='gene_id')

pca <- plotPCA(rldHT, intgroup='condition',  ntop=2000, returnData=TRUE)

return(list(res, pca, dds))


}

#############################################################
#############################################################
importSTARcounts <- function(dir, column="yes") {

require(dplyr)
require(tools)

#V2="no"
#V3="yes"
#V4="reverse"

if(column=='yes'){ 
    column='V3' 
    } else if(column=='reverse') {
    column='V4'
    }

x <- list.files(
    dir, 
    pattern=".tab", 
    recursive=T, 
    full.names=T
    )
names(x) <- list.files( dir )
print(x)

x <- lapply(x, 
            read.delim, header=F)
x <- lapply(x, 
            subset, !V1 %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous") )
x <- lapply(x, 
            dplyr::select, V1, column)

for( i in 1:length(x) ) {
    
colnames(x[[i]]) <- c("gene_id", names(x)[i] )

}
x <- plyr::join_all(x)
x <- head(x, n=-5)
rownames(x) <- x$gene_id
x$gene_id <- NULL

return(x)

}

#############################################################
runGO <- function(test.vector, 
univr_vector, 
ont="BP",
organism='mouse',
simplify=TRUE) {

require(dplyr)
require(ggplot2)
require(clusterProfiler)
require(org.Mm.eg.db)
require(org.Hs.eg.db)

if(organism=='mouse') {
  OrgDb='org.Mm.eg.db'
  } else {
  OrgDb='org.Hs.eg.db'
  }

eg <- read_csv("/project/imoskowitz/shared/annotations/iGenomes//Mus_musculus/mart_cpc.csv")

eg$NCBI <- as.character(eg$NCBI)

univr_vector <- subset(eg, gene_id %in% univr_vector)$NCBI

goPATH <- enrichGO(gene=as.vector(subset(eg, gene_id %in% test.vector)$NCBI), 
OrgDb = OrgDb, 
universe=univr_vector,
ont=ont,
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05)

goPATH <- setReadable(goPATH, OrgDb = OrgDb)         
goPATH <- simplify(goPATH)

goPATH@result <-  goPATH@geneSets %>% 
lapply(. ,length) %>% 
unlist %>% data.frame("GenesInGo" = .) %>%
tibble::rownames_to_column("ID") %>% 
left_join(goPATH@result, .)

return(goPATH)

}

#################################
plotGO <- function(obj, topTerms=5, muhcolor = "grey40", strWrap = 30, barWidth = 0.5, expandY = c(0.1,0)) {

goPATH <- data.frame(obj) %>%
.[order(.$p.adjust, decreasing = F),] %>%
.[c(1:topTerms),] 

goPATH$Desc <- goPATH$Description %>%
stringr::str_wrap(., strWrap)  

goPATH$DEGRatio <- goPATH$Count/goPATH$GenesInGo

ggplot(data=goPATH, 
    aes(x=-log10(p.adjust),y=as.factor(Desc), size = DEGRatio ))+
geom_point(fill=muhcolor, color = muhcolor )+
ylab("")+
xlab("-log10 Adj. Pval")+
scale_size(range = c(0, 20))+
scale_y_discrete(
limits= as.vector(goPATH[order(-log10(goPATH$p.adjust),decreasing=F),])$Desc ,
expand=expandY
)
}



salmon_import <- function(star_salmon_path) {
  pkgs <- c("readr", "stringr", "tximport")#, "testpackage")
  pkg_avail <- sapply(pkgs, \(x) require(x, character.only = TRUE))
  if(!all(pkg_avail)) {
    warning(
      "the following required packages are not installed:",
      paste0("\n - ", pkgs[!pkg_avail])
    )
  }
  ##  Read tx2gene file
  tx2gene <- as.data.frame(readr::read_tsv(
    file.path(star_salmon_path, "tx2gene.tsv"),
    col_names = FALSE
  ))

  ##  Get list of quant.sf files
  files <- list.files(
    path = star_salmon_path,
    pattern = "quant.sf",
    recursive = TRUE,
    full.names = TRUE
    )

  ##  Create names for the files
  toremove <- c(star_salmon_path, "/", "quant.sf")
  names(files) <- stringr::str_remove_all(files, paste(toremove, collapse = "|"))

  ##  Import quantification data
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene)

  return(c(txi, list(tx2gene = tx2gene)))
}

load_metadata <- function(metadata_path) {
  metadata <- readxl::read_excel(metadata_path) %>%
    dplyr::slice(-(1:2)) %>%
    dplyr::slice(1:(n() - 1)) %>%
    data.frame() %>%
    `colnames<-`(gsub(names(.), pattern = '..GEO.sample.field..title..', replacement="")) %>%
    dplyr::select(Experiment.ID, Sample.ID, Sample.name, Age.stage) %>%
    .[complete.cases(.),]

  rownames(metadata) <- paste0(metadata$Experiment.ID, "_", metadata$Sample.ID)

  metadata <- metadata %>%
  mutate(treatment = factor(case_when(grepl("dox|tmx|ko|oe", Sample.name, ignore.case = TRUE) ~ str_extract(tolower(Sample.name), "dox|tmx|ko|oe"),
	TRUE ~ "ctrl"
	 )))
  metadata$condition <- sub("_[^_]+$", "", metadata$Sample.name)
  metadata$time_condition <- paste0(metadata$Age.stage, "_", metadata$treatment)
  return(metadata)
}

get_deseq_results <- function(dds, contrast, treatment, control, alpha = 0.05) {
  contrast_vector <- c(contrast, treatment, control)
  results(dds, contrast = contrast_vector, alpha = alpha) %>%
    lfcShrink(dds = dds,
              res = .,
              contrast = contrast_vector,
              type = 'normal') %>%
    data.frame() %>%
    rownames_to_column("gene_id") %>%
    .[complete.cases(.$padj),]
}


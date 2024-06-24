##this script is to compare SNP and SV heri
##this script is to draw a bar composition plot 

library(reshape2)
library(ggplot2)
library(RColorBrewer)

########
##step01
##s1
##check average h2 of GA
ipt__heri_snp_dt <- read.delim('GA_GD.heritability_snp.txt')
head(ipt__heri_snp_dt)
ipt_heri_snp_dt <- ipt__heri_snp_dt[grepl('.GA',ipt__heri_snp_dt$Type),]
ipt_heri_snp_dt <- ipt__heri_snp_dt[,c('Type','h2')]
colnames(ipt_heri_snp_dt) <- c('Type','h2_snp')

ipt_heri_sv_dt <- read.delim('GA_GD.heritability_sv.txt')
head(ipt_heri_sv_dt)
ipt_heri_sv_dt <- ipt_heri_sv_dt[grepl('.GA',ipt_heri_sv_dt$Type),]
ipt_heri_sv_dt <- ipt_heri_sv_dt[,c('Type','h2')]
colnames(ipt_heri_sv_dt) <- c('Type','h2_sv')

mergedt_dt <- merge(ipt_heri_snp_dt,ipt_heri_sv_dt,by.x = 'Type',by.y = 'Type')
mergedt_reshap_dt <- melt(mergedt_dt)
colnames(mergedt_reshap_dt) <- c('Treat','Cate','Heritability')

##plot the box plot
p <- ggplot(mergedt_reshap_dt, aes(x=Cate, y=Heritability)) + 
  geom_boxplot() + 
  #geom_point(data = points_to_label, color = 'red', size = 3) +
  #geom_text(data = points_to_label, aes(label = rownames(points_to_label)), hjust = -0.2,color='red',cex= 7) + 
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text.x = element_text(colour = "black", size=30,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
    axis.text.x = element_text(size=20,colour = "black"),
    axis.text.y = element_text(size=20,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
#scale_x_continuous(breaks=seq(0, 2, 0.5)) +
#scale_y_continuous(breaks=seq(0, 2, 0.5)) +
#coord_cartesian(xlim = c(0, 2), ylim = c(0,2))
p
pdf('opt_s1_compare_heri_SNP_SV.pdf',width = 8,height = 8)
p
dev.off()


########
##step01 plot which trait the SNPs having higher heri and SVs having higher heri
##s2
##for the GA snp
top_num = 50
ipt_heri_snp_dt <- ipt_heri_snp_dt[order(ipt_heri_snp_dt$h2_snp,decreasing = T),]
ipt_heri_snp_top_dt <- ipt_heri_snp_dt[1:top_num,]
ipt_heri_snp_top_dt$Type <- factor(ipt_heri_snp_top_dt$Type,levels = ipt_heri_snp_top_dt$Type)

p <- ggplot(data=ipt_heri_snp_top_dt, aes(x=Type, y=h2_snp)) +
  geom_bar(stat="identity", width=0.5) + 
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35,),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )

pdf(paste0('./opt_s2_snp_h2_top_barplot.pdf'),width = 20,height = 15)
print(p)
dev.off()

##for the GA SV
ipt_heri_sv_dt <- ipt_heri_sv_dt[order(ipt_heri_sv_dt$h2_sv,decreasing = T),]
ipt_heri_sv_top_dt <- ipt_heri_sv_dt[1:top_num,]
ipt_heri_sv_top_dt$Type <- factor(ipt_heri_sv_top_dt$Type,levels = ipt_heri_sv_top_dt$Type)

p <- ggplot(data=ipt_heri_sv_top_dt, aes(x=Type, y=h2_sv)) +
  geom_bar(stat="identity", width=0.5) + 
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35,),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )

pdf(paste0('./opt_s2_sv_h2_top_barplot.pdf'),width = 20,height = 15)
print(p)
dev.off()

intersect(ipt_heri_snp_top_dt$Type,ipt_heri_sv_top_dt$Type )




########
##step02
input_SNP_top_sample_dt <- read.delim('opt_SNP_top10samplemethodNum.txt',header =F)
head(input_SNP_top_sample_dt)
colnames(input_SNP_top_sample_dt) <- c('TraitSample','Trait','Sample','Number','Mlist')
input_SNP_top_sample_flt_dt <- input_SNP_top_sample_dt[input_SNP_top_sample_dt$Number > 5,]
unique(input_SNP_top_sample_flt_dt$Trait)

input_SV_top_sample_dt <- read.delim('opt_SV_top10samplemethodNum.txt',header = F)
head(input_SV_top_sample_dt)
colnames(input_SV_top_sample_dt) <- c('TraitSample','Trait','Sample','Number','Mlist')
input_SV_top_sample_dt <- input_SV_top_sample_dt[input_SV_top_sample_dt$Number > 5,]
unique(input_SV_top_sample_dt$Trait)

shared_traits <- intersect(unique(input_SNP_top_sample_flt_dt$Trait),unique(input_SV_top_sample_dt$Trait))
all_traits <- unique(c(unique(input_SNP_top_sample_flt_dt$Trait),unique(input_SV_top_sample_dt$Trait)))


##Here we will want to divide the bar plot into three cates
##1) have traits to be shared
##so we care about the shared_traits

outs <- lapply(all_traits,function(x){
  
  
  ipt_snp_target_trait_dt <- input_SNP_top_sample_flt_dt[input_SNP_top_sample_flt_dt$Trait == x,] 
  ipt_sv_target_trait_dt <- input_SV_top_sample_dt[input_SV_top_sample_dt$Trait == x,] 
  
  shared_samples <- intersect(ipt_snp_target_trait_dt$Sample,ipt_sv_target_trait_dt$Sample)
  shared_samples_num <- length(shared_samples)
  
  snp_samples_num <- nrow(ipt_snp_target_trait_dt)
  sv_samples_num <- nrow(ipt_sv_target_trait_dt)
  
  res <- c()
  res$trait <- x
  res$share <- shared_samples_num
  res$snp <- snp_samples_num
  res$sv <- sv_samples_num
  
  return(res)

})

combine_dt <- do.call(rbind,outs)
rownames(combine_dt) <- all_traits
combine_dt <- as.data.frame(combine_dt)
class(combine_dt)
combine_dt$share <- as.numeric(combine_dt$share)
combine_dt$snp <- as.numeric(combine_dt$snp)
combine_dt$sv <- as.numeric(combine_dt$sv)
combine_dt$trait <- as.character(combine_dt$trait)
combine_dt <- combine_dt[order(combine_dt$share,decreasing = T),]
dim(combine_dt)

write.table(combine_dt,'opt_summary_sample_shared_unique_snp_sv.txt',quote = F,sep = '\t')

########
##step02 plot a heatmap to show the number of samples per trait in the top 9
ipt_dt <- read.delim('opt_SV_top10samplemethodNum.txt',header = F)
head(ipt_dt)
ipt_dt <- ipt_dt[ipt_dt$V4 == '9',]
ipt_dt$count <- 1
ipt_dt <- ipt_dt[c('V2','V3','count')]

write.table(ipt_dt,'opt_SV_nineMethod_sparse.txt',sep = '\t',quote = F)
ipt_dt <- read.table('opt_SV_nineMethod_sparse.txt',stringsAsFactors = T)

score_mtx <- sparseMatrix(i=as.numeric(ipt_dt$V2),
                          j=as.numeric(ipt_dt$V3),
                          x=as.numeric(ipt_dt$count),
                          dimnames=list(levels(ipt_dt$V2), levels(ipt_dt$V3)))
score_mtx <- as.matrix(score_mtx)


pdf(paste0('opt_SV_nineMethod_sparse.pdf'),width = 5,height = 8)
pheatmap(score_mtx,
         cluster_cols = F,
         cluster_rows = F,
         #border_color = "black",
         border_color = "white",
         fontsize_col = 12,
         fontsize_row = 12,
         na_col = "white",
         #color=colorRampPalette(rev(c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8',"#91bfdb","#4575b4")))(100),
         color=rev(brewer.pal(n = 8, name = "Set2")),
         fontsize_number = 25
) 
dev.off()


ipt_dt <- read.delim('opt_SNP_top10samplemethodNum.txt',header = F)
head(ipt_dt)
ipt_dt <- ipt_dt[ipt_dt$V4 == '9',]
ipt_dt$count <- 1
ipt_dt <- ipt_dt[c('V2','V3','count')]

write.table(ipt_dt,'opt_SNP_nineMethod_sparse.txt',sep = '\t',quote = F)
ipt_dt <- read.table('opt_SNP_nineMethod_sparse.txt',stringsAsFactors = T)

score_mtx <- sparseMatrix(i=as.numeric(ipt_dt$V2),
                          j=as.numeric(ipt_dt$V3),
                          x=as.numeric(ipt_dt$count),
                          dimnames=list(levels(ipt_dt$V2), levels(ipt_dt$V3)))
score_mtx <- as.matrix(score_mtx)


pdf(paste0('opt_SNP_nineMethod_sparse.pdf'),width =5,height = 5)
pheatmap(score_mtx,
         cluster_cols = F,
         cluster_rows = F,
         #border_color = "black",
         border_color = "white",
         fontsize_col = 12,
         fontsize_row = 12,
         na_col = "white",
         #color=colorRampPalette(rev(c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8',"#91bfdb","#4575b4")))(100),
         color=rev(brewer.pal(n = 8, name = "Set2")),
         fontsize_number = 25
) 
dev.off()



########
##step03 plot the count of samples shared by different samples
ipt_dt_SV <- read.delim('opt_SV_top10samplemethodNum.txt',header = F)
head(ipt_dt_SV)

ipt_sample_number_SV_dt <- as.data.frame(table(ipt_dt_SV$V4))
colnames(ipt_sample_number_SV_dt) <- c('method_num','trait_sample_num')
ipt_sample_number_SV_dt$Cate <- 'SV'

ipt_dt_SNP <- read.delim('opt_SNP_top10samplemethodNum.txt',header = F)
head(ipt_dt_SNP)

ipt_dt_SNP_8 <- ipt_dt_SNP[ipt_dt_SNP$V4 == '8',]
table(ipt_dt_SNP_8$V2)


ipt_sample_number_SNP_dt <- as.data.frame(table(ipt_dt_SNP$V4))
colnames(ipt_sample_number_SNP_dt) <- c('method_num','trait_sample_num')
ipt_sample_number_SNP_dt$Cate <- 'SNP'

ipt_combine_sample_num_dt <- rbind(ipt_sample_number_SV_dt,ipt_sample_number_SNP_dt)


p <- ggplot(data=ipt_combine_sample_num_dt, aes(x=method_num, y=trait_sample_num, fill = Cate)) +
  geom_bar(stat="identity", width=0.5,position=position_dodge()) + 
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35,),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=30,hjust = 0.5 ),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
p

pdf(paste0('./opt_SNP_SV_num_compare.pdf'),width = 25,height = 15)
print(p)
dev.off()



















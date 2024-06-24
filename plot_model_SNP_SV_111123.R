##this script is to compare SNP and SV heri for different machine learning methods
##this script is to draw a bar composition plot 

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(pheatmap)

########
##step01 plot the snp and SV for the different model accuracy
snp_acc_fl <- list.files(path = 'snp_accuracy/')

combine_dt_SNP <- read.delim(paste0('snp_accuracy/',snp_acc_fl[1]),header = F)

for (i in 2:length(snp_acc_fl)){
  
  new_dt <- read.delim(paste0('snp_accuracy/',snp_acc_fl[i]),header = F)
  
  combine_dt_SNP <- rbind(combine_dt_SNP,new_dt)

  
}

dim(combine_dt_SNP)
table(combine_dt_SNP$V2)

median(combine_dt_SNP$V3) ##0.9451241

combine_dt_SNP$cate <- 'SNP'

sv_acc_fl <- list.files(path = 'SV_accuracy//')

combine_dt_SV <- read.delim(paste0('SV_accuracy/',sv_acc_fl[1]),header = F)

for (i in 2:length(sv_acc_fl)){
  
  new_dt <- read.delim(paste0('SV_accuracy/',sv_acc_fl[i]),header = F)
  
  combine_dt_SV <- rbind(combine_dt_SV,new_dt)
  
  
}

combine_dt_SV$cate <- 'SV'

min(combine_dt_SV$V3) 

combine_dt_SNP_SV <- rbind(combine_dt_SNP,combine_dt_SV)
head(combine_dt_SNP_SV)
colnames(combine_dt_SNP_SV) <- c('Trait','Method','Accuracy','cate')

combine_dt_SNP_SV <- combine_dt_SNP_SV[combine_dt_SNP_SV$Accuracy >= 0,]
head(combine_dt_SNP_SV)

p <- ggplot(data=combine_dt_SNP_SV, aes(x=Method, y=Accuracy,fill=cate)) +
  geom_boxplot()+ 
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
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )

pdf(paste0('./opt_compare_acc_different models.pdf'),width = 16,height = 10)
print(p)
dev.off()

##check the wilcxon test
all_methods <- names(table(combine_dt_SNP_SV$Method))
outs <- lapply(all_methods, function(x){
  
  combine_dt_SNP <- combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SNP',]
  combine_dt_SV <-  combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SV',]
  
  combine_dt_SNP_tmethod <- combine_dt_SNP[combine_dt_SNP$Method == x,]
  combine_dt_SV_tmethod <- combine_dt_SV[combine_dt_SV$Method == x,]
  dim(combine_dt_SNP_tmethod)
  dim(combine_dt_SV_tmethod)
  res <- wilcox.test(combine_dt_SNP_tmethod$Acc,combine_dt_SV_tmethod$Acc)
  pval <- res$p.value
    
  opt <- c()
  opt$method <- x
  opt$pval <- pval
  
  return(opt)
  
})

combine_dt <- do.call(rbind,outs)

write.csv(combine_dt,'opt_statistic_test.txt',quote = F)

write.table(combine_dt_SNP_SV, 'opt_all_combine_method_acc.txt',quote = F, sep = '\t')

###########################################
##check the median value of each prediction
##we will first calculate the average and then check the median value
combine_dt_SNP <- combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SNP',]
head(combine_dt_SNP)
colnames(combine_dt_SNP) <- c('Trait','Method','Accuracy','Marker')
aggregate_avg_SNP <- aggregate(Accuracy ~ Trait + Method, data = combine_dt_SNP, mean)
aggregate_avg_median_SNP <- aggregate(Accuracy ~ Method, data = aggregate_avg_SNP, median)
aggregate_avg_median_SNP$Marker <- 'SNP'

aggregate_avg_median_SNP <- aggregate_avg_median_SNP[order(aggregate_avg_median_SNP$Accuracy,decreasing = T),]
target_method_str <- aggregate_avg_median_SNP$Method

combine_dt_SV <- combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SV',]
head(combine_dt_SV)
colnames(combine_dt_SV) <- c('Trait','Method','Accuracy','Marker')
aggregate_avg_SV <- aggregate(Accuracy ~ Trait + Method, data = combine_dt_SV, mean)
aggregate_avg_median_SV <- aggregate(Accuracy ~ Method, data = aggregate_avg_SV, median)
aggregate_avg_median_SV$Marker <- 'SV'

combine_median_SV_SNP <- rbind(aggregate_avg_median_SNP,aggregate_avg_median_SV)

combine_median_SV_SNP$Method <- factor(combine_median_SV_SNP$Method,levels = target_method_str)

p<-ggplot(combine_median_SV_SNP, aes(x=Method, y=Accuracy, group=Marker)) +
  geom_line(aes(color=Marker))+
  geom_point(aes(color=Marker)) +
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35,),
    #axis.text.x = element_blank(),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
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
pdf(paste0('./opt_compare_acc_diff_method.pdf'),width = 8,height = 8)
p
dev.off()








##############################
##check the standard deviation
##for SNP
combine_dt_SNP <- combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SNP',]

all_method_list = unique(combine_dt_SNP$Method)
outs <- lapply(all_method_list,function(x){
  
  combine_dt_SNP_t_method <- combine_dt_SNP[combine_dt_SNP$Method == x,]
  combine_dt_SNP_t_method_sd <- sd(combine_dt_SNP_t_method$Accuracy)
  
  res <- c()
  res$method <- x
  res$sd <- combine_dt_SNP_t_method_sd
  res$cate <- 'SNP'
  
  return(as.data.frame(res))
})

combine_final_SNP_dt <- as.data.frame(do.call(rbind,outs))
combine_final_SNP_order_dt <- combine_final_SNP_dt[order(combine_final_SNP_dt$sd,decreasing = F),]

##for SV
combine_dt_SV <- combine_dt_SNP_SV[combine_dt_SNP_SV$cate == 'SV',]

all_method_list = unique(combine_dt_SV$Method)
outs <- lapply(all_method_list,function(x){
  
  combine_dt_SV_t_method <- combine_dt_SV[combine_dt_SV$Method == x,]
  combine_dt_SV_t_method_sd <- sd(combine_dt_SV_t_method$Accuracy)
  
  res <- c()
  res$method <- x
  res$sd <- combine_dt_SV_t_method_sd
  res$cate <- 'SV'
  return(as.data.frame(res))
})

combine_final_SV_dt <- as.data.frame(do.call(rbind,outs))

combine_final_SNP_SV_dt <- rbind(combine_final_SNP_dt,combine_final_SV_dt)

combine_final_SNP_SV_dt$method <- factor(combine_final_SNP_SV_dt$method,levels = combine_final_SNP_order_dt$method)

p<-ggplot(combine_final_SNP_SV_dt, aes(x=method, y=sd, group=cate)) +
  geom_line(aes(color=cate))+
  geom_point(aes(color=cate)) +
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35,),
    #axis.text.x = element_blank(),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
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
pdf(paste0('./opt_compare_acc_sd_diff_method.pdf'),width = 8,height = 8)
p
dev.off()



##updating 062324
##we will cal the range mean sd for specific traits for all models based on the SV and SNP
ipt_dt <- read.delim('Table S5 model accuracy_adjustnotequalTo1.txt')
head(ipt_dt)

aggregate_trait_mean_dt <- aggregate(Accuracy ~ Cate + Trait, data = ipt_dt, mean)
aggregate_trait_mean_dt$catetrait <- paste0(aggregate_trait_mean_dt$Cate,'__',aggregate_trait_mean_dt$Trait)
rownames(aggregate_trait_mean_dt) <- aggregate_trait_mean_dt$catetrait
aggregate_trait_mean_dt <- aggregate_trait_mean_dt[c('Cate','Trait','Accuracy')]
colnames(aggregate_trait_mean_dt) <- c('Cate','Trait','Mean')

aggregate_trait_sd_dt <- aggregate(Accuracy ~ Cate + Trait, data = ipt_dt, sd)
aggregate_trait_sd_dt$catetrait <- paste0(aggregate_trait_sd_dt$Cate,'__',aggregate_trait_sd_dt$Trait)
rownames(aggregate_trait_sd_dt) <- aggregate_trait_sd_dt$catetrait
aggregate_trait_sd_dt <- aggregate_trait_sd_dt[c('Accuracy')]
colnames(aggregate_trait_sd_dt) <- c('Standard Deviation')

aggregate_trait_min_dt <- aggregate(Accuracy ~ Cate + Trait, data = ipt_dt, min)
aggregate_trait_min_dt$catetrait <- paste0(aggregate_trait_min_dt$Cate,'__',aggregate_trait_min_dt$Trait)
rownames(aggregate_trait_min_dt) <- aggregate_trait_min_dt$catetrait
aggregate_trait_min_dt <- aggregate_trait_min_dt[c('Accuracy')]
colnames(aggregate_trait_min_dt) <- c('Lowest')

aggregate_trait_max_dt <- aggregate(Accuracy ~ Cate + Trait, data = ipt_dt, max)
aggregate_trait_max_dt$catetrait <- paste0(aggregate_trait_max_dt$Cate,'__',aggregate_trait_max_dt$Trait)
rownames(aggregate_trait_max_dt) <- aggregate_trait_max_dt$catetrait
aggregate_trait_max_dt <- aggregate_trait_max_dt[c('Accuracy')]
colnames(aggregate_trait_max_dt) <- c('Highest')

combined_dt <- cbind(aggregate_trait_mean_dt,aggregate_trait_sd_dt,aggregate_trait_min_dt,aggregate_trait_max_dt)
head(combined_dt)

write.table(combined_dt, './opt_final_acc_mean_sd_lowest_highest.txt',quote = FALSE,sep = '\t')





########
##previous
snp_acc_fl <- list.files(path = 'snp_accuracy/')

combine_dt_SNP <- read.delim(paste0('snp_accuracy/',snp_acc_fl[1]),header = F)

for (i in 2:length(snp_acc_fl)){
  
  new_dt <- read.delim(paste0('snp_accuracy/',snp_acc_fl[i]),header = F)
  
  combine_dt_SNP <- rbind(combine_dt_SNP,new_dt)
  
  
}

dim(combine_dt_SNP)
head(combine_dt_SNP)
table(combine_dt_SNP$V1)

##calculate the average
aggregate_avg <- aggregate(V3 ~ V1 + V2, data = combine_dt_SNP, mean)
head(aggregate_avg)
dim(aggregate_avg)

aggregate_sample_median_acc <- aggregate(V3 ~ V1, data = aggregate_avg, median)

aggregate_sample_median_acc <- aggregate_sample_median_acc[order(aggregate_sample_median_acc$V3,decreasing = T),]






write.table(combine_dt_SNP,'opt_threecol_SNP_acc.txt',quote = F,sep = '\t')
ipt_dt <- read.table('opt_threecol_SNP_acc.txt',stringsAsFactors = T)
ipt_mtx <- sparseMatrix(i=as.numeric(ipt_dt$V1),
                         j=as.numeric(ipt_dt$V2),
                         x=as.numeric(ipt_dt$V3),
                         dimnames=list(levels(ipt_dt$V1), levels(ipt_dt$V2)))
ipt_mtx <- as.matrix(ipt_mtx)
dim(ipt_mtx)

pdf(paste0('opt_SNP_sample_acc_heatmap.pdf'),width = 8,height = 15)
pheatmap(ipt_mtx,
         #scale="row",
         #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
         #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
         cluster_cols = F,
         cluster_rows = F,
         border_color='white',
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize_number = 25
) 
dev.off()








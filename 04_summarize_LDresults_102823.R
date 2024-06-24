##updating 102923 we will plot the snp and sv maf distribution
##this script is to draw plot to compare the LD between SNP and TEs
##this script will generate a bar plot for each traits that show the number of TE markers nearby genes
library(ggplot2)
library(reshape2)
library(scater)

opt_cate_fl <- read.table('opt_final_TE_overrank_num_addmaf_cate.txt')
#opt_cate_fl <- read.table('opt_final_TE_overrank_num_addmaf_cate_fltSNPSVmaf01.txt') ##We do not use it as it has very low cases 
head(opt_cate_fl)
colnames(opt_cate_fl) <- c('te','number','maf','category')

##generate part of random number of the 0
#runif(500, min=0, max=100)

opt_cate_fl$category <- factor(opt_cate_fl$category,levels = c('low','high'))

##create LD plot
p <- ggplot(data=opt_cate_fl, aes(x=category, y=maf, fill=category)) +
  geom_boxplot() +
  #geom_bar(stat="identity", color = 'black',position=position_dodge()) +
  labs(x="\nLD", y = paste('MAF\n',sep = ' ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30), ##change all the text size
        axis.text.y = element_text(size=40,colour = "black",face = 'bold'),
        axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1,face = 'bold'),
        axis.title.x = element_text(size=35,colour = "black",face = 'bold'),
        axis.title.y = element_text(size=35,colour = "black",face = 'bold'),
        plot.margin = unit(c(1,1,1,1), "cm")) + ##change the margin more outer
  scale_fill_brewer(palette="Blues") +
  coord_cartesian(ylim = c(0, 0.5))
#p <- p + scale_fill_manual(values=c("lightpink", "indianred1"))
pdf('LD_maf_3cate.pdf',width = 7,height = 8)
#pdf('LD_maf_3cate_fltSNPSVmaf.pdf',width = 7,height = 8)
p
dev.off()

#significant test
opt_cate_low_fl <- opt_cate_fl[opt_cate_fl$category == 'low',]
opt_cate_high_fl <- opt_cate_fl[opt_cate_fl$category == 'high',]
wilcox.test(opt_cate_low_fl$maf,opt_cate_high_fl$maf)
##1.058e-11
##0.003026 for the maf01

##create the bar plot
dt <- read.table('opt_barplot_low_and_high_num.txt')
colnames(dt) <- c('cate','num')
p <- ggplot(data=dt, aes(x=cate, y=num,fill=cate)) +
  geom_bar(stat="identity", width=0.5) +
  #geom_bar(stat="identity", color = 'black',position=position_dodge()) +
  labs(x="\nLD", y = paste('TE number\n',sep = ' ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30), ##change all the text size
        axis.text.y = element_text(size=40,colour = "black",face = 'bold'),
        axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1,face = 'bold'),
        axis.title.x = element_text(size=35,colour = "black",face = 'bold'),
        axis.title.y = element_text(size=35,colour = "black",face = 'bold'),
        plot.margin = unit(c(1,1,1,1), "cm")) + ##change the margin more outer
  scale_fill_brewer(palette="Blues")

pdf('num_maf.pdf',width = 7,height = 8)
p
dev.off()


##create the density plot
t.test(opt_cate_fl[opt_cate_fl$category=='high',]$maf, opt_cate_fl[opt_cate_fl$category=='low',]$maf)
##p-value < 2.2e-16

##draw plot for the number of overank number
p <- ggplot(opt_cate_fl, aes(x=number)) + 
  #geom_histogram(color="darkblue", fill="lightblue")+
  geom_density(color="darkblue", fill="lightblue") + 
  labs(x="\nTE ranks over SNP median", y = 'Density\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=50,hjust = 0.5),
    axis.title = element_text(size =50, face="bold"),
    axis.text.x = element_text(colour = "black", size=50,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=50,colour = "black",face = 'bold'),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  geom_vline(xintercept = 150,linetype="dotted",colour='black',size=1)
  #xlim(0, 300)
#p

pdf('density_112821.pdf',width = 10,height = 10)
p
dev.off()



#################
##udpating 062024
##we will define category to build the density plot
##as the previous plotting with cate did not have enough sv cases
##as if we use the MAF to categorize them, it only returns a small number of caese
opt_cate_fl <- read.table('opt_final_TE_overrank_num_fltSNPSVmaf01.txt')
head(opt_cate_fl)
colnames(opt_cate_fl) <- c('te','number')
opt_cate_fl$category <- ifelse(opt_cate_fl$number >= 150, 'high','low')

ipt_dt <- read.delim('opt_SNP_SV_maf.txt',header = F)
head(ipt_dt)
colnames(ipt_dt) <- c('loc','MAF','cate')
table(ipt_dt$cate)

shared_sv_loc <- intersect(ipt_dt$loc,opt_cate_fl$te)
length(shared_sv_loc)
ipt_shared_sv_dt <- ipt_dt[ipt_dt$loc %in% shared_sv_loc,]
dim(ipt_shared_sv_dt)
head(ipt_shared_sv_dt)

ipt_shared_sv_MAF01_dt <- ipt_shared_sv_dt[ipt_shared_sv_dt$MAF > 0.1,]
dim(ipt_shared_sv_MAF01_dt)



##draw plot for the number of overank number
p <- ggplot(opt_cate_fl, aes(x=number)) + 
  #geom_histogram(color="darkblue", fill="lightblue")+
  geom_density(color="darkblue", fill="lightblue") + 
  labs(x="\nTE ranks over SNP median", y = 'Density\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=50,hjust = 0.5),
    axis.title = element_text(size =50, face="bold"),
    axis.text.x = element_text(colour = "black", size=50,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=50,colour = "black",face = 'bold'),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  geom_vline(xintercept = 150,linetype="dotted",colour='black',size=1)


pdf('density_fltMAF01_062024.pdf',width = 10,height = 10)
p
dev.off()

merged_dt <- merge(opt_cate_fl,ipt_shared_sv_MAF01_dt,by.x = 'te',by.y = 'loc')
head(merged_dt)
merged_dt$category <- factor(merged_dt$category,levels = c('low','high'))

##create LD plot
p <- ggplot(data=merged_dt, aes(x=category, y=MAF, fill=category)) +
  geom_boxplot() +
  #geom_bar(stat="identity", color = 'black',position=position_dodge()) +
  labs(x="\nLD", y = paste('MAF\n',sep = ' ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30), ##change all the text size
        axis.text.y = element_text(size=40,colour = "black",face = 'bold'),
        axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1,face = 'bold'),
        axis.title.x = element_text(size=35,colour = "black",face = 'bold'),
        axis.title.y = element_text(size=35,colour = "black",face = 'bold'),
        plot.margin = unit(c(1,1,1,1), "cm")) + ##change the margin more outer
  scale_fill_brewer(palette="Blues") +
  coord_cartesian(ylim = c(0, 0.5))
#p <- p + scale_fill_manual(values=c("lightpink", "indianred1"))
#pdf('LD_maf_3cate.pdf',width = 7,height = 8)
pdf('LD_maf_3cate_fltSNPSVmaf.pdf',width = 7,height = 8)
p
dev.off()







##updating 102923
##plot the MAF distribution between SNP and TE
ipt_dt <- read.delim('opt_SNP_SV_maf.txt',header = F)
head(ipt_dt)
colnames(ipt_dt) <- c('loc','MAF','cate')
table(ipt_dt$cate)

p <- ggplot(ipt_dt, aes(x=MAF,fill = cate)) + 
  geom_density(alpha=0.4) + 
  #geom_histogram(color="darkblue", fill="lightblue")+
  #geom_density(color="darkblue", fill="lightblue") + 
  labs(x="\nMAF", y = 'Density\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=50,hjust = 0.5),
    axis.title = element_text(size =50, face="bold"),
    axis.text.x = element_text(colour = "black", size=50,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=50,colour = "black",face = 'bold'),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  coord_cartesian(xlim = c(0, 0.5)) +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1))
  #geom_vline(xintercept = 150,linetype="dotted",colour='black',size=1)
#xlim(0, 300)

pdf('opt_distribution_MAF_SNP_SV.pdf',width = 10,height = 10)
p
dev.off()

##updating built the filtred SNP

ipt_snp_dt <- ipt_dt[ipt_dt$cate == 'SNP',]
ipt_snp_mafover0.1_dt <- ipt_snp_dt[ipt_snp_dt$MAF > 0.1,]
dim(ipt_snp_mafover0.1_dt)
ipt_sv_dt <- ipt_dt[ipt_dt$cate == 'SV',]
ipt_sv_mafover0.1_dt <- ipt_sv_dt[ipt_sv_dt$MAF > 0.1,]
dim(ipt_sv_mafover0.1_dt)

mergedt_dt <- rbind(ipt_snp_mafover0.1_dt,ipt_sv_dt)
mergedt_both_dt <- rbind(ipt_snp_mafover0.1_dt,ipt_sv_mafover0.1_dt)


p <- ggplot(mergedt_dt, aes(x=MAF,fill = cate)) + 
  geom_density(alpha=0.4) + 
  #geom_histogram(color="darkblue", fill="lightblue")+
  #geom_density(color="darkblue", fill="lightblue") + 
  labs(x="\nMAF", y = 'Density\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=50,hjust = 0.5),
    axis.title = element_text(size =50, face="bold"),
    axis.text.x = element_text(colour = "black", size=50,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=50,colour = "black",face = 'bold'),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
#geom_vline(xintercept = 150,linetype="dotted",colour='black',size=1)
#xlim(0, 300)

pdf('opt_distribution_MAF_SNP_SV_fltSNPmaf01.pdf',width = 10,height = 10)
p
dev.off()

p <- ggplot(mergedt_both_dt, aes(x=MAF,fill = cate)) + 
  geom_density(alpha=0.4) + 
  #geom_histogram(color="darkblue", fill="lightblue")+
  #geom_density(color="darkblue", fill="lightblue") + 
  labs(x="\nMAF", y = 'Density\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=50,hjust = 0.5),
    axis.title = element_text(size =50, face="bold"),
    axis.text.x = element_text(colour = "black", size=50,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=50,colour = "black",face = 'bold'),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  coord_cartesian(xlim = c(0, 0.5)) +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1))
  
#geom_vline(xintercept = 150,linetype="dotted",colour='black',size=1)
#xlim(0, 300)

pdf('opt_distribution_MAF_SNP_SV_fltSNPSVmaf01.pdf',width = 10,height = 10)
p
dev.off()










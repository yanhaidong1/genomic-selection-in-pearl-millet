##updating 062024 we will compare the GA for the two year comparison

##this script is to compare SNP and SV heri
##

library(reshape2)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

library(VennDiagram)
library(ggplot2)
library(reshape2)
library(scales)

library(dplyr)
library(ggplot2)
library(reshape2)
library(dplyr)

library(ggpubr)
library(cowplot)

library(scater)
library(RColorBrewer)



##step01 check average h2 of GA
ipt_heri_snp_dt <- read.delim('GA_GD.heritability_snp.txt')
head(ipt_heri_snp_dt)
ipt_heri_snp_dt <- ipt_heri_snp_dt[grepl('.GA',ipt_heri_snp_dt$Type),]
ipt_heri_snp_dt <- ipt_heri_snp_dt[,c('Type','h2')]
colnames(ipt_heri_snp_dt) <- c('Type','h2_snp')
dim(ipt_heri_snp_dt)
head(ipt_heri_snp_dt)

ipt_heri_sv_dt <- read.delim('GA_GD.heritability_sv.txt')
head(ipt_heri_sv_dt)
ipt_heri_sv_dt <- ipt_heri_sv_dt[,c('Type','h2')]
ipt_heri_sv_dt <- ipt_heri_sv_dt[grepl('.GA',ipt_heri_sv_dt$Type),]
colnames(ipt_heri_sv_dt) <- c('Type','h2_sv')
dim(ipt_heri_sv_dt)
head(ipt_heri_sv_dt)

mergedt_dt <- merge(ipt_heri_snp_dt,ipt_heri_sv_dt,by.x = 'Type',by.y = 'Type')
head(mergedt_dt)
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

##updating 061724 
##we will compare it with phenotype based heritability
head(mergedt_reshap_dt)
ipt_pheno_herit <- read.table('./wang_pheno_year.txt',header = F)
head(ipt_pheno_herit)
ipt_pheno_herit$Cate <- 'h2_pheno'
colnames(ipt_pheno_herit) <- c('Treat','Heritability','Cate')
ipt_pheno_herit <- ipt_pheno_herit[c('Treat','Cate','Heritability')]
head(ipt_pheno_herit)

merged_pheno_dt <- rbind(mergedt_reshap_dt,ipt_pheno_herit)
head(merged_pheno_dt)

##plot the box plot
p <- ggplot(merged_pheno_dt, aes(x=Cate, y=Heritability)) + 
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

pdf('opt_s1_compare_heri_SNP_SV_pheno.pdf',width = 8,height = 8)
p
dev.off()





##compare the sd
mergedt_reshap_snp_dt <- mergedt_reshap_dt[mergedt_reshap_dt$Cate == 'h2_snp',]
mergedt_reshap_sv_dt <- mergedt_reshap_dt[mergedt_reshap_dt$Cate == 'h2_sv',]
sd(mergedt_reshap_snp_dt$Heritability)
sd(mergedt_reshap_sv_dt$Heritability)

cate <- c("snp_sd",'sv_sd')
sd <- c(sd(mergedt_reshap_snp_dt$Heritability), sd(mergedt_reshap_sv_dt$Heritability))

# Creating a DataFrame using data.frame() function
df <- data.frame(cate = cate, sd = sd)

p <- ggplot(df, aes(x=cate, y=sd)) + 
  geom_bar(stat= 'identity') + 
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
pdf('opt_s1_compare_heri_sd_SNP_SV.pdf',width = 8,height = 8)
p
dev.off()




#############
##step02 plot which trait the SNPs having higher heri and SVs having higher heri
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
##step03
##divide the data into three parts based on the ranking
ipt_snp_dt <- read.delim('GA_GD.heritability_snp.txt',header = T)
head(ipt_snp_dt)

ipt_snp_dt <- ipt_snp_dt[grepl('.GA',ipt_snp_dt$Type),]
head(ipt_snp_dt)
dim(ipt_snp_dt)

ipt_snp_dt <- ipt_snp_dt[order(ipt_snp_dt$h2,decreasing = T),]

ipt_snp_ordered_traits <- ipt_snp_dt$Type

snp_part1_traits <- ipt_snp_ordered_traits[1:60]
snp_part2_traits <- ipt_snp_ordered_traits[61:120]


ipt_sv_dt <- read.delim('GA_GD.heritability_sv.txt',header = T)
head(ipt_sv_dt)

ipt_sv_dt <- ipt_sv_dt[grepl('.GA',ipt_sv_dt$Type),]
head(ipt_sv_dt)
dim(ipt_sv_dt)


ipt_sv_dt <- ipt_sv_dt[order(ipt_sv_dt$h2,decreasing = T),]

ipt_sv_ordered_traits <- ipt_sv_dt$Type

sv_part1_traits <- ipt_sv_ordered_traits[1:60]
sv_part2_traits <- ipt_sv_ordered_traits[61:120]



##Interesct
##98
grid.newpage()
pdf('VN_plot_part1.pdf',width = 6,height = 4)
draw.pairwise.venn(120,120,98,
                   c("SNP", "SV"),
                   col = rep("black", 2), fill = c ("goldenrod1","seagreen3"),cat.fontface = 'bold',
                   cex =2, cat.cex = 2, scaled = FALSE,cat.pos = c(0, 0)) ##12 10
dev.off()

##Intersect part 1 
intersect(snp_part1_traits,sv_part1_traits)
length(intersect(snp_part1_traits,sv_part1_traits))

##Interesct part 2
intersect(snp_part2_traits,sv_part2_traits)
length(intersect(snp_part2_traits,sv_part2_traits))

##Interesct part 1 and part2
intersect(snp_part1_traits,sv_part2_traits)
length(intersect(snp_part1_traits,sv_part2_traits))

##Interesct part 1 and part2
intersect(snp_part2_traits,sv_part1_traits)
length(intersect(snp_part2_traits,sv_part1_traits))

##check if the trait has negative correlation for the SNP and SV
ipt_snp_dt <- read.delim('GA_GD.heritability_snp.txt',header = T)
head(ipt_snp_dt)
ipt_snp_dt <- ipt_snp_dt[grepl('.GA',ipt_snp_dt$Type),]
head(ipt_snp_dt)
dim(ipt_snp_dt)

ipt_sv_dt <- read.delim('GA_GD.heritability_sv.txt',header = T)
head(ipt_sv_dt)
ipt_sv_dt <- ipt_sv_dt[grepl('.GA',ipt_sv_dt$Type),]
head(ipt_sv_dt)
dim(ipt_sv_dt)

mergedt_dt <- merge(ipt_snp_dt,ipt_sv_dt,by.x = 'Type',by.y = 'Type')
head(mergedt_dt)

colnames(mergedt_dt) <- c('Type','h2.snp','h2_SE.snp','h2.sv','h2_SE.sv')
head(mergedt_dt)

correlation <- cor(mergedt_dt$h2.snp, mergedt_dt$h2.sv, method = "spearman")
##-0.09237607


##plot the line plot to show the SNP and SV heritability of the same phenoytpe
mergedt_orderSNP_dt <- mergedt_dt[order(mergedt_dt$h2.snp,decreasing = T),]
mergedt_orderSNP_dt <- mergedt_orderSNP_dt[c('Type','h2.snp','h2.sv')]
mergedt_orderSNP_reshape_dt <- melt(mergedt_orderSNP_dt)

mergedt_orderSNP_reshape_dt$Type <- factor(mergedt_orderSNP_reshape_dt$Type ,levels = mergedt_orderSNP_dt$Type)
p<-ggplot(mergedt_orderSNP_reshape_dt, aes(x=Type, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    #axis.title = element_text(size =35,),
    axis.text.x = element_blank(),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    #axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
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
p
pdf(paste0('./opt_s3_snp_sv_h2_lineplot.pdf'),width = 12,height = 10)
print(p)
dev.off()

##for the sv
mergedt_orderSV_dt <- mergedt_dt[order(mergedt_dt$h2.sv,decreasing = T),]
mergedt_orderSV_dt <- mergedt_orderSV_dt[c('Type','h2.snp','h2.sv')]
mergedt_orderSV_reshape_dt <- melt(mergedt_orderSV_dt)

mergedt_orderSV_reshape_dt$Type <- factor(mergedt_orderSV_reshape_dt$Type ,levels = mergedt_orderSV_dt$Type)
p<-ggplot(mergedt_orderSV_reshape_dt, aes(x=Type, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    #axis.title = element_text(size =35,),
    axis.text.x = element_blank(),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    #axis.text.x = element_text(colour = "black", size=30,angle = 45, hjust = 1),  ##change the text to italic
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
pdf(paste0('./opt_s3_sv_snp_h2_lineplot.pdf'),width = 12,height = 10)
print(p)
dev.off()


########
##step04
##we will compare two years heritability
ipt_heri_snp_dt <- read.delim('GA_GD.heritability_snp.txt')
head(ipt_heri_snp_dt)
ipt_heri_snp_dt <- ipt_heri_snp_dt[grepl('.GA',ipt_heri_snp_dt$Type),]
ipt_heri_snp_dt <- ipt_heri_snp_dt[,c('Type','h2')]
colnames(ipt_heri_snp_dt) <- c('Type','h2')
dim(ipt_heri_snp_dt)
head(ipt_heri_snp_dt)

##for the snp
ipt_heri_snp_dt$year <- gsub('.+_','',ipt_heri_snp_dt$Type)
ipt_heri_snp_dt$year <- gsub('\\.GA','',ipt_heri_snp_dt$year)
ipt_heri_snp_dt$trait <- gsub('_20.+','',ipt_heri_snp_dt$Type)

ipt_heri_snp_dt$trait <- gsub('_Ags','_AgS',ipt_heri_snp_dt$trait)
ipt_heri_snp_dt$trait <- gsub('_AGS','_AgS',ipt_heri_snp_dt$trait)
ipt_heri_snp_dt$trait <- gsub('_Corrected','',ipt_heri_snp_dt$trait)
ipt_heri_snp_dt$trait <- gsub('EarlyS','Early',ipt_heri_snp_dt$trait)
ipt_heri_snp_dt$trait <- gsub('LateS','Late',ipt_heri_snp_dt$trait)
table(ipt_heri_snp_dt$trait)

ipt_heri_snp_dt$marker <- 'SNP'


##for the sv
ipt_heri_sv_dt <- read.delim('GA_GD.heritability_sv.txt')
head(ipt_heri_sv_dt)
ipt_heri_sv_dt <- ipt_heri_sv_dt[,c('Type','h2')]
ipt_heri_sv_dt <- ipt_heri_sv_dt[grepl('.GA',ipt_heri_sv_dt$Type),]
colnames(ipt_heri_sv_dt) <- c('Type','h2')
dim(ipt_heri_sv_dt)
head(ipt_heri_sv_dt)

ipt_heri_sv_dt$year <- gsub('.+_','',ipt_heri_sv_dt$Type)
ipt_heri_sv_dt$year <- gsub('\\.GA','',ipt_heri_sv_dt$year)
ipt_heri_sv_dt$trait <- gsub('_20.+','',ipt_heri_sv_dt$Type)

ipt_heri_sv_dt$trait <- gsub('_Ags','_AgS',ipt_heri_sv_dt$trait)
ipt_heri_sv_dt$trait <- gsub('_AGS','_AgS',ipt_heri_sv_dt$trait)
ipt_heri_sv_dt$trait <- gsub('_Corrected','',ipt_heri_sv_dt$trait)
ipt_heri_sv_dt$trait <- gsub('EarlyS','Early',ipt_heri_sv_dt$trait)
ipt_heri_sv_dt$trait <- gsub('LateS','Late',ipt_heri_sv_dt$trait)
table(ipt_heri_sv_dt$trait)

ipt_heri_sv_dt$marker <- 'SV'

merged_dt <- rbind(ipt_heri_snp_dt,ipt_heri_sv_dt)
merged_dt$markeryear <- paste0(merged_dt$marker,'_',merged_dt$year)
head(merged_dt)

##we will plot the two years information with dot plot

df <- data.frame(
  Group = rep(c("A", "B", "C"), each = 4),
  Condition = rep(c("Condition1", "Condition2", "Condition3", "Condition4"), times = 3),
  Value = c(3, 5, 2, 4, 6, 7, 3, 5, 4, 6, 5, 8)
)

p <- ggplot(merged_dt, aes(x = trait, y = h2, color = markeryear,shape = markeryear)) +
  #geom_point(position = position_dodge(width = 0.8), size = 3, shape = 16) +
  geom_point(size = 3)+
  theme_minimal() +
  labs(title = "", x = "Trait", y = "H2") +
  scale_color_manual(values = c("SNP_2011" = "blue", "SNP_2012" = "red", "SV_2011" = "green", "SV_2012" = "purple"))+
  scale_shape_manual(values = c("SNP_2011" = 16, "SNP_2012" = 16, "SV_2011" = 10, "SV_2012" = 10))+
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    #axis.title = element_text(size =35,),
    #axis.text.x = element_blank(),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=17,angle = 90, vjust = 0.5,hjust = 0.5),  ##change the text to italic
    axis.text.y = element_text(size=17,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )


pdf(paste0('./opt_s4_sv_snp_h2_dotplot_all_trait.pdf'),width = 20,height = 10)
print(p)
dev.off()
















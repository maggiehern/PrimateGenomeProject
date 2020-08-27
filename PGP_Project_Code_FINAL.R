#Packages needed
install.packages("ggplot2")
install.packages("moments")
install.packages("gridExtra")
install.packages("hexbin")
install.packages("readxl")


library(ggplot2)
library(moments)
library(gridExtra)
library(hexbin)
library(readxl)


#Data files needed. Import these into your global environment as excel files. 
#PGP_Data_FINAL
#Codes_and_Quantities

PGP_Data_FINAL<-read_excel("PGP_Data_FINAL.xlsx")
View(PGP_Data_FINAL)

Codes_and_Quantities <- read_excel("Codes_and_Quantities.xlsx")
View(Codes_and_Quantities)

#Log transformations for variables 

PGP_Data_FINAL$log.Mbases = log10(PGP_Data_FINAL$Mbases)
PGP_Data_FINAL[PGP_Data_FINAL$log.Mbases == -Inf,]$log.Mbases = 0

PGP_Data_FINAL$log.RIR_MP = log10(PGP_Data_FINAL$RIR_MP)
PGP_Data_FINAL[PGP_Data_FINAL$log.RIR_MP == -Inf,]$log.RIR_MP = 0

PGP_Data_FINAL$log.Mya_since_LCA = log10(PGP_Data_FINAL$Mya_since_LCA)

PGP_Data_FINAL$log.Importance_in_medical_research = log10(PGP_Data_FINAL$Importance_in_medical_research)
PGP_Data_FINAL[PGP_Data_FINAL$log.Importance_in_medical_research == -Inf,]$log.Importance_in_medical_research = 0

PGP_Data_FINAL$log.Geographical_range_sq_km = log10(as.numeric(PGP_Data_FINAL$Geographical_range_sq_km))

PGP_Data_FINAL$log.Frequency_in_captivity = log10(PGP_Data_FINAL$Frequency_in_captivity)
PGP_Data_FINAL[PGP_Data_FINAL$log.Frequency_in_captivity == -Inf,]$log.Frequency_in_captivity = 0

#Code for independent 2-group Mann-Whitney U Test

wilcox.test(PGP_Data_FINAL$RIR_MP ~PGP_Data_FINAL$Genomic_Data)
tapply(PGP_Data_FINAL$RIR_MP, PGP_Data_FINAL$Genomic_Data, median)
tapply(PGP_Data_FINAL$RIR_MP, PGP_Data_FINAL$Genomic_Data, mean)
tapply(PGP_Data_FINAL$RIR_MP, PGP_Data_FINAL$Genomic_Data, sd)
mean(PGP_Data_FINAL$RIR_MP[PGP_Data_FINAL$Genomic_Data=="Absent"])
mean(PGP_Data_FINAL$RIR_MP[PGP_Data_FINAL$Genomic_Data=="Present"])
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.RIR_MP)) + 
  geom_boxplot(fill="dark grey") +theme_minimal() + ylim(0,5)

wilcox.test(PGP_Data_FINAL$Mya_since_LCA ~PGP_Data_FINAL$Genomic_Data)
tapply(PGP_Data_FINAL$Mya_since_LCA, PGP_Data_FINAL$Genomic_Data, median)
tapply(PGP_Data_FINAL$Mya_since_LCA, PGP_Data_FINAL$Genomic_Data, mean)
tapply(PGP_Data_FINAL$Mya_since_LCA, PGP_Data_FINAL$Genomic_Data, sd)
mean(PGP_Data_FINAL$Mya_since_LCA[PGP_Data_FINAL$Genomic_Data=="Absent"])
mean(PGP_Data_FINAL$Mya_since_LCA[PGP_Data_FINAL$Genomic_Data=="Present"])
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=Mya_since_LCA)) + 
  geom_boxplot(fill="dark grey") +theme_minimal() +ylim(0,80)

wilcox.test(PGP_Data_FINAL$Importance_in_medical_research ~PGP_Data_FINAL$Genomic_Data)
tapply(PGP_Data_FINAL$Importance_in_medical_research, PGP_Data_FINAL$Genomic_Data, median)
tapply(PGP_Data_FINAL$Importance_in_medical_research, PGP_Data_FINAL$Genomic_Data, mean)
tapply(PGP_Data_FINAL$Importance_in_medical_research, PGP_Data_FINAL$Genomic_Data, sd)
mean(PGP_Data_FINAL$Importance_in_medical_research[PGP_Data_FINAL$Genomic_Data=="Absent"])
mean(PGP_Data_FINAL$Importance_in_medical_research[PGP_Data_FINAL$Genomic_Data=="Present"])
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Importance_in_medical_research)) + 
  geom_boxplot(fill="dark grey") +theme_minimal() +ylim(0,5)

PGP_Data_NA_Geo_Removed<- PGP_Data_FINAL[- grep("NA", PGP_Data_FINAL$Geographical_range_sq_km),]
wilcox.test(as.numeric(PGP_Data_FINAL$Geographical_range_sq_km) ~PGP_Data_FINAL$Genomic_Data)
tapply(as.numeric(PGP_Data_NA_Geo_Removed$Geographical_range_sq_km), PGP_Data_NA_Geo_Removed$Genomic_Data, median)
tapply(as.numeric(PGP_Data_NA_Geo_Removed$Geographical_range_sq_km), PGP_Data_NA_Geo_Removed$Genomic_Data, mean)
tapply(as.numeric(PGP_Data_NA_Geo_Removed$Geographical_range_sq_km), PGP_Data_NA_Geo_Removed$Genomic_Data, sd)
mean(PGP_Data_FINAL$Geographical_range_sq_km[PGP_Data_FINAL$Genomic_Data=="Absent"])
mean(PGP_Data_FINAL$Geographical_range_sq_km[PGP_Data_FINAL$Genomic_Data=="Present"])
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Geographical_range_sq_km)) + 
  geom_boxplot(fill="dark grey") +theme_minimal() +ylim(0,5)

wilcox.test(PGP_Data_FINAL$Frequency_in_captivity ~PGP_Data_FINAL$Genomic_Data)
tapply(PGP_Data_FINAL$Frequency_in_captivity, PGP_Data_FINAL$Genomic_Data, median)
tapply(PGP_Data_FINAL$Frequency_in_captivity, PGP_Data_FINAL$Genomic_Data, mean)
tapply(PGP_Data_FINAL$Frequency_in_captivity, PGP_Data_FINAL$Genomic_Data, sd)
mean(PGP_Data_FINAL$Frequency_in_captivity[PGP_Data_FINAL$Genomic_Data=="Absent"])
mean(PGP_Data_FINAL$Frequency_in_captivity[PGP_Data_FINAL$Genomic_Data=="Present"])
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Frequency_in_captivity)) + 
  geom_boxplot(fill="dark grey") +theme_minimal() +ylim(0,5)

#code for ANOVAs
AN1 <- aov(Mbases ~ Red_List_status, data=PGP_Data_FINAL)
summary(AN1)
TukeyHSD(AN1)

PGP_Data_FINAL$Red_List_status <- factor(PGP_Data_FINAL$Red_List_status , levels=c("LC", "NT", "VU", "EN", "CR", "DD"))

ggplot(PGP_Data_FINAL, aes(y=log.Mbases, x=Red_List_status)) + 
  geom_boxplot(fill="dark grey") +theme_minimal()

PGP_Data_Activity_NA_Removed <- PGP_Data_FINAL[- grep("NA", PGP_Data_FINAL$Activity_pattern),]
AN2 <- aov(Mbases ~ Activity_pattern, data=PGP_Data_Activity_NA_Removed)
summary(AN2)
TukeyHSD(AN2)

PGP_Data_Activity_NA_Removed$Activity_pattern <- factor(PGP_Data_Activity_NA_Removed$Activity_pattern , levels=c("Diurnal", "Nocturnal", "Cathemeral"))
ggplot(PGP_Data_Activity_NA_Removed, aes(y=log.Mbases, x=Activity_pattern)) + 
  geom_boxplot(fill="dark grey") +theme_minimal()

# Sum up Mbases by genus for heat map
genus.counts = aggregate(PGP_Data_FINAL$Mbases, by=list(c(PGP_Data_FINAL$Genus_category)), FUN=sum)
names(genus.counts) = c("Genus", "Mbases")

genus.counts$log.Mbases = log10(genus.counts$Mbases)
genus.counts[genus.counts$log.Mbases == -Inf,]$log.Mbases = 0

ggplot(genus.counts, aes(x=0, y=Genus, col=log.Mbases)) + geom_point(pch=15, size=5) +
  scale_color_gradient(low='#e8e6f4', high='#483d8b') + theme_bw()

#heatmap qualtitative data interview codes frequency

ggplot(Codes_and_Quantities, aes(x=0, y=Code, col=Quantity)) + geom_point(pch=15, size=5) +
  scale_color_gradient(low='#e8e6f4', high='#8b483d') + theme_bw()

#Linear regressions with all data  

linearMod1 <- lm(PGP_Data_FINAL$log.Mbases~PGP_Data_FINAL$log.RIR_MP)
summary(linearMod1)
hist(residuals(linearMod1))
qqnorm(residuals(linearMod1))
shapiro.test(residuals(linearMod1))
agostino.test(residuals(linearMod1))

linearMod2 <- lm(PGP_Data_FINAL$log.Mbases~PGP_Data_FINAL$log.Importance_in_medical_research)
summary(linearMod2)
hist(residuals(linearMod2))
qqnorm(residuals(linearMod2))
shapiro.test(residuals(linearMod2))
agostino.test(residuals(linearMod2))

linearMod3 <- lm(PGP_Data_FINAL$log.Mbases~PGP_Data_FINAL$log.Mya_since_LCA)
summary(linearMod3)
hist(residuals(linearMod3))
qqnorm(residuals(linearMod3))
shapiro.test(residuals(linearMod3))
agostino.test(residuals(linearMod3))

linearMod4 <- lm(PGP_Data_FINAL$log.Mbases~PGP_Data_FINAL$log.Geographical_range_sq_km)
summary(linearMod4)
hist(residuals(linearMod4))
qqnorm(residuals(linearMod4))
shapiro.test(residuals(linearMod4))
agostino.test(residuals(linearMod4))

linearMod5 <- lm(PGP_Data_FINAL$log.Mbases~PGP_Data_FINAL$log.Frequency_in_captivity)
summary(linearMod5)
hist(residuals(linearMod5))
qqnorm(residuals(linearMod5))
shapiro.test(residuals(linearMod5))
agostino.test(residuals(linearMod5))

#Data Deficient (DD) variable removed since it could not be placed within the logical order of the different statuses 
PGP_Data_DD_Removed<- subset(PGP_Data_FINAL, Red_List_status != "DD")
PGP_Data_DD_Removed$log.Mbases = log10(PGP_Data_DD_Removed$Mbases)
PGP_Data_DD_Removed[PGP_Data_DD_Removed$log.Mbases == -Inf,]$log.Mbases = 0

linearMod6 <- lm(PGP_Data_DD_Removed$log.Mbases~PGP_Data_DD_Removed$Red_list_code)
summary(linearMod6)
hist(residuals(linearMod6))
qqnorm(residuals(linearMod6))
shapiro.test(residuals(linearMod6))
agostino.test(residuals(linearMod6))

# Linear regressions with subset of species with genomic data present 

PGP_Data_Genomic_Present <- subset(PGP_Data_FINAL, Genomic_Data == "Present")

linearMod7 <- lm(PGP_Data_Genomic_Present$log.Mbases~PGP_Data_Genomic_Present$log.RIR_MP)
summary(linearMod7)
hist(residuals(linearMod7))
qqnorm(residuals(linearMod7))
shapiro.test(residuals(linearMod7))
agostino.test(residuals(linearMod7))

linearMod8 <- lm(PGP_Data_Genomic_Present$log.Mbases~PGP_Data_Genomic_Present$log.Importance_in_medical_research)
summary(linearMod8)
hist(residuals(linearMod8))
qqnorm(residuals(linearMod8))
shapiro.test(residuals(linearMod8))
agostino.test(residuals(linearMod8))

linearMod9 <- lm(PGP_Data_Genomic_Present$log.Mbases~PGP_Data_Genomic_Present$log.Frequency_in_captivity)
summary(linearMod9)
hist(residuals(linearMod9))
qqnorm(residuals(linearMod9))
shapiro.test(residuals(linearMod9))
agostino.test(residuals(linearMod9))

linearMod10 <- lm(PGP_Data_Genomic_Present$log.Mbases~PGP_Data_Genomic_Present$log.Mya_since_LCA)
summary(linearMod10) 
hist(residuals(linearMod10))
qqnorm(residuals(linearMod10))
shapiro.test(residuals(linearMod10))
agostino.test(residuals(linearMod10))

linearMod11 <- lm(formula = PGP_Data_Genomic_Present$log.Mbases ~ as.numeric(PGP_Data_Genomic_Present$log.Geographical_range_sq_km))
summary(linearMod11)
hist(residuals(linearMod11))
qqnorm(residuals(linearMod11))
shapiro.test(residuals(linearMod11))
agostino.test(residuals(linearMod11))

PGP_Data_Genomic_Present_noDD <- PGP_Data_Genomic_Present[!is.na(PGP_Data_Genomic_Present$trimmed.Red_List_status), ]

linearMod12 <- lm(PGP_Data_Genomic_Present_No_DD$log.Mbases~PGP_Data_Genomic_Present_No_DD$Red_list_code)
summary(linearMod12) 
hist(residuals(linearMod12))
qqnorm(residuals(linearMod12))
shapiro.test(residuals(linearMod12))
agostino.test(residuals(linearMod12))


#Figures with regression line through entire dataset and subset with genomic data 
# For each variable, the top code makes a figure with no legend, bottom code with legend

p1<-ggplot(PGP_Data_FINAL, aes(x=log.RIR_MP, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.RIR_MP, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +
  xlab("RIR_MP") + ylab("Megabases of Genomic Data") +ylim(-1,11)+ theme(legend.position="none") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.title.y = element_blank())

ggplot(PGP_Data_FINAL, aes(x=log.RIR_MP, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.RIR_MP, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +
  xlab("RIR_MP") + ylab("Megabases of Genomic Data") +ylim(-1,11)+ theme(legend.position="right")


p2<-ggplot(PGP_Data_FINAL, aes(x=log.Importance_in_medical_research, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  geom_smooth(aes(x=log.Importance_in_medical_research, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+ ylim(-1,11)+
  xlab("# of Medical Papers Published") + ylab("Megabases of Genomic Data") + theme(legend.position="none") + theme(axis.text.y = element_blank(), 
                                                                                                                    axis.ticks.y = element_blank(), 
                                                                                                                    axis.title.y = element_blank())

ggplot(PGP_Data_FINAL, aes(x=log.Importance_in_medical_research, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  geom_smooth(aes(x=log.Importance_in_medical_research, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+ 
  xlab("# of Medical Papers Published") + ylab("Megabases of Genomic Data") + ylim(-1,11)+theme(legend.position="right")



p3<-ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Mya_since_LCA, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +
  xlab("MYA since LCA with Humans") + ylab("Megabases of Genomic Data")+ylim(-1,11) +theme(legend.position="none")

ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Mya_since_LCA, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +
  xlab("MYA since LCA with Humans") + ylab("Megabases of Genomic Data")+ylim(-1,11) +theme(legend.position="right")




p4<-ggplot(PGP_Data_FINAL, aes(x=as.numeric(log.Geographical_range_sq_km), y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Geographical_range_sq_km, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +ylim(-1,11)+
  xlab(bquote('Geographic Range ('*km^2*')')) + ylab("Megabases of Genomic Data") + theme(legend.position="none")+theme(axis.text.y = element_blank(), 
                                                                                                                        axis.ticks.y = element_blank(), 
                                                                                                                        axis.title.y = element_blank())

ggplot(PGP_Data_FINAL, aes(x=as.numeric(log.Geographical_range_sq_km), y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Geographical_range_sq_km, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +ylim(-1,11)+
  xlab(bquote('Geographic Range ('*km^2*')')) + ylab("Megabases of Genomic Data") + theme(legend.position="right")

p5<-ggplot(PGP_Data_FINAL, aes(x=as.numeric(log.Frequency_in_captivity), y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Frequency_in_captivity, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +ylim(-1,11)+
  xlab('Frequency in Zoos') + ylab("Megabases of Genomic Data") + theme(legend.position="none") 

ggplot(PGP_Data_FINAL, aes(x=as.numeric(log.Frequency_in_captivity), y=log.Mbases)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + 
  geom_smooth(aes(x=log.Frequency_in_captivity, y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm", size=2)+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +ylim(-1,11)+
  xlab('Frequency in Zoos') + ylab("Megabases of Genomic Data") + theme(legend.position="right")

p6<-ggplot(PGP_Data_DD_Removed, aes(x= as.character(Red_list_code), y=log.Mbases)) + 
  geom_violin(fill="gray90", color="gray90")+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA)+
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  geom_smooth(aes(x=as.character(Red_list_code), y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm")+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+
  
  theme(legend.position="right") + ylim(-1,11)+
  xlab('IUCN Red List Status') + ylab("Megabases of Genomic Data") +
  scale_x_discrete(labels=c("1" = "LC", "2" = "NT", "3" = "VU", "4" = "EN", "5"= "CR")) + theme(legend.position="none") +theme(axis.text.y = element_blank(), 
                                                                                                                               axis.ticks.y = element_blank(), 
                                                                                                                               axis.title.y = element_blank())

ggplot(PGP_Data_DD_Removed, aes(x= as.character(Red_list_code), y=log.Mbases)) + 
  geom_violin(fill="gray90", color="gray90")+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA)+
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  geom_smooth(aes(x=as.character(Red_list_code), y=log.Mbases, group=Genomic_Data, color=factor(Genomic_Data)), method="lm")+
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+
  theme(legend.position="right") + ylim(-1,11)+
  xlab('IUCN Red List Status') + ylab("Megabases of Genomic Data") +
  scale_x_discrete(labels=c("1" = "LC", "2" = "NT", "3" = "VU", "4" = "EN", "5"= "CR")) + theme(legend.position="right")

#Figure 1 panels 

grid.arrange(p1,p2,ncol=2)

grid.arrange(p3,p4,ncol=2)

grid.arrange(p5,p6,ncol=2)

#Supplementary Figure 

ggplot(PGP_Data_FINAL, aes(x= as.character(PGP_Data_FINAL$Red_list_code), y=PGP_Data_FINAL$log.Mbases)) + 
  geom_violin(fill="dark grey", color="dark grey") + geom_hex(bins = 60) 

PGP_Data_Activity_NA_Removed <- subset(PGP_Data_FINAL, Activity_pattern != "NA")
PGP_Data_Activity_NA_Removed$log.Mbases = log10(PGP_Data_Activity_NA_Removed$Mbases)
PGP_Data_Activity_NA_Removed[PGP_Data_Activity_NA_Removed$log.Mbases == -Inf,]$log.Mbases = 0

ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, y=log.Mbases)) +
  geom_violin(fill="gray90", color="gray90") + stat_binhex(bins = 60) +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme(legend.position="right") + scale_y_continuous(limits = c(-8, 8), breaks = seq(-8,8, by =2))+
  xlab('Activity Pattern') + ylab("Megabases of Genomic Data") + theme(legend.position="none")+theme_minimal() +theme(legend.position="right")

#Figures running the variables against one another 

ggplot(PGP_Data_FINAL, aes(x=log.Importance_in_medical_research, y=log.RIR_MP)) +
  geom_hex(bins = 60) +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="right") +
  theme(legend.position="none")

ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.RIR_MP)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="right")+
  theme(legend.position="none") + xlim(0,2)

ggplot(PGP_Data_FINAL, aes(x=log.Geographical_range_sq_km, y=log.RIR_MP)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="right") +
  theme(legend.position="none") + ylim(-1,5)+xlim(0,7)

ggplot(PGP_Data_FINAL, aes(x=log.Frequency_in_captivity, y=log.RIR_MP)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="right") +
  theme(legend.position="none")

ggplot(PGP_Data_DD_Removed, aes(x=as.character(Red_list_code), y=log.RIR_MP)) +
  geom_violin(fill="gray90", color="gray90")+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="none")

ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, y=log.RIR_MP)) +
  geom_violin(fill="gray90", color="gray90") + stat_binhex(bins = 60) +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme_minimal() +theme(legend.position="none")

ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.Importance_in_medical_research)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none") +
  xlim(0,2)

ggplot(PGP_Data_FINAL, aes(x=log.Geographical_range_sq_km, y=log.Importance_in_medical_research)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none") +
  xlim(0,7)

ggplot(PGP_Data_FINAL, aes(x=log.Frequency_in_captivity, y=log.Importance_in_medical_research)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none")

ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.Geographical_range_sq_km)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none")+
  ylim(0,8)+xlim(0,2)

ggplot(PGP_Data_FINAL, aes(x=log.Frequency_in_captivity, y=log.Geographical_range_sq_km)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none") +ylim(0,8)

ggplot(PGP_Data_DD_Removed, aes(x=as.character(Red_list_code), y=log.Geographical_range_sq_km)) + 
  geom_violin(fill="gray90", color="gray90")+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="none") + ylim(0,8) 
  
ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, log.Geographical_range_sq_km)) +
  geom_violin(fill="gray90", color="gray90") + stat_binhex(bins = 60) +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme_minimal() +theme(legend.position="none")+ylim(0,8)

ggplot(PGP_Data_FINAL, aes(x=log.Mya_since_LCA, y=log.Frequency_in_captivity)) +
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal() +theme(legend.position="none")+ylim(-1,5) +xlim(0,2)

ggplot(PGP_Data_DD_Removed, aes(x=as.character(Red_list_code), y=log.Frequency_in_captivity)) +
  geom_violin(fill="gray90", color="gray90")+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="none")

ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, y=log.Frequency_in_captivity)) +
  geom_violin(fill="gray90", color="gray90") + stat_binhex(bins = 60) +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme_minimal() +theme(legend.position="none")

#used geom jitter here instead

ggplot(PGP_Data_DD_Removed, aes(x=as.character(Red_list_code), y=log.Importance_in_medical_research)) +
  geom_violin(fill="gray90", color="gray90")+geom_jitter(color="gray60")+
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="none")

ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, y=log.Importance_in_medical_research)) +
  geom_violin(fill="gray90", color="gray90") + geom_jitter(color="gray60") +
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme_minimal() +theme(legend.position="none")

ggplot(PGP_Data_DD_Removed, aes(x=as.character(Red_list_code), y=log.Mya_since_LCA)) +
  geom_violin(fill="gray90", color="gray90")+ geom_jitter(color="gray60")+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method = "lm", se=FALSE, color="#808b3d", aes(group=1)) +
  scale_color_manual(values=c("transparent", "#8b3d80")) + theme_minimal()+theme(legend.position="none") +ylim(0,2)

ggplot(PGP_Data_Activity_NA_Removed, aes(x=Activity_pattern, y=log.Mya_since_LCA)) +
  geom_violin(fill="gray90", color="gray90") +  geom_jitter(color="gray60")+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  theme_minimal() +theme(legend.position="none") +ylim(0,2)

  #linear models

linearMod13 <- lm(log.RIR_MP~log.Importance_in_medical_research, data=PGP_Data_FINAL)
summary(linearMod13) 

linearMod14 <- lm(log.RIR_MP~log.Mya_since_LCA, data=PGP_Data_FINAL)
summary(linearMod14) 

linearMod15 <- lm(log.RIR_MP~log.Geographical_range_sq_km, data=PGP_Data_FINAL)
summary(linearMod15) 

linearMod16 <- lm(log.RIR_MP~log.Frequency_in_captivity, data=PGP_Data_FINAL)
summary(linearMod16) 

linearMod17 <- lm(log.RIR_MP~Red_list_code, data=PGP_Data_DD_Removed)
summary(linearMod17)

linearMod18 <- lm(log.Importance_in_medical_research~log.Mya_since_LCA, data=PGP_Data_FINAL)
summary(linearMod18)

linearMod19 <- lm(log.Importance_in_medical_research~log.Geographical_range_sq_km, data=PGP_Data_FINAL)
summary(linearMod19)

linearMod20 <- lm(log.Importance_in_medical_research~log.Frequency_in_captivity, data=PGP_Data_FINAL)
summary(linearMod20)

linearMod21 <- lm(log.Importance_in_medical_research~Red_list_code, data=PGP_Data_DD_Removed)
summary(linearMod21)

linearMod22 <- lm(log.Geographical_range_sq_km~log.Mya_since_LCA, data=PGP_Data_FINAL)
summary(linearMod22)

linearMod23 <- lm(log.Geographical_range_sq_km~log.Frequency_in_captivity, data=PGP_Data_FINAL)
summary(linearMod23)

linearMod24 <- lm(log.Geographical_range_sq_km~Red_list_code, data=PGP_Data_DD_Removed)
summary(linearMod24)

linearMod25 <- lm(log.Frequency_in_captivity ~log.Mya_since_LCA, data=PGP_Data_FINAL)
summary(linearMod25)

linearMod26 <- lm(log.Frequency_in_captivity~Red_list_code, data=PGP_Data_DD_Removed)
summary(linearMod26)

linearMod27 <- lm(log.Mya_since_LCA~Red_list_code, data=PGP_Data_DD_Removed)
summary(linearMod27)



#cumulative percentage plot

PGP_Data_FINAL$Percentage = (PGP_Data_FINAL$Mbases)/sum(PGP_Data_FINAL$Mbases)*100

PGP_Data_max <- subset(PGP_Data_FINAL, Mbases > 5000000)

PGP_Data_min <- subset(PGP_Data_FINAL, Mbases < 5000000)

PGP_Data_max_wOther = add_case(PGP_Data_max, Species= "Other_taxa", Percentage=(sum(PGP_Data_min$Mbases)*100/sum(PGP_Data$Mbases)))

PGP_Data_max_wOther_cumulative = add_column(PGP_Data_max_wOther, Cumulative_percentage=c("34.923639", "50.726223", "62.211129","69.556583", "74.15856", "100"))

PGP_Data_max_wOther_cumulative2 = add_column(PGP_Data_max_wOther_cumulative, Position=c("1", "2", "3","4", "5", "6"))

ggplot(PGP_Data_max_wOther_cumulative2, aes(x=Species, y=Percentage)) + 
  geom_bar(stat="identity", fill="#808b3d") + 
  geom_line(size=1.25, color= "darkgrey", PGP_Data_max_wOther_cumulative2, mapping = aes(x=Species, y=as.numeric(Cumulative_percentage), group=1)) +
  geom_point(size=2, color= "darkgrey", PGP_Data_max_wOther_cumulative2, mapping = aes(x=Species, y=as.numeric(Cumulative_percentage))) + 
  coord_flip() + 
  scale_x_discrete(limits=c("Other_taxa","Macaca_fascicularis","Pan_troglodytes","Chlorocebus_sabaeus","Papio_anubis","Macaca_mulatta")) +
  ylim(0,100) + theme_minimal()

scale_fill_manual(values=c("#483d8b", "#504490", "#534691", "#574a94", "#5b4e97", "lightgrey"))  

#logistical regeressions on entire dataset

mylogit1 <- glm(Genomic_Data_Code ~ log.Importance_in_medical_research + 
                  log.Mya_since_LCA + log.RIR_MP + log.Frequency_in_captivity + 
                  as.numeric(log.Geographical_range_sq_km) +Activity_pattern + Red_List_status, data = PGP_Data_FINAL, family = "binomial")
summary(mylogit1)

#generalized linear models
#entire dataset

GLMGau1<-glm(formula = log.Mbases ~ log.Importance_in_medical_research + 
               log.Mya_since_LCA + log.RIR_MP + log.Frequency_in_captivity + 
               as.numeric(log.Geographical_range_sq_km) +Activity_pattern + Red_List_status, 
             family = "gaussian", data = PGP_Data_FINAL)
summary(GLMGau1)
qqnorm(residuals(GLMGau1))
shapiro.test(residuals(GLMGau1))
agostino.test(residuals(GLMGau1))

#subset with genomic data

GLMGau2<-glm(formula = log.Mbases ~ log.Importance_in_medical_research + 
               log.Mya_since_LCA + log.RIR_MP + log.Frequency_in_captivity + 
               as.numeric(log.Geographical_range_sq_km) +Activity_pattern + Red_List_status, 
             family = "gaussian", data = PGP_Data_Genomic_Present)
summary(GLMGau2)
qqnorm(residuals(GLMGau2))
shapiro.test(residuals(GLMGau2))
agostino.test(residuals(GLMGau2))

AIC(GLMGau1, GLMGau2)

#Chi-square test

table(PGP_Data_FINAL$Genomic_Data, PGP_Data_FINAL$Red_List_status)

chisq_IUCN_Genomic<-chisq.test(PGP_Data_FINAL$Genomic_Data, PGP_Data_FINAL$Red_List_status, correct=FALSE)
chisq_IUCN_Genomic
chisq_IUCN_Genomic$observed
chisq_IUCN_Genomic$residuals

library(corrplot)
corrplot(chisq_IUCN_Genomic$residuals, is.cor = FALSE)

#removing NA values from activity pattern
PGP_Data_NA_Activity_Removed<- PGP_Data_FINAL[- grep("NA", PGP_Data_FINAL$Activity_pattern),]

chisq_Activity_Genomic<-chisq.test(PGP_Data_NA_Activity_Removed$Genomic_Data, PGP_Data_NA_Activity_Removed$Activity_pattern, correct=FALSE)
chisq_Activity_Genomic
chisq_Activity_Genomic$observed
chisq_Activity_Genomic$residuals

#boxplots

ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Research_intensity_strict)) + geom_boxplot(fill="grey")
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Research_intensity_relaxed)) + geom_boxplot(fill="grey")
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Importance_in_medical_research)) + geom_boxplot(fill="grey")
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Frequency_in_captivity)) + geom_boxplot(fill="grey")
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=Mya_since_LCA)) + geom_boxplot(fill="grey")
ggplot(PGP_Data_FINAL, aes(x=Genomic_Data, y=log.Geographical_range_sq_km)) + geom_boxplot(fill="grey")

#Leave One Out Validation Approach 

library(lattice)
library(caret)

data_ctrl <- trainControl(method = "LOOCV", savePredictions = "final")

  #Entire dataset 

LOOCV1 <- train(log.Mbases ~ log.Importance_in_medical_research + 
                       log.Mya_since_LCA + log.RIR_MP + log.Frequency_in_captivity + 
                       as.numeric(log.Geographical_range_sq_km) +Activity_pattern + Red_List_status,  
                     data = PGP_Data_FINAL,                        
                     trControl = data_ctrl,             
                     method = "lm",
                     na.action = na.pass)

write.csv(LOOCV1$pred, "LOOCV1_Data.csv")
LOOCV1_Data <- read.csv("C:/Users/maggi/Box/Penn State/Master's/Primate Genome R Project/LOOCV1_Data.csv")

ggplot(LOOCV1_Data, aes(x=pred, y=obs))+geom_point()

LOOCV_p1 <- ggplot(LOOCV1_Data, aes(x=obs, y=pred))+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + theme_minimal()+xlim(-1,11)+
  theme(legend.position="right") +expand_limits(x = 0, y = 0) +theme(legend.position="none") +ylim(-1,11)

ggplot(LOOCV1_Data, aes(x=obs, y=pred))+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + theme_minimal()+xlim(-1,11)+
  theme(legend.position="right") +expand_limits(x = 0, y = 0) +theme(legend.position="right") +ylim(-1,11)

linearModLOOCV1 <- lm(formula = LOOCV1_Data$pred ~ LOOCV1_Data$obs)
summary(linearModLOOCV1)


  #Only species with genomic data present

LOOCV2 <- train(log.Mbases ~ log.Importance_in_medical_research + 
                  log.Mya_since_LCA + log.RIR_MP + log.Frequency_in_captivity + 
                  as.numeric(log.Geographical_range_sq_km) +Activity_pattern + Red_List_status,  
                data = PGP_Data_Genomic_Present,                        
                trControl = data_ctrl,             
                method = "lm",
                na.action = na.pass)

write.csv(LOOCV2$pred, "LOOCV2_Data.csv")
LOOCV2_Data <- read.csv("C:/Users/maggi/Box/Penn State/Master's/Primate Genome R Project/LOOCV2_Data.csv")


LOOCV_p2 <- ggplot(LOOCV2_Data, aes(x=obs, y=pred))+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + theme_minimal()+ylim(-1,11)+ xlim(-1,11)+
  theme(legend.position="right")+expand_limits(x = 0, y = 0) +theme(legend.position="none") +theme(axis.text.y = element_blank(), 
                                                                                                   axis.ticks.y = element_blank(), 
                                                                                                   axis.title.y = element_blank())

ggplot(LOOCV2_Data, aes(x=obs, y=pred))+
  geom_hex(bins = 60)+
  scale_fill_gradientn(colours=c("gray60","#8b483d"),name = "Number of Species",na.value=NA) +
  geom_smooth(method=lm, color="#808b3d", size=2) + theme_minimal()+ylim(-1,11)+ xlim(-1,11)+
  theme(legend.position="right")+expand_limits(x = 0, y = 0) +theme(legend.position="right") +theme(axis.text.y = element_blank(), 
                                                                                                   axis.ticks.y = element_blank(), 
                                                                                                   axis.title.y = element_blank())

linearModLOOCV2 <- lm(formula = LOOCV2_Data$pred ~ LOOCV2_Data$obs)
summary(linearModLOOCV2)


grid.arrange(LOOCV_p1,LOOCV_p2,ncol=2)


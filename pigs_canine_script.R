### Script used for the analysis of canine morphology ###
### Cherin et al. submitted ###

library(dplyr);library(ggplot2);library(ggsignif);library(ggrepel)
library(psych);library(gridExtra);library(ggbiplot);library(compositions);library(ggord)
setwd("~/Dropbox/Marco_Crotti/Maiali/Maiali/")
setwd("C:/Users/mario/Dropbox/Marco_Crotti/Maiali/Maiali")

pig <- read.csv("canine_data_final.csv", header=TRUE)

pig2 <- pig[-c(55,75,78:81,82,93,94),] # remove weird c.f. strozzi, 2 juveniles, 4 females



pig3 <- mutate(pig2, ratio_distal_labial = DISTAL/LABIAL,
               ratio_lingual_labial = LINGUAL/LABIAL,
               ratio_lingual_distal = LINGUAL/DISTAL)





### comparison of the labial/posterior ration between canine types

wilcox.test(pig3$ratio_distal_labial ~ pig3$Type)

m1 <- aov(pig3$ratio_distal_labial ~ pig3$Type)
summary(m1)


m2 <- aov(pig3$ratio_lingual_labial ~ pig3$Type)
summary(m2)

distal_labial  <- ggplot() + geom_boxplot(aes(x=Type,y=ratio_distal_labial,fill=Type), data = pig3) + ylim(0.5,2.5) +
  geom_point(size=3, aes(x=Type,y=ratio_distal_labial), data = pig3) +
  scale_fill_brewer(palette = "Dark2") + theme_bw() + 
  geom_point(size=3, aes(x=Type,y=ratio_distal_labial), data = subset(pig3, Locality == "Vallparadis"),shape=21,fill="white",color="black") +
  stat_summary(aes(x=Type,y=ratio_distal_labial,fill=Type), data = pig3,fun.y=mean, shape=23,fill="#666666",color="white", geom="point", 
               shape=18, size=6,show.legend = FALSE) + 
  geom_signif(aes(x=Type,y=ratio_distal_labial), data = pig3,
              comparisons = list(c("scrofic","verrucosic")), y_position = 1.7, map_signif_level = TRUE, textsize = 10, size = 1) +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.position = "none") + 
  geom_label_repel(data = subset(pig3, Locality == "Vallparadis"),aes(x=Type,y=ratio_distal_labial,label = Code),force = 3,ylim = 1,show.legend  = F ) + 
  labs(x="Canine type", y="distal/labial ratio") 
distal_labial



lingual_labial <- ggplot() + geom_boxplot(aes(x=Type,y=ratio_lingual_labial,fill=Type), data = pig3) + ylim(0.5,2.5) +
  geom_point(size=3, aes(x=Type,y=ratio_lingual_labial), data = pig3) +
  scale_fill_brewer(palette = "Dark2") + theme_bw() + 
  geom_point(size=3, aes(x=Type,y=ratio_lingual_labial), data = subset(pig3, Locality == "Vallparadis"),shape=21,fill="white",color="black") +
  stat_summary(aes(x=Type,y=ratio_lingual_labial,fill=Type), data = pig3,fun.y=mean, shape=23,fill="#666666",color="white", geom="point", 
               shape=18, size=6,show.legend = FALSE) + 
  geom_signif(aes(x=Type,y=ratio_lingual_labial), data = pig3,
              comparisons = list(c("scrofic","verrucosic")), y_position = 2.3, map_signif_level = TRUE, textsize = 10, size = 1) +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.position = "none") + 
  geom_label_repel(data = subset(pig3, Locality == "Vallparadis"),aes(x=Type,y=ratio_lingual_labial,label = Code),force = 3,ylim = 0.85,show.legend  = F ) + 
  labs(x="Canine type", y="lingual/labial ratio") 
lingual_labial

grid.arrange(distal_labial, lingual_labial, nrow = 1)




### bivariate plot of lingual/labial and distal/labial between species



figure7 <- ggplot() + theme_bw() + geom_point(aes(x=ratio_distal_labial,y=ratio_lingual_labial,shape=Extant_Fossil,
                          colour=Type), data = subset(pig3, Locality!= "Vallparadis"), size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(aes(x=ratio_distal_labial,y=ratio_lingual_labial,shape=Extant_Fossil), 
             data = subset(pig3, Locality == "Vallparadis"), size=4,shape=24,fill="#D95F02",color="black") +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 20)) + 
  geom_label_repel(aes(x=ratio_distal_labial,y=ratio_lingual_labial,shape=Extant_Fossil,label = Code),
                   data = subset(pig3, Locality == "Vallparadis"),force = 10,xlim = 0.9,show.legend  = F ) + labs(x="distal/labial ratio",
                                                                                                                  y="lingual/labial ratio")
figure7



#######################
#######################   PCA analysis with untransformed data

pca_data <- pig3[,c(1:2,13:16)] ## select data for the pca analysis
pca_data <- na.omit(pca_data) ## remove missing data (1 observation lost)

res.pca <- prcomp(pca_data[,c(3:5)],scale=TRUE)

pc_scores <- as.data.frame(res.pca$x[,1:3])
pca_results <- data.frame(pca_data$Species,pca_data$Extant_Fossil,pca_data$Type,
                          pc_scores)

eig <- (res.pca$sdev)^2  # Eigenvalues
variance <- eig*100/sum(eig) # Variances in percentage
cumvar <- cumsum(variance) # Cumulative variances
eig.pig1 <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
head(eig.pig1)


pca1 <- ggplot(pca_results, aes(x=PC1,y=PC2,colour=pca_data.Type,shape=pca_data.Extant_Fossil)) +
  geom_point(size=4) +
  scale_color_brewer(palette = "Dark2") + theme_bw() +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 10)) +
  labs(x="PC1 80%", y="PC2 18%")

pca1 <- ggbiplot(res.pca, varname.size = 4) + theme_bw() + 
  geom_point(aes(colour= pca_data$Type, shape = pca_data$Extant_Fossil),size=4) + 
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.title = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 10)) + labs(x="PC1 80%",y="PC2 18%")
pca1
  

############# PCA analysis with transformed data

pca_data <- pig3[,c(1:2,3,4,13:16)] ## select data for the pca analysis
pca_data <- na.omit(pca_data)

geo_mean <- geometricmeanRow(pca_data[,5:7]) ### calculate geometric mean
pca_data$GeoMean <- geo_mean

pca_data2 <- mutate(pca_data, adj.LABIAL = LABIAL/GeoMean,
                    adj.LINGUAL = LINGUAL/GeoMean,
                    adj.DISTAL = DISTAL/GeoMean)


res.pca2 <- prcomp(pca_data2[,c(10:12)],scale=TRUE)

pc_scores2 <- as.data.frame(res.pca2$x[,1:3])
pca_results2 <- data.frame(pca_data2$Species,pca_data2$Extant_Fossil,pca_data2$Type,pca_data2$Code,
                          pc_scores2,pca_data2$Locality)

eig2 <- (res.pca2$sdev)^2  # Eigenvalues
variance2 <- eig2*100/sum(eig2) # Variances in percentage
cumvar2 <- cumsum(variance2) # Cumulative variances
eig.pig2 <- data.frame(eig = eig2, variance = variance2,cumvariance = cumvar2)
head(eig.pig2)



ggplot() + geom_point(data = subset(pca_results2, pca_data2.Locality != "Vallparadis"), aes(x=PC1,y=PC2,colour=pca_data2.Type,shape=pca_data2.Extant_Fossil), size = 4) + 
  scale_color_brewer(palette = "Dark2") + theme_bw() + 
  geom_point(data = subset(pca_results2, pca_data2.Locality == "Vallparadis"), 
             aes(x=PC1,y=PC2,colour=pca_data2.Type,shape=pca_data2.Extant_Fossil), size = 4,shape=24,fill="#D95F02",color="black") + 
  geom_label_repel(data = subset(pca_results2, pca_data2.Locality == "Vallparadis"), 
                   aes(x=PC1,y=PC2,colour=pca_data2.Type,shape=pca_data2.Extant_Fossil,label=pca_data2.Code),force = 10,xlim = -1.3,show.legend  = F ) +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 10)) +
  labs(x="PC1 86%", y="PC2 14%") 



ggbiplot(res.pca2)

pca2 <- ggbiplot(res.pca2, varname.size = 4) + theme_bw() + 
  geom_point(aes(colour= pca_data2$Type, shape = pca_data2$Extant_Fossil),size=4)  +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) + 
  theme(legend.title = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 10)) + labs(x="PC1 86%",y="PC2 14%") + 
  xlim(-2.5,2.5)
pca2



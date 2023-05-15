library(RColorBrewer)
library(forcats)
library(ggplot2)
library(ggpubr)
library(Hmisc)

monument_read_abun <- read.csv('monumentome_read_abundance.csv',check.names = F)
monument_metadata <- read.csv('Metadata_Monumentome_2022-04-25_remove_comma.csv',check.names = F)


# read filter (100 reads)

read_posi <- c()
for (i in 1:length(monument_read_abun$name)){
  sample <- monument_read_abun$name[i]
  
  sample_id <- strsplit(sample,'.',fixed = T)[[1]][1]
  bracken <- strsplit(sample,'.',fixed = T)[[1]][2]
  
  if (bracken == 'bracken_num'){
    read_posi <- c(read_posi,i)
    rownames(monument_read_abun)[i] <- sample_id
  }
}

monument_read <- monument_read_abun[read_posi,]
monument_read <- monument_read[,-1]
monument_read[monument_read < 100] <- 0


# Just keep bacterial
taxa_anno_file <- read.csv('monument_taxonomic_level-Adjust_final.csv',check.names = F,header = T)
bacteria_list <- taxa_anno_file$taxa[which(taxa_anno_file$Kingdom == 'Bacteria')]

monument_read <- monument_read[,bacteria_list]


# Calculate RSA
data_normalise <- function(df) {
  return(df/rowSums(df))
}

monument_abun <- monument_read/rowSums(monument_read) # relative abundance data

monument_abun[is.na(monument_abun)] <- 0

monument_abun$sample_id <- rownames(monument_abun)


# abundance data table
monument_data <- merge(monument_metadata,monument_abun,by.x="sample_id",by.y="sample_id")

for (i in 1:length(monument_data$city)){
  if (is.na(monument_data$latitude[i])){
    city <- monument_data$city[i]
    if (city == 'Tel Megiddo'){
      monument_data[i,'latitude'] <- 32.5806
      monument_data[i,'longitude'] <- 35.1795
    }else if (city == 'Ulaanbaatar'){
      monument_data[i,'latitude'] <- 47.8923
      monument_data[i,'longitude'] <- 106.9090
    }else if (city == 'Seoul'){
      monument_data[i,'latitude'] <- 37.5480
      monument_data[i,'longitude'] <- 126.9887
    }
  }
}

# write.csv(monument_data,'monument_data_filter_with_control_b.csv',row.names = F)


## read data table
monument_read$sample_id <- rownames(monument_read)

monument_read_data <- merge(monument_metadata,monument_read,by.x="sample_id",by.y="sample_id")

for (i in 1:length(monument_read_data$city)){
  if (is.na(monument_read_data$latitude[i])){
    city <- monument_read_data$city[i]
    if (city == 'Tel Megiddo'){
      monument_read_data[i,'latitude'] <- 32.5806
      monument_read_data[i,'longitude'] <- 35.1795
    }else if (city == 'Ulaanbaatar'){
      monument_read_data[i,'latitude'] <- 47.8923
      monument_read_data[i,'longitude'] <- 106.9090
    }else if (city == 'Seoul'){
      monument_read_data[i,'latitude'] <- 37.5480
      monument_read_data[i,'longitude'] <- 126.9887
    }
  }
}

# write.csv(monument_read_data,'monument_read_data_with_control_b.csv',row.names = F)

# ------------------------------------------------------------------------------
# Extract Tel Megiddo samples

telm_file <- monument_data[monument_data$city == "Tel Megiddo",]
monument_file_filter <- monument_data[-c(which(rowSums(monument_data[,12:length(colnames(monument_data))]) == 0)),-c(which(colSums(monument_data[,12:length(colnames(monument_data))]) == 0) + 11)]

# remove taxa with zero abundance
taxa_all_list <- colnames(telm_file)[12:length(colnames(telm_file))]
taxa_abun_sum <- colSums(telm_file[12:length(colnames(telm_file))])

posi <- which(taxa_abun_sum!=0)

telm_file <- telm_file[,c(colnames(telm_file)[1:11],taxa_all_list[posi])]


#-------------------------------------------------------------------------------

# taxa count per city

# Create a data form to store the number of times that taxa appears in each city
taxa_all_list <- colnames(monument_file_filter)[12:length(colnames(monument_file_filter))]

monument_file_filter$city <- factor(monument_file_filter$city)

taxa_count_per_city <- data.frame(row.names = c(levels(monument_file_filter$city)))

for (i in 1:length(levels(monument_file_filter$city))){
  for (j in 1: length(taxa_all_list)){
    city <- levels(monument_file_filter$city)[i]
    taxa <- taxa_all_list[j]
    count <- sum(monument_file_filter[monument_file_filter$city == city,][,taxa] != 0)
    taxa_count_per_city[city,taxa] <- count
  }
}

# write.csv(taxa_count_per_city,'monument_taxa_conunt_per_city.csv')

telm_uniq_taxa <- c()
for (i in colnames(taxa_count_per_city)){
    if (sum(taxa_count_per_city[,i])==taxa_count_per_city["Tel Megiddo",i]){
      telm_uniq_taxa <- c(telm_uniq_taxa,i)
    }
}

telm_uniq_taxa_file <- data.frame(telm_uniq_taxa)
# write.csv(telm_uniq_taxa_file,'TelMegiddo_uniq_taxa.csv')



## Calculate average relative abundance in city and global 
monument_abun <- data.frame(row.names = c(levels(monument_file_filter$city)))

for (i in 1:length(levels(monument_file_filter$city))){
  for (j in 1: length(taxa_all_list)){
    city <- levels(monument_file_filter$city)[i]
    taxa <- taxa_all_list[j]  
    abun <- sum(monument_file_filter[monument_file_filter$city == city,][,taxa])/length(monument_file_filter[monument_file_filter$city == city,][,taxa])
    monument_abun[city,taxa] <- abun
  }
}

# Calculate global average relative abundance
global_proportion <- c()

for (i in 1:(length(taxa_all_list))){
  abun <- sum(monument_abun[,i])/length(monument_abun[,i])
  global_proportion <- c(global_proportion,abun)
}

monument_abun <- rbind.data.frame(global_proportion,monument_abun)
rownames(monument_abun)[1] <- 'Global'

# write.csv(monument_abun,'monument_abun_per_city.csv')

taxa_anno_filter <- taxa_anno_file[which(taxa_anno_file$taxa %in% colnames(monument_abun)),]
taxa_anno_filter[taxa_anno_filter$Phylum =='','Phylum'] <- 'None'

# taxa_anno_filter <- droplevels(taxa_anno_filter)
taxa_anno_filter$Phylum <- factor(taxa_anno_filter$Phylum)

phylum_df <- data.frame(row.names = rownames(monument_abun)[2:dim(monument_abun)[1]])

for (i in 2:length(rownames(monument_abun))){
  
  city <- rownames(monument_abun)[i]
  
  for (j in 1:length(levels(taxa_anno_filter$Phylum))){
    phylum <- levels(taxa_anno_filter$Phylum)[j]
    taxa_list <- taxa_anno_filter[taxa_anno_filter$Phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    percentage <- sum(as.numeric(monument_abun[city,taxa_list]))
    phylum_df[city,phylum] <- percentage
  }
  phylum_df[city,'No-annotation'] <- 1 - sum(phylum_df[city,levels(taxa_anno_filter$Phylum)])
}

# write.csv(phylum_df,'20221128monument_phylum_abundance.csv')
# write.csv(taxa_anno_filter,'monument_taxonomy_annotation.csv',row.names = F)


# ------------------------------------------------------------------------------
# Tel Megiddo vs. other monument cites
city_abun <- as.data.frame(t(monument_abun))
city_abun <- city_abun[order(city_abun$`Tel Megiddo`,decreasing = T),]

data <- city_abun[1:25,]
data2 <- t(data)


for (i in colnames(data2)){
  data2[,i] <- (data2[,i] - min(data2[,i]))/(max(data2[,i]) - min(data2[,i]))
}

data2 <- t(data2)
data2 <- as.data.frame(data2)

library(reshape2)
data2$Taxa <- rownames(data2)

data2$Taxa <- fct_inorder(data2$Taxa)
data2$Taxa <- fct_relevel(data2$Taxa,rev(levels(data2$Taxa)))                            

city_name <- colnames(data2)

for (i in 1:length(city_name)){
  city_name[i] <- gsub('_',', ',city_name[i])
}

colnames(data2) <- city_name

data_m <- melt(data2, id.vars=c("Taxa"))
head(data_m)

p = ggplot(data_m, aes(x=variable,y=Taxa))+
  geom_tile(aes(fill=value))+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1))+
  labs(fill="RSA\n(normalized)")+
  xlab('City')
# ggsave('top25_TelM_RSA_normalized.png',p,
#        width = 16,
#        height = 14,
#        units = "cm",
#        dpi = 600)

data$Taxa <- rownames(data)

data$Taxa <- fct_inorder(data$Taxa)
data$Taxa <- fct_relevel(data$Taxa,rev(levels(data$Taxa)))    
colnames(data) <- city_name

data_m <- melt(data, id.vars=c("Taxa"))
head(data_m)

p = ggplot(data_m, aes(x=variable,y=Taxa))+
  geom_tile(aes(fill=value))+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1))+
  labs(fill="RSA")+
  xlab('City')
# ggsave('top25_TelM_RSA_new.png',p,
#        width = 15,
#        height = 14,
#        units = "cm",
#        dpi = 600)


# ------------------------------------------------------------------------------
# Tel Megiddo taxa analysis

# sum(monument_abun['Tel Megiddo',]) 
taxa_anno_file <- read.csv('monument_taxonomic_level-Adjust_final.csv',check.names = F,header = T)


# 1

library(readxl)
telm_anno_file_2 <- read_excel('TelMegiddo - sampling_final.xlsx',sheet = 1, na = 'NA')

telm_file_backup <- telm_file
for (i in 1:length(telm_file$pangea_barcode)){
  s_id <- telm_file$pangea_barcode[i]
  s_id <- strsplit(s_id,'-')[[1]][2]
  s_id <- substr(s_id,2,length(strsplit(s_id,'')[[1]]))
  telm_file[i,'Site'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Site_merged']
  telm_file[i,'Material'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Material']
  telm_file[i,'Access'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Access']
}

telm_file2 <- telm_file[,c(1:11,3430,3431,3432,12:3429)]
rownames(telm_file2) <- telm_file2$sample_id


# Filter of read data
monument_read_file_filter <- monument_read_data[-c(which(rowSums(monument_read_data[,12:length(colnames(monument_read_data))]) == 0)),-c(which(colSums(monument_read_data[,12:length(colnames(monument_read_data))]) == 0)+11)]
telm_read_file <- monument_read_file_filter[monument_read_file_filter$city == "Tel Megiddo",]

# remove taxa with zero abundance
taxa_all_list2 <- colnames(telm_read_file)[12:length(colnames(telm_read_file))]
taxa_read_sum <- colSums(telm_read_file[12:length(colnames(telm_read_file))])

posi2 <- which(taxa_read_sum!=0)

telm_read_file <- telm_read_file[,c(colnames(telm_read_file)[1:11],taxa_all_list2[posi2])]


for (i in 1:length(telm_read_file$pangea_barcode)){
  s_id <- telm_read_file$pangea_barcode[i]
  s_id <- strsplit(s_id,'-')[[1]][2]
  s_id <- substr(s_id,2,length(strsplit(s_id,'')[[1]]))
  telm_read_file[i,'Site'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Site_merged']
  telm_read_file[i,'Material'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Material']
  telm_read_file[i,'Access'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Access']
}

telm_read_file2 <- telm_read_file[,c(1:11,3430,3431,3432,12:3429)]
rownames(telm_read_file2) <- telm_read_file2$sample_id

telm_read_data <- telm_read_file2[,15:dim(telm_read_file2)[2]] # 44 samples
rownames(telm_read_data) <- telm_read_file2$sample_id

telm_taxa_list <- colnames(telm_file2)[15:dim(telm_file2)[2]]



## diversity analysis ----------------------------------------------------------

# a diversity
library(vegan)

alpha_diversity <- function(otu){
  Read_number <- rowSums(otu)
  Observed_species <- estimateR(otu)[1, ] # Richness index
  Chao1  <- estimateR(otu)[2, ] # Chao 1 index
  ACE  <- estimateR(otu)[4, ] # ACE index
  Shannon <- diversity(otu,'shannon') # shannon index
  Simpson <- diversity(otu,'simpson') # Gini-Simpson index
  Goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)
  output <- data.frame(Read_number, Observed_species, Chao1, ACE, Shannon, Simpson, Goods_coverage)
  return(output)
}

telm_alpha_diversity <- alpha_diversity(telm_read_data)

telm_alpha_diversity <- telm_alpha_diversity[-c(which(telm_alpha_diversity$Read_number<1000)),]
telm_alpha_diversity <- telm_alpha_diversity[-c(which(telm_alpha_diversity$Observed_species<10)),]

for (i in 1:length(rownames(telm_alpha_diversity))){
  id <- rownames(telm_alpha_diversity)[i]
  telm_alpha_diversity[i,'Site'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Site']
  telm_alpha_diversity[i,'Material'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Material']
  telm_alpha_diversity[i,'Access'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Access']
}
# write.csv(telm_alpha_diversity,'telm_alpha_diversity_final.csv')

sample_id_list <- rownames(telm_alpha_diversity)

telm_file2 <- telm_file2[sample_id_list,]
telm_read_file2 <- telm_read_file2[sample_id_list,]

# write.csv(telm_file2,'telm_abun_40samples.csv')
# write.csv(telm_read_file2,'telm_reads_40samples.csv')

telm_read_data <- telm_read_file2[,15:dim(telm_read_file2)[2]] # just contain read data
telm_abun_data <- telm_file2[,15:dim(telm_file2)[2]] # just contain abundance data


## Access-Material (stone access sample) ---------------------------------------
# just run when need to control material when analysis the access
telm_file2 <- telm_file2[telm_file2$Material == 'stone',]

telm_read_file2 <- telm_read_file2[telm_read_file2$Material == 'stone',]

telm_read_data <- telm_read_file2[,15:dim(telm_read_file2)[2]] # just contain read data
telm_abun_data <- telm_file2[,15:dim(telm_file2)[2]] # just contain abundance data

# telm_alpha_diversity <- read.csv('telm_alpha_diversity_final.csv',check.names = T,header = T)
telm_alpha_diversity <- telm_alpha_diversity[telm_alpha_diversity$X %in% telm_file2$sample_id,]
rownames(telm_alpha_diversity) <- telm_alpha_diversity$X
telm_alpha_diversity <- telm_alpha_diversity[,-1]

##------------------------------------------------------------------------------


## telm taxonomic annotation

anno_num <- which(colnames(telm_file2) %in% taxa_anno_file$taxa)
sum(colSums(telm_file2[,which(colnames(telm_file2) %in% taxa_anno_file$taxa)])/dim(telm_file2)[1]) #  0.7695763

taxa_avg_abun <- colSums(telm_file2[,15:dim(telm_file2)[2]])/dim(telm_file2)[1]
sum(taxa_avg_abun) 
taxa_avg_abun <- as.data.frame(taxa_avg_abun)


taxa_anno_file_telm <- taxa_anno_file[which(taxa_anno_file$taxa %in% telm_taxa_list),]

taxa_anno_file_telm$Phylum[which(taxa_anno_file_telm$Phylum=="")] <- 'None'
taxa_anno_file_telm$Order[which(taxa_anno_file_telm$Order=="")] <- 'None'

taxa_anno_file_telm$Phylum <- factor(taxa_anno_file_telm$Phylum)
taxa_anno_file_telm$Order <- factor(taxa_anno_file_telm$Order)

# write.csv(taxa_anno_file_telm,'taxa_anno_file_telm.csv')

## group analysis ____________________________________________________________


variable <- 'Site'
# variable <- 'Material'
# variable <- 'Access'

telm_file2[,variable] <- factor(telm_file2[,variable])

for (i in 1:length(rownames(telm_alpha_diversity))){
  id <- rownames(telm_alpha_diversity)[i]
  telm_alpha_diversity[i,'group'] <- telm_read_file2[telm_read_file2$sample_id == id, variable]
}

sample_n <- table(telm_alpha_diversity$group)

mycompare <- list()
for (i in 1:(length(levels(factor(telm_alpha_diversity$group)))-1)){
  k1 <- levels(factor(telm_alpha_diversity$group))[i]
  for (j in (i+1):length(levels(factor(telm_alpha_diversity$group)))){
    k2 <- levels(factor(telm_alpha_diversity$group))[j]
    k <- c(k1,k2)
    mycompare <- c(mycompare,list(k))
  }
}

mycompare <- list(c("Two-chambered City Gate","Entrance"),c("Two-chambered City Gate","Area K"),c("Two-chambered City Gate","Area S" ),
                  c("Entrance",'Area K'))
mycompare <- list(c("Two-chambered City Gate","Area J"),c("Two-chambered City Gate","Area K"),c("Two-chambered City Gate","Area S" ),
                  c("Entrance",'Area K'),
                  c("Area J",'Area K'))
mycompare <- list(c("Two-chambered City Gate","Area K"),c("Two-chambered City Gate","Area S" ),
                  c("Entrance",'Area S'),
                  c("Area J",'Area K'),c("Area J",'Area S'))

mycompare <- list(c('cement','soil'),c('cement','stone'),
                  c('soil','stone'))

# mycompare <- list(c('non-public','public'))

kruskal.test(Observed_species~group,data=telm_alpha_diversity)

p = ggplot(telm_alpha_diversity,aes(x=group,y=Observed_species,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  labs(x='',y='Observed species')+
  ylim(0,2600)+#ylim(0,3500)+#ylim(0,4200)+# ylim(0,4600)+
  geom_signif(comparisons = mycompare,
              step_increase = 0.12,
              # color = 'black',
              map_signif_level = T,
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))


# ggsave(paste0('telm_',variable,'_obseved_species_new.png'),
#        p,
#        width = 6,
#        height = 8.5,
#        units = "cm",
#        dpi = 600)


kruskal.test(Shannon~group,data=telm_alpha_diversity)
p = ggplot(telm_alpha_diversity,aes(x=group,y=Shannon,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  ylim(2,7.2)+# ylim(2,8.3)+# ylim(2,9.5)+#ylim(1,9.1)+#ylim(1,12.9)+
  labs(x='',y='Shannon ')+
  geom_signif(comparisons = mycompare,
              step_increase = 0.13,
              map_signif_level = T,
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))

# ggsave(paste0('telm_',variable,'_shannon_new.png'),
#        p,
#        width = 6,
#        height = 8.5,
#        units = "cm",
#        dpi = 600)

kruskal.test(Simpson~group,data=telm_alpha_diversity)
p = ggplot(telm_alpha_diversity,aes(x=group,y=Simpson,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  ylim(0.77,1.02)+# ylim(0.77,1.08)+# ylim(0.77,1.14)+#ylim(0.4,1.2)+#ylim(0.4,1.62)+
  labs(x='',y='Simpson')+
  geom_signif(comparisons = mycompare,
              step_increase = 0.12,
              map_signif_level = T,
              # color = 'black',
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))
# 
# ggsave(paste0('telm_',variable,'_simpson_new.png'),
#        p,
#        width = 6,
#        height = 8.5,
#        units = "cm",
#        dpi = 600)


# beta diversity -------------------------------------------------------

library(amplicon)

sample_group <- as.data.frame(telm_alpha_diversity$group)
colnames(sample_group) <- 'Group'
rownames(sample_group) <- rownames(telm_alpha_diversity)


sample_n <- table(telm_alpha_diversity$group)


# PCoA
pcoa_result <- BetaDiv(otu=t(telm_read_data),map=sample_group,group='Group',
                       dist='bray',method='PCoA',Micromet='adonis')

pcoa_points <- pcoa_result[[2]]
pcoa_method = 'PCoA'
pcoa_adonis <- pcoa_result[[5]]

library(dplyr)
library(tidyr)
hull_data <- 
  pcoa_points %>%
  drop_na() %>%
  group_by(Group) %>% 
  slice(chull(x, y))


ggplot(pcoa_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste0(pcoa_method," 1 (36.82%)"),
       y=paste0(pcoa_method," 2 (12.28%)"),
       # title=pcoa_adonis)+
       title=expression(paste('adonis:',italic(r),'=0.028,',italic(p),'=0.324'))) +
  # stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))+
  geom_polygon(data = hull_data,
               aes(fill = Group,
                   colour = Group),
               alpha = 0.3,
               show.legend = FALSE)+
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-110,hjust=1.85,size = 10),
        text=element_text(size=10),
        axis.text =element_text(colour="black",size=9),
        # axis.text.y=element_text(colour="black",size=7),
        legend.text=element_text(size=9),
        legend.background = element_rect(
          size=0.3, linetype="solid",
          colour ="black"),
        legend.position= 'right')+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(pcoa_points$Group))),' (',sample_n,')')))
# ggsave(paste0('telm_',variable,'_PCoA_Bray_Crutis.png'),
#        p,
#        width = 16,
#        height = 10,
#        units = "cm",
#        dpi = 600)


##Access----

p = ggplot(pcoa_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste0(pcoa_method," 1 (32.01%)"),
       y=paste0(pcoa_method," 2 (15.78%)"),
       # title=pcoa_adonis)+
       title=expression(paste('adonis:',italic(r),'=0.031,',italic(p),'=0.645'))) +
  # stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))+
  geom_polygon(data = hull_data,
               aes(fill = Group,
                   colour = Group),
               alpha = 0.3,
               show.legend = FALSE)+
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-120,hjust=0.95,size = 10),
        text=element_text(size=10),
        axis.text =element_text(colour="black",size=9),
        # axis.text.y=element_text(colour="black",size=7),
        legend.text=element_text(size=9),
        legend.background = element_rect(
          size=0.3, linetype="solid",
          colour ="black"),
        legend.position= "top")+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(pcoa_points$Group))),' (',sample_n,')')))
# ggsave(paste0('telm_',variable,'_PCoA_Bray_Crutis.png'),
#        p,
#        width = 10,
#        height = 11,
#        units = "cm",
#        dpi = 600)

##----

# NMDS

nmds_result <- BetaDiv(otu=t(telm_read_data),map=sample_group,group='Group',
                       dist='bray',method='NMDS',Micromet='adonis')

nmds_points <- nmds_result[[2]]
nmds_method = 'NMDS'
nmds_adonis <- nmds_result[[5]]

hull_data <- 
  nmds_points %>%
  drop_na() %>%
  group_by(Group) %>% 
  slice(chull(x, y))

p = ggplot(nmds_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste(nmds_method,"1", sep=""),
       y=paste(nmds_method,"2",sep=""),
       title=expression(paste('adonis:',italic(r),'=0.288,',italic(p),'=0.002')))+
  annotate("text", x= 0.4, y= 1.1,label='Stress=0.13',size=3.5 )+
  # stat_ellipse( linetype=2,level=0.68,aes(group  =Group, colour= Group))+
  # Plot_ConvexHull(nmds_points, col_g)+
  geom_polygon(data = hull_data,
               aes(fill = Group,
                   colour = Group),
               alpha = 0.3,
               show.legend = FALSE)+
  
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-110,hjust=1.85,size = 10),
        text=element_text(size=10),
        # axis.title.y =element_text(size=7,face="bold",colour="black"),
        # axis.title.x =element_text(size=7,face="bold",colour="black"),
        # axis.text=element_text(size=7,face="bold"),
        axis.text =element_text(colour="black",size=9),
        # axis.text.y=element_text(colour="black",size=7),
        legend.text=element_text(size=9),
        legend.background = element_rect(
          size=0.3, linetype="solid",
          colour ="black"),
        legend.position= 'right')+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(nmds_points$Group))),' (',sample_n,')')))

# ggsave(paste0('telm_',variable,'_NMDS_Bray_Crutis.png'),
#        p,
#        width = 16,
#        height = 10,
#        units = "cm",
#        dpi = 600)

## access -----

p = ggplot(nmds_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste(nmds_method,"1", sep=""),
       y=paste(nmds_method,"2",sep=""),
       title=expression(paste('adonis:',italic(r),'=0.028,',italic(p),'=0.326')))+
  annotate("text", x= 0.4, y= 1,label='Stress=0.13',size=3.5 )+
  # stat_ellipse( linetype=2,level=0.68,aes(group  =Group, colour= Group))+
  # Plot_ConvexHull(nmds_points, col_g)+
  geom_polygon(data = hull_data,
               aes(fill = Group,
                   colour = Group),
               alpha = 0.3,
               show.legend = FALSE)+
  
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-28,hjust=0.01,size = 10),
        text=element_text(size=10),
        # axis.title.y =element_text(size=7,face="bold",colour="black"),
        # axis.title.x =element_text(size=7,face="bold",colour="black"),
        # axis.text=element_text(size=7,face="bold"),
        axis.text =element_text(colour="black",size=9),
        # axis.text.y=element_text(colour="black",size=7),
        legend.text=element_text(size=9),
        legend.background = element_rect(
          size=0.3, linetype="solid",
          colour ="black"),
        legend.position= 'top')+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(nmds_points$Group))),' (',sample_n,')')))

# ggsave(paste0('telm_',variable,'_NMDS_Bray_Crutis.png'),
#        p,
#        width = 10,
#        height = 11,
#        units = "cm",
#        dpi = 600)
### ------


## std for adonis in PCoA and NMDS (1000 times)

##delate ----------------------
# variable <- 'Site'
variable <- 'Material'
# variable <- 'Access'


# telm_alpha_diversity <- read.csv("telm_alpha_diversity_final.csv",header = T)
rownames(telm_alpha_diversity) <- telm_alpha_diversity[,1]
telm_alpha_diversity <- telm_alpha_diversity[,-1]

telm_file2[,variable] <- factor(telm_file2[,variable])

for (i in 1:length(rownames(telm_alpha_diversity))){
  id <- rownames(telm_alpha_diversity)[i]
  telm_alpha_diversity[i,'group'] <- telm_read_file2[telm_read_file2$sample_id == id, variable]
}

##delate ----------------------

r_1000 <- c()
p_1000 <- c()

for (i in 1:1000){
  telm_alpha_diversity2 <- telm_alpha_diversity[sample(dim(telm_alpha_diversity)[1],20.8),] # 80%
  telm_read_data2 <- telm_read_data[rownames(telm_alpha_diversity2),]
  
  sample_group <- as.data.frame(telm_alpha_diversity2$group)
  colnames(sample_group) <- 'Group'
  rownames(sample_group) <- rownames(telm_alpha_diversity2)
  
  
  sample_n <- table(telm_alpha_diversity2$group)
  
  pcoa_result <- BetaDiv(otu=t(telm_read_data2),map=sample_group,group='Group',
                         dist='bray',method='PCoA',Micromet='adonis')
  
  pcoa_adonis <- pcoa_result[[5]]
  
  r_1000 <- c(r_1000, as.numeric(strsplit(pcoa_adonis,' ')[[1]][2]))
  p_1000 <- c(p_1000, as.numeric(strsplit(pcoa_adonis,' ')[[1]][4]))
  
  print(i)
}

avg_r <- mean(r_1000) 
avg_p <- mean(p_1000) 

sum(r_1000>0)
sum(p_1000<0.05) 

r_1000_df <- as.data.frame(r_1000)

# p = ggplot(r_1000_df, aes(x = r_1000)) +
#   geom_density(alpha = 0.3)+
#   xlab("adonis:r")+
#   theme_bw()
# # ggsave('material_PCoA_adonis_r.png',p,dpi = 600)


p_1000_df <- as.data.frame(p_1000)

# p = ggplot(p_1000_df, aes(x = p_1000)) +
#   geom_density(alpha = 0.3)+
#   xlab("adonis:p")+
#   theme_bw()
# # ggsave('material_PCoA_adonis_p.png',p,dpi = 600)

r_p_adonis_1000 <- cbind.data.frame(r_1000,p_1000)
# write.csv(r_p_adonis_1000,'PCoA_1000_adonis_access-stone.csv')

#####

# sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))} 
# sd.p(r_1000)
# sd(p_1000)
# mean(r_1000)

###### repeat 1000 times NMDS

r_1000 <- c()
p_1000 <- c()

for (i in 1:1000){
  telm_alpha_diversity2 <- telm_alpha_diversity[sample(dim(telm_alpha_diversity)[1],20.8),] # 80%
  telm_read_data2 <- telm_read_data[rownames(telm_alpha_diversity2),]
  
  sample_group <- as.data.frame(telm_alpha_diversity2$group)
  colnames(sample_group) <- 'Group'
  rownames(sample_group) <- rownames(telm_alpha_diversity2)
  
  sample_n <- table(telm_alpha_diversity2$group)
  
  nmds_result <- BetaDiv(otu=t(telm_read_data2),map=sample_group,group='Group',
                         dist='bray',method='NMDS',Micromet='adonis')
  nmds_adonis <- nmds_result[[5]]
  
  r_1000 <- c(r_1000, as.numeric(strsplit(nmds_adonis,' ')[[1]][2]))
  p_1000 <- c(p_1000, as.numeric(strsplit(nmds_adonis,' ')[[1]][4]))
}

avg_r <- mean(r_1000) 
avg_p <- mean(p_1000) 

sum(r_1000>0) 
sum(p_1000<0.05) 

# r_1000_df <- as.data.frame(r_1000)

# p = ggplot(r_1000_df, aes(x = r_1000)) +
#   geom_density(alpha = 0.3)+
#   xlab("adonis:r")+
#   theme_bw()
# ggsave('material_PCoA_adonis_r.png',p,dpi = 600)


# p_1000_df <- as.data.frame(p_1000)
# 
# p = ggplot(p_1000_df, aes(x = p_1000)) +
#   geom_density(alpha = 0.3)+
#   xlab("adonis:p")+
#   theme_bw()
# ggsave('material_PCoA_adonis_p.png',p,dpi = 600)

r_p_adonis_1000 <- cbind.data.frame(r_1000,p_1000)
# write.csv(r_p_adonis_1000,'NMDS_1000_adonis_access-stone.csv')


# analysis microbiome composition   ----------------------------------

telm_file2[,variable] <- factor(telm_file2[,variable])
abun_df <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  for (j in 1: length(telm_taxa_list)){
    v <- levels(telm_file2[,variable])[i]
    taxa <- telm_taxa_list[j]  
    abun <- sum(telm_file2[telm_file2[,variable] == v,][,taxa])/length(telm_file2[telm_file2[,variable] == v,][,taxa])
    abun_df[v,taxa] <- abun
  }
}

sum(abun_df[1,])

# write.csv(abun_df,'site/Tel_Megiddo_site_taxa_abun.csv')
# write.csv(abun_df,'material/Tel_Megiddo_material_taxa_abun.csv')
# write.csv(abun_df,'access/Tel_Megiddo_access_taxa_abun.csv')
# write.csv(abun_df,'access_stone/Tel_Megiddo_access_stone_taxa_abun.csv')

rownames(abun_df) <- abun_df[,1]
abun_df <- abun_df[,-1]


## Phylum relative abundance
phylum_df <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  
  city <- levels(telm_file2[,variable])[i]
  
  for (j in 1:length(levels(taxa_anno_file_telm$Phylum))){
    phylum <- levels(taxa_anno_file_telm$Phylum)[j]
    taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    percentage <- sum(as.numeric(abun_df[city,taxa_list]))
    phylum_df[city,phylum] <- percentage
  }
}

# write.csv(phylum_df,'site/Tel_Megiddo_site_phylum_abundance.csv')
# write.csv(phylum_df,'material/Tel_Megiddo_material_phylum_abundance.csv')
# write.csv(phylum_df,'access/Tel_Megiddo_access_phylum_abundance.csv')
# write.csv(phylum_df,'access_stone/Tel_Megiddo_access_stone_phylum_abundance.csv')

rownames(phylum_df) <- phylum_df[,1]
phylum_df <- phylum_df[,-1]

## order RSA table

order_df <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  
  city <- levels(telm_file2[,variable])[i]
  
  for (j in 1:length(levels(taxa_anno_file_telm$Order))){
    Order <- levels(taxa_anno_file_telm$Order)[j]
    taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Order == Order,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    percentage <- sum(as.numeric(abun_df[city,taxa_list]))
    order_df[city,Order] <- percentage
  }
}

# write.csv(order_df,'site/Tel_Megiddo_site_order_RSA.csv')
# write.csv(order_df,'material/Tel_Megiddo_material_order_abundance.csv')
# write.csv(order_df,'access/Tel_Megiddo_access_order_abundance.csv')
# write.csv(order_df,'access_stone/Tel_Megiddo_access_stone_order_abundance.csv')



# unique_phylum
# phylum_df2 <- phylum_df
# phylum_df2[phylum_df2 > 0] <- 1
# phylum_occupancy <- colSums(phylum_df2)
# phylum_occupancy <- phylum_occupancy[order(phylum_occupancy)]
# 
# uniq_p_city <- phylum_df2[,c(names(phylum_occupancy[1:8]))]

phylum_average_abun <- colSums(phylum_df)/dim(phylum_df)[1]
phylum_average_abun <- phylum_average_abun[order(phylum_average_abun,decreasing = T)]

# phylum name list
# select_phylum <- names(phylum_average_abun[1:8])

select_phylum_df <- phylum_df[,select_phylum]

for (i in 1:dim(select_phylum_df)[1]){
  select_phylum_df[i,'Other'] <- sum(phylum_df[i,-which(colnames(phylum_df) %in% select_phylum)])
}

select_phylum_df <- select_phylum_df[order(select_phylum_df$Actinobacteria),]
# write.csv(select_phylum_df,'material/select_phylum_df_plot.csv')

sample_n <- table(telm_file2[,variable])
sample_n2 <- sample_n[rownames(select_phylum_df)]


palette <- colorRampPalette(brewer.pal(10, "Set3"))(10)[1:9] # 9 7

# png("TelM_material_phylum_abundance.png", width = 6,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(6.5,0,0,8.7))
bp <- barplot(t(select_phylum_df*100), col=palette,
              # names.arg= ,
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F, axes = T, space =0,
              ylim = c(0,100))
mtext(text = paste0(capitalize(rownames(select_phylum_df)),' (',sample_n2,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Relative abundance % (Phylum level)",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-0.7,0), colnames(select_phylum_df),  fill = palette , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

## access-----------------
# png("TelM_access_phylum_abundance.png", width = 4,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(3.5,0,0,9))
bp <- barplot(t(select_phylum_df*100), col=palette,
              # names.arg= ,
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F, axes = T, space =0,
              ylim = c(0,100))
mtext(text = paste0(capitalize(rownames(select_phylum_df)),' (',sample_n2,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Relative abundance % (Phylum level)",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-2.3,0), colnames(select_phylum_df),  fill = palette , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()
##

########

## taxa number in each phylum appeared in each site

phylum_taxa <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  
  city <- levels(telm_file2[,variable])[i]
  
  for (j in 1:length(levels(taxa_anno_file_telm$Phylum))){
    phylum <- levels(taxa_anno_file_telm$Phylum)[j]
    taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    num <- sum(abun_df[city,taxa_list]>0)
    
    phylum_taxa[city,phylum] <- num
  }
  
  phylum_taxa[city,'No-annotation'] <- sum(abun_df[city,]>0) - sum(phylum_taxa[city,levels(taxa_anno_file_telm$Phylum)])
  phylum_taxa[city,'Total_taxa_number'] <- sum(abun_df[city,]>0)
}

# write.csv(phylum_taxa,'site/Tel_Megiddo_site_phylum_count.csv')
# write.csv(phylum_taxa,'material/Tel_Megiddo_material_phylum_count.csv')
# write.csv(phylum_taxa,'access/Tel_Megiddo_access_phylum_count.csv')
# write.csv(phylum_taxa,'access_stone/Tel_Megiddo_access_stone_phylum_count.csv')



# Pathogen analysis ------------------------------------------------------------

taxa_anno_file_telm2 <- taxa_anno_file_telm
taxa_anno_file_telm2[is.na(taxa_anno_file_telm2)] <- 0

animal_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 1,'taxa']
plant_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 2,'taxa']
both_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 3,'taxa']

total_pathogen <- c(animal_pathogen,plant_pathogen,both_pathogen)
abun_df_pathogen <- abun_df[,total_pathogen]

pathogen_abun <- data.frame(row.names = rownames(abun_df))

for (i in 1:dim(pathogen_abun)[1]){
  pathogen_abun[i,'animal_pathogen_n'] <- sum(abun_df[i,animal_pathogen] > 0)
  pathogen_abun[i,'animal_pathogen_p'] <- sum(as.numeric(abun_df[i,animal_pathogen]))
  pathogen_abun[i,'plant_pathogen_n'] <- sum(abun_df[i,plant_pathogen] > 0)
  pathogen_abun[i,'plant_pathogen_p'] <- sum(as.numeric(abun_df[i,plant_pathogen]))
  pathogen_abun[i,'both_pathogen_n'] <- sum(abun_df[i,both_pathogen] > 0)
  pathogen_abun[i,'both_pathogen_p'] <- sum(as.numeric(abun_df[i,both_pathogen]))
}

abun_df_pathogen <- cbind(abun_df_pathogen,pathogen_abun)
# write.csv(abun_df_pathogen,'site_pathogen_RSA_summary.csv')


taxa_p <- data.frame(row.names = 'type')
for (i in total_pathogen){
  taxa_p['type',i] <- taxa_anno_file_telm2[taxa_anno_file_telm2$taxa == i,'pathogen_type']
}
# write.csv(taxa_p,'site/site_pathogen_RSA_summary2.csv')
# write.csv(taxa_p,'material/material_pathogen_RSA_summary2.csv')
# write.csv(taxa_p,'access/access_pathogen_RSA_summary2.csv')
# write.csv(taxa_p,'access_stone/access_stone_pathogen_RSA_summary2.csv')



bar_df2 <- pathogen_abun[,c('plant_pathogen_p','both_pathogen_p','animal_pathogen_p')]

sample_n2 <- sample_n[rownames(bar_df2)]

png("pathogen_abundance_material.png", width = 6,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(6.5,0,0,8.2))
bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D"), 
              ylim = c(0,25),
              names.arg= paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')') ,
              args.legend = list(x = "right", inset=c(-0.5,0)), las =2, cex.names=.7,
              ylab = "",axisnames = F, axes = T)
mtext(text = paste0(capitalize((rownames(pathogen_abun))),' (',sample_n2,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Pathogen relative abundance %",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-0.65,0),
       c('Animal Pathogen','Dual-kingdom\nPathogen','Plant Pathogen'),
       fill = rev(c("#619CFF","#00BA38","#F8766D")) , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

##access

# access
png("pathogen_abundance_access.png", width = 4,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(3.5,0,0,8.5))
bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D"), 
              ylim = c(0,15),
              names.arg= paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')') ,
              args.legend = list(x = "right", inset=c(-0.5,0)), las =2, cex.names=.7,
              ylab = "",axisnames = F, axes = T)
mtext(text = paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Pathogen relative abundance %",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-1.95,0),
       c('Animal Pathogen','Dual-kingdom\nPathogen','Plant Pathogen'),
       fill = rev(c("#619CFF","#00BA38","#F8766D")) , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
##


# Whole dataset sign -----------------------------------------------------------

## total pathogen abundance

taxa_avg_abun <- data.frame(colSums(telm_abun_data)/dim(telm_abun_data)[1])

pathogen_abun <- data.frame(row.names = 'total')

pathogen_abun[1,'animal_pathogen_n'] <- sum(taxa_avg_abun[animal_pathogen,1] > 0)
pathogen_abun[1,'animal_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[animal_pathogen,1]))
pathogen_abun[1,'plant_pathogen_n'] <- sum(taxa_avg_abun[plant_pathogen,1] > 0)
pathogen_abun[1,'plant_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[plant_pathogen,1]))
pathogen_abun[1,'both_pathogen_n'] <- sum(taxa_avg_abun[both_pathogen,1] > 0)
pathogen_abun[1,'both_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[both_pathogen,1]))

for (i in 1:dim(taxa_avg_abun)[1]){
  taxa <- rownames(taxa_avg_abun)[i]
  if (taxa %in% taxa_anno_file_telm2$taxa){
    taxa_avg_abun[i,'pathogen_type'] <- taxa_anno_file_telm2[taxa_anno_file_telm2$taxa==taxa,'pathogen_type']
  }else{
    taxa_avg_abun[i,'pathogen_type'] <- 'No_annotation'
  }
}

# write.csv(taxa_avg_abun,'Telm_taxa_avg_abundance_and_pathogen_annotation.csv')



## alpha diversity data for all samples
alpha_file <- telm_alpha_diversity
alpha_file <- cbind.data.frame(rownames(alpha_file),alpha_file)
colnames(alpha_file)[1] <- 'id'

alpha_file <- alpha_file[order(alpha_file$Observed_species),]

alpha_file$id <- factor(alpha_file$id)
alpha_file$id <- fct_inorder(alpha_file$id)
alpha_file$id <- fct_relevel(alpha_file$id,rev(levels(alpha_file$id))) 

alpha_file$Simpson

p = ggplot(alpha_file,aes(x=id, y=Simpson))+
  geom_bar(stat="identity",fill = 'lightblue')+
  ylab('Simpson')+
  xlab('Sample ID')+
  theme_bw()+
  # theme(axis.text.x = element_text( vjust = 0.5, hjust = 1, angle = 90))+
  coord_flip()
# ggsave('Simpson - all.png',
#        p,
#        width = 10,
#        height = 15,
#        units = "cm",
#        dpi = 600)


## taxa distribution in area ---------------------------------------------------
# Create a data form to store the number of times that taxa appears in each city

bar_df <- data.frame(row.names = c(levels(telm_file2$Site)))

telm_taxa_list 

for (i in 1:length(levels(telm_file2$Site))){
  for (j in 1: length(telm_taxa_list)){

    city <- levels(telm_file2$Site)[i]
    taxa <- telm_taxa_list[j]

    count <- sum(telm_file2[telm_file2$Site == city,][,taxa] != 0)

    bar_df[city,taxa] <- count
  }
}

Overall <- c()
for (i in telm_taxa_list){
  count_sum <- sum(bar_df[,i])
  Overall <- c(Overall,count_sum)
}

bar_df <- rbind(bar_df,Overall)
row.names(bar_df) <- c(levels((telm_file2$Site)),"Overall")

# write.csv(bar_df,"taxa_count_in_site.csv",row.names = T)

bar_df_t <- data.frame(t(bar_df))
bar_df_t <- data.frame(taxa_list = colnames(bar_df),bar_df_t)

for (i in 1:dim(bar_df_t)[1]){
  bar_df_t[i,'Percentage_overall'] <- round(bar_df_t[i,'Overall']/40*100,2)
}

for ( i in 1:dim(bar_df_t)[1]){
  bar_df_t[i,'Overall_city'] <- sum(bar_df_t[i,][,2:11] != 0)
  bar_df_t[i,'Percentage_city'] <- round(bar_df_t[i,'Overall_city']/10*100,2)
}


count_sample2 <- bar_df_t[order(bar_df_t$Percentage_city,decreasing=T),]

# write.csv(count_sample2,'count_sample_site_with_occupancy.csv')

percentage2 <- matrix(ncol=2,nrow=10)

percentage2[1,1] <- '100%'
percentage2[2,1] <- '90%'
percentage2[3,1] <- '80%'
percentage2[4,1] <- '70%'
percentage2[5,1] <- '60%'
percentage2[6,1] <- '50%'
percentage2[7,1] <- '40%'
percentage2[8,1] <- '30%'
percentage2[9,1] <- '20%'
percentage2[10,1] <- '10%\n1 Area'

percentage2[1,2] <- sum(count_sample2[,'Overall_city'] == 10)
percentage2[2,2] <- sum(count_sample2[,'Overall_city'] == 9)
percentage2[3,2] <- sum(count_sample2[,'Overall_city'] == 8)
percentage2[4,2] <- sum(count_sample2[,'Overall_city'] == 7)
percentage2[5,2] <- sum(count_sample2[,'Overall_city'] == 6)
percentage2[6,2] <- sum(count_sample2[,'Overall_city'] == 5)
percentage2[7,2] <- sum(count_sample2[,'Overall_city'] == 4)
percentage2[8,2] <- sum(count_sample2[,'Overall_city'] == 3)
percentage2[9,2] <- sum(count_sample2[,'Overall_city'] == 2)
percentage2[10,2] <- sum(count_sample2[,'Overall_city'] == 1)

percentage2 <- data.frame(percentage2)
percentage2$X2 <- as.numeric(percentage2$X2)
percentage2$X1 <- fct_inorder(percentage2$X1)

p = ggplot(percentage2,aes(x=X1,y=X2))+
  geom_bar(stat='identity',fill = "lightblue")+
  ylim(0,610)+
  # geom_text(aes(label=X2,vjust = 0),
  #           position=position_dodge(width=1))+
  geom_text(aes(label= c(paste0(X2,' (',round(X2/3418*100),'%)')),vjust = -0.5,y = c(X2) + 0.5),position=position_dodge(0))+
  ggtitle('Taxa distribution over areas (3418 taxa, 10 areas)')+
  labs(x = 'Percentage of areas (%)', y = 'Number of taxa')+
  theme_bw()+
  theme(axis.title =  element_text(size=10.5,face = "bold"),
        axis.text =   element_text(size=10.5)) # adjust size of label

# ggsave('taxa site occupancy - TelM.png',
#        p,
#        width = 19,
#        height = 10,
#        units = "cm",
#        dpi = 600)



### compare within monumentome -------------------------------------------------

#### overlap in tel megiddo and other city
metasub <- monument_abun

taxa_overlap <- c()
for (i in telm_taxa_list){
    if (sum(monument_abun[2:13,i])!=monument_abun["Tel Megiddo",i]){
      taxa_overlap <- c(taxa_overlap,i)
    }
} # 3393 taxa overlap

## abundance comparison
taxa_overlap_access <- data.frame(row.names = c('Public','Non-public'))
for (i in taxa_overlap){
  taxa_overlap_access['Public',i] <- mean(telm_file2[telm_file2$Access == 'public',i]) # 18 samples
  taxa_overlap_access['Non-public',i] <- mean(telm_file2[telm_file2$Access == 'non-public',i]) # 22 samples
}
taxa_overlap_access['Public','Sum'] <- sum(taxa_overlap_access['Public',taxa_overlap])
taxa_overlap_access['Non-public','Sum'] <- sum(taxa_overlap_access['Non-public',taxa_overlap])
taxa_overlap_access$Sum # public:0.9999619   non-public:0.9999730

## count comparison
sum(taxa_overlap_access['Public',taxa_overlap] != 0) # 2606
sum(taxa_overlap_access['Non-public',taxa_overlap] != 0) # 3298



# overlap in metasub and tel megiddo -------------------------------------------
metasub <- read.csv('Metasub/MetaSub_taxa_city_proportion.csv')
rownames(metasub) <- metasub[,1]
metasub <- metasub[,-1]

m_Overall <- c()
for (i in colnames(metasub)){
  count_sum <- sum(metasub[,i])
  m_Overall <- c(m_Overall,count_sum)
}

metasub <- rbind.data.frame(metasub,m_Overall)
rownames(metasub)[54] <- 'Overall'

which(m_Overall==0)
metasub <- metasub[,-which(m_Overall==0)]

taxa_metasub <- colnames(metasub)
taxa_telm <- make.names(telm_taxa_list)

taxa_overlap <- c()
for (i in taxa_telm){
  if(i %in% taxa_metasub){
    taxa_overlap <- c(taxa_overlap,i)
  }
} 

### uniq taxa overlap with metasub
taxa_not_overlap <- c()
for (i in taxa_telm){
  if(i %in% taxa_metasub){
    next
  }else{
    taxa_not_overlap <- c(taxa_not_overlap,i)
  }
} 

uniq_both <- c()
for (i in telm_uniq_taxa){
  i2 <- make.names(i)
  if (i2 %in% taxa_not_overlap){
    uniq_both <- c(uniq_both,i)
  }
}  
# write.csv(uniq_both,'telm_uniq_monument_and_metasub.csv',row.names = F)


## abundance comparison
taxa_overlap_access <- data.frame(row.names = c('Public','Non-public'))
for (i in taxa_overlap){
  taxa_overlap_access['Public',i] <- mean(telm_file3[telm_file3$Access == 'public',i]) # 18 samples
  taxa_overlap_access['Non-public',i] <- mean(telm_file3[telm_file3$Access == 'non-public',i]) # 22 samples
}
taxa_overlap_access['Public','Sum'] <- sum(taxa_overlap_access['Public',taxa_overlap])
taxa_overlap_access['Non-public','Sum'] <- sum(taxa_overlap_access['Non-public',taxa_overlap])
taxa_overlap_access$Sum 

## count comparison
sum(taxa_overlap_access['Public',taxa_overlap] != 0) # 1622 -> 1566
sum(taxa_overlap_access['Non-public',taxa_overlap] != 0) # 2023 -> 1949


###### pathogen in each sample
taxa_anno_file <- read.csv('monument_taxonomic_level-Adjust_final.csv',header = T)

taxa_anno_file_telm <- taxa_anno_file[which(make.names(taxa_anno_file$taxa) %in% taxa_telm),]

taxa_anno_file_telm <- droplevels(taxa_anno_file_telm)

taxa_anno_file_telm[is.na(taxa_anno_file_telm)] <- 0

animal_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 1,'taxa'] 
plant_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 2,'taxa'] 
both_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 3,'taxa'] 


# 
# abun_df <- read.csv('Tel_Megiddo_access_taxa_abun.csv',header = T,check.names = F)
# rownames(abun_df) <- abun_df[,1]
# abun_df <- abun_df[,-1]
pathogen_abun <- data.frame(row.names = rownames(abun_df))

for (i in 1:dim(pathogen_abun)[1]){
  pathogen_abun[i,'animal_pathogen_n'] <- sum(abun_df[i,animal_pathogen] > 0)
  pathogen_abun[i,'animal_pathogen_p'] <- sum(as.numeric(abun_df[i,animal_pathogen]))
  pathogen_abun[i,'plant_pathogen_n'] <- sum(abun_df[i,plant_pathogen] > 0)
  pathogen_abun[i,'plant_pathogen_p'] <- sum(as.numeric(abun_df[i,plant_pathogen]))
  pathogen_abun[i,'both_pathogen_n'] <- sum(abun_df[i,both_pathogen] > 0)
  pathogen_abun[i,'both_pathogen_p'] <- sum(as.numeric(abun_df[i,both_pathogen]))
}
# write.csv(pathogen_abun,'access/221123pathogen_abun_access.csv')

sample_pathogen <- data.frame(row.names = c(animal_pathogen,plant_pathogen,both_pathogen))
sample_pathogen[,'pathogen_type'] <- c(rep('animal_pathogen',length(animal_pathogen)),rep('plant_pathogen',length(plant_pathogen)),rep('both_pathogen',length(both_pathogen)))

for(i in telm_file2$sample_id){
  for(j in rownames(sample_pathogen)){
    if (telm_file2[i,j]!=0){
      sample_pathogen[j,i] <- 1
    }else if(telm_file2[i,j]==0){
      sample_pathogen[j,i] <- 0
    }
  }
}

animal_pathogen_n <- colSums(sample_pathogen[sample_pathogen$pathogen_type=='animal_pathogen',2:41])
plant_pathogen_n <- colSums(sample_pathogen[sample_pathogen$pathogen_type=='plant_pathogen',2:41])
both_pathogen_n <- colSums(sample_pathogen[sample_pathogen$pathogen_type=='both_pathogen',2:41])

sample_pathogen <- rbind.data.frame(c('animal_pathogen_sum',animal_pathogen_n),c('plant_pathogen_sum',plant_pathogen_n),c('both_pathogen_sum',both_pathogen_n),sample_pathogen)
rownames(sample_pathogen)[1:3] <- c('Number of animal pathogen','Number of plant pathogen','Number of both pathogen')

# write.csv(sample_pathogen,'Telm_sample_pathogen.csv')


## uniq taxa in TelM table -----------------------------------------------------

telm_uniq_abun <- monument_abun[,c(telm_uniq_taxa)]
telm_uniq_abun <- telm_uniq_abun[-1,]
telm_uniq_abun <- t(telm_uniq_abun)
# write.csv(telm_uniq_abun,"TelM_uniq_taxa.csv")

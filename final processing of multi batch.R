library(readxl)
library(tidyverse)
library(magrittr)

Height_0_20211020122<-read.delim('UCDavis\\Research\\Madison temp\\Height_0_20211020122.txt',header = T, sep = '\t')
labels_pos<-read_excel('UCDavis\\Research\\Madison temp\\Labels.xlsx','in')
height_filter<-Height_0_20211020122[4:nrow(Height_0_20211020122),]
colnames(height_filter)<-Height_0_20211020122[3,]

#S/N filter
height_filter<-height_filter%>%
  mutate(`S/N average` = as.numeric(`S/N average`))%>%
  filter(`S/N average`>10.)
#Blank filter
height_filter[,23:133]<-height_filter%>%
  select(23:133)%>%
  sapply(as.numeric)

height_filter_t<-height_transpose(height_filter,labels_pos,23,133)

blank_df <- filter(height_filter_t,Type == 'blank')%>%
  select(10:ncol(height_filter_t))%>%
  apply(2,mean)%>%
  data.frame()%>%
  set_colnames(c('Blanks'))
blank_df<- rownames_to_column(blank_df, var = 'Alignment ID')

blank_df2 <- filter(height_filter_t,!is.na(Month))%>%
  select(10:ncol(height_filter_t))%>%
  apply(2,max)%>%
  data.frame()%>%
  set_colnames(c('Samples'))%>%
  rownames_to_column(var = 'Alignment ID')

blank_df <- full_join(blank_df,blank_df2,by = 'Alignment ID')
blank_df <- mutate(blank_df,Sample_over_blank = Samples/Blanks)

height_filter<- filter(height_filter,`Alignment ID`%in%(blank_df$`Alignment ID`[blank_df$Sample_over_blank>10]))

#Retention time filter
height_filter<-height_filter%>%
  mutate(`Average Rt(min)` = as.numeric(`Average Rt(min)`))%>%
  filter(`Average Rt(min)`>=4.5)

colnames(height_filter)<-colnames(Height_0_20211020122)
height_filter[,c(2,19,23:133)]<-height_filter%>%
  select(c(2,19,23:133))%>%
  sapply(as.character)
Height_0_20211020122<-Height_0_20211020122[1:3,]%>%
  bind_rows(height_filter)

#Write to file to input to MSFLO
write_tsv(Height_0_20211020122,'Height_0_20211020122_filtered.txt',col_names = F, quote = "none")
########################
#  Put into MSFLO      #
########################
multi_batchNH4<-read_excel('UCDavis\\Research\\Madison temp\\Height_0_20211020122_filtered_NH3_MSFLO\\Height_0_20211020122_filtered_processed.xlsx')
multi_batchNa<-read_excel('UCDavis\\Research\\Madison temp\\Height_0_20211020122_filtered_Na\\Height_0_20211020122_filtered_processed_Na.xlsx')

#Remove unnecessary metadata columns
multi_batchNa<-multi_batchNa%>%
  select(c(1:7,21,23:135))
multi_batchNH4<-multi_batchNH4%>%
  select(c(1:7,21,23:135))

#Merging the two MS-FLO outputs is tricky; some [M+H]+ might have been merged with an [M+Na]+ and [M+NH4]+
NH4_unmatched<-filter(multi_batchNH4,str_detect(adduct_flag,'Matched',negate = T)|is.na(adduct_flag))
Na_unmatched<-filter(multi_batchNa,str_detect(adduct_flag,'Matched',negate = T)|is.na(adduct_flag))
unmatched<-filter(Na_unmatched,`Alignment ID`%in%NH4_unmatched$`Alignment ID`)
NH4_matched<-filter(multi_batchNH4,str_detect(adduct_flag,'Matched'))
NH4_matched$`Adduct ion name` <-'[M+H]+'
Na_matched<-filter(multi_batchNa,str_detect(adduct_flag,'Matched'))
Na_matched$`Adduct ion name`<-<-'[M+H]+'
split_AID<-function(x){
  x = str_split(x,"_")[[1]]
  # x = as.data.frame(c(AID1 = x[1],AID2 = x[2]))
}

NH4_AID<-NH4_matched%>%
  select(`Alignment ID`)%>%
  apply(1,split_AID)%>%
  t()%>%
  as.data.frame()
Na_AID<-Na_matched%>%
  select(`Alignment ID`)%>%
  apply(1,split_AID)%>%
  t()%>%
  as.data.frame()

na_matched_list<-sapply(Na_matched$`Alignment ID`,str_split,"_")%>%
  unlist(recursive = F)%>%
  unname()
nh4_matched_list<-sapply(NH4_matched$`Alignment ID`,str_split,"_")%>%
  unlist(recursive = F)%>%
  unname()  
overlap<-nh4_matched_list[nh4_matched_list%in%na_matched_list]
overlap2<-na_matched_list[na_matched_list%in%nh4_matched_list]

NH4_matched$AID_1<-NH4_AID$V1
NH4_matched$AID_2<-NH4_AID$V2
Na_matched$AID_1<-Na_AID$V1
Na_matched$AID_2<-Na_AID$V2

overlap_Na<-filter(Na_matched,AID_1%in%overlap|AID_2%in%overlap)
overlap_Na$overlap = 0
overlap_Na$overlap[overlap_Na$AID_1%in%overlap] = 1
overlap_Na$overlap[overlap_Na$AID_2%in%overlap] = 2
overlap_Na<-filter(overlap_Na,overlap==1)

overlap_NH4<-filter(NH4_matched,AID_1%in%overlap|AID_2%in%overlap)
overlap_NH4$overlap = 0
overlap_NH4$overlap[overlap_NH4$AID_1%in%overlap] = 1
overlap_NH4$overlap[overlap_NH4$AID_2%in%overlap] = 2
overlap_NH4<-filter(overlap_NH4,overlap==1)

# Fortunately AID_1 (presumably M+H) matches! so we can join by this
# Merging text and adding peak heights
overlap_merge<-bind_rows(overlap_Na,overlap_NH4)

text_merge<-overlap_merge%>%
  select(c(1:10,122:123))%>%
  gather(column,value,c(1:10,12))%>%
  group_by(AID_1,column)%>%
  summarise(sum = paste(value,collapse = "_"))%>%
  spread(column,sum)

height_merge<-overlap_merge%>%
  select(c(11:122))%>%
  gather(sample_name,value,1:111)%>%
  group_by(AID_1,sample_name)%>%
  summarise(sum_height = sum(value))%>%
  spread(sample_name,sum_height)

merged_df<-full_join(text_merge,height_merge,by = "AID_1")
merged_df$`Adduct ion name`<-'[M+H]+'

Na_matched<-filter(Na_matched,!(AID_1%in%overlap_Na$AID_1|AID_2%in%overlap_Na$AID_2))
NH4_matched<-filter(NH4_matched,!(AID_1%in%overlap_NH4$AID_1|AID_2%in%overlap_NH4$AID_2))

final_df2<-bind_rows(unmatched,merged_df)%>%
  bind_rows(Na_matched)%>%
  bind_rows(NH4_matched)

# Need to fix up m/z and Rts so that we can use the split-feature joining script:
fun0 <- function(df){
  df = sapply(df,as.numeric)
  df = min(df)
  return(df)
}

fun1 <- function(df){
  df = sapply(df,as.numeric)
  df = mean(df)
  return(df)
}

final_df2$`Average Mz`<- final_df2%>%
  select(`Average Mz`)%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun0)

final_df2$`Average Rt(min)`<-final_df2%>%
  select(`Average Rt(min)`)%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun1)

final_df2$`MS1 isotopic spectrum`<-final_df2%>%
  select(9)%>%
  apply(1,str_split, "_")%>%
  data.frame()%>%
  t()%>%
  as.data.frame()%>%
  set_rownames(NULL)%>%
  use_series('V1')

write_csv(final_df2,'final_df_WW_20220114.csv')
############################
# to split_join.R          #
############################

pos_df<-read_csv('multi_split_joined_df_20220114_02.csv')%>%
  filter(!is.na(`Alignment ID`))
pos_t<-height_transpose(pos_df,labels_pos,10,119)
pos_t<-filter(pos_t,!is.na(new_name))

# Detection counting in samples for detection frequency filtering
sample_count2 <- filter(pos_t, !is.na(new_name))%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/56)%>%
  data.frame()%>%
  rownames_to_column(var = 'Alignment ID')%>%
  set_colnames(c('Alignment ID','Counts'))

sample_count2$May <- filter(pos_t, (Month == 'May'))%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$June <- filter(pos_t, Month == 'June')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$July <- filter(pos_t, (Month == 'July'& Type != 'special'))%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$August <- filter(pos_t, Month == 'August')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$September <- filter(pos_t, (Month == 'September'& Type != 'special'))%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$November <- filter(pos_t, Month == 'November')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$January <- filter(pos_t, Month == 'January')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$MV2 <- filter(pos_t, Location == 'MV 7-2')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$DM <- filter(pos_t, Location == 'DM 9-1')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA4 <- filter(pos_t, Location == 'CPA4 4-1')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$EPA <- filter(pos_t, Location == 'EPA 5-1')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA18 <- filter(pos_t, Location == 'CPA 18-1')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA2 <- filter(pos_t, Location == 'CPA 2-1')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$In <- filter(pos_t, Location == 'Influent')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$Eff <- filter(pos_t, Location == 'Effluent')%>%
  select(12:ncol(pos_t))%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$num_groups_high_detect_6<- select(sample_count2,"May":"Eff")%>%
  apply(1,function(c)sum(c>=0.6))

filter(sample_count2,num_groups_high_detect_6>=1)%>%
  nrow()
pos_df <- full_join(pos_df,sample_count2,by = 'Alignment ID')

pos_df<-filter(pos_df,num_groups_high_detect_6>=1)

pos_t<-height_transpose(pos_df,labels_pos,10,119)%>%
  filter(!is.na(new_name))

##########################################################

effs<-filter(labels_pos,Type =='effluent')%>%use_series('Sample name')
pos_df[,effs]<-pos_df[,effs]/5

#I want a file to have for posterity:
write_csv(pos_df,'pos_df_20220114.csv')
pos_df<-read_csv("pos_df_20220114.csv")

pos_df<-pos_df%>%
  select(1:120)%>%
  select(c(colnames(pos_df)[1:10],labels_pos$`Sample name`[!is.na(labels_pos$new_name)]))

pos_df<-pos_df[,c(1:5,7:117,119)]
pos_log_t<-pos_df%>%
  column_to_rownames(var = 'Alignment ID')%>%
  select(7:116)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(labels_pos,by ='Sample name')%>%
  filter(!Type%in%c('blank','Spike','Standard'))%>%
  #filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')

n_feat<-nrow(pos_df)  
pos_log_t[,1:n_feat]<-log2(pos_log_t[,1:n_feat]+1)
pos_log<-t(pos_log_t)
library(preprocessCore)
pos_log_quant_norm<-pos_log_t%>%
  #column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  t()%>%
  normalize.quantiles(copy = T)%>%
  as.data.frame()%>%
  set_colnames(labels_pos$`Sample name`[!labels_pos$Type%in%c("blank","Spike","Standard")])%>%
  set_rownames(pos_df$`Alignment ID`)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(labels_pos,by ='Sample name')

pos_log_quant_norm<-pos_log_quant_norm%>%
  filter(!is.na(Month))

# Let's look at box plots to see what we're dealing with:
pos_log_quant_norm%>%
  #  t()%>%
  #  as.data.frame()%>%
  #  rownames_to_column(var = 'Sample name')%>%
  #  full_join(pos_t[,c('Sample name','Month')],by = 'Sample name')%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  #  filter(!is.na(Month))%>%
  gather(feature,abundance,2:n_feat)%>%
  ggplot(aes(x = Month, y = abundance))+
  geom_boxplot()

pos_log_quant_norm%>%write_csv('pos_log_quant_norm_20220114.csv')
pos_pca<- pos_log_quant_norm%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  prcomp()

library(factoextra)
fviz_eig(pos_pca)
pos_ind_coord <- pos_pca$x %>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(labels_pos[,c('Sample name','Month','new_name')],by = 'Sample name')%>%
  filter(!is.na(`PC1`))%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  mutate(Month_location = paste(Month,new_name,sep = " "))

combat_ind_coord%>%
  ggplot(aes(x = `PC1`,y = `PC2`,color = Month,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,14,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

pos_ind_coord$Batch<-3
pos_ind_coord$Batch[pos_ind_coord$Month=='May']<-1
pos_ind_coord$Batch[pos_ind_coord$Month=='June']<-2
pos_ind_coord$Batch[pos_ind_coord$Month=='November']<-4
pos_ind_coord$Batch[pos_ind_coord$Month=='January']<-4

#To make a PC-pair plot:
pairs(pos_ind_coord[,2:6],
      pch = c(3,4,8,1)[pos_ind_coord$Batch],
      col = pos_ind_coord$Month,
      cex = 1.5,
      oma=c(3,3,3,15),
      font.labels = 2)
par(xpd = TRUE)
legend("bottomright",fill = unique(pos_ind_coord$Month), legend = c(levels(pos_ind_coord$Month)),cex = 1,xpd = TRUE)
legend("right",pch = c(3,4,8,1)[unique(pos_ind_coord$Batch)],legend = c('Batch 1','Batch 2','Batch 3','Batch 4'),cex = 1,xpd = TRUE)

#Preparing to make a 3-d PCA that will never see the light of day
colors <- brewer.pal(n = 8, name = 'Dark2')
pos_ind_coord$color<-colors[1]
pos_ind_coord$color[pos_ind_coord$Month=="June"]<-colors[2]
pos_ind_coord$color[pos_ind_coord$Month=="July"]<-colors[3]
pos_ind_coord$color[pos_ind_coord$Month=="August"]<-colors[4]
pos_ind_coord$color[pos_ind_coord$Month=="September"]<-colors[5]
pos_ind_coord$color[pos_ind_coord$Month=="November"]<-colors[6]
pos_ind_coord$color[pos_ind_coord$Month=="January"]<-colors[7]

pos_ind_coord$shape<-4
pos_ind_coord$shape[pos_ind_coord$new_name=="B"]<-15
pos_ind_coord$shape[pos_ind_coord$new_name=="C"]<-14
pos_ind_coord$shape[pos_ind_coord$new_name=="D"]<-17
pos_ind_coord$shape[pos_ind_coord$new_name=="E"]<-18
pos_ind_coord$shape[pos_ind_coord$new_name=="G"]<-19
pos_ind_coord$shape[pos_ind_coord$new_name=="Influent"]<-8
pos_ind_coord$shape[pos_ind_coord$new_name=="Effluent"]<-9


s3d <- scatterplot3d(pos_ind_coord[,c('PC1','PC2','PC3')], 
                     pch = pos_ind_coord$shape, 
                     color=pos_ind_coord$color,
                     cex.symbols = 1.5)
legend(s3d$xyz.convert(500,600,500),legend=c('May','June','July','August','September','November','January'),
       col = colors, pch = 16,cex = 1.25)
legend(s3d$xyz.convert(500,600,-400),legend=c('A','B','C','D','E','G','Influent','Effluent'), 
       pch = pos_ind_coord$shape,cex = 1.25)

########################
# Go to "pvca final.R" #
########################

pos_ind_coord$Month_short<-'May'
pos_ind_coord$Month_short[pos_ind_coord$Month=='June']<-'Jun'
pos_ind_coord$Month_short[pos_ind_coord$Month=='July']<-'Jul'
pos_ind_coord$Month_short[pos_ind_coord$Month=='August']<-'Aug'
pos_ind_coord$Month_short[pos_ind_coord$Month=='September']<-'Sep'
pos_ind_coord$Month_short[pos_ind_coord$Month=='November']<-'Nov'
pos_ind_coord$Month_short[pos_ind_coord$Month=='January']<-'Jan'
pos_ind_coord<-mutate(pos_ind_coord,Month_loc_short = paste(Month_short,new_name,sep = " "))

pos_ind_coord_clust<-pos_ind_coord%>%
  column_to_rownames(var = 'Month_loc_short')%>%
  select(2:7)

pos_distance <- get_dist(pos_ind_coord_clust)
fviz_dist(pos_distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

res.hc <- pos_ind_coord_clust%>%
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering


# Visualize using factoextra
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 7, # Cut in four groups
          cex = 0.75, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
##################
batch<-labels_pos%>%
  filter(!is.na(new_name))%>%
  use_series(Batch)
library(sva)
combat_corrected = pos_log_quant_norm%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  t()%>%
  as.data.frame()%>%
  ComBat(batch=batch, mod=NULL, par.prior=TRUE, prior.plots=T)
combat_corrected<-t(combat_corrected)%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(labels_pos,by = 'Sample name')%>%
  filter(!is.na(Month))

n_feat = nrow(pos_df)
combat_corrected%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  gather(feature,abundance,2:n_feat+1)%>%
  ggplot(aes(x = Month, y = abundance))+
  geom_boxplot()

combat_pca<-combat_corrected%>%
  filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  prcomp()
fviz_eig(combat_pca)
combat_ind_coord <- combat_pca$x %>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(labels_pos[,c('Sample name','Month','new_name')],by = 'Sample name')%>%
  filter(!is.na(`PC1`))%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  mutate(Month_location = paste(Month,new_name,sep = " "))

combat_ind_coord$Month<-as.factor(combat_ind_coord$Month)

combat_ind_coord$Batch<-3
combat_ind_coord$Batch[combat_ind_coord$Month=='May']<-1
combat_ind_coord$Batch[combat_ind_coord$Month=='June']<-2
combat_ind_coord$Batch[combat_ind_coord$Month=='November']<-4
combat_ind_coord$Batch[combat_ind_coord$Month=='January']<-4


pairs(combat_ind_coord[,2:6],
      pch = c(3,4,8,1)[combat_ind_coord$Batch],
      col = combat_ind_coord$Month,
      cex = 1.5,
      oma=c(3,3,3,15),
      font.labels = 2)
par(xpd = TRUE)
legend("bottomright",fill = unique(combat_ind_coord$Month), legend = c(levels(combat_ind_coord$Month)),cex = 0.75,xpd = TRUE)
legend("right",pch = unique(combat_ind_coord$Batch),legend = c('Batch 1','Batch 2','Batch 3','Batch 4'),cex = 0.75,xpd = TRUE)

colors <- brewer.pal(n = 8, name = 'Dark2')
combat_ind_coord$color<-colors[1]
combat_ind_coord$color[combat_ind_coord$Month=="June"]<-colors[2]
combat_ind_coord$color[combat_ind_coord$Month=="July"]<-colors[3]
combat_ind_coord$color[combat_ind_coord$Month=="August"]<-colors[4]
combat_ind_coord$color[combat_ind_coord$Month=="September"]<-colors[5]
combat_ind_coord$color[combat_ind_coord$Month=="November"]<-colors[6]
combat_ind_coord$color[combat_ind_coord$Month=="January"]<-colors[7]

combat_ind_coord$shape<-4
combat_ind_coord$shape[combat_ind_coord$new_name=="B"]<-15
combat_ind_coord$shape[combat_ind_coord$new_name=="C"]<-14
combat_ind_coord$shape[combat_ind_coord$new_name=="D"]<-17
combat_ind_coord$shape[combat_ind_coord$new_name=="E"]<-18
combat_ind_coord$shape[combat_ind_coord$new_name=="G"]<-19
combat_ind_coord$shape[combat_ind_coord$new_name=="Effluent"]<-8
combat_ind_coord$shape[combat_ind_coord$new_name=="Influent"]<-9

s3d <- scatterplot3d(combat_ind_coord[,c('PC1','PC2','PC3')], 
                     pch = combat_ind_coord$shape, 
                     color=combat_ind_coord$color,
                     cex.symbols = 1.5)
legend(s3d$xyz.convert(300,600,275),legend=c('May','June','July','August','September','November','January'),
       col = colors, pch = 16,cex = 1.25)
legend(s3d$xyz.convert(300,600,-650),legend=c('A','B','C','D','E','G','Influent','Effluent'), 
       pch = combat_ind_coord$shape,cex = 1.25)

combat_ind_coord_clust<-combat_ind_coord%>%
  full_join(pos_ind_coord[,c('Month_location','Month_loc_short')],by = 'Month_location')%>%
  column_to_rownames(var = 'Month_loc_short')%>%
  select(2:14)

combat_distance <- get_dist(combat_ind_coord_clust)
fviz_dist(combat_distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

res.hc2 <-combat_ind_coord_clust%>%
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

fviz_dend(res.hc2, k = 7, # Cut in seven groups
          cex = 0.75, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

write_csv(combat_corrected,'combat_corrected_20220118.csv')

pos_log_quant_norm$batch_col<-"blue"
pos_log_quant_norm$batch_col[pos_log_quant_norm$Batch==2]<-"green"
pos_log_quant_norm$batch_col[pos_log_quant_norm$Batch==3]<-"orange"
pos_log_quant_norm$batch_col[pos_log_quant_norm$Batch==4]<-"red"
# colored bar according to batch https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend.html
par(mar=c(1,1,1,10))
par(mar=c(1,1,1,10))
as.dendrogram(res.hc) %>%
  raise.dendrogram(100)%>%
  set("labels_cex",value = 0.9)%>%
  #set("labels_col", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k=7) %>%
  #set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k = 7) %>%
  set("labels_color",pos_ind_coord$Batch[res.hc$order])%>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.01) %>% 
  plot(horiz=TRUE,axes = FALSE)

# Add the colored bar
colored_bars(colors = pos_log_quant_norm$batch_col, dend = as.dendrogram(res.hc), rowLabels = NULL,horiz = TRUE)

#tanglegram
dendlist <- dendlist(
  res.hc %>% 
    as.dendrogram()%>%
    set("labels_color",pos_ind_coord$Batch[res.hc$order]) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k = 7),
  res.hc2 %>% 
    as.dendrogram()%>% 
    set("labels_color", combat_ind_coord$Batch[res.hc2$order]) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k = 7)
)

tanglegram(dendlist, 
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=7,
           lwd=2
)

# Graphs for number of clusters
fviz_nbclust(pos_ind_coord_clust,FUN = hcut, method = "wss")
fviz_nbclust(pos_ind_coord_clust,FUN = hcut, method = "sil")
fviz_nbclust(pos_ind_coord_clust,FUN = hcut, method = "gap_stat")

fviz_nbclust(combat_ind_coord_clust,FUN = hcut, method = "wss")
fviz_nbclust(combat_ind_coord_clust,FUN = hcut, method = "sil")
fviz_nbclust(combat_ind_coord_clust,FUN = hcut, method = "gap_stat")

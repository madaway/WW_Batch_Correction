#filtering, feature joining, for SB
Height_0_202162933<-read.delim('UCDavis\\Research\\Madison temp\\mega data WW analysis\\reruns\\Height_0_202162933.txt',header = T, sep = '\t')
Height_0_202162933<-Height_0_202162933[,1:147]
height_filter<-Height_0_202162933[5:nrow(Height_0_202162933),]
colnames(height_filter)<-Height_0_202162933[4,]

# S/N and RT filter
height_filter<-height_filter%>%
  mutate(`S/N average` = as.numeric(`S/N average`))%>%
  filter(`S/N average`>10.)%>%
  mutate(`Average Rt(min)` = as.numeric(`Average Rt(min)`))%>%
  filter(`Average Rt(min)`>=4.5)

# Blank filtering
height_filter[,33:147]<-height_filter%>%
  select(33:147)%>%
  sapply(as.numeric)

height_filter_t<-height_transpose(height_filter,single_batch_labels,33,147)

blank_df <- filter(height_filter_t,Type == 'Blank')%>%
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

colnames(height_filter)<-colnames(Height_0_202162933)
height_filter[,c(2,29,33:147)]<-height_filter%>%
  select(c(2,29,33:147))%>%
  sapply(as.character)
Height_0_202162933<-Height_0_202162933[1:4,]%>%
  bind_rows(height_filter)

# write to file to input to MSFLO
write_tsv(Height_0_202162933,'Height_0_202162933_filtered.txt',col_names = F, quote = "none")
###################
# Submit to MSFLO #
###################
single_batch_labels<-read_excel('UCDavis\\Research\\Madison temp\\mega data WW analysis\\reruns\\Pos reruns.xlsx','Labels')
single_batchNH4<-read_excel('UCDavis\\Research\\Madison temp\\Height_0_202162933_filtered_NH3_MSFLO\\Height_0_202162933_filtered_processed.xlsx')
single_batchNa<-read_excel('UCDavis\\Research\\Madison temp\\Height_0_202162933_filtered_Na\\Height_0_202162933_filtered_processed_Na.xlsx')

# Get rid of unwanted metadata rows
single_batchNa<-single_batchNa%>%
  select(c(2:7,31,33,35:149))

single_batchNH4<-single_batchNH4%>%
  select(c(2:7,31,33,35:149))


NH4_unmatched<-filter(single_batchNH4,str_detect(adduct_flag,'Matched',negate = T)|is.na(adduct_flag))
Na_unmatched<-filter(single_batchNa,str_detect(adduct_flag,'Matched',negate = T)|is.na(adduct_flag))
unmatched<-filter(Na_unmatched,`Alignment ID`%in%NH4_unmatched$`Alignment ID`)
NH4_matched<-filter(single_batchNH4,str_detect(adduct_flag,'Matched'))
Na_matched<-filter(single_batchNa,str_detect(adduct_flag,'Matched'))

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

NH4_matched$AID_1<-NH4_AID$V1
NH4_matched$AID_2<-NH4_AID$V2
Na_matched$AID_1<-Na_AID$V1
Na_matched$AID_2<-Na_AID$V2

overlap_Na<-filter(Na_matched,AID_1%in%overlap|AID_2%in%overlap)
overlap_NH4<-filter(NH4_matched,AID_1%in%overlap|AID_2%in%overlap)
# Fortunately AID_1 (presumably M+H) matches! so we can join by this
overlap_NH4$overlap = 0
overlap_NH4$overlap[overlap_NH4$AID_2%in%overlap] = 2

overlap_merge<-bind_rows(overlap_Na,overlap_NH4)

text_merge<-overlap_merge%>%
  select(c(1:8,124:125))%>%
  gather(column,value,c(1:8,10))%>%
  group_by(AID_1,column)%>%
  summarise(sum = paste(value,collapse = "_"))%>%
  spread(column,sum)

height_merge<-overlap_merge%>%
  select(c(9:124))%>%
  gather(sample_name,value,1:115)%>%
  group_by(AID_1,sample_name)%>%
  summarise(sum_height = sum(value))%>%
  spread(sample_name,sum_height)

merged_df<-full_join(text_merge,height_merge,by = "AID_1")
  
Na_matched<-filter(Na_matched,!(AID_1%in%overlap|AID_2%in%overlap))
NH4_matched<-filter(NH4_matched,!(AID_1%in%overlap|AID_2%in%overlap))

final_df<-bind_rows(unmatched,merged_df)%>%
  bind_rows(Na_matched)%>%
  bind_rows(NH4_matched)

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

final_df$`Average Mz`<- final_df%>%
  select(`Average Mz`)%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun0)

final_df$`Average Rt(min)`<-final_df%>%
  select(`Average Rt(min)`)%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun1)

final_df$`MS1 isotopic spectrum`<-final_df%>%
  select(8)%>%
  apply(1,str_split, "_")%>%
  data.frame()%>%
  t()%>%
  as.data.frame()%>%
  set_rownames(NULL)%>%
  use_series('V1')
final_df$`Adduct type`[str_detect(final_df$adduct_flag,'Matched')]<-'[M+H]+'

write_csv(final_df,'single adduct joined 20211223.csv')  
################################
# take final_df to split_join.R#
################################
single_batch_height<-read_csv("single_split_joined_df_20211223.csv")
colnames(single_batch_height)
single_batch_t<-height_transpose(single_batch_height,single_batch_labels,9,123)
colnames(single_batch_labels)[7]<-'new_name'
single_batch_t<-single_batch_t%>%filter(!is.na(new_name))

# Detection frequency filter
sample_count2 <- filter(single_batch_t, !is.na(new_name))%>%
  select(8:ncol(single_batch_t))%>%
  apply(2,function(c)sum(c>3000)/56)%>%
  data.frame()%>%
  rownames_to_column(var = 'Alignment ID')%>%
  set_colnames(c('Alignment ID','Counts'))

sample_count2$May <- filter(single_batch_t, Month == 'May')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$June <- filter(single_batch_t, Month == 'June')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$July <- filter(single_batch_t, (Month == 'July'& is.na(Specialty)&Type=='Sample'))%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$August <- filter(single_batch_t, Month == 'August')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$September <- filter(single_batch_t, (Month == 'September'& is.na(Specialty)))%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$November <- filter(single_batch_t, Month == 'November')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$January <- filter(single_batch_t, Month == 'January')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/8)

sample_count2$MV2 <- filter(single_batch_t, Location == 'MV 7-2')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$DM <- filter(single_batch_t, Location == 'DM 9-1')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA4 <- filter(single_batch_t, Location == 'CPA4 4-1')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$EPA <- filter(single_batch_t, Location == 'EPA 5-1')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA18 <- filter(single_batch_t, Location == 'CPA 18-1')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$CPA2 <- filter(single_batch_t, Location == 'CPA 2-1')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$In <- filter(single_batch_t, Location == 'Influent')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$Eff <- filter(single_batch_t, Location == 'Effluent')%>%
  select(colnames(single_batch_t)[!colnames(single_batch_t)%in%colnames(single_batch_labels)])%>%
  apply(2,function(c)sum(c>3000)/7)

sample_count2$num_groups_high_detect<- select(sample_count2,"May":"Eff")%>%
  apply(1,function(c)sum(c>=0.6))
filter(sample_count2,num_groups_high_detect>=1)%>%
  nrow()

single_batch_height <- full_join(single_batch_height,sample_count2,by = 'Alignment ID')

single_batch_height<-filter(single_batch_height,num_groups_high_detect>=1)
single_batch_t <- height_transpose(single_batch_height,single_batch_labels,9,123)%>%filter(!is.na(new_name))

#Divide effluent samples by 5 since effluent sample volume was 5x influent/trunkline volumne
effs_s<-filter(single_batch_labels,Location == 'Effluent')%>%use_series('Sample name')
single_batch_height[,effs_s]<-single_batch_height[,effs_s]/5
single_batch_t<-height_transpose(single_batch_height,single_batch_labels,9,123)%>%filter(!is.na(new_name))

# log2 transform
single_batch_log_t<-filter(single_batch_t,!is.na(Month))%>%
  column_to_rownames('Sample name')%>%
  select(7:(ncol(single_batch_t)-1))
single_batch_log_t<-log2(single_batch_log_t+1)
single_batch_log<-t(single_batch_log_t)

# quantile normalization
library(preprocessCore)
single_log_quant_norm<-normalize.quantiles(single_batch_log,copy = T)%>%
  as.data.frame()%>%
  set_colnames(rownames(single_batch_log_t))%>%
  set_rownames(colnames(single_batch_log_t))

n_feat<-nrow(single_batch_height)

# Let's look at box plots to see what we're dealing with:
single_log_quant_norm%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(single_batch_labels[,c('Sample name','Month')],by = 'Sample name')%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  filter(!is.na(Month))%>%
  gather(feature,abundance,2:n_feat)%>%
  ggplot(aes(x = Month, y = abundance))+
  geom_boxplot()

# PCA
single_pca<- single_log_quant_norm%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(single_batch_labels,by = 'Sample name')%>%
  filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  prcomp()

library(factoextra)
fviz_eig(single_pca)
single_ind_coord <- single_pca$x %>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')%>%
  full_join(single_batch_labels,by = 'Sample name')%>%
  filter(!is.na(`PC1`))%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  mutate(Month_location = paste(Month,new_name,sep = " "))

single_ind_coord$Batch<-3
single_ind_coord$Batch[single_ind_coord$Month=='May']<-1
single_ind_coord$Batch[single_ind_coord$Month=='June']<-2
single_ind_coord$Batch[single_ind_coord$Month=='November']<-4
single_ind_coord$Batch[single_ind_coord$Month=='January']<-4

single_ind_coord$Month<-as.factor(single_ind_coord$Month)

# PC Pair plot
pairs(single_ind_coord[,2:6],
      pch = c(3,4,8,1)[single_ind_coord$Batch],
      col = single_ind_coord$Month,
      cex = 1.5,
      oma=c(3,3,3,15),
      font.labels = 2)
par(xpd = TRUE)
legend("bottomright",fill = unique(single_ind_coord$Month), legend = c(levels(single_ind_coord$Month)),cex = 1,xpd = TRUE)
legend("right",pch = unique(single_ind_coord$Batch),legend = c('Batch 1','Batch 2','Batch 3','Batch 4'),cex = 1,xpd = TRUE)

#pairs(single_ind_coord[,2:6],col = single_ind_coord$Month,oma=c(3,3,3,15))
#par(xpd = TRUE)
#legend("bottomright", fill = unique(single_ind_coord$Month), legend = c(levels(single_ind_coord$Month)))


single_ind_coord%>%
  ggplot(aes(x = `PC1`,y = `PC2`,color = Month,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,16,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

library(RColorBrewer)
colors <- brewer.pal(n = 8, name = 'Dark2')
single_ind_coord$color<-colors[1]
single_ind_coord$color[single_ind_coord$Month=="June"]<-colors[2]
single_ind_coord$color[single_ind_coord$Month=="July"]<-colors[3]
single_ind_coord$color[single_ind_coord$Month=="August"]<-colors[4]
single_ind_coord$color[single_ind_coord$Month=="September"]<-colors[5]
single_ind_coord$color[single_ind_coord$Month=="November"]<-colors[6]
single_ind_coord$color[single_ind_coord$Month=="January"]<-colors[7]

single_ind_coord$shape<-4
single_ind_coord$shape[single_ind_coord$new_name=="B"]<-15
single_ind_coord$shape[single_ind_coord$new_name=="C"]<-16
single_ind_coord$shape[single_ind_coord$new_name=="D"]<-17
single_ind_coord$shape[single_ind_coord$new_name=="E"]<-18
single_ind_coord$shape[single_ind_coord$new_name=="G"]<-19
single_ind_coord$shape[single_ind_coord$new_name=="Effluent"]<-8
single_ind_coord$shape[single_ind_coord$new_name=="Effluent"]<-9

library(scatterplot3d)
s3d <- scatterplot3d(single_ind_coord[,c('PC1','PC2','PC3')], 
                     pch = single_ind_coord$shape, 
                     color=single_ind_coord$color)
legend(s3d$xyz.convert(500,50,220),legend=c('May','June','July','August','September','November','January'),
       col = colors, pch = 16,cex = 1.25)
legend(s3d$xyz.convert(500,50,-10),legend=c('A','B','C','D','E','G','Influent','Effluent'), 
       pch = single_ind_coord$shape,cex = 1.25)
single_ind_coord$month_short<-"May"
single_ind_coord$month_short[single_ind_coord$Month=="June"] <- "Jun"
single_ind_coord$month_short[single_ind_coord$Month=="July"] <- "Jul"
single_ind_coord$month_short[single_ind_coord$Month=="August"] <- "Aug"
single_ind_coord$month_short[single_ind_coord$Month=="September"] <- "Sep"
single_ind_coord$month_short[single_ind_coord$Month=="November"] <- "Nov"
single_ind_coord$month_short[single_ind_coord$Month=="January"] <- "Jan"
single_ind_coord<-mutate(single_ind_coord,Month_location2 = paste(month_short,new_name,sep = " "))
ind_coord_clust<-single_ind_coord%>%
  column_to_rownames(var = 'Month_location2')%>%
  select(2:17)

distance <- get_dist(ind_coord_clust)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
res.hc <- ind_coord_clust%>%
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

  # Visualize using factoextra
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

par(mar=c(1,1,1,10))
as.dendrogram(res.hc) %>%
  set("labels_cex",value = 0.9)%>%
  set("labels_col", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k=7) %>%
  set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k = 7) %>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.01) %>% 
  plot(horiz=TRUE,axes = FALSE)

single_ind_coord%>%
  ggplot(aes(x = `PC1`,y = `PC2`,color = Month,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,14,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

write_csv(single_log_quant_norm,'final single log quant 20211224.csv')
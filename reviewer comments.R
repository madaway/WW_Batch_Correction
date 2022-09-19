library(tidyverse)
library(magrittr)
library(readxl)
library(ggrepel)

istds_df<-read_excel("UCDavis\\Research\\Madison temp\\mega data WW analysis\\ISTD maybe.xlsx","improved MB")
istds_df$`Average Rt(min)`
pos_df%>%
  ggplot(aes(x = `Average Rt(min)`))+
  stat_ecdf(geom = "point")+
  #geom_histogram()+
  geom_point(data = istds_df,aes(x = `Average Rt(min)`,y = c(0.25,0.25,0.25,0.25,0.25,0.25,0.25)),color = "red",size = 3)+
  geom_text_repel(data = istds_df,aes(x = `Average Rt(min)`,y = c(0.25,0.25,0.25,0.25,0.25,0.25,0.25),label = name))+
  theme(axis.line = element_line(colour = "darkgrey"),
          panel.background = element_blank(),
          text = element_text(size = 14),
          panel.grid = element_line(colour = "lightgrey"))+
  ylab("Cumulative percent features")
  
  istds_df_t<-istds_df%>%
  column_to_rownames(var = "name")%>%
  select(33:143)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(labels_pos,by = "Sample name")



istds_df_t%>%
  gather(ISTD,height,2:8)%>%
  filter(!is.na(new_name))%>%
  mutate(new_name = factor(new_name,levels = c("A","B","C","D","E","G","Influent","Effluent")))%>%
  mutate(Month = factor(Month,levels = c("May","June","July","August","September","November","January")))%>%
  mutate(Batch = factor(Batch,levels = c(1,2,3,4)))%>%
  ggplot(aes(x = new_name,y = height,fill = ISTD))+
  geom_boxplot()+
  theme(axis.line = element_line(colour = "darkgrey"),
      panel.background = element_blank(),
      text = element_text(size = 14),
      panel.grid = element_line(colour = "lightgrey"))+
  # scale_shape_manual(values = c(16,8))+
  scale_y_continuous(trans = "log10")
 
istds_df_t$median_height<-istds_df_t%>%
  select(2:8)%>%
  apply(1,median)
istds_df_t<-filter(istds_df_t,!is.na(new_name))
filter(istds_df_t,median_height<3000)
#order: 65, 74, 89, 92

istds_df_t$median_height[istds_df_t$order==65]<-mean(istds_df_t$median_height[between(istds_df_t$order,64,66)])
istds_df_t$median_height[istds_df_t$order==74]<-mean(istds_df_t$median_height[between(istds_df_t$order,73,75)])
istds_df_t$median_height[istds_df_t$order==89]<-mean(istds_df_t$median_height[between(istds_df_t$order,88,90)])
istds_df_t$median_height[istds_df_t$order==92]<-mean(istds_df_t$median_height[between(istds_df_t$order,91,93)])


pos_log_quant_norm_w_istd<-pos_df%>%
  column_to_rownames(var = "Alignment ID")%>%
  select(9:65)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(istds_df_t,by = "Sample name")%>%
  na.omit()

pos_log_quant_norm_w_istd[,2:25830]<-log2(pos_log_quant_norm_w_istd[,2:25830]+1)
pos_log_quant_norm_w_istd<-pos_log_quant_norm_w_istd%>%
  #column_to_rownames(var = 'Sample name')%>%
  select(2:25830)%>%
  t()%>%
  normalize.quantiles(copy = T)%>%
  as.data.frame()%>%
  set_colnames(labels_pos$`Sample name`[!is.na(labels_pos$new_name)])%>%
  set_rownames(colnames(pos_log_quant_norm_w_istd)[2:25830])%>%
  t()%>%
  as.data.frame()


pos_t[,c(12:25833)]<-pos_t[,c(12:25833)]/pos_t$median_height

pos_df$pen_d<-abs(pos_df$`Average Rt(min)` - istds_df$`Average Rt(min)`[istds_df$name=="Pendimethalin-D5"])
min_diff_df<-pos_df%>%
  select(c("Alignment ID","met_d":"pen_d"))%>%
  gather(istd,diff,2:8)%>%
  group_by(`Alignment ID`)%>%
  summarise(diff = min(diff))%>%
  mutate(is_min = 1)
  
min_diff_df<-pos_df%>%
  select(c("Alignment ID","met_d":"pen_d"))%>%
  gather(istd,diff,2:8)%>%
  full_join(min_diff_df,by = c("diff","Alignment ID"))%>%
  filter(!is.na(is_min))

min_diff_df$istd%>%table()
pos_t_scaled_2<-pos_df%>%
  column_to_rownames(var = "Alignment ID")%>%
  select(9:65)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(istds_df_t,by = "Sample name")
pos_t_scaled_2<-na.omit(pos_t_scaled_2)

#low_istd_indx<-which(pos_t_scaled_2$`Pendimethalin-D5`<3000) 
#Pendimethalin in may samples is just really low...
#pos_t_scaled_2$`Pendimethalin-D5`[1:5] = quantile(pos_t_scaled_2$`Pendimethalin-D5`[pos_t_scaled_2$Month=="May"])[4]*rnorm(5,mean = 1)
#istds_df_t<-pos_t_scaled_2[,c(1,25824:25841)]

#for(i in low_istd_indx){
#  pos_t_scaled_2$`Pendimethalin-D5`[i]<-mean(pos_t_scaled_2$`Pendimethalin-D5`[(i-1):(i+1)])
#}

pos_log_quant_norm_w_istd_RT_match<-pos_log_quant_norm_w_istd
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="met_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="met_d"]]/pos_log_quant_norm_w_istd_RT_match$`Methomyl-D3`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="imi_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="imi_d"]]/pos_log_quant_norm_w_istd_RT_match$`Imidacloprid-D4`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="dim_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="dim_d"]]/pos_log_quant_norm_w_istd_RT_match$`Dimethoate-D4`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="sim_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="sim_d"]]/pos_log_quant_norm_w_istd_RT_match$`Simazine-D10`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="diu_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="diu_d"]]/pos_log_quant_norm_w_istd_RT_match$`Diuron-D6`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="bos_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="bos_d"]]/pos_log_quant_norm_w_istd_RT_match$`Boscalid-D4`
pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="pen_d"]]<-pos_log_quant_norm_w_istd_RT_match[,min_diff_df$`Alignment ID`[min_diff_df$istd=="pen_d"]]/pos_log_quant_norm_w_istd_RT_match$`Pendimethalin-D5`

pos_log_quant_norm_w_istd$istd_median<-select(pos_log_quant_norm_w_istd,25823:25829)%>%
  apply(1,median)

pos_log_quant_norm_w_istd_med_scaled<-pos_log_quant_norm_w_istd
pos_log_quant_norm_w_istd_med_scaled[,1:25823]<-pos_log_quant_norm_w_istd_med_scaled[,1:25823]/pos_log_quant_norm_w_istd_med_scaled$istd_median


library(factoextra)

istd_scaled_pca<- pos_log_quant_norm_w_istd_med_scaled%>%  
  #column_to_rownames(var = 'Sample name')%>%
  #select(1:25822)%>%
  select(1:25822)%>%
   prcomp()

fviz_eig(istd_scaled_pca)
istd_scaled_ind_coord <- istd_scaled_pca$x %>%
  as.data.frame()%>%
  rownames_to_column(var = 'Sample name')
View(istd_scaled_ind_coord)

istd_scaled_ind_coord<-istd_scaled_ind_coord%>%
  full_join(labels_pos,by= "Sample name")%>%
  filter(!is.na(`PC1`))%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))%>%
  mutate(Month_location = paste(Month,new_name,sep = " "))

istd_scaled_ind_coord%>%
  ggplot(aes(x = `PC1`,y = `PC2`,color = Month,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,14,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

#To make a PC-pair plot:
pairs(istd_scaled_ind_coord[,2:6],
      pch = c(3,4,8,1)[istd_scaled_ind_coord$Batch],
      col = istd_scaled_ind_coord$Month,
      cex = 1.5,
      oma=c(3,3,3,15),
      font.labels = 2)
par(xpd = TRUE)
legend("bottomright",fill = unique(istd_scaled_ind_coord$Month), legend = c(levels(istd_scaled_ind_coord$Month)),cex = 1,xpd = TRUE)
legend("right",pch = c(3,4,8,1)[unique(istd_scaled_ind_coord$Batch)],legend = c('Batch 1','Batch 2','Batch 3','Batch 4'),cex = 1,xpd = TRUE)

#HCA
labels_pos$Month_short<-labels_pos$Month
labels_pos$Month_short[labels_pos$Month=="June"]<-"Jun"
labels_pos$Month_short[labels_pos$Month=="July"]<-"Jul"
labels_pos$Month_short[labels_pos$Month=="August"]<-"Aug"
labels_pos$Month_short[labels_pos$Month=="September"]<-"Sep"
labels_pos$Month_short[labels_pos$Month=="November"]<-"Nov"
labels_pos$Month_short[labels_pos$Month=="January"]<-"Jan"


istd_scaled_ind_coord_clust<-istd_scaled_ind_coord%>%
  full_join(labels_pos[!is.na(labels_pos$new_name),c("Month_short","Month","new_name")],by = c("Month","new_name"))%>%
  mutate(Month_loc_short = paste(Month_short,new_name,sep = "_"))%>%
  column_to_rownames(var = 'Month_loc_short')

View(istd_scaled_ind_coord_clust)
istd_scaled_ind_coord_clust<-select(istd_scaled_ind_coord_clust,2:(pc_n+1))

res.hc <-istd_scaled_ind_coord_clust%>%
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering
labels_pos$batch_col<-"blue"
labels_pos$batch_col[labels_pos$Batch==2]<-"green"
labels_pos$batch_col[labels_pos$Batch==3]<-"orange"
labels_pos$batch_col[labels_pos$Batch==4]<-"red"
# colored bar according to batch https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend.html
par(mar=c(1,1,1,10))
par(mar=c(1,1,1,10))
library(dendextend)
as.dendrogram(res.hc) %>%
  raise.dendrogram(100)%>%
  set("labels_cex",value = 0.9)%>%
  set("labels_col", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k=7) %>%
  set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid","aquamarine4","deeppink","darkslategray"), k = 7) %>%
  set("labels_color",pos_ind_coord$Batch[res.hc$order])%>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.01) %>% 
  plot(horiz=TRUE,axes = FALSE)

# Add the colored bar
colored_bars(colors = labels_pos$batch_col[!is.na(labels_pos$new_name)], dend = as.dendrogram(res.hc), rowLabels = NULL,horiz = TRUE)


##################################
pct_threshold = .6 # Amount of variability desired to be explained by the principal components.  User can adjust this to a higher (>= 0.8) number but < 1.0

exp_design_multi <- read_excel('UCDavis//Research//Madison temp//Labels.xlsx','exp_design')
exp_design<-exp_design_multi[,c(1:4,6)]
effectsNames_multi = c('Location','Sampling date','Batch','Residual')

## In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model
########## Load data ##########
# theDataMatrix has (presumably) samples as the columns and probes as the rows, so I would want my matrix to be
# in the "height" format, without metadata, and Alignment IDs as rownames

# Considering data set "MB-unC/C"
batch_correction<-'multi'
n_feat<-nrow(pos_df)
# MB-unC
theDataMatrix<-pos_log_quant_norm_w_istd_med_scaled%>%
 # set_rownames(NULL)%>%
#  column_to_rownames(var = 'Sample name')%>%
  select(1:25822)%>%
  t()%>%
  as.data.frame()

##########
theDataMatrix[]<-lapply(theDataMatrix, as.numeric)

dataRowN <- nrow(theDataMatrix)
dataColN <- ncol(theDataMatrix)

########## Center the data (center rows) ##########
theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
#scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
theDataMatrixCentered = t(theDataMatrixCentered_transposed)


expDesignRowN <- nrow(exp_design)
expDesignColN <- ncol(exp_design)
myColNames <- names(exp_design)

theDataMatrixCentered<-na.omit(theDataMatrixCentered)
########## Compute correlation matrix ##########

theDataCor <- cor(theDataMatrixCentered)

########## Obtain eigenvalues ##########

eigenData <- eigen(theDataCor)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum 

########## Merge experimental file and eigenvectors for n components ##########

my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
  my_sum_2  = my_sum_2 - percents_PCs[i]
  if ((my_sum_2) <= pct_threshold ){
    my_counter_2 = my_counter_2 + 1
  }
  
}
if (my_counter_2 < 3){
  pc_n  = 3
  
}else {
  pc_n = my_counter_2 
}

# pc_n is the number of principal components to model
pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
  for (j in 1:expDesignRowN){
    mycounter <- mycounter + 1
    pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
    
  }
}

AAA <- exp_design[rep(1:expDesignRowN,pc_n),]

Data <- cbind(AAA,pc_data_matrix)

####### Edit these variables according to your factors #######
Data$Time <- as.factor(Data$Time)
Data$Treatment<-as.factor(Data$Treatment)
Data$Batch <- as.factor(Data$Batch)
########## Mixed linear model ##########
op <- options(warn = (-1)) 
effects_n = 4
randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
library(lme4)

for (i in 1:pc_n){
  y = (((i-1)*expDesignRowN)+1)
  randomEffects <- as.data.frame(summary(Rm1ML <- lmer(pc_data_matrix ~ (1|Time)+(1|Treatment)+(1|Batch), Data[y:(i*expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit))$varcor)
  for (j in 1:effects_n){
    randomEffectsMatrix[i,j] = as.numeric(randomEffects[j,4])
  }
  
}


effectsNames <- randomEffects[,1]

########## Standardize Variance ##########

randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
  mySum = sum(randomEffectsMatrix[i,])
  for (j in 1:effects_n){
    randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
  }
}

########## Compute Weighted Proportions ##########

randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
  weight = eigenValues[i]/eigenValuesSum
  for (j in 1:effects_n){
    randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
  }
}

########## Compute Weighted Ave Proportions ##########

randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
totalSum = sum(randomEffectsSums)
randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

for (j in 1:effects_n){
  randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
  
}


if(str_detect(batch_correction,'single')){
  randomEffectsMatrixWtAveProp<- t(randomEffectsMatrixWtAveProp)%>%
    as.data.frame()%>%
    set_colnames(c('wt_avg'))%>%
    mutate(Effects = effectsNames_single)
  position = effectsNames_single
} else{
  randomEffectsMatrixWtAveProp<- t(randomEffectsMatrixWtAveProp)%>%
    as.data.frame()%>%
    set_colnames(c('wt_avg'))%>%
    mutate(Effects = effectsNames_multi)  
  position = effectsNames_multi
}

new_values = round(randomEffectsMatrixWtAveProp$wt_avg , 3)

# Generates 1 bar chart per dataset.
# To create a stacked bar chart as in Fig. 1, use weighted average proportion variance for each factor
# as reported in this chart:
ggplot(randomEffectsMatrixWtAveProp,aes(x = Effects,y = wt_avg))+
  geom_col()+
  xlab('Effects')+
  ylab('Weighted average proportion variance')+
  geom_text(label = (new_values),nudge_y = 0.02)+
  scale_x_discrete(limits = position)+
  theme_bw()

pos_df%>%
  column_to_rownames(var = "Alignment ID")%>%
  select(9:65)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(istds_df_t,by = "Sample name")%>%
  filter(!is.na(new_name))%>%
  gather(feature,height,2:(n_feat+1))%>%
  mutate(log_height = log2(height+1))%>%
 # group_by(feature)%>%
  ggplot(aes(x = order,y = log_height,group = feature))+
  geom_line()+
  theme(legend.position = NULL)

pos_t%>%
  gather(feature,height,12:25833)%>%
  ggplot(aes(x = order,y = height,color = feature))+
  geom_line()+
  theme(legend.position = NULL)+
  theme(axis.line = element_line(colour = "darkgrey"),
                                      panel.background = element_blank(),
                                      text = element_text(size = 14),
                                      panel.grid = element_line(colour = "lightgrey"))+
  scale_y_continuous(trans = "log10")

pos_t_scaled_2%>%
  gather(feature,height,2:25823)%>%
  ggplot(aes(x = order,y = height,color = feature))+
  geom_line()+
  theme(legend.position = NULL)+
    theme(axis.line = element_line(colour = "darkgrey"),
          panel.background = element_blank(),
          text = element_text(size = 14),
          panel.grid = element_line(colour = "lightgrey"))+
    scale_y_continuous(trans = "log10")

combat_corrected%>%
  gather(feature,height,2:25823)%>%
  ggplot(aes(x = order,y = height,color = feature))+
  geom_line()+
  theme(legend.position = NULL)+
  theme(axis.line = element_line(colour = "darkgrey"),
        panel.background = element_blank(),
        text = element_text(size = 14),
        panel.grid = element_line(colour = "lightgrey"))+
  scale_y_continuous(trans = "log10")

pvca_stacked<-read_excel("UCDavis//Research//publications//Batchin//PVCA stacked bar.xlsx")%>%select(2:4)

pvca_stacked%>%
  mutate(Factor = fct_relevel(Factor,"Residual" ,"Location", "Sampling date","Analytical batch"))%>%
  mutate(Dataset = fct_relevel(Dataset,"SB","MB-unC","MB-IS","MB-C"))%>%
  ggplot(aes(x = Dataset,y = `Weighted average proportion variance`,fill = Factor))+
  geom_col()+
  theme(axis.line = element_line(colour = "darkgrey"),
        panel.background = element_blank(),
        text = element_text(size = 14),
        panel.grid = element_line(colour = "lightgrey"))


istds_SB<-istds_df<-read_excel("UCDavis\\Research\\Madison temp\\mega data WW analysis\\ISTD maybe.xlsx","improved SB")
single_batch_labels<-read_excel('UCDavis\\Research\\Madison temp\\mega data WW analysis\\reruns\\Pos reruns.xlsx','Labels')
istds_SB<-istds_SB%>%
  column_to_rownames(var = "name")%>%
  select(33:147)%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(single_batch_labels,by = "Sample name")
istds_SB%>%
  gather(ISTD,height,2:8)%>%
  filter(!is.na(`Injection number`))%>%
  ggplot(aes(x = `Injection number`,y = height,shape = ISTD,color = Month))+
  geom_point()+
  scale_y_continuous(trans = "log10")

istds_SB%>%
  full_join(labels_pos[,c("Month","Location","new_name")],by = c("Month","Location"))%>%
  gather(ISTD,height,2:8)%>%
  filter(!is.na(new_name))%>%
  mutate(new_name = factor(new_name,levels = c("A","B","C","D","E","G","Influent","Effluent")))%>%
  mutate(Month = factor(Month,levels = c("May","June","July","August","September","November","January")))%>%
 # mutate(Batch = factor(Batch,levels = c(1,2,3,4)))%>%
  ggplot(aes(x = Month,y = height,fill = ISTD))+
  geom_boxplot()+
  theme(axis.line = element_line(colour = "darkgrey"),
        panel.background = element_blank(),
        text = element_text(size = 14),
        panel.grid = element_line(colour = "lightgrey"))+
  # scale_shape_manual(values = c(16,8))+
  scale_y_continuous(trans = "log10")
#####################
library(Rtsne)
tsne<-Rtsne(as.matrix(combat_corrected[,2:25823]),perplexity = 10)

tsne<-tsne$Y%>%as.data.frame()
colnames(tsne)<-c("X","Y")  
tsne<-bind_cols(tsne,combat_corrected[,25824:25833])

tsne%>%
  mutate(Batch = as.factor(Batch))%>%
  ggplot(aes(x = X,y = Y,color = Batch,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,14,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 14))+
  labs(x = "tSNE1",y = "tSNE2")

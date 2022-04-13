labels_WW <- read_excel('UCDavis//Research//Madison temp//MEGA DATASET.xlsx','Parameters',skip = 1)

labels_WW<-labels_WW%>%
  filter(!is.na(Month))%>%
  filter(!Location%in%c('Groomer','PCO','Laundry','72-I 14-4','MV 7-1','SM 6-1'))

write_csv(labels_WW,'WW_parameters_for_PCA.csv')
labels_WW<-read_excel('WW_parameters_for_PCA.xlsx','WW_parameters_for_PCA')
labels_WW<-labels_WW[1:56,]
labels_WW$short_month<-'May'
labels_WW$short_month[labels_WW$Month=='June']<-'Jun'
labels_WW$short_month[labels_WW$Month=='July']<-'Jul'
labels_WW$short_month[labels_WW$Month=='August']<-'Aug'
labels_WW$short_month[labels_WW$Month=='September']<-'Sep'
labels_WW$short_month[labels_WW$Month=='November']<-'Nov'
labels_WW$short_month[labels_WW$Month=='January']<-'Jan'

labels_WW<-mutate(labels_WW,month_location = paste(short_month,new_name,sep = " "))

labels_WW_PCA<-labels_WW%>%
  mutate(month_location = paste(short_month, new_name, sep = " "))%>%
  column_to_rownames(var = 'month_location')%>%
  select(c("Bifenthrin":"Imidacloprid"))


labels_WW_PCA[is.na(labels_WW_PCA)]<-0

labels_WW_PCA<-labels_WW_PCA%>%
  t()%>%
  as.data.frame()

# Z-scale before PCA
labels_WW_PCA$mean<-apply(labels_WW_PCA,1,mean)
labels_WW_PCA$sd<-apply(labels_WW_PCA,1,sd)

labels_WW_PCA[,1:56]<-(labels_WW_PCA[,1:56]-labels_WW_PCA$mean)/labels_WW_PCA$sd
labels_WW_PCA<-labels_WW_PCA%>%
  select(1:56)%>%
  t()%>%
  as.data.frame()

# PCA
pca<-prcomp(labels_WW_PCA)
fviz_eig(pca)
fviz_pca_biplot(pca)
fviz_pca_var(pca)

ind_coord <- pca$x%>%
  as.data.frame()%>%
  rownames_to_column(var = 'month_location')%>%
  full_join(labels_WW[,c('month_location','Month','short_month','new_name')],by = 'month_location')%>%
  filter(!is.na(`PC1`))%>%
  mutate(Month = fct_relevel(Month, "May",'June','July','August','September','November','January'))
ind_coord$Batch<-3
ind_coord$Batch[ind_coord$Month=='May']<-1
ind_coord$Batch[ind_coord$Month=='June']<-2
ind_coord$Batch[ind_coord$Month=='November']<-4
ind_coord$Batch[ind_coord$Month=='January']<-4

pairs(ind_coord[,2:6],
     pch = c(3,4,8,1)[ind_coord$Batch],
     col = ind_coord$Month,
     cex = 1.5,
     oma=c(3,3,3,15),
     font.labels = 2)
par(xpd = TRUE)
legend("bottomright",fill = unique(ind_coord$Month), legend = c(levels(ind_coord$Month)),cex = 1,xpd = TRUE)
legend("right",pch = c(3,4,8,1)[unique(ind_coord$Batch)],legend = c('Batch 1','Batch 2','Batch 3','Batch 4'),cex = 1,xpd = TRUE)

on.exit(par(opar))

colors <- brewer.pal(n = 8, name = 'Dark2')
ind_coord$color<-colors[1]
ind_coord$color[ind_coord$Month=="June"]<-colors[2]
ind_coord$color[ind_coord$Month=="July"]<-colors[3]
ind_coord$color[ind_coord$Month=="August"]<-colors[4]
ind_coord$color[ind_coord$Month=="September"]<-colors[5]
ind_coord$color[ind_coord$Month=="November"]<-colors[6]
ind_coord$color[ind_coord$Month=="January"]<-colors[7]

s3d <- scatterplot3d(ind_coord[,c('PC1','PC2','PC3')], pch = 16, color=ind_coord$color)
legend(s3d$xyz.convert(10,0,0),legend=c('May','June','July','August','September','November','January'),col = colors, pch = 16)
ind_coord_clust<-ind_coord%>%
  column_to_rownames(var = 'month_location')%>%
  select(2:6)
fviz_pca_var(pca,col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red") +
  theme_minimal()

ind_coord%>%
  ggplot(aes(x = `PC1`,y = `PC2`,color = Month,shape = new_name))+
  geom_point(size = 4)+
  scale_shape_manual(name = 'Site',values = c(4,15,16,17,18,19,8,9))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

pos_distance <- get_dist(ind_coord_clust)
fviz_dist(pos_distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
#res.hc<-ind_coord%>%
#  column_to_rownames(var = 'month_location')%>%
#  select(1:3)%>%
res.hc <- ind_coord_clust%>%
  #res.hc<- pos_log_quant_norm%>%
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
  set("labels_col", value = c("deepskyblue4", "darkorange", "coral4","darkorchid"), k=4) %>%
  set("branches_k_color", value = c("deepskyblue4", "darkorange", "coral4","darkorchid"), k = 4) %>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.01) %>% 
  plot(horiz=TRUE,axes = FALSE)

labels_WW$batch_col<-"orange"
labels_WW$batch_col[labels_WW$Month=='June']<-"green"
labels_WW$batch_col[labels_WW$Month=='May']<-"blue"
labels_WW$batch_col[labels_WW$Month=='January']<-"red"
labels_WW$batch_col[labels_WW$Month=='November']<-"red"

colored_bars(colors = labels_WW$batch_col, dend = as.dendrogram(res.hc), rowLabels = NULL,horiz = TRUE)


# Contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 10)+
  theme(text = element_text(size = 14))
        
# Contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 10)+
  theme(text = element_text(size = 14))
fviz_contrib(pca, choice = "var", axes = 3, top = 10)+
  theme(text = element_text(size = 14))

var <- get_pca_var(pca)
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)   

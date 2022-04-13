library(limma)
#This is necessary for m/z, RT, and other such feature "metadata"
pos_df<-read_csv('pos_df_20220114.csv')

labels_pos<-read_xlsx('UCDavis//Research//Madison temp//Labels.xlsx','in')
# MB-C
combat_corrected<-read_csv("combat_corrected_20220118.csv")
# MB-unC
pos_log_quant_norm<-read_csv("pos_log_quant_norm_20220114.csv")

design_month<-read_excel('UCDavis//Research//Madison temp//Labels.xlsx','design matrix month')
design_site<-read_excel('UCDavis//Research//Madison temp//Labels.xlsx','design matrix site')
colnames(design_site)<-c('single_batch','multi_batch','D','G','C','E','A','B','In','Eff')

design_cluster_uncorrect<-read_excel('UCDavis//Research//Madison temp//Labels.xlsx','design matrix cluster uncorrect')

design_cluster_combat<-read_excel('UCDavis//Research//Madison temp//Labels.xlsx','design matrix cluster combat')
colnames(design_cluster_combat)[1]<-'multi_batch'
################
# Differential Expression in uncorrected, multi-batch dataset
n_feat<-nrow(pos_df)
uncorrected_dat<-pos_log_quant_norm%>%
  #t()%>%
  #as.data.frame()%>%
  #rownames_to_column(var = 'Sample name')%>%
  #full_join(labels_pos,by = 'Sample name')%>%
  #na.omit()%>%
  #filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  t()

uncorrect_month_fit<-lmFit(uncorrected_dat,as.matrix(design_month[,3:9]))
month_contrast = makeContrasts(June-May,July-May,August-May,September-May,November-May,January-May,
                               July-June,August-June,September-June,November-June,January-June,
                               August-July,September-July,November-July, January-July,
                               September-August, November-August,January-August,
                               November-September, January-September,
                               January-November,levels = as.matrix(design_month[,3:9]))
month_fit_contrasts=contrasts.fit(uncorrect_month_fit,month_contrast)
efit_month_contrasts=eBayes(month_fit_contrasts)
par(mfrow=c(3,7))
for (i in 1:ncol(efit_month_contrasts$p.value)){
 hist(efit_month_contrasts$p.value[,i],main=colnames(efit_month_contrasts$p.value)[i],xlab = 'p value')
  }


colnames(pos_df)[1]<-'ID'
path<-'UCDavis//Research//Madison temp//Diff Exp 20220118//Multi//contrast '

for(i in 1:ncol(efit_month_contrasts)){
  table<-topTable(efit_month_contrasts,coef =i,adjust.method = 'BH',n = Inf,genelist = rownames(uncorrected_dat))%>%
    as.data.frame()%>%
    full_join(pos_df[,1:5],by = 'ID')%>%
    na.omit()
  
  if(nrow(table)==0){
    next
  } else{
    table<-table%>% filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(month_fit_contrasts)[i],'.csv',sep = "")
    write_csv(table,filepath)
  }  
  
  
}



uncorrect_site_fit<-lmFit(uncorrected_dat,as.matrix(design_site[,3:10]))
site_contrast = makeContrasts(In-Eff,A-Eff,B-Eff,C-Eff,D-Eff,E-Eff,G-Eff,
                               In-G,A-G,B-G,C-G,D-G,E-G,
                               In-E,A-E,B-E,C-E,D-E,
                               In-D,A-D,B-D,C-D,
                               In-C,A-C,B-C,
                               In-B,A-B,
                               In-A,levels = as.matrix(design_site[,3:10]))
site_fit_contrasts=contrasts.fit(uncorrect_site_fit,site_contrast)
efit_site_contrasts=eBayes(site_fit_contrasts)
par(mfrow=c(4,7))
for (i in 1:ncol(efit_site_contrasts$p.value)){
  hist(efit_site_contrasts$p.value[,i],main=colnames(efit_site_contrasts$p.value)[i],xlab = 'p value')
}

for(i in 1:ncol(efit_site_contrasts)){
  print(i)
  table<-topTable(efit_site_contrasts,coef =i,adjust.method = 'BH',n =Inf,genelist = rownames(uncorrected_dat))
  if(nrow(table)==0){
    next
  } else{
    table<-as.data.frame(table)%>%
    full_join(pos_df[,1:5],by = 'ID')%>%
    na.omit()%>%
    filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(site_fit_contrasts)[i],'.csv',sep = "")
    write_csv(table,filepath)
  }  
}
###############
#clusters ok
uncorrect_cluster_fit<-lmFit(uncorrected_dat,as.matrix(design_cluster_uncorrect[,c(2:8)]))
cluster_contrast = makeContrasts(one-two,one-three,one-four,one-five,one-six,one-seven,
                                 two-three,two-four,two-five,two-six,two-seven,
                                 three-four,three-five,three-six,three-seven,
                                 four-five,four-six,four-seven,
                                 five-six,five-seven,six-seven,
                                 levels = as.matrix(design_cluster_uncorrect[,c(2:8)]))
cluster_fit_contrasts=contrasts.fit(uncorrect_cluster_fit,cluster_contrast)
efit_cluster_contrasts=eBayes(cluster_fit_contrasts)
par(mfrow=c(3,7))
for (i in 1:ncol(efit_cluster_contrasts$p.value)){
  hist(efit_cluster_contrasts$p.value[,i],main=colnames(efit_cluster_contrasts$p.value)[i],xlab = 'p value')
}

for(i in 1:ncol(efit_cluster_contrasts)){
  print(i)
  table<-topTable(efit_cluster_contrasts,coef =i,adjust.method = 'BH',n = Inf,genelist = rownames(uncorrected_dat))
  if(nrow(table)==0){
    next
  } else{
    table<-as.data.frame(table)%>%
      full_join(pos_df[,1:5],by = 'ID')%>%
      na.omit()%>%
      filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(cluster_fit_contrasts)[i],'.csv',sep = "")
    write_csv(table,filepath)
  }  
  
  
  
}

####################
# Differential expression in combat-corrected multi-batch dataset
corrected_dat<-combat_corrected%>%
  filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:nrow(pos_df))%>%
  t()

correct_month_fit<-lmFit(corrected_dat,as.matrix(design_month[,3:9]))
month_fit_contrasts2=contrasts.fit(correct_month_fit,month_contrast)
efit_month_contrasts2=eBayes(month_fit_contrasts2)
par(mfrow=c(3,7))

for (i in 1:ncol(efit_month_contrasts2$p.value)){
  hist(efit_month_contrasts2$p.value[,i],main=colnames(efit_month_contrasts2$p.value)[i],xlab='p value')
}

path<-'UCDavis//Research//Madison temp//Diff Exp 20220118//Combat//contrast '
#colnames(pos_df)[1]<-'ID'
for(i in 1:ncol(efit_month_contrasts2)){
  table<-topTable(efit_month_contrasts2,coef =i,adjust.method = 'BH',n = Inf,genelist = rownames(corrected_dat))
  if(nrow(table)==0){
    next
  } else{
    table<-as.data.frame(table)%>%
      full_join(pos_df[,1:5],by = 'ID')%>%
      na.omit()%>%
      filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(month_fit_contrasts2)[i],'.csv',sep = "")
    write_csv(table,filepath)
  } 
  
}



correct_site_fit<-lmFit(corrected_dat,as.matrix(design_site[,3:10]))
site_fit_contrasts2=contrasts.fit(correct_site_fit,site_contrast)
efit_site_contrasts2=eBayes(site_fit_contrasts2)
par(mfrow=c(4,7))
for (i in 1:ncol(efit_site_contrasts2$p.value)){
  hist(efit_site_contrasts2$p.value[,i],main=colnames(efit_site_contrasts2$p.value)[i],xlab = 'p value')
}

for(i in 1:ncol(efit_site_contrasts2)){
  print(i)
  table<-topTable(efit_site_contrasts2,coef =i,adjust.method = 'BH',n = Inf,genelist = rownames(corrected_dat))
  if(nrow(table)==0){
    next
  } else{
    table<-as.data.frame(table)%>%
      full_join(pos_df[,1:5],by = 'ID')%>%
      na.omit()%>%
      filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(site_fit_contrasts2)[i],'.csv',sep = "")
    write_csv(table,filepath)
  }  
  
}

correct_cluster_fit<-lmFit(corrected_dat,as.matrix(design_cluster_combat[,2:8]))

cluster_fit_contrasts2=contrasts.fit(correct_cluster_fit,cluster_contrast)
efit_cluster_contrasts2=eBayes(cluster_fit_contrasts2)
par(mfrow=c(3,7))
for (i in 1:ncol(efit_cluster_contrasts2$p.value)){
  hist(efit_cluster_contrasts2$p.value[,i],main=colnames(efit_cluster_contrasts2$p.value)[i],xlab = 'p value')
}

for(i in 1:ncol(efit_cluster_contrasts2)){
  print(i)
  table<-topTable(efit_cluster_contrasts2,coef =i,adjust.method = 'BH',n = Inf,genelist = rownames(corrected_dat))
  if(nrow(table)==0){
    next
  } else{
    table<-as.data.frame(table)%>%
      full_join(pos_df[,1:5],by = 'ID')%>%
      na.omit()%>%
      filter(adj.P.Val<0.05)
    filepath<-paste(path,colnames(cluster_fit_contrasts2)[i],'.csv',sep = "")
    write_csv(table,filepath)
  }
}

##############
# Made a sheet of all the file-paths for contrasts that had significant results
file_df<-read_excel('UCDavis//Research//Madison temp//Diff Exp 20220118//file paths.xlsx','Sheet1')

file_df<-file_df%>%
  mutate(file_path = paste('UCDavis//Research//Madison temp//Diff Exp 20220118//',Folder,'//contrast ',Contrast,'.csv',sep = ""))

contrast_table<-file_df%>%
  select(2)%>%
  table()%>%
  as.matrix()%>%
  as.data.frame()%>%
  set_colnames("instances")%>%
  rownames_to_column(var = 'contrast')

max(contrast_table$instances)

contrast_table_1<-filter(contrast_table,instances == 1)
contrast_table_2<-filter(contrast_table,instances == 2)
#############

for(i in 1:nrow(contrast_table)){
  print(i)
  set<-filter(file_df,Contrast ==contrast_table$contrast[1])
  contrast<-set$Contrast[1]
  
  results1<-read_csv(set$file_path[1])%>%
    select(1:11)%>%
    set_colnames(c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B","Average Mz",
                   "Average Rt(min)","adduct_flag","Metabolite name"))%>%
    mutate(dataset = set$Folder[1])%>%
    mutate(contrast = set$Contrast[1])%>%
    mutate(ID = as.character(ID))
    
  results2<-read_csv(set$file_path[2])%>%
    select(1:11)%>%
    set_colnames(c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B","Average Mz",
                   "Average Rt(min)","adduct_flag","Metabolite name"))%>%
    mutate(dataset = set$Folder[2])%>%
    mutate(contrast = set$Contrast[2])%>%
    mutate(ID = as.character(ID))#%>%
  
  results = bind_rows(results1,results2)%>%arrange(ID)
  
  
   
  write_csv(results,paste('UCDavis//Research//Madison temp//Diff Exp 20211208//results ',contrast,'.csv'))
}


#looking at significant differential abundance/expression results
file_df$Folder%>%table()
multi_de <-54
combat_de <- 17

results_multi$contrast%>%table()
results_all<-results1[FALSE,]

results_all<-file_df%>%
  #filter(Folder == 'Multi')%>%
  apply(1,file_combine_fun, results_df = results_all)%>%
  map_dfr(bind_rows)

file_combine_fun<-function(df,results_df){
  dataset_df<-df[1]
  contrast_df = df[2]
  filepath = df[4]
  file<-read_csv(filepath)%>%select(1:11)%>%
    set_colnames(c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B","Average Mz",
                   "Average Rt(min)","adduct_flag","Metabolite name"))%>%
    mutate(ID = as.character(ID))%>%
    mutate(dataset = dataset_df)%>%
    mutate(contrast = contrast_df)
  results_df<-bind_rows(results_df,file)
  return(results_df)
}

results_combat<-results[F,]
results_combat<-file_df%>%
  filter(Folder == 'Combat2')%>%
  apply(1,file_combine_fun,results_df = results_combat)%>%
  map_dfr(bind_rows)

results_single<-results[F,]
results_single<-file_df%>%
  filter(Folder == 'Single')%>%
  apply(1,file_combine_fun,results_df = results_single)%>%
  map_dfr(bind_rows)

results_all<-bind_rows(results_multi,results_combat)

write_csv(results_all,'UCDavis//Research//Madison temp//Diff Exp 20220118//all results.csv')
results_all<-read_csv('UCDavis//Research//Madison temp//Diff Exp 20220118//all results.csv')

results_all%>%
  #filter(contrast%in%contrast_table_2$contrast)%>%
  group_by(dataset)%>%
  select("contrast")%>%
  table()%>%as.data.frame()%>%write_csv('UCDavis//Research//Madison temp//Diff Exp 20211224//cluster results2.csv')
  filter(Freq==2)%>%
  #select(ID)%>%unique()%>%nrow()
  select(contrast)%>%
  table()
results_all$dataset[results_all$dataset=='Combat']<-'MB-C'

results_all$FC_direction<-'increase'
results_all$FC_direction[results_all$logFC<0]<-'decrease'

results_all%>%
  filter(str_detect(dataset,'MB-unC'))%>%
  select('contrast')%>%
  table()%>%
  as.data.frame()%>%
  write_csv('UCDavis//Research//Madison temp//Diff Exp 20220118//MB-unC contras table.csv')

overlap_contrast<-contrast_table%>%
  filter(instances==2)%>%
  filter(str_detect(contrast,'one',negate = T))%>%
  filter(str_detect(contrast,'two',negate = T))%>%
  filter(str_detect(contrast,'three',negate = T))%>%
  filter(str_detect(contrast,'four',negate = T))%>%
  filter(str_detect(contrast,'five',negate = T))%>%
  filter(str_detect(contrast,'six',negate = T))%>%
  select('contrast')%>%unique()
  

for (i in 2:7){
  #i = 1
  month_overlap<-filter(results_all,contrast==overlap_contrast$contrast[i])%>%
    select('ID')%>%
    table()%>%
    as.data.frame()%>%
    filter(Freq>1)%>%
    select(1)%>%
    set_colnames('ID')%>%
    mutate(contrast = overlap_contrast$contrast[i])%>%
    bind_rows(month_overlap)
}


month_overlap%>%
  mutate(keep = 1)%>%
  full_join(results_all,by = c('ID','contrast'))%>%
  na.omit()%>%
  #filter(str_detect(contrast,'July'))%>%
  #mutate(kendrick_mass = `Average Mz`*(32/32.0116))%>% #PEG repeating unit (?)
  #mutate(kendrick_mass = `Average Mz`*(14/14.01565))%>% #CH2
  #mutate(kendrick_mass = `Average Mz`*(44/43.96))%>% #CO2
  #mutate(kendrick_mass_defect = round(kendrick_mass)-kendrick_mass)%>%
  #mutate(mass_defect = `Average Mz` - round(`Average Mz`))%>%
  #ggplot(aes(x = kendrick_mass,y = kendrick_mass_defect,color = contrast))+
  #ggplot(aes(x = `Average Rt(min)`,y = mass_defect,color = contrast,shape = FC_bin))+
  ggplot(aes(x = `Average Rt(min)`,y = `Average Mz`,color = contrast,shape = FC_direction))+
  geom_point(alpha = 0.5,cex = 3)+
  scale_shape_manual(name = 'FC_bin',values = c(17,18))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))

clusters_list<-file_df%>%
  filter(Type == 'Cluster')%>%
  select("Contrast")%>%
  unique()%>%
  use_series(Contrast)

for (i in clusters_list){
  x = filter(results_all,contrast == i)%>%
    select('ID')%>%
    table()%>%
    as.data.frame()%>%
    filter(Freq>1)%>%
    nrow()
  print(c(i,x))
}

list<-c('one','two','three','four','five','six','seven')
#months<-c('May','June','July','August','September','November','January')
const<-read_excel('cluster_constant.xlsx')
le_dataset = 'MB-C'
i = 'one'
COR_cluster_sig_plus =  
filter(results_all,str_detect(contrast,i))%>%
  filter(dataset == le_dataset)%>%
  full_join(const,by = 'contrast')%>%
  na.omit()%>%
  gather(cluster,constant,15:21)%>%
  filter(cluster ==i)%>%
  mutate(FC_new = constant*logFC)%>%
  filter(FC_new>0)%>%
  select('ID')%>%
  table()%>%
  as.data.frame()%>%
  filter(Freq==6)%>%
  set_colnames(c('ID','Freq'))%>%
  mutate(ID = as.character(ID))%>%
  mutate(cluster = i)

for (i in list[2:7]){
  COR_cluster_sig_plus =  filter(results_all,str_detect(contrast,i))%>%
    filter(dataset == le_dataset)%>%
    full_join(const,by = 'contrast')%>%
    na.omit()%>%
    gather(cluster,constant,15:21)%>%
    filter(cluster ==i)%>%
    mutate(FC_new = constant*logFC)%>%
    filter(FC_new > 0)%>%
    select('ID')%>%
    table()%>%
    as.data.frame()%>%
    filter(Freq==6)%>%
    set_colnames(c('ID','Freq'))%>%
    mutate(ID = as.character(ID))%>%
    mutate(cluster = i)%>%
    bind_rows(COR_cluster_sig_plus)
}
#############
for (i in months){
  x1 =  filter(month_sites,str_detect(contrast,i))%>%
    filter(dataset == le_dataset)%>%
    select('ID')%>%
    table()%>%
    as.data.frame()%>%
    filter(Freq==6)%>%
    set_colnames(c('ID','Freq'))%>%
    mutate(ID = as.character(ID))%>%
    mutate(cluster = i)%>%
    bind_rows(x1)
} 

x2<-x1%>%filter(cluster%in%months)
x1<-x1%>%filter(cluster%in%list)%>%full_join(x2[c('ID','Month')],by = 'ID')%>%na.omit

x1<-x1%>%
  mutate(cluster_month = paste(cluster,Month,sep = " "))
x1%>%
  select(cluster_month)%>%
  table()%>%as.data.frame()
##########

COR_cluster_sig_plus$FC_dir <- 'increase'
COR_cluster_sig_minus$FC_dir<-'decrease'

i = list[4]

bind_rows(COR_cluster_sig_minus,COR_cluster_sig_plus)%>%
  full_join(results_all,by = 'ID')%>%
  filter(str_detect(contrast,cluster))%>%
  filter(dataset == 'MB-C')%>%
  filter(FC_dir == "increase")%>%
  #filter(cluster == i)%>%
  #mutate(cluster = factor(cluster,levels = c('one','two','three','four','five','six','seven')))%>%#write_csv('UCDavis//Research//Madison temp//Diff Exp 20220118//clusters significant features cor.csv')
  select(c('ID','Average Mz','Average Rt(min)','cluster','FC_dir'))%>%
  distinct()%>%
  ggplot(aes(x = `Average Rt(min)`,y = `Average Mz`,color =  cluster))+
  #ggplot(aes(x = `Average Mz`,y=mass_defect,color = cluster))+
  geom_point(alpha = 0.5,cex = 3)+
  #scale_shape_manual(name = 'Fold change direction',values = c("square","triangle"))+
  theme_bw()+
  theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13),axis.text = element_text(size = 12),axis.title = element_text(size = 12))+
  #ggtitle(paste('Features significant to cluster',i,sep = " "))+
  scale_x_continuous(limits = c(4.5,25))
  

i = list[7]

bind_rows(UNC_cluster_sig_minus,UNC_cluster_sig_plus)%>%
  full_join(results_all,by = 'ID')%>%
  filter(str_detect(contrast,cluster))%>%
  filter(dataset == 'MB-unC')%>%
  filter(FC_dir == "decrease")%>%
  #mutate(cluster = factor(cluster,levels = c('one','two','three','four','five','six','seven')))%>%#write_csv('UCDavis//Research//Madison temp//Diff Exp 20220118//clusters significant features cor.csv')
  select(c('ID','Average Mz','Average Rt(min)','cluster','FC_dir'))%>%
  distinct()%>%
  ggplot(aes(x = `Average Rt(min)`,y = `Average Mz`,color =  cluster))+
  #ggplot(aes(x = `Average Mz`,y=mass_defect,color = cluster))+
  geom_point(alpha = 0.5,cex = 3)+
  #scale_shape_manual(name = 'Fold change direction',values = c('circle','triangle'))+
  theme_bw()+
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12),axis.text = element_text(size = 12),axis.title = element_text(size = 12))+
  #ggtitle(paste('Features significant to cluster',i,sep = " "))+
  scale_x_continuous(limits = c(4.5,25))
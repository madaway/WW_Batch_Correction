packages<-c('tidyverse','readxl','magrittr','ggplot2','useful','beepr')
lapply(packages,library,character.only = TRUE)

# height_test is our MS-DIAL output 
height_test <- read_csv("final_df_WW_20220114.csv")
colnames(height_test)


labels <-read_excel("UCDavis\\Research\\Madison temp\\Labels.xlsx","in")

# For computation time, remove features that are not in the 56-sample subset that will be retained for analysis
height_test_sample_count<-height_test%>%
  select(c("Alignment ID",first_col:last_col))%>%
  #distinct()%>%
  column_to_rownames(var = "Alignment ID")%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample name")%>%
  full_join(labels,by = "Sample name")%>%
  filter(!is.na(new_name))%>%
  select(2:82178)%>%
  apply(2,function(c)sum(c>3000))%>%
  as.data.frame()%>%
  set_colnames("count")%>%
  rownames_to_column(var = "Alignment ID")

height_test<-height_test%>%
  filter(`Alignment ID`%in%height_test_sample_count$`Alignment ID`[height_test_sample_count$count>0])

#Define tolerances and scoring weights here:
# d_mz and d_RT are used in the initial filtering step that is used to search for neighbor features within a certain window
d_mz<- 0.015 #Da
d_RT<- 15    #s

# RTT and mass_error_tol are used in scoring the groups of neighbor features. Could use more stringent tolerances here than in the initial filtering step if you want to consider the more "marginal" split features
RTT <-15            #s
mass_error_tol<-10  #ppm

#Define minimum height for isotopic peak for scoring MS1 isotopic patterns
min_peak <- 3000

# scoring weights used to create the combined score of RT score, mz score, and MS1 isotopic score
MS1_weight<-1/3
mz_weight<-1/3
rt_weight<-1/3

# define: first_col (first sample name), last_col (last sample name), and n_sample = number of sample columns
# this is important for joining of split features
first_col <-"160522_10_pos_Pest-STD-100"
last_col <- "170303_pos_41_WW_Jan_Spike"
n_sample <- 110

#get_MS1_fun_2 is used to extract the MS1 isotopic spectra from the data frame for MS1 scoring
get_MS1_fun_2 <- function(compound_ID,df){
  MS1_mat<-filter(df,`Alignment ID`==compound_ID)%>%
    magrittr::extract(1,)%>%
    select('MS1 isotopic spectrum')%>%    
    sapply(str_split, " ") %>%
    data.frame() %>%
    t() %>%
    unname() %>%
    data.frame() %>%
    sapply(str_split,":") %>%
    data.frame() %>%
    sapply(as.numeric) %>%
    t()%>%
    as.data.frame()%>%
    set_colnames(c('mz','abundance'))
  
  MS1_mat$rel_abundance = MS1_mat$abundance/MS1_mat$abundance[1]
  
  rownames(MS1_mat)<-c(1:nrow(MS1_mat))
  
  return(MS1_mat)
}


#0. convert RT from min>s (*60)
height_test$`Average Rt(min)`<-as.numeric(height_test$`Average Rt(min)`)*60
height_test$`Average Mz`<-as.numeric(height_test$`Average Mz`)

# 1. Define "neighbor pairs" by
#    a. RT within +/- d_RT
#    b. mz within +/- d_mz

height_test$neighbors<-NA
height_test$`neighbors2`<-NA
height_test$`neighbors3`<-NA
height_test$`neighbors4`<-NA
height_test$`neighbors5`<-NA

for (i in 1:nrow(height_test)){
  mz_i<-height_test$`Average Mz`[i]
  rt_i<-height_test$`Average Rt(min)`[i]
  id_i<-height_test$`Alignment ID`[i]
  adduct_i<-height_test$`Adduct ion name`[i]
  #Check to see if a feature has already been flagged as a neighbor of another feature
  if(id_i%in%height_test$neighbors | id_i%in%height_test$neighbors2 | id_i%in%height_test$neighbors3 | id_i%in%height_test$neighbors4 | id_i%in%height_test$neighbors5){
    next
  }
  #filter all the features to include those within +/- 0.015 Da of mz, +/-0.2 min, with the same adduct type
  neighbor_list<-filter(height_test,(between(`Average Mz`,mz_i-d_mz,mz_i+d_mz) & between(`Average Rt(min)`,rt_i-d_RT,rt_i+d_RT) & `Adduct ion name`==adduct_i))%>%
    #obviously don't want to match it with itself  
    filter(`Alignment ID`!=id_i)%>%
    use_series(`Alignment ID`)
  
  # this is just because R gets mad about saving lists of different lengths in the same column, hence the columns "neighbors",
  # "neighbors2", and "neighbors3"
  if(length(neighbor_list)>4){
    height_test$neighbors[i]<-neighbor_list[1]
    height_test$`neighbors2`[i]<-neighbor_list[2]
    height_test$`neighbors3`[i]<-neighbor_list[3]
    height_test$`neighbors4`[i]<-neighbor_list[4]
    height_test$`neighbors5`[i]<-neighbor_list[5]
  }
  else if(length(neighbor_list)>3){
    height_test$neighbors[i]<-neighbor_list[1]
    height_test$`neighbors2`[i]<-neighbor_list[2]
    height_test$`neighbors3`[i]<-neighbor_list[3]    
    height_test$`neighbors3`[i]<-neighbor_list[4]
    
  }
  else if(length(neighbor_list)>2){
    height_test$neighbors[i]<-neighbor_list[1]
    height_test$`neighbors2`[i]<-neighbor_list[2]
    height_test$`neighbors3`[i]<-neighbor_list[3]
  }
  else if(length(neighbor_list)>1){
    height_test$neighbors[i]<-neighbor_list[1]
    height_test$`neighbors2`[i]<-neighbor_list[2]
  } 
  else if(length(neighbor_list)>0){
    height_test$neighbors[i]<-neighbor_list
  } 
}

neighbor_list<-filter(height_test,!is.na(neighbors))%>%
  select(c('Alignment ID','neighbors','neighbors2','neighbors3','neighbors4','neighbors5'))

height_test<-height_test%>%
  select(1:(ncol(height_test)-5))
nrow_neighbor <-nrow(neighbor_list)
neighbor_list$group_no<-c(1:nrow_neighbor)
neighbor_list<-column_to_rownames(neighbor_list,var = 'group_no')%>%
  t()%>%
  as.data.frame()%>%
  gather(group,`Alignment ID`,1:nrow_neighbor)%>%
  filter(!is.na(`Alignment ID`))

height_test<-full_join(height_test,neighbor_list,by = 'Alignment ID')

# 2. Score groups (aka neighbor pairs) based on RT and m/z
#    a. RT_score = 1.001-(RT_delta/RTT), where RTT is the retention time "tolerance"
#    b. mz_score = 1.001-((max(mz)-min(mz)/mean(mz))*1e6)/mass_error_tol, where mass_error_tol is the maximum tolerated mass error

neighbors_only_df<-filter(height_test,!is.na(group))
group_scores<-group_by(neighbors_only_df,group)%>%
  summarise(RT_score = max(`Average Rt(min)`)-min(`Average Rt(min)`),mz_score = (max(`Average Mz`)-min(`Average Mz`))/mean(`Average Mz`))    

group_scores$RT_score<-1.001-(group_scores$RT_score/RTT)

group_scores$mz_score<-1.001-(group_scores$mz_score*1e6/mass_error_tol)

# 3. Isotope score
#    a. Are the M+1 and M+2 peaks > min_peak (by feature)
#    b. Compute the CV for M+1/M and M+2/M (by group)

#Ignore the warnings generated by this chunk if they say "Warning in FUN(newX[, i], ...) :
#  no non-missing arguments to min; returning Inf", R is just mad about dividing by 0
MS1_df<-data.frame(row.names = c(1:nrow(neighbors_only_df)))
MS1_df$`Alignment ID`<-neighbors_only_df$`Alignment ID`
MS1_df<-distinct(MS1_df)
MS1_df[,c('M','M+1','M+2','abund','abund+1','abund+2','rel abund','rel abund+1','rel abund+2')]<-NA
for(i in 1:nrow(MS1_df)){
  id<-MS1_df$`Alignment ID`[i]
  ms1<-get_MS1_fun_2(id,neighbors_only_df)
  MS1_df[i,2:10]<-flatten(ms1)
}

MS1_df[,c('flag raw M+1','flag raw M+2','flag rel M+1','flag rel M+2')]<-NA
MS1_df$`flag raw M+1`[MS1_df$`abund+1`<min_peak]<-'M+1 less than 3000'
MS1_df$`flag raw M+2`[MS1_df$`abund+2`<min_peak]<-'M+2 less than 3000'
MS1_df$`flag rel M+1`[MS1_df$`rel abund+1`>=1.07]<-'possibly fluorinated'
MS1_df$`flag rel M+2`[MS1_df$`rel abund+2`>=0.2]<-'use M+2'

MS_CV_df<-MS1_df%>%
  full_join(neighbor_list,by = 'Alignment ID')%>%
  #mutate(group = as.factor(group))%>%
  group_by(`group`)%>%
  summarise(`CV_M+1` = sd(`rel abund+1`)/mean(`rel abund+1`),`CV_M+2` = sd(`rel abund+2`)/mean(`rel abund+2`))


MS_CV_df$CV_score<-MS_CV_df%>%
  select(2:3)%>%
  apply(1,min, na.rm = T)

MS_CV_df$CV_score = 1.001 - (MS_CV_df$CV_score/0.2)

# Use all three scores to determine whether features should be joined
# Basically, count all the thresholds that the group meets: 
# is the RT difference less than RTT? is mass error less than mass error tolerance? is the CV of the 
# MS1 isotopic spectrum 20% or less? 
# If yes to all 3, merge; if yes to 2, flag for inspection; if yes to 1, ignore
group_scores<-full_join(group_scores,MS_CV_df,by = 'group')

group_scores$place_holder_rt<-group_scores$RT_score
group_scores$place_holder_mz<-group_scores$mz_score
group_scores$place_holder_cv<-group_scores$CV_score

group_scores$place_holder_rt[group_scores$place_holder_rt<0]<-0
group_scores$place_holder_mz[group_scores$place_holder_mz<0]<-0
group_scores$place_holder_cv[group_scores$place_holder_cv<0]<-0
group_scores<-mutate(group_scores,combined_score = (MS1_weight*place_holder_cv)+(mz_weight*place_holder_mz)+(rt_weight*place_holder_rt))

group_scores$count<-group_scores%>%
  select(7:9)%>%
  apply(1,function(c)sum(c>0))

group_scores$Flag<-NA
group_scores$Flag[group_scores$count==1]<-'Unlikely'
group_scores$Flag[group_scores$count==2]<-'Manual inspect'
group_scores$Flag[group_scores$count==3]<-'Merged'

group_scores<-group_scores[,c(1:6,10,12)]

#Create merge df for feature-groups that will be joined
height_test<-full_join(height_test,MS1_df,by = 'Alignment ID')
height_test<-full_join(height_test,group_scores,by='group')
merge_df<-filter(height_test,Flag=='Merged'&combined_score>0.5)

#Somehow I had an issue with some features getting counted twice, so I wanted to make sure they were only
#counted once
dbl_dippers<-merge_df$`Alignment ID`%>%
  table()%>%
  as.data.frame()%>%
  filter(Freq>1)%>%select(".")
dbl_dippers<-filter(merge_df,`Alignment ID`%in%dbl_dippers$.)
dbl_dippers<-dbl_dippers[,c(1:2,135:154)]
dbl_dippers<-dbl_dippers%>%
  group_by(`Alignment ID`)%>%
  summarise(combined_score = max(combined_score))%>%
  mutate(keep = 1)%>%
  full_join(dbl_dippers,by = c("Alignment ID","combined_score"))
dbl_dippers$keep[is.na(dbl_dippers$keep)]<-0

merge_df<-dbl_dippers%>%
  select(c("Alignment ID","combined_score","keep"))%>%
  full_join(merge_df,by = c("Alignment ID","combined_score"))%>%
  filter((is.na(keep)|keep ==1))

colnames(merge_df)
merge_df<-merge_df[,c(1,4:7,23:150)]
dont_merge_df<-filter(height_test,!`Alignment ID`%in%merge_df$`Alignment ID`)%>%
  select("Alignment ID":last_col)%>%
  distinct()

text_merge<-merge_df%>%
  select(c('Alignment ID':'MS/MS spectrum','neutral_mass':'flag rel M+2'))%>%
  gather(column,value,c(1:8,10:22))%>%
  group_by(group,column)%>%
  summarise(sum = paste(value,collapse = "_"))%>%
  spread(column,sum)

height_merge<-merge_df%>%
  select(c('group',first_col:last_col))%>%
  gather(sample_name,value,2:(n_sample+1))%>%
  group_by(group,sample_name)%>%
  summarise(max_height = max(value))%>%
  spread(sample_name,max_height)

merged_df<-full_join(text_merge,height_merge,by = "group")%>%
  full_join(group_scores,by = 'group')%>%
  filter(Flag == 'Merged')%>%
  mutate(group = as.numeric(group))%>%
  filter(!is.na(abund))

merged_df$`Average Mz`<- ungroup(merged_df)%>%
  select("Average Mz")%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun1)

merged_df$`Average Rt(min)`<-ungroup(merged_df)%>%
  select(`Average Rt(min)`)%>%
  sapply(str_split,"_")%>%
  data.frame()%>%
  apply(1,fun2)


merged_df_t<-merged_df%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "column_names")

dont_merge_df_t<-dont_merge_df%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "column_names")%>%
  mutate(order = c(1:ncol(dont_merge_df)))

split_joined_df<-full_join(merged_df_t,dont_merge_df_t,by = "column_names")%>%
  arrange(order)%>%
  column_to_rownames(var = "column_names")%>%
  t()%>%
  as.data.frame()%>%
  set_rownames(NULL)

split_joined_df<-split_joined_df[1:nrow(split_joined_df)-1,]
split_joined_df$`Average Mz`<-as.numeric(split_joined_df$`Average Mz`)
split_joined_df$`Average Rt(min)`<-as.numeric(split_joined_df$`Average Rt(min)`)
split_joined_df$`Average Rt(min)`[split_joined_df$`Average Rt(min)`>25]=split_joined_df$`Average Rt(min)`[split_joined_df$`Average Rt(min)`>25]/60

#Change file paths here if desired
write_csv(split_joined_df,'multi_split_joined_df_20220114_02.csv')
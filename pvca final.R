packages<-c('tidyverse','readxl','magrittr','ggplot2','useful','beepr')
lapply(packages,library,character.only = TRUE)

########## load libraries ##########
library(lme4)
####### Edit these variables according to user defined parameters and the path to your data and data files names ##########
pct_threshold = .6 # Amount of variability desired to be explained by the principal components.  User can adjust this to a higher (>= 0.8) number but < 1.0

exp_design_multi <- read_excel('UCDavis//Research//Madison temp//Labels.xlsx','exp_design')
exp_design_single<-read_excel('UCDavis\\Research\\Madison temp\\mega data WW analysis\\reruns\\Pos Reruns.xlsx','exp_design')
colnames(exp_design_single)<-colnames(exp_design_multi)

exp_design<-exp_design_multi[,c(1:4,6)]
effectsNames_multi = c('Location','Sampling date','Batch','Residual')

effectsNames_single = c('Location','Sampling date','Batch','Residual')
## In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model
########## Load data ##########
# theDataMatrix has (presumably) samples as the columns and probes as the rows, so I would want my matrix to be
# in the "height" format, without metadata, and Alignment IDs as rownames

# Considering Data set "SB"
batch_correction<-'single'
n_feat<-3108
single_log_quant_norm<-read_csv('final single log quant 20211224.csv')
exp_design<-exp_design_single[,c(1:4,6)]%>%
  arrange(sample)
theDataMatrix<-single_log_quant_norm

# Considering data set "MB-unC/C"
batch_correction<-'multi'
n_feat<-nrow(pos_df)
# MB-unC
theDataMatrix<-pos_log_quant_norm%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  t()%>%
  as.data.frame()
# MB-C
theDataMatrix<-combat_corrected%>%
  filter(!is.na(new_name))%>%
  column_to_rownames(var = 'Sample name')%>%
  select(1:n_feat)%>%
  t()%>%
  as.data.frame

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


library(readxl)
library(dplyr)


my_path='/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/behavior_data/merged_MWM_hsm_090123.csv'

my_path='/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/behavior_data/MWM_080124.csv'
behave=read.csv(my_path ) # %>% select(Animal, )
behave = as.data.frame(behave)
behave$NormSWDist<-behave$SW_Distance/behave$Distance
behave$DistTot<-behave$NE_Distance+behave$NW_Distance+behave$SE_Distance+behave$SW_Distance
behave$SW.Dist.Norm<-behave$SW_Distance/behave$DistTot
behave = behave[behave$Stage=="Probe_D8",]
# bheave
# 
# behave = as.data.frame(t(na.omit(t(behave))))

# 
master_path = '/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/MasterSheet_Experiments2021.xlsx'
#%master_df = read_xlsx(master_path, sheet = "18ABB11_readable02.22.22_BJ_Cor") %>%dplyr::select( DWI,BadeaID, Diet, Age_Handling )#subselect
master_df = read_xlsx(master_path, sheet = "18ABB11_readable02.22.22_BJ_Cor") %>%dplyr::select( DWI,CIVMID, BadeaID, Diet, Age_Months )#subselect
#master_df = read_xlsx(master_path, sheet = "18ABB11_readable02.22.22_BJ_Cor") %>%dplyr::select( DWI)
master_df <- master_df[!is.na(master_df$DWI), ]

master_df$BadeaID = gsub("-","_",  master_df$BadeaID)
master_df$CIVMID = gsub("-","_",  master_df$CIVMID)
master_df$DWI = substr( master_df$DWI , 1, 6)
sum(!is.na(master_df$DWI ))
#master_df = na.omit(master_df)
#
# index_keep_master = which( master_df$BadeaID %in% gsub("-", "_",behave$Animal ) )
# master_df = master_df[index_keep_master,]

behave$DWI =NA
behave$Age_correct = NA 
behave$Diet_correct = NA  


for (i in 1:dim(behave)[1]) {
  # index = which(behave$Animal[i] == master_df$BadeaID )
  # index = which(behave$Animal[i] == master_df$CIVMID )
  # index <- which(behave$Animal[i] == master_df$CIVMID |  master_df$BadeaID == behave$Animal[i] )
  
  #index in master
  index <- which(behave$Animal[i] == master_df$CIVMID |  master_df$BadeaID == behave$Animal[i] )
  
  #index in behave
  #matches <- behave$Animal %in% master_df$CIVMID | behave$Animal %in% master_df$BadeaID
  #index <- which(matches)
  
  if(length(index)>0 ) {
    behave$DWI[i]   = master_df$DWI[index]
    behave$Age_correct[i]   = master_df$Age_Months[index]
    behave$Diet_correct[i]   = master_df$Diet[index]
    
    
  }
  
}

genotype_counts <- table(behave$Genotype_mastersheet)
print(genotype_counts)

behave = behave[!is.na(behave$DWI),]

behave = behave [substr( behave$DWI , 1, 1) =="N", ]
table(behave$Genotype_mastersheet)

behave = behave[behave$Genotype_mastersheet!="HN",]
#HN
# index_match = match( gsub("-","_",master_df$Cardiac_ID)   , gsub("-","_", cardiac$ID )   )
#
# cardiac$Arunno = NA
# cardiac$Arunno [ index_match ] = master_df$ARunno
# cardiac = cardiac[!is.na(cardiac$Arunno),]
# 

path_connec="/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/connectome_data/"
file_list=list.files(path_connec)
plain_index = grep("_conn_sift.csv", file_list)
file_list = file_list [plain_index]









temp_conn= read.csv( paste0(path_connec,file_list[1]) , header = F )
# connectivity=array( NA ,dim=c(length(including),length(including),dim(cardiac)[1]))
connectivity=array( NA ,dim=c(dim(temp_conn)[1],dim(temp_conn)[2],dim(behave)[1]))
dim(connectivity)
notfound = 0 

for (i in 1:dim(behave)[1] ) {
  index= which(behave$DWI[i] == substr(file_list,1,6) )
  if (length(index) > 0 ) {
    #print(i)
    temp = read.csv( paste0(path_connec,file_list[index]) , header = F )
   # temp[is.na(temp)] = 0
    #temp = temp[,1:90]
    # temp = as.matrix(temp[index_include,index_include])
    temp = as.matrix(temp)
    
    #temp=(temp - mean(temp)) /sd(temp)
    connectivity[,,i]=as.matrix(temp)
    
    
  }
  else
  { notfound = c(notfound, i) }
}

# connectivity = connectivity[index_include, index_include, ]

# connectivity = connectivity [ , ,-notfound[2:length(notfound)] ]

aaa = which(is.na(connectivity) , arr.ind = T)

response=behave
save(response, file="response.rda")
save(connectivity, file="connectivity.rda")


# library(R.matlab)

# writeMat(con ="fmri_cardiac.mat" ,response=cardiac, connectivity=connectivity )


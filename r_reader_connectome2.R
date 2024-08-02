library(readxl)
library(dplyr)


cardiac_path='Full_Cardiac_Metrics_02082024_Exercise.csv'
cardiac=read.csv(cardiac_path )
cardiac = as.data.frame(cardiac)
cardiac = as.data.frame(t(na.omit(t(cardiac))))


master_path = 'MasterSheet_Experiments2021.xlsx'
master_df = read_xlsx(master_path, sheet = "18ABB11_readable02.22.22_BJ_Cor") %>%dplyr::select( ARunno, Cardiac_ID)#subselect
master_df = na.omit(master_df)



index_keep_master = which( master_df$Cardiac_ID %in% gsub("-", "_",cardiac$ID ) )
master_df = master_df[index_keep_master,]

index_match = match( gsub("-","_",master_df$Cardiac_ID)   , gsub("-","_", cardiac$ID )   )

cardiac$Arunno = NA
cardiac$Arunno [ index_match ] = master_df$ARunno 
cardiac = cardiac[!is.na(cardiac$Arunno),]


path_connec="/Users/ali/Desktop/Mar24/fmri_all_to_batch9_conn/time_ser/"
file_list=list.files(path_connec)
plain_index = grep("FC", file_list)
file_list = file_list [plain_index]







#### 324 to 332

# chass324_path = "Chass324index2.csv"
# chass324 = read.csv(chass324_path, header = F)
# index2_1000 = chass324[chass324>166]
# chass324[chass324>166] = index2_1000 - 1000 + 166

# noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab
# 
# including =c(57
#              ,91
#              ,51
#              ,64
#              ,121
#              ,65
#              ,59
#              ,19
#              ,81
#              ,25
#              ,20
#              ,72
#              ,15
#              ,16
#              ,147
#              ,62
#              ,82
#              ,33
#              ,60
#              ,42
#              ,89
#              ,66
#              ,27
#              ,116
#              ,104
#              ,17
#              ,37
#              ,32
#              ,8
#              ,69
#              ,127
#              ,73
#              ,148
#              ,28
#              ,53
#              ,34
#              ,55
#              ,44
#              ,85
#              ,29
#              ,99
#              ,40
#              ,9
#              ,119
#              ,63
#              ,114
#              ,80
#              ,78
#              ,49
#              ,38
#              ,68
#              ,120
#              ,1
#              ,41
#              ,7
#              ,12
#              ,77
#              ,26
#              ,133
#              ,3
#              ,166
#              ,113
#              ,107
#              ,11
#              ,39
#              ,36
#              ,135
#              ,21
#              ,43
#              ,10
#              ,46
#              ,54
#              ,115
#              ,124
#              ,151
#              ,100
#              ,125
#              ,144
#              ,140)
# including = setdiff(including, noreadcsf)


# 
# including = c(including, including + 166)
# 
# index_include = which( chass324$V1 %in% including)


# index_removal = which( chass324$V1 %in% noreadcsf)
# connectivity [index_removal,,]=0
# connectivity [,index_removal,]=0





temp_conn= read.csv( paste0(path_connec,file_list[1]) , header = F )
# connectivity=array( NA ,dim=c(length(including),length(including),dim(cardiac)[1]))
connectivity=array( NA ,dim=c(dim(temp_conn)[1],dim(temp_conn)[2],dim(cardiac)[1]))
dim(connectivity)
notfound = 0 

for (i in 1:dim(cardiac)[1] ) {
  index= which(cardiac$Arunno[i] == substr(file_list,4,12) )
  if (length(index) > 0 ) {
    #print(i)
    temp = read.csv( paste0(path_connec,file_list[index]) , header = F )
    temp[is.na(temp)] = 0
    #temp = temp[,1:90]
    # temp = as.matrix(temp[index_include,index_include])
    temp = as.matrix(temp)
    
    temp=(temp - mean(temp)) /sd(temp)
    connectivity[,,i]=as.matrix(temp)
    
    
  }
  else
  { notfound = c(notfound, i) }
}

# connectivity = connectivity[index_include, index_include, ]

connectivity = connectivity [ , ,-notfound[2:length(notfound)] ]
cardiac = cardiac[-notfound[2:length(notfound)],]

aaa = which(is.na(connectivity) , arr.ind = T)

response=cardiac
save(response, file="response.rda")
save(connectivity, file="connectivity.rda")


# library(R.matlab)

# writeMat(con ="fmri_cardiac.mat" ,response=cardiac, connectivity=connectivity )


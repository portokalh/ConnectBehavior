library(dplyr)

noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab

load(file='connectivity.rda')
load(file='response.rda')

# temp=connectivity[-noreadcsf,-noreadcsf,1]


for (i in 1:dim(connectivity)[3]) {
  connectivity[noreadcsf, noreadcsf,i] = 0
}



temp=connectivity[,,1]
indexlower=lower.tri(temp, diag=FALSE)
indexlowertrue=which(indexlower==TRUE)
temp=temp[indexlower]
len=sum(indexlower)  



#riskfactors=matrix(NA,  dim(response)[1], (dim(response)[2]-1))
riskfactors=response %>% select( Genotype_mastersheet , Sex_mastersheet, Age_correct, Diet_correct ,NormSWDist)


riskfactors =t(na.omit(t(riskfactors) ))
riskfactors  = as.data.frame(riskfactors)
riskfactors$HN = NA
riskfactors$HN[ grep("HN",riskfactors$Genotype_mastersheet) ] = "1"
riskfactors$HN[is.na(riskfactors$HN)] = 0 
# table(response$Genotype_mastersheet)
# table(response$Diet)


##########no HN
# indeceis_apoe_no_hn= riskfactors$Genotype %in% c("APOE22", "APOE33", "APOE44", "APOE22HN", "APOE33HN", "APOE44HN")
# riskfactors=riskfactors[indeceis_apoe_no_hn,]
# riskfactors=riskfactors%>%dplyr::select(Sex, Age_Months)
dim(connectivity)
# connectivity = connectivity [ , , indeceis_apoe_no_hn , drop =T] 


Gene=riskfactors$Genotype_mastersheet
Gene[Gene=="APOE44" |  Gene== "APOE44HN" ]=1
Gene[Gene!=1]=-1
# Gene[Gene=="E33"]=0
# Gene[Gene=="E44"]=4
# Gene[Gene=="E2HN"]=2
# Gene[Gene=="E3HN"]=3
# Gene[Gene=="E4HN"]=4
# Gene[Gene=="KO"]=0
riskfactors$Genotype_mastersheet=as.numeric(Gene)

Sex=riskfactors$Sex_mastersheet
Sex[Sex=="male"]=-1
Sex[Sex=="female"]=1
riskfactors$Sex_mastersheet=as.numeric(Sex)

Diet=riskfactors$Diet_correct
Diet[Diet=="Control"]=-1
Diet[Diet=="HFD"]=1
riskfactors$Diet_correct=as.numeric(Diet)


# 
# Exercise=riskfactors$Exercise
# Exercise[Exercise=="NO"]=-1
# Exercise[Exercise=="YES"]=1
# riskfactors$Exercise=as.numeric(Exercise)


riskfactors$Age_correct=as.numeric(riskfactors$Age_correct)

riskfactors=as.data.frame(riskfactors)
riskfactors[] <- lapply(riskfactors, as.numeric)
class(riskfactors)


riskfactors_orig = riskfactors



lm=lm(NormSWDist~Age_correct*Sex_mastersheet*Diet_correct*Genotype_mastersheet*HN, data =riskfactors )
write.csv2(file ="anova_NormSWDist_probday8.csv" ,anova(lm) )


library(ggplot2)
library(ggpubr)

my_comparisons=list(c(0,1))


HN_factor = riskfactors$HN
HN_factor[riskfactors$HN == 0 ] = "Non-HN"
HN_factor[riskfactors$HN == 1 ] = "HN"
NormSWDist = as.numeric(riskfactors$NormSWDist)


ggplot(riskfactors, aes(x = HN_factor, y = NormSWDist, color = HN_factor)) +
  geom_violin(trim = FALSE, aes(fill = HN_factor), alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("High" = "purple", "Low" = "green")) +
  scale_fill_manual(values = c("High" = "purple", "Low" = "green")) +
  labs(title = "Violin Plot of NormSWDist by HN_factor",
       x = "HN Factor",
       y = "Normalized SW Distance") +
  theme_bw()

ggsave( "HN.png" , plot= last_plot(), device='png', 
       scale=1, width=8, 
       height=8, unit=c("in") )






# riskfactors1=riskfactors%>%select(Sex, Age_Months, Genotype, Diet)
#riskfactors2=riskfactors%>%select(NormSWDist, Distance, Winding)
# riskfactors2=riskfactors%>%select( NormSWDist, Winding)

inddz1=0
for (i in 1:dim(riskfactors)[2]) if(sd(riskfactors[,i])==0 ) {inddz1=rbind(inddz1,i);  cat ( i , sd(riskfactors[,i]), "\n" );}
if (length(inddz1)>1){
  inddz1=inddz1[2:dim(inddz1)[1]]
  riskfactors=riskfactors[,-inddz1]
}



# 
# inddz2=0
# for (i in 1:dim(riskfactors2)[2]) if(sd(riskfactors2[,i])==0 ) {inddz2=rbind(inddz2,i);  cat ( i , sd(riskfactors2[,i]), "\n" );}
# if (length(inddz2)>1){
#   inddz2=inddz2[2:dim(inddz2)[1]]
#   riskfactors2=riskfactors2[,-inddz2]
# }

# temp=connectivity[-noreadcsf,-noreadcsf,1]
temp=connectivity[,,1]
indexlower=lower.tri(temp, diag=FALSE)
indexlowertrue=which(indexlower==TRUE)
temp=temp[indexlower]
len=sum(indexlower)  

image=matrix(NA,  dim(connectivity)[3], len) # -6 becasue of cfs removal

for (i in 1:dim(connectivity)[3]){
  # temp=connectivity[-noreadcsf,-noreadcsf,i]
  temp=connectivity[,,i]
  indexlower=lower.tri(temp, diag=FALSE)
  temp=temp[indexlower]
  image[i,]=temp
}
dim(image)





indd=0
for (i in 1:dim(image)[2]) if(sd(image[,i])==0 ) {indd=rbind(indd,i);  cat ( i , sd(image[,i]), "\n" );}
if (length(indd)>1){
  indd=indd[2:dim(indd)[1]]
  image=image[,-indd] }









#lets run
## Not run:
#install.packages("PMA")
#install.packages("https://gitlab.oit.duke.edu/am983/PMA2/-/archive/master/PMA2-master.tar.gz", repos = NULL, type="source")
library(PMA)

set.seed(3189) #for reproductivity

# Can run CCA with default settings, and can get e.g. 3 components
#??glmnet

#errors in this loop
# this is selecing a lambda for traits so it would penalize the most without sparisity 
for (i in 100:1) {
  zlamb=i/100
  out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb)
  numzerv=sum(out$v==0)
  #if (numzerv>0.6*length(out$v) ) { i=i+1;  zlamb=i/100;  break }
  if (numzerv>0 ) { i=i+1;  zlamb=i/100;  break }
}

numzeru=sum(out$u!=0)
xlamb=out$penaltyx
sparsity=0.95
persnonsparse=floor((1-sparsity)*length(out$u))
while (numzeru>persnonsparse) {
  xlamb=0.7*xlamb
  out2 <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
  numzeru=sum(out2$u!=0); 
}

for (i in 100:1) {
  zlamb=i/100
  out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
  numzerv=sum(out$v==0)
  if (numzerv>0 ) { i=i+1;  zlamb=i/100;  break }
}

# 
# 
# numzerv=sum(out$v!=0)
# nonsparsv=0.1
# persnonsparse=floor((1-nonsparsv)*length(out$v))
# while (numzerv>persnonsparse) {
#   zlamb=0.9*zlamb
#   out2 <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
#   numzerv=sum(out2$v!=0); 
# }





# zlamb = 0.85
# xlamb = 0.0172944

# out2 <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
out2 <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
# out2 <- CCA(x=image,z=riskfactors,typex="standard",typez="standard")
out2
sum(out2$u!=0)
out2$v

# perm.out <- CCA.permute(x=image,z=riskfactors,typex="standard",typez="standard",nperms=101, standardize=TRUE, SD=TRUE, penaltyxs = xlamb, penaltyzs = zlamb)

# if you dont want to run the 1000 permuttion again just load the 1000 permutation R data and run from here:



# print(perm.out)

# plot(perm.out)
#out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard",
#           penaltyx=perm.out$bestpenaltyx/20,penaltyz=perm.out$bestpenaltyz,
#           v=perm.out$v.init)
# out <- CCA(x=image,z=riskfactors1,typex="standard",typez="standard",
           # penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, UVperms= perm.out$UVperms, allpenaltyxs = perm.out$penaltyxs)
# out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb)

out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard", penaltyz = zlamb, penaltyx = xlamb)
#out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard")




print(out) # could do print(out,verbose=TRUE)
#print(image[out$u!=0]) 

u=out$u


#########
# u[abs(u)<quantile(u,0.999)]  = 0
nonzerou = u[u!=0]
nonzerou[abs(nonzerou)<quantile(nonzerou,0.999)]   = 0
u[u!=0] = nonzerou
sum(u!=0)
######
v=out$v


sum(u==0)
#len=length(u)
sum(u!=0)
u[u!=0]
sum(u==0)/len #sparsity 

uout=matrix(NA, dim(u)[1],1 )
#put those zeros back
# uout[indd]=0
uout=u
#uout[-indd]=coef



# indd
# indexlowertrue=which(indexlower==TRUE)
# temp[indd]
# indexlowertrue[indd]

connectivityexample=connectivity[,,1]

connectivityexample[indexlowertrue[indd]] ##yes the're them
connectivityexample[indexlowertrue[indd]]="zeros" # lest make them known for a special word
indexofzeros=which(connectivityexample=="zeros", arr.ind=TRUE)

# indexofzeros[,1]
# indexofzeros

# #lest check really quick:
# for (j in 1:dim(indexofzeros)[1]) {
# for (i in 1:dim(connectivity)[3]) {  cat(  "subject", i, "at position", indexofzeros[j,] , connectivity[matrix(c(indexofzeros[j,],i),1,3)] , "\n")
# }
# } ## yes theyre all odly zeros

#results of connectivities that matter:
nonzeroindex=which(uout!=0)

connectivityexample=connectivity[,,1]
connectivityexample[]=0
connectivitvals=connectivityexample
nonzerouout=uout[uout!=0]
for (i in 1:length(nonzeroindex)) {
  connectivityexample[indexlowertrue[nonzeroindex[i]]]=c("nonzero") # lest make them known for a special word
  connectivitvals[indexlowertrue[nonzeroindex[i]]]=nonzerouout[i] #store their coefitient values
}



library('igraph');
connectivitvalsones=connectivitvals
#####
diag(connectivitvals)=0
#####
t=which(connectivitvalsones!=0, arr.ind=TRUE)
t <- cbind(t, connectivitvals[which(connectivitvals!=0,arr.ind=TRUE)]) 
t.graph=graph.data.frame(t,directed=F)
E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'blue','red') 
#t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
minC <- rep(-Inf, vcount(t.graph))
maxC <- rep(Inf, vcount(t.graph))
minC[1] <- maxC[1] <- 0
l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                    miny=minC, maxy=maxC)      

pathnames='mouse_anatomy.csv'
datanmes=read.csv(pathnames, header = TRUE, sep = ",", quote = "")
datanmes$ROI

# noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab

#datanmes=datanmes[-noreadcsf]

# datanmess=datanmes$ROI[-noreadcsf] # remove csf
datanmess=datanmes$ROI



par(mfrow=c(1,1))

#set.vertex.attribute(t.graph, "name", value=datanmes$ROI   )


 jpeg("nets.jpeg", units="in", width=10, height=5, res=300)  

plot(t.graph, layout=l, 
     rescale=T,
     asp=0,
     edge.arrow.size=0.1, 
     vertex.label.cex=0.8, 
     vertex.label.family="Helvetica",
     vertex.label.font=4,
     #vertex.label=t.names,
     vertex.shape="circle", 
     vertex.size=5, 
     vertex.color="white",
     vertex.label.color="black", 
     #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
     edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))

 dev.off()
connectivitvals=connectivitvals+t(connectivitvals) #symetric


# nonzeroposition=which(connectivityexample=="nonzero", arr.ind=TRUE)
getwd()
filename=paste(getwd(), "/", "valandpos.mat", sep = "")
#writeMat(filename, nonzeroposition = nonzeroposition, connectivitvals = connectivitvals , oddzeroposition=indexofzeros)




subnets=groups(components(t.graph))
subnetsresults=vector(mode = "list", length = length(subnets))
colsumabs=colSums(abs(connectivitvals))
colsum=colSums(connectivitvals)

leftright=datanmes$Hemisphere




####################3
for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,8,length(temp) )
  net[1,]=as.numeric(temp)
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  #indofleftright=tt>=164
  #net[5,][indofleftright]="Right"
  #net[5,][!indofleftright]="Left"
  
  
  net[2,]=datanmess[temp]
  net[5,]=leftright[temp]
  net[1,]=paste(net[1,],net[5,])
  net[3,]= as.numeric( colsum[temp]   )
  net[4,]= as.numeric( colsumabs[temp]   )
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  subnetsresults[[i]]=net 
}

##############new net

net_new=matrix(NA, length(subnetsresults),4)


for (j in 1:dim(net_new)[1]) {
  temps=subnetsresults[[j]]
  net_new[j,1]=j
  net_new[j,2]= paste(temps[8,], collapse = ", ")
  net_new[j,3] = paste(paste(temps[5,],temps[2,]), collapse = ", ")
  net_new[j,4] = paste(temps[7,1])
  
}
colnames(net_new)=c("Sub-Network", "Region Number", "Region Name", "Sub-Network Weight")


#install.packages("xlsx")
library(xlsx)


for (i in 1:length(subnetsresults)){
  net=t(subnetsresults[[i]])
  write.xlsx2(net, "nets.xlsx", sheetName =  paste0(i), append=TRUE )
}

write.xlsx2(net_new, "net_new.xlsx" )


# install.packages("vioplot")
library("vioplot")


for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,8,length(temp) )
  net[2,]=datanmess[temp]
  net[1,]=as.numeric(temp)
  net[3,]= as.numeric( colsumabs[temp]   )
  net[4,]= as.numeric( colsum[temp]   )
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  indofleftright=tt>=164
  net[5,][indofleftright]="Right"
  net[5,][!indofleftright]="Left"
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  subnetsresults[[i]]=net
}





for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,8,length(temp) )
  net[2,]=datanmess[temp]
  net[1,]=as.numeric(temp)
  net[3,]= as.numeric( colsumabs[temp]   )
  net[4,]= as.numeric( colsum[temp]   )
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  indofleftright=tt>=164
  net[5,][indofleftright]="Right"
  net[5,][!indofleftright]="Left"
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  for (j in 1:length( net[8,])) {
    tempindex=which(datanmes$ROI %in% net[2,j]  )
    if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
  }
  subnetsresults[[i]]=net 
}



#for (i in 1:length(subnetsresults)) {
#  net=subnetsresults[i]
#  print(net[[1]][1:2,])
#}


for (i in 1:length(subnetsresults)) {
  net=subnetsresults[i]
  cat( i,'th sub-net: the summation of all edges in this sub-net is' ,sum(as.numeric(net[[1]][4,])), 'and summation of absolut values of all edges in this subnet is', sum(as.numeric(net[[1]][3,])),'\n')
  cat(  'the fsirst row is the Region #, second row is the name of Region, the third row is the sum of absulote values of the edges of each region, and the last row is the sum of edges of each region \n')
  print(net)
  cat( '\n \n \n')
}


capture.output(subnetsresults, file = "subnet.txt")


write.csv(subnetsresults, row.names = T)




#install.packages("xlsx")
library(xlsx)

# 
# for (i in 1:length(subnetsresults)){
#   net=subnetsresults[[i]]
#   write.xlsx2(net, "ssir.xlsx", sheetName =  paste0(i), append=TRUE )
# }


# install.packages("vioplot")
library("vioplot")





# path3="/Users/ali/Desktop/mar/mice/apoe234/resultsresponse_arrayraw.mat"
# data3=readMat(path3)
# responseraw=data3$response.arrayraw
# 
# riskfactors=matrix(NA,  dim(responseraw)[1], (dim(responseraw)[2]-1))
#sum(riskfactors[,2]==3)
riskfactors=response

# #subjnameofconnectivity=data$subjlist
# 
# for (i in 1:dim(riskfactors)[1]) {
#   ind=which(response[i,1]==subjnameofconnectivity)
#   #if (i!=ind) cat("here", i)
#   riskfactors[ind,]=responseraw[ind, 2:(dim(response)[2])     ]
# }




####################
#############
########3
# 
# ### histograms of nets
# histdata=matrix(0,length(subnetsresults),dim(connectivity)[3])
# #t
# 
# for (j in 1:length(subnetsresults)){
#   net=subnetsresults[[j]]
#   subnetsuperset=as.numeric(net[1,])
#   for (i in 1:dim(t)[1]){
#     if ( t[i,][1]%in%subnetsuperset){
#       for (k in 1:dim(connectivity)[3]) {
#         temp=connectivity[,,k]
#         histdata[j,k]=histdata[j,k]+ temp[t[i,][1],t[i,][2]]+temp[t[i,][2],t[i,][1]]
#       }
#     }
#   }
# }
# histdata=cbind(seq(1,length(subnetsresults)),histdata)
# 
# 
# 
# ##split plots.
# histdatasplit=histdata[,2:dim(histdata)[2]]
# #uniqgeneafterreg=unique(riskfactors[,5])
# apoe2=histdatasplit[,riskfactors$Genotype==2] 
# apoe3=histdatasplit[,riskfactors$Genotype==3]
# apoe4=histdatasplit[,riskfactors$Genotype==4]
# 
# library(plotrix)
# GenoTypes=as.factor(riskfactors$Genotype)
# #GenoTypes[GenoTypes==2]="APOE22";GenoTypes[GenoTypes==3]="APOE33";GenoTypes[GenoTypes==4]="APOE44";
# Sex=riskfactors$Sex;
# Sexname=Sex; Sexname[Sex=="female"]="2 Female"; Sexname[Sex=="male"]="1 Male" 
# Diet=as.factor(riskfactors$Diet)
# #Dietname=Diet; Diet[Diet==1]="Control"; Diet[Diet==2]="HFD";
# #Sexname=Sex; Sexname[Sexname=="1M"]="Male" ;  Sexname[Sexname=="2F"]="Female";
# 
# ###### here specify   
# xaxisvar=Sexname; 
# xaxis="Sex" ;
# xaxisvarnames=Sexname;
# brightnessvar=GenoTypes; 
# brightness="GenoTypes";
# #######
# 
# 
# names(xaxisvar)=xaxisvarnames
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# jpeg("violin.jpeg", units="in", width=20, height=10, res=300)  
# 
# sqrt=sqrt(length(subnetsresults))
# if (sqrt>4) sqrt=3
# par(mfrow = c(ceiling(sqrt), ceiling(length(subnetsresults)/sqrt)))
# for (j in 1:length(subnetsresults)){
#   #cols <- brewer.pal(8,'Set2')[6:8]
#   cols=c( 'chartreuse1','blueviolet')
#   colsorig=adjustcolor(cols, alpha.f = 0.5)
#   colsorig2=cols
#   # cols[length(unique(xaxisvar))]=rgb(1,0,0,)
#   cols=adjustcolor(cols, alpha.f = 0.25)
#   legend=paste0("Net ",j, ". Networks medians: ( ");
#   for (m in 1:length(unique(xaxisvar))) {
#     if (m==1) {legend=paste(legend,  round(median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[m]]),2), sep =""  )}
#     else {legend=paste(legend,  round(median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[m]]),2), sep =","  )}
#   }
#   legend=paste(legend,")")
#   vioplot(histdatasplit[j,]~ xaxisvar , plotCentre = "dot", col =cols,  ylab="net weight" ,
#           main = legend, xlab ="" , cex.axis=1.5, cex.lab=2,
#           rectCol=colsorig, lineCol=colsorig)
#   a=subnetsresults[[j]]; a=as.table(a);a= as.numeric(a[1,]); 
#   #mtext(paste("R",a, collapse=', '), side=4, cex=0.5)
#   
#   for (l in 1:length(sort(unique(brightnessvar), decreasing=T))) {
#     cols=c( 'chartreuse1','blueviolet')
#     cols[length(unique(xaxisvar))]=rgb(1,0,0,)
#     tempnum=length(sort(unique(brightnessvar), decreasing=T))
#     #cols=adjustcolor(cols, alpha.f = l/length(sort(unique(brightnessvar), decreasing=T)))
#     cols=c( 'chartreuse1','blueviolet')
#     stripchart(histdatasplit[j,brightnessvar==sort(unique(brightnessvar), decreasing=T)[l]]~xaxisvar[brightnessvar==sort(unique(brightnessvar), decreasing=T)[l]], vertical = TRUE, method = "jitter", points=50,
#                pch = (17:length(sort(unique(brightnessvar), decreasing=T))), add = TRUE, col =cols , offset=0, cex = 1.2)
#     cat(sort(unique(brightnessvar), decreasing=T)[l],"  "  ,l/length(unique(brightnessvar, "\n") )  )
#   }
#   
#   # mtext(paste("Dark Jitt=", sort(unique(brightnessvar), decreasing=T)[l], brightness) , side=1, cex=0.6)
#   
#   for (k in 1:length(unique(xaxisvar))) {
#     ablineclip(h=median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[k]]), col=colsorig[k], lwd = 2, x1=0.4, x2=k, lty="dotted", cex=2)
#     
#   }
# }
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# #####3 only some violin 
# c=c(2,7,9,11,15,19)
# 
# jpeg("violin_selected.jpeg", units="in", width=20, height=10, res=300)  
# 
# par(mfrow = c(2,3))
# for (j in c){
#   #cols <- brewer.pal(8,'Set2')[6:8]
#   cols=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0,))
#   cols[length(unique(xaxisvar))]=rgb(1,0,0,)
#   cols=adjustcolor(cols, alpha.f = 0.1)
#   legend=paste0("Net ",j, ". Networks medians: ( ");
#   for (m in 1:length(unique(xaxisvar))) {
#     if (m==1) {legend=paste(legend,  round(median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[m]]),2), sep =""  )}
#     else {legend=paste(legend,  round(median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[m]]),2), sep =","  )}
#   }
#   legend=paste(legend,")")
#   vioplot(histdatasplit[j,]~ xaxisvar , plotCentre = "dot", col =cols,  ylab="net weight" ,
#           main = legend, xlab ="" )
#   a=subnetsresults[[j]]; a=as.table(a);a= as.numeric(a[1,]); 
#   #mtext(paste("R",a, collapse=', '), side=4, cex=0.5)
#   
#   for (l in 1:length(sort(unique(brightnessvar), decreasing=T))) {
#     cols=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0,))
#     cols[length(unique(xaxisvar))]=rgb(1,0,0,)
#     tempnum=length(sort(unique(brightnessvar), decreasing=T))
#     cols=adjustcolor(cols, alpha.f = l/length(sort(unique(brightnessvar), decreasing=T)))
#     stripchart(histdatasplit[j,brightnessvar==sort(unique(brightnessvar), decreasing=T)[l]]~xaxisvar[brightnessvar==sort(unique(brightnessvar), decreasing=T)[l]], vertical = TRUE, method = "jitter", points=50,
#                pch = (17:length(sort(unique(brightnessvar), decreasing=T))), add = TRUE, col =cols , offset=0, cex = 1.2)
#     cat(sort(unique(brightnessvar), decreasing=T)[l],"  "  ,l/length(unique(brightnessvar, "\n") )  )
#   }
#   
#   #mtext(paste("Dark Jitt=", sort(unique(brightnessvar), decreasing=T)[l], brightness) , side=1, cex=0.6)
#   
#   for (k in 1:length(unique(xaxisvar))) {
#     ablineclip(h=median(histdatasplit[j,xaxisvar==sort(unique(xaxisvar))[k]]), col=cols[k], lwd = 1, x1=0.4, x2=k, lty="dotted")
#     
#   }
# }
# 
# dev.off()
# 


##################################3
##############################3
############################3
####### making graph match the original atlas with csf

# 
# library('igraph');
# #connectivitvalsones=connectivitvals
# t=which(connectivitvalsones!=0, arr.ind=TRUE)
# 
# for (i in 1:(dim(t)[1]*dim(t)[2])) {
#   for (j in noreadcsf) {
#     if (j <= t[i]) {
#       t[i] = t[i] + 1
#     }
#   }
#   #t[i] = t[i] - 1
# }
# 
# t <- cbind(t, connectivitvals[which(connectivitvals!=0,arr.ind=TRUE)]) 
# t.graph=graph.data.frame(t,directed=F)
# E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'blue','red') 
# #t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
# minC <- rep(-Inf, vcount(t.graph))
# maxC <- rep(Inf, vcount(t.graph))
# minC[1] <- maxC[1] <- 0
# l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
#                     miny=minC, maxy=maxC)      
# 
# pathnames='/Users/ali/Desktop/Jul/apoe/mouse_anatomy.csv'
# datanmes=read.csv(pathnames, header = TRUE, sep = ",", quote = "")
# datanmes$ROI
# 
# noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab
# 
# #datanmes=datanmes[-noreadcsf]
# 
# datanmess=datanmes$ROI[-noreadcsf] # remove csf
# #datanmess=datanmes$ROI
# 
# 
# 
# par(mfrow=c(1,1))
# 
# #set.vertex.attribute(t.graph, "name", value=datanmes$ROI   )
# 
# 
# jpeg("nets2.jpeg", units="in", width=15, height=7, res=300)  
# 
# plot(t.graph, layout=l, 
#      rescale=T,
#      asp=0,
#      edge.arrow.size=0.1, 
#      vertex.label.cex=0.8, 
#      vertex.label.family="Helvetica",
#      vertex.label.font=3.5,
#      #vertex.label=t.names,
#      vertex.shape="circle", 
#      vertex.size=4.5, 
#      vertex.color="white",
#      vertex.label.color="black", 
#      #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
#      edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))
# 
# dev.off()
# 
# 


###compare behavior traits with top 50 edges


#Top 50 edges attaches as a 332x332x7 array with slices corresponding to age,  sex,  diet, 22->44, 22->33, 33->44, and HN, respectively.

# risk_nets=readRDS('/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/steven_results/mouse/tensor_top50.rds')
# 
# 
# # intersect the risk nets with connectivity 
# 
# for (i in 1:dim(risk_nets)[3]) {
#   
#   risk_net_one = risk_nets[,,i]
#   risk_net_one_1_0 =risk_net_one
#   risk_net_one_1_0[risk_net_one_1_0!=0] =1
#   connectivitvals_1_0 = connectivitvals
#   connectivitvals_1_0[connectivitvals_1_0!=0]=1
#   result = risk_net_one_1_0*connectivitvals_1_0
#   print(which(result!=0))
#   
}

# "age"  "sex"  "diet" "2244" "2233" "3344" "hn"  
res=readRDS('/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/steven_results/risk_full_networks/features_correct.rds')

res2=readRDS('/Users/alex/AlexBadea_MyCodes/DTI2Behavior_Mouse/steven_results/risk_full_networks/tensor_full_combined.rds')
 
for (ii in 1:dim(res2)[3]) {
  
  risk_net_one = res2[,,ii]
  risk_net_one_1_0 =risk_net_one
  risk_net_one_1_0[risk_net_one_1_0!=0] =1
  connectivitvals_1_0 = connectivitvals
  connectivitvals_1_0[connectivitvals_1_0!=0]=1
  result = risk_net_one_1_0*connectivitvals_1_0
   print(sum(result!=0))
# print( min(risk_net_one*result) )
  
  weighted_result = risk_net_one*result
  
  
  
  
  t=which(weighted_result!=0, arr.ind=TRUE)
  t <- cbind(t, weighted_result[which(weighted_result!=0,arr.ind=TRUE)]) 
  t.graph=graph.data.frame(t,directed=F)
  E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'blue','red') 
  #t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
  minC <- rep(-Inf, vcount(t.graph))
  maxC <- rep(Inf, vcount(t.graph))
  minC[1] <- maxC[1] <- 0
  l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                      miny=minC, maxy=maxC)      
  
  pathnames='mouse_anatomy.csv'
  datanmes=read.csv(pathnames, header = TRUE, sep = ",", quote = "")

  
  subnets=groups(components(t.graph))
  subnetsresults=vector(mode = "list", length = length(subnets))
  colsumabs=colSums(abs(weighted_result))
  colsum=colSums(weighted_result)
  
  leftright=datanmes$Hemisphere
  
  
  
  ####################3
  for (i in 1:length(subnets)) {
    temp=subnets[[i]]
    temp=as.numeric(temp)
    net=matrix(NA,8,length(temp) )
    net[1,]=as.numeric(temp)
    tt=as.numeric(net[1,])
    #tt=c(1,200)
    #indofleftright=tt>=164
    #net[5,][indofleftright]="Right"
    #net[5,][!indofleftright]="Left"
    
    
    net[2,]=datanmess[temp]
    net[5,]=leftright[temp]
    net[1,]=paste(net[1,],net[5,])
    net[3,]= as.numeric( colsum[temp]   )
    net[4,]= as.numeric( colsumabs[temp]   )
    net[6,]=sum(as.numeric(net[4,]))
    net[7,]=sum(as.numeric(net[3,]))
    for (j in 1:length( net[8,])) {
      tempindex=which(datanmes$ROI %in% net[2,j]  )
      if (net[5,j]=="Right" ) {net[8,j]= max(tempindex) } else { net[8,j]=min(tempindex) }
    }
    subnetsresults[[i]]=net 
  }
  
  ##############new net
  
  net_new=matrix(NA, length(subnetsresults),4)
  
  
  for (j in 1:dim(net_new)[1]) {
    temps=subnetsresults[[j]]
    net_new[j,1]=j
    net_new[j,2]= paste(temps[8,], collapse = ", ")
    net_new[j,3] = paste(paste(temps[5,],temps[2,]), collapse = ", ")
    net_new[j,4] = paste(temps[7,1])
    
  }
  colnames(net_new)=c("Sub-Network", "Region Number", "Region Name", "Sub-Network Weight")
  
  
  #install.packages("xlsx")
  library(xlsx)
  
  
  
  write.xlsx2(net_new, paste0("common_with_",res[ii] ,".xlsx") )
  
  
  
  
}



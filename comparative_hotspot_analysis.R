####################################################################
####################################################################
#
#   HOTSPOT ANAYLSIS
#     2 parametric models (persistence, conditional on presence)
#     1 spatial model (Getis-Ord Gi*)
#       (2nd spatial model KDE done in arcgis)
#
#   COMPARATIVE ANALYSIS
#     Percent congruence
#     Pearson's correlation
#
#   Runtime: ~16 minutes
#
#   Authors: Sussman A.L., Gardner B., Adams E.M., Salas L., Kenow K.P., 
#   Luukkonen D.R., Monfils M.J., Mueller W.P., Williams K.A., Leduc-Lapierre M., and Zipkin E.F
#
####################################################################
####################################################################

#set working directory
setwd("C:/Users/")

#read in standardized data
blocks_surv <- read.csv("data_blocks.csv")                                          #distinct block list
bird_recs_persist <- read.csv("data_birds_persistence.csv")                         #bird data for persistence model
bird_recs <- read.csv("data_birds.csv")                                             #bird data for conditional on presence & Gi* models (also used for kde in ArcGIS)

#set species code for selected species/group
species = 'ALLSP'

####### KERNEL DENSITY ESTIMATION MODEL (Spatial Approach #1) ######
#
#   Sources: 
#     Wilson L.J., McSorley, C.A., Gray C.M., Dean B.J., Dunn T.E., Webb A., and Reid J.B. 2009. 
#       Radio-telemetry as a tool to define protected areas for seabirds in the marine environment. Biological Conservation 142:1808-1817.
#     O'Brien S.H., Webb A., Brewer M.J., and Reid J.B. 2012. Use of kernel density estimation and maximum curvature to set 
#       Marine Protected Area boundaries: Identifying a Special Protection Area for wintering red-throated divers in the UK. 
#       Biological Conservation 156:15-21.
#     Suryan R.M., Santora J.A., and Veit R.R. 2012. New approach for using remotely sensed chlorophyll a to 
#       identify seabird hotspots. Marine Ecology Progress Series 451:213-225.
#     Wong S.N.P., Gjerdrum C., Morgan K.H., and Mallory M.L. 2014. Hotspots in cold seas: The composition, distribution, 
#       and abundance of marine birds in the North American Arctic. Journal of Geophysical Research: Oceans 119:1691-1705.
#     
####################################################################

# this model is done in ArcGIS using the kernel density tool, as at the time of this  
# analysis R did not have the capability of conducting KDE with weighting (e.g., using abundance or density values)


#### GETIS-ORD Gi* HOTSPOT ANALYSIS MODEL (Spatial Approach #2) ####
#
#   Sources:  
#     Getis A. and Boots B. 1978. Models of spatial processes: an approach to the study of point, line, and area patterns. 
#       Cambridge University Press, Cambridge, England. 198pp.
#     Getis A. and Ord J.K. 1992. The analysis of spatial association by use of distance statistics. Geographical Analysis 24(3):189-206.
#     Santora J.A., Reiss C.S., Loeb V.J., and Veit R.R. 2010. Spatial association between hotspots of baleen whales and demographic 
#       patterns of Antarctic krill Euphasia superba suggests size-dependent predation. Marine Ecology Progress Series 405:255-269.
#     Kuletz K.J., Ferguson M.C., Hurley B., Gall A.E., Labunski E.A., and Morgan T.C. 2015. Seasonal spatial patterns in seabird and 
#       marine mammal distribution in the eastern Chukchi and western Beaufort seas: Identifying biologically important areas. 
#       Progress in Oceanography 136:175-200.
#
####################################################################

#load libraries
library(maptools)
library(rgdal)
library(spdep)
library(plyr)
library(data.table)

#import shapefile of grid cells with standardized effort-corrected counts 
      #(could also use grid cell layer joined to bird data and convernt to spatial feature...)
abund_shp<-readShapePoly("C:/Users/ALLSP_poly.shp")

#set projection (important for maintaining km in transects)
gl_albers<-CRS("+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-84.455955 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
proj4string(abund_shp) <- gl_albers

#merge shapefile with block list to eliminate all blocks with transects <1km
g.data<-merge(x = abund_shp, y = blocks_surv, by = "grid_id" )
g.data<-g.data[which(g.data$trans_len>=1),]

#calculate contiguity based neighbors (i.e., sharing boundary)
data_nbr<-poly2nb(g.data,queen=FALSE)                                                   #rook method
data_nbr<-include.self(data_nbr)                                                        #include self 

#create row standardized spatial weights matrix, "listw" object (default style = "W")
nbr_w<- nb2listw(data_nbr)

#run Getis-Ord G-test (uses the data column and the spatial weights matrix)
G <- localG(g.data$abun_std, nbr_w)

#convert localG object to dataframe
t.df <- do.call("rbind", lapply(G, as.data.frame))
row.names(t.df)<-g.data$grid_id

#do some calculations on zscores to determine which level of confidence you are: 90 , 95 , 99 (zscore: 2.575, 1.96, 1.645); positive = hotspot, negative = coldspot
t.df$bin<- ifelse(t.df$'X[[i]]' <=-2.575, "-3",
                  ifelse( (t.df$'X[[i]]' >-2.575 & t.df$'X[[i]]' <= -1.96) , "-2",
                          ifelse( (t.df$'X[[i]]' >-1.96 & t.df$'X[[i]]' <= -1.645) , "-1",
                                  ifelse( (t.df$'X[[i]]' >-1.645 & t.df$'X[[i]]' < 1.645) , "0",
                                          ifelse( (t.df$'X[[i]]' >=1.645 & t.df$'X[[i]]' <1.96) , "1",
                                                  ifelse( (t.df$'X[[i]]' >=1.96 & t.df$'X[[i]]' <2.575) , "2",
                                                          ifelse(t.df$'X[[i]]'>=2.575, "3",
                                                                 NA  ))))))) # all other values map to NA

#convert back to dataframe and merge with orig data (data.merge) to get block counts and num times surveyed included in output
setDT(t.df, keep.rownames = TRUE)[]
t.df$grid_id <- as.numeric(t.df$rn)
colnames(t.df) <- c("blockid","zscore","bin","grid_id")
t.merg<-merge(x = t.df, y = g.data, by= "grid_id")
gstat.hs<-as.data.frame(t.merg)
row.names(gstat.hs)<-gstat.hs$grid_id

#select specific columns for output
gstat.hs<-gstat.hs[,c(1,3:5,15,18:20,25)]

#export model output
gstat_out<-paste(species,"_gstat.csv",sep="", collapse = "")
write.csv(gstat.hs, (gstat_out))  

############ PERSISTENCE MODEL (Parametric Approach #1) ############
#
#     Sources:       
#       Suryan R.M., Santora J.A., and Veit R.R. 2012. 
#         New approach for using remotely sensed chlorophyll a to identify seabird hotspots. Marine Ecology Progress Series 451:213-225.
#       Santora J.A. and Veit R.R. 2013. Spatio-temporal persistence of top predator hotspots near the Antarctic Peninsula. 
#         Marine Ecology Progress Series 487:287-304.
#       Johnson S.M., Connelly E.E., Williams K.A., Adams E.M., Stenhouse I.J., and Gilbert A.T. 2015. 
#         Integrating data across survey methods to identify spatial and temporal patterns in wildlife distributions in: 
#         Wildlife densities and habitat use across temporal and spatial scales on the Mid- Atlantic Outer Continental Shelf: 
#         Final report to the Department of Energy Efficiency and Renewable Energy Wind & Water Power Technologies Office. 
#         Williams K.A., Connelly E.E., Johnson S.M., and Stenhouse I.J. (eds.) Award Number: DE-EE0005362. 
#         Report BRI 2015-11, Biodiversity Research Institute, Portland, Maine. 56 pp.
#
####################################################################

#create function that loops through all survey data
persist.hs.function<-function(data){
  require(fitdistrplus)                                                                 #pull in library to fit gamma distribution to data
  all<-data                                                                             #read in data
  survey_prob<-matrix(NA,nrow=nrow(all))                                                #make blank matrix to store survey specific hotspots
  all_prob<-matrix(NA,nrow=nrow(all))                                                   #make blank matrix to store all survey hotspots
  
  #get probability of being a hotspot for a unique survey
  for (i in 1:nrow(all)){
    surv<-all[i,]$surveyid                                                              #this is a factor of distinct survey IDs
    tmp<-all[which(all$surveyid==surv),]                                                #pull out all the data from a given survey for this point
    if(nrow(tmp)>1){
      tmp.g<-fitdist(tmp$block_survey_count,"gamma",method="mme")                       #fit a gamma distribution to the survey data
      survey_prob[i]<-pgamma(all[i,9],shape=tmp.g$estimate[1],rate=tmp.g$estimate[2])   #assign a quantile to the data point given the gamma distribution we just fit
    } 
  }
  
  #get probability of being a hotspot for ALL surveys
  for (i in 1:nrow(all)){
    all.g<-fitdist(all$block_survey_count,"gamma",method="mme")                         #fit a gamma distribution for all data
    all_prob[i]<-pgamma(all[i,9],shape=all.g$estimate[1],rate=all.g$estimate[2])        #assign a quantile to the data point given the gamma distribution for all data
  }
  
  combo_prob<-(survey_prob+all_prob)/2                                                  #combine the two probabilities together, essentially overall hotspots plus survey-specific hotspots (averages the two)
  all.final<-data.frame(all,survey_prob,all_prob,combo_prob)                            #put all new calculated fields into dataframe with original data
  plot(all.g)                                                                           #plot the cumulative density frequencies, QQ plot, PP plot
  dev.off()                                                                             #turn off plot
  return(all.final)
}

#merge bird data with surveyed block list and remove all blocks that have less than 1km transects
persist.data.blocks<-merge(x = bird_recs_persist, y = blocks_surv, by = "grid_id" )
persist.data.blocks<-persist.data.blocks[which(persist.data.blocks$trans_len>=1),]

#run function on bird data
persist.hs<-persist.hs.function(data=persist.data.blocks)
persist.hs[is.na(persist.hs)] <- 0                                                      #remove NAs

#export model output                                                                  
persist_out<-paste(species,"_persistence.csv",sep="", collapse = "")
write.csv(persist.hs, (persist_out)) 

###### CONDITIONAL ON PRESENCE MODEL (Parametric Approach #2) ######
#
#   Sources:
#     Kinlan et al. 2012, Zipkin et al. 2015
#     Kinlan, B.P., E.F. Zipkin, A.F. O'Connell, and C. Caldow. 2012. 
#       Statistical analyses to support guidelines for marine avian sampling: 
#       final report. U.S. Department of the Interior, Bureau of Ocean Energy Management, 
#       Office of Renewable Energy Programs, Herndon, VA. OCS Study BOEM 2012-101. NOAA Technical Memorandum NOS NCCOS 158. xiv+77 pp.
#     Zipkin E.F., Kinlan B.P., Sussman A., Rypkema D., Wimer and M., O'Connell A.F. 2015. 
#       Statistical guidelines for assessing marine avian hotspots and coldspots: 
#       A case study on wind energy development in the U.S. Atlantic Ocean. Biological Conservation 191: 216-223.
#
####################################################################

#load libraries
library(fitdistrplus)

#aggregate on unique blocks using the total standardized effort-corrected count
cond.blocks<-with(bird_recs, aggregate(bird_recs[,12], data.frame(grid_id), sum))       
colnames(cond.blocks) <- c("grid_id",species)
data.merge<-merge(x = cond.blocks, y = blocks_surv, by = "grid_id" )                    #merge with surveyed blocks
data.merge$stdabund<-data.merge[,species]/data.merge$num_times_surveyed                 #add column to standardize for effort (number of times surveyed)
data.merge<-data.merge[which(data.merge$trans_len>=1),]                                 #cutout all blocks that have less than 1km transect

#pull out all records for the species
x=data.merge$stdabund
names(x)=data.merge$grid_id
n=length(x)
#pull out all the nonzero counts
x1=x[which(x>0)]
n1=length(x1)

#create a vector with the number of times that a block was surveyed, including times when no individuals were observed
times.surv = data.merge$num_times_surveyed
names(times.surv)=data.merge$grid_id

#create a vector with the mean count of each block
species.mean = rep(0,length(times.surv))
names(species.mean)=names(times.surv)
#fill in the vector with the mean values for each block
for(k in 1:length(species.mean)) {
  a=which(names(x)==names(times.surv[k]))
  species.mean[k]=data.merge$stdabund[a]
}

#find the prevalence
prev=n1/n

#num simulations
nsim=10000

#fit lognormal distribution to nonzero counts
fit.lnorm <- fitdist(x1, "lnorm")

#pull out the parameter estimates
ln.mean = fit.lnorm$estimate[1]
ln.sd = fit.lnorm$estimate[2]

#simulate data for lnorm distribution given lognormal distribution mean and sd on count data
rn.binom <- rbinom(n*nsim,1,prev)                                                       #binom = prob of seeing a bird (yes/no)
binmat<-matrix(rn.binom,nrow=n,ncol=nsim)
rn.lnorm <- rlnorm(n*nsim,ln.mean,ln.sd)                                                #lnorm = okay, we see a bird, what's our species count?
lnmat<-matrix(rn.lnorm, nrow=n, ncol=nsim)
#combine into matrix: 1767 rows (num blocks surveyed) by 10,000 cols (num simulations)
y=binmat*lnmat

#estimate the hotspot value 
ln.hotpvalue.3xmean=rep(NA,length(times.surv))                                          #create empty vector
for(j in 1:length(times.surv)) {
  ln.hotpvalue.3xmean[j]=length(which((3*apply(y,2,mean))<species.mean[j]))/dim(y)[2]
}                                                                                       #loop through all surveyed blocks (~15min)
#get the 3x mean value for each of the 10,000 columns in the simulated data y
#how many of those y values are smaller than the true mean? (divide by num simulations to return to scale)
    #a value that is close to 0 indicates that the location is not a hotspot
    #a value that is close to 1 indicates that the location is significantly a hotspot

#combine results & export model output
conditional.hs=cbind(as.numeric(names(times.surv)),times.surv,species.mean, ln.hotpvalue.3xmean)
rownames(conditional.hs)=rep(names(times.surv))
colnames(conditional.hs)=c("grid_id","times_sampled","mean_count", "hotpvalue_3xmean")
conditional.hs<-as.data.frame(conditional.hs)
conditional.hs_out<-paste(species,"_conditional.csv",sep="", collapse = "")
write.csv(conditional.hs, (conditional.hs_out))  

####################### COMPARATIVE ANALYSIS #######################
#
#   Format and combine results, look at percent congruence
#   and run pearson's correlations on all four hotspot models
#
####################################################################

#load libraries
library(limma)
library(psych)

################### PREPARE HOTSPOT MODEL RESULTS ##################
#read in kernel density estimation hotspot analysis data (from ArcGIS)
kde.hs<-read.table("ALLSP_kde.csv", sep=",",header=TRUE)
names(kde.hs) <- tolower(names(kde.hs))                                                 #change column names to lower case
names(kde.hs)[names(kde.hs)=="gid"] <- "grid_id"                                        #set row names = blocks
kde.hs<-merge(x=blocks_surv,y=kde.hs,by="grid_id",all.x=TRUE)                           #merge with blocks list
kde.hs[is.na(kde.hs)] <- 0                                                              #set NAs to zero

#Gi*
names(gstat.hs)[names(gstat.hs)=="X"] <- "grid_id"                                      #change grid col name to grid_id to match others
gstat.hs[is.na(gstat.hs)] <- 0                                                          #set NAs to zero

#read in hotspot persistence (persistence calculated in PostgreSQL)
persist.hs.tot <- read.csv("ALLSP_persistence_totals.csv",header=TRUE)
persist.hs.tot[is.na(persist.hs.tot)] <- 0                                              #set NAs to zero

#conditional on presence
conditional.hs.b<-merge(x=blocks_surv,y=conditional.hs,by="grid_id",all=TRUE)           #merge with blocks list
conditional.hs.b[is.na(conditional.hs.b)] <- 0                                          #set NAs to zero

################ COMBINE RESULTS INTO ONE DATAFRAME ################
#create empty df, and join each model's hotspot data 
all.values<-data.frame(matrix("NA", ncol = 6, nrow = 1767))
colnames(all.values)<-c("grid_id","num_times_surv","kde","gstat","persistence","conditional")
all.values$grid_id<-blocks_surv$grid_id
all.values$num_times_surv<-blocks_surv$num_times_surveyed[match(all.values$grid_id, blocks_surv$grid_id)]
row.names(all.values)<-all.values$grid_id
#populate models' columns
all.values$persistence<-persist.hs.tot$prop_hotspot_75[match(all.values$grid_id, persist.hs.tot$grid_id)]
all.values$conditional<-conditional.hs.b$hotpvalue_3xmean[match(all.values$grid_id, conditional.hs.b$grid_id)]
all.values$kde<-kde.hs$mean[match(all.values$grid_id, kde.hs$grid_id)]
all.values$gstat<-gstat.hs$zscore[match(all.values$grid_id, gstat.hs$grid_id)]

#we only want to run rank/correlations on blocks with 4 or more sampling events
all.values<-all.values[which(all.values$num_times_surv>=4),]

####################### LOOK AT % CONGRUENCE #######################
#determine if value/grid cell falls within the top 25% (1) or not (0)
all.values$kde.bin<-ifelse(all.values$kde <= quantile(all.values$kde,c(0.75)), 0,1)
all.values$gstat.bin<-ifelse(all.values$gstat <= quantile(all.values$gstat,c(0.75)), 0,1)
all.values$persistence.bin<-ifelse(all.values$persistence <= quantile(all.values$persistence,c(0.75)), 0,1)
all.values$conditional.bin<-ifelse(all.values$conditional <= quantile(all.values$conditional,c(0.75)), 0,1)
c<-all.values[,7:10]    #combine into one df

#create venn diagram to look at percent congruence
meth.venn<- vennCounts(c)
vennDiagram(meth.venn,names=c("kde","getis-ord gi*", "persistence","conditional"),circle.col =c("skyblue", "pink1", "mediumorchid", "orange") )

##################### RUN PEARSON'S CORRELATION ####################
#create new df for each model and add model.75 column --> selects all raw values above 75th percentile, everything else gets 0 (i.e., taking the top 25%)

#kde
kde_top25<-transform(all.values, kde.75 = ifelse(all.values$kde <= quantile(all.values$kde,0.75), 0, all.values$kde))
#select only relevant cols
kde_top25<-kde_top25[,c(2,1,3,7,11)]
names(kde_top25)[names(kde_top25)=="grid_id"] <- "kde_grid_id"

#Gi*
gstat_top25<-transform(all.values, gstat.75 = ifelse(all.values$gstat <= quantile(all.values$gstat,0.75),0, all.values$gstat))
#select only relevant cols
gstat_top25<-gstat_top25[,c(1,4,8,11)]
names(gstat_top25)[names(gstat_top25)=="grid_id"] <- "gstat_grid_id"

#persistence
persistence_top25<-transform(all.values, persistence.75 = ifelse(all.values$persistence <= quantile(all.values$persistence,0.75), 0, all.values$persistence))
#select only relevant cols
persistence_top25<-persistence_top25[,c(1,5,9,11)]
names(persistence_top25)[names(persistence_top25)=="grid_id"] <- "persistence_grid_id"

#conditional
conditional_top25<-transform(all.values, conditional.75 = ifelse(all.values$conditional <= quantile(all.values$conditional,0.75), 0, all.values$conditional))
#select only relevant cols
conditional_top25<-conditional_top25[,c(1,6,10,11)]
names(conditional_top25)[names(conditional_top25)=="grid_id"] <- "conditional_grid_id"

#combine all models together
allmethods_top25=cbind(kde_top25,gstat_top25,persistence_top25,conditional_top25)

#pull out four 75th percentile cols for correlations
allmethods_top25_corr_rawvals<-allmethods_top25[,c(3,7,11,15)]
#run correlations
cor.four.top25.raw<-corr.test(allmethods_top25_corr_rawvals,use = "complete",method="pearson",adjust="bonferroni", alpha=.05,ci=TRUE)
#export correlations & all model values for mapping
four.top25.raw.out<-paste(species,"_vals.csv",sep="", collapse = "")
write.csv(allmethods_top25, (four.top25.raw.out))  
cor.four.top25.raw.out<-paste(species,"_corr.csv",sep="", collapse = "")
write.csv(cor.four.top25.raw$r, (cor.four.top25.raw.out))  

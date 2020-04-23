## ---------------------------
##
## Script name: Main Script To View Network Rep of Patients
##
## Purpose of script: 
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 22-04-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

################################################################################

### Load Libs

################################################################################

library(igraph)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plyr)

################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# F to Add Column of Positive Sample Date if no collection date in range then shift date to last

addPosSpecDateCol = function(collections, pathways){
  pathways$PosSampleDT = rep(0,nrow(pathways))
  #Put a 1 in collections row for each patient when positive collect DT
  for (i in 1:nrow(collections)){
    name = collections$`Anonymised names`[i]
    pos.date = collections$`Specimen collection date`[i]
    index = intersect(which(pathways$`Anonymised names` == name),which(pathways$WardDate == pos.date))
    
    #Correct DT Date if not in range as last date in hospital
    if((length(index) != 1)){
      print(paste("Collection Date Not in Range",name))
      last.index = tail(which(pathways$`Anonymised names` == name), n=1)
      pathways$PosSampleDT[last.index] = 1
    }else{
      pathways$PosSampleDT[index] = 1
    }
  }
  return(pathways)
}

# ------------------------------------------------------------------------------

## F to Subset a patient dataset from the prior n days from positive sample

# x = dataframe of patient hisotry
# n.p = number of days prior
# l = lowerbound number of days
# c = TRUE/FALSE Depending on looking at confirmed vs possible HOCI

subset_stay = function(x,n.p = 14,l = 7,c = TRUE){
  
  # Get index of positve sample date
  pos.sample.index = which(x$PosSampleDT == 1)
  
  # Get index minus n.prior (or 1 if smaller than size of dataframe)
  proir.dat.index = max(1, (pos.sample.index - n.p))
  
  if(c){ #If Confirmed HOCI
    
    #If position less than 14 days return NULL
    if(pos.sample.index<n.p){
      return(NULL)
    }else{
      #Subset Data
      x = x[proir.dat.index:pos.sample.index,]
      return(x)
    }
  }else{ # asume possible HOCI
    
    #If less than the lower bound or not looking at COCI also return NULL
    if(pos.sample.index<l | pos.sample.index > (n.p-1)){
      return(NULL)
    }else{
      #Subset Data
      x = x[proir.dat.index:pos.sample.index,]
      return(x)
    }
  }
}


# ------------------------------------------------------------------------------

# F to get edges from a supplied dataframe subsetted by LoS

getAllEdges = function(patient.df){
  #Dataframe for edges
  edge.list <-list()
  # Loop through and look for change in ward name
  for(i in 1:(nrow(patient.df)-1)){
    ward.now = patient.df$CurrentLastWard[i]
    ward.next = patient.df$CurrentLastWard[i+1]
    if(ward.now != ward.next){
      #print(paste("Move",ward.now,"to",ward.next))
      df = data.frame(ward1=ward.now,
                      ward2=ward.next,
                      #ward1date=patient.df$WardDate[i],
                      #ward2date=patient.df$WardDate[i+1],
                      weight=1)
      edge.list[[(length(edge.list)+1)]] = df
    }
  }
  return(do.call(rbind, edge.list))
}



# ------------------------------------------------------------------------------

# F wrapper to get weighted edge df based on movement if hosptial stay longer than 14 days or more

# x = dataframe of patient hisotry
# n.p = number of days prior
# l = lowerbound number of days
# c = TRUE/FALSE Depending on looking at confirmed vs possible HOCI

get_wrapper_edge = function(x,n.p = 14,l = 7,c = TRUE){
  
  # Get index of positve sample date
  pos.sample.index = which(x$PosSampleDT == 1)
  
  # Get index minus n.prior (or 1 if smaller than size of dataframe)
  proir.dat.index = max(1, (pos.sample.index - n.p))
  
  if(c){ #If Confirmed HOCI
    
    #If position less than 14 days return NULL
    if(pos.sample.index<n.p){
      return(NULL)
    }else{
      #Subset Data
      x = x[proir.dat.index:pos.sample.index,]
      
      #get edges from sub df
      return(getAllEdges(x))
    }
  }else{ # asume possible HOCI
    
    #If less than the lower bound or not looking at COCI also return NULL
    if(pos.sample.index<l | pos.sample.index > (n.p-1)){
      return(NULL)
    }else{
      #Subset Data
      x = x[proir.dat.index:pos.sample.index,]
      
      #get edges from sub df
      return(getAllEdges(x))
    }
  }
}




################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


# ------------------------------------------------------------------------------

##
##  0. Preprocess Data
##

# 0.1 - Patient Pathways
path.raw <- read_excel("HOCIPath90PT_20200422.xlsx",
                       col_types = c("numeric", "text", "text", "date", "text", "text", "text", "text", "date"))

# Reset date columns
path.raw$AdmitDateTime = as.Date(path.raw$AdmitDateTime,format='%d/%m/%Y')
path.raw$WardDate = as.Date(path.raw$WardDate,format='%d/%m/%Y')


# 0.2 Patient Positive Collection Dates
col.dt <- read_excel("CollectionDates.xlsx", 
                     col_types = c("text", "numeric", "date", "text"))
col.dt$`Specimen collection date` = as.Date(col.dt$`Specimen collection date`,format='%d/%m/%Y')



# 0.3 Add column to pathways of psotive specimen date
path.raw = addPosSpecDateCol(col.dt,path.raw)

# 0.4 Order Date!!
path.raw = path.raw[order(as.Date(path.raw$WardDate, format="%d/%m/%Y")),]

# Split df by patient ID
patient.list <- split( path.raw , f = path.raw$`Anonymised names` )



# ------------------------------------------------------------------------------

##
##  1. Histomgrams and patients stays
##

day.seq = c(1,5,10,15,20,25)

count.list = lapply(day.seq, function(x){
  # Subset list for n days prior then combine into single dataframe
  sub.df = bind_rows(lapply(patient.list,subset_stay,x))
  # Create a coutn table 
  count.tab = as.data.frame(table(sub.df$CurrentLastWard))
  colnames(count.tab) = c("Ward",paste("Count",x,sep = ""))
  return(count.tab)
  
})



full.join.df = count.list %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Ward"), .)

full.join.df[is.na(full.join.df)] <- 0
full.join.df = cbind(full.join.df[,1:2],(full.join.df[,3:ncol(full.join.df)] - full.join.df[,2:(ncol(full.join.df)-1)]))

full.join.df.long = full.join.df %>% pivot_longer(-Ward, names_to = "Days_prior", values_to = "Count")

full.join.df.long$Days_prior = factor(full.join.df.long$Days_prior, levels = rev(paste("Count",day.seq,sep = "")))


# Make plot
ggplot(full.join.df.long, aes(x=Ward, y=Count,fill = Days_prior)) + 
  geom_bar(position="dodge",stat="identity")+
  theme_minimal()+
  labs(title = "Aggregate HOCI Covid Patient Days on Wards grouped by days prior (Patients w/ Los >= 14 days)",
       x = "", 
       y = "Aggregate Patient Ward Days")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))










# ------------------------------------------------------------------------------

##
##  2. Get Edge List 
##


# 2.1 Get Node Df
node.df = count.list %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Ward"), .)
node.df[is.na(node.df)] <- 0
#Take only total
node.df = node.df[,c(1,ncol(node.df))]

#Add column which is hosptial
map = unique(path.raw[,5:6])
colnames(map)[2] = "Ward"

node.df = left_join(node.df , map , by = "Ward")
colnames(node.df) = c("Ward","Size","Facility")

write.csv(node.df,"node.df.csv")

# Edge weights come from Movement Prior to Postive Sample Collection Date

# 2.2 First Method weights edges by 1 and only take edges within lookback window

# Subset list for n days prior then combine into single dataframe
meta.edges = bind_rows(lapply(patient.list,get_edge,25))

# Count duplicates assuming directed``````````````````````````````
meta.edges = count(meta.edges, vars = c("ward1","ward2"))
colnames(meta.edges) = c("source","target","weight")
write.csv(meta.edges,"edge.df.csv")


# 2.3 Second Method weights edges by a kernal based on distance from onset 







# ------------------------------------------------------------------------------

##
##  2. Plot Network
##










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
##  1.  Get Node Sizes
##

### Definintions

## ------------ Confirmed HOCI

# Confirmed HOCI was defined as any patient who had a COVID-19 positive sample 
# sent >14 days after admission. As the sample was taken outside the incubation 
# period from community-acquisition the most plausible explanation is healthcare-acquisition, 
# irrespective of symptoms at admission. 

## ------------ Possible HOCI 

# Possible HOCI was defined as any patient who had a COVID-19 positive sample 
# sent >7 days but < 14 days after admission AND did not have any symptoms of 
# COVID-19 (fever, cough, SOB, malaise) on admission. This adjusts for the 
# incubation period is up to 14 days, as it is still plausible that these 
# patients could have possibly acquired their infection in the community.


# 1.1 - Confirmed HOCI Covid Patient Data Movement
con.hoci.df = bind_rows(lapply(patient.list,subset_stay))
count.con.hoci = as.data.frame(table(con.hoci.df$CurrentLastWard))
colnames(count.con.hoci) = c("Ward",paste("Count"))

# 1.2 - Possible HOCI Covid Patient Data Movement
pos.hoci.df = bind_rows(lapply(patient.list,subset_stay,14,7,FALSE))
count.pos.hoci = as.data.frame(table(pos.hoci.df$CurrentLastWard))
colnames(count.pos.hoci) = c("Ward",paste("Count"))

# 1.3 - Join Two counts
count.con.pos = list(count.con.hoci,count.pos.hoci) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Ward"), .)
colnames(count.con.pos) = c("Ward","con_HOCI","pos_HOCI")


#Add column which is hosptial
map = unique(path.raw[,5:6])
colnames(map)[2] = "Ward"

count.con.hoci = left_join(count.con.hoci , map , by = "Ward")
colnames(count.con.hoci) = c("Ward","Size","Facility")

count.pos.hoci = left_join(count.pos.hoci , map , by = "Ward")
colnames(count.pos.hoci) = c("Ward","Size","Facility")

count.con.pos = left_join(count.con.pos , map , by = "Ward")
colnames(count.con.pos) = c("Ward","con_HOCI","pos_HOCI","Facility")




# 1.4 - Save as Node CSV's
write.csv(count.con.hoci,"confirmed_nodes.csv",row.names = F)
write.csv(count.pos.hoci,"possible_nodes.csv",row.names = F)
write.csv(count.con.pos,"combined_nodes.csv",row.names = F)



# ------------------------------------------------------------------------------

##
##  2.  Get Edges
##


# 2.1 - Confirmed HOCI Covid Patient Data Movement
con.hoci.edges = bind_rows(lapply(patient.list,get_wrapper_edge))
con.hoci.edges = count(con.hoci.edges, vars = c("ward1","ward2"))
colnames(con.hoci.edges) = c("source","target","weight")
con.hoci.edges$hoci = rep("Confirmed",nrow(con.hoci.edges))


# 2.2 - Possible HOCI Covid Patient Data Movement
pos.hoci.edges = bind_rows(lapply(patient.list,get_wrapper_edge,14,7,FALSE))
pos.hoci.edges = count(pos.hoci.edges, vars = c("ward1","ward2"))
colnames(pos.hoci.edges) = c("source","target","weight")
pos.hoci.edges$hoci = rep("Possible",nrow(pos.hoci.edges))


# 2.3 - Combine
combined.edges = rbind(con.hoci.edges,pos.hoci.edges)


# 2.4 - Save as Node CSV's
write.csv(con.hoci.edges,"confirmed_edges.csv",row.names = F)
write.csv(pos.hoci.edges,"possible_edges.csv",row.names = F)
write.csv(combined.edges,"combined_edges.csv",row.names = F)




# ------------------------------------------------------------------------------

##
##  3.  Create Node and Edge especially for analysis
##

networkx.nodes = count.con.pos
networkx.nodes$TotalSize = networkx.nodes$con_HOCI + networkx.nodes$pos_HOCI
write.csv(networkx.nodes,"networkx.nodes.csv",row.names = F)





# ------------------------------------------------------------------------------

##
##  4.
##











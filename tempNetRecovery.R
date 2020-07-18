## ---------------------------
##
## Script name: Temporal Network Analysis
##
## Purpose of script: 
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 04-06-2020
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

library(sna)
library(tsna)
library(ndtv)
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

## Function to check shared wards then call function to look for time intersect

# pateint.lists - list of patient dfs
# index_a - index of first patient to compare
# index_b - index of second patient to compare
# t - time difference for overlap
# mech = whether patients need to have same resitance mechanism

checkWards = function(patient.list,index_a,index_b,t=7,mech = F){
  
  # - Get wards which overlap
  wards = unique(patient.list[[index_a]]$CurrentLastWard[which(patient.list[[index_a]]$CurrentLastWard %in% patient.list[[index_b]]$CurrentLastWard)])
  
  # - For each overlapping ward run ward intersect function 
  cross.overs <- list()
  count = 1
  for(i in wards){
    sect = wardTimeIntersect(filter(patient.list[[index_a]], CurrentLastWard %in% i),
                             filter(patient.list[[index_b]], CurrentLastWard %in% i),
                             dif = t)
    
    # If cross over found save to list of df
    if(sect$Bool == TRUE){
      
      #Make df with connection 
      l_df = data.frame("Ward.cross" = rep(i,nrow(sect$Crosses)), 
                        "t_dif" = sect$Crosses$Time_dif, 
                        "date" = sect$Crosses$date) 
      cross.overs[[count]] = l_df
      count = count + 1
    }
  }
  
  # Check list of DF's is non empty 
  if(length(cross.overs) == 0){
    return(list("Bool" =  FALSE,"t"= NA))
  }else{
    return(list("Bool" =  TRUE,"df"= bind_rows(cross.overs)))
  }
}


# ------------------------------------------------------------------------------

## Function to check intersecting pathways of two patients for a single ward

# p_a - first patient
# p_b - second patient
# t - time difference 

wardTimeIntersect = function(p_a,p_b,dif=7){
  
  # For each row compare the days difference (only need to cycle through one patient and compare)
  df.list <- list()
  count = 1
  
  for(t in 1:nrow(p_a)){
    t_difs = as.numeric(p_a$WardDate[t] - p_b$WardDate)
    if(any(t_difs <= dif & t_difs >= dif)){
      t_diff = min(abs(t_difs))
      l_df = data.frame("Time_dif" = t_diff,
                        "date" = p_a$WardDate[t]+t_diff/1) # Date of a plus half the number of days between patients
      df.list[[count]] = l_df
      count = count + 1
    }
  }
  
  # Check list of DF's is non empty 
  if(length(df.list) == 0){
    return(list("Bool" =  FALSE,"Crosses"= NA))
  }else{
    return(list("Bool" = TRUE,"Crosses"= bind_rows(df.list)))
  }
}



################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


# ------------------------------------------------------------------------------

##
##  1. Recover Imp & Cleaned Pathways 
##

path.raw <- read.csv("Data/IMP_pathways10072020.csv")

#Remove redund cols
path.raw = path.raw[,1:13]

colnames(path.raw)[1] = "CollectDT"

#Clean colums
path.raw$case.ID = gsub(" ","_",path.raw$case.ID)

path.raw$cluster = gsub("cluster c","C",path.raw$cluster)
path.raw$cluster = gsub("c","C",path.raw$cluster)
path.raw$cluster = gsub("a","A",path.raw$cluster)
path.raw$cluster = gsub("b","B",path.raw$cluster)

#Replace missing columns with above if missing 
for(i in 1:nrow(path.raw)){
  if(path.raw[i,1] == ""){
    path.raw[i,1:6] = path.raw[i-1,1:6]
  }
}


# Reset date columns
path.raw$CollectDT = as.Date(path.raw$CollectDT,format='%d/%m/%Y')
path.raw$AdmDT = as.Date(path.raw$AdmDT,format='%d/%m/%Y')
path.raw$DisDT = as.Date(path.raw$DisDT,format='%d/%m/%Y')

# Temp fix to remove anons
path.raw = filter(path.raw, AdmDT <= DisDT)

# 0.5 Split df by patient ID
patient.list <- split( path.raw , f = path.raw$case.ID )

# Strange issue so remove missing element 
patient.list[[1]] = NULL


## Transform Patient Df into right long format

patient.list <- lapply(patient.list,function(x){
  if(nrow(x)==0){return(NA)}
  #Take each row and expand into df (lists)
  stay.list = lapply(1:nrow(x), function(x,df){
    #Select row of df
    row = df[x,]
    dates = seq.Date(row$AdmDT,row$DisDT,by=1)
    l = length(dates)
    df.row = row[rep(1,l),]
    #Drop Admit and leave date col
    df.row = df.row[,-c(9,10)]
    # Add properly formated dates
    df.row$WardDate = dates
    return(df.row)
  },df=x)
  df.stay = bind_rows(stay.list, .id = "Ward_stay")
  
  
  
  
  return(df.stay)
})


# Get csv with dates for when after PosT a patient moved to a side room
library(readr)
side_room <- read_csv("Data/Move_side_room.csv", 
                           col_types = cols(Comment = col_skip()))
# Remove NA's
side_room = side_room[-which(side_room$Last_date == "na"),]

#Fix Dates
side_room$Last_date = as.Date(side_room$Last_date,format='%d/%m/%Y')

patient.list.clean = lapply(patient.list, function(x,side_room){
  
  df = x
  df = x[order(as.Date(df$WardDate)),]
  
  #Check dupliacted dates
  if(any(duplicated(df$WardDate))){
    #print(paste("Dupliacted Ward Dates in",df$case.ID[1]))
  }
  
  # Add column to pathways of pos specimen date
  df$PosSampleDT = rep(0,nrow(df))
  df$PosSampleDT[which(df$WardDate == df$CollectDT)] = 1
  
  # If in set where patients were'nt moved to side room straight away move CollectDT down
  if(df$case.ID[1] %in% side_room$Case_no){
    side_idx = which(side_room$Case_no == df$case.ID[1])
    df$PosSampleDT[which(df$WardDate == side_room$Last_date[side_idx])] = 1
  }
  
  #If found then exit
  if(sum(df$PosSampleDT) > 0){
    return(df)
  }
  
  #Check it's not after the end
  if(df$CollectDT[1]>df$WardDate[nrow(df)]){
    #Add to end
    df$PosSampleDT[nrow(df)] = 1
    #If found then exit
    if(sum(df$PosSampleDT) == 1){
      return(df)
    }
  }
  
  # Check it's not in the middle but patient wasnt on a ward 
  if((df$CollectDT[1] > min(df$WardDate))&((df$CollectDT[1] < max(df$WardDate)))){
    #Set PosSample DT to least negtive index from this vector (will be included in pathway analysis but anything positive i.e after will be excluded)
    df$PosSampleDT = df$WardDate - df$CollectDT[1]
    index = which(df$PosSampleDT == max((filter(df, PosSampleDT < 0))$PosSampleDT))
    df$PosSampleDT = rep(0,nrow(df))
    df$PosSampleDT[index] = 1
    #return(df)
    #If found then exit
    if(sum(df$PosSampleDT) > 0){
      return(df)
    }
  }
  
  #Check it's not before start of pathways - if so then return NA
  if(df$CollectDT[1]<df$WardDate[1]){
    print(paste("PosSampleDT before pathway",df$case.ID[1]))
    return(NA)
  }
  
  if(sum(df$PosSampleDT) > 0){
    print(paste("Missing PosSampleDT",df$case.ID[1]))
  }
},side_room=side_room)


# Remove Case_66 - CollectionDT before start of pathway history
patient.list.clean[["Case_65"]] = NULL


# Remove Data after PosSampleDT
patient.list.clean = lapply(patient.list.clean, function(x){
  if(is.na(x)){return(NULL)}
  index = which(x$PosSampleDT == 1)
  x = x[1:max(index),] # <------------------------- Max because some overlaps mean mutliple wards can hposave the pos
  return(x)
})



# ------------------------------------------------------------------------------

##
##  Extract Temporal Edges
##



insec_m = matrix(rep(0, len=(length(patient.list.clean))^2), nrow = length(patient.list.clean))
row.names(insec_m) = names(patient.list.clean)
colnames(insec_m) = names(patient.list.clean)
insec_edge_l <- list()
count = 1

for (i in 1:length(patient.list.clean)) {
  for (j in 1:length(patient.list.clean)) {
    if(i != j){
      if(any(patient.list.clean[[i]]$CurrentLastWard %in% patient.list[[j]]$CurrentLastWard)){
        isec = checkWards(patient.list.clean,i,j,t=0,F)
        
        if(isec$Bool){
          print(paste(row.names(insec_m)[i],row.names(insec_m)[j]))
          

          l_df = data.frame("source"=rep(row.names(insec_m)[i],nrow(isec$df)),
                            "target"=rep(row.names(insec_m)[j],nrow(isec$df)))
          
          l_df = cbind(l_df,isec$df)
          
          
          insec_edge_l [[count]] = l_df
          count = count + 1
        }
      }
    }
  }
}


# Change dates to number 
date.map = data.frame("date"=seq(as.Date(min(bind_rows(patient.list.clean)$WardDate)),as.Date(max(bind_rows(patient.list.clean)$WardDate)),"day"),
                      "t" = 1:length(seq(as.Date(min(bind_rows(patient.list.clean)$WardDate)),as.Date(max(bind_rows(patient.list.clean)$WardDate)),"day")))


# Pull out dynamic edges indexed by date numebr
dyn.edges = bind_rows(lapply(insec_edge_l, function(x,d.mapping){
  temp = x
  
  # Join extra col of date index
  temp = left_join(temp, d.mapping, by = c('date'))
  #Remove duplicated rows --- could be from data entry errors
  temp = temp[!duplicated(temp$t), ]
  #Sort edges by time
  temp = temp[order(temp$t),]
  # Aggregate edges (if any duplicates)
  temp$agg = c(0,abs(diff(temp$t)))
  temp$agg_group = rep(0,nrow(temp))
  mark = 0
  for (i in 1:nrow(temp)){
    if(temp$agg[i] != 1){
      mark = mark + 1
    }
    temp$agg_group[i] = mark
  }
  
  edge.groups <- split(temp, f = temp$agg_group)
  edges = bind_rows(lapply(edge.groups, function(x){
    #start
    s = min(x$t)
    #finish
    f = max(x$t)
    #duraion
    d = f - s
    edge = data.frame("source"=x$source[1],"target"=x$target[1],"ward"=x$Ward.cross[1],
                      "start"=s,"finish"=f,"duration"=d)
    return(edge)
  }))
  
  return(edges)
  
},d.mapping = date.map))



## Filter out dupliacted edges

# Duplaicted have an edge which is connecting the same wards and starts at the same time

dyn.edges = dyn.edges[!duplicated(dyn.edges[,1:4]),]




## Map across by case ID to index

# Map cases to numbers
nodes = bind_rows(lapply(patient.list.clean, function(x){return(x[1,3:7])}))
nodes$vertex.id = 1:nrow(nodes)

nodes$source = nodes$case.ID
nodes$target = nodes$case.ID

dyn.edges = left_join(dyn.edges, nodes[,c(6,7)], by = c('source'))
dyn.edges = dyn.edges[,-1]
colnames(dyn.edges)[ncol(dyn.edges)] = "source"

dyn.edges = left_join(dyn.edges, nodes[,c(6,8)], by = c('target'))
dyn.edges = dyn.edges[,-1]
colnames(dyn.edges)[ncol(dyn.edges)] = "target"

nodes= nodes[,-c(7,8)]



## Get dynamic vertex

patient.list.clean=patient.list.clean[-(which(sapply(patient.list.clean,is.null),arr.ind=TRUE))]

dyn.nodes = bind_rows(lapply(patient.list.clean, function(x){
  df = data.frame("case.ID" = (x$case.ID[1]),
                  "start" = min(x$WardDate),
                  "end" = max(x$WardDate))
  return(df)
}))


#map dates
dyn.nodes = left_join(dyn.nodes, nodes[,c(1,6)], by = c('case.ID'))
dyn.nodes = dyn.nodes[,-1]

date.map$start = date.map$date
date.map$end = date.map$date

dyn.nodes = left_join(dyn.nodes, date.map[,2:3], by = c('start'))
dyn.nodes = dyn.nodes[,-1]

dyn.nodes = left_join(dyn.nodes, date.map[,c(2,4)], by = c('end'))
dyn.nodes = dyn.nodes[,-1]

colnames(dyn.nodes)[2:3] = c("start","end")
dyn.nodes$duration = dyn.nodes$end - dyn.nodes$start

##
##    Format Data for dynamic network analysis
##

## Edges
dynamic.edges = data.frame("onset" = dyn.edges$start, #start time
                           "terminus" = dyn.edges$finish, # end time
                           "tail" = dyn.edges$source, # source
                           "head" = dyn.edges$target, # target
                           "onset.censored" = rep(FALSE,nrow(dyn.edges)), # ignore edge onset
                           "terminus.censored" = rep(FALSE,nrow(dyn.edges)), # ignore edge terminys
                           "duration" = dyn.edges$duration,
                           "edge.id" = 1:nrow(dyn.edges))

## Dynamic nodes
dynamic.nodes = data.frame("onset" = dyn.nodes$start,
                           "terminus" = dyn.nodes$end,
                           "vertex.id" = dyn.nodes$vertex.id,
                           "onset.censored" = rep(FALSE,nrow(dyn.nodes)),
                           "terminus.censored"= rep(FALSE,nrow(dyn.nodes)),
                           "duration" = dyn.nodes$duration)



## Vertex attributes
vert.at = data.frame("vertex.id" = nodes$vertex.id,"case.ID"=nodes$case.ID,
                     "cluster" =nodes$cluster,"Plasmid" = nodes$Plasmid)

## Static edges
stat.edges = data.frame("tail"= dyn.edges$source,"head"=dyn.edges$target)
stat.edges = stat.edges[!duplicated(stat.edges), ]








write.csv(dynamic.edges,"dynamic.edges.csv",row.names = F)
write.csv(dyn.nodes,"dyn.nodes.csv",row.names = F)
write.csv(stat.edges,"static.edges.csv",row.names = F)
write.csv(vert.at,"vert.attributes.csv",row.names = F)








# ------------------------------------------------------------------------------

# Analysis of intersections

ward.insecs = as.data.frame(table(dyn.edges$ward))

library(ggplot2)
p <- ggplot(ward.insecs , aes(x = Var1, y = Freq))+
  geom_col( width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  ggtitle("CPE Imp Pathway intersection locations") +
  xlab("Ward") + ylab("Pathway intersections")
p



